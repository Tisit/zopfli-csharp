using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Runtime.CompilerServices;
using System.Runtime.InteropServices;
using ZopfliCSharp;

/*
Streaming output sink that mirrors the subset of List<byte> operations the
deflate encoder uses (Add, AddRange, Count, and indexing the last byte while it
is being bit-packed), but flushes completed bytes to an underlying Stream after
each master block instead of holding the whole compressed result in memory.

The bit writers only ever mutate the most recently added byte (index Count-1);
once a byte is complete and a new one is appended, the old one is never touched
again. That invariant is what makes it safe to flush everything except the last
(possibly still-being-packed) byte.
*/
sealed class OutputSink
{
    readonly Stream _out;
    readonly List<byte> _buf;   /* Bytes not yet written to _out. */
    long _flushed;              /* Count of bytes already written to _out. */

    public OutputSink(Stream outStream)
    {
        _out = outStream;
        _buf = new List<byte>();
        _flushed = 0;
    }

    /* Total bytes logically written so far (flushed + buffered). Callers use this
       as a running offset to measure the emitted size of a block. */
    public int Count => (int)(_flushed + _buf.Count);

    public void Add(byte b) => _buf.Add(b);

    public void AddRange(IEnumerable<byte> items) => _buf.AddRange(items);

    /* Only ever indexed at Count-1, the byte currently being bit-packed, which is
       always still in the buffer because FlushCompleted keeps the last byte. */
    public byte this[int index]
    {
        get => _buf[index - (int)_flushed];
        set => _buf[index - (int)_flushed] = value;
    }

    /* Write out every completed byte, keeping only the final byte (which may still
       receive more bits in the next block) in memory. */
    public void FlushCompleted()
    {
        int n = _buf.Count - 1;   /* Keep the last byte. */
        if (n <= 0) return;
        _out.Write(CollectionsMarshal.AsSpan(_buf).Slice(0, n));
        _buf.RemoveRange(0, n);
        _flushed += n;
    }

    /* Write everything remaining and flush the underlying stream. */
    public void Finish()
    {
        if (_buf.Count > 0)
        {
            _out.Write(CollectionsMarshal.AsSpan(_buf));
            _flushed += _buf.Count;
            _buf.Clear();
        }
        _out.Flush();
    }
}

class ZopfliLZ77Store
{
    /* The per-symbol arrays below are parallel and have logical length `size`.
       Their allocated capacity may be larger; grow via EnsureMainCapacity. */
    public ushort[] litlens;  /* Lit or len. */
    public ushort[] dists;  /* If 0: indicates literal in corresponding litlens,
      if > 0: length in corresponding litlens, this is the distance. */
    public ulong size;

    public byte[] data;  /* original data */
    public ulong[] pos;  /* position in data where this LZ77 command begins */

    public ushort[] ll_symbol;
    public ushort[] d_symbol;

    /* Cumulative histograms wrapping around per chunk. Each chunk has the amount
    of distinct symbols as length, so using 1 value per LZ77 symbol, we have a
    precise histogram at every N symbols, and the rest can be calculated by
    looping through the actual symbols of this chunk. Logical lengths are
    ll_counts_size / d_counts_size. */
    public ulong[] ll_counts;
    public ulong[] d_counts;
    public int ll_counts_size;
    public int d_counts_size;

    //ZopfliInitLZ77Store() from original
    public ZopfliLZ77Store(byte[] indata)
    {
        SetDefaults(indata);
    }
    public void ResetStore(byte[] indata)
    {
        SetDefaults(indata);
    }
    void SetDefaults(byte[] indata)
    {
        const int cap = 16;
        litlens = new ushort[cap];
        dists = new ushort[cap];
        size = 0;
        data = indata;
        pos = new ulong[cap];
        ll_symbol = new ushort[cap];
        d_symbol = new ushort[cap];
        ll_counts = Array.Empty<ulong>();
        d_counts = Array.Empty<ulong>();
        ll_counts_size = 0;
        d_counts_size = 0;
    }

    /* Ensures the parallel per-symbol arrays hold at least `needed` elements. */
    public void EnsureMainCapacity(int needed)
    {
        if (litlens.Length >= needed) return;
        int cap = litlens.Length * 2;
        if (cap < needed) cap = needed;
        Array.Resize(ref litlens, cap);
        Array.Resize(ref dists, cap);
        Array.Resize(ref pos, cap);
        Array.Resize(ref ll_symbol, cap);
        Array.Resize(ref d_symbol, cap);
    }

    /* Grows the histogram arrays with doubling (amortized O(1)); a fixed-step
       Array.Resize per chunk would be quadratic over the whole store. */
    public void EnsureLLCountsCapacity(int needed)
    {
        if (ll_counts.Length >= needed) return;
        int cap = ll_counts.Length == 0 ? needed : ll_counts.Length * 2;
        if (cap < needed) cap = needed;
        Array.Resize(ref ll_counts, cap);
    }

    public void EnsureDCountsCapacity(int needed)
    {
        if (d_counts.Length >= needed) return;
        int cap = d_counts.Length == 0 ? needed : d_counts.Length * 2;
        if (cap < needed) cap = needed;
        Array.Resize(ref d_counts, cap);
    }
}

/*
A longest match for some position in run-length form: for each length in
(end[r - 1], end[r]] (starting at ZOPFLI_MIN_MATCH for r = 0), dist[r] is the
smallest distance that reaches that length. This replaces the per-length "sublen"
array of the C original. The match finder visits candidates in order of increasing
distance, so each time it finds a longer match that starts a new run, and runs are
never adjacent with equal distances.
*/
class MatchRuns
{
    public readonly int[] end = new int[Compress.ZOPFLI_MAX_MATCH + 1];
    public readonly ushort[] dist = new ushort[Compress.ZOPFLI_MAX_MATCH + 1];
    public int count;
}

class ZopfliLongestMatchCache
{
    public ushort[] length;
    public ushort[] dist;
    public byte[] runs;  /* Per position, the runs of the match (see ZopfliRunsToCache). */
}
class ZopfliBlockState
{
    /*
    For longest match cache. max 256. Uses huge amounts of memory but makes it
    faster. Uses this many times three bytes per single byte of the input data.
    This is so because longest match finding has to find the exact distance
    that belongs to each length for the best lz77 strategy.
    Good values: e.g. 5, 8.
    */
    const int ZOPFLI_CACHE_LENGTH = 8;

    /* Cache for length/distance pairs found so far. */
    public ZopfliLongestMatchCache lmc;

    /* The start (inclusive) and end (not inclusive) of the current block. */
    public int blockstart;
    public int blockend;

    // does ZopfliInitBlockState
    public ZopfliBlockState(int istart, int iend, int add_lmc)
    {
        blockstart = istart;
        blockend = iend;

        int blocksize = blockend - blockstart;

        if (add_lmc > 0)
        {
            lmc = new ZopfliLongestMatchCache();
            lmc.length = new ushort[blocksize];
            for (int i = 0; i < blocksize; i++)
            {
                lmc.length[i] = 1;
            }
            lmc.dist = new ushort[blocksize];
            lmc.runs = new byte[ZOPFLI_CACHE_LENGTH * 3 * blocksize];
        }
        else
        {
            lmc = new ZopfliLongestMatchCache();
        }
    }

    /*
    Stores the runs of a longest match found for pos, whose length is length, in the
    cache. Holds the first ZOPFLI_CACHE_LENGTH runs as (last length - 3, distance)
    triples; if there are fewer, the last slot's length records the longest cached
    length.
    */
    public void ZopfliRunsToCache(MatchRuns runs, int pos, int length)
    {
        byte[] cs = lmc.runs;
        int cachebase = ZOPFLI_CACHE_LENGTH * pos * 3;

        if (length < 3) return;
        int n = Math.Min(runs.count, ZOPFLI_CACHE_LENGTH);
        for (int j = 0; j < n; j++)
        {
            int end = runs.end[j];
            ushort dist = runs.dist[j];
            cs[cachebase + j * 3] = (byte)(end - 3);
            cs[cachebase + j * 3 + 1] = (byte)(dist % 256);
            cs[cachebase + j * 3 + 2] = (byte)((dist >> 8) % 256);
        }
        if (n < ZOPFLI_CACHE_LENGTH)
        {
            Debug.Assert(runs.end[n - 1] == length);
            cs[cachebase + (ZOPFLI_CACHE_LENGTH - 1) * 3] = (byte)(length - 3);
        }
        Debug.Assert(runs.end[n - 1] <= length);
        Debug.Assert(runs.end[n - 1] == ZopfliMaxCachedLength(pos));
    }

    /*
    Loads the cached runs for pos into runs (see ZopfliRunsToCache). Only valid when
    runs are cached for pos, i.e. ZopfliMaxCachedLength(pos) > 0.
    */
    public void ZopfliCacheToRuns(int pos, MatchRuns runs)
    {
        int maxlength = ZopfliMaxCachedLength(pos);
        int cachebase = ZOPFLI_CACHE_LENGTH * pos * 3;
        byte[] cs = lmc.runs;
        int n = 0;
        for (int j = 0; j < ZOPFLI_CACHE_LENGTH; j++)
        {
            int o = cachebase + j * 3;
            int length2 = cs[o] + 3;
            runs.end[n] = length2;
            runs.dist[n] = (ushort)(cs[o + 1] + 256 * cs[o + 2]);
            n++;
            if (length2 == maxlength) break;
        }
        runs.count = n;
    }

    /*
    The longest match cache lookup of the optimal parser's hot loop (limit ==
    ZOPFLI_MAX_MATCH, runs wanted), with ZopfliMaxCachedLength and ZopfliCacheToRuns
    folded in. Hits exactly when Compress.TryGetFromLongestMatchCache would. Returns
    whether it hit; if so, sets length and runs.
    The duplicated decoding is deliberate: sharing one helper with ZopfliCacheToRuns
    (inlined or not) measured ~1.7% slower overall.
    */
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public bool TryGetCachedRuns(int pos, MatchRuns runs, ref ushort length)
    {
        ushort[] lmcLength = lmc.length;
        if (lmcLength == null) return false;  /* No cache (add_lmc == 0). */

        int lmcpos = pos - blockstart;
        ushort cachedLen = lmcLength[lmcpos];
        /* Length > 0 and dist 0 is invalid combination, which indicates on purpose
           that this cache value is not filled in yet. */
        if (cachedLen != 0 && lmc.dist[lmcpos] == 0) return false;

        /* One bounds check for the position's whole cache entry. */
        ReadOnlySpan<byte> entry = lmc.runs.AsSpan(ZOPFLI_CACHE_LENGTH * 3 * lmcpos,
                                                     ZOPFLI_CACHE_LENGTH * 3);
        int maxlength = entry[1] == 0 && entry[2] == 0
            ? 0 : entry[(ZOPFLI_CACHE_LENGTH - 1) * 3] + 3;
        if (cachedLen > maxlength) return false;

        length = cachedLen;
        int n = 0;
        if (cachedLen >= 3)
        {
            int[] end = runs.end;
            ushort[] dist = runs.dist;
            for (int j = 0; j < ZOPFLI_CACHE_LENGTH; j++)
            {
                int length2 = entry[j * 3] + 3;
                end[n] = length2;
                dist[n] = (ushort)(entry[j * 3 + 1] + 256 * entry[j * 3 + 2]);
                n++;
                if (length2 == maxlength) break;
            }
        }
        runs.count = n;
        return true;
    }

    /*
    Returns the length up to which could be stored in the cache.
    */
    public int ZopfliMaxCachedLength(int pos)
    {
        byte[] cs = lmc.runs;
        int cachebase = ZOPFLI_CACHE_LENGTH * pos * 3;
        if (cs[cachebase + 1] == 0 && cs[cachebase + 2] == 0) return 0;  /* No runs cached. */
        return cs[cachebase + (ZOPFLI_CACHE_LENGTH - 1) * 3] + 3;
    }
}


class ZopfliHash
{
    /*
The window size for deflate. Must be a power of two. This should be 32768, the
maximum possible by the deflate spec. Anything less hurts compression more than
speed.
*/
    public const int ZOPFLI_WINDOW_SIZE = 32768;
    const int HASH_SHIFT = 5;
    const int HASH_MASK = 32767;
    public const int ZOPFLI_WINDOW_MASK = ZOPFLI_WINDOW_SIZE - 1;
    const int ZOPFLI_MIN_MATCH = 3;

    /* Hash values are masked to 15 bits (HASH_MASK), and positions to the window, so
       the head tables need only HASH_MASK + 1 short entries (64 KB each), which keeps
       both within the L2 cache. */
    public short[] head;  /* Hash value to index of its most recent occurrence. */
    public short[] prev;  /* Index to index of prev. occurrence of same hash. */
    public short[] hashval;  /* Index to hash value at this index. */
    public int val;  /* Current hash value. */

    /* Fields with similar purpose as the above hash, but for the second hash with
    a value that is calculated differently.  */
    public short[] head2;  /* Hash value to index of its most recent occurrence. */
    public short[] prev2;  /* Index to index of prev. occurrence of same hash. */
    public short[] hashval2;  /* Index to hash value at this index. */
    public int val2;  /* Current hash value. */

    public ushort[] same;  /* Amount of repetitions of same byte after this .*/

    public ZopfliHash()
    {
        head = new short[HASH_MASK + 1];
        prev = new short[ZOPFLI_WINDOW_SIZE];
        hashval = new short[ZOPFLI_WINDOW_SIZE];

        same = new ushort[ZOPFLI_WINDOW_SIZE];

        head2 = new short[HASH_MASK + 1];
        prev2 = new short[ZOPFLI_WINDOW_SIZE];
        hashval2 = new short[ZOPFLI_WINDOW_SIZE];
    }

    public void ZopfliResetHash()
    {
        int i;

        Array.Fill<short>(head, -1 /* -1 indicates no head so far. */);

        for (i = 0; i < ZOPFLI_WINDOW_SIZE; i++)
        {
            prev[i] = (short)i;  /* If prev[j] == j, then prev[j] is uninitialized. */
            prev2[i] = (short)i;
        }

        Array.Fill<short>(hashval, -1 /* -1 indicates no head so far. */);
        Array.Clear(same, 0, same.Length);

        val = 0;
        val2 = 0;
        Array.Fill<short>(head2, -1 /* -1 indicates no head so far. */);
        Array.Fill<short>(hashval2, -1 /* -1 indicates no head so far. */);
    }
    /*
    Update the sliding hash value with the given byte. All calls to this function
    must be made on consecutive input characters. Since the hash value exists out
    of multiple input bytes, a few warmups with this function are needed initially.
    */
    void UpdateHashValue(byte c)
    {
        val = (((val) << HASH_SHIFT) ^ (c)) & HASH_MASK;
    }

    public void ZopfliUpdateHash(byte[] array, int pos, int end)
    {
        int hpos = pos & ZOPFLI_WINDOW_MASK;
        int amount = 0;

        byte t;
        if (pos + ZOPFLI_MIN_MATCH <= end)
        {
            t = array[pos + ZOPFLI_MIN_MATCH - 1];
        } else
            t = 0;
        UpdateHashValue(t);
        int v = val;
        hashval[hpos] = (short)v;
        int h = head[v];
        prev[hpos] = (short)(h != -1 && hashval[h] == v ? h : hpos);
        head[v] = (short)hpos;

        /* Update "same". */
        int prevsame = same[(pos - 1) & ZOPFLI_WINDOW_MASK];
        if (prevsame > 1)
        {
            amount = prevsame - 1;
        }
        while (pos + amount + 1 < end && array[pos] == array[pos + amount + 1] && amount < ushort.MaxValue) {
            amount++;
        }
        same[hpos] = (ushort)amount;

        int v2 = ((amount - ZOPFLI_MIN_MATCH) & 255) ^ v;
        val2 = v2;
        hashval2[hpos] = (short)v2;
        int h2 = head2[v2];
        prev2[hpos] = (short)(h2 != -1 && hashval2[h2] == v2 ? h2 : hpos);
        head2[v2] = (short)hpos;
    }
    public void ZopfliWarmupHash(byte[] InFile, int pos, int end)
    {
        UpdateHashValue(InFile[pos + 0]);
        if (pos + 1 < end) UpdateHashValue(InFile[pos + 1]);
    }
}

/*
The chains of a ZopfliHash, precomputed for every position of a range. The rolling
hash is deterministic, so for a given range the chains seen at a position are the
same in every pass over it, and the optimal parser makes two passes per iteration
(~30 per block). Building the hash once and recording, per position, everything the
longest match search reads from it replaces all that hash maintenance with lookups.
The arrays are indexed by (absolute position - start), where start is the window
start of the range, i.e. where the hash would have been reset and warmed up.
*/
class ZopfliHashChains
{
    public readonly int start;

    /* Distance back to the previous position in this position's chain of the first
       hash, or 0 if it is the end of its chain (ZopfliHash's prev[i] == i). */
    public readonly ushort[] step;
    /* Same, for the second hash. */
    public readonly ushort[] step2;
    /* Amount of repetitions of same byte after this. */
    public readonly ushort[] same;
    /* Value of the second hash at this position. */
    public readonly ushort[] val2;

    /* Runs ZopfliHash over [instart - window, inend) exactly as a pass that searches
       matches in [instart, inend) would. */
    public ZopfliHashChains(byte[] array, int instart, int inend)
    {
        const int mask = ZopfliHash.ZOPFLI_WINDOW_MASK;
        start = instart > ZopfliHash.ZOPFLI_WINDOW_SIZE
            ? instart - ZopfliHash.ZOPFLI_WINDOW_SIZE : 0;
        int n = Math.Max(inend - start, 0);
        step = new ushort[n];
        step2 = new ushort[n];
        same = new ushort[n];
        val2 = new ushort[n];
        if (instart >= inend) return;

        ZopfliHash h = new ZopfliHash();
        h.ZopfliResetHash();
        h.ZopfliWarmupHash(array, start, inend);
        for (int i = start; i < inend; i++)
        {
            h.ZopfliUpdateHash(array, i, inend);
            int hpos = i & mask;
            int k = i - start;
            step[k] = (ushort)((hpos - h.prev[hpos]) & mask);
            step2[k] = (ushort)((hpos - h.prev2[hpos]) & mask);
            same[k] = h.same[hpos];
            val2[k] = (ushort)h.val2;
        }
    }
}

class SplitCostContext
{
    public ZopfliLZ77Store lz77;
    public ulong start;
    public ulong end;
}

class SymbolStats
{
    /* The literal and length symbols. */
    public ulong[] litlens = new ulong[Compress.ZOPFLI_NUM_LL];
    /* The 32 unique dist symbols, not the 32768 possible dists. */
    public ulong[] dists = new ulong[Compress.ZOPFLI_NUM_D];

    /* Length of each lit/len symbol in bits. */
    public double[] ll_symbols = new double[Compress.ZOPFLI_NUM_LL];
    /* Length of each dist symbol in bits. */
    public double[] d_symbols = new double[Compress.ZOPFLI_NUM_D];

    void ZopfliCalculateEntropy(ulong[] count, int n, double[] bitlengths)
    {
        const double kInvLog2 = 1.4426950408889;  /* 1.0 / log(2.0) */
        uint sum = 0;
        uint i;
        double log2sum;
        for (i = 0; i < n; ++i)
        {
            sum += (uint)count[i];
        }
        log2sum = (sum == 0 ? Math.Log(n) : Math.Log(sum)) * kInvLog2;
        for (i = 0; i < n; ++i)
        {
            /* When the count of the symbol is 0, but its cost is requested anyway, it
            means the symbol will appear at least once anyway, so give it the cost as if
            its count is 1.*/
            if (count[i] == 0) bitlengths[i] = log2sum;
            else bitlengths[i] = log2sum - Math.Log(count[i]) * kInvLog2;
            /* Depending on compiler and architecture, the above subtraction of two
            floating point numbers may give a negative result very close to zero
            instead of zero (e.g. -5.973954e-17 with gcc 4.1.2 on Ubuntu 11.4). Clamp
            it to zero. These floating point imprecisions do not affect the cost model
            significantly so this is ok. */
            if (bitlengths[i] < 0 && bitlengths[i] > -1e-5) bitlengths[i] = 0;
            Debug.Assert(bitlengths[i] >= 0);
        }
    }

    /* Calculates the entropy of the statistics */
    public void CalculateStatistics()
    {
        ZopfliCalculateEntropy(litlens, Compress.ZOPFLI_NUM_LL, ll_symbols);
        ZopfliCalculateEntropy(dists, Compress.ZOPFLI_NUM_D, d_symbols);
    }
}

class RanState
{
    public uint m_w = 1;
    public uint m_z = 2;
}