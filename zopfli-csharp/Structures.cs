using System;
using System.Collections.Generic;
using System.Diagnostics;
using ZopfliCSharp;

class ZopfliLZ77Store
{
    public List<ushort> litlens;  /* Lit or len. */
    public List<ushort> dists;  /* If 0: indicates literal in corresponding litlens,
      if > 0: length in corresponding litlens, this is the distance. */
    public ulong size;

    public byte[] data;  /* original data */
    public List<ulong> pos;  /* position in data where this LZ77 command begins */

    public List<ushort> ll_symbol;
    public List<ushort> d_symbol;

    /* Cumulative histograms wrapping around per chunk. Each chunk has the amount
    of distinct symbols as length, so using 1 value per LZ77 symbol, we have a
    precise histogram at every N symbols, and the rest can be calculated by
    looping through the actual symbols of this chunk. */
    public List<ulong> ll_counts;
    public List<ulong> d_counts;

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
        litlens = new List<ushort>();
        dists = new List<ushort>();
        size = 0;
        data = indata;
        pos = new List<ulong>();
        ll_symbol = new List<ushort>();
        d_symbol = new List<ushort>();
        ll_counts = new List<ulong>();
        d_counts = new List<ulong>();
    }
}

class ZopfliLongestMatchCache
{
    public ushort[] length;
    public ushort[] dist;
    public byte[] sublen;
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
            lmc.sublen = new byte[ZOPFLI_CACHE_LENGTH * 3 * blocksize];
        }
        else
        {
            lmc = new ZopfliLongestMatchCache();
        }
    }

    public void ZopfliSublenToCache(ushort[] sublen, int pos, int length)
    {
        int i;
        int j = 0;
        uint bestlength = 0;
        
        if (length < 3) return;
        for (i = 3; i <= length; i++)
        {
            if (i == length || sublen[i] != sublen[i + 1])
            {
                lmc.sublen[ZOPFLI_CACHE_LENGTH * pos * 3 + j * 3] = (byte)(i - 3);
                lmc.sublen[ZOPFLI_CACHE_LENGTH * pos * 3 + j * 3 + 1] = (byte)(sublen[i] % 256);
                lmc.sublen[ZOPFLI_CACHE_LENGTH * pos * 3 + j * 3 + 2] = (byte)((sublen[i] >> 8) % 256);
                bestlength = (uint)i;
                j++;
                if (j >= ZOPFLI_CACHE_LENGTH) break;
            }
        }
        if (j < ZOPFLI_CACHE_LENGTH)
        {
            Debug.Assert(bestlength == length);
            lmc.sublen[ZOPFLI_CACHE_LENGTH * pos * 3 + (ZOPFLI_CACHE_LENGTH - 1) * 3] = (byte)(bestlength - 3);
        }
        else
        {
            Debug.Assert(bestlength <= length);
        }
        Debug.Assert(bestlength == ZopfliMaxCachedSublen((ulong)pos));
        
    }

    public void ZopfliCacheToSublen(
                         ulong pos, ulong length,
                         ushort[] sublen)
    {
        ulong i, j;
        uint maxlength = ZopfliMaxCachedSublen(pos);
        int prevlength = 0;

        if (length < 3) return;
        for (j = 0; j < ZOPFLI_CACHE_LENGTH; j++)
        {
            int length2 = lmc.sublen[ZOPFLI_CACHE_LENGTH * pos * 3 + j * 3] + 3;
            int dist = lmc.sublen[ZOPFLI_CACHE_LENGTH * pos * 3 + j * 3 + 1] + 256 * lmc.sublen[ZOPFLI_CACHE_LENGTH * pos * 3 + j * 3 + 2];
            for (i = (ulong)prevlength; i <= (ulong)length2; i++)
            {
                sublen[i] = (ushort)dist;
            }
            if (length2 == maxlength) break;
            prevlength = length2 + 1;
        }
    }

    /*
    Returns the length up to which could be stored in the cache.
    */
    public uint ZopfliMaxCachedSublen(ulong pos)
    {
        if (lmc.sublen[ZOPFLI_CACHE_LENGTH * pos * 3 + 1] == 0 && lmc.sublen[ZOPFLI_CACHE_LENGTH * pos * 3 + 2] == 0) return 0;  /* No sublen cached. */
        return (uint)(lmc.sublen[ZOPFLI_CACHE_LENGTH * pos * 3 + (ZOPFLI_CACHE_LENGTH - 1) * 3] + 3);
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

    public int[] head;  /* Hash value to index of its most recent occurrence. */
    public short[] prev;  /* Index to index of prev. occurrence of same hash. */
    public short[] hashval;  /* Index to hash value at this index. */
    public int val;  /* Current hash value. */

    /* Fields with similar purpose as the above hash, but for the second hash with
    a value that is calculated differently.  */
    public int[] head2;  /* Hash value to index of its most recent occurrence. */
    public short[] prev2;  /* Index to index of prev. occurrence of same hash. */
    public short[] hashval2;  /* Index to hash value at this index. */
    public int val2;  /* Current hash value. */

    public ushort[] same;  /* Amount of repetitions of same byte after this .*/

    public ZopfliHash()
    {
        head = new int[65536];
        prev = new short[ZOPFLI_WINDOW_SIZE];
        hashval = new short[ZOPFLI_WINDOW_SIZE];

        same = new ushort[ZOPFLI_WINDOW_SIZE];

        head2 = new int[65536];
        prev2 = new short[ZOPFLI_WINDOW_SIZE];
        hashval2 = new short[ZOPFLI_WINDOW_SIZE];
    }

    public void ZopfliResetHash()
    {
        int i;

        Array.Fill(head, -1 /* -1 indicates no head so far. */);

        for (i = 0; i < ZOPFLI_WINDOW_SIZE; i++)
        {
            prev[i] = (short)i;  /* If prev[j] == j, then prev[j] is uninitialized. */
            prev2[i] = (short)i;
        }

        Array.Fill<short>(hashval, -1 /* -1 indicates no head so far. */);
        same.Initialize();

        val = 0;
        val2 = 0;
        Array.Fill(head2, -1 /* -1 indicates no head so far. */);
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
        short hpos = (short)(pos & ZOPFLI_WINDOW_MASK);
        int amount = 0;

        byte t;
        if (pos + ZOPFLI_MIN_MATCH <= end)
        {
            t = array[pos + ZOPFLI_MIN_MATCH - 1];
        } else
            t = 0;
        UpdateHashValue(t);
        hashval[hpos] = (short)val;
        if (head[val] != -1 && hashval[head[val]] == val)
        {
            prev[hpos] = (short)head[val];
            Debug.Assert(head[val] <= short.MaxValue);
        }
        else prev[hpos] = hpos;
        head[val] = hpos;

        /* Update "same". */
        if (same[(pos - 1) & ZOPFLI_WINDOW_MASK] > 1)
        {
            amount = same[(pos - 1) & ZOPFLI_WINDOW_MASK] - 1;
        }
        while (pos + amount + 1 < end && array[pos] == array[pos + amount + 1] && amount < ushort.MaxValue) {
            amount++;
        }
        same[hpos] = (ushort)amount;

        val2 = ((same[hpos] - ZOPFLI_MIN_MATCH) & 255) ^ val;
        hashval2[hpos] = (short)val2;
        if (head2[val2] != -1 && hashval2[head2[val2]] == val2)
        {
            prev2[hpos] = (short)head2[val2];
            Debug.Assert(head[val2] <= short.MaxValue);

        }
        else prev2[hpos] = hpos;
        head2[val2] = hpos;
    }
    public void ZopfliWarmupHash(byte[] InFile, int pos, int end)
    {
        UpdateHashValue(InFile[pos + 0]);
        if (pos + 1 < end) UpdateHashValue(InFile[pos + 1]);
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