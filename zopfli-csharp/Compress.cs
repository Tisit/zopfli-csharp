using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using zopfli_csharp;

namespace ZopfliCSharp
{
    public partial class Compress
    {
        /*
        A block structure of huge, non-smart, blocks to divide the input into, to allow
        operating on huge files without exceeding memory, such as the 1GB wiki9 corpus.
        The whole compression algorithm, including the smarter block splitting, will
        be executed independently on each huge block.
        Dividing into huge blocks hurts compression, but not much relative to the size.
        Set it to 0 to disable master blocks.
        */
        const int ZOPFLI_MASTER_BLOCK_SIZE = 1000000;

        /* Minimum and maximum length that can be encoded in deflate. */
        const int ZOPFLI_MAX_MATCH = 258;
        const int ZOPFLI_MIN_MATCH = 3;

        /* Number of distinct literal/length and distance symbols in DEFLATE */
        public const int ZOPFLI_NUM_LL = 288;
        public const int ZOPFLI_NUM_D = 32;

        /*
        limit the max hash chain hits for this hash value. This has an effect only
        on files where the hash value is the same very often. On these files, this
        gives worse compression (the value should ideally be 32768, which is the
        ZOPFLI_WINDOW_SIZE, while zlib uses 4096 even for best level), but makes it
        faster on some specific files.
        Good value: e.g. 8192.
        */
        const int ZOPFLI_MAX_CHAIN_HITS = 8192;

        const int NUM = 9;  /* Good value: 9. */

        /*
        Used to initialize costs for example
        */
        const double ZOPFLI_LARGE_FLOAT = 1e30;

        public static byte[] ZopfliCompress(byte[] InFile)
        {

            if (Globals.output_type == ZopfliFormat.ZOPFLI_FORMAT_GZIP)
            {
                return ZopfliGzipCompress(InFile);
            }
            else if (Globals.output_type == ZopfliFormat.ZOPFLI_FORMAT_ZLIB)
            {
                return ZopfliZlibCompress(InFile);
            }
            else //(Globals.output_type == ZopfliFormat.ZOPFLI_FORMAT_DEFLATE)
            {
                return ZopfliDeflate(InFile, 2 /* Dynamic block */, 1).ToArray();
            }
        }
        /* CRC polynomial: 0xedb88320 */
        static readonly ulong [] crc32_table = {
           0, 1996959894, 3993919788, 2567524794,  124634137, 1886057615,
          3915621685, 2657392035,  249268274, 2044508324, 3772115230, 2547177864,
           162941995, 2125561021, 3887607047, 2428444049,  498536548, 1789927666,
          4089016648, 2227061214,  450548861, 1843258603, 4107580753, 2211677639,
           325883990, 1684777152, 4251122042, 2321926636,  335633487, 1661365465,
          4195302755, 2366115317,  997073096, 1281953886, 3579855332, 2724688242,
          1006888145, 1258607687, 3524101629, 2768942443,  901097722, 1119000684,
          3686517206, 2898065728,  853044451, 1172266101, 3705015759, 2882616665,
           651767980, 1373503546, 3369554304, 3218104598,  565507253, 1454621731,
          3485111705, 3099436303,  671266974, 1594198024, 3322730930, 2970347812,
           795835527, 1483230225, 3244367275, 3060149565, 1994146192,   31158534,
          2563907772, 4023717930, 1907459465,  112637215, 2680153253, 3904427059,
          2013776290,  251722036, 2517215374, 3775830040, 2137656763,  141376813,
          2439277719, 3865271297, 1802195444,  476864866, 2238001368, 4066508878,
          1812370925,  453092731, 2181625025, 4111451223, 1706088902,  314042704,
          2344532202, 4240017532, 1658658271,  366619977, 2362670323, 4224994405,
          1303535960,  984961486, 2747007092, 3569037538, 1256170817, 1037604311,
          2765210733, 3554079995, 1131014506,  879679996, 2909243462, 3663771856,
          1141124467,  855842277, 2852801631, 3708648649, 1342533948,  654459306,
          3188396048, 3373015174, 1466479909,  544179635, 3110523913, 3462522015,
          1591671054,  702138776, 2966460450, 3352799412, 1504918807,  783551873,
          3082640443, 3233442989, 3988292384, 2596254646,   62317068, 1957810842,
          3939845945, 2647816111,   81470997, 1943803523, 3814918930, 2489596804,
           225274430, 2053790376, 3826175755, 2466906013,  167816743, 2097651377,
          4027552580, 2265490386,  503444072, 1762050814, 4150417245, 2154129355,
           426522225, 1852507879, 4275313526, 2312317920,  282753626, 1742555852,
          4189708143, 2394877945,  397917763, 1622183637, 3604390888, 2714866558,
           953729732, 1340076626, 3518719985, 2797360999, 1068828381, 1219638859,
          3624741850, 2936675148,  906185462, 1090812512, 3747672003, 2825379669,
           829329135, 1181335161, 3412177804, 3160834842,  628085408, 1382605366,
          3423369109, 3138078467,  570562233, 1426400815, 3317316542, 2998733608,
           733239954, 1555261956, 3268935591, 3050360625,  752459403, 1541320221,
          2607071920, 3965973030, 1969922972,   40735498, 2617837225, 3943577151,
          1913087877,   83908371, 2512341634, 3803740692, 2075208622,  213261112,
          2463272603, 3855990285, 2094854071,  198958881, 2262029012, 4057260610,
          1759359992,  534414190, 2176718541, 4139329115, 1873836001,  414664567,
          2282248934, 4279200368, 1711684554,  285281116, 2405801727, 4167216745,
          1634467795,  376229701, 2685067896, 3608007406, 1308918612,  956543938,
          2808555105, 3495958263, 1231636301, 1047427035, 2932959818, 3654703836,
          1088359270,  936918000, 2847714899, 3736837829, 1202900863,  817233897,
          3183342108, 3401237130, 1404277552,  615818150, 3134207493, 3453421203,
          1423857449,  601450431, 3009837614, 3294710456, 1567103746,  711928724,
          3020668471, 3272380065, 1510334235,  755167117
        };

        /* Returns the CRC32 */
        static ulong CRC(byte[] data, int size) 
        {
            ulong result = 0xffffffff;
            for (int i = 0; size > 0; size--, i++) 
            {
                result = crc32_table[(result ^ data[i]) & 0xff] ^ (result >> 8);
            }
            return result ^ 0xffffffff;
        }

        /* Compresses the data according to the gzip specification, RFC 1952. */
        static byte[] ZopfliGzipCompress(byte[] InFile)
        {
            ulong crcvalue = CRC(InFile, InFile.Length);
            List<byte> OutFile = new List<byte>();

            OutFile.Add(31);  /* ID1 */
            OutFile.Add(139);  /* ID2 */
            OutFile.Add(8);  /* CM */
            OutFile.Add(0);  /* FLG */
            /* MTIME */
            OutFile.Add(0);
            OutFile.Add(0);
            OutFile.Add(0);
            OutFile.Add(0);

            OutFile.Add(2);  /* XFL, 2 indicates best compression. */
            OutFile.Add(3);  /* OS follows Unix conventions. */

            OutFile.AddRange(ZopfliDeflate(InFile, 2 /* Dynamic block */, 1));

            /* CRC */
            OutFile.Add((byte)(crcvalue % 256));
            OutFile.Add((byte)((crcvalue >> 8) % 256));
            OutFile.Add((byte)((crcvalue >> 16) % 256));
            OutFile.Add((byte)((crcvalue >> 24) % 256));

            /* ISIZE */
            OutFile.Add((byte)(InFile.Length % 256));
            OutFile.Add((byte)((InFile.Length >> 8) % 256));
            OutFile.Add((byte)((InFile.Length >> 16) % 256));
            OutFile.Add((byte)((InFile.Length >> 24) % 256));

            if (Globals.verbose == 1)
            {
               Console.WriteLine("Original Size: "+ InFile.Length + ", Gzip: " + OutFile.Count + 
                   " , Compression: " + 100.0 * (double)(InFile.Length - OutFile.Count) / (double)InFile.Length + " Removed");
            }
            return OutFile.ToArray();
        }

        /* Calculates the adler32 checksum of the data */
        static uint adler32(byte[] InFile)
        {
            const int mod = 65521;
            uint a = 1, b = 0;
            foreach (byte c in InFile)
            {
                a = (a + c) % mod;
                b = (b + a) % mod;
            }
            return (b << 16) | a;
        }

        static byte[] ZopfliZlibCompress(byte[] InFile)
        {
            uint checksum = adler32(InFile);
            uint cmf = 120;  /* CM 8, CINFO 7. See zlib spec.*/
            uint flevel = 3;
            uint fdict = 0;
            uint cmfflg = 256 * cmf + fdict * 32 + flevel * 64;
            uint fcheck = 31 - cmfflg % 31;
            List<byte> OutFile = new List<byte>();


            cmfflg += fcheck;

            OutFile.Add((byte)(cmfflg / 256));
            OutFile.Add((byte)(cmfflg % 256));

            OutFile.AddRange(ZopfliDeflate(InFile, 2 /* dynamic block */, 1 /* final */));

            OutFile.Add((byte)((checksum >> 24) % 256));
            OutFile.Add((byte)((checksum >> 16) % 256));
            OutFile.Add((byte)((checksum >> 8) % 256));
            OutFile.Add((byte)(checksum % 256));

            if (Globals.verbose > 0)
            {
                Console.WriteLine(
                        "Original Size: " + InFile.Length + ", Zlib: " + OutFile.Count + ", Compression: "
                        + 100.0 * (double)(InFile.Length - OutFile.Count) / (double)InFile.Length + "% Removed");
            }
            return OutFile.ToArray();
        }

        static List<byte> ZopfliDeflate(byte[] InFile, int btype, int final)
        {
            int i = 0;
            List<byte> OutFile = new List<byte>();

            byte[] bp = new byte[1];
            bp[0] = 0;

            do
            {
                bool masterfinal = (i + ZOPFLI_MASTER_BLOCK_SIZE >= InFile.Length);
                bool final2 = (final > 0) && masterfinal;
                int size = masterfinal ? InFile.Length - i : ZOPFLI_MASTER_BLOCK_SIZE;
                ZopfliDeflatePart(btype, final2, InFile, i, i + size, bp, OutFile);
                i += size;
            } while (i < InFile.Length);

            return OutFile;
        }

        /*
        Finds how long the match of scan and match is. Can be used to find how many
        bytes starting from scan, and from match, are equal. Returns the last byte
        after scan, which is still equal to the correspondinb byte after match.
        scan is the position to compare
        match is the earlier position to compare.
        end is the last possible byte, beyond which to stop looking.
        safe_end is a few (8) bytes before end, for comparing multiple bytes at once.
        */
        static int GetMatch(int scan, int match, int end, byte[] array)
        {

            /* just do it the naive way */
            while (scan < end && array[scan] == array[match])
            {
                scan++; match++;
            }

            return scan;
        }

        /*
        Gets distance, length and sublen values from the cache if possible.
        Returns 1 if it got the values from the cache, 0 if not.
        Updates the limit value to a smaller one if possible with more limited
        information from the cache.
        */
        static int TryGetFromLongestMatchCache(ZopfliBlockState s,
            ulong pos, ref int limit,
            ushort[] sublen, ref ushort distance, ref ushort length)
        {
            ulong lmcpos = pos - (ulong)s.blockstart;

            /* Length > 0 and dist 0 is invalid combination, which indicates on purpose
               that this cache value is not filled in yet. */
            bool cache_available = s.lmc.length != null  && (s.lmc.length[lmcpos] == 0 ||
                s.lmc.dist[lmcpos] != 0);
            bool limit_ok_for_cache = cache_available && 
                (limit == ZOPFLI_MAX_MATCH || s.lmc.length[lmcpos] <= limit ||
                (sublen != null && s.ZopfliMaxCachedSublen(lmcpos) >= limit));

            if (limit_ok_for_cache && s.lmc != null)
            {
                if (sublen == null || s.lmc.length[lmcpos]
                    <= s.ZopfliMaxCachedSublen(lmcpos))
                {
                    length = s.lmc.length[lmcpos];
                    if (length > limit) length = (ushort)limit;
                    if (sublen != null)
                    {
                        s.ZopfliCacheToSublen(lmcpos, length, sublen);
                        distance = sublen[length];
                        if (limit == ZOPFLI_MAX_MATCH && length >= ZOPFLI_MIN_MATCH)
                        {
                            Debug.Assert(sublen[length] == s.lmc.dist[lmcpos]);
                        }
                    }
                    else
                    {
                        distance = s.lmc.dist[lmcpos];
                    }
                    return 1;
                }
                /* Can't use much of the cache, since the "sublens" need to be calculated,
                   but at  least we already know when to stop. */
                limit = s.lmc.length[lmcpos];
            }

            return 0;
        }

        static void StoreInLongestMatchCache(ZopfliBlockState s, int pos, int limit, ushort[] sublen, ushort distance, ushort length)
        {
            /* The LMC cache starts at the beginning of the block rather than the
             beginning of the whole array. */
            int lmcpos = pos - s.blockstart;

            /* Length > 0 and dist 0 is invalid combination, which indicates on purpose
            that this cache value is not filled in yet. */
            bool cache_available = (s.lmc.length != null &&(s.lmc.length[lmcpos] == 0 ||
                s.lmc.dist[lmcpos] != 0));

            if (s.lmc.length != null && sublen != null && limit == ZOPFLI_MAX_MATCH && sublen.Length > 0 && !cache_available)
            {
                Debug.Assert(s.lmc.length[lmcpos] == 1 && s.lmc.dist[lmcpos] == 0);
                s.lmc.dist[lmcpos] = (ushort)(length < ZOPFLI_MIN_MATCH ? 0 : distance);
                s.lmc.length[lmcpos] = (ushort)(length < ZOPFLI_MIN_MATCH ? 0 : length);
                Debug.Assert(!(s.lmc.length[lmcpos] == 1 && s.lmc.dist[lmcpos] == 0));
                s.ZopfliSublenToCache(sublen, lmcpos, length);
            }
        }

        static void ZopfliFindLongestMatch(ZopfliBlockState s, ZopfliHash h, byte[] array, int pos, int size, int limit, ushort[] sublen,
                          ref ushort distance, ref ushort length)
        {
            ushort hpos = (ushort)(pos & ZopfliHash.ZOPFLI_WINDOW_MASK), p, pp;
            ushort bestdist = 0;
            ushort bestlength = 1;
            int scan;
            int match;
            
            int chain_counter = ZOPFLI_MAX_CHAIN_HITS;  /* For quitting early. */

            int[] hhead = h.head;
            short[] hprev = h.prev;
            short[] hhashval = h.hashval;
            int hval = h.val;

            if (TryGetFromLongestMatchCache(s, (ulong)pos, ref limit, sublen, ref distance, ref length) > 0) {
                Debug.Assert(pos + length <= size);
                return;
            }

            Debug.Assert(limit <= ZOPFLI_MAX_MATCH);
            Debug.Assert(limit >= ZOPFLI_MIN_MATCH);
            Debug.Assert(pos < size);

            if (size - pos < ZOPFLI_MIN_MATCH)
            {
                /* The rest of the code assumes there are at least ZOPFLI_MIN_MATCH bytes to
                   try. */
                length = 0;
                distance = 0;
                return;
            }

            if (pos + limit > size)
            {
                limit = size - pos;
            }

            int arrayend = pos + limit;

            Debug.Assert(hval < 65536);

            pp = (ushort)hhead[hval];  /* During the whole loop, p == hprev[pp]. */
            p = (ushort)hprev[pp];

            Debug.Assert(pp == hpos);

            ushort dist = (ushort)(p < pp ? pp - p : ZopfliHash.ZOPFLI_WINDOW_SIZE - p + pp);

            /* Go through all distances. */
            while (dist < ZopfliHash.ZOPFLI_WINDOW_SIZE)
            {
                ushort currentlength = 0;

                Debug.Assert(p < ZopfliHash.ZOPFLI_WINDOW_SIZE);
                Debug.Assert(p == hprev[pp]);
                Debug.Assert(hhashval[p] == hval);

                if (dist > 0)
                {
                    Debug.Assert(pos < size);
                    Debug.Assert(dist <= pos);
                    scan = pos;
                    match = pos - dist;

                    /* Testing the byte at position bestlength first, goes slightly faster. */
                    if (pos + bestlength >= size
                        || array[scan + bestlength] == array[match + bestlength])
                    {

                        ushort same0 = h.same[hpos];
                        if (same0 > 2 && array[scan] == array[match])
                        {
                            ushort same1 = h.same[(pos - dist) & ZopfliHash.ZOPFLI_WINDOW_MASK];
                            ushort same = same0 < same1 ? same0 : same1;
                            if (same > limit) same = (ushort)limit;
                            scan += same;
                            match += same;
                        }
                        scan = GetMatch(scan, match, arrayend, array);
                        currentlength = (ushort)(scan - pos);  /* The found length. */
                    }

                    if (currentlength > bestlength)
                    {
                        if (sublen != null && sublen.Length > 0)
                        {
                            ushort j;
                            for (j = (ushort)(bestlength + 1); j <= currentlength; j++)
                            {
                                sublen[j] = dist;
                            }
                        }
                        bestdist = dist;
                        bestlength = currentlength;
                        if (currentlength >= limit) break;
                    }
                }


                /* Switch to the other hash once this will be more efficient. */
                if (hhead != h.head2 && bestlength >= h.same[hpos] &&
                    h.val2 == h.hashval2[p])
                {
                    /* Now use the hash that encodes the length and first byte. */
                    hhead = h.head2;
                    hprev = h.prev2;
                    hhashval = h.hashval2;
                    hval = h.val2;
                }

                pp = p;
                p = (ushort)hprev[p];
                if (p == pp) break;  /* Uninited prev value. */

                dist += p < pp ? (ushort)(pp - p) : (ushort)((ZopfliHash.ZOPFLI_WINDOW_SIZE - p) + pp);

                chain_counter--;
                if (chain_counter <= 0) break;
            }

            StoreInLongestMatchCache(s, pos, limit, sublen, bestdist, bestlength);

            Debug.Assert(bestlength <= limit);

            distance = bestdist;
            length = bestlength;
            Debug.Assert(pos + length <= size);
        }

        /*
        Gets a score of the length given the distance. Typically, the score of the
        length is the length itself, but if the distance is very long, decrease the
        score of the length a bit to make up for the fact that long distances use large
        amounts of extra bits.

        This is not an accurate score, it is a heuristic only for the greedy LZ77
        implementation. More accurate cost models are employed later. Making this
        heuristic more accurate may hurt rather than improve compression.

        The two direct uses of this heuristic are:
        -avoid using a length of 3 in combination with a long distance. This only has
         an effect if length == 3.
        -make a slightly better choice between the two options of the lazy matching.

        Indirectly, this affects:
        -the block split points if the default of block splitting first is used, in a
         rather unpredictable way
        -the first zopfli run, so it affects the chance of the first run being closer
         to the optimal output
        */
        static int GetLengthScore(int length, int distance)
        {
            /*
            At 1024, the distance uses 9+ extra bits and this seems to be the sweet spot
            on tested files.
            */
            return distance > 1024 ? length - 1 : length;
        }

        static public void ZopfliVerifyLenDist(byte[] data, int datasize, int pos,
                         ushort dist, ushort length) 
        {

            /* TODO(lode): make this only run in a debug compile, it's for assert only. */
            int i;

            Debug.Assert(pos + length <= datasize);
            for (i = 0; i<length; i++) 
            {
                if (data[pos - dist + i] != data[pos + i]) 
                {
                    Debug.Assert(data[pos - dist + i] == data[pos + i]);
                    break;
                }
            }
        }

        /*
        Appends the length and distance to the LZ77 arrays of the ZopfliLZ77Store.
        context must be a ZopfliLZ77Store*.
        */
        static void ZopfliStoreLitLenDist(ushort length, ushort dist, int pos, ZopfliLZ77Store store)
        {
            int i;
            /* Needed for using ZOPFLI_APPEND_DATA multiple times. */
            int origsize = (int)store.size;
            int llstart = ZOPFLI_NUM_LL * (origsize / ZOPFLI_NUM_LL);
            int dstart = ZOPFLI_NUM_D * (origsize / ZOPFLI_NUM_D);

            /* Everytime the index wraps around, a new cumulative histogram is made: we're
            keeping one histogram value per LZ77 symbol rather than a full histogram for
            each to save memory. */
            if (origsize % ZOPFLI_NUM_LL == 0)
            {
                for (i = 0; i < ZOPFLI_NUM_LL; i++)
                {
                    store.ll_counts.Add(origsize == 0 ? 0 : store.ll_counts[origsize - ZOPFLI_NUM_LL + i]);
                }
            }

            if (origsize % ZOPFLI_NUM_D == 0)
            {
                for (i = 0; i < ZOPFLI_NUM_D; i++)
                {
                    store.d_counts.Add(origsize == 0 ? 0 : store.d_counts[origsize - ZOPFLI_NUM_D + i]);
                }
            }

            store.litlens.Add(length);
            store.dists.Add(dist);
            store.pos.Add((ulong)pos);
            Debug.Assert(length < 259);
            store.size++;

            if (dist == 0)
            {
                store.ll_symbol.Add(length);
                store.d_symbol.Add(0);
                store.ll_counts[llstart + length]++;
            }
            else
            {
                store.ll_symbol.Add(Symbols.ZopfliGetLengthSymbol(length));
                store.d_symbol.Add(Symbols.ZopfliGetDistSymbol(dist));
                store.ll_counts[llstart + Symbols.ZopfliGetLengthSymbol(length)]++;
                store.d_counts[dstart + Symbols.ZopfliGetDistSymbol(dist)]++;
            }

        }

        //TODO: doesn't work entirely correctly. WHYYYY???
        static void ZopfliLZ77Greedy(ZopfliBlockState s, byte[] InFile, int instart, int inend, ZopfliLZ77Store store, ZopfliHash h)
        {
            int i, j;
            ushort leng = 0;
            ushort dist = 0;
            int lengthscore;
            int windowstart = instart > ZopfliHash.ZOPFLI_WINDOW_SIZE
                ? instart - ZopfliHash.ZOPFLI_WINDOW_SIZE : 0;
            ushort[] dummysublen = new ushort[259];

            /* Lazy matching. */
            uint prev_length = 0;
            uint prev_match = 0;
            int prevlengthscore;
            int match_available = 0;

            if (instart == inend) return;

            h.ZopfliResetHash();
            h.ZopfliWarmupHash(InFile, windowstart, inend);
            for (i = windowstart; i < instart; i++)
            {
                h.ZopfliUpdateHash(InFile, i, inend);
            }

            for (i = instart; i < inend; i++)
            {
                h.ZopfliUpdateHash(InFile, i, inend);
                ZopfliFindLongestMatch(s, h, InFile, i, inend, ZOPFLI_MAX_MATCH, dummysublen,
                           ref dist, ref leng);
                lengthscore = GetLengthScore(leng, dist);
                /* Lazy matching. */
                prevlengthscore = GetLengthScore((int)prev_length, (int)prev_match);
                if (match_available > 0)
                {
                    match_available = 0;
                    if (lengthscore > prevlengthscore + 1)
                    {
                        ZopfliStoreLitLenDist(InFile[i - 1], 0, i - 1, store);
                        if (lengthscore >= ZOPFLI_MIN_MATCH && leng < ZOPFLI_MAX_MATCH)
                        {
                            match_available = 1;
                            prev_length = leng;
                            prev_match = dist;
                            continue;
                        }
                    }
                    else
                    {
                        /* Add previous to output. */
                        leng = (ushort)prev_length;
                        dist = (ushort)prev_match;
                        lengthscore = prevlengthscore;
                        /* Add to output. */
                        ZopfliVerifyLenDist(InFile, inend, i - 1, dist, leng);
                        ZopfliStoreLitLenDist(leng, dist, i - 1, store);
                        for (j = 2; j < leng; j++)
                        {
                            Debug.Assert(i < inend);
                            i++;
                            h.ZopfliUpdateHash(InFile, i, inend);
                        }
                        continue;
                    }
                }
                else if (lengthscore >= ZOPFLI_MIN_MATCH && leng < ZOPFLI_MAX_MATCH)
                {
                    match_available = 1;
                    prev_length = leng;
                    prev_match = dist;
                    continue;
                }
                /* End of lazy matching. */

                /* Add to output. */
                if (lengthscore >= ZOPFLI_MIN_MATCH)
                {
                    ZopfliVerifyLenDist(InFile, inend, i, dist, leng);
                    ZopfliStoreLitLenDist(leng, dist, i, store);
                }
                else
                {
                    leng = 1;
                    ZopfliStoreLitLenDist(InFile[i], 0, i, store);
                }
                for (j = 1; j < leng; j++)
                {
                    Debug.Assert(i < inend);
                    i++;
                    h.ZopfliUpdateHash(InFile, i, inend);
                }
            }
        }

        static ulong ZopfliLZ77GetByteRange(ZopfliLZ77Store lz77,
                                      ulong lstart, ulong lend) {
            int l = (int)lend - 1;
            if (lstart == lend) return 0;
            return lz77.pos[l] + ((lz77.dists[l] == 0) ?
                (ulong)1 : lz77.litlens[l]) - lz77.pos[(int)lstart];
        }

        static void GetFixedTree(uint[] ll_lengths, uint[] d_lengths)
        {
            ulong i;
            for (i = 0; i < 144; i++) ll_lengths[i] = 8;
            for (i = 144; i < 256; i++) ll_lengths[i] = 9;
            for (i = 256; i < 280; i++) ll_lengths[i] = 7;
            for (i = 280; i < 288; i++) ll_lengths[i] = 8;
            for (i = 0; i < 32; i++) d_lengths[i] = 5;
        }

        static ulong CalculateBlockSymbolSizeSmall( uint[] ll_lengths,
                                                    uint[] d_lengths,
                                                    ZopfliLZ77Store lz77,
                                                    ulong lstart, ulong lend) 
        {
            ulong result = 0;
            ulong i;
            for (i = lstart; i < lend; i++) 
            {
                Debug.Assert(i < lz77.size);
                Debug.Assert(lz77.litlens[(int)i] < 259);
                if (lz77.dists[(int)i] == 0) 
                {
                    result += ll_lengths[lz77.litlens[(int)i]];
                } else 
                {
                    int ll_symbol = Symbols.ZopfliGetLengthSymbol(lz77.litlens[(int)i]);
                    int d_symbol = Symbols.ZopfliGetDistSymbol(lz77.dists[(int)i]);
                    result += ll_lengths[ll_symbol];
                    result += d_lengths[d_symbol];
                    result += (ulong)Symbols.ZopfliGetLengthSymbolExtraBits(ll_symbol);
                    result += (ulong)Symbols.ZopfliGetDistSymbolExtraBits(d_symbol);
                }
            }
            result += ll_lengths[256]; /*end symbol*/
            return result;
        }

        static void ZopfliLZ77GetHistogramAt(ZopfliLZ77Store lz77, ulong lpos,
                                     uint[] ll_counts, uint[] d_counts)
        {
            /* The real histogram is created by using the histogram for this chunk, but
            all superfluous values of this chunk subtracted. */
            ulong llpos = ZOPFLI_NUM_LL * (lpos / ZOPFLI_NUM_LL);
            ulong dpos = ZOPFLI_NUM_D * (lpos / ZOPFLI_NUM_D);
            ulong i;
            for (i = 0; i < ZOPFLI_NUM_LL; i++)
            {
                ll_counts[i] = (uint)lz77.ll_counts[(int)(llpos + i)];
            }
            for (i = lpos + 1; i < llpos + ZOPFLI_NUM_LL && i < lz77.size; i++)
            {
                ll_counts[lz77.ll_symbol[(int)i]]--;
            }
            for (i = 0; i < ZOPFLI_NUM_D; i++)
            {
                d_counts[i] = (uint)lz77.d_counts[(int)(dpos + i)];
            }
            for (i = lpos + 1; i < dpos + ZOPFLI_NUM_D && i < lz77.size; i++)
            {
                if (lz77.dists[(int)i] != 0) d_counts[lz77.d_symbol[(int)i]]--;
            }
        }

        static void ZopfliLZ77GetHistogram(ZopfliLZ77Store lz77,
                           ulong lstart, ulong lend,
                           uint[] ll_counts, uint[] d_counts)
        {
            ulong i;
            if (lstart + ZOPFLI_NUM_LL * 3 > lend)
            {
                for (i = lstart; i < lend; i++)
                {
                    ll_counts[lz77.ll_symbol[(int)i]]++;
                    if (lz77.dists[(int)i] != 0) d_counts[lz77.d_symbol[(int)i]]++;
                }
            }
            else
            {
                /* Subtract the cumulative histograms at the end and the start to get the
                histogram for this range. */
                ZopfliLZ77GetHistogramAt(lz77, lend - 1, ll_counts, d_counts);
                if (lstart > 0)
                {
                    uint[] ll_counts2 = new uint[ZOPFLI_NUM_LL];
                    uint[] d_counts2 = new uint[ZOPFLI_NUM_D];
                    ZopfliLZ77GetHistogramAt(lz77, lstart - 1, ll_counts2, d_counts2);

                    for (i = 0; i < ZOPFLI_NUM_LL; i++)
                    {
                        ll_counts[i] -= ll_counts2[i];
                    }
                    for (i = 0; i < ZOPFLI_NUM_D; i++)
                    {
                        d_counts[i] -= d_counts2[i];
                    }
                }
            }
        }

        /*
        Same as CalculateBlockSymbolSize, but with the histogram provided by the caller.
        */
        static ulong CalculateBlockSymbolSizeGivenCounts(uint[] ll_counts,
                                                  uint[] d_counts,
                                                  uint[] ll_lengths,
                                                  uint[] d_lengths,
                                                  ZopfliLZ77Store lz77,
                                                  ulong lstart, ulong lend)
        {
            ulong result = 0;
            int i;
            if (lstart + ZOPFLI_NUM_LL * 3 > lend)
            {
                return CalculateBlockSymbolSizeSmall(
                    ll_lengths, d_lengths, lz77, lstart, lend);
            }
            else
            {
                for (i = 0; i < 256; i++)
                {
                    result += ll_lengths[i] * ll_counts[i];
                }
                for (i = 257; i < 286; i++)
                {
                    result += ll_lengths[i] * ll_counts[i];
                    result += (ulong)Symbols.ZopfliGetLengthSymbolExtraBits(i) * ll_counts[i];
                }
                for (i = 0; i < 30; i++)
                {
                    result += d_lengths[i] * d_counts[i];
                    result += (ulong)Symbols.ZopfliGetDistSymbolExtraBits(i) * d_counts[i];
                }
                result += ll_lengths[256]; /*end symbol*/
                return result;
            }
        }

        /*
        Calculates size of the part after the header and tree of an LZ77 block, in bits.
        */
        static ulong CalculateBlockSymbolSize( uint[] ll_lengths,
                                       uint[] d_lengths,
                                       ZopfliLZ77Store lz77,
                                       ulong lstart, ulong lend) 
        {
            if (lstart + ZOPFLI_NUM_LL* 3 > lend) 
            {
                return CalculateBlockSymbolSizeSmall(
                    ll_lengths, d_lengths, lz77, lstart, lend);
            } else 
            {
                uint[] ll_counts = new uint[ZOPFLI_NUM_LL];
                uint[] d_counts = new uint[ZOPFLI_NUM_D];
                ZopfliLZ77GetHistogram(lz77, lstart, lend, ll_counts, d_counts);
                return CalculateBlockSymbolSizeGivenCounts(
                    ll_counts, d_counts, ll_lengths, d_lengths, lz77, lstart, lend);
            }
        }

        static void ZopfliCalculateBitLengths(uint[] count, ulong n, int maxbits,
                               uint[] bitlengths)
        {
            int error = Katajainen.ZopfliLengthLimitedCodeLengths(count, (int)n, maxbits, bitlengths);
            Debug.Assert(error == 0);
        }

        /*
        Ensures there are at least 2 distance codes to support buggy decoders.
        Zlib 1.2.1 and below have a bug where it fails if there isn't at least 1
        distance code (with length > 0), even though it's valid according to the
        deflate spec to have 0 distance codes. On top of that, some mobile phones
        require at least two distance codes. To support these decoders too (but
        potentially at the cost of a few bytes), add dummy code lengths of 1.
        References to this bug can be found in the changelog of
        Zlib 1.2.2 and here: http://www.jonof.id.au/forum/index.php?topic=515.0.

        d_lengths: the 32 lengths of the distance codes.
        */
        static void PatchDistanceCodesForBuggyDecoders(uint[] d_lengths)
        {
            int num_dist_codes = 0; /* Amount of non-zero distance codes */
            int i;
            for (i = 0; i < 30 /* Ignore the two unused codes from the spec */; i++)
            {
                if (d_lengths[i] > 0) num_dist_codes++;
                if (num_dist_codes >= 2) return; /* Two or more codes is fine. */
            }

            if (num_dist_codes == 0)
            {
                d_lengths[0] = d_lengths[1] = 1;
            }
            else if (num_dist_codes == 1)
            {
                d_lengths[d_lengths[0] > 0 ? 1 : 0] = 1;
            }
        }

        /*
        Changes the population counts in a way that the consequent Huffman tree
        compression, especially its rle-part, will be more likely to compress this data
        more efficiently. length contains the size of the histogram.
        */
        static void OptimizeHuffmanForRle(int length, uint[] counts)
        {
            int i, k, stride;
            ulong symbol, sum, limit;

            /* 1) We don't want to touch the trailing zeros. We may break the
            rules of the format by adding more data in the distance codes. */
            for (; length >= 0; --length)
            {
                if (length == 0)
                {
                    return;
                }
                if (counts[length - 1] != 0)
                {
                    /* Now counts[0..length - 1] does not have trailing zeros. */
                    break;
                }
            }
            /* 2) Let's mark all population counts that already can be encoded
            with an rle code.*/
            int[] good_for_rle = new int[length];
            for (i = 0; i < length; ++i) good_for_rle[i] = 0;

            /* Let's not spoil any of the existing good rle codes.
            Mark any seq of 0's that is longer than 5 as a good_for_rle.
            Mark any seq of non-0's that is longer than 7 as a good_for_rle.*/
            symbol = counts[0];
            stride = 0;
            for (i = 0; i < length + 1; ++i)
            {
                if (i == length || counts[i] != symbol)
                {
                    if ((symbol == 0 && stride >= 5) || (symbol != 0 && stride >= 7))
                    {
                        for (k = 0; k < stride; ++k)
                        {
                            good_for_rle[i - k - 1] = 1;
                        }
                    }
                    stride = 1;
                    if (i != length)
                    {
                        symbol = counts[i];
                    }
                }
                else
                {
                    ++stride;
                }
            }

            /* 3) Let's replace those population counts that lead to more rle codes. */
            stride = 0;
            limit = counts[0];
            sum = 0;
            for (i = 0; i < length + 1; ++i)
            {
                if (i == length || good_for_rle[i] > 0
                    /* Heuristic for selecting the stride ranges to collapse. */
                    || Math.Abs((int)(counts[i] - limit)) >= 4)
                {
                    if (stride >= 4 || (stride >= 3 && sum == 0))
                    {
                        /* The stride must end, collapse what we have, if we have enough (4). */
                        int count = ((int)sum + stride / 2) / stride;
                        if (count < 1) count = 1;
                        if (sum == 0)
                        {
                            /* Don't make an all zeros stride to be upgraded to ones. */
                            count = 0;
                        }
                        for (k = 0; k < stride; ++k)
                        {
                            /* We don't want to change value at counts[i],
                            that is already belonging to the next stride. Thus - 1. */
                            counts[i - k - 1] = (uint)count;
                        }
                    }
                    stride = 0;
                    sum = 0;
                    if (i < length - 3)
                    {
                        /* All interesting strides have a count of at least 4,
                        at least when non-zeros. */
                        limit = (counts[i] + counts[i + 1] +
                                 counts[i + 2] + counts[i + 3] + 2) / 4;
                    }
                    else if (i < length)
                    {
                        limit = counts[i];
                    }
                    else
                    {
                        limit = 0;
                    }
                }
                ++stride;
                if (i != length)
                {
                    sum += counts[i];
                }
            }

        }

        /*
        Tries out OptimizeHuffmanForRle for this block, if the result is smaller,
        uses it, otherwise keeps the original. Returns size of encoded tree and data in
        bits, not including the 3-bit block header.
        */
        static double TryOptimizeHuffmanForRle(
            ZopfliLZ77Store lz77, ulong lstart, ulong lend,
            uint[] ll_counts, uint[] d_counts,
            uint[] ll_lengths, uint[] d_lengths)
        {
            uint[] ll_counts2 = new uint[ZOPFLI_NUM_LL];
            uint[] d_counts2 = new uint[ZOPFLI_NUM_D];
            uint[] ll_lengths2 = new uint[ZOPFLI_NUM_LL];
            uint[] d_lengths2 = new uint[ZOPFLI_NUM_D];
            double treesize;
            double datasize;
            double treesize2;
            double datasize2;

            treesize = CalculateTreeSize(ll_lengths, d_lengths);
            datasize = CalculateBlockSymbolSizeGivenCounts(ll_counts, d_counts,
                ll_lengths, d_lengths, lz77, lstart, lend);

            ll_counts.CopyTo(ll_counts2, 0);
            d_counts.CopyTo(d_counts2, 0);
            OptimizeHuffmanForRle(ZOPFLI_NUM_LL, ll_counts2);
            OptimizeHuffmanForRle(ZOPFLI_NUM_D, d_counts2);
            ZopfliCalculateBitLengths(ll_counts2, ZOPFLI_NUM_LL, 15, ll_lengths2);
            ZopfliCalculateBitLengths(d_counts2, ZOPFLI_NUM_D, 15, d_lengths2);
            PatchDistanceCodesForBuggyDecoders(d_lengths2);

            treesize2 = CalculateTreeSize(ll_lengths2, d_lengths2);
            datasize2 = CalculateBlockSymbolSizeGivenCounts(ll_counts, d_counts,
                ll_lengths2, d_lengths2, lz77, lstart, lend);

            if (treesize2 + datasize2 < treesize + datasize)
            {
                ll_lengths2.CopyTo(ll_lengths, 0);
                d_lengths2.CopyTo(d_lengths, 0);
                return treesize2 + datasize2;
            }
            return treesize + datasize;
        }

        /*
        Calculates the bit lengths for the symbols for dynamic blocks. Chooses bit
        lengths that give the smallest size of tree encoding + encoding of all the
        symbols to have smallest output size. This are not necessarily the ideal Huffman
        bit lengths. Returns size of encoded tree and data in bits, not including the
        3-bit block header.
        */
        static double GetDynamicLengths(ZopfliLZ77Store lz77,
                                        ulong lstart, ulong lend,
                                        uint[] ll_lengths, uint[] d_lengths)
        {
            uint[] ll_counts = new uint[ZOPFLI_NUM_LL];
            uint[] d_counts = new uint[ZOPFLI_NUM_D];

            ZopfliLZ77GetHistogram(lz77, lstart, lend, ll_counts, d_counts);
            ll_counts[256] = 1;  /* End symbol. */
            ZopfliCalculateBitLengths(ll_counts, ZOPFLI_NUM_LL, 15, ll_lengths);
            ZopfliCalculateBitLengths(d_counts, ZOPFLI_NUM_D, 15, d_lengths);
            PatchDistanceCodesForBuggyDecoders(d_lengths);
            return TryOptimizeHuffmanForRle(
                lz77, lstart, lend, ll_counts, d_counts, ll_lengths, d_lengths);
        }

        static double ZopfliCalculateBlockSize(ZopfliLZ77Store lz77,
                                        ulong lstart, ulong lend, int btype) 
        {
            uint[] ll_lengths = new uint[ZOPFLI_NUM_LL];
            uint[] d_lengths = new uint[ZOPFLI_NUM_D];

            double result = 3; /* bfinal and btype bits */

            if (btype == 0) 
            {
                ulong length = ZopfliLZ77GetByteRange(lz77, lstart, lend);
                ulong rem = length % 65535;
                ulong blocks = length / 65535 + (ulong)(rem > 0 ? 1 : 0);
                /* An uncompressed block must actually be split into multiple blocks if it's
                   larger than 65535 bytes long. Eeach block header is 5 bytes: 3 bits,
                   padding, LEN and NLEN (potential less padding for first one ignored). */
                return blocks* 5 * 8 + length* 8;
            } 
            
            if (btype == 1) 
            {
                GetFixedTree(ll_lengths, d_lengths);
                result += CalculateBlockSymbolSize(
                    ll_lengths, d_lengths, lz77, lstart, lend);
            } else 
            {
                result += GetDynamicLengths(lz77, lstart, lend, ll_lengths, d_lengths);
            }

            return result;
        }

        static double ZopfliCalculateBlockSizeAutoType(ZopfliLZ77Store lz77,
                                        ulong lstart, ulong lend) 
        {
            double uncompressedcost = ZopfliCalculateBlockSize(lz77, lstart, lend, 0);
            /* Don't do the expensive fixed cost calculation for larger blocks that are
                unlikely to use it. */
            double fixedcost = (lz77.size > 1000) ?
                uncompressedcost : ZopfliCalculateBlockSize(lz77, lstart, lend, 1);
            double dyncost = ZopfliCalculateBlockSize(lz77, lstart, lend, 2);
            return (uncompressedcost<fixedcost && uncompressedcost<dyncost)
              ? uncompressedcost
              : (fixedcost < dyncost ? fixedcost : dyncost);
        }

        /*
        Returns estimated cost of a block in bits.  It includes the size to encode the
        tree and the size to encode all literal, length and distance symbols and their
        extra bits.

        litlens: lz77 lit/lengths
        dists: ll77 distances
        lstart: start of block
        lend: end of block (not inclusive)
        */
        static double EstimateCost(ZopfliLZ77Store lz77, ulong lstart, ulong lend) 
         {
            return ZopfliCalculateBlockSizeAutoType(lz77, lstart, lend);
        }

        /*
        Gets the cost which is the sum of the cost of the left and the right section
        of the data.
        type: FindMinimumFun
        */
        static double SplitCost(ulong i, SplitCostContext context)
        {
            return EstimateCost(context.lz77, context.start, i) + EstimateCost(context.lz77, i, context.end);
        }

        static void AddSorted(ulong value, List<ulong> out1, ref ulong npoints)
        {
            ulong i;
            out1.Add(value);
            npoints++;
            for (i = 0; i + 1 < npoints; i++)
            {
                if (out1[(int)i] > value)
                {
                    ulong j;
                    for (j = npoints - 1; j > i; j--)
                    {
                        out1[(int)j] = out1[(int)j - 1];
                    }
                    out1[(int)i] = value;
                    break;
                }
            }
        }

        static ulong FindMinimum(SplitCostContext context, ulong start, ulong end, out double smallest)
        {
            if (end - start < 1024)
            {
                double best = ZOPFLI_LARGE_FLOAT;
                ulong result = start;
                ulong i;
                for (i = start; i < end; i++)
                {
                    double v = SplitCost(i, context);
                    if (v < best)
                    {
                        best = v;
                        result = i;
                    }
                }
                smallest = best;
                return result;
            }
            else
            {
                /* Try to find minimum faster by recursively checking multiple points. */
                ulong i;
                ulong[] p = new ulong[NUM];
                double[] vp = new double[NUM];
                ulong besti;
                double best;
                double lastbest = ZOPFLI_LARGE_FLOAT;
                ulong pos = start;

                for (; ; )
                {
                    if (end - start <= NUM) break;

                    for (i = 0; i < NUM; i++)
                    {
                        p[i] = start + (i + 1) * ((end - start) / (NUM + 1));
                        vp[i] = SplitCost(p[i], context);
                    }
                    besti = 0;
                    best = vp[0];
                    for (i = 1; i < NUM; i++)
                    {
                        if (vp[i] < best)
                        {
                            best = vp[i];
                            besti = i;
                        }
                    }
                    if (best > lastbest) break;

                    start = besti == 0 ? start : p[besti - 1];
                    end = besti == NUM - 1 ? end : p[besti + 1];

                    pos = p[besti];
                    lastbest = best;
                }
                smallest = lastbest;
                return pos;
            }
        }

        /*
        Prints the block split points as decimal and hex values in the terminal.
        */
        static void PrintBlockSplitPoints(ZopfliLZ77Store lz77,
                                  List<ulong> lz77splitpoints,
                                  ulong nlz77points)
        {
            List<ulong> splitpoints  = new List<ulong>();
            ulong npoints = 0;
            ulong i;
            /* The input is given as lz77 indices, but we want to see the uncompressed
            index values. */
            ulong pos = 0;
            if (nlz77points > 0)
            {
                for (i = 0; i < lz77.size; i++)
                {
                    ulong length = lz77.dists[(int)i] == 0 ? (ulong)1 : lz77.litlens[(int)i];
                    if (lz77splitpoints[(int)npoints] == i)
                    {
                        splitpoints.Add(pos);
                        if (npoints == nlz77points) break;
                    }
                    pos += length;
                }
            }
            Debug.Assert(npoints == nlz77points);

            Console.Write("block split points: ");
            for (i = 0; i < npoints; i++)
            {
                Console.Write("%d ", (int)splitpoints[(int)i]);
            }
            Console.Write("(hex:");
            for (i = 0; i < npoints; i++)
            {
                Console.Write(" %x", (int)splitpoints[(int)i]);
            }
            Console.Write(")\n");
        }

        /*
        Finds next block to try to split, the largest of the available ones.
        The largest is chosen to make sure that if only a limited amount of blocks is
        requested, their sizes are spread evenly.
        lz77size: the size of the LL77 data, which is the size of the done array here.
        done: array indicating which blocks starting at that position are no longer
            splittable (splitting them increases rather than decreases cost).
        splitpoints: the splitpoints found so far.
        npoints: the amount of splitpoints found so far.
        lstart: output variable, giving start of block.
        lend: output variable, giving end of block.
        returns 1 if a block was found, 0 if no block found (all are done).
        */
        static int FindLargestSplittableBlock(
            ulong lz77size, byte[] done,
            List<ulong> splitpoints, ulong npoints,
            ref ulong lstart, ref ulong lend)
        {
            ulong longest = 0;
            int found = 0;
            ulong i;
            for (i = 0; i <= npoints; i++)
            {
                ulong start = i == 0 ? 0 : splitpoints[(int)i - 1];
                ulong end = i == npoints ? lz77size - 1 : splitpoints[(int)i];
                if (done[start] == 0 && end - start > longest)
                {
                    lstart = start;
                    lend = end;
                    found = 1;
                    longest = end - start;
                }
            }
            return found;
        }

        static void ZopfliAppendLZ77Store(ZopfliLZ77Store store,
                           ZopfliLZ77Store target)
        {
            int i;
            for (i = 0; i < (int)store.size; i++)
            {
                ZopfliStoreLitLenDist(store.litlens[i], store.dists[i],
                                      (int)store.pos[i], target);
            }
        }

        /*
        bp = bitpointer, always in range [0, 7].
        The outsize is number of necessary bytes to encode the bits.
        Given the value of bp and the amount of bytes, the amount of bits represented
        is not simply bytesize * 8 + bp because even representing one bit requires a
        whole byte. It is: (bp == 0) ? (bytesize * 8) : ((bytesize - 1) * 8 + bp)
        */
        static void AddBit(int bit,
                            byte[] bp, List<byte> OutFile)
        {
            if (bp[0] == 0) OutFile.Add(0);
            OutFile[OutFile.Count - 1] |= (byte)(bit << bp[0]);
            bp[0] = (byte)((bp[0] + 1) & 7);
        }


        static void AddBits(uint symbol, uint length,
                    byte[] bp, List<byte> OutFile)
        {
            /* TODO(lode): make more efficient (add more bits at once). */
            ushort i;
            for (i = 0; i < length; i++)
            {
                uint bit = (symbol >> i) & 1;
                if (bp[0] == 0) OutFile.Add(0);
                OutFile[OutFile.Count - 1] |= (byte)(bit << bp[0]);
                bp[0] = (byte)((bp[0] + 1) & 7);
            }
        }

        /*
        Adds bits, like AddBits, but the order is inverted. The deflate specification
        uses both orders in one standard.
        */
        static void AddHuffmanBits(uint symbol, uint length,
                                   byte[] bp, List<byte> OutFile)
        {
            /* TODO(lode): make more efficient (add more bits at once). */
            uint i;
            for (i = 0; i < length; i++)
            {
                uint bit = (symbol >> (ushort)(length - i - 1)) & 1;
                if (bp[0] == 0) OutFile.Add(0);
                OutFile[OutFile.Count - 1] |= (byte)(bit << bp[0]);
                bp[0] = (byte)((bp[0] + 1) & 7);
            }
        }

        static void ZopfliLengthsToSymbols(uint[] lengths, ulong n, uint maxbits,
                            uint[] symbols)
        {
            ulong[] bl_count = new ulong[maxbits + 1];
            uint[] next_code = new uint[maxbits + 1];
            uint bits, i;
            uint code;

            for (i = 0; i < n; i++)
            {
                symbols[i] = 0;
            }

            /* 1) Count the number of codes for each code length. Let bl_count[N] be the
            number of codes of length N, N >= 1. */
            for (bits = 0; bits <= maxbits; bits++)
            {
                bl_count[bits] = 0;
            }
            for (i = 0; i < n; i++)
            {
                Debug.Assert(lengths[i] <= maxbits);
                bl_count[lengths[i]]++;
            }
            /* 2) Find the numerical value of the smallest code for each code length. */
            code = 0;
            bl_count[0] = 0;
            for (bits = 1; bits <= maxbits; bits++)
            {
                code = (uint)(code + bl_count[bits - 1]) << 1;
                next_code[bits] = code;
            }
            /* 3) Assign numerical values to all codes, using consecutive values for all
            codes of the same length with the base values determined at step 2. */
            for (i = 0; i < n; i++)
            {
                uint len = lengths[i];
                if (len != 0)
                {
                    symbols[i] = next_code[len];
                    next_code[len]++;
                }
            }

        }
        static readonly uint[] EncodeTreeOrder = { 16, 17, 18, 0, 8, 7, 9, 6, 10, 5, 11, 4, 12, 3, 13, 2, 14, 1, 15 };
        /*
        Encodes the Huffman tree and returns how many bits its encoding takes. If out
        is a null pointer, only returns the size and runs faster.
        */
        static ulong EncodeTree(uint[] ll_lengths,
                         uint[] d_lengths,
                         bool use_16, bool use_17, bool use_18,
                         byte[] bp,
                         List<byte> OutFile, bool size_only)
        {
            uint lld_total;  /* Total amount of literal, length, distance codes. */
            /* Runlength encoded version of lengths of litlen and dist trees. */
            List<uint> rle = new List<uint>();
            List<uint> rle_bits = new List<uint>();  /* Extra bits for rle values 16, 17 and 18. */
            uint hlit = 29;  /* 286 - 257 */
            uint hdist = 29;  /* 32 - 1, but gzip does not like hdist > 29.*/
            uint hclen;
            uint hlit2;
            ulong i, j;
            uint[] clcounts= new uint[19];
            uint[] clcl = new uint[19];  /* Code length code lengths. */
            uint[] clsymbols = new uint[19];
            /* The order in which code length code lengths are encoded as per deflate. */
            
            ulong result_size = 0;

            /* Trim zeros. */
            while (hlit > 0 && ll_lengths[257 + hlit - 1] == 0) hlit--;
            while (hdist > 0 && d_lengths[1 + hdist - 1] == 0) hdist--;
            hlit2 = hlit + 257;

            lld_total = hlit2 + hdist + 1;

            for (i = 0; i < lld_total; i++)
            {
                /* This is an encoding of a huffman tree, so now the length is a symbol */
                byte symbol = i < hlit2 ? (byte) ll_lengths[i] : (byte) d_lengths[i - hlit2];
                uint count = 1;
                if (use_16 || (symbol == 0 && (use_17 || use_18)))
                {
                    for (j = i + 1; j < lld_total && symbol ==
                        (j < hlit2 ? ll_lengths[j] : d_lengths[j - hlit2]); j++)
                    {
                        count++;
                    }
                }
                i += count - 1;

                /* Repetitions of zeroes */
                if (symbol == 0 && count >= 3)
                {
                    if (use_18)
                    {
                        while (count >= 11)
                        {
                            uint count2 = count > 138 ? 138 : count;
                            if (!size_only)
                            {
                                rle.Add(18);
                                rle_bits.Add(count2 - 11);
                            }
                            clcounts[18]++;
                            count -= count2;
                        }
                    }
                    if (use_17)
                    {
                        while (count >= 3)
                        {
                            uint count2 = count > 10 ? 10 : count;
                            if (!size_only)
                            {
                                rle.Add(17);
                                rle_bits.Add(count2 - 3);
                            }
                            clcounts[17]++;
                            count -= count2;
                        }
                    }
                }

                /* Repetitions of any symbol */
                if (use_16 && count >= 4)
                {
                    count--;  /* Since the first one is hardcoded. */
                    clcounts[symbol]++;
                    if (!size_only)
                    {
                        rle.Add(symbol);
                        rle_bits.Add(0);
                    }
                    while (count >= 3)
                    {
                        uint count2 = count > 6 ? 6 : count;
                        if (!size_only)
                        {
                            rle.Add(16);
                            rle_bits.Add(count2 - 3);
                        }
                        clcounts[16]++;
                        count -= count2;
                    }
                }

                /* No or insufficient repetition */
                clcounts[symbol] += count;
                while (count > 0)
                {
                    if (!size_only)
                    {
                        rle.Add(symbol);
                        rle_bits.Add(0);
                    }
                    count--;
                }
            }

            ZopfliCalculateBitLengths(clcounts, 19, 7, clcl);
            if (!size_only) ZopfliLengthsToSymbols(clcl, 19, 7, clsymbols);

            hclen = 15;
            /* Trim zeros. */
            while (hclen > 0 && clcounts[EncodeTreeOrder[hclen + 4 - 1]] == 0) hclen--;

            if (!size_only)
            {
                AddBits(hlit, 5, bp, OutFile);
                AddBits(hdist, 5, bp, OutFile);
                AddBits(hclen, 4, bp, OutFile);

                for (i = 0; i < hclen + 4; i++)
                {
                    AddBits(clcl[EncodeTreeOrder[i]], 3, bp, OutFile);
                }

                for (int ii = 0; ii < (int)rle.Count; ii++)
                {
                    uint symbol = clsymbols[rle[ii]];
                    AddHuffmanBits(symbol, clcl[rle[ii]], bp, OutFile);
                    /* Extra bits. */
                    if (rle[ii] == 16) AddBits(rle_bits[ii], 2, bp, OutFile);
                    else if (rle[ii] == 17) AddBits(rle_bits[ii], 3, bp, OutFile);
                    else if (rle[ii] == 18) AddBits(rle_bits[ii], 7, bp, OutFile);
                }
            }

            result_size += 14;  /* hlit, hdist, hclen bits */
            result_size += (hclen + 4) * 3;  /* clcl bits */
            for (i = 0; i < 19; i++)
            {
                result_size += clcl[i] * clcounts[i];
            }
            /* Extra bits. */
            result_size += clcounts[16] * 2;
            result_size += clcounts[17] * 3;
            result_size += clcounts[18] * 7;

            /* Note: in case of "size_only" these are null pointers so no effect. */

            return result_size;
        }

        static void AddDynamicTree(uint[] ll_lengths,
                           uint[] d_lengths,
                           byte[] bp,
                           List<byte> OutFile)
        {
            int i;
            int best = 0;
            ulong bestsize = 0;

            for (i = 0; i < 8; i++)
            {
                ulong size = EncodeTree(ll_lengths, d_lengths,
                                         Convert.ToBoolean(i & 1), Convert.ToBoolean(i & 2), Convert.ToBoolean(i & 4),
                                         bp, OutFile, true);
                if (bestsize == 0 || size < bestsize)
                {
                    bestsize = size;
                    best = i;
                }
            }

            EncodeTree(ll_lengths, d_lengths,
                       Convert.ToBoolean(best & 1), Convert.ToBoolean(best & 2), Convert.ToBoolean(best & 4),
                       bp, OutFile, false);
        }

        /*
        Gives the exact size of the tree, in bits, as it will be encoded in DEFLATE.
        */
        static ulong CalculateTreeSize(uint[] ll_lengths,
                                uint[] d_lengths)
        {
            ulong result = 0;
            int i;
            //these are only used to have all parameters for EncodeTree
            byte[] DummyBP = new byte[1];
            List<byte> DummyOutFIle = new List<byte>();

            for (i = 0; i < 8; i++)
            {
                ulong size = EncodeTree(ll_lengths, d_lengths,
                                         Convert.ToBoolean(i & 1), Convert.ToBoolean(i & 2), Convert.ToBoolean(i & 4),
                                         DummyBP, DummyOutFIle, true);
                if (result == 0 || size < result) result = size;
            }

            return result;
        }

        /*
        Adds all lit/len and dist codes from the lists as huffman symbols. Does not add
        end code 256. expected_data_size is the uncompressed block size, used for
        assert, but you can set it to 0 to not do the assertion.
        */
        static void AddLZ77Data(ZopfliLZ77Store lz77,
                                ulong lstart, ulong lend,
                                ulong expected_data_size,
                        uint[] ll_symbols, uint[] ll_lengths,
                        uint[] d_symbols, uint[] d_lengths,
                        byte[] bp,
                        List<byte> OutFile)
        {
            ulong testlength = 0;
            ulong i;

            for (i = lstart; i < lend; i++)
            {
                uint dist = lz77.dists[(int)i];
                uint litlen = lz77.litlens[(int)i];
                if (dist == 0)
                {
                    Debug.Assert(litlen < 256);
                    Debug.Assert(ll_lengths[litlen] > 0);
                    AddHuffmanBits(ll_symbols[litlen], ll_lengths[litlen], bp, OutFile);
                    testlength++;
                }
                else
                {
                    uint lls = Symbols.ZopfliGetLengthSymbol((ushort)litlen);
                    uint ds = Symbols.ZopfliGetDistSymbol((ushort)dist);
                    Debug.Assert(litlen >= 3 && litlen <= 288);
                    Debug.Assert(ll_lengths[lls] > 0);
                    Debug.Assert(d_lengths[ds] > 0);
                    AddHuffmanBits(ll_symbols[lls], ll_lengths[lls], bp, OutFile);
                    AddBits(Symbols.ZopfliGetLengthExtraBitsValue((ushort)litlen),
                            (uint)Symbols.ZopfliGetLengthExtraBits((int)litlen),
                            bp, OutFile);
                    AddHuffmanBits(d_symbols[ds], d_lengths[ds], bp, OutFile);
                    AddBits((uint)Symbols.ZopfliGetDistExtraBitsValue((int)dist),
                            (uint)Symbols.ZopfliGetDistExtraBits((int)dist),
                            bp, OutFile);
                    testlength += litlen;
                }
            }
            Debug.Assert(expected_data_size == 0 || testlength == expected_data_size);
        }

        /* Since an uncompressed block can be max 65535 in size, it actually adds
        multible blocks if needed. */
        static void AddNonCompressedBlock(bool final,
                                          byte[] InFile, ulong instart,
                                          ulong inend,
                                          byte[] bp,
                                          List<byte> OutFile)
        {
            ulong pos = instart;
            for (; ; )
            {
                ulong i;
                ushort blocksize = 65535;
                ushort nlen;
                bool currentfinal;

                if (pos + blocksize > inend) blocksize = (ushort)(inend - pos);
                currentfinal = pos + blocksize >= inend;

                nlen = (ushort)~blocksize;

                AddBit(Convert.ToInt32(final && currentfinal), bp, OutFile);
                /* BTYPE 00 */
                AddBit(0, bp, OutFile);
                AddBit(0, bp, OutFile);

                /* Any bits of input up to the next byte boundary are ignored. */
                bp[0] = 0;

                OutFile.Add((byte)(blocksize % 256));
                OutFile.Add((byte)(blocksize / 256));
                OutFile.Add((byte)(nlen % 256));
                OutFile.Add((byte)((nlen / 256) % 256));

                for (i = 0; i < blocksize; i++)
                {
                    OutFile.Add(InFile[pos + i]);
                }

                if (currentfinal) break;
                pos += blocksize;
            }
        }

        /*
        Adds a deflate block with the given LZ77 data to the output.
        options: global program options
        btype: the block type, must be 1 or 2
        final: whether to set the "final" bit on this block, must be the last block
        litlens: literal/length array of the LZ77 data, in the same format as in
            ZopfliLZ77Store.
        dists: distance array of the LZ77 data, in the same format as in
            ZopfliLZ77Store.
        lstart: where to start in the LZ77 data
        lend: where to end in the LZ77 data (not inclusive)
        expected_data_size: the uncompressed block size, used for assert, but you can
          set it to 0 to not do the assertion.
        bp: output bit pointer
        out: dynamic output array to append to
        outsize: dynamic output array size
        */
        static void AddLZ77Block(int btype, bool final,
                         ZopfliLZ77Store lz77,
                         ulong lstart, ulong lend,
                         ulong expected_data_size,
                         byte[] bp,
                         List<byte> OutFile)
        {
            uint[] ll_lengths = new uint[ZOPFLI_NUM_LL];
            uint[] d_lengths = new uint[ZOPFLI_NUM_D];
            uint[] ll_symbols = new uint[ZOPFLI_NUM_LL];
            uint[] d_symbols = new uint[ZOPFLI_NUM_D];
            int detect_block_size /* header size*/;
            ulong compressed_size;
            ulong uncompressed_size = 0;
            ulong i;
            if (btype == 0)
            {
                ulong length = ZopfliLZ77GetByteRange(lz77, lstart, lend);
                ulong pos = lstart == lend ? 0 : lz77.pos[(int)lstart];
                ulong end = pos + length;
                AddNonCompressedBlock(final,
                                      lz77.data, pos, end, bp, OutFile);
                return;
            }

            AddBit(Convert.ToInt32(final), bp, OutFile);
            AddBit(btype & 1, bp, OutFile);
            AddBit((btype & 2) >> 1, bp, OutFile);

            if (btype == 1)
            {
                /* Fixed block. */
                GetFixedTree(ll_lengths, d_lengths);
            }
            else
            {
                /* Dynamic block. */
                uint detect_tree_size;
                Debug.Assert(btype == 2);

                GetDynamicLengths(lz77, lstart, lend, ll_lengths, d_lengths);

                detect_tree_size = (uint)OutFile.Count + 10 /* header size*/;
                AddDynamicTree(ll_lengths, d_lengths, bp, OutFile);
                if (Globals.verbose > 0)
                {
                    Console.WriteLine ("treesize: " + (int)(OutFile.Count - detect_tree_size));
                }
            }

            ZopfliLengthsToSymbols(ll_lengths, ZOPFLI_NUM_LL, 15, ll_symbols);
            ZopfliLengthsToSymbols(d_lengths, ZOPFLI_NUM_D, 15, d_symbols);

            detect_block_size = OutFile.Count + 10 /* header size*/;
            AddLZ77Data(lz77, lstart, lend, expected_data_size,
                        ll_symbols, ll_lengths, d_symbols, d_lengths,
                        bp, OutFile);
            /* End symbol. */
            AddHuffmanBits(ll_symbols[256], ll_lengths[256], bp, OutFile);

            for (i = lstart; i < lend; i++)
            {
                uncompressed_size += (ulong)(lz77.dists[(int)i] == 0 ? 1 : lz77.litlens[(int)i]);
            }
            compressed_size = (ulong)(OutFile.Count + 10 /* header size*/ - detect_block_size);
            if (Globals.verbose > 0)
            {
                Console.WriteLine("compressed block size: %d (%dk) (unc: %d)" + compressed_size + "(" + compressed_size / 1024 + "k) (unc: "
                       + uncompressed_size + ")");
            }
        }

        static void AddLZ77BlockAutoType(bool final,
                                 ZopfliLZ77Store lz77,
                                 ulong lstart, ulong lend,
                                 ulong expected_data_size,
                                 byte[] bp,
                                 List<byte> OutFilePart)
        {

            double uncompressedcost = ZopfliCalculateBlockSize(lz77, lstart, lend, 0);
            double fixedcost = ZopfliCalculateBlockSize(lz77, lstart, lend, 1);
            double dyncost = ZopfliCalculateBlockSize(lz77, lstart, lend, 2);

            /* Whether to perform the expensive calculation of creating an optimal block
            with fixed huffman tree to check if smaller. Only do this for small blocks or
            blocks which already are pretty good with fixed huffman tree. */
            bool expensivefixed = (lz77.size < 1000) || fixedcost <= dyncost * 1.1;

            ZopfliLZ77Store fixedstore = new ZopfliLZ77Store(lz77.data);
            if (lstart == lend)
            {
                /* Smallest empty block is represented by fixed block */
                AddBits(Convert.ToUInt32(final), 1, bp, OutFilePart);
                AddBits(1, 2, bp, OutFilePart);  /* btype 01 */
                AddBits(0, 7, bp, OutFilePart);  /* end symbol has code 0000000 */
                return;
            }
            if (expensivefixed)
            {
                /* Recalculate the LZ77 with ZopfliLZ77OptimalFixed */
                ulong instart = lz77.pos[(int)lstart];
                ulong inend = instart + ZopfliLZ77GetByteRange(lz77, lstart, lend);

                ZopfliBlockState s = new ZopfliBlockState((int)instart, (int)inend, 1);
                ZopfliLZ77OptimalFixed(s, lz77.data, instart, inend, fixedstore);
                fixedcost = ZopfliCalculateBlockSize(fixedstore, 0, fixedstore.size, 1);
            }

            if (uncompressedcost < fixedcost && uncompressedcost < dyncost)
            {
                AddLZ77Block(0, final, lz77, lstart, lend,
                             expected_data_size, bp, OutFilePart);
            }
            else if (fixedcost < dyncost)
            {
                if (expensivefixed)
                {
                    AddLZ77Block(1, final, fixedstore, 0, fixedstore.size,
                                 expected_data_size, bp, OutFilePart);
                }
                else
                {
                    AddLZ77Block(1, final, lz77, lstart, lend,
                                 expected_data_size, bp, OutFilePart);
                }
            }
            else
            {
                AddLZ77Block(2, final, lz77, lstart, lend,
                             expected_data_size, bp, OutFilePart);
            }

        }

        static void ZopfliBlockSplitLZ77(ZopfliLZ77Store lz77, List<ulong> splitpoints, out ulong npoints)
        {
            ulong lstart, lend;
            ulong llpos;
            int numblocks = 1;
            byte[] done;
            double splitcost, origcost;

            npoints = 0;
            if (lz77.size < 10)
            {
                return;  /* This code fails on tiny files. */
            }

            done = new byte[lz77.size];

            lstart = 0;
            lend = lz77.size;
            for (; ; )
            {
                SplitCostContext c = new SplitCostContext();

                if (Globals.blocksplittingmax > 0 && numblocks >= Globals.blocksplittingmax)
                {
                    break;
                }

                c.lz77 = lz77;
                c.start = lstart;
                c.end = lend;
                Debug.Assert(lstart < lend);
                llpos = FindMinimum(c, lstart + 1, lend, out splitcost);

                Debug.Assert(llpos > lstart);
                Debug.Assert(llpos < lend);

                origcost = EstimateCost(lz77, lstart, lend);

                if (splitcost > origcost || llpos == lstart + 1 || llpos == lend)
                {
                    done[lstart] = 1;
                }
                else
                {
                    AddSorted(llpos, splitpoints, ref npoints);
                    numblocks++;
                }
                
                if (FindLargestSplittableBlock(
                    lz77.size, done, splitpoints, npoints, ref lstart, ref lend) == 0)
                {
                    break;  /* No further split will probably reduce compression. */
                }

                if (lend - lstart < 10)
                {
                    break;
                }
            }

            if (Globals.verbose> 0)
            {
                PrintBlockSplitPoints(lz77, splitpoints, npoints);
            }
        }

        static void ZopfliBlockSplit(byte[] InFile, int instart, int inend, List<ulong> splitpoints)
        {
            int pos;
            ulong i;
            ZopfliBlockState s = new ZopfliBlockState(instart, inend, 0);
            List<ulong> lz77splitpoints = new List<ulong>();
            ulong nlz77points;
            ZopfliLZ77Store store = new ZopfliLZ77Store(InFile);
            ZopfliHash h = new ZopfliHash();

            /* Unintuitively, Using a simple LZ77 method here instead of ZopfliLZ77Optimal
            results in better blocks. */
            ZopfliLZ77Greedy(s, InFile, instart, inend, store, h);

            ZopfliBlockSplitLZ77(store, lz77splitpoints, out nlz77points );

            /* Convert LZ77 positions to positions in the uncompressed input. */
            pos = instart;
            if (nlz77points > 0)
            {
                for (i = 0; i < store.size; i++)
                {
                    ulong length = (ulong)(store.dists[(int)i] == 0 ? 1 : store.litlens[(int)i]);
                    if (lz77splitpoints[splitpoints.Count] == i)
                    {
                        splitpoints.Add((ulong)pos);
                        if (splitpoints.Count == (int)nlz77points) break;
                    }
                    pos += (int)length;
                }
            }
            Debug.Assert(lz77splitpoints.Count == (int)nlz77points);

        }

        /* Since an uncompressed block can be max 65535 in size, it actually adds
        multible blocks if needed. */
        static void AddNonCompressedBlock(bool final,
                                  byte[] InFile, int instart,
                                  int inend,
                                  byte[] bp,
                                  List<byte> OutFilePart)
        {
            int pos = instart;
            for (; ; )
            {
                int i;
                ushort blocksize = 65535;
                ushort nlen;
                bool currentfinal;

                if (pos + blocksize > inend) blocksize = (ushort)(inend - pos);
                currentfinal = pos + blocksize >= inend;

                nlen = (ushort)~blocksize;

                AddBit(Convert.ToInt32(final && currentfinal), bp, OutFilePart);
                /* BTYPE 00 */
                AddBit(0, bp, OutFilePart);
                AddBit(0, bp, OutFilePart);

                /* Any bits of input up to the next byte boundary are ignored. */
                bp[0] = 0;

                OutFilePart.Add((byte)(blocksize % 256));
                OutFilePart.Add((byte)((blocksize / 256) % 256));
                OutFilePart.Add((byte)(nlen % 256));
                OutFilePart.Add((byte)((nlen / 256) % 256));

                for (i = 0; i < blocksize; i++)
                {
                    OutFilePart.Add(InFile[pos + i]);
                }

                if (currentfinal) break;
                pos += blocksize;
            }
        }

        /*
        Deflate a part, to allow ZopfliDeflate() to use multiple master blocks if
        needed.
        It is possible to call this function multiple times in a row, shifting
        instart and inend to next bytes of the data. If instart is larger than 0, then
        previous bytes are used as the initial dictionary for LZ77.
        This function will usually output multiple deflate blocks. If final is 1, then
        the final bit will be set on the last block.
        */
        static void ZopfliDeflatePart(int btype, bool final, byte[] InFile, int instart, int inend , byte[] bp, List<byte> OutFilePart)
        {
            int i;
            /* byte coordinates rather than lz77 index */
            List<ulong> splitpoints_uncompressed = new List<ulong>();
            ulong npoints = 0;
            List<ulong> splitpoints = new List<ulong>();
            double totalcost = 0;
            ZopfliLZ77Store lz77 = new ZopfliLZ77Store(InFile);

            /* If btype=2 is specified, it tries all block types. If a lesser btype is
            given, then however it forces that one. Neither of the lesser types needs
            block splitting as they have no dynamic huffman trees. */
            if (btype == 0)
            {
                AddNonCompressedBlock(final, InFile, instart, inend, bp, OutFilePart);
                return;
            }
            else if (btype == 1)
            {
                ZopfliBlockState s = new ZopfliBlockState(instart, inend, 1);
                ZopfliLZ77Store store = new ZopfliLZ77Store(InFile);

                ZopfliLZ77OptimalFixed(s, store.data, (ulong)instart, (ulong)inend, store);
                AddLZ77Block( btype, final, store, 0, store.size, 0,
                             bp, OutFilePart);

                return;
            }

            if (Globals.blocksplitting == 1)
            {
                ZopfliBlockSplit( InFile, instart, inend, splitpoints_uncompressed);
                npoints = (ulong)splitpoints_uncompressed.Count;

                for (i = 0; i <= (int)npoints; i++)
                {
                    int start = i == 0 ? instart : (int)splitpoints_uncompressed[i - 1];
                    int end = i == (int)npoints ? inend : (int)splitpoints_uncompressed[i];
                    ZopfliBlockState s = new ZopfliBlockState(start, end, 1);
                    ZopfliLZ77Store store = new ZopfliLZ77Store(InFile);
                    ZopfliLZ77Optimal(s, InFile, start, end, Globals.numiterations, store);
                    totalcost += ZopfliCalculateBlockSizeAutoType(store, 0, store.size);

                    ZopfliAppendLZ77Store(store, lz77);
                    if (i < (int)npoints) splitpoints.Add(lz77.size);

                }
            }

            /* Second block splitting attempt */
            if (Globals.blocksplitting == 1 && npoints > 1)
            {
                List<ulong> splitpoints2 = new List<ulong>();
                ulong npoints2;
                double totalcost2 = 0;

                ZopfliBlockSplitLZ77(lz77, splitpoints2, out npoints2);

                for (i = 0; i <= (int)npoints2; i++)
                {
                    ulong start = i == 0 ? 0 : splitpoints2[i - 1];
                    ulong end = i == (int)npoints2 ? lz77.size : splitpoints2[i];
                    totalcost2 += ZopfliCalculateBlockSizeAutoType(lz77, start, end);
                }

                if (totalcost2 < totalcost)
                {
                    splitpoints = new List<ulong>(splitpoints2);
                    npoints = npoints2;
                }
            }

            for (i = 0; i <= (int)npoints; i++)
            {
                ulong start = i == 0 ? 0 : splitpoints[i - 1];
                ulong end = i == (int)npoints ? lz77.size : splitpoints[i];
                AddLZ77BlockAutoType(i == (int)npoints && final == true,
                                     lz77, start, end, 0, bp, OutFilePart);
            }

        }
    }
}
