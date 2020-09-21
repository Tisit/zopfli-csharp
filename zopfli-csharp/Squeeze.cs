using System;
using System.Collections.Generic;
using System.Diagnostics;
using zopfli_csharp;

namespace ZopfliCSharp
{
    public partial class Compress
    {
        static void CopyStats(SymbolStats source, ref SymbolStats dest)
        {
            Array.Copy(source.dists, dest.dists, source.dists.Length);
            Array.Copy(source.litlens, dest.litlens, source.litlens.Length);
            Array.Copy(source.ll_symbols, dest.ll_symbols, source.ll_symbols.Length);
            Array.Copy(source.d_symbols, dest.d_symbols, source.d_symbols.Length);
        }

        /* Adds the bit lengths. */
        static void AddWeighedStatFreqs( SymbolStats stats1, double w1,
                                 SymbolStats stats2, double w2,
                                SymbolStats result)
        {
            ulong i;
            for (i = 0; i < ZOPFLI_NUM_LL; i++)
            {
                result.litlens[i] = (ulong)(stats1.litlens[i] * w1 + stats2.litlens[i] * w2);
            }
            for (i = 0; i < ZOPFLI_NUM_D; i++)
            {
                result.dists[i] = (ulong)(stats1.dists[i] * w1 + stats2.dists[i] * w2);
            }
            result.litlens[256] = 1;  /* End symbol. */
        }

        static void ClearStatFreqs(SymbolStats stats)
        {
            ulong i;
            for (i = 0; i < ZOPFLI_NUM_LL; i++) stats.litlens[i] = 0;
            for (i = 0; i < ZOPFLI_NUM_D; i++) stats.dists[i] = 0;
        }

        /* Get random number: "Multiply-With-Carry" generator of G. Marsaglia */
        static uint Ran(RanState state)
        {
            state.m_z = 36969 * (state.m_z & 65535) + (state.m_z >> 16);
            state.m_w = 18000 * (state.m_w & 65535) + (state.m_w >> 16);
            return (state.m_z << 16) + state.m_w;  /* 32-bit result. */
        }

        static void RandomizeFreqs(RanState state, ulong[] freqs, int n)
        {
            int i;
            for (i = 0; i < n; i++)
            {
                if ((Ran(state) >> 4) % 3 == 0) freqs[i] = freqs[Ran(state) % n];
            }
        }

        static void RandomizeStatFreqs(RanState state, SymbolStats stats)
        {
            RandomizeFreqs(state, stats.litlens, ZOPFLI_NUM_LL);
            RandomizeFreqs(state, stats.dists, ZOPFLI_NUM_D);
            stats.litlens[256] = 1;  /* End symbol. */
        }

        /* Appends the symbol statistics from the store. */
        static void GetStatistics(ZopfliLZ77Store store, SymbolStats stats)
        {
            ulong i;
            for (i = 0; i < store.size; i++)
            {
                if (store.dists[(int)i] == 0)
                {
                    stats.litlens[store.litlens[(int)i]]++;
                }
                else
                {
                    stats.litlens[Symbols.ZopfliGetLengthSymbol(store.litlens[(int)i])]++;
                    stats.dists[Symbols.ZopfliGetDistSymbol(store.dists[(int)i])]++;
                }
            }
            stats.litlens[256] = 1;  /* End symbol. */

            stats.CalculateStatistics();
        }

        static double GetCostFixed(int litlen, int dist)
        {
            if (dist == 0)
            {
                if (litlen <= 143) return 8;
                else return 9;
            }
            else
            {
                int dbits = Symbols.ZopfliGetDistExtraBits(dist);
                int lbits = Symbols.ZopfliGetLengthExtraBits(litlen);
                int lsym = Symbols.ZopfliGetLengthSymbol((ushort)litlen);
                int cost = 0;
                if (lsym <= 279) cost += 7;
                else cost += 8;
                cost += 5;  /* Every dist symbol has length 5. */
                return dbits + lbits + cost;
            }
        }

        /*
        Cost model based on symbol statistics.
        type: CostModelFun
        */
        static double GetCostStat(int litlen, int dist, SymbolStats stats)
        {
            if (dist == 0)
            {
                return stats.ll_symbols[litlen];
            }
            else
            {
                int dbits = Symbols.ZopfliGetDistExtraBits(dist);
                int lbits = Symbols.ZopfliGetLengthExtraBits(litlen);
                int lsym = Symbols.ZopfliGetLengthSymbol((ushort)litlen);
                int dsym = Symbols.ZopfliGetDistSymbol((ushort)dist);
                return dbits + lbits + stats.ll_symbols[lsym] + stats.d_symbols[dsym];
            }
        }

        /*
        Finds the minimum possible cost this cost model can return for valid length and
        distance symbols.
        */
        static double GetCostModelMinCost(SymbolStats stats)
        {
            double mincost;
            int bestlength = 0; /* length that has lowest cost in the cost model */
            int bestdist = 0; /* distance that has lowest cost in the cost model */
            int i;
            /*
            Table of distances that have a different distance symbol in the deflate
            specification. Each value is the first distance that has a new symbol. Only
            different symbols affect the cost model so only these need to be checked.
            See RFC 1951 section 3.2.5. Compressed blocks (length and distance codes).
            */
            int[] dsymbols = {
                    1, 2, 3, 4, 5, 7, 9, 13, 17, 25, 33, 49, 65, 97, 129, 193, 257, 385, 513,
                    769, 1025, 1537, 2049, 3073, 4097, 6145, 8193, 12289, 16385, 24577
                  };

            mincost = ZOPFLI_LARGE_FLOAT;
            for (i = 3; i < 259; i++)
            {
                double c = GetCostStat(i, 1, stats);
                if (c < mincost)
                {
                    bestlength = i;
                    mincost = c;
                }
            }

            mincost = ZOPFLI_LARGE_FLOAT;
            for (i = 0; i < 30; i++)
            {
                double c = GetCostStat(3, dsymbols[i], stats);
                if (c < mincost)
                {
                    bestdist = dsymbols[i];
                    mincost = c;
                }
            }

            return GetCostStat(bestlength, bestdist, stats);
        }

        static ulong zopfli_min(ulong a, ulong b)
        {
            return a < b ? a : b;
        }

        /*
        Performs the forward pass for "squeeze". Gets the most optimal length to reach
        every byte from a previous byte, using cost calculations.
        s: the ZopfliBlockState
        in: the input data array
        instart: where to start
        inend: where to stop (not inclusive)
        costmodel: function to calculate the cost of some lit/len/dist pair.
        costcontext: abstract context for the costmodel function
        length_array: output array of size (inend - instart) which will receive the best
            length to reach this byte from a previous byte.
        returns the cost that was, according to the costmodel, needed to get to the end.
        */
        static double GetBestLengths(ZopfliBlockState s,
                             byte[] InFile,
                             int instart, int inend,
                             SymbolStats stats,
                             ushort[] length_array,
                             ZopfliHash h, float[] costs, bool fixedcosts)
        {
            /* Best cost to get here so far. */
            ulong blocksize = (ulong)(inend - instart);
            ulong i, k, kend;
            ushort leng = 0; //bogus value
            ushort dist = 0; //bogusvalue
            ushort[] sublen = new ushort[259];
            ulong windowstart = (ulong)(instart > ZopfliHash.ZOPFLI_WINDOW_SIZE
                ? instart - ZopfliHash.ZOPFLI_WINDOW_SIZE : 0);
            double result;
            double mincost = GetCostModelMinCost(stats);
            double mincostaddcostj;

            if (instart == inend) return 0;

            h.ZopfliResetHash();
            h.ZopfliWarmupHash(InFile, (int)windowstart, inend);
            for (i = windowstart; i < (ulong)instart; i++)
            {
                h.ZopfliUpdateHash(InFile, (int)i, inend);
            }

            Array.Fill(costs, (float)Compress.ZOPFLI_LARGE_FLOAT);
            costs[0] = 0;  /* Because it's the start. */
            length_array[0] = 0;

            for (i = (ulong)instart; i < (ulong)inend; i++)
            {
                ulong j = i - (ulong)instart;  /* Index in the costs array and length_array. */
                h.ZopfliUpdateHash(InFile, (int)i, inend);

                /* If we're in a long repetition of the same character and have more than
                ZOPFLI_MAX_MATCH characters before and after our position. */
                if (h.same[i & ZopfliHash.ZOPFLI_WINDOW_MASK] > ZOPFLI_MAX_MATCH * 2
                    && (int)i > instart + ZOPFLI_MAX_MATCH + 1
                    && i + ZOPFLI_MAX_MATCH * 2 + 1 < (ulong)inend
                    && h.same[(i - ZOPFLI_MAX_MATCH) & ZopfliHash.ZOPFLI_WINDOW_MASK]
                        > ZOPFLI_MAX_MATCH)
                {
                    double symbolcost;

                    if (fixedcosts)
                    {
                        symbolcost = GetCostFixed(ZOPFLI_MAX_MATCH, 1);
                    } else
                    {
                        symbolcost = GetCostStat(ZOPFLI_MAX_MATCH, 1, stats);
                    }
                    /* Set the length to reach each one to ZOPFLI_MAX_MATCH, and the cost to
                    the cost corresponding to that length. Doing this, we skip
                    ZOPFLI_MAX_MATCH values to avoid calling ZopfliFindLongestMatch. */
                    for (k = 0; k < ZOPFLI_MAX_MATCH; k++)
                    {
                        costs[j + ZOPFLI_MAX_MATCH] = costs[j] + (float)symbolcost;
                        length_array[j + ZOPFLI_MAX_MATCH] = ZOPFLI_MAX_MATCH;
                        i++;
                        j++;
                        h.ZopfliUpdateHash(InFile, (int)i, inend);
                    }
                }

                ZopfliFindLongestMatch(s, h, InFile, (int)i, inend, ZOPFLI_MAX_MATCH, sublen,
                                        ref dist,  ref leng);

                /* Literal. */
                if (i + 1 <= (ulong)inend)
                {
                    double newCost;
                    if (fixedcosts)
                    {
                        newCost = GetCostFixed(InFile[i], 0) + costs[j];
                    }
                    else
                    {
                        newCost = GetCostStat(InFile[i], 0, stats) + costs[j];
                    }
                    Debug.Assert(newCost >= 0);
                    if (newCost < costs[j + 1])
                    {
                        costs[j + 1] = (float)newCost;
                        length_array[j + 1] = 1;
                    }
                }
                /* Lengths. */
                kend = zopfli_min(leng, (ulong)inend - i);
                mincostaddcostj = mincost + costs[j];
                for (k = 3; k <= kend; k++)
                {
                    double newCost;

                    /* Calling the cost model is expensive, avoid this if we are already at
                    the minimum possible cost that it can return. */
                    if (costs[j + k] <= mincostaddcostj) continue;

                    if (fixedcosts)
                    {
                        newCost = GetCostFixed((int)k, sublen[k]) + costs[j];
                    }
                    else
                    {
                        newCost = GetCostStat((int)k, sublen[k], stats) + costs[j];
                    }
                    Debug.Assert(newCost >= 0);
                    if (newCost < costs[j + k])
                    {
                        Debug.Assert(k <= ZOPFLI_MAX_MATCH);
                        costs[j + k] = (float)newCost;
                        length_array[j + k] = (ushort)k;
                    }
                }
            }

            Debug.Assert(costs[blocksize] >= 0);
            result = costs[blocksize];

            return result;
        }

        /*
        Calculates the optimal path of lz77 lengths to use, from the calculated
        length_array. The length_array must contain the optimal length to reach that
        byte. The path will be filled with the lengths to use, so its data size will be
        the amount of lz77 symbols.
        */
        static void TraceBackwards(ulong size, ushort[] length_array,
                                   List<ushort> path, ref ulong  pathsize)
        {
            int index = (int)size;
            if (size == 0) return;
            for (; ; )
            {
                path.Add(length_array[index]);
                pathsize++;
                Debug.Assert(length_array[index] <= index);
                Debug.Assert(length_array[index] <= ZOPFLI_MAX_MATCH);
                Debug.Assert(length_array[index] != 0);
                index -= length_array[index];
                if (index == 0) break;
            }

            /* Mirror result. */
            for (index = 0; index < (int)pathsize / 2; index++)
            {
                ushort temp = path[index];
                path[index] = path[(int)pathsize - index - 1];
                path[(int)pathsize - index - 1] = temp;
            }
        }

        static void FollowPath(ZopfliBlockState s,
                       byte[] InFile, int instart, int inend,
                       List<ushort> path, ulong pathsize,
                       ZopfliLZ77Store store, ZopfliHash h)
        {
            int i, j, pos;
            int windowstart = instart > ZopfliHash.ZOPFLI_WINDOW_SIZE
                ? instart - ZopfliHash.ZOPFLI_WINDOW_SIZE : 0;

            ulong total_length_test = 0;

            if (instart == inend) return;

            h.ZopfliResetHash();
            h.ZopfliWarmupHash(InFile, windowstart, inend);
            for (i = windowstart; i < instart; i++)
            {
                h.ZopfliUpdateHash(InFile, i, inend);
            }

            pos = instart;
            for (i = 0; i < (int)pathsize; i++)
            {
                ushort length = path[i];
                ushort dummy_length = 0;
                ushort dist = 0;
                Debug.Assert(pos < inend);

                h.ZopfliUpdateHash(InFile, pos, inend);

                /* Add to output. */
                if (length >= ZOPFLI_MIN_MATCH)
                {
                    /* Get the distance by recalculating longest match. The found length
                    should match the length from the path. */
                    ZopfliFindLongestMatch(s, h, InFile, pos, inend, length, null,
                                            ref dist,  ref dummy_length);
                    Debug.Assert(!(dummy_length != length && length > 2 && dummy_length > 2));
                    ZopfliVerifyLenDist(InFile, inend, pos, dist, length);
                    ZopfliStoreLitLenDist(length, dist, pos, store);
                    total_length_test += length;
                }
                else
                {
                    length = 1;
                    ZopfliStoreLitLenDist(InFile[pos], 0, pos, store);
                    total_length_test++;
                }


                Debug.Assert(pos + length <= inend);
                for (j = 1; j < length; j++)
                {
                    h.ZopfliUpdateHash(InFile, pos + j, inend);
                }

                pos += length;
            }
        }

        static void ZopfliCopyLZ77Store(
            ZopfliLZ77Store source, ZopfliLZ77Store dest)
        {
            dest.ResetStore(source.data);
            dest.litlens = new List<ushort>(source.litlens);
            dest.ll_counts = new List<ulong>(source.ll_counts);
            dest.d_counts = new List<ulong>(source.d_counts);
            dest.ll_symbol = new List<ushort>(source.ll_symbol);
            dest.d_symbol = new List<ushort>(source.d_symbol);
            dest.pos = new List<ulong>(source.pos);
            dest.dists = new List<ushort>(source.dists);
            dest.size = source.size;
        }

        /*
        Does a single run for ZopfliLZ77Optimal. For good compression, repeated runs
        with updated statistics should be performed.
        s: the block state
        in: the input data array
        instart: where to start
        inend: where to stop (not inclusive)
        path: pointer to dynamically allocated memory to store the path
        pathsize: pointer to the size of the dynamic path array
        length_array: array of size (inend - instart) used to store lengths
        costmodel: function to use as the cost model for this squeeze run
        costcontext: abstract context for the costmodel function
        store: place to output the LZ77 data
        returns the cost that was, according to the costmodel, needed to get to the end.
            This is not the actual cost.
        */
        static double LZ77OptimalRun(ZopfliBlockState s,
            byte[] InFile, int instart, int inend, List<ushort> path,
            ushort[] length_array, SymbolStats stats, ZopfliLZ77Store store,
            ZopfliHash h, float[] costs, bool fixedcosts)
        {
            double cost = GetBestLengths(s, InFile, instart, inend, stats, length_array, h, costs, fixedcosts);
            ulong pathsize = 0;
            path.Clear();
            TraceBackwards((ulong)(inend - instart), length_array, path, ref pathsize);
            FollowPath(s, InFile, instart, inend, path, pathsize, store, h);
            Debug.Assert(cost < ZOPFLI_LARGE_FLOAT);
            return cost;
        }

        static void ZopfliLZ77Optimal(ZopfliBlockState s,
                       byte[] InFile, int instart, int inend,
                       int numiterations,
                       ZopfliLZ77Store store)
        {
            /* Dist to get to here with smallest cost. */
            int blocksize = inend - instart;
            ushort[] length_array = new ushort[blocksize + 1];
            List<ushort> path = new List<ushort>();
            ZopfliLZ77Store currentstore = new ZopfliLZ77Store(InFile);
            ZopfliHash h = new ZopfliHash();
            SymbolStats stats = new SymbolStats();
            SymbolStats beststats = new SymbolStats();
            SymbolStats laststats = new SymbolStats();
            int i;
            float[] costs = new float [blocksize + 1];
            double cost;
            double bestcost = ZOPFLI_LARGE_FLOAT;
            double lastcost = 0;
            /* Try randomizing the costs a bit once the size stabilizes. */
            RanState ran_state = new RanState();
            int lastrandomstep = -1;

            /* Do regular deflate, then loop multiple shortest path runs, each time using
            the statistics of the previous run. */

            /* Initial run. */
            ZopfliLZ77Greedy(s, InFile, instart, inend, currentstore, h);
            GetStatistics(currentstore, stats);

            /* Repeat statistics with each time the cost model from the previous stat
            run. */
            for (i = 0; i < numiterations; i++)
            {
                currentstore.ResetStore(InFile);
                LZ77OptimalRun(s, InFile, instart, inend, path, length_array, stats,
                               currentstore, h, costs, false);
                cost = ZopfliCalculateBlockSize(currentstore, 0, currentstore.size, 2);
                if (Globals.verbose_more > 0 || (Globals.verbose > 0 && cost < bestcost))
                {
                    Console.WriteLine("Iteration " + i + ": " + cost + " bit");
                }
                if (cost < bestcost)
                {
                    /* Copy to the output store. */
                    ZopfliCopyLZ77Store(currentstore, store);
                    CopyStats(stats, ref beststats);
                    bestcost = cost;
                }
                CopyStats(stats, ref laststats);
                ClearStatFreqs(stats);
                GetStatistics(currentstore, stats);
                if (lastrandomstep != -1)
                {
                    /* This makes it converge slower but better. Do it only once the
                    randomness kicks in so that if the user does few iterations, it gives a
                    better result sooner. */
                    AddWeighedStatFreqs(stats, 1.0, laststats, 0.5, stats);
                    stats.CalculateStatistics();
                }
                if (i > 5 && cost == lastcost)
                {
                    CopyStats(beststats, ref stats);
                    RandomizeStatFreqs(ran_state, stats);
                    stats.CalculateStatistics();
                    lastrandomstep = i;
                }
                lastcost = cost;
            }

        }

        static void ZopfliLZ77OptimalFixed(ZopfliBlockState s,
                            byte[] InFile,
                            ulong instart, ulong inend,
                            ZopfliLZ77Store store)
        {
            /* Dist to get to here with smallest cost. */
            ulong blocksize = inend - instart;
            ushort[] length_array = new ushort[blocksize + 1];
            List<ushort> path = new List<ushort>();
            path.Add(0);
            ZopfliHash h = new ZopfliHash();
            float[] costs = new float[blocksize + 1];
            SymbolStats stats = new SymbolStats();

            s.blockstart = (int)instart;
            s.blockend = (int)inend;

            /* Shortest path for fixed tree This one should give the shortest possible
            result for fixed tree, no repeated runs are needed since the tree is known. */
            LZ77OptimalRun(s, InFile, (int)instart, (int)inend, path,
                           length_array, stats, store, h, costs, true);


        }

    }
}
