using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Numerics;
using System.Runtime.CompilerServices;
using System.Runtime.InteropServices;
using System.Runtime.Intrinsics;
using System.Runtime.Intrinsics.X86;
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
            Array.Clear(stats.litlens);
            Array.Clear(stats.dists);
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
        static double GetCostModelMinCost(SymbolStats stats, bool fixedcosts)
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
                double c = fixedcosts ? GetCostFixed(i, 1) : GetCostStat(i, 1, stats);
                if (c < mincost)
                {
                    bestlength = i;
                    mincost = c;
                }
            }

            mincost = ZOPFLI_LARGE_FLOAT;
            for (i = 0; i < 30; i++)
            {
                double c = fixedcosts ? GetCostFixed(3, dsymbols[i]) : GetCostStat(3, dsymbols[i], stats);
                if (c < mincost)
                {
                    bestdist = dsymbols[i];
                    mincost = c;
                }
            }

            return fixedcosts ? GetCostFixed(bestlength, bestdist) : GetCostStat(bestlength, bestdist, stats);
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
                             ZopfliHashChains chains, float[] costs, bool fixedcosts)
        {
            /* Best cost to get here so far. */
            ulong blocksize = (ulong)(inend - instart);
            ulong i, k, kend;
            ushort leng = 0; //bogus value
            ushort dist = 0; //bogusvalue
            double result;
            double mincost = GetCostModelMinCost(stats, fixedcosts);
            double mincostaddcostj;

            if (instart == inend) return 0;

            /* Precompute the length- and distance-dependent parts of the cost model.
               For dist > 0 the cost is separable: cost(len, dist) = costLen[len] +
               costDist[dist]. This turns the inner cost loop into two array lookups
               instead of repeated symbol/extra-bit function calls per position. */
            double[] costLen = new double[ZOPFLI_MAX_MATCH + 1];
            double[] costDist = new double[ZopfliHash.ZOPFLI_WINDOW_SIZE + 1];
            for (int len = ZOPFLI_MIN_MATCH; len <= ZOPFLI_MAX_MATCH; len++)
            {
                int lbits = Symbols.ZopfliGetLengthExtraBits(len);
                int lsym = Symbols.ZopfliGetLengthSymbol((ushort)len);
                costLen[len] = fixedcosts
                    ? lbits + (lsym <= 279 ? 7 : 8)
                    : lbits + stats.ll_symbols[lsym];
            }
            /* Single-precision copy of costLen for the pre-filter in UpdateCostsForRun. */
            float[] costLenF = new float[ZOPFLI_MAX_MATCH + 1];
            for (int len = 0; len <= ZOPFLI_MAX_MATCH; len++)
            {
                costLenF[len] = (float)costLen[len];
            }
            for (int d = 1; d <= ZopfliHash.ZOPFLI_WINDOW_SIZE; d++)
            {
                int dbits = Symbols.ZopfliGetDistExtraBits(d);
                costDist[d] = fixedcosts
                    ? dbits + 5
                    : dbits + stats.d_symbols[Symbols.ZopfliGetDistSymbol((ushort)d)];
            }

            ushort[] same = chains.same;
            int cstart = chains.start;

            Array.Fill(costs, (float)Compress.ZOPFLI_LARGE_FLOAT);
            costs[0] = 0;  /* Because it's the start. */
            length_array[0] = 0;

            /* Element-0 refs for the inner length loop (UpdateCostsForRun). Every index
               there is provably in range under the loop invariants (j + k <= j + kend <=
               blocksize < costs/length_array length; k <= ZOPFLI_MAX_MATCH < costLen
               length), but the JIT can't prove it (kend is a runtime min), so it would
               emit a bounds check per access in this, the hottest loop in the program.
               Indexing through these refs elides them; the Debug.Asserts enforce each
               invariant in Debug builds (compiled out in Release, same pattern as
               GetMatch). These arrays are never reallocated, and managed byrefs are
               GC-tracked, so hoisting is safe. */
            ref float costs0 = ref MemoryMarshal.GetArrayDataReference(costs);
            ref ushort la0 = ref MemoryMarshal.GetArrayDataReference(length_array);
            ref double costLen0 = ref MemoryMarshal.GetArrayDataReference(costLen);
            ref float costLenF0 = ref MemoryMarshal.GetArrayDataReference(costLenF);

            /* The current match as runs of equal distance. */
            MatchRuns runs = new MatchRuns();
            int[] runEnd = runs.end;
            ushort[] runDist = runs.dist;

            for (i = (ulong)instart; i < (ulong)inend; i++)
            {
                ulong j = i - (ulong)instart;  /* Index in the costs array and length_array. */

                /* If we're in a long repetition of the same character and have more than
                ZOPFLI_MAX_MATCH characters before and after our position. */
                if (same[(int)i - cstart] > ZOPFLI_MAX_MATCH * 2
                    && (int)i > instart + ZOPFLI_MAX_MATCH + 1
                    && i + ZOPFLI_MAX_MATCH * 2 + 1 < (ulong)inend
                    && same[(int)i - ZOPFLI_MAX_MATCH - cstart] > ZOPFLI_MAX_MATCH)
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
                    }
                }

                /* Get the match as runs of equal distance: straight from the longest
                   match cache when possible (the common case after the first pass),
                   otherwise by searching. */
                if (!s.TryGetCachedRuns((int)i, runs, ref leng))
                {
                    ZopfliFindLongestMatch(s, chains, InFile, (int)i, inend, ZOPFLI_MAX_MATCH, runs,
                                            ref dist, ref leng);
                }
                int nruns = runs.count;

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
                kend = Math.Min(leng, (ulong)inend - i);
                float costsj = costs[j];  /* Loop-invariant; also feeds each newCost. */
                mincostaddcostj = mincost + costsj;
                /* Within a run the distance, and so its cost, is constant, which lets
                   the per-length update be vectorized. */
                Debug.Assert((int)j + (int)kend < costs.Length);
                int lo = ZOPFLI_MIN_MATCH;
                for (int r = 0; r < nruns && lo <= (int)kend; r++)
                {
                    int hi = Math.Min(runEnd[r], (int)kend);
                    Debug.Assert(runDist[r] <= ZopfliHash.ZOPFLI_WINDOW_SIZE);
                    UpdateCostsForRun(ref Unsafe.Add(ref costs0, (nint)j),
                                      ref Unsafe.Add(ref la0, (nint)j), ref costLen0, ref costLenF0,
                                      lo, hi, costDist[runDist[r]], costsj, mincostaddcostj);
                    lo = runEnd[r] + 1;
                }
            }

            Debug.Assert(costs[blocksize] >= 0);
            result = costs[blocksize];

            return result;
        }

        /*
        Relaxes costs[k] / length_array[k] for every match length k in [lo, hi] that
        uses the same distance, whose cost is costDist. costs and lengths point at
        element j (the current position), so index k means position j + k.
        Bit-identical to the scalar form: every update is decided and computed by the
        exact double-precision code (Relax4 or the tail loop); the vector pre-filter
        only skips groups of lengths that provably cannot pass that test.
        */
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        static void UpdateCostsForRun(ref float costs, ref ushort lengths, ref double costLen,
                                      ref float costLenF, int lo, int hi, double costDist,
                                      float costsj, double mincostaddcostj)
        {
            Debug.Assert(lo >= ZOPFLI_MIN_MATCH && hi <= ZOPFLI_MAX_MATCH);
            int k = lo;
            if (Avx.IsSupported)
            {
                Vector256<double> vDist = Vector256.Create(costDist);
                Vector256<double> vCostj = Vector256.Create((double)costsj);
                Vector256<double> vMin = Vector256.Create(mincostaddcostj);

                /* Pre-filter, 8 lengths at a time in single precision. With L = costLen[k],
                   D = costDist and S = costsj (all >= 0), the exact test takes a length
                   only if fl64(fl64(L + D) + S) < cur, so L + D + S < cur * (1 + 2^-51).
                   The estimate est = fl32(fl32(L) + fl32(fl32(D) + S)) has 4 roundings
                   of relative error u = 2^-24 each, so est <= (L + D + S) * (1 + u)^3
                   < cur * (1 + 3.01u) for every length the exact test would take. The
                   threshold fl32(fl32(cur * (1 + 8u)) + 1e-30) is at least
                   cur * (1 + 7u - 8u^2) (the 1e-30 absorbs subnormal rounding, which
                   is absolute rather than relative), so no such length is filtered out.
                   Groups where no lane passes (the vast majority) are skipped. */
                Vector256<float> vC2 = Vector256.Create((float)costDist + costsj);
                Vector256<float> vScale = Vector256.Create(1f + 1f / (1 << 21));
                Vector256<float> vTiny = Vector256.Create(1e-30f);
                for (; k + 7 <= hi; k += 8)
                {
                    Vector256<float> cur = Vector256.LoadUnsafe(ref costs, (nuint)k);
                    Vector256<float> est = Avx.Add(Vector256.LoadUnsafe(ref costLenF, (nuint)k), vC2);
                    Vector256<float> threshold = Avx.Add(Avx.Multiply(cur, vScale), vTiny);
                    int maybe = Avx.MoveMask(Avx.CompareLessThanOrEqual(est, threshold));
                    if (maybe == 0) continue;
                    if ((maybe & 0x0F) != 0)
                        Relax4(ref costs, ref lengths, ref costLen, k, vDist, vCostj, vMin);
                    if ((maybe & 0xF0) != 0)
                        Relax4(ref costs, ref lengths, ref costLen, k + 4, vDist, vCostj, vMin);
                }
                for (; k + 3 <= hi; k += 4)
                {
                    Relax4(ref costs, ref lengths, ref costLen, k, vDist, vCostj, vMin);
                }
            }
            for (; k <= hi; k++)
            {
                float cur = Unsafe.Add(ref costs, k);
                /* Calling the cost model is expensive, avoid this if we are already at
                the minimum possible cost that it can return. */
                if (cur <= mincostaddcostj) continue;
                double newCost = Unsafe.Add(ref costLen, k) + costDist + costsj;
                Debug.Assert(newCost >= 0);
                if (newCost < cur)
                {
                    Unsafe.Add(ref costs, k) = (float)newCost;
                    Unsafe.Add(ref lengths, k) = (ushort)k;
                }
            }
        }

        /*
        Exact relaxation of the 4 lengths k..k+3, vectorized. Bit-identical to the
        scalar tail loop of UpdateCostsForRun: each lane computes
        (costLen[k] + costDist) + costsj in double, compares against the float cost
        widened to double, and rounds to float on store.
        */
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        static void Relax4(ref float costs, ref ushort lengths, ref double costLen, int k,
                           Vector256<double> vDist, Vector256<double> vCostj,
                           Vector256<double> vMin)
        {
            ref float cp = ref Unsafe.Add(ref costs, k);
            Vector128<float> curf = Vector128.LoadUnsafe(ref cp);
            Vector256<double> cur = Avx.ConvertToVector256Double(curf);
            Vector256<double> newCost = Avx.Add(
                Avx.Add(Vector256.LoadUnsafe(ref costLen, (nuint)k), vDist), vCostj);
            /* Same early-out as the scalar loop (cur <= mincostaddcostj skips). */
            Vector256<double> take = Avx.And(Avx.CompareLessThan(newCost, cur),
                                             Avx.CompareGreaterThan(cur, vMin));
            int bits = Avx.MoveMask(take);
            if (bits == 0) return;

            /* Narrow the 64-bit lane mask to 32-bit lanes (each lane is all-ones or
               zero, so either half will do) and blend the rounded costs in. */
            Vector128<float> take32 = Sse.Shuffle(take.GetLower().AsSingle(),
                                                  take.GetUpper().AsSingle(), 0b10_00_10_00);
            Vector128<float> newf = Avx.ConvertToVector128Single(newCost);
            Sse41.BlendVariable(curf, newf, take32).StoreUnsafe(ref cp);
            do
            {
                int b = BitOperations.TrailingZeroCount(bits);
                Unsafe.Add(ref lengths, k + b) = (ushort)(k + b);
                bits &= bits - 1;
            } while (bits != 0);
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
                       ZopfliLZ77Store store, ZopfliHashChains chains)
        {
            int i, pos;

            if (instart == inend) return;

            pos = instart;
            for (i = 0; i < (int)pathsize; i++)
            {
                ushort length = path[i];
                ushort dummy_length = 0;
                ushort dist = 0;
                Debug.Assert(pos < inend);

                /* Add to output. */
                if (length >= ZOPFLI_MIN_MATCH)
                {
                    /* Get the distance by recalculating longest match. The found length
                    should match the length from the path. */
                    ZopfliFindLongestMatch(s, chains, InFile, pos, inend, length, null,
                                            ref dist,  ref dummy_length);
                    Debug.Assert(!(dummy_length != length && length > 2 && dummy_length > 2));
                    ZopfliVerifyLenDist(InFile, inend, pos, dist, length);
                    ZopfliStoreLitLenDist(length, dist, pos, store);
                }
                else
                {
                    length = 1;
                    ZopfliStoreLitLenDist(InFile[pos], 0, pos, store);
                }

                Debug.Assert(pos + length <= inend);
                pos += length;
            }
        }

        static void ZopfliCopyLZ77Store(
            ZopfliLZ77Store source, ZopfliLZ77Store dest)
        {
            dest.data = source.data;
            dest.size = source.size;

            int n = (int)source.size;
            dest.EnsureMainCapacity(n);
            Array.Copy(source.litlens, dest.litlens, n);
            Array.Copy(source.dists, dest.dists, n);
            Array.Copy(source.pos, dest.pos, n);
            Array.Copy(source.ll_symbol, dest.ll_symbol, n);
            Array.Copy(source.d_symbol, dest.d_symbol, n);

            if (dest.ll_counts.Length < source.ll_counts_size)
                dest.ll_counts = new ulong[source.ll_counts_size];
            Array.Copy(source.ll_counts, dest.ll_counts, source.ll_counts_size);
            dest.ll_counts_size = source.ll_counts_size;

            if (dest.d_counts.Length < source.d_counts_size)
                dest.d_counts = new ulong[source.d_counts_size];
            Array.Copy(source.d_counts, dest.d_counts, source.d_counts_size);
            dest.d_counts_size = source.d_counts_size;
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
            ZopfliHashChains chains, float[] costs, bool fixedcosts)
        {
            double cost = GetBestLengths(s, InFile, instart, inend, stats, length_array, chains, costs, fixedcosts);
            ulong pathsize = 0;
            path.Clear();
            TraceBackwards((ulong)(inend - instart), length_array, path, ref pathsize);
            FollowPath(s, InFile, instart, inend, path, pathsize, store, chains);
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
            /* Every pass over this block sees the same hash chains; build them once. */
            ZopfliHashChains chains = new ZopfliHashChains(InFile, instart, inend);
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
            ZopfliLZ77Greedy(s, InFile, instart, inend, currentstore, chains);
            GetStatistics(currentstore, stats);

            /* Repeat statistics with each time the cost model from the previous stat
            run. */
            for (i = 0; i < numiterations; i++)
            {
                currentstore.ResetStore(InFile);
                LZ77OptimalRun(s, InFile, instart, inend, path, length_array, stats,
                               currentstore, chains, costs, false);
                cost = ZopfliCalculateBlockSize(currentstore, 0, currentstore.size, 2);
                if (Globals.verbose_more > 0 || (Globals.verbose > 0 && cost < bestcost))
                {
                    Console.Error.WriteLine("Iteration " + i + ": " + cost + " bit");
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
            ZopfliHashChains chains = new ZopfliHashChains(InFile, (int)instart, (int)inend);
            float[] costs = new float[blocksize + 1];
            SymbolStats stats = new SymbolStats();

            s.blockstart = (int)instart;
            s.blockend = (int)inend;

            /* Shortest path for fixed tree This one should give the shortest possible
            result for fixed tree, no repeated runs are needed since the tree is known. */
            LZ77OptimalRun(s, InFile, (int)instart, (int)inend, path,
                           length_array, stats, store, chains, costs, true);


        }

    }
}
