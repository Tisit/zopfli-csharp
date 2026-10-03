using System;
using System.Collections.Generic;
using System.Text;
using ZopfliCSharp;

namespace zopfli_csharp
{
    static class Katajainen
    {
        const int CHAR_BIT = 8;

        /*
        Nodes forming chains, stored as parallel arrays (structure of arrays) and
        referenced by index, with -1 as the null chain. Compared to a pool of Node
        objects this avoids pointer chasing and the covariant type check on every
        reference store into a Node[].
        weight: total weight (symbol count) of this chain.
        tail: previous node(s) of this chain, or -1 if none.
        count: number of leaves before this chain.
        The pool persists across calls and only grows, so a node can still hold data
        from an earlier call when it is handed out; the result does not depend on it
        (the C original uses uninitialized malloc memory here).
        */
        static uint[] _weight = Array.Empty<uint>();
        static int[] _tail = Array.Empty<int>();
        static int[] _count = Array.Empty<int>();

        static void EnsurePool(int length)
        {
            if (_weight.Length >= length) return;
            Array.Resize(ref _weight, length);
            Array.Resize(ref _tail, length);
            Array.Resize(ref _count, length);
        }

        static void InitNode(uint weight, int count, int tail, int node)
        {
            _weight[node] = weight;
            _count[node] = count;
            _tail[node] = tail;
        }

        /*
        Performs a Boundary Package-Merge step. Puts a new chain in the given list. The
        new chain is, depending on the weights, a leaf or a combination of two chains
        from the previous list.
        lists: The lists of chains.
        leafWeights: The weights of the leaves, one per symbol, lightest first.
        numsymbols: Number of leaves.
        PoolNext: the next free node of the node memory pool.
        index: The index of the list in which a new chain or leaf is required.
        */
        /* `lists` is a flattened [maxbits, 2] array of node indices: element [i, j]
           lives at [i * 2 + j]. */
        static void BoundaryPM(int[] lists, uint[] leafWeights, int numsymbols,
                               ref int PoolNext, int index)
        {
            int oldchain = lists[index * 2 + 1];
            int lastcount = _count[oldchain];  /* Count of last chain of list. */
            int newchain = PoolNext++;

            lists[index * 2] = oldchain;
            lists[index * 2 + 1] = newchain;

            if (index == 0)
            {
                if (lastcount >= numsymbols) return;
                /* New leaf node in list 0. */
                InitNode(leafWeights[lastcount], lastcount + 1, -1, newchain);
            }
            else
            {
                uint sum = _weight[lists[(index - 1) * 2]] + _weight[lists[(index - 1) * 2 + 1]];
                if (lastcount < numsymbols && sum > leafWeights[lastcount])
                {
                    /* New leaf inserted in list, so count is incremented. */
                    InitNode(leafWeights[lastcount], lastcount + 1, _tail[oldchain], newchain);
                }
                else
                {
                    InitNode(sum, lastcount, lists[(index - 1) * 2 + 1], newchain);
                    /* Two lookahead chains of previous list used up, create new ones. */
                    BoundaryPM(lists, leafWeights, numsymbols, ref PoolNext, index - 1);
                    BoundaryPM(lists, leafWeights, numsymbols, ref PoolNext, index - 1);
                }
            }
        }

        static void BoundaryPMFinal(int[] lists, uint[] leafWeights, int numsymbols,
                               ref int PoolNext, int index)
        {
            int lastcount = _count[lists[index * 2 + 1]];  /* Count of last chain of list. */

            ulong sum = _weight[lists[(index - 1) * 2]] + _weight[lists[(index - 1) * 2 + 1]];

            if (lastcount < numsymbols && sum > leafWeights[lastcount])
            {
                int newchain = PoolNext;
                int oldchain = _tail[lists[index * 2 + 1]];

                lists[index * 2 + 1] = newchain;
                _count[newchain] = lastcount + 1;
                _tail[newchain] = oldchain;
            }
            else
            {
                _tail[lists[index * 2 + 1]] = lists[(index - 1) * 2 + 1];
            }
        }

        /*
        Initializes each list with as lookahead chains the two leaves with lowest
        weights.
        */
        static void InitLists(ref int PoolNext, uint[] leafWeights, int maxbits, int[] lists)
        {
            int i;
            int node0 = PoolNext++;
            int node1 = PoolNext++;
            InitNode(leafWeights[0], 1, -1, node0);
            InitNode(leafWeights[1], 2, -1, node1);
            for (i = 0; i < maxbits; i++)
            {
                lists[i * 2] = node0;
                lists[i * 2 + 1] = node1;
            }
        }

        /*
        Converts result of boundary package-merge to the bitlengths. The result in the
        last chain of the last list contains the amount of active leaves in each list.
        chain: Chain to extract the bit length from (last chain from last list).
        leafSymbols: The symbol each leaf represents, in leaf order.
        */
        static void ExtractBitLengths(int chain, int[] leafSymbols, uint[] bitlengths)
        {
            Span<int> counts = stackalloc int[16];

            int end = 16;
            int ptr = 15;
            uint value = 1;
            int val;

            for (int node = chain; node != -1; node = _tail[node])
            {
                counts[--end] = _count[node];
            }

            val = counts[15];
            while (ptr >= end)
            {
                for (; val > counts[ptr - 1]; val--)
                {
                    bitlengths[leafSymbols[val - 1]] = value;
                }
                ptr--;
                value++;
            }
        }

        public static int ZopfliLengthLimitedCodeLengths(
            uint[] frequencies, int n, int maxbits, uint[] bitlengths)
        {
            int PoolNext = 0;
            int i;
            int numsymbols = 0;  /* Amount of symbols with frequency > 0. */
            int numBoundaryPMRuns;

            /* Array of lists of chains. Each list requires only two lookahead chains at
            a time, so each list is a array of two node indices. Flattened to a 1-D int[]
            of length maxbits*2 (element [i, j] at [i * 2 + j]). */
            int[] lists;

            /* Initialize all bitlengths at 0. */
            bitlengths.Initialize();

            /* Count used symbols. */
            for (i = 0; i < n; i++)
            {
                if (frequencies[i] > 0)
                {
                    numsymbols++;
                }
            }

            /* Check special cases and error conditions. */
            if ((1 << maxbits) < numsymbols)
            {
                return 1;  /* Error, too few maxbits to represent symbols. */
            }
            if (numsymbols == 0)
            {
                return 0;  /* No symbols at all. OK. */
            }

            /* Place the used symbols in the leaves, in symbol order. */
            int[] leafSymbols = new int[numsymbols];
            int SymbolNumber = 0;
            for (i = 0; i < n; i++)
            {
                if (frequencies[i] > 0)
                {
                    leafSymbols[SymbolNumber++] = i;
                }
            }

            if (numsymbols == 1)
            {
                bitlengths[leafSymbols[0]] = 1;
                return 0;  /* Only one symbol, give it bitlength 1, not 0. OK. */
            }
            if (numsymbols == 2)
            {
                bitlengths[leafSymbols[0]]++;
                bitlengths[leafSymbols[1]]++;
                return 0;
            }

            /* Sort the leaves from lightest to heaviest. Add count into the same
            variable for stable sorting. */
            /* weight is a 32-bit uint and the low 9 bits are used to pack the count for
               stable sorting, so a weight must fit in the remaining 23 bits. (The C
               original uses size_t weights and thus a 55-bit limit; here it is 23.)
               The packed keys are unique, so sorting them as plain integers yields
               exactly the order of the original comparison sort. */
            const uint MAX_SORT_WEIGHT = 1u << (32 - 9);
            uint[] keys = new uint[numsymbols];
            for (i = 0; i < numsymbols; i++)
            {
                uint weight = frequencies[leafSymbols[i]];
                if (weight >= MAX_SORT_WEIGHT)
                {
                    return 1;  /* Error, we need 9 bits for the count. */
                }
                keys[i] = (weight << 9) | (uint)leafSymbols[i];
            }
            Array.Sort(keys);
            uint[] leafWeights = new uint[numsymbols];
            for (i = 0; i < numsymbols; i++)
            {
                leafWeights[i] = keys[i] >> 9;
                leafSymbols[i] = (int)(keys[i] & 511);
            }

            if (numsymbols - 1 < maxbits)
            {
                maxbits = numsymbols - 1;
            }

            /* Initialize node memory pool. */
            EnsurePool(maxbits * 2 * numsymbols);

            lists = new int[maxbits * 2];
            InitLists(ref PoolNext, leafWeights, maxbits, lists);

            /* In the last list, 2 * numsymbols - 2 active chains need to be created. Two
            are already created in the initialization. Each BoundaryPM run creates one. */
            numBoundaryPMRuns = 2 * numsymbols - 4;
            for (i = 0; i < numBoundaryPMRuns - 1; i++)
            {
                BoundaryPM(lists, leafWeights, numsymbols, ref PoolNext, maxbits - 1);
            }
            BoundaryPMFinal(lists, leafWeights, numsymbols, ref PoolNext, maxbits - 1);

            ExtractBitLengths(lists[(maxbits - 1) * 2 + 1], leafSymbols, bitlengths);

            return 0;  /* OK. */

        }

    }
}
