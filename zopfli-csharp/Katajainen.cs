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
        Nodes forming chains. Also used to represent leaves.
        */
        class Node
        {
            public uint weight;  /* Total weight (symbol count) of this chain. */
            public Node tail;  /* Previous node(s) of this chain, or 0 if none. */
            public int count;  /* Leaf symbol index, or number of leaves before this chain. */
        };

        static void InitNode(uint weight, int count, Node tail, Node node)
        {
            node.weight = weight;
            node.count = count;
            node.tail = tail;
        }

        /*
        Performs a Boundary Package-Merge step. Puts a new chain in the given list. The
        new chain is, depending on the weights, a leaf or a combination of two chains
        from the previous list.
        lists: The lists of chains.
        maxbits: Number of lists.
        leaves: The leaves, one per symbol.
        numsymbols: Number of leaves.
        pool: the node memory pool.
        index: The index of the list in which a new chain or leaf is required.
        */
        static void BoundaryPM(Node[ , ]lists, Node[] leaves, int numsymbols,
                               Span<Node> pool, ref int PoolNext, int index)
        {
            Node newchain;
            Node oldchain;
            int lastcount = lists[index, 1].count;  /* Count of last chain of list. */


            newchain = pool[PoolNext++];
            oldchain = lists[index, 1];

            /* These are set up before the recursive calls below, so that there is a list
            pointing to the new node, to let the garbage collection know it's in use. */
            lists[index, 0] = oldchain;
            lists[index, 1] = newchain;

            if (index == 0)
            {
                if (lastcount >= numsymbols) return;
                /* New leaf node in list 0. */
                InitNode(leaves[lastcount].weight, lastcount + 1, null, newchain);
            }
            else
            {
                uint sum = lists[index - 1, 0].weight + lists[index - 1, 1].weight;
                if (lastcount < numsymbols && sum > leaves[lastcount].weight)
                {
                    /* New leaf inserted in list, so count is incremented. */
                    InitNode(leaves[lastcount].weight, lastcount + 1, oldchain.tail,
                        newchain);
                }
                else
                {
                    InitNode(sum, lastcount, lists[index - 1, 1], newchain);
                    /* Two lookahead chains of previous list used up, create new ones. */
                    BoundaryPM(lists, leaves, numsymbols, pool, ref PoolNext, index - 1);
                    BoundaryPM(lists, leaves, numsymbols, pool, ref PoolNext, index - 1);
                }
            }

        }

        static void BoundaryPMFinal(Node[,] lists, Node[] leaves, int numsymbols,
                               Span<Node> pool, ref int PoolNext, int index)
        {
            int lastcount = lists[index, 1].count;  /* Count of last chain of list. */

            ulong sum = lists[index - 1, 0].weight + lists[index - 1, 1].weight;

            if (lastcount < numsymbols && sum > leaves[lastcount].weight)
            {
                Node newchain = pool[PoolNext];
                Node oldchain = lists[index, 1].tail;

                lists[index, 1] = newchain;
                newchain.count = lastcount + 1;
                newchain.tail = oldchain;
            }
            else
            {
                lists[index, 1].tail = lists[index - 1, 1];
            }
        }

        /*
        Initializes each list with as lookahead chains the two leaves with lowest
        weights.
        */
        static void InitLists(
            Span<Node> pool, ref int PoolNext, Node[] leaves, int maxbits, Node[,] lists)
        {
            int i;
            Node node0 = pool[PoolNext++];
            Node node1 = pool[PoolNext++];
            InitNode(leaves[0].weight, 1, null, node0);
            InitNode(leaves[1].weight, 2, null, node1);
            for (i = 0; i < maxbits; i++)
            {
                lists[i, 0] = node0;
                lists[i, 1] = node1;
            }
        }

        /*
        Converts result of boundary package-merge to the bitlengths. The result in the
        last chain of the last list contains the amount of active leaves in each list.
        chain: Chain to extract the bit length from (last chain from last list).
        */
        static void ExtractBitLengths(Node chain, Node[] leaves, uint[] bitlengths)
        {
            int[] counts = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

            uint end = 16;
            uint ptr = 15;
            uint value = 1;
            Node node;
            int val;

            for (node = chain; node != null; node = node.tail)
            {
                counts[--end] = node.count;
            }

            val = counts[15];
            while (ptr >= end)
            {
                for (; val > counts[ptr - 1]; val--)
                {
                    bitlengths[leaves[val - 1].count] = value;
                }
                ptr--;
                value++;
            }
        }

        static Node[] _nodes;

        static Span<Node> InitializeNodes(int length)
        {
            if (_nodes is null || _nodes.Length < length)
            {
                Node[] temp_nodes;
                temp_nodes = new Node[length];
                for (int i = 0; i < length; i++)
                {
                    if(!(_nodes is null) && _nodes.Length > i)
                    {
                        temp_nodes[i] = _nodes[i];
                    }
                    else
                    {
                        temp_nodes[i] = new Node();
                    }
                 }

                _nodes = temp_nodes;

            }
            return new Span<Node>(_nodes, 0, length);
        }

        static Node[] _leaves;

        static Span<Node> InitializeLeaves(int length)
        {
            if (_leaves is null || _leaves.Length < length)
            {
                Node[] temp_nodes;
                temp_nodes = new Node[length];
                for (int i = 0; i < length; i++)
                {
                    if (!(_leaves is null) && _leaves.Length > i)
                    {
                        temp_nodes[i] = _leaves[i];
                    }
                    else
                    {
                        temp_nodes[i] = new Node();
                    }
                }

                _leaves = temp_nodes;

            }
            return new Span<Node>(_leaves, 0, length);
        }

        public static int ZopfliLengthLimitedCodeLengths(
            uint[] frequencies, int n, int maxbits, uint[] bitlengths)
        {
            int PoolNext = 0;
            int i;
            int numsymbols = 0;  /* Amount of symbols with frequency > 0. */
            int numBoundaryPMRuns;
            Span<Node> nodes;

            /* Array of lists of chains. Each list requires only two lookahead chains at
            a time, so each list is a array of two Node*'s. */
            Node[,] lists;

            /* One leaf per symbol. Only numsymbols leaves will be used. */
            //Node[] leaves = Helpers.InitializeArray<Node>(n);

            /* Initialize all bitlengths at 0. */
            bitlengths.Initialize();

            /* Count used symbols and place them in the leaves. */
            for (i = 0; i < n; i++)
            {
                if (frequencies[i] > 0)
                {
                    numsymbols++;
                }
            }

            Helpers.DebugCounter++;

            Node[] leaves = InitializeLeaves(numsymbols).ToArray();
            int SymbolNumber = 0;

            /* Count used symbols and place them in the leaves. */
            for (i = 0; i < n; i++)
            {

                if (frequencies[i] > 0)
                {
                    leaves[SymbolNumber].weight = frequencies[i];
                    leaves[SymbolNumber].count = i;  /* Index of symbol this leaf represents. */
                    SymbolNumber++;
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
            if (numsymbols == 1)
            {
                bitlengths[leaves[0].count] = 1;
                return 0;  /* Only one symbol, give it bitlength 1, not 0. OK. */
            }
            if (numsymbols == 2)
            {
                bitlengths[leaves[0].count]++;
                bitlengths[leaves[1].count]++;
                return 0;
            }

            /* Sort the leaves from lightest to heaviest. Add count into the same
            variable for stable sorting. */
            for (i = 0; i < numsymbols; i++)
            {
                if (leaves[i].weight >=
                    (ulong)(1 << (int)(sizeof(ulong) * leaves[0].weight * CHAR_BIT - 9))) 
                {
                    return 1;  /* Error, we need 9 bits for the count. */
                }
                leaves[i].weight = (leaves[i].weight << 9) | (uint)leaves[i].count;
            }
            Array.Sort(leaves, delegate (Node x, Node y) { return x.weight.CompareTo(y.weight); }) ;
            for (i = 0; i < numsymbols; i++)
            {
                leaves[i].weight >>= 9;
            }

            if (numsymbols - 1 < maxbits)
            {
                maxbits = numsymbols - 1;
            }

            /* Initialize node memory pool. */
            nodes = InitializeNodes(maxbits * 2 * numsymbols);

            lists = new Node[maxbits, 2];
            InitLists(nodes, ref PoolNext, leaves, maxbits, lists);

            /* In the last list, 2 * numsymbols - 2 active chains need to be created. Two
            are already created in the initialization. Each BoundaryPM run creates one. */
            numBoundaryPMRuns = 2 * numsymbols - 4;
            for (i = 0; i < numBoundaryPMRuns - 1; i++)
            {
                BoundaryPM(lists, leaves, numsymbols, nodes, ref PoolNext, maxbits - 1);
            }
            BoundaryPMFinal(lists, leaves, numsymbols, nodes, ref PoolNext, maxbits - 1);

            ExtractBitLengths(lists[maxbits - 1, 1], leaves, bitlengths);

            return 0;  /* OK. */

        }

    }
}
