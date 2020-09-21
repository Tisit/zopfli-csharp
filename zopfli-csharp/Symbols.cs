using System;
using System.Collections.Generic;
using System.Text;

namespace zopfli_csharp
{
    /*
    Utilities for using the lz77 symbols of the deflate spec.
    */
    static class Symbols
    {
        static readonly int[] ExtraBits = {
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1,
                2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
                3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
                3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
                4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
                4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
                4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
                4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
                5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
                5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
                5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
                5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
                5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
                5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
                5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
                5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 0
            };

        static readonly ushort[] ExtraBitsValue = {
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 2, 3, 0,
                    1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 4, 5, 6, 7, 0, 1, 2, 3, 4, 5,
                    6, 7, 0, 1, 2, 3, 4, 5, 6, 7, 0, 1, 2, 3, 4, 5, 6, 7, 0, 1, 2, 3, 4, 5, 6,
                    7, 8, 9, 10, 11, 12, 13, 14, 15, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12,
                    13, 14, 15, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 0, 1, 2,
                    3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,
                    10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28,
                    29, 30, 31, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17,
                    18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 0, 1, 2, 3, 4, 5, 6,
                    7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26,
                    27, 28, 29, 30, 31, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15,
                    16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 0
                };

        static readonly ushort[] LengthSymbol = {
                    0, 0, 0, 257, 258, 259, 260, 261, 262, 263, 264,
                    265, 265, 266, 266, 267, 267, 268, 268,
                    269, 269, 269, 269, 270, 270, 270, 270,
                    271, 271, 271, 271, 272, 272, 272, 272,
                    273, 273, 273, 273, 273, 273, 273, 273,
                    274, 274, 274, 274, 274, 274, 274, 274,
                    275, 275, 275, 275, 275, 275, 275, 275,
                    276, 276, 276, 276, 276, 276, 276, 276,
                    277, 277, 277, 277, 277, 277, 277, 277,
                    277, 277, 277, 277, 277, 277, 277, 277,
                    278, 278, 278, 278, 278, 278, 278, 278,
                    278, 278, 278, 278, 278, 278, 278, 278,
                    279, 279, 279, 279, 279, 279, 279, 279,
                    279, 279, 279, 279, 279, 279, 279, 279,
                    280, 280, 280, 280, 280, 280, 280, 280,
                    280, 280, 280, 280, 280, 280, 280, 280,
                    281, 281, 281, 281, 281, 281, 281, 281,
                    281, 281, 281, 281, 281, 281, 281, 281,
                    281, 281, 281, 281, 281, 281, 281, 281,
                    281, 281, 281, 281, 281, 281, 281, 281,
                    282, 282, 282, 282, 282, 282, 282, 282,
                    282, 282, 282, 282, 282, 282, 282, 282,
                    282, 282, 282, 282, 282, 282, 282, 282,
                    282, 282, 282, 282, 282, 282, 282, 282,
                    283, 283, 283, 283, 283, 283, 283, 283,
                    283, 283, 283, 283, 283, 283, 283, 283,
                    283, 283, 283, 283, 283, 283, 283, 283,
                    283, 283, 283, 283, 283, 283, 283, 283,
                    284, 284, 284, 284, 284, 284, 284, 284,
                    284, 284, 284, 284, 284, 284, 284, 284,
                    284, 284, 284, 284, 284, 284, 284, 284,
                    284, 284, 284, 284, 284, 284, 284, 285
                };

        static readonly int[] LengthSymbolExtraBits = {
                0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2,
                3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 0
            };

        static readonly int[] DistSymbolExtraBits = {
                0, 0, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8,
                9, 9, 10, 10, 11, 11, 12, 12, 13, 13
            };

        /* Gets the amount of extra bits for the given dist, cfr. the DEFLATE spec. */
        public static int ZopfliGetDistExtraBits(int dist)
        {
            if (dist < 5) return 0;
            else if (dist < 9) return 1;
            else if (dist < 17) return 2;
            else if (dist < 33) return 3;
            else if (dist < 65) return 4;
            else if (dist < 129) return 5;
            else if (dist < 257) return 6;
            else if (dist < 513) return 7;
            else if (dist < 1025) return 8;
            else if (dist < 2049) return 9;
            else if (dist < 4097) return 10;
            else if (dist < 8193) return 11;
            else if (dist < 16385) return 12;
            else return 13;
        }

        /* Gets value of the extra bits for the given dist, cfr. the DEFLATE spec. */
        public static int ZopfliGetDistExtraBitsValue(int dist)
        {
            if (dist < 5) return 0;
            else if (dist < 9) return (dist - 5) & 1;
            else if (dist < 17) return (dist - 9) & 3;
            else if (dist < 33) return (dist - 17) & 7;
            else if (dist < 65) return (dist - 33) & 15;
            else if (dist < 129) return (dist - 65) & 31;
            else if (dist < 257) return (dist - 129) & 63;
            else if (dist < 513) return (dist - 257) & 127;
            else if (dist < 1025) return (dist - 513) & 255;
            else if (dist < 2049) return (dist - 1025) & 511;
            else if (dist < 4097) return (dist - 2049) & 1023;
            else if (dist < 8193) return (dist - 4097) & 2047;
            else if (dist < 16385) return (dist - 8193) & 4095;
            else return (dist - 16385) & 8191;
        }

        /* Gets the symbol for the given dist, cfr. the DEFLATE spec. */
        public static ushort ZopfliGetDistSymbol(ushort dist)
        {
            if (dist < 193)
            {
                if (dist < 13)
                {  /* dist 0..13. */
                    if (dist < 5) return (ushort)(dist - 1);
                    else if (dist < 7) return 4;
                    else if (dist < 9) return 5;
                    else return 6;
                }
                else
                {  /* dist 13..193. */
                    if (dist < 17) return 7;
                    else if (dist < 25) return 8;
                    else if (dist < 33) return 9;
                    else if (dist < 49) return 10;
                    else if (dist < 65) return 11;
                    else if (dist < 97) return 12;
                    else if (dist < 129) return 13;
                    else return 14;
                }
            }
            else
            {
                if (dist < 2049)
                {  /* dist 193..2049. */
                    if (dist < 257) return 15;
                    else if (dist < 385) return 16;
                    else if (dist < 513) return 17;
                    else if (dist < 769) return 18;
                    else if (dist < 1025) return 19;
                    else if (dist < 1537) return 20;
                    else return 21;
                }
                else
                {  /* dist 2049..32768. */
                    if (dist < 3073) return 22;
                    else if (dist < 4097) return 23;
                    else if (dist < 6145) return 24;
                    else if (dist < 8193) return 25;
                    else if (dist < 12289) return 26;
                    else if (dist < 16385) return 27;
                    else if (dist < 24577) return 28;
                    else return 29;
                }
            }
        }

        /* Gets the amount of extra bits for the given length, cfr. the DEFLATE spec. */
        public static int ZopfliGetLengthExtraBits(int l)
        {
            return ExtraBits[l];
        }

        /* Gets value of the extra bits for the given length, cfr. the DEFLATE spec. */
        public static ushort ZopfliGetLengthExtraBitsValue(ushort l)
        {
           
            return ExtraBitsValue[l];
        }
        /*
        Gets the symbol for the given length, cfr. the DEFLATE spec.
        Returns the symbol in the range [257-285] (inclusive)
        */
        public static ushort ZopfliGetLengthSymbol(ushort l)
        {
           

            return LengthSymbol[l];
        }

        /* Gets the amount of extra bits for the given length symbol. */
        public static int ZopfliGetLengthSymbolExtraBits(int s) 
        {
            
            return LengthSymbolExtraBits[s - 257];
        }

        /* Gets the amount of extra bits for the given distance symbol. */
        public static int ZopfliGetDistSymbolExtraBits(int s) {
            
            return DistSymbolExtraBits[s];
        }
    }
}
