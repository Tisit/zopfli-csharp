using System;
using System.Diagnostics;
using System.IO;

namespace ZopfliCSharp
{
    static class Globals
    {
        /* Whether to print output */
        static public int verbose = 0;

        /* Whether to print more detailed output */
        static public int verbose_more = 0;

        /*
        Maximum amount of times to rerun forward and backward pass to optimize LZ77
        compression cost. Good values: 10, 15 for small files, 5 for files over
        several MB in size or it will be too slow.
        */
        static public int numiterations = 15;

        /*
        If true, splits the data in multiple deflate blocks with optimal choice
        for the block boundaries. Block splitting gives better compression. Default:
        true (1).
        */
        static public int blocksplitting = 1;

        /*
        No longer used, left for compatibility.
        */
        static public int blocksplittinglast = 0;

        /*
        Maximum amount of blocks to split into (0 for unlimited, but this can give
        extreme results that hurt compression on some files). Default value: 15.
        */
        static public int blocksplittingmax = 15;

        static public bool output_to_stdout = false;
        static public string filename = "";
        static public string outfilename;
        static public ZopfliFormat output_type = ZopfliFormat.ZOPFLI_FORMAT_GZIP;
    }
    
    static class Program
    {

        static void CompressFile()
        {
            int testsize;
            long InSize;
            byte[] InFile;

            using (FileStream fs = File.Open(Globals.filename, FileMode.Open))
            {
                InSize = fs.Length;
                if (InSize > 2147483647)
                {
                    Console.WriteLine("Files larger than 2GB are not supported.");
                    Environment.Exit(1);
                }

                InFile = new byte[InSize];
                testsize = fs.Read(InFile, 0, (int)InSize);

                if(testsize != (int) InSize)
                {
                    Console.WriteLine("Invalid filename: " + Globals.filename);
                    Environment.Exit(1);
                }
            }

            byte[] OutFile;
            OutFile = Compress.ZopfliCompress(InFile);

            if(Globals.outfilename != "")
            {
                File.WriteAllBytes(Globals.outfilename, OutFile);
            } else
            {
                Console.Write(OutFile);
            }

        }

        static int Main(string[] args)
        {
            

            foreach (string arg in args)
            {
                switch(arg)
                {
                    case "-v":
                        Globals.verbose = 1;
                        break;
                    case "-c":
                        Globals.output_to_stdout = true;
                        break;
                    case "--deflate":
                        Globals.output_type = ZopfliFormat.ZOPFLI_FORMAT_DEFLATE;
                        break;
                    case "--zlib":
                        Globals.output_type = ZopfliFormat.ZOPFLI_FORMAT_ZLIB;
                        break;
                    case "--gzip":
                        Globals.output_type = ZopfliFormat.ZOPFLI_FORMAT_GZIP;
                        break;
                    case "--splitlast":
                        /* ignore*/
                        break;
                    case "-h":
                        Console.WriteLine("Usage: zopfli[OPTION]... FILE...\n" +
                                          "  -h    gives this help\n" +
                                          "  -c    write the result on standard output, instead of disk" +
                                          " filename + '.gz'\n" +
                                          "  -v    verbose mode\n" +
                                          "  --i#  perform # iterations (default 15). More gives" +
                                          " more compression but is slower." +
                                          " Examples: --i10, --i50, --i1000\n!" +
                                          "  --gzip        output to gzip format (default)\n" +
                                          "  --zlib        output to zlib format instead of gzip\n" + 
                                          "  --deflate     output to deflate format instead of gzip\n" +
                                          "  --splitlast   ignored, left for backwards compatibility\n"
                        );
                        return 0;
                    default:
                        if (arg[0] == '-' && arg[1] == '-' && arg[2] == 'i' && arg[3] >= '0' && arg[3] <= '9')
                        {
                            Globals.numiterations = Int32.Parse(arg.Substring(3));
                        }
                        else if (arg[0] != '-')
                        {
                            Globals.filename = arg;
                            if (Globals.output_to_stdout)
                            {
                                Globals.outfilename = "";
                            }
                            else if (Globals.output_type == ZopfliFormat.ZOPFLI_FORMAT_GZIP)
                            {
                                Globals.outfilename = Globals.filename + ".gz";
                            }
                            else if (Globals.output_type == ZopfliFormat.ZOPFLI_FORMAT_ZLIB)
                            {
                                Globals.outfilename = Globals.filename + ".zlib";
                            }
                            else
                            {
                                Debug.Assert(Globals.output_type == ZopfliFormat.ZOPFLI_FORMAT_DEFLATE);
                                Globals.outfilename = Globals.filename + ".deflate";
                            }
                            if (Globals.verbose == 1 && Globals.outfilename != "")
                            {
                                Console.WriteLine("Saving to: " + Globals.outfilename);
                            }
                        }
                        else
                            Console.WriteLine("Unknown parameter!");
                        break;
                }
            }

            if (Globals.filename == "" || Globals.filename is null)
            {
                string AppName = System.Diagnostics.Process.GetCurrentProcess().ProcessName;
                Console.WriteLine("Please provide filename\nFor help, type: " + AppName + " -h");
                return 0;
            }

            if (Globals.numiterations < 1)
            {
                Console.WriteLine("Error: must have 1 or more iterations");
                return 0;
            }

            CompressFile();

            return 0;
        }
    }

    static class Helpers
    {
        //taken from
        //https://stackoverflow.com/questions/3301678/how-to-declare-an-array-of-objects-in-c-sharp
        static public T[] InitializeArray<T>(int length) where T : new()
        {
            T[] array = new T[length];
            for (int i = 0; i < length; ++i)
            {
                array[i] = new T();
            }

            return array;
        }

        //used for debugging
        public static ulong DebugCounter = 0;
    }
}
