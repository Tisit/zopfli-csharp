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
            /* Stream both input and output so the whole file and the whole compressed
               result are never held in memory at once. File I/O errors (missing file,
               permission denied, unwritable target) are reported cleanly rather than
               surfacing as an unhandled exception with a stack trace. */
            try
            {
                using FileStream inStream = new FileStream(Globals.filename, FileMode.Open,
                    FileAccess.Read, FileShare.Read, 1 << 16);

                if (inStream.Length > 2147483647)
                {
                    Console.WriteLine("Files larger than 2GB are not supported.");
                    Environment.Exit(1);
                }

                if (Globals.outfilename == "")
                {
                    /* Streaming to stdout: there is no file to make atomic. */
                    using Stream stdout = Console.OpenStandardOutput();
                    using BufferedStream sbs = new BufferedStream(stdout, 1 << 16);
                    Compress.ZopfliCompress(inStream, sbs);
                    return;
                }

                /* Write to a temporary file in the same directory, force its bytes to
                   physical media, then atomically replace the target. A crash, Ctrl+C,
                   or power loss therefore leaves either the previous output (if any) or
                   the complete new file on disk -- never a partial or corrupt one. The
                   temp is removed on any failure. This uses only portable BCL APIs, so
                   it behaves the same on Windows (FlushFileBuffers + MoveFileEx) and
                   Linux (fsync + rename). */
                string finalPath = Globals.outfilename;
                string tempPath = finalPath + ".tmp" + Environment.ProcessId;

                ConsoleCancelEventHandler onCancel = (s, e) => TryDelete(tempPath);
                Console.CancelKeyPress += onCancel;
                try
                {
                    /* FileShare.Delete lets the cancel handler (or any concurrent
                       cleanup) remove the temp file while this handle is still open.
                       On Windows that marks it delete-pending, so the OS erases it as
                       soon as the handle closes during process teardown -- otherwise a
                       Ctrl+C landing mid-write would leave the temp behind. */
                    using (FileStream outFile = new FileStream(tempPath, FileMode.Create, FileAccess.Write, FileShare.Read | FileShare.Delete))
                    using (BufferedStream bs = new BufferedStream(outFile, 1 << 16))
                    {
                        Compress.ZopfliCompress(inStream, bs);
                        bs.Flush();
                        outFile.Flush(flushToDisk: true);  /* fsync / FlushFileBuffers */
                    }
                    File.Move(tempPath, finalPath, overwrite: true);  /* atomic replace */
                }
                catch
                {
                    TryDelete(tempPath);
                    throw;
                }
                finally
                {
                    Console.CancelKeyPress -= onCancel;
                }
            }
            catch (Exception ex) when (ex is IOException || ex is UnauthorizedAccessException)
            {
                Console.WriteLine("Error: " + ex.Message);
                Environment.Exit(1);
            }
        }

        /* Best-effort removal of a temp file; a leftover temp is harmless clutter. */
        static void TryDelete(string path)
        {
            try { File.Delete(path); }
            catch { /* ignore */ }
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
                        if (arg.Length >= 4 && arg[0] == '-' && arg[1] == '-' && arg[2] == 'i' && arg[3] >= '0' && arg[3] <= '9')
                        {
                            if (!Int32.TryParse(arg.Substring(3), out Globals.numiterations))
                            {
                                Console.WriteLine("Error: invalid iteration count in '" + arg + "'");
                                return 1;
                            }
                        }
                        else if (arg.Length > 0 && arg[0] != '-')
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
}
