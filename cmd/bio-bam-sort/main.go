package main

// bio-bam-sort sorts a BAM file in increasing coordinate order.
//
// Usage: bio-bam-sort input.bam output.bam

import (
	"flag"
	"io"
	"os"
	"runtime"

	"github.com/biogo/hts/bam"
	"github.com/biogo/hts/sam"
	"github.com/grailbio/base/file"
	"github.com/grailbio/base/grail"
	"github.com/grailbio/base/log"
	"github.com/grailbio/base/vcontext"
	"github.com/grailbio/bio/cmd/bio-bam-sort/sorter"
)

var (
	samInputFlag           = flag.Bool("sam", true, "Specify that the inputs are in SAM format")
	shardIndexFlag         = flag.Int("shard-index", 0, "Value of bam.SorterOptions.ShardIndex")
	bamFlag                = flag.String("bam", "", "Merge multiple sortshard files into one BAM file specified by this flag")
	pamFlag                = flag.String("pam", "", "Merge multiple sortshard files into one PAM file specified by this flag")
	parallelismFlag        = flag.Int("parallelism", 64, "Parallelism during PAM generation.")
	recordsPerPAMShardFlag = flag.Int64("records-per-pam-shard", 128<<20,
		"Approx. size of each PAM shard, in number of reads.")
)

// recordReader is implemented by both biogo sam.Reader and biogo bam.Reader.
type recordReader interface {
	Header() *sam.Header
	Read() (*sam.Record, error)
}

// openInput creates a BAM or SAM reader.
func openInput(inPath string) recordReader {
	var in io.Reader
	if inPath == "-" {
		in = os.Stdin
	} else {
		ctx := vcontext.Background()
		f, err := file.Open(ctx, inPath) // Note: f is leaked.
		if err != nil {
			log.Panicf("open %v: %v", inPath, err)
		}
		in = f.Reader(ctx)
	}

	var err error
	var reader recordReader
	if *samInputFlag {
		reader, err = sam.NewReader(in)
		if err != nil {
			log.Panicf("open %v: failed to open SAM: %v", inPath, err)
		}
	} else {
		reader, err = bam.NewReader(in, runtime.NumCPU())
		if err != nil {
			log.Panicf("open %v: failed to open BAM: %v", inPath, err)
		}
	}
	return reader
}

// sort sorts a sequence of sam.Records in inPath to a sortshard file outPath.
func sort(inPath, outPath string) {
	in := openInput(inPath)
	sorter := sorter.NewSorter(outPath, in.Header(), sorter.SortOptions{ShardIndex: uint32(*shardIndexFlag)})
	for nRecs := 0; ; nRecs++ {
		rec, err := in.Read()
		if rec == nil {
			if err != io.EOF {
				log.Panicf("%v: failed to read %dth record: %v", inPath, nRecs, err)
			}
			break
		}
		sorter.AddRecord(rec)
	}
	if err := sorter.Close(); err != nil {
		log.Panicf("Sorter close failed: %v", err)
	}
}

func main() {
	log.SetFlags(log.Ldate | log.Ltime | log.Lmicroseconds | log.Lshortfile)
	flag.Usage = func() {
		os.Stderr.WriteString(`Usage:
This command can either sort a list of BAM/SAM records into a sortshard file, or
merge multiple sortshard files into a BAM/PAM file. It is invoked one of the
following three ways.

1. bio-bam-sort [-sam] <input> <sortshard>

   The command reads a sequence of bam or sam records from input, sorts them,
   and produces file <sortshard>. If <input> is '-', records are read from
   stdin.  With -sam flag, the records are assumed to be in SAM format. Else, it
   is assumed to be in BAM format.

2. bio-bam-sort -bam <foo.bam> <sortshard...>

   The command reads a list of sortshard files and merges them into foo.bam.
   Existing contents of foo.bam, if any, are destroyed.

3. bio-bam-sort -pam <foo.pam> <sortshard...>

   The command reads a list of sortshard files and merges them into foo.pam.
   Existing contents of foo.pam, if any, are destroyed.
`)
		flag.PrintDefaults()
	}
	shutdown := grail.Init()
	defer shutdown()

	args := flag.Args()
	if *bamFlag != "" {
		if len(args) < 1 {
			flag.Usage()
			os.Exit(1)
		}
		err := sorter.BAMFromSortShards(args, *bamFlag)
		if err != nil {
			log.Panicf("merge %v to %v: %v", args, *bamFlag, err)
		}
	} else if *pamFlag != "" {
		if len(args) < 1 {
			flag.Usage()
			os.Exit(1)
		}
		err := sorter.PAMFromSortShards(args, *pamFlag, *recordsPerPAMShardFlag, *parallelismFlag)
		if err != nil {
			log.Panicf("merge %v to %v: %v", args, *pamFlag, err)
		}
	} else {
		if len(args) != 2 {
			flag.Usage()
			os.Exit(1)
		}
		sort(args[0], args[1])
	}
}
