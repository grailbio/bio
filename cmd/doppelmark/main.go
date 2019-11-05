package main

/*
  doppelmark is a tool for marking and removing PCR and optical
  duplicates. For more information, see
  github.com/grailbio/bio/markduplicates/doc.go
*/

import (
	"flag"
	"runtime"
	"strings"

	"github.com/grailbio/base/grail"
	"github.com/grailbio/base/log"
	"github.com/grailbio/base/vcontext"
	gbam "github.com/grailbio/bio/encoding/bam"
	"github.com/grailbio/bio/encoding/bamprovider"
	md "github.com/grailbio/bio/markduplicates"
)

var (
	bamFile              = flag.String("bam", "", "Input BAM filename")
	indexFile            = flag.String("index", "", "Input BAM index filename. By default, set to input BAM filename + .bai")
	outputPath           = flag.String("output", "", "Output filename")
	format               = flag.String("format", "bam", "Output format. Value is either 'bam' or 'pam'.")
	metricsFile          = flag.String("metrics", "", "Output metrics file")
	scratchDir           = flag.String("scratch-dir", "/tmp", "Directory to put scratch files")
	parallelism          = flag.Int("parallelism", runtime.NumCPU(), "Number of parallel computations to run during the markdup phase")
	queueLength          = flag.Int("queue-length", runtime.NumCPU()*5, "Number shards to queue while waiting for flush")
	shardSize            = flag.Int("shard-size", 5000000, "approx shard size in bytes")
	minBases             = flag.Int("min-bases", 5000, "minimum number of bases per shard")
	padding              = flag.Int("clip-padding", 143, "padding in bp, this must be larger than the largest per-read clipping distance")
	clearExisting        = flag.Bool("clear-existing", false, "clear existing duplicate flag before marking")
	removeDups           = flag.Bool("remove-dups", false, "remove duplicates instead of flagging them")
	tagDups              = flag.Bool("tag-duplicates", false, "tag duplicates as DT:Z:SQ (optical) or DT:Z:LB (pcr), and include DI and DS tags")
	useUmis              = flag.Bool("use-umis", false, "use Umi information in read names for grouping duplicates")
	umiFile              = flag.String("umi-file", "", "perform UMI error correction with the known UMIs in this file")
	scavengeUmis         = flag.Int("scavenge-umis", -1, "scavenge UMIs with at most this edit distance")
	separateSingletons   = flag.Bool("separate-singletons", false, "keep singletons separate from pairs, don't bag them together")
	intDI                = flag.Bool("int-di", false, "use integer formatting for DI tags, sets the maximum number of reads to 2147483647 (use for testing only)")
	opticalDistance      = flag.Int("optical-distance", 2500, "pixel distance threshold for optical duplicates, use -1 to disable")
	diskMateShards       = flag.Int("disk-mate-shards", 0, "number of disk shards to use for distant mate storage, use 0 to keep mates in memory.  A value of 1000 is a reasonable choice when using disk, but will require an increase in file descriptor limit, e.g. 'ulimit -n 2000'.")
	emitUnmodifiedFields = flag.Bool("emit-unmodified-fields", false, "Write fields that are not modified. This flag is meaningful only when --format=pam.")
	strandSpecific       = flag.Bool("strand-specific", false, "mark reads only if their r1 strands match")
	opticalHistogram     = flag.String("optical-histogram", "", "path to optical distance histogram output file")
	// The default opticalHistogramMax is set to 2000. Experimentally, the runtimes with 2000 seem reasonable, and it will still consider many duplicate pairs.
	// The histograms looked the same between the full set of duplicate pairs and when capped at 2000.
	opticalHistogramMax = flag.Int("optical-histogram-max", 2000, "maximum number of bag entries to compare when computing optical histogram. Setting to -1 reports for all bag entries.")
)

func main() {
	shutdown := grail.Init()
	defer shutdown()

	// Validate parameters.
	if flag.NArg() > 0 {
		a := flag.Args()
		log.Fatalf("unparsed flags, please check flag syntax: '%s'", strings.Join(a[len(a)-flag.NArg():], " "))
	}

	opts := md.Opts{
		BamFile:              *bamFile,
		IndexFile:            *indexFile,
		MetricsFile:          *metricsFile,
		Format:               *format,
		ShardSize:            *shardSize,
		MinBases:             *minBases,
		Padding:              *padding,
		DiskMateShards:       *diskMateShards,
		ScratchDir:           *scratchDir,
		Parallelism:          *parallelism,
		QueueLength:          *queueLength,
		ClearExisting:        *clearExisting,
		RemoveDups:           *removeDups,
		TagDups:              *tagDups,
		IntDI:                *intDI,
		UseUmis:              *useUmis,
		UmiFile:              *umiFile,
		ScavengeUmis:         *scavengeUmis,
		EmitUnmodifiedFields: *emitUnmodifiedFields,
		SeparateSingletons:   *separateSingletons,
		OutputPath:           *outputPath,
		StrandSpecific:       *strandSpecific,
		OpticalHistogram:     *opticalHistogram,
		OpticalHistogramMax:  *opticalHistogramMax,
	}

	// Create the provider.
	bamOpts := bamprovider.ProviderOpts{Index: opts.IndexFile}
	if !opts.EmitUnmodifiedFields {
		bamOpts.DropFields = []gbam.FieldType{
			gbam.FieldMapq,
			gbam.FieldTempLen,
		}
	}
	provider := bamprovider.NewProvider(*bamFile, bamOpts)

	// Create optical duplicate detector if necessary.
	if *opticalDistance >= 0 {
		opts.OpticalDetector = &md.TileOpticalDetector{
			OpticalDistance: *opticalDistance,
		}
	}

	ctx := vcontext.Background()
	if err := md.SetupAndMark(ctx, provider, &opts); err != nil {
		log.Fatalf(err.Error())
	}
	log.Debug.Printf("exiting")
}
