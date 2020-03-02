package bamprovider_test

import (
	"compress/gzip"
	"flag"
	"os"
	"path/filepath"
	"runtime"
	"strings"
	"sync"
	"testing"

	gbam "github.com/grailbio/bio/encoding/bam"
	"github.com/grailbio/bio/encoding/bamprovider"
	"github.com/grailbio/hts/bam"
	"github.com/grailbio/hts/sam"
	"github.com/grailbio/testutil"
	"github.com/grailbio/testutil/assert"
	"github.com/grailbio/testutil/expect"
	"v.io/x/lib/vlog"
)

var (
	// Flags for BenchmarkWrite
	inFile = flag.String("in",
		"//go/src/grail.com/bio/encoding/bam/testdata/170614_WGS_LOD_Pre_Library_B3_27961B_05.merged.10000.bam",
		"Input BAM filename. If the path starts with '//', it is assumed relative to the relative of the workspace")
	outFile             = flag.String("out", "", "Output BAM filename. IF empty, writes to a temporary file")
	useShardedBAMWriter = flag.Bool("useshardedbamwriter", false, "use ShardedBAMWriter")
	shardSize           = flag.Int("shard-size", 1000000, "shard size")
	gzLevel             = flag.Int("gz-level", gzip.DefaultCompression, "gz compression level")
	parallelism         = flag.Int("parallelism", 2*runtime.NumCPU(), "parallelism")
	queueLength         = flag.Int("queue-length", 4*runtime.NumCPU(), "queue length")
)

func processShards(b *testing.B, provider bamprovider.Provider, worker int, channel chan gbam.Shard,
	shardedwriter *gbam.ShardedBAMWriter, biogoout chan []*sam.Record) {

	var compressor *gbam.ShardedBAMCompressor
	if *useShardedBAMWriter {
		compressor = shardedwriter.GetCompressor()
	}

	for {
		shard, ok := <-channel
		if !ok {
			vlog.VI(1).Infof("worker %d done", worker)
			break
		}

		iter := provider.NewIterator(shard)
		vlog.VI(1).Infof("starting shard (%s,%d,%d,%d)", shard.StartRef.Name(), shard.Start, shard.End, shard.ShardIdx)
		if *useShardedBAMWriter {
			assert.NoError(b, compressor.StartShard(shard.ShardIdx))
		}

		outlist := make([]*sam.Record, 0)
		for iter.Scan() {
			record := iter.Record()
			if *useShardedBAMWriter {
				compressor.AddRecord(record)
			} else {
				outlist = append(outlist, record)
			}
		}

		if *useShardedBAMWriter {
			err := compressor.CloseShard()
			if err != nil {
				b.Fatalf("Error closing shard %v", err)
			}
		} else {
			biogoout <- outlist
		}
		vlog.VI(1).Infof("finished shard (%s,%d,%d,%d)", shard.StartRef.Name(), shard.Start, shard.End, shard.ShardIdx)
		assert.NoError(b, iter.Close())
	}
}

func biogowriter(b *testing.B, biogoout chan []*sam.Record, bamwriter *bam.Writer) {
	for {
		outlist, ok := <-biogoout
		if !ok {
			break
		}
		for _, r := range outlist {
			err := bamwriter.Write(r)
			if err != nil {
				b.Fatalf("Error writing shard %v", err)
			}
		}
	}
}

func shardedCopy(b *testing.B, inFile, outFile string) {
	// Prepare inputs.
	provider := bamprovider.NewProvider(inFile, bamprovider.ProviderOpts{})
	header, err := provider.GetHeader()
	if err != nil {
		b.Fatalf("Could not read header from file %s: %s", inFile, err)
	}

	// Prepare outputs
	var shardedwriter *gbam.ShardedBAMWriter
	var biogoout chan []*sam.Record
	var outGroup sync.WaitGroup

	out, err := os.Create(outFile)
	if err != nil {
		b.Fatalf("error creating output file %s", outFile)
	}
	if *useShardedBAMWriter {
		// Write the header
		shardedwriter, err = gbam.NewShardedBAMWriter(out, *gzLevel, *queueLength, header)
		if err != nil {
			b.Fatalf("Error initializing ShardedBAM writer: %v", err)
		}
	} else {
		writer, err := bam.NewWriterLevel(out, header, *gzLevel, *parallelism)
		if err != nil {
			b.Fatalf("Error initializing BAM writer: %v", err)
		}

		// start biogo receiver
		biogoout = make(chan []*sam.Record, 100)
		outGroup.Add(1)
		go func() {
			defer outGroup.Done()
			biogowriter(b, biogoout, writer)
			writer.Close()
		}()
	}
	defer out.Close()

	// start workers
	var workerGroup sync.WaitGroup
	shardList, err := gbam.GetPositionBasedShards(header, *shardSize, 0, true)
	expect.Nil(b, err)
	shardChannel := gbam.NewShardChannel(shardList)
	for i := 0; i < *parallelism; i++ {
		vlog.VI(1).Infof("Creating worker %d", i)
		workerGroup.Add(1)
		go func(worker int) {
			defer workerGroup.Done()
			processShards(b, provider, worker, shardChannel, shardedwriter, biogoout)
		}(i)
	}
	workerGroup.Wait()

	if *useShardedBAMWriter {
		if err := shardedwriter.Close(); err != nil {
			b.Fatalf("error in close: %v", err)
		}
	} else {
		close(biogoout)
	}
	outGroup.Wait()
}

// This benchmark allows us to compare the performance of
// ShardedBAMWriter to biogo's writer.  Currently, the biogo output is
// out of order, so it's not a completely fair comparison.
//
// TODO(josh): When this benchmark is located in github.com/grailbio/bio/encoding/bam with the code
// it exercises, the Bazel go_default_test build fails with a package height error that may be
// similar to https://github.com/bazelbuild/rules_go/issues/1877. Consider moving this back to
// package bam when that's fixed.
func BenchmarkWrite(b *testing.B) {
	in := *inFile
	if strings.HasPrefix(in, "//") {
		in = testutil.GetFilePath(in)
	}
	out := *outFile
	if out == "" {
		tmpDir, cleanup := testutil.TempDir(b, "", "")
		defer cleanup()
		out = filepath.Join(tmpDir, "out.bam")
	}
	for i := 0; i < b.N; i++ {
		shardedCopy(b, in, out)
	}
}
