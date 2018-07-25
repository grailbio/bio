package bam_test

import (
	"bytes"
	"flag"
	"io"
	"os"
	"path/filepath"
	"runtime"
	"strings"
	"sync"
	"testing"

	"github.com/biogo/hts/bam"
	"github.com/biogo/hts/sam"
	gbam "github.com/grailbio/bio/encoding/bam"
	"github.com/grailbio/bio/encoding/bamprovider"
	"github.com/grailbio/testutil"
	"github.com/grailbio/testutil/assert"
	"github.com/grailbio/testutil/expect"
	"github.com/klauspost/compress/gzip"
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

func verifyBAM(t *testing.T, records []*sam.Record, bamBuffer *bytes.Buffer) {
	reader, err := bam.NewReader(bamBuffer, 1)
	expect.Nil(t, err)
	i := 0
	for {
		r, err := reader.Read()
		if err == io.EOF {
			expect.EQ(t, i, len(records), "not enough records in bam output %d vs %d", len(records), i)
			break
		}
		expected, err := records[i].MarshalText()
		expect.Nil(t, err)
		actual, err := r.MarshalText()
		expect.Nil(t, err)

		expect.EQ(t, actual, expected, "record[%d] does not match %v vs %v", i, records[i], r)
		i++
	}
}

func writeAndVerify(t *testing.T, header *sam.Header, records []*sam.Record, compressors, shards int, forward bool) {
	var bamBuffer bytes.Buffer
	w, err := gbam.NewShardedBAMWriter(&bamBuffer, gzip.DefaultCompression, 10, header)
	if err != nil {
		t.Errorf("error creating ShardedBAMWriter: %v", err)
	}

	shardCompressors := make([]*gbam.ShardedBAMCompressor, compressors)
	for i := 0; i < compressors; i++ {
		shardCompressors[i] = w.GetCompressor()
	}

	for i := 0; i < shards; i++ {
		var shardNum int
		if forward {
			shardNum = i
		} else {
			shardNum = shards - 1 - i
		}

		c := shardNum % compressors
		err := shardCompressors[c].StartShard(shardNum)
		assert.Nil(t, err)

		for _, r := range records[shardNum*(len(records)/shards) : (shardNum+1)*(len(records)/shards)] {
			err := shardCompressors[c].AddRecord(r)
			expect.Nil(t, err)
		}
		// If there are remainders from uneven division, then add them to the last shard.
		if shardNum == shards-1 {
			for _, r := range records[(shardNum+1)*(len(records)/shards):] {
				err := shardCompressors[c].AddRecord(r)
				expect.Nil(t, err)
			}
		}

		err = shardCompressors[c].CloseShard()
		expect.Nil(t, err)
	}
	err = w.Close()
	expect.Nil(t, err)
	verifyBAM(t, records, &bamBuffer)
}

func TestShardedBAMSmall(t *testing.T) {
	chr1, err := sam.NewReference("chr1", "", "", 1000, nil, nil)
	expect.Nil(t, err)
	chr2, err := sam.NewReference("chr2", "", "", 2000, nil, nil)
	expect.Nil(t, err)
	header, err := sam.NewHeader(nil, []*sam.Reference{chr1, chr2})
	expect.Nil(t, err)
	cigar := []sam.CigarOp{
		sam.NewCigarOp(sam.CigarSoftClipped, 1),
		sam.NewCigarOp(sam.CigarMatch, 8),
		sam.NewCigarOp(sam.CigarSoftClipped, 1),
	}

	r1F := sam.Paired | sam.Read1
	r2R := sam.Paired | sam.Read2 | sam.Reverse

	records := []*sam.Record{
		&sam.Record{Name: "A::::10:1:1", Ref: chr1, Pos: 0, Flags: r1F, MatePos: 10, MateRef: chr1, Cigar: cigar},
		&sam.Record{Name: "A::::10:1:1", Ref: chr1, Pos: 10, Flags: r2R, MatePos: 0, MateRef: chr1, Cigar: cigar},
		&sam.Record{Name: "B::::10:1:1", Ref: chr1, Pos: 20, Flags: r1F, MatePos: 30, MateRef: chr1, Cigar: cigar},
		&sam.Record{Name: "B::::10:1:1", Ref: chr1, Pos: 30, Flags: r2R, MatePos: 20, MateRef: chr1, Cigar: cigar},
		&sam.Record{Name: "C::::10:1:1", Ref: chr1, Pos: 40, Flags: r1F, MatePos: 50, MateRef: chr1, Cigar: cigar},
		&sam.Record{Name: "C::::10:1:1", Ref: chr1, Pos: 50, Flags: r2R, MatePos: 40, MateRef: chr1, Cigar: cigar},
		&sam.Record{Name: "D::::10:1:1", Ref: chr2, Pos: 60, Flags: r1F, MatePos: 70, MateRef: chr2, Cigar: cigar},
		&sam.Record{Name: "D::::10:1:1", Ref: chr2, Pos: 70, Flags: r2R, MatePos: 60, MateRef: chr2, Cigar: cigar},
		&sam.Record{Name: "E::::10:1:1", Ref: chr2, Pos: 80, Flags: r1F, MatePos: 90, MateRef: chr2, Cigar: cigar},
		&sam.Record{Name: "E::::10:1:1", Ref: chr2, Pos: 90, Flags: r2R, MatePos: 80, MateRef: chr2, Cigar: cigar},
		&sam.Record{Name: "F::::10:1:1", Ref: chr2, Pos: 100, Flags: r1F, MatePos: 110, MateRef: chr2, Cigar: cigar},
		&sam.Record{Name: "F::::10:1:1", Ref: chr2, Pos: 110, Flags: r2R, MatePos: 100, MateRef: chr2, Cigar: cigar},
	}
	writeAndVerify(t, header, records, 2, 3, true)
	writeAndVerify(t, header, records, 2, 3, false)
}

func TestShardedBAMLarge(t *testing.T) {
	filename := testutil.GetFilePath("//go/src/grail.com/bio/encoding/bam/testdata/170614_WGS_LOD_Pre_Library_B3_27961B_05.merged.10000.bam")
	f, err := os.Open(filename)
	expect.Nil(t, err)
	reader, err := bam.NewReader(f, 1)
	expect.Nil(t, err)

	records := []*sam.Record{}
	for {
		record, err := reader.Read()
		if err == io.EOF {
			break
		}
		records = append(records, record)
	}

	writeAndVerify(t, reader.Header(), records, 3, 6, true)
	writeAndVerify(t, reader.Header(), records, 3, 6, false)
}

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
