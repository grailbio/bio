// Copyright 2018 GRAIL, Inc. All rights reserved.
// Use of this source code is governed by the Apache 2.0
// license that can be found in the LICENSE file.

package pam_test

import (
	"flag"
	"fmt"
	"io"
	"io/ioutil"
	"math"
	"os"
	"path/filepath"
	"runtime"
	"strings"
	"sync/atomic"
	"testing"
	"time"

	"github.com/grailbio/base/grail"
	"github.com/grailbio/base/traverse"
	"github.com/grailbio/base/vcontext"
	"github.com/grailbio/bio/biopb"
	gbam "github.com/grailbio/bio/encoding/bam"
	"github.com/grailbio/bio/encoding/bamprovider"
	"github.com/grailbio/bio/encoding/converter"
	"github.com/grailbio/bio/encoding/pam"
	"github.com/grailbio/bio/encoding/pam/pamutil"
	"github.com/grailbio/hts/bam"
	"github.com/grailbio/hts/sam"
	"github.com/grailbio/testutil"
	"github.com/grailbio/testutil/assert"
	"github.com/grailbio/testutil/expect"
	"v.io/x/lib/vlog"
)

func mustOpenBAM(t testing.TB, bamPath string) *bam.Reader {
	in, err := os.Open(bamPath)
	assert.NoError(t, err)
	r, err := bam.NewReader(in, runtime.NumCPU())
	assert.NoError(t, err)
	// Note: file descriptor for "in" leaks here.
	return r
}

func generatePAM(t testing.TB, opts pam.WriteOpts, pamPath, bamPath string) {
	assert.NoError(t, pamutil.ValidateCoordRange(&opts.Range))
	rbam := mustOpenBAM(t, bamPath)
	w := pam.NewWriter(opts, rbam.Header(), pamPath)
	n := 0
	for {
		rec, err := rbam.Read()
		if err != nil {
			assert.EQ(t, io.EOF, err)
			break
		}
		vlog.VI(1).Infof("Org: %v", rec)
		recAddr := gbam.CoordFromSAMRecord(rec, 0)
		if !opts.Range.Contains(recAddr) {
			continue
		}
		w.Write(rec)
		assert.NoError(t, w.Err())
		sam.PutInFreePool(rec)
		n++
	}
	assert.NoError(t, w.Close())
	vlog.Infof("Converted %v -> %v (%+v), %d records", bamPath, pamPath, opts, n)
}

func newPAMPath(bamPath string, tempDir string) string {
	return filepath.Join(tempDir, filepath.Base(bamPath))
}

func verifyPAM(t *testing.T, opts pam.ReadOpts, pamPath, bamPath string) {
	verifyPAMWithShardedReader(t, opts, pamPath, bamPath,
		[]biopb.CoordRange{gbam.UniversalRange})
}

func verifyPAMWithShardedReader(t *testing.T, opts pam.ReadOpts, pamPath, bamPath string, shards []biopb.CoordRange) {
	assert.NoError(t, pamutil.ValidateCoordRange(&opts.Range))
	in, err := os.Open(bamPath)
	assert.NoError(t, err)
	defer in.Close()
	rbam, err := bam.NewReader(in, 1)
	assert.NoError(t, err)
	bamAddr := gbam.NewCoordGenerator()
	readBAM := func() *sam.Record {
		for {
			rec, err := rbam.Read()
			if err != nil {
				if err == io.EOF {
					return nil
				}
				t.Fatal(err)
			}
			recAddr := bamAddr.GenerateFromRecord(rec)
			if opts.Range.Contains(recAddr) {
				return rec
			}
			vlog.VI(1).Infof("Skip %v %v %v", recAddr, opts.Range, rec)
			sam.PutInFreePool(rec)
			continue
		}
	}

	vlog.VI(1).Infof("Comparing %v and %v with %d shards %+v", pamPath, bamPath, len(shards), shards)
	n := 0
	for _, bound := range shards {
		vlog.VI(1).Infof("Start reading shard %+v, n=%d", bound, n)
		localOpts := opts
		localOpts.Range = bound
		rpam := pam.NewReader(localOpts, pamPath)
		for rpam.Scan() {
			recPAM := rpam.Record()
			recBAM := readBAM()
			assert.NotNil(t, recBAM, "%d: missing BAM record for %v, with opts %+v", n, recPAM, localOpts)
			assert.EQ(t, recPAM.String(), recBAM.String())
			sam.PutInFreePool(recPAM)
			sam.PutInFreePool(recBAM)
			n++
		}
		assert.NoError(t, rpam.Close())
	}
	if rec := readBAM(); rec != nil {
		t.Fatalf("%d: Excess record in BAM: %v", n, rec)
	}
}

func TestReadWriteMultipleBlocks(t *testing.T) {
	tempDir, cleanup := testutil.TempDir(t, "", "")
	defer cleanup()
	bamPath := testutil.GetFilePath("//go/src/grail.com/bio/encoding/bam/testdata/test.bam")
	pamPath := filepath.Join(tempDir, "test")

	assert.NoError(t, converter.ConvertToPAM(pam.WriteOpts{MaxBufSize: 150}, pamPath, bamPath, "", math.MaxInt64))
	verifyPAM(t, pam.ReadOpts{}, pamPath, bamPath)
}

func TestWriteEmptyFile(t *testing.T) {
	tempDir, cleanup := testutil.TempDir(t, "", "")
	defer cleanup()
	rbam := mustOpenBAM(t, testutil.GetFilePath("//go/src/grail.com/bio/encoding/bam/testdata/test.bam"))
	pamPath := filepath.Join(tempDir, "test")
	w := pam.NewWriter(pam.WriteOpts{}, rbam.Header(), pamPath)
	assert.NoError(t, w.Close())
	r := pam.NewReader(pam.ReadOpts{}, pamPath)
	assert.False(t, r.Scan(), "There should be no record")
	assert.NoError(t, r.Close())
}

func TestNewWriterError(t *testing.T) {
	rbam := mustOpenBAM(t, testutil.GetFilePath("//go/src/grail.com/bio/encoding/bam/testdata/test.bam"))
	rec, err := rbam.Read()
	assert.NoError(t, err)
	assert.NoError(t, rbam.Close())

	w := pam.NewWriter(pam.WriteOpts{}, rbam.Header(), "/non/existing")
	w.Write(rec)
	err = w.Close()
	assert.NotNil(t, err)
	assert.Regexp(t, err.Error(), "no such file or directory")
}

func TestNewReaderError0(t *testing.T) {
	r := pam.NewReader(pam.ReadOpts{}, "/non/existing")
	assert.False(t, r.Scan(), "No record is expected")
	err := r.Close()
	assert.NotNil(t, err)
	assert.Regexp(t, err, ".*no pam file found.*")
}

func TestReadSubsetColumns(t *testing.T) {
	tempDir, cleanup := testutil.TempDir(t, "", "")
	defer cleanup()
	bamPath := testutil.GetFilePath("//go/src/grail.com/bio/encoding/bam/testdata/test.bam")
	readModel := func() []*sam.Record {
		model := []*sam.Record{}
		rbam := mustOpenBAM(t, bamPath)
		for {
			rec, err := rbam.Read()
			if err == io.EOF {
				break
			}
			assert.NoError(t, err)
			model = append(model, rec)
		}
		return model
	}

	pamPath := newPAMPath(bamPath, tempDir)
	assert.NoError(t, converter.ConvertToPAM(pam.WriteOpts{}, pamPath, bamPath, "", math.MaxInt64))
	r := pam.NewReader(pam.ReadOpts{DropFields: []gbam.FieldType{gbam.FieldQual, gbam.FieldName}}, pamPath)
	n := 0
	model := readModel()
	for r.Scan() {
		rec := r.Record()
		m := model[n]
		m.Qual = pam.GetDummyQual(rec.Seq.Length)
		m.Name = ""
		assert.EQ(t, m.String(), rec.String())
		n++
	}
	assert.EQ(t, len(model), n)
	assert.NoError(t, r.Close())

	r = pam.NewReader(pam.ReadOpts{DropFields: []gbam.FieldType{gbam.FieldQual, gbam.FieldSeq}}, pamPath)
	model = readModel()
	n = 0
	for r.Scan() {
		rec := r.Record()
		m := model[n]
		m.Qual = nil
		m.Seq.Length = 0
		m.Seq.Seq = nil
		assert.EQ(t, m.String(), rec.String())
		n++
	}
	assert.EQ(t, len(model), n)
	assert.NoError(t, r.Close())

	r = pam.NewReader(pam.ReadOpts{DropFields: []gbam.FieldType{gbam.FieldSeq, gbam.FieldAux}}, pamPath)
	model = readModel()
	n = 0
	for r.Scan() {
		rec := r.Record()
		m := model[n]
		m.Seq = pam.GetDummySeq(len(m.Qual))
		m.AuxFields = nil
		assert.EQ(t, m.String(), rec.String())
		n++
	}
	assert.EQ(t, len(model), n)
	assert.NoError(t, r.Close())
}

func TestReadWriteUnmapped(t *testing.T) {
	tempDir, cleanup := testutil.TempDir(t, "", "")
	defer cleanup()

	bamPath := testutil.GetFilePath("//go/src/grail.com/bio/encoding/bam/testdata/test-unmapped.bam")
	pamPath := newPAMPath(bamPath, tempDir)
	assert.NoError(t, converter.ConvertToPAM(pam.WriteOpts{}, pamPath, bamPath, "", math.MaxInt64))
	verifyPAM(t, pam.ReadOpts{}, pamPath, bamPath)
	// Test subrange reads.
	testCases := []struct {
		opts     pam.ReadOpts
		expected []string // Expected names of records.
	}{
		// Exclude unmapped segments.
		{
			opts: pam.ReadOpts{Range: biopb.CoordRange{
				Start: biopb.Coord{0, 0, 0},
				Limit: biopb.Coord{biopb.LimitValidRefID, biopb.InfinityPos, 0}}},
			expected: []string{"read1", "read2", "read3"},
		},
		// Read only unmapped segments.
		{
			opts: pam.ReadOpts{Range: biopb.CoordRange{
				Start: biopb.Coord{biopb.UnmappedRefID, 0, 0},
				Limit: biopb.Coord{biopb.UnmappedRefID, biopb.InfinityPos, 0}}},
			expected: []string{"read10", "read10"},
		},
	}
	for _, tc := range testCases {
		vlog.VI(1).Infof("Start test %+v", tc)
		r := pam.NewReader(tc.opts, pamPath)
		for _, name := range tc.expected {
			assert.True(t, r.Scan(), tc)
			rec := r.Record()
			assert.EQ(t, name, rec.Name, tc)
		}
		assert.False(t, r.Scan(), "extra rec", tc)
		assert.NoError(t, r.Close(), tc)
	}
	// Do GC frequently and make sure we haven't screwed up unsafe arena
	// management.
	{
		r := pam.NewReader(pam.ReadOpts{}, pamPath)
		for r.Scan() {
			s0 := r.Record().String()
			runtime.GC()
			s1 := r.Record().String()
			assert.EQ(t, s1, s0)
		}
		assert.NoError(t, r.Close())
	}
}

func TestReadWriteLarge(t *testing.T) {
	tempDir, cleanup := testutil.TempDir(t, "", "")
	defer cleanup()

	bamPath := testutil.GetFilePath("//go/src/grail.com/bio/encoding/bam/testdata/170614_WGS_LOD_Pre_Library_B3_27961B_05.merged.10000.bam")
	pamPath := newPAMPath(bamPath, tempDir)

	assert.NoError(t, converter.ConvertToPAM(pam.WriteOpts{}, pamPath, bamPath, "", math.MaxInt64))
	verifyPAM(t, pam.ReadOpts{}, pamPath, bamPath)
}

func mustGenerateReadShards(t *testing.T, opts pamutil.GenerateReadShardsOpts, pamPath string) []biopb.CoordRange {
	shards, err := pamutil.GenerateReadShards(vcontext.Background(), opts, pamPath, gbam.FieldNames)
	assert.NoError(t, err)
	return shards
}

func newRange(ref0, pos0, ref1, pos1 int) biopb.CoordRange {
	return biopb.CoordRange{
		biopb.Coord{int32(ref0), int32(pos0), 0},
		biopb.Coord{int32(ref1), int32(pos1), 0}}
}

func TestSharder0(t *testing.T) {
	tempDir, cleanup := testutil.TempDir(t, "", "")
	defer cleanup()
	bamPath := testutil.GetFilePath("//go/src/grail.com/bio/encoding/bam/testdata/test-unmapped.bam")
	pamPath := newPAMPath(bamPath, tempDir)
	assert.NoError(t, converter.ConvertToPAM(pam.WriteOpts{}, pamPath, bamPath, "", math.MaxInt64))

	ranges := mustGenerateReadShards(t, pamutil.GenerateReadShardsOpts{NumShards: 1}, pamPath)
	expect.EQ(t, boundString(ranges), "0:0,-:-")
	ranges = mustGenerateReadShards(t, pamutil.GenerateReadShardsOpts{
		NumShards:                          1,
		AlwaysSplitMappedAndUnmappedCoords: true,
	}, pamPath)
	expect.EQ(t, boundString(ranges), "0:0,-:0 -:0,-:-")
	ranges = mustGenerateReadShards(t, pamutil.GenerateReadShardsOpts{
		Range:     newRange(1, 2, 2, 100),
		NumShards: 1,
	}, pamPath)
	expect.EQ(t, boundString(ranges), "1:2,2:100")
	ranges = mustGenerateReadShards(t, pamutil.GenerateReadShardsOpts{
		Range:                              newRange(1, 2, 2, 100),
		NumShards:                          1,
		AlwaysSplitMappedAndUnmappedCoords: true,
	}, pamPath)
	expect.EQ(t, boundString(ranges), "1:2,2:100")
}

func boundString(bounds []biopb.CoordRange) string {
	s := make([]string, len(bounds))
	for i, b := range bounds {
		s[i] = pamutil.CoordRangePathString(b)
	}
	return strings.Join(s, " ")
}

func TestConvert(t *testing.T) {
	tempDir, cleanup := testutil.TempDir(t, "", "")
	defer cleanup()
	ctx := vcontext.Background()
	bamPath := testutil.GetFilePath("//go/src/grail.com/bio/encoding/bam/testdata/170614_WGS_LOD_Pre_Library_B3_27961B_05.merged.10000.bam")
	pamPath := newPAMPath(bamPath, tempDir)

	// The bam file is 2.8MB, so with 1MB shard size, we expect three PAM
	// shard files.
	assert.NoError(t, converter.ConvertToPAM(pam.WriteOpts{}, pamPath, bamPath, "", 1<<20))
	verifyPAM(t, pam.ReadOpts{}, pamPath, bamPath)
	indexes, err := pamutil.ListIndexes(ctx, pamPath)
	assert.NoError(t, err)
	assert.EQ(t, 3, len(indexes), "Index:", indexes)

	// Try 256KB shard size.
	assert.NoError(t, converter.ConvertToPAM(pam.WriteOpts{}, pamPath, bamPath, "", 1<<18))
	verifyPAM(t, pam.ReadOpts{}, pamPath, bamPath)
	indexes, err = pamutil.ListIndexes(ctx, pamPath)
	assert.NoError(t, err)
	assert.EQ(t, 11, len(indexes), "Index:", indexes)
}

func TestSharder1(t *testing.T) {
	tempDir, cleanup := testutil.TempDir(t, "", "")
	defer cleanup()

	// Create three PAM rowshards.
	bamPath := testutil.GetFilePath("//go/src/grail.com/bio/encoding/bam/testdata/170614_WGS_LOD_Pre_Library_B3_27961B_05.merged.10000.bam")
	pamPath := newPAMPath(bamPath, tempDir)

	for _, shardRange := range []biopb.CoordRange{
		biopb.CoordRange{biopb.Coord{0, 0, 0}, biopb.Coord{1, 0, 0}},
		biopb.CoordRange{biopb.Coord{1, 0, 0}, biopb.Coord{3, 0, 0}},
		biopb.CoordRange{biopb.Coord{3, 0, 0}, biopb.Coord{biopb.InfinityRefID, biopb.InfinityPos, 0}},
	} {
		generatePAM(t, pam.WriteOpts{Range: shardRange}, pamPath, bamPath)
	}
	verifyPAM(t, pam.ReadOpts{}, pamPath, bamPath)

	// Try creating just one shard. There will be one shard for each file.
	shards := mustGenerateReadShards(t, pamutil.GenerateReadShardsOpts{NumShards: 1}, pamPath)
	assert.EQ(t, len(shards), 3)
	verifyPAMWithShardedReader(t, pam.ReadOpts{}, pamPath, bamPath, shards)

	// The same test, but specify the params via BytesPerShard.
	shards = mustGenerateReadShards(t, pamutil.GenerateReadShardsOpts{BytesPerShard: math.MaxInt64}, pamPath)
	assert.EQ(t, len(shards), 3)
	verifyPAMWithShardedReader(t, pam.ReadOpts{}, pamPath, bamPath, shards)
}

type syntheticTester struct {
	t               *testing.T
	tmpDir          string
	header          *sam.Header
	seq             int
	cleanupCallback func()

	// PAM generated in generatePAM.
	pamPath string
	// Records generated in generatePAM.
	recs []*sam.Record
}

func newSyntheticTester(t *testing.T) *syntheticTester {
	in := mustOpenBAM(t, testutil.GetFilePath("//go/src/grail.com/bio/encoding/bam/testdata/test.bam"))
	tmpDir, cleanup := testutil.TempDir(t, "", "")
	return &syntheticTester{
		t:               t,
		tmpDir:          tmpDir,
		header:          in.Header(),
		cleanupCallback: cleanup,
	}
}

func (st *syntheticTester) cleanup() {
	st.cleanupCallback()
}

// Create a new PAM file. Returns the pam path prefix.
func (st *syntheticTester) generatePAM(opts pam.WriteOpts, nRecords int,
	posCallback func(index int) (*sam.Reference, int)) string {
	path := filepath.Join(st.tmpDir, fmt.Sprintf("test%03d.pam", st.seq))
	st.seq++
	w := pam.NewWriter(opts, st.header, path)

	seq := []byte{}
	qual := []byte{}
	for i := 0; i < 100; i++ {
		seq = append(seq, 'C')
		qual = append(qual, 50)
	}
	cigar := sam.Cigar{sam.NewCigarOp(sam.CigarMatch, len(seq))}
	st.recs = nil
	st.pamPath = path
	for i := 0; i < nRecords; i++ {
		ref, pos := posCallback(i)
		rec, err := sam.NewRecord(fmt.Sprintf("seq%06d", i), ref, st.header.Refs()[1],
			pos /*pos*/, pos+100 /*matepos*/, 10 /*templen*/, 60 /*mapq*/, cigar, seq, qual, nil)
		assert.NoError(st.t, err, "ref=%v, pos=%d", ref, pos)
		w.Write(rec)
		st.recs = append(st.recs, rec)
	}
	assert.NoError(st.t, w.Close())
	return path
}

// Read PAM generated by the last call to generatePAM, using sharding strategy
// "shards". Verify that records read match those produced by generatePAM.
func (st *syntheticTester) verifyPAM(shards []biopb.CoordRange) {
	vlog.VI(1).Infof("Start verify using shards %+v", shards)
	n := 0
	for _, bound := range shards {
		vlog.VI(1).Infof("Start reading shard %+v, n=%d", bound, n)
		opts := pam.ReadOpts{Range: bound}
		r := pam.NewReader(opts, st.pamPath)
		for r.Scan() {
			rec := r.Record()
			assert.True(st.t, n < len(st.recs), "n=%d, len=%d", n, len(st.recs))
			assert.EQ(st.t, rec.String(), st.recs[n].String(), "n=%d", n)
			n++
		}
		assert.NoError(st.t, r.Close())
	}
	assert.EQ(st.t, n, len(st.recs))
}

// Create a synthetic PAM file where reads have unique coordinates.
func TestSyntheticUniquePositions(t *testing.T) {
	st := newSyntheticTester(t)
	defer st.cleanup()

	const nRecords = 10000
	writeOpts := pam.WriteOpts{MaxBufSize: 1024}
	ref := st.header.Refs()[0]
	tmpPAMPath := st.generatePAM(writeOpts, nRecords, func(index int) (*sam.Reference, int) { return ref, index })
	shards := mustGenerateReadShards(t, pamutil.GenerateReadShardsOpts{NumShards: 16}, tmpPAMPath)
	expect.EQ(t, len(shards), 16, "Shards: %+v", shards)
	for _, shard := range shards {
		assert.True(t, shard.Start.Seq == 0 && shard.Limit.Seq == 0, "Shard: %+v", shard)
	}
	st.verifyPAM(shards)
}

// Create a synthetic PAM file where all reads are at coordinate (0,0)
func TestSyntheticAllReadsAtZero(t *testing.T) {
	st := newSyntheticTester(t)
	defer st.cleanup()

	const nRecords = 10000
	writeOpts := pam.WriteOpts{MaxBufSize: 1024}
	tmpPAMPath := st.generatePAM(writeOpts, nRecords, func(index int) (*sam.Reference, int) { return st.header.Refs()[0], 0 })
	// With the default sharder, there will be only one shard.
	shards := mustGenerateReadShards(t, pamutil.GenerateReadShardsOpts{NumShards: 16}, tmpPAMPath)
	expect.EQ(t, len(shards), 1, "Shards: %+v", shards)
	assert.True(t, shards[0].Start.Seq == 0 && shards[0].Limit.Seq == 0, "Shards: %+v", shards)
	st.verifyPAM(shards)

	// Allow splitting positions.
	const numShards = 16
	shards = mustGenerateReadShards(t, pamutil.GenerateReadShardsOpts{
		SplitMappedCoords:   true,
		SplitUnmappedCoords: true,
		NumShards:           numShards}, tmpPAMPath)
	expect.EQ(t, len(shards), numShards)
	expectedRecordsPerShard := nRecords / numShards
	for i, shard := range shards {
		if shard.Limit.RefId == 0 {
			expect.EQ(t, 0, int(shard.Start.RefId), "shard=%v i=%v shards=%v", shard, i, shards)
			expect.EQ(t, 0, int(shard.Limit.Pos), "shard=%v i=%v shards=%v", shard, i, shards)
			expect.True(t, shard.Limit.Seq > 0, "shard=%v i=%v shards=%v", shard, i, shards)
			nRecords := float64(shard.Limit.Seq - shard.Start.Seq)
			expect.True(t,
				nRecords > float64(expectedRecordsPerShard)*0.8 && nRecords < float64(expectedRecordsPerShard)*1.2,
				"shard=%v nrecords=%v, expected=%v", shard, nRecords, expectedRecordsPerShard)

		}
	}
	st.verifyPAM(shards)
}

// Create a synthetic PAM file where half the records are at (0,0), other half are unmapped.
func TestSyntheticHalfUnmapped(t *testing.T) {
	st := newSyntheticTester(t)
	defer st.cleanup()

	const nRecords = 10000
	writeOpts := pam.WriteOpts{MaxBufSize: 1024}
	tmpPAMPath := st.generatePAM(writeOpts, nRecords, func(index int) (*sam.Reference, int) {
		if index < nRecords/2 {
			return st.header.Refs()[0], 0 // mapped
		}
		return nil, -1 // unmapped
	})
	shards := mustGenerateReadShards(t, pamutil.GenerateReadShardsOpts{
		SplitMappedCoords:   false,
		SplitUnmappedCoords: true,
		NumShards:           16}, tmpPAMPath)
	expect.True(t, len(shards) > 1 && len(shards) < 16, shards)
	for i, shard := range shards {
		expect.EQ(t,

			shard.Start.RefId == biopb.UnmappedRefID, shard.Start.Seq != 0,

			"Shard %d: %+v", i, shard)

	}
	st.verifyPAM(shards)

	shards = mustGenerateReadShards(t, pamutil.GenerateReadShardsOpts{
		SplitMappedCoords:   true,
		SplitUnmappedCoords: false,
		NumShards:           16}, tmpPAMPath)
	expect.True(t, len(shards) > 1 && len(shards) < 16, shards)
	t.Logf("Found %v shards, %+v", len(shards), shards)
	for i, shard := range shards {
		expect.EQ(t,

			shard.Limit.RefId != biopb.UnmappedRefID, shard.Limit.Seq != 0,

			"Shard %d: %+v", i, shard)

	}
	st.verifyPAM(shards)
}

func TestShardedUnmappedReads(t *testing.T) {
	tempDir, cleanup := testutil.TempDir(t, "", "")
	defer cleanup()
	bamPath := testutil.GetFilePath("//go/src/grail.com/bio/encoding/bam/testdata/170614_WGS_LOD_Pre_Library_B3_27961B_05.merged.10000.bam")
	pamPath := filepath.Join(tempDir, "test")
	generatePAM(t, pam.WriteOpts{MaxBufSize: 10000}, pamPath, bamPath)

	// shards0 won't split unmapped sequences.
	shards0 := mustGenerateReadShards(t, pamutil.GenerateReadShardsOpts{NumShards: 2000}, pamPath)
	// shards1 will split unmapped sequences.
	shards1 := mustGenerateReadShards(t, pamutil.GenerateReadShardsOpts{SplitUnmappedCoords: true, NumShards: 2000}, pamPath)

	// shards1 contain the same set of shards for mapped reads, but it also
	// splits unmapped reads into multiple shards.
	if len(shards1) <= len(shards0) {
		t.Errorf("Wrong shard sizes: %v %v", len(shards1), len(shards0))
	}
	for _, s := range shards0 {
		if s.Start.Seq != 0 || s.Limit.Seq != 0 {
			t.Error(s)
		}
	}
	// The last few shards of shards1 must be for unmapped reads.
	nOpen := 0
	nClosed := 0
	for _, s := range shards1 {
		if s.Start.Seq != 0 || s.Limit.Seq != 0 {
			nOpen++
		} else {
			nClosed++
		}
	}
	if nOpen < 2 || nClosed < 2 {
		t.Fatal(nOpen, nClosed, shards1)
	}
	verifyPAMWithShardedReader(t, pam.ReadOpts{}, pamPath, bamPath, shards0)
	verifyPAMWithShardedReader(t, pam.ReadOpts{}, pamPath, bamPath, shards1)
}

var (
	bamFlag        = flag.String("bam", "/scratch-nvme/bam/CNVS-NORM-110033752-cfDNA-WGBS-Rep1.bam", "File to generate in benchmark")
	pamFlag        = flag.String("pam", "", "PAM file to produce. If empty, write in a temp dir.")
	dropFieldsFlag = flag.String("drop-fields", "", "Comma-separated fields to drop during PAM benchmarks")
	unmappedFlag   = flag.Bool("unmapped", true, "If true, read unmapped sequences as well as mapped ones")
	tmpdirFlag     = flag.String("tmpdir", "", "Temp dir used in benchmarks")
)

func BenchmarkConvert(b *testing.B) {
	pamPath := *pamFlag
	if pamPath == "" {
		tempDir, cleanup := testutil.TempDir(b, "", "")
		defer cleanup()
		pamPath = filepath.Join(tempDir, "bench.pam")
	}
	for n := 0; n < b.N; n++ {
		assert.NoError(b, converter.ConvertToPAM(pam.WriteOpts{}, pamPath, *bamFlag, "", math.MaxInt64))
	}
}

func BenchmarkReadPAM(b *testing.B) {
	b.StopTimer()
	pamPath := *pamFlag
	if pamPath == "" {
		if *bamFlag == "" {
			vlog.Fatal("No input specified")
		}
		tempDir, cleanup := testutil.TempDir(b, *tmpdirFlag, "")
		defer cleanup()
		pamPath = filepath.Join(tempDir, *bamFlag)
	}
	if files, err := filepath.Glob(pamPath + "*.index"); err != nil || len(files) == 0 {
		vlog.Infof("Generating PAM files %v", pamPath)
		assert.NoError(b, converter.ConvertToPAM(pam.WriteOpts{}, pamPath, *bamFlag, "", math.MaxInt64))
	} else {
		vlog.Infof("Reusing PAM files %v", pamPath)
	}
	b.StartTimer()
	for n := 0; n < b.N; n++ {
		opts := pam.ReadOpts{}
		if *dropFieldsFlag != "" {
			for _, fieldName := range strings.Split(*dropFieldsFlag, ",") {
				f, err := gbam.ParseFieldType(fieldName)
				assert.NoError(b, err)
				opts.DropFields = append(opts.DropFields, f)
			}
		}
		if !*unmappedFlag {
			vlog.Infof("Skipping unmapped reads")
			opts.Range = gbam.MappedRange
		}
		ctx := vcontext.Background()
		bounds, err := pamutil.GenerateReadShards(ctx, pamutil.GenerateReadShardsOpts{Range: opts.Range}, pamPath, gbam.FieldNames)
		assert.NoError(b, err)
		boundCh := make(chan biopb.CoordRange, len(bounds))
		for _, r := range bounds {
			boundCh <- r
		}
		close(boundCh)
		totalRecs := int64(0)
		traverse.CPU(func() error { // nolint: errcheck
			for bound := range boundCh {
				start := time.Now()
				vlog.Infof("Start read %+v", bound)
				r := pam.NewReader(pam.ReadOpts{DropFields: opts.DropFields, Range: bound}, pamPath)
				nRecs := int64(0)
				for r.Scan() {
					rec := r.Record()
					sam.PutInFreePool(rec)
					nRecs++
				}
				assert.NoError(b, r.Close())
				atomic.AddInt64(&totalRecs, nRecs)
				end := time.Now()
				vlog.Infof("Finish read %+v, %d recs %d ms", bound, nRecs, end.Sub(start)/time.Millisecond)
			}
			return nil
		})
		vlog.Infof("%v: read %d records", pamPath, totalRecs)
	}
}

func BenchmarkReadBAM(b *testing.B) {
	if *bamFlag == "" {
		vlog.Fatal("No input specified")
	}
	for n := 0; n < b.N; n++ {
		provider := bamprovider.NewProvider(*bamFlag)
		header, err := provider.GetHeader()
		assert.NoError(b, err)
		shardList, err := gbam.GetPositionBasedShards(header, 100000, 0, *unmappedFlag)
		assert.NoError(b, err)
		shardCh := gbam.NewShardChannel(shardList)
		assert.NoError(b, err)
		parallelism := runtime.NumCPU()
		totalRecs := int64(0)

		traverse.Each(parallelism, func(_ int) error { // nolint: errcheck
			for {
				shard, ok := <-shardCh
				if !ok {
					break
				}
				nRecs := int64(0)
				iter := provider.NewIterator(shard)
				for iter.Scan() {
					record := iter.Record()
					nRecs++
					sam.PutInFreePool(record)
				}
				atomic.AddInt64(&totalRecs, nRecs)
				assert.NoError(b, iter.Close())
			}
			return nil
		})
		assert.NoError(b, provider.Close())
		vlog.Infof("%v: read %d records", *bamFlag, totalRecs)
	}
}

func mustCreate(t *testing.T, path string) {
	fd, err := os.Create(path)
	assert.NoError(t, err)
	fd.WriteString("foo")
	assert.NoError(t, fd.Close())
}

func TestRemove(t *testing.T) {
	tempDir, cleanup := testutil.TempDir(t, "", "")
	defer cleanup()

	// Removing non-existing file should get no error.
	assert.NoError(t, pamutil.Remove(filepath.Join(tempDir, "gg")))

	path := filepath.Join(tempDir, "f")
	assert.NoError(t, os.MkdirAll(path, 0777))
	mustCreate(t, filepath.Join(path, "0:0,2:123.index"))
	mustCreate(t, filepath.Join(path, "0:0,2:123.aux"))

	_, err := os.Stat(path)
	assert.NoError(t, err)
	assert.NoError(t, pamutil.Remove(path))
	_, err = os.Stat(path)
	assert.True(t, os.IsNotExist(err))
}

func TestListIndexes(t *testing.T) {
	tempDir, cleanup := testutil.TempDir(t, "", "")
	defer cleanup()
	path := filepath.Join(tempDir, "f")
	assert.NoError(t, os.MkdirAll(path, 0777))
	mustCreate(t, path+"/0:0,2:123.index")
	mustCreate(t, path+"/2:123,10:200.index")
	mustCreate(t, path+"/10:200,-:-.index")

	ctx := vcontext.Background()
	indexes, err := pamutil.ListIndexes(ctx, path)
	assert.NoError(t, err)
	expected := []biopb.CoordRange{
		biopb.CoordRange{biopb.Coord{0, 0, 0}, biopb.Coord{2, 123, 0}},
		biopb.CoordRange{biopb.Coord{2, 123, 0}, biopb.Coord{10, 200, 0}},
		biopb.CoordRange{biopb.Coord{10, 200, 0}, biopb.Coord{biopb.InfinityRefID, biopb.InfinityPos, 0}}}
	expect.EQ(t, len(expected), len(indexes))
	for i, e := range expected {
		expect.EQ(t, indexes[i].Range, e)
	}

	indexes, err = pamutil.ListIndexes(ctx, path+"/")
	assert.NoError(t, err)
	expect.EQ(t, len(expected), len(indexes))
	for i, e := range expected {
		expect.EQ(t, indexes[i].Range, e)
	}
}

func TestMain(m *testing.M) {
	shutdown := grail.Init()
	if *tmpdirFlag == "" {
		tryDir := func(path string) bool {
			tempDir, err := ioutil.TempDir(path, "pam")
			if err != nil {
				return false
			}
			vlog.Infof("Using tempdir %v", tempDir)
			os.RemoveAll(tempDir)
			return true
		}
		if !tryDir("/scratch-nvme/tmp") {
			*tmpdirFlag = "/scratch-nvme/tmp"
		}
	}
	status := m.Run()
	shutdown()
	os.Exit(status)
}
