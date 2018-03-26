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
	"sync"
	"sync/atomic"
	"testing"
	"time"

	"github.com/biogo/hts/bam"
	"github.com/biogo/hts/sam"
	"github.com/grailbio/base/grail"
	"github.com/grailbio/bio/biopb"
	gbam "github.com/grailbio/bio/encoding/bam"
	"github.com/grailbio/bio/encoding/bamprovider"
	"github.com/grailbio/bio/encoding/converter"
	"github.com/grailbio/bio/encoding/pam"
	"github.com/grailbio/internal/testutil"
	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
	"v.io/x/lib/vlog"
)

func mustOpenBAM(t require.TestingT, bamPath string) *bam.Reader {
	in, err := os.Open(bamPath)
	require.NoError(t, err)
	r, err := bam.NewReader(in, runtime.NumCPU())
	require.NoError(t, err)
	// Note: file descriptor for "in" leaks here.
	return r
}

func generatePAM(t require.TestingT, opts pam.WriteOpts, pamPath, bamPath string) {
	require.NoError(t, pam.ValidateCoordRange(&opts.Range))
	rbam := mustOpenBAM(t, bamPath)
	w := pam.NewWriter(opts, rbam.Header(), pamPath)
	n := 0
	for {
		rec, err := rbam.Read()
		if err != nil {
			require.Equal(t, err, io.EOF)
			break
		}
		vlog.VI(1).Infof("Org: %v", rec)
		recAddr := gbam.CoordFromSAMRecord(rec, 0)
		if !opts.Range.Contains(recAddr) {
			continue
		}
		w.Write(rec)
		require.NoError(t, w.Err())
		sam.PutInFreePool(rec)
		n++
	}
	require.NoError(t, w.Close())
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
	require.NoError(t, pam.ValidateCoordRange(&opts.Range))
	in, err := os.Open(bamPath)
	require.NoError(t, err)
	defer in.Close()
	rbam, err := bam.NewReader(in, 1)
	require.NoError(t, err)
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
			if recBAM == nil {
				t.Fatalf("%d: missing BAM record for %v, with opts %+v", n, recPAM, localOpts)
			}
			if recBAM.String() != recPAM.String() {
				t.Fatalf("%d: Records differ:\nbam:\n %v\npam:\n %v",
					n, recBAM, recPAM)
			}
			sam.PutInFreePool(recPAM)
			sam.PutInFreePool(recBAM)
			n++
		}
		require.NoError(t, rpam.Close())
	}
	if rec := readBAM(); rec != nil {
		t.Fatalf("%d: Excess record in BAM: %v", n, rec)
	}
}

func TestReadWriteMultipleBlocks(t *testing.T) {
	tempDir, cleanup := testutil.TempDir(t, "", "")
	defer cleanup()
	bamPath := testutil.GetFilePath("@grailgo//bio/encoding/bam/testdata/test.bam")
	pamPath := filepath.Join(tempDir, "test")

	require.NoError(t, converter.ConvertToPAM(pam.WriteOpts{MaxBufSize: 150}, pamPath, bamPath, "", math.MaxInt64))
	verifyPAM(t, pam.ReadOpts{}, pamPath, bamPath)
}

func TestWriteEmptyFile(t *testing.T) {
	tempDir, cleanup := testutil.TempDir(t, "", "")
	defer cleanup()
	rbam := mustOpenBAM(t, testutil.GetFilePath("@grailgo//bio/encoding/bam/testdata/test.bam"))
	pamPath := filepath.Join(tempDir, "test")
	w := pam.NewWriter(pam.WriteOpts{}, rbam.Header(), pamPath)
	require.NoError(t, w.Close())
	r := pam.NewReader(pam.ReadOpts{}, pamPath)
	require.False(t, r.Scan(), "There should be no record")
	require.NoError(t, r.Close())
}

func TestNewWriterError(t *testing.T) {
	rbam := mustOpenBAM(t, testutil.GetFilePath("@grailgo//bio/encoding/bam/testdata/test.bam"))
	rec, err := rbam.Read()
	require.NoError(t, err)
	require.NoError(t, rbam.Close())

	w := pam.NewWriter(pam.WriteOpts{}, rbam.Header(), "/non/existing")
	w.Write(rec)
	err = w.Close()
	require.Error(t, err)
	require.Regexp(t, "no such file or directory", err.Error())
}

func TestNewReaderError0(t *testing.T) {
	r := pam.NewReader(pam.ReadOpts{}, "/non/existing")
	require.False(t, r.Scan(), "No record is expected")
	err := r.Close()
	require.Error(t, err)
	require.Regexp(t, ".*no index files found.*", err)
}

func TestReadSubsetColumns(t *testing.T) {
	tempDir, cleanup := testutil.TempDir(t, "", "")
	defer cleanup()
	bamPath := testutil.GetFilePath("@grailgo//bio/encoding/bam/testdata/test.bam")
	pamPath := newPAMPath(bamPath, tempDir)

	require.NoError(t, converter.ConvertToPAM(pam.WriteOpts{}, pamPath, bamPath, "", math.MaxInt64))
	r := pam.NewReader(pam.ReadOpts{DropFields: []gbam.FieldType{gbam.FieldQual, gbam.FieldName}}, pamPath)
	n := 0
	for r.Scan() {
		rec := r.Record()
		require.True(t, len(rec.Qual) > 0, "Rec:", rec)
		require.Equal(t, len(rec.Qual), rec.Seq.Length, "Rec:", rec)
		require.Equal(t, rec.Name, "", "Rec:", rec)
		n++
	}
	require.Equal(t, n, 3)
	require.NoError(t, r.Close())

	r = pam.NewReader(pam.ReadOpts{DropFields: []gbam.FieldType{gbam.FieldQual, gbam.FieldSeq}}, pamPath)
	n = 0
	for r.Scan() {
		rec := r.Record()
		require.Equal(t, 0, len(rec.Qual), "Rec:", rec)
		require.Equal(t, 0, rec.Seq.Length, "Rec:", rec)
		require.Equal(t, 0, len(rec.Seq.Seq), "Rec:", rec)
		n++
	}
	require.Equal(t, n, 3)
	require.NoError(t, r.Close())
}

func TestReadWriteUnmapped(t *testing.T) {
	tempDir, cleanup := testutil.TempDir(t, "", "")
	defer cleanup()

	bamPath := testutil.GetFilePath("@grailgo//bio/encoding/bam/testdata/test-unmapped.bam")
	pamPath := newPAMPath(bamPath, tempDir)
	require.NoError(t, converter.ConvertToPAM(pam.WriteOpts{}, pamPath, bamPath, "", math.MaxInt64))
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
			require.True(t, r.Scan(), tc)
			rec := r.Record()
			require.Equal(t, rec.Name, name, tc)
		}
		require.False(t, r.Scan(), "extra rec", tc)
		require.NoError(t, r.Close(), tc)
	}
	// Do GC frequently and make sure we haven't screwed up unsafe arena
	// management.
	{
		r := pam.NewReader(pam.ReadOpts{}, pamPath)
		for r.Scan() {
			s0 := r.Record().String()
			runtime.GC()
			s1 := r.Record().String()
			require.Equal(t, s0, s1)
		}
		require.NoError(t, r.Close())
	}
}

func TestReadWriteLarge(t *testing.T) {
	tempDir, cleanup := testutil.TempDir(t, "", "")
	defer cleanup()

	bamPath := testutil.GetFilePath("@grailgo//bio/encoding/bam/testdata/170614_WGS_LOD_Pre_Library_B3_27961B_05.merged.10000.bam")
	pamPath := newPAMPath(bamPath, tempDir)

	require.NoError(t, converter.ConvertToPAM(pam.WriteOpts{}, pamPath, bamPath, "", math.MaxInt64))
	verifyPAM(t, pam.ReadOpts{}, pamPath, bamPath)
}

func mustGenerateReadShards(t *testing.T, opts pam.GenerateReadShardsOpts, pamPath string) []biopb.CoordRange {
	shards, err := pam.GenerateReadShards(opts, pamPath)
	require.NoError(t, err)
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
	bamPath := testutil.GetFilePath("@grailgo//bio/encoding/bam/testdata/test-unmapped.bam")
	pamPath := newPAMPath(bamPath, tempDir)
	require.NoError(t, converter.ConvertToPAM(pam.WriteOpts{}, pamPath, bamPath, "", math.MaxInt64))

	ranges := mustGenerateReadShards(t, pam.GenerateReadShardsOpts{NumShards: 1}, pamPath)
	assert.Equal(t, "0:0,-:-", boundString(ranges))
	ranges = mustGenerateReadShards(t, pam.GenerateReadShardsOpts{
		NumShards:                          1,
		AlwaysSplitMappedAndUnmappedCoords: true,
	}, pamPath)
	assert.Equal(t, "0:0,-:0 -:0,-:-", boundString(ranges))
	ranges = mustGenerateReadShards(t, pam.GenerateReadShardsOpts{
		Range:     newRange(1, 2, 2, 100),
		NumShards: 1,
	}, pamPath)
	assert.Equal(t, "1:2,2:100", boundString(ranges))
	ranges = mustGenerateReadShards(t, pam.GenerateReadShardsOpts{
		Range:                              newRange(1, 2, 2, 100),
		NumShards:                          1,
		AlwaysSplitMappedAndUnmappedCoords: true,
	}, pamPath)
	assert.Equal(t, "1:2,2:100", boundString(ranges))
}

func boundString(bounds []biopb.CoordRange) string {
	s := make([]string, len(bounds))
	for i, b := range bounds {
		s[i] = pam.CoordRangePathString(b)
	}
	return strings.Join(s, " ")
}

func TestConvert(t *testing.T) {
	tempDir, cleanup := testutil.TempDir(t, "", "")
	defer cleanup()

	bamPath := testutil.GetFilePath("@grailgo//bio/encoding/bam/testdata/170614_WGS_LOD_Pre_Library_B3_27961B_05.merged.10000.bam")
	pamPath := newPAMPath(bamPath, tempDir)

	// The bam file is 2.8MB, so with 1MB shard size, we expect three PAM
	// shard files.
	require.NoError(t, converter.ConvertToPAM(pam.WriteOpts{}, pamPath, bamPath, "", 1<<20))
	verifyPAM(t, pam.ReadOpts{}, pamPath, bamPath)
	indexes, err := pam.ListIndexes(pamPath)
	require.NoError(t, err)
	require.Equal(t, len(indexes), 3, "Index:", indexes)

	// Try 256KB shard size.
	require.NoError(t, converter.ConvertToPAM(pam.WriteOpts{}, pamPath, bamPath, "", 1<<18))
	verifyPAM(t, pam.ReadOpts{}, pamPath, bamPath)
	indexes, err = pam.ListIndexes(pamPath)
	require.NoError(t, err)
	require.Equal(t, len(indexes), 11, "Index:", indexes)
}

func TestSharder1(t *testing.T) {
	tempDir, cleanup := testutil.TempDir(t, "", "")
	defer cleanup()

	// Create three PAM rowshards.
	bamPath := testutil.GetFilePath("@grailgo//bio/encoding/bam/testdata/170614_WGS_LOD_Pre_Library_B3_27961B_05.merged.10000.bam")
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
	shards := mustGenerateReadShards(t, pam.GenerateReadShardsOpts{NumShards: 1}, pamPath)
	require.Equal(t, 3, len(shards))
	verifyPAMWithShardedReader(t, pam.ReadOpts{}, pamPath, bamPath, shards)

	// The same test, but specify the params via BytesPerShard.
	shards = mustGenerateReadShards(t, pam.GenerateReadShardsOpts{BytesPerShard: math.MaxInt64}, pamPath)
	require.Equal(t, 3, len(shards))
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
	in := mustOpenBAM(t, testutil.GetFilePath("@grailgo//bio/encoding/bam/testdata/test.bam"))
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
		require.NoErrorf(st.t, err, "ref=%v, pos=%d", ref, pos)
		w.Write(rec)
		st.recs = append(st.recs, rec)
	}
	require.NoError(st.t, w.Close())
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
			require.Truef(st.t, n < len(st.recs), "n=%d, len=%d", n, len(st.recs))
			require.Equalf(st.t, st.recs[n].String(), rec.String(), "n=%d", n)
			n++
		}
		require.NoError(st.t, r.Close())
	}
	require.Equal(st.t, len(st.recs), n)
}

// Create a synthetic PAM file where reads have unique coordinates.
func TestSyntheticUniquePositions(t *testing.T) {
	st := newSyntheticTester(t)
	defer st.cleanup()

	const nRecords = 10000
	writeOpts := pam.WriteOpts{MaxBufSize: 1024}
	ref := st.header.Refs()[0]
	tmpPAMPath := st.generatePAM(writeOpts, nRecords, func(index int) (*sam.Reference, int) { return ref, index })
	shards := mustGenerateReadShards(t, pam.GenerateReadShardsOpts{NumShards: 16}, tmpPAMPath)
	assert.Equalf(t, 16, len(shards), "Shards: %+v", shards)
	for _, shard := range shards {
		require.Truef(t, shard.Start.Seq == 0 && shard.Limit.Seq == 0, "Shard: %+v", shard)
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
	shards := mustGenerateReadShards(t, pam.GenerateReadShardsOpts{NumShards: 16}, tmpPAMPath)
	assert.Equalf(t, 1, len(shards), "Shards: %+v", shards)
	require.Truef(t, shards[0].Start.Seq == 0 && shards[0].Limit.Seq == 0, "Shards: %+v", shards)
	st.verifyPAM(shards)

	// Allow splitting positions.
	const numShards = 16
	shards = mustGenerateReadShards(t, pam.GenerateReadShardsOpts{
		SplitMappedCoords:   true,
		SplitUnmappedCoords: true,
		NumShards:           numShards}, tmpPAMPath)
	assert.Equal(t, numShards, len(shards))
	expectedRecordsPerShard := nRecords / numShards
	for i, shard := range shards {
		if shard.Limit.RefId == 0 {
			assert.Equalf(t, int(shard.Start.RefId), 0, "shard=%v i=%v shards=%v", shard, i, shards)
			assert.Equalf(t, int(shard.Limit.Pos), 0, "shard=%v i=%v shards=%v", shard, i, shards)
			assert.True(t, shard.Limit.Seq > 0, "shard=%v i=%v shards=%v", shard, i, shards)
			nRecords := float64(shard.Limit.Seq - shard.Start.Seq)
			assert.True(t,
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
	shards := mustGenerateReadShards(t, pam.GenerateReadShardsOpts{
		SplitMappedCoords:   false,
		SplitUnmappedCoords: true,
		NumShards:           16}, tmpPAMPath)
	assert.True(t, len(shards) > 1 && len(shards) < 16, shards)
	for i, shard := range shards {
		assert.Equalf(t,
			shard.Start.Seq != 0,
			shard.Start.RefId == biopb.UnmappedRefID,
			"Shard %d: %+v", i, shard)
	}
	st.verifyPAM(shards)

	shards = mustGenerateReadShards(t, pam.GenerateReadShardsOpts{
		SplitMappedCoords:   true,
		SplitUnmappedCoords: false,
		NumShards:           16}, tmpPAMPath)
	assert.True(t, len(shards) > 1 && len(shards) < 16, shards)
	t.Logf("Found %v shards, %+v", len(shards), shards)
	for i, shard := range shards {
		assert.Equalf(t,
			shard.Limit.Seq != 0,
			shard.Limit.RefId != biopb.UnmappedRefID,
			"Shard %d: %+v", i, shard)
	}
	st.verifyPAM(shards)
}

func TestShardedUnmappedReads(t *testing.T) {
	tempDir, cleanup := testutil.TempDir(t, "", "")
	defer cleanup()
	bamPath := testutil.GetFilePath("@grailgo//bio/encoding/bam/testdata/170614_WGS_LOD_Pre_Library_B3_27961B_05.merged.10000.bam")
	pamPath := filepath.Join(tempDir, "test")
	generatePAM(t, pam.WriteOpts{MaxBufSize: 10000}, pamPath, bamPath)

	// shards0 won't split unmapped sequences.
	shards0 := mustGenerateReadShards(t, pam.GenerateReadShardsOpts{NumShards: 2000}, pamPath)
	// shards1 will split unmapped sequences.
	shards1 := mustGenerateReadShards(t, pam.GenerateReadShardsOpts{SplitUnmappedCoords: true, NumShards: 2000}, pamPath)

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
		require.NoError(b, converter.ConvertToPAM(pam.WriteOpts{}, pamPath, *bamFlag, "", math.MaxInt64))
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
		require.NoError(b, converter.ConvertToPAM(pam.WriteOpts{}, pamPath, *bamFlag, "", math.MaxInt64))
	} else {
		vlog.Infof("Reusing PAM files %v", pamPath)
	}
	b.StartTimer()
	for n := 0; n < b.N; n++ {
		opts := pam.ReadOpts{}
		if *dropFieldsFlag != "" {
			for _, fieldName := range strings.Split(*dropFieldsFlag, ",") {
				f, err := gbam.ParseFieldType(fieldName)
				require.NoError(b, err)
				opts.DropFields = append(opts.DropFields, f)
			}
		}
		if !*unmappedFlag {
			vlog.Infof("Skipping unmapped reads")
			opts.Range = gbam.MappedRange
		}
		bounds, err := pam.GenerateReadShards(pam.GenerateReadShardsOpts{Range: opts.Range}, pamPath)
		require.NoError(b, err)
		boundCh := make(chan biopb.CoordRange, len(bounds))
		for _, r := range bounds {
			boundCh <- r
		}
		close(boundCh)
		totalRecs := int64(0)
		wg := sync.WaitGroup{}
		for i := 0; i < runtime.NumCPU(); i++ {
			wg.Add(1)
			go func() {
				defer wg.Done()
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
					require.NoError(b, r.Close())
					atomic.AddInt64(&totalRecs, nRecs)
					end := time.Now()
					vlog.Infof("Finish read %+v, %d recs %d ms", bound, nRecs, end.Sub(start)/time.Millisecond)
				}
			}()
		}
		wg.Wait()
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
		require.NoError(b, err)
		shardList, err := gbam.GetPositionBasedShards(header, 100000, 0, *unmappedFlag)
		require.NoError(b, err)
		shardCh := gbam.NewShardChannel(shardList)
		require.NoError(b, err)
		parallelism := runtime.NumCPU()
		totalRecs := int64(0)
		wg := sync.WaitGroup{}
		for i := 0; i < parallelism; i++ {
			wg.Add(1)
			go func() {
				defer wg.Done()
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
					require.NoError(b, iter.Close())
				}
			}()
		}
		wg.Wait()
		require.NoError(b, provider.Close())
		vlog.Infof("%v: read %d records", *bamFlag, totalRecs)
	}
}

func mustCreate(t *testing.T, path string) {
	fd, err := os.Create(path)
	require.NoError(t, err)
	fd.WriteString("foo")
	require.NoError(t, fd.Close())
}

func TestRemove(t *testing.T) {
	tempDir, cleanup := testutil.TempDir(t, "", "")
	defer cleanup()

	// Removing non-existing file should get no error.
	require.NoError(t, pam.Remove(filepath.Join(tempDir, "gg")))

	path := filepath.Join(tempDir, "f")
	require.NoError(t, os.MkdirAll(path, 0777))
	mustCreate(t, filepath.Join(path, "0:0,2:123.index"))
	mustCreate(t, filepath.Join(path, "0:0,2:123.aux"))

	_, err := os.Stat(path)
	require.NoError(t, err)
	require.NoError(t, pam.Remove(path))
	_, err = os.Stat(path)
	require.True(t, os.IsNotExist(err))
}

func TestListIndexes(t *testing.T) {
	tempDir, cleanup := testutil.TempDir(t, "", "")
	defer cleanup()
	path := filepath.Join(tempDir, "f")
	require.NoError(t, os.MkdirAll(path, 0777))
	mustCreate(t, path+"/0:0,2:123.index")
	mustCreate(t, path+"/2:123,10:200.index")
	mustCreate(t, path+"/10:200,-:-.index")

	indexes, err := pam.ListIndexes(path)
	require.NoError(t, err)
	expected := []biopb.CoordRange{
		biopb.CoordRange{biopb.Coord{0, 0, 0}, biopb.Coord{2, 123, 0}},
		biopb.CoordRange{biopb.Coord{2, 123, 0}, biopb.Coord{10, 200, 0}},
		biopb.CoordRange{biopb.Coord{10, 200, 0}, biopb.Coord{biopb.InfinityRefID, biopb.InfinityPos, 0}}}
	assert.Equal(t, len(indexes), len(expected))
	for i, e := range expected {
		assert.Equal(t, e, indexes[i].Range)
	}

	indexes, err = pam.ListIndexes(path + "/")
	require.NoError(t, err)
	assert.Equal(t, len(indexes), len(expected))
	for i, e := range expected {
		assert.Equal(t, e, indexes[i].Range)
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
