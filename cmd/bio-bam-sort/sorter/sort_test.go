package sorter

import (
	"bytes"
	"flag"
	"fmt"
	"io"
	"log"
	"math"
	"math/rand"
	"os"
	"path/filepath"
	"runtime"
	"sort"
	"strings"
	"sync"
	"sync/atomic"
	"syscall"
	"testing"

	"github.com/biogo/hts/bam"
	"github.com/biogo/hts/bgzf"
	"github.com/biogo/hts/sam"
	"github.com/grailbio/base/errorreporter"
	"github.com/grailbio/base/grail"
	gbam "github.com/grailbio/bio/encoding/bam"
	"github.com/grailbio/bio/encoding/bamprovider"
	"github.com/grailbio/testutil"
	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
	"v.io/x/lib/gosh"
	"v.io/x/lib/lookpath"
)

func headerString(t *testing.T, header *sam.Header) string {
	s, err := header.MarshalText()
	require.NoError(t, err)
	return "[" + string(s) + "]"
}

func TestCompareKeys(t *testing.T) {
	ref, err := sam.NewReference("chr1", "0", "human", 10000, nil, nil)
	require.NoError(t, err)
	_, err = sam.NewHeader(nil, []*sam.Reference{ref})
	require.NoError(t, err)

	testRecord := func(name string, alignPos int, flags sam.Flags) *sam.Record {
		seq := []byte{'A'}
		qual := []byte{30}
		rec, err := sam.NewRecord(name, ref, ref, alignPos, alignPos+223 /*matepos*/, 0 /*tlen*/, 1 /*mapq*/, nil /*cigar*/, seq, qual, nil)
		require.NoError(t, err)
		rec.Flags = flags
		return rec
	}

	type testCase struct {
		r0   *sam.Record
		seq0 uint64
		r1   *sam.Record
		seq1 uint64
		e    int // expected result of sortKey.compare.
	}

	for _, tc := range []testCase{
		{
			testRecord("r01", 1, sam.Flags(0)), 123,
			testRecord("r01", 1, sam.Flags(0)), 123,
			0},
		{
			testRecord("r01", 1, sam.Flags(0)), 234,
			testRecord("r01", 1, sam.Flags(0)), 123,
			1},
		{
			testRecord("r01", 1, sam.Flags(0)), 123,
			testRecord("r01", 2, sam.Flags(0)), 123,
			-1},
		{
			testRecord("r01", 1, sam.Reverse), 123,
			testRecord("r01", 2, sam.Flags(0)), 123,
			-1},
		{
			testRecord("r01", 1, sam.Flags(0)), 123,
			testRecord("r02", 1, sam.Reverse), 123,
			-1},
	} {
		r0Key := makeSortKey(tc.r0, tc.seq0)
		r1Key := makeSortKey(tc.r1, tc.seq1)
		// Test trivial comparison ops.
		assert.Equalf(t, 0, r0Key.compare(r0Key), "R0: %+v", r0Key)
		assert.Equalf(t, 0, r1Key.compare(r1Key), "R0: %+v", r1Key)
		assert.Equalf(t, r0Key.compare(r1Key), tc.e,
			"R0: %+v, r1: %+v, compare %d, expected %d",
			r0Key, r1Key, r0Key.compare(r1Key), tc.e)
	}
}

const invalidRefID = -99 // -1 is taken by unmapped reads, so use -99.
type recordPos struct {
	refID int // ID of the reference sequence
	pos   int // aligned position
}

// Read the header and list of records from a BAM or PAM file.
func readRecords(t *testing.T, path string) (*sam.Header, []*sam.Record) {
	records := []*sam.Record{}
	if strings.HasSuffix(path, ".bam") {
		// BAM used in this test is unindexed, so we can't use a Provider.
		in, err := os.Open(path)
		require.NoErrorf(t, err, "Open %v", path)
		defer in.Close()
		r, err := bam.NewReader(in, 16)
		require.NoErrorf(t, err, "NewReader %v", path)
		for {
			rec, err := r.Read()
			if err != nil {
				require.True(t, err == io.EOF)
				break
			}
			records = append(records, rec)
		}
		require.NoError(t, r.Close())
		return r.Header(), records
	}
	provider := bamprovider.NewProvider(path)
	header, err := provider.GetHeader()
	require.NoError(t, err)
	iter := provider.NewIterator(gbam.UniversalShard(header))
	for iter.Scan() {
		records = append(records, iter.Record())
	}
	require.NoError(t, iter.Close())
	require.NoError(t, provider.Close())
	return header, records
}

// Hash-split recs into "n" sortshards.
func createSortShards(t *testing.T, header *sam.Header, recs []*sam.Record, tempDir string, n int) []string {
	shards := make([]string, n)
	for i := 0; i < len(shards); i++ {
		shards[i] = fmt.Sprintf("%s/shard%d", tempDir, i)
		log.Printf("Creating sort shard %v", shards[i])
		sorter := NewSorter(shards[i], header, SortOptions{ShardIndex: uint32(i)})
		for j, rec := range recs {
			if j%len(shards) == i {
				sorter.AddRecord(rec)
			}
		}
		require.NoError(t, sorter.Close())
		log.Printf("Finished creating sort shard %v", shards[i])
	}
	return shards
}

// Verify that the two files store the same set of BAM records in the same
// order. Both files can either be a BAM or a PAM file.
func compareFiles(t *testing.T, orgPath, newPath string) {
	orgHeader, orgRecs := readRecords(t, orgPath)
	newHeader, newRecs := readRecords(t, newPath)

	require.Equal(t, headerString(t, orgHeader), headerString(t, newHeader))
	require.Equal(t, len(orgRecs), len(newRecs))
	// Read the result, make sure that records are actually sorted, and that the
	// file contains the same records as the input.
	orgNames := make(map[string]*sam.Record)
	newNames := make(map[string]*sam.Record)
	for i := 0; i < len(orgRecs); i++ {
		r0 := orgRecs[i]
		r1 := newRecs[i]
		orgNames[r0.Name] = r0
		newNames[r1.Name] = r1
		require.Equalf(t, r0.Ref.Name(), r1.Ref.Name(), "%d: Records: org:\n%v, new\n%v", i, r0, r1)
		require.Equalf(t, r0.Pos, r1.Pos, "%d: Records: org:\n%v, new\n%v", i, r0, r1)
	}

	// Make sure that two files srtore the same set of records.
	require.Equal(t, len(orgNames), len(newNames))
	for name := range orgNames {
		_, ok := newNames[name]
		assert.Truef(t, ok, "Name %d not found in new", name)
	}
}

var testBAMFiles = []string{
	testutil.GetFilePath("@grailgo//bio/encoding/bam/testdata/test.bam"),
	testutil.GetFilePath("@grailgo//bio/encoding/bam/testdata/test-unmapped.bam"),
	testutil.GetFilePath("@grailgo//bio/encoding/bam/testdata/170614_WGS_LOD_Pre_Library_B3_27961B_05.merged.10000.bam"),
	// Sad. This file is outside the Bazel workspace.
	// testutil.GetFilePath("/tbr/pecan/tests/data/stages/HD753_titrated_fusions.bam")
}

func TestSortBAMEndToEnd(t *testing.T) {
	tempDir, cleanup := testutil.TempDir(t, "", "")
	defer testutil.NoCleanupOnError(t, cleanup)

	// Feed the records in "inPath" in random order to the sorter+merger.  Verify
	// that the results match the original bamPath.
	for _, bamPath := range testBAMFiles {
		newBAMPath := filepath.Join(tempDir, "test.bam")
		log.Printf("Sorting %v to %v", bamPath, newBAMPath)

		header, recs := shuffleRecords(t, bamPath, "")
		shards := createSortShards(t, header, recs, tempDir, 3)
		require.NoError(t, BAMFromSortShards(shards[:], newBAMPath))
		compareFiles(t, bamPath, newBAMPath)
	}
}

func TestSortPAMEndToEnd(t *testing.T) {
	tempDir, cleanup := testutil.TempDir(t, "", "")
	defer testutil.NoCleanupOnError(t, cleanup)

	// Feed the records in "inPath" in random order to the sorter+merger.  Verify
	// that the results match the original bamPath.
	for _, bamPath := range testBAMFiles {
		newPAMPath := filepath.Join(tempDir, "test.pam")
		log.Printf("Sorting %v to %v", bamPath, newPAMPath)
		header, recs := shuffleRecords(t, bamPath, "")
		shards := createSortShards(t, header, recs, tempDir, 3)
		require.NoError(t, PAMFromSortShards(shards[:], newPAMPath, math.MaxInt64))
		compareFiles(t, bamPath, newPAMPath)

		// Test PAM sharding using the large test file
		if !strings.Contains(bamPath, "170614") {
			continue
		}
		newPAMPath = filepath.Join(tempDir, "sharded.pam")

		// Try a few different sharding params.
		require.NoError(t, PAMFromSortShards(shards[:], newPAMPath, 1024))
		compareFiles(t, bamPath, newPAMPath)
		require.NoError(t, PAMFromSortShards(shards[:], newPAMPath, 3000))
		compareFiles(t, bamPath, newPAMPath)
		require.NoError(t, PAMFromSortShards(shards[:], newPAMPath, 500))
		compareFiles(t, bamPath, newPAMPath)
	}
}

func TestEmpty(t *testing.T) {
	tempDir, cleanup := testutil.TempDir(t, "", "")
	defer testutil.NoCleanupOnError(t, cleanup)
	shardPath := filepath.Join(tempDir, "shard0")
	bamPath := filepath.Join(tempDir, "test.bam")

	header, err := sam.NewHeader(nil, nil)
	require.NoError(t, err)
	sorter := NewSorter(shardPath, header)
	require.NoError(t, sorter.Close())
	require.NoError(t, BAMFromSortShards([]string{shardPath}, bamPath))

	_, recs := readRecords(t, bamPath)
	require.Truef(t, len(recs) == 0, "Recs: %+v", recs)
}

func sortSAM(t *testing.T, opts SortOptions, shardPath, text string) {
	r, err := sam.NewReader(bytes.NewBuffer([]byte(text)))
	require.NoError(t, err)
	sorter := NewSorter(shardPath, r.Header(), opts)
	for {
		rec, err := r.Read()
		if rec == nil {
			require.Truef(t, err == io.EOF, "err: %v", err)
			break
		}
		require.NoError(t, err)
		sorter.AddRecord(rec)
	}
	require.NoError(t, sorter.Close())
}

func TestUnmapped(t *testing.T) {
	tempDir, cleanup := testutil.TempDir(t, "", "")
	defer testutil.NoCleanupOnError(t, cleanup)

	shard0 := fmt.Sprintf("%s/shard0", tempDir)
	sortSAM(t, SortOptions{}, shard0, `@HD	VN:1.3	SO:coordinate
read10	12	*	0	0	10M	*	0	0	ACGTACGTAC	ABCDEFGHIJ	RG:Z:NA12878
read11	12	*	0	0	10M	*	0	0	ACGTACGTAC	ABCDEFGHIJ	RG:Z:NA12879
`)
	errReporter := errorreporter.T{}
	pool := newSortShardBlockPool()
	shardReader := newSortShardReader(shard0, pool, &errReporter)

	header, err := mergeHeader([]*sortShardReader{shardReader})
	require.NoError(t, err)
	pamPath := filepath.Join(tempDir, "test.pam")
	generatePAMShard([]*sortShardReader{shardReader},
		pamPath, header, unmappedCoord, infinityCoord, pool, &errReporter)
	require.NoError(t, errReporter.Err())
	_, recs := readRecords(t, pamPath)
	log.Printf("Read recs: %v", recs)
	expected := []string{"read10", "read11"}
	n := 0
	for _, rec := range recs {
		log.Printf("Read rec: %v", rec)
		assert.Equalf(t, expected[n], rec.Name, "rec=%+v", rec)
		n++
	}
	assert.Equal(t, n, len(expected))
}

func TestMergeDifferentRefs(t *testing.T) {
	tempDir, cleanup := testutil.TempDir(t, "", "")
	defer testutil.NoCleanupOnError(t, cleanup)

	shard0 := fmt.Sprintf("%s/shard0", tempDir)
	shard1 := fmt.Sprintf("%s/shard1", tempDir)

	sortSAM(t, SortOptions{ShardIndex: 1}, shard0, `@HD	VN:1.3	SO:coordinate
@SQ	SN:chr1	LN:10000
read1	0	chr1	123	60	10M	=	456	20	AAAAAAAAAA	ABCDEFGHIJ	NM:i:1
`)

	sortSAM(t, SortOptions{ShardIndex: 2}, shard1, `@HD	VN:1.3	SO:coordinate
@SQ	SN:chr1	LN:10000
@SQ	SN:chr2	LN:9999
read2	0	chr1	100	60	10M	=	456	20	CCCCCCCCCC	ABCDEFGHIJ	NM:i:1
read3	0	chr2	123	60	10M	=	456	20	GGGGGGGGGG	ABCDEFGHIJ	NM:i:1
`)

	bamPath := filepath.Join(tempDir, "test.bam")
	require.NoError(t, BAMFromSortShards([]string{shard0, shard1}, bamPath))
	header, recs := readRecords(t, bamPath)
	log.Printf("Read recs: %v", recs)
	require.Contains(t, headerString(t, header), "@SQ	SN:chr1	LN:10000")
	require.Contains(t, headerString(t, header), "@SQ	SN:chr2	LN:9999")

	expected := []string{"read2", "read1", "read3"}
	n := 0
	for _, rec := range recs {
		log.Printf("Read rec: %v", rec)
		assert.Equalf(t, expected[n], rec.Name, "rec=%+v", rec)
		n++
	}
	assert.Equal(t, n, len(expected))
}

// Merge two shards with the same shard index.  The resulting BAM should contain
// all the records in the input shards, although the sort order of reads at the
// same coordinate is unspecified.
func TestMergeSameShardIndex(t *testing.T) {
	tempDir, cleanup := testutil.TempDir(t, "", "")
	defer testutil.NoCleanupOnError(t, cleanup)

	samText := `@HD	VN:1.3	SO:coordinate
@SQ	SN:chr1	LN:10000
read1	0	chr1	123	60	10M	=	456	20	AAAAAAAAAA	ABCDEFGHIJ	NM:i:1
`
	shard0 := fmt.Sprintf("%s/shard0", tempDir)
	shard1 := fmt.Sprintf("%s/shard1", tempDir)
	sortSAM(t, SortOptions{ShardIndex: 1}, shard0, samText)
	sortSAM(t, SortOptions{ShardIndex: 1}, shard1, samText)

	bamPath := filepath.Join(tempDir, "test.bam")
	require.NoError(t, BAMFromSortShards([]string{shard0, shard1}, bamPath))
	header, recs := readRecords(t, bamPath)
	log.Printf("Read recs: %v", recs)
	require.Contains(t, headerString(t, header), "@SQ	SN:chr1	LN:10000")

	expected := []string{"read1", "read1"}
	n := 0
	for _, rec := range recs {
		log.Printf("Read rec: %v", rec)
		assert.Equalf(t, expected[n], rec.Name, "rec=%+v", rec)
		n++
	}
	assert.Equal(t, n, len(expected))
}

func runCmd(t *testing.T, sh *gosh.Shell, arg0 string, args ...string) {
	cmd := sh.Cmd(arg0, args...)
	cmd.Run()
	if cmd.Err != nil {
		t.Fatal(cmd.Err)
	}
}

// Sort the given file using sambamba and Sort(). Make sure that the results
// match.
func testAgainstSambamba(t *testing.T, bamPath string) {
	sh := gosh.NewShell(t)
	defer sh.Cleanup()
	if _, err := lookpath.Look(sh.Vars, "sambamba"); err != nil {
		t.Skipf("sambamba not found on the machine. Skipping the test for %v", bamPath)
		return
	}
	log.Printf("Comparing outputs of sambamba and Sort() for %s", bamPath)
	unsortedFile := sh.MakeTempFile()

	// Create a file that's sorted by querynames.
	runCmd(t, sh, "sambamba", "sort", "-N", "-o", unsortedFile.Name(), bamPath)

	// Sort by position using sambamba, then convert the output into SAM.
	smbSortedFile := sh.MakeTempFile()
	runCmd(t, sh, "sambamba", "sort", "-o", smbSortedFile.Name(), unsortedFile.Name())
	smbSortedTextFile := sh.MakeTempFile()
	runCmd(t, sh, "sambamba", "view", "-h", "-f", "sam", "-o", smbSortedTextFile.Name(), smbSortedFile.Name())
	t.Logf("Sorted %v into %v using sambamba", unsortedFile.Name(), smbSortedFile.Name())

	bamReader, err := bam.NewReader(unsortedFile, 1)
	require.NoError(t, err)

	// Sort by position using Sort(), then convert the output into SAM.
	shardFile := sh.MakeTempFile()
	sorter := NewSorter(shardFile.Name(), bamReader.Header())
	for {
		rec, err := bamReader.Read()
		if rec == nil {
			break
		}
		require.NoError(t, err)
		sorter.AddRecord(rec)
	}
	require.NoError(t, bamReader.Close())
	require.NoError(t, sorter.Close())

	sortedFile := sh.MakeTempFile()
	require.NoError(t, BAMFromSortShards([]string{shardFile.Name()}, sortedFile.Name()))
	t.Logf("Sorted %v into %v using Sort", unsortedFile.Name(), sortedFile.Name())

	sortedTextFile := sh.MakeTempFile()
	runCmd(t, sh, "sambamba", "view", "-h", "-f", "sam", "-o", sortedTextFile.Name(), sortedFile.Name())

	diff := sh.Cmd("diff", "-I", "^@PG", smbSortedTextFile.Name(), sortedTextFile.Name()).Stdout()
	if diff != "" {
		t.Errorf("Diff: %v", diff)
	}
}

func TestAgainstSamBamba(t *testing.T) {
	for _, path := range testBAMFiles {
		testAgainstSambamba(t, path)
	}
}

func TestSortKey(t *testing.T) {
	ref0, err := sam.NewReference("chr1", "", "", 1000, nil, nil)
	require.NoError(t, err)
	ref1, err := sam.NewReference("chr2", "", "", 1000, nil, nil)
	require.NoError(t, err)
	_, err = sam.NewHeader(nil, []*sam.Reference{ref0, ref1})
	require.NoError(t, err)
	refid, pos, reverse := parseCoord(coordFromRecord(&sam.Record{Ref: ref0, Pos: 0}))
	assert.True(t, refid == 0 && pos == 0 && !reverse, "rec", refid, pos, reverse)

	refid, pos, reverse = parseCoord(coordFromRecord(&sam.Record{Ref: ref1, Pos: 1}))
	assert.True(t, refid == 1 && pos == 1 && !reverse, "rec", refid, pos, reverse)

	refid, pos, reverse = parseCoord(coordFromRecord(&sam.Record{Ref: nil, Pos: -1}))
	assert.True(t, refid == -1 && pos == -1 && !reverse, "rec", refid, pos, reverse)
}

func increaseRlimit() error {
	var l syscall.Rlimit
	err := syscall.Getrlimit(syscall.RLIMIT_NOFILE, &l)
	if err != nil {
		return err
	}
	l.Cur = l.Max
	return syscall.Setrlimit(syscall.RLIMIT_NOFILE, &l)
}

var (
	bamFlag             = flag.String("bam", "/scratch-nvme/bam/small.bam", "File to sort in benchmark")
	baiFlag             = flag.String("bai", "", "BAM index file to use during benchmarking. May be empty")
	sortShardDirFlag    = flag.String("sortsharddir", "/scratch-nvme/sortshard", "Directory to store outputs during benchmarking")
	tmpDirFlag          = flag.String("tmpdir", "/scratch-nvme/tmp", "Tmp dir to use")
	compressFlag        = flag.Bool("compress", true, "Compress tmp files using snappy")
	splitFlag           = flag.Bool("split", false, "If set, split the BAM file into multiple shards")
	bytesPerShardFlag   = flag.Int64("bytes-per-shard", 8<<30, "BAM/PAM read shard size")
	recsPerShardFlag    = flag.Int64("recs-per-shard", 32<<20, "Target sortshard size")
	sortBatchSizeFlag   = flag.Int("sort-batch-size", 1<<20, "Sort batch size during benchmark")
	sortParallelismFlag = flag.Int("parallelism", 16, "Number of sortshards to create in parallel")
)

// generateShards compute Shard boundarie for the given BAM file, where each
// shard is close to bytesPerShard , ignoring reference (chromosome) boundaries.
func generateShards(b *testing.B, bamPath, baiPath string, bytesPerShard int64) []gbam.Shard {
	shards := []gbam.Shard{}

	in, err := os.Open(bamPath)
	require.NoError(b, err)
	defer in.Close()
	bamr, err := bam.NewReader(in, 1)
	defer bamr.Close()
	header := bamr.Header()
	require.NoError(b, err)
	if baiPath == "" {
		baiPath = bamPath + ".bai"
	}
	indexIn, err := os.Open(baiPath)
	require.NoError(b, err)
	defer indexIn.Close()
	index, err := gbam.ReadIndex(indexIn)
	require.NoError(b, err)

	// Get chunks
	chunks := []bgzf.Offset{}
	for _, chunksForRef := range index.AllOffsets() {
		for _, c := range chunksForRef {
			chunks = append(chunks, c)
		}
	}
	sort.SliceStable(chunks,
		func(a, b int) bool {
			return chunks[a].File < chunks[b].File
		})
	require.True(b, len(chunks) > 0)
	startChunk := chunks[0]
	for _, chunk := range chunks[1:] {
		if chunk.File-startChunk.File >= bytesPerShard {
			startCoord, err := gbam.GetCoordAtOffset(bamr, startChunk)
			require.NoError(b, err)
			limitCoord, err := gbam.GetCoordAtOffset(bamr, chunk)
			require.NoError(b, err)
			shard := gbam.Shard{
				StartRef: header.Refs()[startCoord.RefId],
				Start:    int(startCoord.Pos),
				EndRef:   header.Refs()[limitCoord.RefId],
				End:      int(limitCoord.Pos),
				ShardIdx: len(shards),
			}
			startChunk = chunk
			shards = append(shards, shard)
		}
	}
	startCoord, err := gbam.GetCoordAtOffset(bamr, startChunk)
	require.NoError(b, err)
	shard := gbam.Shard{
		StartRef: header.Refs()[startCoord.RefId],
		Start:    int(startCoord.Pos),
		EndRef:   nil,
		End:      math.MaxInt32,
		ShardIdx: len(shards),
	}
	return append(shards, shard)
}

// Read a BAM file then random-shuffle the records within.
func shuffleRecords(t require.TestingT, bamPath, baiPath string) (*sam.Header, []*sam.Record) {
	provider := bamprovider.NewProvider(bamPath, bamprovider.ProviderOpts{Index: baiPath})
	header, err := provider.GetHeader()
	require.NoError(t, err)

	iter := provider.NewIterator(gbam.UniversalShard(header))
	recs := []*sam.Record{}
	for iter.Scan() {
		recs = append(recs, iter.Record())
	}
	require.NoError(t, iter.Close())
	require.NoError(t, provider.Close())
	r := rand.New(rand.NewSource(0))
	for i := range recs {
		j := r.Intn(i + 1)
		recs[i], recs[j] = recs[j], recs[i]
	}
	return header, recs
}

// Create a sortshard file from a set of records.
func benchmarkSortShard(b *testing.B, recs []*sam.Record, header *sam.Header, shardPath string, options SortOptions) {
	if _, err := os.Stat(shardPath); err == nil {
		return
	}
	log.Printf("Creating sortshard %v. This may take minutes", shardPath)
	s := NewSorter(shardPath, header, options)
	for i := range recs {
		s.AddRecord(recs[i])
		recs[i] = nil
	}
	require.NoError(b, s.Close())
}

// Example of sorting & merging a large BAM file for benchmarking:
//
// Example of splitting a BAM file ~60 ways and creating a sortshard file for each.
//
//    bazel-bin/cmd/bio-bam-sort/sorter/go_default_test --test.run xxxx -test.v -test.bench Split -logtostderr --split -bam /scratch-nvme/bam/CNVS-NORM-110033752-cfDNA-WGBS-Rep1.bam
//
// Example of merging the sortshard files created above.
//
//    bazel-bin/cmd/bio-bam-sort/bio-bam-sort -profile-interval-s 60 -heap-profile heap -cpu-profile cpu --bam /scratch-nvme/tmp/foo.bam /scratch-nvme/sortshard/*.sortshard
//
func BenchmarkSplitBAM(b *testing.B) {
	if !*splitFlag {
		b.Skipf("--split not set")
		return
	}
	require.NoError(b, increaseRlimit())
	p := bamprovider.NewProvider(*bamFlag, bamprovider.ProviderOpts{Index: *baiFlag})
	header, err := p.GetHeader()
	require.NoError(b, err)
	shards := generateShards(b, *bamFlag, *baiFlag, *bytesPerShardFlag)

	sortCh := make(chan *sam.Record, 1<<20)
	sortWg := sync.WaitGroup{}
	shardSeq := int32(0)
	for i := 0; i < *sortParallelismFlag; i++ {
		sortWg.Add(1)
		go func() {
			var sorter *Sorter
			var nRecsInShard int64
			for rec := range sortCh {
				if sorter == nil {
					seq := atomic.AddInt32(&shardSeq, 1)
					shardPath := fmt.Sprintf("%s/test-%04d.sortshard", *sortShardDirFlag, seq)
					log.Printf("Creating %v", shardPath)
					sorter = NewSorter(shardPath, header, SortOptions{
						TmpDir:             *tmpDirFlag,
						NoCompressTmpFiles: !*compressFlag,
						ShardIndex:         uint32(seq),
					})
				}
				sorter.AddRecord(rec)
				nRecsInShard++
				if nRecsInShard >= *recsPerShardFlag {
					require.NoError(b, sorter.Close())
					sorter = nil
					nRecsInShard = 0
				}
			}
			require.NoError(b, sorter.Close())
			sortWg.Done()
		}()
	}

	// Do N-way sharded reads. Each reader puts records in random sortshards.
	ch := gbam.NewShardChannel(shards)
	readWg := sync.WaitGroup{}
	totalNumRecs := int64(0)
	for i := 0; i < runtime.NumCPU(); i++ {
		readWg.Add(1)
		go func() {
			n := int64(0)
			for shard := range ch {
				iter := p.NewIterator(shard)
				for iter.Scan() {
					sortCh <- iter.Record()
					n++
					if n == 1<<20 {
						if total := atomic.AddInt64(&totalNumRecs, n); total%(16<<20) == 0 {
							log.Printf("Read %d records", total)
						}
						n = 0
					}
				}
				require.NoError(b, iter.Close())
			}
			readWg.Done()
		}()
	}
	readWg.Wait()
	close(sortCh)
	sortWg.Wait()
	require.NoError(b, p.Close())
}

func BenchmarkSort(b *testing.B) {
	require.NoError(b, increaseRlimit())
	outPath := fmt.Sprintf("%s/sorttest.sortshard", *sortShardDirFlag)
	recs, header := shuffleRecords(b, *bamFlag, *baiFlag)
	for n := 0; n < b.N; n++ {
		_ = os.Remove(outPath) // ignore errors
		benchmarkSortShard(b, header, recs, outPath, SortOptions{
			TmpDir:             *tmpDirFlag,
			NoCompressTmpFiles: !*compressFlag,
			SortBatchSize:      *sortBatchSizeFlag,
		})
	}
}

func BenchmarkMerge(b *testing.B) {
	require.NoError(b, increaseRlimit())
	recs, header := shuffleRecords(b, *bamFlag, *baiFlag)
	shardPaths := []string{
		fmt.Sprintf("%s/test-0000.sortshard", *sortShardDirFlag),
		fmt.Sprintf("%s/test-0001.sortshard", *sortShardDirFlag),
	}
	wg := sync.WaitGroup{}
	for _, shardPath := range shardPaths {
		wg.Add(1)
		go func(shardPath string) {
			benchmarkSortShard(b, header, recs, shardPath, SortOptions{
				TmpDir:             *tmpDirFlag,
				NoCompressTmpFiles: !*compressFlag,
			})
			wg.Done()
		}(shardPath)
	}
	wg.Wait()

	outPath := filepath.Join(filepath.Dir(*bamFlag), "out.bam")
	for n := 0; n < b.N; n++ {
		require.NoError(b, BAMFromSortShards(shardPaths, outPath))
	}
}

func TestMain(m *testing.M) {
	// Enable the profile handlers.
	shutdown := grail.Init()
	defer shutdown()
	os.Exit(m.Run())
}
