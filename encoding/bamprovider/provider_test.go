package bamprovider_test

import (
	"flag"
	"fmt"
	"io"
	"math"
	"math/rand"
	"os"
	"path/filepath"
	"runtime"
	"sync"
	"testing"

	"github.com/biogo/hts/bam"
	"github.com/biogo/hts/sam"
	"github.com/grailbio/base/grail"
	"github.com/grailbio/bio/biopb"
	gbam "github.com/grailbio/bio/encoding/bam"
	"github.com/grailbio/bio/encoding/bamprovider"
	"github.com/grailbio/bio/encoding/converter"
	"github.com/grailbio/bio/encoding/pam"
	"github.com/grailbio/internal/testutil"
	"github.com/stretchr/testify/require"
	"v.io/x/lib/vlog"
)

func TestMain(m *testing.M) {
	shutdown := grail.Init()
	status := m.Run()
	shutdown()
	os.Exit(status)
}

func doRead(t *testing.T, path string) []string {
	p := bamprovider.NewProvider(path)
	shards, err := p.GenerateShards(bamprovider.GenerateShardsOpts{
		IncludeUnmapped: true})
	require.NoError(t, err)

	var names []string
	// Repeat the test to test iterator-reuse code path.
	for i := 0; i < 3; i++ {
		names = []string{}
		for _, shard := range shards {
			i := p.NewIterator(shard)
			for i.Scan() {
				names = append(names, i.Record().Name)
			}
			require.NoError(t, i.Err())
			require.NoError(t, i.Close())
		}
		require.NoError(t, p.Close())
	}
	return names
}

func TestPAMSmall(t *testing.T) {
	tmpDir, cleanup := testutil.TempDir(t, "", "")
	defer cleanup()
	bamPath := testutil.GetFilePath("@grailgo//bio/encoding/bam/testdata/test-unmapped.bam")
	pamPath := filepath.Join(tmpDir, "test-unmapped.pam")
	require.NoError(t, converter.ConvertToPAM(pam.WriteOpts{}, pamPath, bamPath, "", math.MaxInt64))
	require.Equal(t, doRead(t, pamPath),
		[]string{"read1", "read2", "read3", "read10", "read10"})
}

func TestPAMLarge(t *testing.T) {
	tmpDir, cleanup := testutil.TempDir(t, "", "")
	defer cleanup()
	bamPath := testutil.GetFilePath("@grailgo//bio/encoding/bam/testdata/170614_WGS_LOD_Pre_Library_B3_27961B_05.merged.10000.bam")
	pamPath := filepath.Join(tmpDir, "large.pam")
	// The bam file is 2.8MB, so with 1MB shard size, we expect 3 shards.
	require.NoError(t, converter.ConvertToPAM(pam.WriteOpts{}, pamPath, bamPath, "", 1<<20))

	p := bamprovider.NewProvider(pamPath)
	shards, err := p.GetFileShards()
	require.NoError(t, err)
	header, err := p.GetHeader()
	require.NoError(t, err)
	const n = 3
	require.Equal(t, n, len(shards), "Shards:", shards)
	require.Equal(t, shards[0].StartRef, header.Refs()[0], "Shards:", shards)
	require.Equal(t, shards[0].Start, 0, "Shards:", shards)
	require.True(t, shards[n-1].EndRef == nil, "Shards:", shards) // the last shard is always unmapped
	require.Equal(t, shards[n-1].End, math.MaxInt32, "Shards:", shards)
}

func TestError(t *testing.T) {
	for _, test := range []struct{ path, errRe string }{
		{"nonexistent.bam", "no such file"},
		{"nonexistint.pam", "no index files found"},
	} {
		p := bamprovider.NewProvider(test.path)
		_, err := p.GenerateShards(bamprovider.GenerateShardsOpts{IncludeUnmapped: true})
		require.Regexp(t, test.errRe, err.Error())

		iter := p.NewIterator(gbam.Shard{StartRef: nil, EndRef: nil, Start: 0, End: 1})
		require.Regexp(t, test.errRe, iter.Close())
		require.Regexp(t, test.errRe, p.Close().Error())
	}
}

func TestBAM(t *testing.T) {
	require.Equal(t,
		doRead(t, testutil.GetFilePath("@grailgo//bio/encoding/bam/testdata/test.bam")),
		[]string{"read1", "read2", "read3"})
}

func TestBAMUnmapped(t *testing.T) {
	require.Equal(t,
		doRead(t, testutil.GetFilePath("@grailgo//bio/encoding/bam/testdata/test-unmapped.bam")),
		[]string{"read1", "read2", "read3", "read10", "read10"})
}

func TestBAMUnmappedOnly(t *testing.T) {
	require.Equal(t,
		doRead(t, testutil.GetFilePath("@grailgo//bio/encoding/bam/testdata/test-unmapped-only.bam")),
		[]string{"read10", "read10"})
}

// Test reading random ranges.
func testRandom(t *testing.T, randomSeed int64) {
	tmpDir, cleanup := testutil.TempDir(t, "", "")
	defer cleanup()

	bamPath := testutil.GetFilePath("@grailgo//bio/encoding/bam/testdata/170614_WGS_LOD_Pre_Library_B3_27961B_05.merged.10000.bam")
	tester := newRandomTester(t, bamPath, randomSeed)
	bamProvider := bamprovider.NewProvider(bamPath)

	pamPath := filepath.Join(tmpDir, "test.pam")
	require.NoError(t, converter.ConvertToPAM(pam.WriteOpts{}, pamPath, bamPath, "", math.MaxInt64))
	pamProvider := bamprovider.NewProvider(pamPath)

	for i := 0; i < 20; i++ {
		start, limit := tester.newRecAddr(), tester.newRecAddr()
		if start.EQ(limit) {
			continue
		}
		if limit.LT(start) {
			start, limit = limit, start
		}
		tester.testOnce(pamProvider, start, limit)

		// BAM doesn't support reading from a middle of a coord, so clear RecAddr.Seq.
		start.Seq = 0
		limit.Seq = 0
		if start.EQ(limit) {
			continue
		}
		tester.testOnce(bamProvider, start, limit)
	}
}

func TestRandom0(t *testing.T) { testRandom(t, 0) }
func TestRandom1(t *testing.T) { testRandom(t, 1) }
func TestRandom2(t *testing.T) { testRandom(t, 2) }
func TestRandom3(t *testing.T) { testRandom(t, 3) }

type randomTester struct {
	t         *testing.T
	r         *rand.Rand
	header    *sam.Header
	allRecs   []*sam.Record
	allCoords []biopb.Coord // coordinate of allRecs elems.
}

// Generate a random coordinate.
func (r *randomTester) newRecAddr() biopb.Coord {
	return r.allCoords[r.r.Int()%len(r.allCoords)]
}

// Return subset of r.allRecs that are in [start,limit).
func (r *randomTester) findRecordsInRange(start, limit biopb.Coord) []*sam.Record {
	recs := []*sam.Record{}
	gen := gbam.NewCoordGenerator()
	for _, r := range r.allRecs {
		addr := gen.GenerateFromRecord(r)
		if addr.GE(start) && addr.LT(limit) {
			recs = append(recs, r)
		}
	}
	return recs
}

// Read [start,limit) from provider, and check that the result matches the
// brute-forced list from findRecordsInRange.
func (r *randomTester) testOnce(provider bamprovider.Provider, start, limit biopb.Coord) {
	getRef := func(refID int32) *sam.Reference {
		if refID == -1 {
			return nil
		}
		return r.header.Refs()[refID]
	}

	vlog.Infof("Test: [%+v, %+v)", start, limit)
	shard := gbam.Shard{
		StartRef: getRef(start.RefId),
		Start:    int(start.Pos),
		StartSeq: int(start.Seq),
		EndRef:   getRef(limit.RefId),
		End:      int(limit.Pos),
		EndSeq:   int(limit.Seq),
	}
	iter := provider.NewIterator(shard)
	expected := r.findRecordsInRange(start, limit)
	n := 0
	for iter.Scan() {
		require.Equalf(r.t, expected[n].String(), iter.Record().String(),
			`Record %d mismatch, for range [%+v,%+v),
expected:%v
found:   %v`,
			n, start, limit, expected[n], iter.Record())
		n++
	}
	if len(expected) != n {
		vlog.Infof("Missed record, %v", expected[n])
		vlog.Infof("Missed record(2), %v", expected[n+1])
	}
	require.Equal(r.t, len(expected), n)
	require.NoError(r.t, iter.Close())
}

func newRandomTester(t *testing.T, bamPath string, randomSeed int64) *randomTester {
	tester := &randomTester{
		t: t,
		r: rand.New(rand.NewSource(randomSeed)),
	}
	in, err := os.Open(bamPath)
	require.NoError(t, err)
	bamr, err := bam.NewReader(in, 64)
	require.NoError(t, err)
	tester.header = bamr.Header()
	in, err = os.Open(bamPath)
	require.NoError(t, err)
	gen := gbam.NewCoordGenerator()
	for {
		rec, err := bamr.Read()
		if err == io.EOF {
			break
		}
		require.NoError(t, err)
		tester.allRecs = append(tester.allRecs, rec)
		tester.allCoords = append(tester.allCoords, gen.GenerateFromRecord(rec))
	}
	require.NoError(t, bamr.Close())
	require.NoError(t, in.Close())
	return tester
}

var (
	inputFlag = flag.String("input", "go/src/grail.com/bio/encoding/bam/testdata/170614_WGS_LOD_Pre_Library_B3_27961B_05.merged.10000.bam", "File to read in benchmark")
)

// BenchmarkSequentialRead reads the whole bam/pam file sequentially.
func BenchmarkSequentialRead(b *testing.B) {
	provider := bamprovider.NewProvider(*inputFlag)
	header, err := provider.GetHeader()
	require.NoError(b, err)
	for n := 0; n < b.N; n++ {
		shard := gbam.Shard{
			StartRef: header.Refs()[0],
			Start:    0,
			EndRef:   nil,
			End:      math.MaxInt32,
		}
		iter := provider.NewIterator(shard)
		nRecs := 0
		for iter.Scan() {
			record := iter.Record()
			nRecs++
			sam.PutInFreePool(record)
		}
		require.NoError(b, iter.Close())
		require.True(b, nRecs > 0)
	}
	require.NoError(b, provider.Close())
}

// Example of reading a BAM file in parallel.
func Example_shardedread() {
	path := testutil.GetFilePath("@grailgo//bio/encoding/bam/testdata/170614_WGS_LOD_Pre_Library_B3_27961B_05.merged.10000.bam")
	provider := bamprovider.NewProvider(path)
	shards, err := provider.GenerateShards(bamprovider.GenerateShardsOpts{})
	if err != nil {
		panic(err)
	}
	shardCh := gbam.NewShardChannel(shards)

	wg := sync.WaitGroup{}
	for i := 0; i < runtime.NumCPU(); i++ {
		wg.Add(1)
		go func() {
			defer wg.Done()
			for shard := range shardCh {
				iter := provider.NewIterator(shard)
				for iter.Scan() {
					// use iter.Record
				}
				if err := iter.Close(); iter != nil {
					panic(err)
				}
			}
		}()
	}
	wg.Wait()
	if err := provider.Close(); err != nil {
		panic(err)
	}
}

func ExampleParseFileType() {
	fmt.Printf("%d\n", bamprovider.ParseFileType("bam"))
	fmt.Printf("%d\n", bamprovider.ParseFileType("pam"))
	fmt.Printf("%d\n", bamprovider.ParseFileType("invalid"))
	// Output:
	// 1
	// 2
	// 0
}
