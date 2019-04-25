package bamprovider_test

import (
	"flag"
	"fmt"
	"io"
	"log"
	"math"
	"math/rand"
	"os"
	"path/filepath"
	"runtime"
	"sync"
	"testing"

	"github.com/aws/aws-sdk-go/aws/session"
	"github.com/grailbio/base/file"
	"github.com/grailbio/base/file/s3file"
	"github.com/grailbio/base/grail"
	"github.com/grailbio/bio/biopb"
	gbam "github.com/grailbio/bio/encoding/bam"
	"github.com/grailbio/bio/encoding/bamprovider"
	"github.com/grailbio/bio/encoding/converter"
	"github.com/grailbio/bio/encoding/pam"
	"github.com/grailbio/hts/bam"
	"github.com/grailbio/hts/sam"
	"github.com/grailbio/testutil"
	"github.com/grailbio/testutil/assert"
	"github.com/grailbio/testutil/expect"
	"github.com/grailbio/testutil/h"
)

func TestMain(m *testing.M) {
	shutdown := grail.Init()
	file.RegisterImplementation("s3", func() file.Implementation {
		return s3file.NewImplementation(s3file.NewDefaultProvider(session.Options{}), s3file.Options{})
	})
	status := m.Run()
	shutdown()
	os.Exit(status)
}

func readIterator(i bamprovider.Iterator) []string {
	var names []string
	for i.Scan() {
		names = append(names, i.Record().Name)
	}
	return names
}

func doRead(t *testing.T, path string) []string {
	p := bamprovider.NewProvider(path)
	shards, err := p.GenerateShards(bamprovider.GenerateShardsOpts{
		IncludeUnmapped: true})
	assert.NoError(t, err)

	var names []string
	// Repeat the test to test iterator-reuse code path.
	for i := 0; i < 3; i++ {
		names = []string{}
		for _, shard := range shards {
			i := p.NewIterator(shard)
			names = append(names, readIterator(i)...)
			assert.NoError(t, i.Err())
			assert.NoError(t, i.Close())
		}
		assert.NoError(t, p.Close())
	}
	return names
}

func TestPAMSmall(t *testing.T) {
	tmpDir, cleanup := testutil.TempDir(t, "", "")
	defer cleanup()
	bamPath := testutil.GetFilePath("//go/src/grail.com/bio/encoding/bam/testdata/test-unmapped.bam")
	pamPath := filepath.Join(tmpDir, "test-unmapped.pam")
	assert.NoError(t, converter.ConvertToPAM(pam.WriteOpts{}, pamPath, bamPath, "", math.MaxInt64))
	assert.EQ(t, []string{"read1", "read2", "read3", "read10", "read10"}, doRead(t, pamPath))
}

func TestPAMLarge(t *testing.T) {
	tmpDir, cleanup := testutil.TempDir(t, "", "")
	defer cleanup()
	bamPath := testutil.GetFilePath("//go/src/grail.com/bio/encoding/bam/testdata/170614_WGS_LOD_Pre_Library_B3_27961B_05.merged.10000.bam")
	pamPath := filepath.Join(tmpDir, "large.pam")
	// The bam file is 2.8MB, so with 1MB shard size, we expect 3 shards.
	assert.NoError(t, converter.ConvertToPAM(pam.WriteOpts{}, pamPath, bamPath, "", 1<<20))

	p := bamprovider.NewProvider(pamPath)
	shards, err := p.GetFileShards()
	assert.NoError(t, err)
	header, err := p.GetHeader()
	assert.NoError(t, err)
	const n = 3
	assert.EQ(t, len(shards), n, "Shards:", shards)
	assert.EQ(t, header.Refs()[0], shards[0].StartRef, "Shards:", shards)
	assert.EQ(t, 0, shards[0].Start, "Shards:", shards)
	assert.True(t, shards[n-1].EndRef == nil, "Shards:", shards) // the last shard is always unmapped
	assert.EQ(t, math.MaxInt32, shards[n-1].End, "Shards:", shards)
}

func TestError(t *testing.T) {
	for _, test := range []struct{ path, errRe string }{
		{"nonexistent.bam", "no such file"},
		{"nonexistent.pam", "no pam file found"},
	} {
		p := bamprovider.NewProvider(test.path)
		_, err := p.GenerateShards(bamprovider.GenerateShardsOpts{IncludeUnmapped: true})
		assert.Regexp(t, err.Error(), test.errRe)

		iter := p.NewIterator(gbam.Shard{StartRef: nil, EndRef: nil, Start: 0, End: 1})
		assert.Regexp(t, iter.Close(), test.errRe)
		assert.Regexp(t, p.Close().Error(), test.errRe)
	}
}

func TestBAM(t *testing.T) {
	assert.That(t, doRead(t, testutil.GetFilePath("//go/src/grail.com/bio/encoding/bam/testdata/test.bam")),
		h.ElementsAre("read1", "read2", "read3"))
}

func TestBAMUnmapped(t *testing.T) {
	assert.That(t, doRead(t, testutil.GetFilePath("//go/src/grail.com/bio/encoding/bam/testdata/test-unmapped.bam")),
		h.ElementsAre("read1", "read2", "read3", "read10", "read10"))
}

func TestRefByName(t *testing.T) {
	p := bamprovider.NewProvider(testutil.GetFilePath("//go/src/grail.com/bio/encoding/bam/testdata/test.bam"))
	h, err := p.GetHeader()
	assert.NoError(t, err)
	assert.EQ(t, bamprovider.RefByName(h, "chr1").Name(), "chr1")
	assert.Nil(t, bamprovider.RefByName(h, "chr999"))
	assert.NoError(t, p.Close())
}

func TestNewRefIterator(t *testing.T) {
	p := bamprovider.NewProvider(testutil.GetFilePath("//go/src/grail.com/bio/encoding/bam/testdata/170614_WGS_LOD_Pre_Library_B3_27961B_05.merged.10000.bam"))

	// Read the first few reads in the file.
	iter := bamprovider.NewRefIterator(p, "chr1", 709300, 709306)
	expect.That(t, readIterator(iter), h.ElementsAre(
		"E00587:46:HK2FFALXX:1:1101:2798:35660:CCATTT+TAACGA"))
	assert.NoError(t, iter.Close())

	iter = bamprovider.NewRefIterator(p, "chr1", 709300, 709356)
	expect.That(t, readIterator(iter), h.ElementsAre(
		"E00587:46:HK2FFALXX:1:1101:2798:35660:CCATTT+TAACGA",
	))
	assert.NoError(t, iter.Close())

	iter = bamprovider.NewRefIterator(p, "chr1", 709300, 709357)
	expect.That(t, readIterator(iter), h.ElementsAre(
		"E00587:46:HK2FFALXX:1:1101:2798:35660:CCATTT+TAACGA",
		"E00587:46:HK2FFALXX:1:1101:2798:35660:CCATTT+TAACGA"))
	assert.NoError(t, iter.Close())

	assert.NoError(t, p.Close())
}

func TestBAMUnmappedOnly(t *testing.T) {
	assert.That(t, doRead(t, testutil.GetFilePath("//go/src/grail.com/bio/encoding/bam/testdata/test-unmapped-only.bam")),
		h.ElementsAre("read10", "read10"))
}

// Test reading random ranges.
func testRandom(t *testing.T, randomSeed int64) {
	tmpDir, cleanup := testutil.TempDir(t, "", "")
	defer cleanup()

	bamPath := testutil.GetFilePath("//go/src/grail.com/bio/encoding/bam/testdata/170614_WGS_LOD_Pre_Library_B3_27961B_05.merged.10000.bam")
	tester := newRandomTester(t, bamPath, randomSeed)
	bamProvider := bamprovider.NewProvider(bamPath)

	pamPath := filepath.Join(tmpDir, "test.pam")
	assert.NoError(t, converter.ConvertToPAM(pam.WriteOpts{}, pamPath, bamPath, "", math.MaxInt64))
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

	log.Printf("Test: [%+v, %+v)", start, limit)
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
		assert.EQ(r.t, iter.Record().String(), expected[n].String(),
			`Record %d mismatch, for range [%+v,%+v),
expected:%v
found:   %v`,
			n, start, limit, expected[n], iter.Record())

		n++
	}
	if len(expected) != n {
		log.Printf("Missed record, %v", expected[n])
		log.Printf("Missed record(2), %v", expected[n+1])
	}
	assert.EQ(r.t, n, len(expected))
	assert.NoError(r.t, iter.Close())
}

func newRandomTester(t *testing.T, bamPath string, randomSeed int64) *randomTester {
	tester := &randomTester{
		t: t,
		r: rand.New(rand.NewSource(randomSeed)),
	}
	in, err := os.Open(bamPath)
	assert.NoError(t, err)
	bamr, err := bam.NewReader(in, 64)
	assert.NoError(t, err)
	tester.header = bamr.Header()
	in, err = os.Open(bamPath)
	assert.NoError(t, err)
	gen := gbam.NewCoordGenerator()
	for {
		rec, err := bamr.Read()
		if err == io.EOF {
			break
		}
		assert.NoError(t, err)
		tester.allRecs = append(tester.allRecs, rec)
		tester.allCoords = append(tester.allCoords, gen.GenerateFromRecord(rec))
	}
	assert.NoError(t, bamr.Close())
	assert.NoError(t, in.Close())
	return tester
}

var (
	inputFlag = flag.String("input", "go/src/grail.com/bio/encoding/bam/testdata/170614_WGS_LOD_Pre_Library_B3_27961B_05.merged.10000.bam", "File to read in benchmark")
)

// BenchmarkSequentialRead reads the whole bam/pam file sequentially.
func BenchmarkSequentialRead(b *testing.B) {
	provider := bamprovider.NewProvider(*inputFlag)
	header, err := provider.GetHeader()
	assert.NoError(b, err)
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
		assert.NoError(b, iter.Close())
		assert.True(b, nRecs > 0)
	}
	assert.NoError(b, provider.Close())
}

// BenchmarkGenerateShards benchmarks the GenerateShards function.
func BenchmarkGenerateShards(b *testing.B) {
	for i := 0; i < b.N; i++ {
		provider := bamprovider.NewProvider(*inputFlag, bamprovider.ProviderOpts{Index: *inputFlag + ".bai"})
		shards, err := provider.GenerateShards(bamprovider.GenerateShardsOpts{
			NumShards: 16,
			Strategy:  bamprovider.ByteBased,
		})
		assert.NoError(b, err)
		assert.NoError(b, provider.Close())
		b.Logf("Created %d shards: %+v", len(shards), shards)
	}
}

// Example of reading a BAM file in parallel.
func Example_shardedread() {
	path := testutil.GetFilePath("//go/src/grail.com/bio/encoding/bam/testdata/170614_WGS_LOD_Pre_Library_B3_27961B_05.merged.10000.bam")
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

func getReadNames(t *testing.T, provider bamprovider.Provider) []string {
	opts := bamprovider.GenerateShardsOpts{
		Strategy:        bamprovider.ByteBased,
		IncludeUnmapped: true,
		BytesPerShard:   128 * 1024,
	}

	shards, err := provider.GenerateShards(opts)
	assert.NoError(t, err)
	names := make([]string, 0)
	for _, shard := range shards {
		iter := provider.NewIterator(shard)
		n := 0
		for iter.Scan() {
			r := iter.Record()
			names = append(names, r.Name)
			n++
		}
		assert.NoError(t, iter.Close())
	}
	return names
}

func TestGIndex(t *testing.T) {
	// Test that sharded reading using both .bai and .gbai return the
	// same set of read names.
	tmpDir, cleanup := testutil.TempDir(t, "", "")
	defer cleanup()
	bamPath := testutil.GetFilePath("//go/src/grail.com/bio/encoding/bam/testdata/170614_WGS_LOD_Pre_Library_B3_27961B_05.merged.10000.bam")
	baiPath := bamPath + ".bai"

	// Generate .gbai index file.
	r, err := os.Open(bamPath)
	assert.NoError(t, err)
	gbaiPath := filepath.Join(tmpDir, "foo.gbai")
	w, err := os.Create(gbaiPath)
	assert.NoError(t, err)
	err = gbam.WriteGIndex(w, r, 1024, 4)
	assert.NoError(t, err)
	assert.NoError(t, w.Close())
	assert.NoError(t, r.Close())

	// Get all the read names using .bai and .gbai, then compare the
	// read names.
	provider := bamprovider.NewProvider(bamPath, bamprovider.ProviderOpts{Index: baiPath})
	expected := getReadNames(t, provider)
	provider = bamprovider.NewProvider(bamPath, bamprovider.ProviderOpts{Index: gbaiPath})
	actual := getReadNames(t, provider)

	assert.EQ(t, len(expected), len(actual))
	for i := range expected {
		assert.EQ(t, expected[i], actual[i])
	}
}
