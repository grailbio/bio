package bamprovider_test

import (
	"sort"
	"sync"
	"testing"

	"github.com/biogo/hts/sam"
	gbam "github.com/grailbio/bio/encoding/bam"
	"github.com/grailbio/bio/encoding/bamprovider"
	"github.com/grailbio/testutil"
	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
	"v.io/x/lib/vlog"
)

func newRecord(name string, ref *sam.Reference, pos int, mateRef *sam.Reference, matePos int, flags sam.Flags) *sam.Record {
	r := gbam.CastUp(gbam.GetFromFreePool())
	r.Name = name
	r.Ref = ref
	r.Pos = pos
	r.MateRef = mateRef
	r.MatePos = matePos
	r.Flags = flags
	return r
}

var (
	chr8, _          = sam.NewReference("chr8", "", "", 2000000, nil, nil)
	chr9, _          = sam.NewReference("chr9", "", "", 3000000, nil, nil)
	processHeader, _ = sam.NewHeader(nil, []*sam.Reference{chr8, chr9})
	read4            = newRecord("ABCDEFG", chr8, 123, nil, -1, sam.Read2)
	read5            = newRecord("ABCDEFG", chr8, 456, nil, -1, sam.Read1)
	read6            = newRecord("XYZ", chr8, 1024, nil, -1, sam.Read1)
	read7            = newRecord("foo", chr9, 777, nil, -1, sam.Read2)
	read8            = newRecord("foo", chr9, 1000001, nil, -1, sam.Read1)
	read9            = newRecord("XYZ", chr9, 2000000, nil, -1, sam.Read2)
	read9Secondary   = newRecord("XYZ", chr9, 2000002, nil, -1, sam.Read2|sam.Secondary)
	unmapped00       = newRecord("unmapped0", nil, -1, nil, -1, sam.Read1|sam.Unmapped|sam.MateUnmapped)
	unmapped01       = newRecord("unmapped0", nil, -1, nil, -1, sam.Read2|sam.Unmapped|sam.MateUnmapped)
	unmapped10       = newRecord("unmapped1", nil, -1, nil, -1, sam.Read1|sam.Unmapped|sam.MateUnmapped)
	unmapped11       = newRecord("unmapped1", nil, -1, nil, -1, sam.Read2|sam.Unmapped|sam.MateUnmapped)
	unmapped20       = newRecord("unmapped2", nil, -1, nil, -1, sam.Read1|sam.Unmapped|sam.MateUnmapped)
	unmapped21       = newRecord("unmapped2", nil, -1, nil, -1, sam.Read2|sam.Unmapped|sam.MateUnmapped)
)

type pair struct {
	r1 *sam.Record
	r2 *sam.Record
}

type pairByR1Pos []pair

func (a pairByR1Pos) Len() int      { return len(a) }
func (a pairByR1Pos) Swap(i, j int) { a[i], a[j] = a[j], a[i] }
func (a pairByR1Pos) Less(i, j int) bool {
	if a[i].r1.Ref.Name() != a[j].r1.Ref.Name() {
		return a[i].r1.Ref.Name() < a[j].r1.Ref.Name()
	}
	if a[i].r1.Pos != a[j].r1.Pos {
		return a[i].r1.Pos < a[j].r1.Pos
	}
	// for sorting unmapped records
	return a[i].r1.Flags < a[j].r1.Flags
}

func pairsEqualAnyOrder(t *testing.T, testName string, expected, actual []pair) {
	canonicalize := func(p *pair) {
		// Canonicalize the order of r1 and r2.
		if p.r2.Pos < p.r1.Pos || (p.r2.Pos == p.r1.Pos && p.r2.Flags < p.r1.Flags) {
			p.r1, p.r2 = p.r2, p.r1
		}
	}
	for i := range expected {
		canonicalize(&expected[i])
	}
	for i := range actual {
		canonicalize(&actual[i])
	}
	sort.Sort(pairByR1Pos(expected))
	sort.Sort(pairByR1Pos(actual))
	assert.Equal(t, expected, actual, "test %v", testName)
}

func TestGetPairs(t *testing.T) {
	tests := []struct {
		name          string
		records       []*sam.Record
		pairs         []pair
		unmapped      bool
		expectedPairs []pair
	}{
		{
			"basic",
			[]*sam.Record{read4, read5, read6, read7, read8, read9, read9Secondary},
			[]pair{pair{read5, read4}, pair{read6, read9}, pair{read8, read7}},
			false,
			nil,
		},
		{
			"basic2",
			[]*sam.Record{read4, read5, read6, read7, read8, read9, read9Secondary},
			[]pair{pair{read5, read4}, pair{read6, read9}, pair{read8, read7}},
			false,
			nil,
		},
		{
			"has unmapped reads, but does not read them",
			[]*sam.Record{read4, read5, read6, read7, read8, read9, read9Secondary, unmapped00, unmapped01},
			[]pair{pair{read5, read4}, pair{read6, read9}, pair{read8, read7}, pair{unmapped00, unmapped01}},
			false,
			[]pair{pair{read5, read4}, pair{read6, read9}, pair{read8, read7}},
		},
		{
			"has unmapped reads, and reads them",
			[]*sam.Record{read4, read5, read6, read7, read8, read9, read9Secondary, unmapped00, unmapped01},
			[]pair{pair{read5, read4}, pair{read6, read9}, pair{read8, read7}, pair{unmapped00, unmapped01}},
			true,
			nil,
		},
		{
			">1 unmapped reads across shards",
			[]*sam.Record{unmapped00, unmapped10, unmapped11, unmapped20, unmapped21, unmapped01},
			[]pair{pair{unmapped10, unmapped11}, pair{unmapped20, unmapped21}, pair{unmapped00, unmapped01}},
			true,
			nil,
		},
	}
	for i, test := range tests {
		vlog.Infof("Start test %v %+v", i, test)
		provider := bamprovider.NewFakeProvider(processHeader, test.records)
		iters, err := bamprovider.NewPairIterators(provider, test.unmapped)
		require.NoError(t, err)

		var pairs []pair
		for _, iter := range iters {
			for iter.Scan() {
				pairOrError := iter.Record()
				require.NoError(t, pairOrError.Err)
				pairs = append(pairs, pair{pairOrError.R1, pairOrError.R2})
			}
		}

		var expected []pair
		if test.expectedPairs != nil {
			expected = test.expectedPairs
		} else {
			expected = test.pairs
		}
		pairsEqualAnyOrder(t, test.name, expected, pairs)
	}
}

// Example_pairiterators is an example of NewPairIterator
func ExampleNewPairIterators() {
	bamPath := testutil.GetFilePath("//go/src/grail.com/bio/encoding/bam/testdata/170614_WGS_LOD_Pre_Library_B3_27961B_05.merged.10000.bam")
	provider := bamprovider.NewProvider(bamPath)
	iters, err := bamprovider.NewPairIterators(provider, true)
	if err != nil {
		panic(err)
	}

	wg := sync.WaitGroup{}
	for _, iter := range iters {
		wg.Add(1)
		go func(iter *bamprovider.PairIterator) {
			defer wg.Done()
			for iter.Scan() {
				p := iter.Record()
				if p.Err != nil {
					panic(p.Err)
				}
				// use p.R1 and p.R2
			}
		}(iter)
	}
	wg.Wait()
	if err := bamprovider.FinishPairIterators(iters); err != nil {
		panic(err)
	}
}
