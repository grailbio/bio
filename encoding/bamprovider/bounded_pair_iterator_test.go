package bamprovider_test

import (
	"testing"

	"github.com/grailbio/bio/encoding/bamprovider"
	"github.com/grailbio/hts/sam"
	"github.com/grailbio/testutil/assert"
	"v.io/x/lib/vlog"
)

var (
	read10 = newRecord("dist999", chr9, 1000, chr9, 1999, sam.Read1)
	read11 = newRecord("dist999", chr9, 1999, chr9, 1000, sam.Read2)
	read12 = newRecord("dist1000", chr9, 2000, chr9, 3000, sam.Read1)
	read13 = newRecord("dist1000", chr9, 3000, chr9, 2000, sam.Read2)
	read14 = newRecord("dist1001", chr9, 3000, chr9, 4001, sam.Read1)
	read15 = newRecord("dist1001", chr9, 4001, chr9, 3000, sam.Read2)
)

func TestGetBoundedPairs(t *testing.T) {
	tests := []struct {
		name          string
		records       []*sam.Record
		pairs         []pair
		expectedPairs []pair
	}{
		{
			"basic",
			[]*sam.Record{read4, read5, read6, read7, read8, read9, read9Secondary},
			[]pair{pair{read5, read4}, pair{read6, read9}, pair{read8, read7}},
			[]pair{pair{read5, read4}},
		},
		{
			"has unmapped reads, but does not read them",
			[]*sam.Record{read4, read5, read6, read7, read8, read9, read9Secondary, unmapped00, unmapped01},
			[]pair{pair{read5, read4}, pair{read6, read9}, pair{read8, read7}, pair{unmapped00, unmapped01}},
			[]pair{pair{read5, read4}},
		},
		{
			"padding test 1",
			[]*sam.Record{read10, read11, read12, read13},
			[]pair{pair{read10, read11}, pair{read12, read13}},
			[]pair{pair{read10, read11}, pair{read12, read13}},
		},
		{
			"padding test 2",
			[]*sam.Record{read10, read11, read14, read15},
			[]pair{pair{read10, read11}, pair{read14, read15}},
			[]pair{pair{read10, read11}},
		},
	}
	for i, test := range tests {
		vlog.Infof("Start test %v %+v", i, test)
		provider := bamprovider.NewFakeProvider(processHeader, test.records)
		iters, err := bamprovider.NewBoundedPairIterators(provider, bamprovider.BoundedPairIteratorOpts{})
		assert.NoError(t, err)

		var pairs []pair
		for _, iter := range iters {
			for iter.Scan() {
				pairOrError := iter.Record()
				assert.NoError(t, pairOrError.Err)
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
