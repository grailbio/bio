package pamutil_test

import (
	"testing"

	"github.com/grailbio/bio/biopb"
	"github.com/grailbio/bio/encoding/pam/pamutil"
	"github.com/grailbio/testutil/expect"
)

func TestBlockRangeIntersects(t *testing.T) {
	tests := []struct {
		blockStartRefid, blockStartPos, blockEndRefid, blockEndPos,
		rangeStartRefid, rangeStartPos, rangeLimitRefid, rangeLimitPos int32
		intersect bool
	}{
		{3, 2, 10, 5, 10, 5, 11, 0, true},
		{3, 2, 10, 5, 10, 6, 11, 0, false},
		{3, 2, 10, 5, 3, 0, 3, 2, false},
		{3, 2, 10, 5, 3, 0, 3, 3, true},
	}
	for _, test := range tests {
		blockStart := biopb.Coord{test.blockStartRefid, test.blockStartPos, 0}
		blockEnd := biopb.Coord{test.blockEndRefid, test.blockEndPos, 0}
		r := biopb.CoordRange{biopb.Coord{test.rangeStartRefid, test.rangeStartPos, 0},
			biopb.Coord{test.rangeLimitRefid, test.rangeLimitPos, 0}}
		expect.EQ(t, pamutil.BlockIntersectsRange(blockStart, blockEnd, r), test.intersect, test)
	}
}
