package sorter

import (
	"fmt"
	"testing"

	"github.com/grailbio/bio/biopb"
	"github.com/stretchr/testify/assert"
)

func block(key, seq, offset, numRecords int) biopb.SortShardBlockIndex {
	return biopb.SortShardBlockIndex{StartKey: uint64(key), StartSeq: uint64(seq), FileOffset: uint64(offset), NumRecords: uint32(numRecords)}
}

func keysString(coords []recCoord) string {
	s := ""
	for _, k := range coords {
		if s != "" {
			s += ","
		}
		s += fmt.Sprintf("%d", k)
	}
	return s
}
func TestShardBounds(t *testing.T) {
	var blocks = []biopb.SortShardBlockIndex{
		block(1, 0, 10, 100),
		block(2, 0, 11, 100),
		block(2, 0, 12, 100),
		block(2, 0, 12, 100),
		block(5, 0, 13, 100),
		block(6, 0, 13, 100),
		block(8, 0, 13, 100),
	}
	assert.Equal(t, "", keysString(computePAMShardBounds([][]biopb.SortShardBlockIndex{blocks}, 10000)))
	assert.Equal(t, "2,4,8", keysString(computePAMShardBounds([][]biopb.SortShardBlockIndex{blocks}, 200)))
	assert.Equal(t, "2,4,6,8", keysString(computePAMShardBounds([][]biopb.SortShardBlockIndex{blocks}, 50)))
}

func TestShardFileOffset(t *testing.T) {
	var blocks = []biopb.SortShardBlockIndex{
		block(1, 10, 10, 99),
		block(2, 0, 11, 99),
		block(2, 0, 12, 99),
		block(2, 10, 13, 99),
		block(8, 0, 14, 99),
	}

	assert.Equal(t, int64(0), startFileOffset(blocks, 0))
	assert.Equal(t, int64(0), startFileOffset(blocks, 1))
	assert.Equal(t, int64(11), startFileOffset(blocks, 2))
	assert.Equal(t, int64(13), startFileOffset(blocks, 5))
	assert.Equal(t, int64(14), startFileOffset(blocks, 8))
	assert.Equal(t, int64(14), startFileOffset(blocks, 1000))

	assert.Equal(t, int64(10), limitFileOffset(blocks, 0))
	assert.Equal(t, int64(10), limitFileOffset(blocks, 1))
	assert.Equal(t, int64(11), limitFileOffset(blocks, 2))
	assert.Equal(t, int64(14), limitFileOffset(blocks, 5))
	assert.Equal(t, int64(14), limitFileOffset(blocks, 8))
}
