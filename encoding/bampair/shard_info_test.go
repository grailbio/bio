package bampair

import (
	"testing"

	"github.com/grailbio/hts/sam"
	"github.com/stretchr/testify/assert"
	grailbam "grail.com/bio/encoding/bam"
)

func newShard(ref *sam.Reference, start, end, padding, index int) *grailbam.Shard {
	return &grailbam.Shard{ref, ref, start, end, 0, 0, padding, index}
}

func TestShardInfo(t *testing.T) {
	chr1, _ := sam.NewReference("chr1", "", "", 1000, nil, nil)
	chr2, _ := sam.NewReference("chr2", "", "", 2000, nil, nil)

	shardInfo := newShardInfo()
	shardInfo.add(newShard(chr1, 0, 100, 10, 0))
	shardInfo.add(newShard(chr1, 100, 200, 10, 1))
	shardInfo.add(newShard(chr1, 200, 300, 10, 2))
	shardInfo.add(newShard(chr2, 0, 100, 0, 3))

	checkInfoByRecord := func(r *sam.Record, expectedRef *sam.Reference, expectedShardStart int) {
		info := shardInfo.getInfoByRecord(r)
		assert.Equal(t, expectedRef.ID(), info.Shard.StartRef.ID())
		assert.Equal(t, expectedShardStart, info.Shard.Start)
	}
	checkInfoByRecord(newRecord("foo", chr1, 0, 0, 0, nil, nil), chr1, 0)
	checkInfoByRecord(newRecord("foo", chr1, 99, 0, 0, nil, nil), chr1, 0)
	checkInfoByRecord(newRecord("foo", chr1, 100, 0, 0, nil, nil), chr1, 100)
	checkInfoByRecord(newRecord("foo", chr1, 101, 0, 0, nil, nil), chr1, 100)
	checkInfoByRecord(newRecord("foo", chr2, 0, 0, 0, nil, nil), chr2, 0)

	checkInfoByShard := func(s *grailbam.Shard, expectedRef *sam.Reference, expectedShardStart int) {
		info := shardInfo.GetInfoByShard(s)
		assert.Equal(t, expectedRef.ID(), info.Shard.StartRef.ID())
		assert.Equal(t, expectedShardStart, info.Shard.Start)
	}
	checkInfoByShard(newShard(chr1, 0, 100, 10, 0), chr1, 0)
	checkInfoByShard(newShard(chr1, 100, 200, 10, 2), chr1, 100)

	checkMateShard := func(r *sam.Record, expectedRef *sam.Reference, expectedShardStart int) {
		shard := shardInfo.getMateShard(r)
		assert.Equal(t, expectedRef.ID(), shard.StartRef.ID())
		assert.Equal(t, expectedShardStart, shard.Start)
	}
	checkMateShard(newRecord("foo", nil, 0, 0, 0, chr1, nil), chr1, 0)
	checkMateShard(newRecord("foo", nil, 0, 0, 99, chr1, nil), chr1, 0)
	checkMateShard(newRecord("foo", nil, 0, 0, 100, chr1, nil), chr1, 100)
	checkMateShard(newRecord("foo", nil, 0, 0, 101, chr1, nil), chr1, 100)
	checkMateShard(newRecord("foo", nil, 0, 0, 0, chr2, nil), chr2, 0)
}
