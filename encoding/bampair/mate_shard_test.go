package bampair

import (
	"testing"

	"github.com/grailbio/hts/sam"
	"github.com/grailbio/testutil"
	"github.com/stretchr/testify/assert"
	"grail.com/bio/encoding/bam"
)

func createShardInfo() *ShardInfo {
	shardInfo := newShardInfo()
	shard1 := &bam.Shard{
		StartRef: chr1,
		EndRef:   chr1,
		Start:    0,
		End:      5,
		ShardIdx: 0,
	}
	shard2 := &bam.Shard{
		StartRef: chr1,
		EndRef:   nil,
		Start:    5,
		End:      100,
		ShardIdx: 1,
	}
	shardInfo.add(shard1)
	shardInfo.add(shard2)
	shardInfo.updateInfoByShard(shard1, 0, 1000)
	shardInfo.updateInfoByShard(shard2, 0, 2000)
	shardInfo.computeFileIndexes()
	return shardInfo
}

func TestMemMateShard(t *testing.T) {
	r1 := &sam.Record{Name: "foo", Ref: chr1, Pos: 0, MateRef: chr1, MatePos: 10}
	r2 := &sam.Record{Name: "foo", Ref: chr1, Pos: 10, MateRef: chr1, MatePos: 0}
	shardInfo := createShardInfo()

	shard := newMemMateShard()
	shard.add(r1, 123)
	shard.add(r2, 456)
	shard.closeWriter()

	// Check that mate is correct, and fileIdx is computed correctly.
	shard.openReader()
	mate, fileIdx := shard.getMate(shardInfo, r1)
	assert.True(t, recordsEqual(r2, mate))
	assert.Equal(t, uint64(1456), fileIdx)
}

func TestDiskMateShard(t *testing.T) {
	r1 := &sam.Record{Name: "foo", Ref: chr1, Pos: 0, MateRef: chr1, MatePos: 10}
	r2 := &sam.Record{Name: "foo", Ref: chr1, Pos: 10, MateRef: chr1, MatePos: 0}
	shardInfo := createShardInfo()

	tempDir, cleanup := testutil.TempDir(t, "", "")
	defer cleanup()

	shard, err := newDiskMateShard(header, tempDir, 0, 1000)
	assert.NoError(t, err)
	shard.add(r1, 123)
	shard.add(r2, 456)
	shard.closeWriter()

	// Check that mate is correct, and fileIdx is computed correctly.
	shard.openReader()
	mate, fileIdx := shard.getMate(shardInfo, r1)
	assert.True(t, recordsEqual(r2, mate))
	assert.Equal(t, uint64(1456), fileIdx)
}
