package bam_test

import (
	"testing"

	"github.com/biogo/hts/sam"
	"github.com/grailbio/bio/encoding/bam"
	"github.com/grailbio/testutil"
	"github.com/stretchr/testify/assert"
)

func TestShard(t *testing.T) {
	ref1, err := sam.NewReference("chr1", "", "", 100, nil, nil)
	assert.NoError(t, err)
	s := bam.Shard{StartRef: ref1, Start: 20, EndRef: ref1, End: 90, Padding: 3}
	assert.Equal(t, 17, s.PaddedStart())
	assert.Equal(t, 93, s.PaddedEnd())
	assert.Equal(t, 12, s.PadStart(8))
	assert.Equal(t, 0, s.PadStart(21))
	assert.Equal(t, 100, s.PadEnd(11))

	ref2, err := sam.NewReference("chr2", "", "", 200, nil, nil)
	assert.NoError(t, err)
	s = bam.Shard{StartRef: ref1, Start: 20, EndRef: ref2, End: 0, Padding: 3}
	assert.Equal(t, 17, s.PaddedStart())
	assert.Equal(t, 0, s.PaddedEnd())

	s = bam.Shard{StartRef: ref1, Start: 20, EndRef: ref2, End: 0, Padding: 3, EndSeq: 1}
	assert.Equal(t, 17, s.PaddedStart())
	assert.Equal(t, 3, s.PaddedEnd())
}

func TestNewShardChannel(t *testing.T) {
	ref1, err := sam.NewReference("chr1", "", "", 100, nil, nil)
	assert.NoError(t, err)
	ref2, err := sam.NewReference("chr2", "", "", 101, nil, nil)
	assert.NoError(t, err)
	ref3, err := sam.NewReference("chr3", "", "", 1, nil, nil)
	assert.NoError(t, err)
	header, _ := sam.NewHeader(nil, []*sam.Reference{ref1, ref2, ref3})
	shardList, err := bam.GetPositionBasedShards(header, 50, 10, false)
	assert.NoError(t, err)
	shardChan := bam.NewShardChannel(shardList)

	shards := []bam.Shard{}
	for s := range shardChan {
		shards = append(shards, s)
	}

	assert.Equal(t, 6, len(shards))
	assert.Equal(t, shards[0], bam.Shard{ref1, ref1, 0, 50, 0, 0, 10, 0})
	assert.Equal(t, shards[1], bam.Shard{ref1, ref1, 50, 100, 0, 0, 10, 1})
	assert.Equal(t, shards[2], bam.Shard{ref2, ref2, 0, 50, 0, 0, 10, 2})
	assert.Equal(t, shards[3], bam.Shard{ref2, ref2, 50, 100, 0, 0, 10, 3})
	assert.Equal(t, shards[4], bam.Shard{ref2, ref2, 100, 101, 0, 0, 10, 4})
	assert.Equal(t, shards[5], bam.Shard{ref3, ref3, 0, 1, 0, 0, 10, 5})
}

func TestGetByteBasedShards(t *testing.T) {
	bamPath := testutil.GetFilePath("//go/src/grail.com/bio/encoding/bam/testdata/170614_WGS_LOD_Pre_Library_B3_27961B_05.merged.10000.bam")
	baiPath := testutil.GetFilePath("//go/src/grail.com/bio/encoding/bam/testdata/170614_WGS_LOD_Pre_Library_B3_27961B_05.merged.10000.bam.bai")

	shardList, err := bam.GetByteBasedShards(bamPath, baiPath, 100000, 5000, 400, true)
	assert.Nil(t, err)

	for i, shard := range shardList {
		t.Logf("shard[%d]: %v", i, shard)
	}

	bam.ValidateShardList(shardList, 400)
}
