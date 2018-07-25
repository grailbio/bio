package bam_test

import (
	"testing"

	"github.com/biogo/hts/sam"
	"github.com/grailbio/bio/encoding/bam"
	"github.com/grailbio/testutil"
	"github.com/grailbio/testutil/expect"
)

func TestShard(t *testing.T) {
	ref1, err := sam.NewReference("chr1", "", "", 100, nil, nil)
	expect.NoError(t, err)
	s := bam.Shard{StartRef: ref1, Start: 20, EndRef: ref1, End: 90, Padding: 3}
	expect.EQ(t, s.PaddedStart(), 17)
	expect.EQ(t, s.PaddedEnd(), 93)
	expect.EQ(t, s.PadStart(8), 12)
	expect.EQ(t, s.PadStart(21), 0)
	expect.EQ(t, s.PadEnd(11), 100)

	ref2, err := sam.NewReference("chr2", "", "", 200, nil, nil)
	expect.NoError(t, err)
	s = bam.Shard{StartRef: ref1, Start: 20, EndRef: ref2, End: 0, Padding: 3}
	expect.EQ(t, s.PaddedStart(), 17)
	expect.EQ(t, s.PaddedEnd(), 0)

	s = bam.Shard{StartRef: ref1, Start: 20, EndRef: ref2, End: 0, Padding: 3, EndSeq: 1}
	expect.EQ(t, s.PaddedStart(), 17)
	expect.EQ(t, s.PaddedEnd(), 3)
}

func TestNewShardChannel(t *testing.T) {
	ref1, err := sam.NewReference("chr1", "", "", 100, nil, nil)
	expect.NoError(t, err)
	ref2, err := sam.NewReference("chr2", "", "", 101, nil, nil)
	expect.NoError(t, err)
	ref3, err := sam.NewReference("chr3", "", "", 1, nil, nil)
	expect.NoError(t, err)
	header, _ := sam.NewHeader(nil, []*sam.Reference{ref1, ref2, ref3})
	shardList, err := bam.GetPositionBasedShards(header, 50, 10, false)
	expect.NoError(t, err)
	shardChan := bam.NewShardChannel(shardList)

	shards := []bam.Shard{}
	for s := range shardChan {
		shards = append(shards, s)
	}

	expect.EQ(t, len(shards), 6)
	expect.EQ(t, bam.Shard{ref1, ref1, 0, 50, 0, 0, 10, 0}, shards[0])
	expect.EQ(t, bam.Shard{ref1, ref1, 50, 100, 0, 0, 10, 1}, shards[1])
	expect.EQ(t, bam.Shard{ref2, ref2, 0, 50, 0, 0, 10, 2}, shards[2])
	expect.EQ(t, bam.Shard{ref2, ref2, 50, 100, 0, 0, 10, 3}, shards[3])
	expect.EQ(t, bam.Shard{ref2, ref2, 100, 101, 0, 0, 10, 4}, shards[4])
	expect.EQ(t, bam.Shard{ref3, ref3, 0, 1, 0, 0, 10, 5}, shards[5])
}

func TestGetByteBasedShards(t *testing.T) {
	bamPath := testutil.GetFilePath("//go/src/grail.com/bio/encoding/bam/testdata/170614_WGS_LOD_Pre_Library_B3_27961B_05.merged.10000.bam")
	baiPath := testutil.GetFilePath("//go/src/grail.com/bio/encoding/bam/testdata/170614_WGS_LOD_Pre_Library_B3_27961B_05.merged.10000.bam.bai")

	shardList, err := bam.GetByteBasedShards(bamPath, baiPath, 100000, 5000, 400, true)
	expect.Nil(t, err)

	for i, shard := range shardList {
		t.Logf("shard[%d]: %v", i, shard)
	}

	bam.ValidateShardList(shardList, 400)
}
