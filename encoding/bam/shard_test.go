package bam_test

import (
	"fmt"
	"testing"

	"github.com/grailbio/bio/biopb"
	"github.com/grailbio/bio/encoding/bam"
	"github.com/grailbio/bio/encoding/bamprovider"
	"github.com/grailbio/hts/sam"
	"github.com/grailbio/testutil"
	"github.com/grailbio/testutil/assert"
	"github.com/grailbio/testutil/expect"
	"github.com/grailbio/testutil/h"
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

	expect.That(t, shards, h.UnorderedElementsAre(
		bam.Shard{ref1, ref1, 0, 50, 0, 0, 10, 0},
		bam.Shard{ref1, ref1, 50, 100, 0, 0, 10, 1},
		bam.Shard{ref2, ref2, 0, 50, 0, 0, 10, 2},
		bam.Shard{ref2, ref2, 50, 100, 0, 0, 10, 3},
		bam.Shard{ref2, ref2, 100, 101, 0, 0, 10, 4},
		bam.Shard{ref3, ref3, 0, 1, 0, 0, 10, 5}))
}

func validateShards(t *testing.T, p bamprovider.Provider, shards []bam.Shard, includeUnmapped bool) {
	header, err := p.GetHeader()
	assert.NoError(t, err)
	unionShard := bam.UniversalShard(header)
	if !includeUnmapped {
		lastRef := header.Refs()[len(header.Refs())-1]
		unionShard.EndRef = lastRef
		unionShard.End = lastRef.Len()
	}
	iter0 := p.NewIterator(unionShard)
	for _, shard := range shards {
		t.Logf("Reading shard %+v", shard)
		n := 0
		iter1 := p.NewIterator(shard)
		shardRange := bam.ShardToCoordRange(shard)
		for iter1.Scan() {
			rec := iter1.Record()
			coord := biopb.Coord{RefId: int32(rec.Ref.ID()), Pos: int32(rec.Pos)}
			if rec.Ref != nil {
				assert.True(t, coord.GE(shardRange.Start), "shardrange=%+v, coord=%+v", shardRange, coord)
				assert.True(t, coord.LT(shardRange.Limit), "shardrange=%+v, coord=%+v rec=%s", shardRange, coord, rec.String())
			} else {
				// Unmapped reads should be in its own shard.
				assert.True(t, shard.StartRef == nil)
				assert.True(t, shard.EndRef == nil)
			}
			assert.True(t, iter0.Scan())
			assert.EQ(t, rec.String(), iter0.Record().String(), "n=%d", n)
			n++
		}
		iter1.Close()
	}
	assert.False(t, iter0.Scan())
	iter0.Close()
}

func TestGenerateShards(t *testing.T) {
	bamPath := testutil.GetFilePath("//go/src/grail.com/bio/encoding/bam/testdata/170614_WGS_LOD_Pre_Library_B3_27961B_05.merged.10000.bam")
	baiPath := testutil.GetFilePath("//go/src/grail.com/bio/encoding/bam/testdata/170614_WGS_LOD_Pre_Library_B3_27961B_05.merged.10000.bam.bai")
	p := bamprovider.NewProvider(bamPath, bamprovider.ProviderOpts{Index: baiPath})

	n := 0
	test := func(opts bamprovider.GenerateShardsOpts, cb func(*testing.T, []bam.Shard)) {
		t.Run(fmt.Sprint(n), func(t *testing.T) {
			t.Parallel()
			shardList, err := p.GenerateShards(opts)
			assert.NoError(t, err)
			cb(t, shardList)
		})
		n++
	}

	test(bamprovider.GenerateShardsOpts{
		Strategy:         bamprovider.ByteBased,
		IncludeUnmapped:  true,
		BytesPerShard:    50000,
		MinBasesPerShard: 5000,
	}, func(t *testing.T, shardList []bam.Shard) {
		assert.EQ(t, len(shardList), 51)
		validateShards(t, p, shardList, true)
	})

	test(bamprovider.GenerateShardsOpts{
		Strategy:        bamprovider.ByteBased,
		IncludeUnmapped: true,
		NumShards:       1,
	}, func(t *testing.T, shardList []bam.Shard) {
		assert.EQ(t, len(shardList), 2)
		validateShards(t, p, shardList, true)
	})
	test(bamprovider.GenerateShardsOpts{
		Strategy:        bamprovider.ByteBased,
		IncludeUnmapped: false,
		NumShards:       1,
	}, func(t *testing.T, shardList []bam.Shard) {
		assert.EQ(t, len(shardList), 1)
		validateShards(t, p, shardList, false)
	})
	assert.NoError(t, p.Close())
}
