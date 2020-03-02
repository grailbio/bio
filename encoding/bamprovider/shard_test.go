package bamprovider_test

import (
	"fmt"
	"testing"

	"github.com/grailbio/bio/biopb"
	"github.com/grailbio/bio/encoding/bam"
	"github.com/grailbio/bio/encoding/bamprovider"
	"github.com/grailbio/testutil"
	"github.com/grailbio/testutil/assert"
)

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

// TODO(josh): When this test is located in github.com/grailbio/bio/encoding/bam with the code
// it exercises, the Bazel go_default_test build fails with a package height error that may be
// similar to https://github.com/bazelbuild/rules_go/issues/1877. Consider moving this back to
// package bam when that's fixed.
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
