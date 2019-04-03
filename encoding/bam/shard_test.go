package bam

import (
	"fmt"
	"math"
	"os"
	"path/filepath"
	"reflect"
	"testing"

	biogobam "github.com/grailbio/hts/bam"
	"github.com/grailbio/hts/sam"
	"github.com/grailbio/testutil"
	"github.com/grailbio/testutil/expect"
	"github.com/grailbio/testutil/h"
	"github.com/stretchr/testify/require"
)

func TestShard(t *testing.T) {
	ref1, err := sam.NewReference("chr1", "", "", 100, nil, nil)
	expect.NoError(t, err)
	s := Shard{StartRef: ref1, Start: 20, EndRef: ref1, End: 90, Padding: 3}
	expect.EQ(t, s.PaddedStart(), 17)
	expect.EQ(t, s.PaddedEnd(), 93)
	expect.EQ(t, s.PadStart(8), 12)
	expect.EQ(t, s.PadStart(21), 0)
	expect.EQ(t, s.PadEnd(11), 100)

	ref2, err := sam.NewReference("chr2", "", "", 200, nil, nil)
	expect.NoError(t, err)
	s = Shard{StartRef: ref1, Start: 20, EndRef: ref2, End: 0, Padding: 3}
	expect.EQ(t, s.PaddedStart(), 17)
	expect.EQ(t, s.PaddedEnd(), 0)

	s = Shard{StartRef: ref1, Start: 20, EndRef: ref2, End: 0, Padding: 3, EndSeq: 1}
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
	shardList, err := GetPositionBasedShards(header, 50, 10, false)
	expect.NoError(t, err)
	shardChan := NewShardChannel(shardList)

	shards := []Shard{}
	for s := range shardChan {
		shards = append(shards, s)
	}

	expect.That(t, shards, h.UnorderedElementsAre(
		Shard{ref1, ref1, 0, 50, 0, 0, 10, 0},
		Shard{ref1, ref1, 50, 100, 0, 0, 10, 1},
		Shard{ref2, ref2, 0, 50, 0, 0, 10, 2},
		Shard{ref2, ref2, 50, 100, 0, 0, 10, 3},
		Shard{ref2, ref2, 100, 101, 0, 0, 10, 4},
		Shard{ref3, ref3, 0, 1, 0, 0, 10, 5}))
}

func TestGetByteBasedShards(t *testing.T) {
	bamPath := testutil.GetFilePath("//go/src/grail.com/bio/encoding/bam/testdata/170614_WGS_LOD_Pre_Library_B3_27961B_05.merged.10000.bam")
	baiPath := testutil.GetFilePath("//go/src/grail.com/bio/encoding/bam/testdata/170614_WGS_LOD_Pre_Library_B3_27961B_05.merged.10000.bam.bai")

	shardList, err := GetByteBasedShards(bamPath, baiPath, 100000, 5000, 400, true)
	expect.Nil(t, err)

	for i, shard := range shardList {
		t.Logf("shard[%d]: %v", i, shard)
	}

	bamIn, err := os.Open(bamPath)
	expect.NoError(t, err)
	defer func() {
		require.NoError(t, bamIn.Close())
	}()
	bamr, err := biogobam.NewReader(bamIn, 1)
	expect.NoError(t, err)
	header := bamr.Header()
	ValidateShardList(header, shardList, 400)
}

func TestGIndexShards(t *testing.T) {
	ref1, err := sam.NewReference("chr1", "", "", 100, nil, nil)
	expect.NoError(t, err)
	ref2, err := sam.NewReference("chr2", "", "", 101, nil, nil)
	expect.NoError(t, err)
	ref3, err := sam.NewReference("chr3", "", "", 102, nil, nil)
	expect.NoError(t, err)
	header1, _ := sam.NewHeader(nil, []*sam.Reference{ref1})
	header2, _ := sam.NewHeader(nil, []*sam.Reference{ref2, ref3})

	tests := []struct {
		header          *sam.Header
		padding         int
		includeUnmapped bool
		indexEntries    []GIndexEntry
		expectedShards  []Shard
	}{
		{ // Test case 0, one index entry per ref, one ref.  Excludes unmapped.
			header:          header1,
			padding:         44,
			includeUnmapped: false,
			indexEntries: []GIndexEntry{
				GIndexEntry{
					RefID:   0,
					Pos:     0,
					Seq:     0,
					VOffset: 1<<16 | 2,
				},
			},
			expectedShards: []Shard{
				Shard{
					StartRef: ref1,
					EndRef:   ref1,
					Start:    0,
					End:      100,
					Padding:  44,
					ShardIdx: 0,
				},
			},
		},
		{ // Test case 1, one index entry per ref, two refs.  Excludes unmapped.
			header:          header2,
			padding:         44,
			includeUnmapped: false,
			indexEntries: []GIndexEntry{
				GIndexEntry{
					RefID:   0,
					Pos:     0,
					Seq:     0,
					VOffset: 1<<16 | 2,
				},
				GIndexEntry{
					RefID:   1,
					Pos:     0,
					Seq:     0,
					VOffset: 1<<16 | 12345,
				},
			},
			expectedShards: []Shard{
				Shard{
					StartRef: ref2,
					EndRef:   ref2,
					Start:    0,
					End:      101,
					Padding:  44,
					ShardIdx: 0,
				},
				Shard{
					StartRef: ref3,
					EndRef:   ref3,
					Start:    0,
					End:      102,
					Padding:  44,
					ShardIdx: 1,
				},
			},
		},
		{ // Test case 2, one index entry per ref.  Includes unmapped
			// even though there is no index entry for unmapped reads.
			header:          header2,
			padding:         44,
			includeUnmapped: true,
			indexEntries: []GIndexEntry{
				GIndexEntry{
					RefID:   0,
					Pos:     0,
					Seq:     0,
					VOffset: 1<<16 | 2,
				},
				GIndexEntry{
					RefID:   1,
					Pos:     0,
					Seq:     0,
					VOffset: 1<<16 | 12345,
				},
			},
			expectedShards: []Shard{
				Shard{
					StartRef: ref2,
					EndRef:   ref2,
					Start:    0,
					End:      101,
					Padding:  44,
					ShardIdx: 0,
				},
				Shard{
					StartRef: ref3,
					EndRef:   ref3,
					Start:    0,
					End:      102,
					Padding:  44,
					ShardIdx: 1,
				},
				Shard{
					StartRef: nil,
					EndRef:   nil,
					Start:    0,
					End:      math.MaxInt32,
					Padding:  0,
					ShardIdx: 2,
				},
			},
		},
		{ // Test case 3, two index entries, far enough apart to
			// create a new shard.
			header:          header1,
			padding:         44,
			includeUnmapped: false,
			indexEntries: []GIndexEntry{
				GIndexEntry{
					RefID:   0,
					Pos:     0,
					Seq:     0,
					VOffset: 1<<16 | 2,
				},
				GIndexEntry{
					RefID:   0,
					Pos:     50,
					Seq:     0,
					VOffset: 1001<<16 | 12345,
				},
			},
			expectedShards: []Shard{
				Shard{
					StartRef: ref1,
					EndRef:   ref1,
					Start:    0,
					End:      50,
					Padding:  44,
					ShardIdx: 0,
				},
				Shard{
					StartRef: ref1,
					EndRef:   ref1,
					Start:    50,
					End:      100,
					Padding:  44,
					ShardIdx: 1,
				},
			},
		},
		{ // Test case 4, two index entries, too close to create a new
			// shard.
			header:          header1,
			padding:         44,
			includeUnmapped: false,
			indexEntries: []GIndexEntry{
				GIndexEntry{
					RefID:   0,
					Pos:     0,
					Seq:     0,
					VOffset: 1<<16 | 2,
				},
				GIndexEntry{
					RefID:   0,
					Pos:     50,
					Seq:     0,
					VOffset: 1000<<16 | 12345,
				},
			},
			expectedShards: []Shard{
				Shard{
					StartRef: ref1,
					EndRef:   ref1,
					Start:    0,
					End:      100,
					Padding:  44,
					ShardIdx: 0,
				},
			},
		},
		{ // Test case 5, two index entries, far enough apart to
			// create a new shard, but the genomic positions are too
			// close to exceed min_bases.
			header:          header1,
			padding:         44,
			includeUnmapped: false,
			indexEntries: []GIndexEntry{
				GIndexEntry{
					RefID:   0,
					Pos:     0,
					Seq:     0,
					VOffset: 1<<16 | 2,
				},
				GIndexEntry{
					RefID:   0,
					Pos:     20,
					Seq:     0,
					VOffset: 1001<<16 | 12345,
				},
			},
			expectedShards: []Shard{
				Shard{
					StartRef: ref1,
					EndRef:   ref1,
					Start:    0,
					End:      100,
					Padding:  44,
					ShardIdx: 0,
				},
			},
		},
		{ // Test case 6, one index entry, that points to pos 1.
			// Shard should start at position 0.
			header:          header1,
			padding:         44,
			includeUnmapped: false,
			indexEntries: []GIndexEntry{
				GIndexEntry{
					RefID:   0,
					Pos:     1,
					Seq:     0,
					VOffset: 1<<16 | 2,
				},
			},
			expectedShards: []Shard{
				Shard{
					StartRef: ref1,
					EndRef:   ref1,
					Start:    0,
					End:      100,
					Padding:  44,
					ShardIdx: 0,
				},
			},
		},
	}

	for testIdx, test := range tests {
		t.Run(fmt.Sprintf("test case %d", testIdx), func(t *testing.T) {
			tempDir, cleanup := testutil.TempDir(t, "", "")
			defer cleanup()

			// Write out a gindex file using indexEntries.
			gbaiPath := filepath.Join(tempDir, "foo.gbai")
			f, err := os.Create(gbaiPath)
			require.NoError(t, err)
			writer := newGIndexWriter(f)
			writer.writeHeader()
			for _, e := range test.indexEntries {
				t.Logf("appending %v", e)
				require.NoError(t, writer.append(&e))
			}
			require.NoError(t, writer.close())
			require.NoError(t, f.Close())

			// Calculate shards
			actual, err := gbaiByteBasedShards(test.header, gbaiPath, 1000, 25,
				test.padding, test.includeUnmapped)
			require.NoError(t, err)
			expect.True(t, reflect.DeepEqual(test.expectedShards, actual), "%v --- %v",
				test.expectedShards, actual)
		})
	}
}
