package bam

import (
	"bytes"
	"context"
	"fmt"
	"io"
	"math"
	"os"
	"path/filepath"
	"testing"

	biogobam "github.com/grailbio/hts/bam"
	"github.com/grailbio/hts/bgzf"
	"github.com/grailbio/hts/sam"
	"github.com/grailbio/testutil"
	"github.com/grailbio/testutil/expect"
	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

const (
	offset0 = 1<<16 | 2
	offset1 = 3<<16 | 4
	offset2 = 5<<16 | 6
	offset3 = 7<<16 | 8
)

func bgzfOffset(offset uint64) bgzf.Offset {
	return bgzf.Offset{int64(offset >> 16), uint16(offset & 0xffff)}
}

func TestRecordOffset(t *testing.T) {
	index := make(GIndex, 0)
	index = append(index, GIndexEntry{0, 0, 1, offset0})
	index = append(index, GIndexEntry{0, 5, 0, offset1})
	index = append(index, GIndexEntry{1, 0, 3, offset2})
	index = append(index, GIndexEntry{-1, 0, 0, offset3})

	offset := index.RecordOffset(0, 0, 1)
	assert.Equal(t, bgzfOffset(offset0), offset)

	offset = index.RecordOffset(0, 6, 0)
	assert.Equal(t, bgzfOffset(offset1), offset)

	offset = index.RecordOffset(1, 0, 0)
	assert.Equal(t, bgzfOffset(offset1), offset)

	offset = index.RecordOffset(1, 0, 3)
	assert.Equal(t, bgzfOffset(offset2), offset)

	offset = index.UnmappedOffset()
	assert.Equal(t, bgzfOffset(offset3), offset)
}

func TestUnmapped(t *testing.T) {
	index := GIndex{GIndexEntry{-1, 0, 0, offset0}}
	offset := index.UnmappedOffset()
	assert.Equal(t, bgzfOffset(offset0), offset)
}

func TestWriteGIndex(t *testing.T) {
	bamPath := testutil.GetFilePath("//go/src/grail.com/bio/encoding/bam/testdata/170614_WGS_LOD_Pre_Library_B3_27961B_05.merged.10000.bam")
	f, err := os.Open(bamPath)
	require.NoError(t, err)
	defer func() {
		require.NoError(t, f.Close())
	}()

	// Write a .gbai index and read it back.
	var indexBuf bytes.Buffer
	require.NoError(t, WriteGIndex(&indexBuf, f, 1024, 1))
	index, err := ReadGIndex(&indexBuf)
	require.NoError(t, err)

	// Open the bam file
	f2, err := os.Open(bamPath)
	require.NoError(t, err)
	defer func() {
		require.NoError(t, f2.Close())
	}()
	reader, err := biogobam.NewReader(f2, 1)
	require.NoError(t, err)

	bamPosCounts := make(map[string]int)
	indexedPosCounts := make(map[string]int)

	// Seek to each entry in the index, and check that the record's
	// (ref, pos) is equal to the index entry's (ref, pos).
	for _, e := range *index {
		err := reader.Seek(ToBGZFOffset(e.VOffset))
		require.NoError(t, err)
		record, err := reader.Read()
		require.NoError(t, err)
		assert.Equal(t, int(e.RefID), record.Ref.ID())
		assert.Equal(t, int(e.Pos), record.Pos)

		// For now, we expect the Seq number to always be zero.
		assert.Equal(t, int(e.Seq), 0)

		// Count how many records have the same (ref, pos).
		indexedPosCounts[fmt.Sprintf("%d:%d", record.Ref.ID(), record.Pos)]++
		for {
			record, err := reader.Read()
			if err == io.EOF {
				break
			}
			assert.NoError(t, err)
			if int(e.RefID) != record.Ref.ID() || int(e.Pos) != record.Pos {
				break
			}
			indexedPosCounts[fmt.Sprintf("%d:%d", record.Ref.ID(), record.Pos)]++
		}
	}
	assert.NoError(t, reader.Close())

	// Open the bam file again, to read all position counts.  Check
	// that the counts match the bam.  This confirms that since e.Seq
	// is zero, then the number of times the above reader sees a
	// particular (ref, position) pair is equal to the number of
	// records in the bam at (ref, position).
	f3, err := os.Open(bamPath)
	require.NoError(t, err)
	defer func() {
		require.NoError(t, f3.Close())
	}()
	reader, err = biogobam.NewReader(f3, 1)
	require.NoError(t, err)

	for {
		record, err := reader.Read()
		if err == io.EOF {
			break
		}
		assert.NoError(t, err)
		bamPosCounts[fmt.Sprintf("%d:%d", record.Ref.ID(), record.Pos)]++
	}

	for k, v := range indexedPosCounts {
		assert.Equal(t, bamPosCounts[k], v, "key: %s", k)
	}
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
					EndRef:   ref3,
					Start:    0,
					End:      102,
					Padding:  44,
					ShardIdx: 0,
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
					EndRef:   ref3,
					Start:    0,
					End:      102,
					Padding:  44,
					ShardIdx: 0,
				},
				Shard{
					StartRef: nil,
					EndRef:   nil,
					Start:    0,
					End:      math.MaxInt32,
					ShardIdx: 1,
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
		t.Run(fmt.Sprint(testIdx), func(t *testing.T) {
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
			actual, err := gbaiByteBasedShards(context.Background(), test.header, gbaiPath, 1000, 25,
				test.padding, test.includeUnmapped)
			require.NoError(t, err)
			expect.EQ(t, actual, test.expectedShards, "%v --- %v",
				actual, test.expectedShards)
		})
	}
}
