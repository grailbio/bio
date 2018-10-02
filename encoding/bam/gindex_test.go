package bam

import (
	"bytes"
	"fmt"
	"io"
	"os"
	"testing"

	"github.com/biogo/hts/bam"
	"github.com/biogo/hts/bgzf"
	"github.com/grailbio/testutil"
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
	reader, err := bam.NewReader(f2, 1)
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
	reader, err = bam.NewReader(f3, 1)
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
