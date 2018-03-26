package bam

import (
	"bytes"
	"encoding/binary"
	"io"
	"reflect"
	"sort"
	"strconv"
	"strings"
	"testing"

	"github.com/biogo/hts/bgzf"
	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

func toInt(t *testing.T, s string) int {
	i, err := strconv.Atoi(s)
	require.Nil(t, err)
	return i
}

func writeBin(t *testing.T, w io.Writer, s string) {
	bins := strings.Split(s, ":")
	// Write the number of bins
	err := binary.Write(w, binary.LittleEndian, int32(len(bins)))
	assert.Nil(t, err)

	for _, bin := range bins {
		binContent := strings.Split(bin, ",")

		// Write the bin number
		err = binary.Write(w, binary.LittleEndian, uint32(toInt(t, binContent[0])))
		assert.Nil(t, err)
		binContent = binContent[1:]

		// Write the number of chunks
		err = binary.Write(w, binary.LittleEndian, int32(len(binContent)/2))
		assert.Nil(t, err)

		// Write the chunks
		for _, voffset := range binContent {
			err = binary.Write(w, binary.LittleEndian, uint64(toInt(t, voffset)))
			assert.Nil(t, err)
		}
	}
}

func writeIntervals(t *testing.T, w io.Writer, s string) {
	intervals := strings.Split(s, ",")

	// Write the number of intervals
	err := binary.Write(w, binary.LittleEndian, int32(len(intervals)))
	assert.Nil(t, err)

	for _, voffset := range intervals {
		err = binary.Write(w, binary.LittleEndian, uint64(toInt(t, voffset)))
		assert.Nil(t, err)
	}
}

func writeIndex(t *testing.T, bins, intervals []string, unmapped int) *bytes.Buffer {
	var buf bytes.Buffer
	magic := []byte{'B', 'A', 'I', 0x1}
	_, err := buf.Write(magic)
	assert.Nil(t, err)

	// Two references
	err = binary.Write(&buf, binary.LittleEndian, int32(len(bins)))
	assert.Nil(t, err)

	for i := range bins {
		writeBin(t, &buf, bins[i])
		writeIntervals(t, &buf, intervals[i])
	}

	// Write unmapped count
	if unmapped >= 0 {
		err = binary.Write(&buf, binary.LittleEndian, uint64(unmapped))
		assert.Nil(t, err)
	}

	return &buf
}

func TestReadIndex(t *testing.T) {
	tests := []struct {
		bins      []string
		intervals []string
		unmapped  int
	}{
		{
			bins: []string{
				"100,1,2:200,3,4:37450,5,6,7,8",
				"100,10,22",
				"37450,5,6,7,8",
				"200,100002,200003", // Use a voffset larger than 16 bits to check that .File is parsed.
			},
			intervals: []string{
				"1000,1001",
				"2000,2001",
				"103000,103001", // Use a voffset larger than 16 bits to check that .File is parsed.
				"4000,4001",
			},
			unmapped: 999,
		},
		{
			bins: []string{
				"100,1,2:200,3,4:37450,5,6,7,8",
			},
			intervals: []string{
				"1000,1001",
			},
			unmapped: -1,
		},
	}

	for _, test := range tests {
		expectedOffsets := map[int][]bgzf.Offset{}

		// Write an index to buf.
		buf := writeIndex(t, test.bins, test.intervals, test.unmapped)
		index, err := ReadIndex(buf)
		require.Nil(t, err)

		assert.Equal(t, [4]byte{'B', 'A', 'I', 0x1}, index.Magic)
		assert.Equal(t, len(test.bins), len(index.Refs))
		for refId := range test.bins {
			expectedOffsets[refId] = make([]bgzf.Offset, 0)

			// Check the bins
			bins := strings.Split(test.bins[refId], ":")
			binCount := 0
			for binNum, binString := range bins {
				binInfo := strings.Split(binString, ",")
				if toInt(t, binInfo[0]) == 37450 {
					// Check meta chunk info
					assert.Equal(t, uint64(toInt(t, binInfo[1])), index.Refs[refId].Meta.UnmappedBegin)
					assert.Equal(t, uint64(toInt(t, binInfo[2])), index.Refs[refId].Meta.UnmappedEnd)
					assert.Equal(t, uint64(toInt(t, binInfo[3])), index.Refs[refId].Meta.MappedCount)
					assert.Equal(t, uint64(toInt(t, binInfo[4])), index.Refs[refId].Meta.UnmappedCount)
				} else {
					// Check regular chunk offsets
					assert.Equal(t, uint32(toInt(t, binInfo[0])), index.Refs[refId].Bins[binNum].BinNum)
					for i := 1; i < len(binInfo)-1; i += 2 {
						beginFile := int64(toInt(t, binInfo[i])) >> 16
						beginBlock := uint16(toInt(t, binInfo[i]))
						endFile := int64(toInt(t, binInfo[i+1])) >> 16
						endBlock := uint16(toInt(t, binInfo[i+1]))

						assert.Equal(t, beginFile, index.Refs[refId].Bins[binNum].Chunks[(i-1)/2].Begin.File)
						assert.Equal(t, beginBlock, index.Refs[refId].Bins[binNum].Chunks[(i-1)/2].Begin.Block)
						assert.Equal(t, endFile, index.Refs[refId].Bins[binNum].Chunks[(i-1)/2].End.File)
						assert.Equal(t, endBlock, index.Refs[refId].Bins[binNum].Chunks[(i-1)/2].End.Block)

						expectedOffsets[refId] = append(expectedOffsets[refId], bgzf.Offset{beginFile, beginBlock})
					}
					binCount++
				}
			}
			assert.Equal(t, binCount, len(index.Refs[refId].Bins))

			// Check the intervals
			intervals := strings.Split(test.intervals[refId], ",")
			assert.Equal(t, len(intervals), len(index.Refs[refId].Intervals))
			for i, intervalStr := range intervals {
				intervalFile := int64(toInt(t, intervalStr)) >> 16
				intervalBlock := uint16(toInt(t, intervalStr))

				assert.Equal(t, intervalFile, index.Refs[refId].Intervals[i].File)
				assert.Equal(t, intervalBlock, index.Refs[refId].Intervals[i].Block)

				expectedOffsets[refId] = append(expectedOffsets[refId], bgzf.Offset{intervalFile, intervalBlock})
			}

			// Sort and unique the expected offsets
			sort.SliceStable(expectedOffsets[refId], func(i, j int) bool {
				c0 := expectedOffsets[refId][i]
				c1 := expectedOffsets[refId][j]
				if c0.File != c1.File {
					return c0.File < c1.File
				}
				return c0.Block < c1.Block
			})
			uniq := make([]bgzf.Offset, 0)
			previous := bgzf.Offset{-1, 0}
			for _, offset := range expectedOffsets[refId] {
				if offset != previous {
					uniq = append(uniq, offset)
					previous = offset
				}
			}
			expectedOffsets[refId] = uniq

		}
		if test.unmapped >= 0 {
			assert.Equal(t, uint64(test.unmapped), *index.UnmappedCount)
		} else {
			assert.Nil(t, index.UnmappedCount)
		}

		// Verify AllOffsets output.
		actualOffsets := index.AllOffsets()
		assert.Equal(t, len(expectedOffsets), len(actualOffsets))
		for refId, expected := range expectedOffsets {
			actual, found := actualOffsets[refId]
			require.True(t, found)
			assert.True(t, reflect.DeepEqual(expected, actual), "%v %v", expected, actual)
		}
	}
}
