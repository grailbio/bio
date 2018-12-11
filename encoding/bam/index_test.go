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

	"github.com/grailbio/hts/bgzf"
	"github.com/grailbio/testutil/assert"
	"github.com/grailbio/testutil/expect"
)

func toInt(t *testing.T, s string) int {
	i, err := strconv.Atoi(s)
	assert.Nil(t, err)
	return i
}

func writeBin(t *testing.T, w io.Writer, s string) {
	bins := strings.Split(s, ":")
	// Write the number of bins
	err := binary.Write(w, binary.LittleEndian, int32(len(bins)))
	expect.Nil(t, err)

	for _, bin := range bins {
		binContent := strings.Split(bin, ",")

		// Write the bin number
		err = binary.Write(w, binary.LittleEndian, uint32(toInt(t, binContent[0])))
		expect.Nil(t, err)
		binContent = binContent[1:]

		// Write the number of chunks
		err = binary.Write(w, binary.LittleEndian, int32(len(binContent)/2))
		expect.Nil(t, err)

		// Write the chunks
		for _, voffset := range binContent {
			err = binary.Write(w, binary.LittleEndian, uint64(toInt(t, voffset)))
			expect.Nil(t, err)
		}
	}
}

func writeIntervals(t *testing.T, w io.Writer, s string) {
	intervals := strings.Split(s, ",")

	// Write the number of intervals
	err := binary.Write(w, binary.LittleEndian, int32(len(intervals)))
	expect.Nil(t, err)

	for _, voffset := range intervals {
		err = binary.Write(w, binary.LittleEndian, uint64(toInt(t, voffset)))
		expect.Nil(t, err)
	}
}

func writeIndex(t *testing.T, bins, intervals []string, unmapped int) *bytes.Buffer {
	var buf bytes.Buffer
	magic := []byte{'B', 'A', 'I', 0x1}
	_, err := buf.Write(magic)
	expect.Nil(t, err)

	// Two references
	err = binary.Write(&buf, binary.LittleEndian, int32(len(bins)))
	expect.Nil(t, err)

	for i := range bins {
		writeBin(t, &buf, bins[i])
		writeIntervals(t, &buf, intervals[i])
	}

	// Write unmapped count
	if unmapped >= 0 {
		err = binary.Write(&buf, binary.LittleEndian, uint64(unmapped))
		expect.Nil(t, err)
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
		assert.Nil(t, err)

		expect.EQ(t, index.Magic, [4]byte{'B', 'A', 'I', 0x1})
		expect.EQ(t, len(index.Refs), len(test.bins))
		for refId := range test.bins {
			expectedOffsets[refId] = make([]bgzf.Offset, 0)

			// Check the bins
			bins := strings.Split(test.bins[refId], ":")
			binCount := 0
			for binNum, binString := range bins {
				binInfo := strings.Split(binString, ",")
				if toInt(t, binInfo[0]) == 37450 {
					// Check meta chunk info
					expect.EQ(t, index.Refs[refId].Meta.UnmappedBegin, uint64(toInt(t, binInfo[1])))
					expect.EQ(t, index.Refs[refId].Meta.UnmappedEnd, uint64(toInt(t, binInfo[2])))
					expect.EQ(t, index.Refs[refId].Meta.MappedCount, uint64(toInt(t, binInfo[3])))
					expect.EQ(t, index.Refs[refId].Meta.UnmappedCount, uint64(toInt(t, binInfo[4])))
				} else {
					// Check regular chunk offsets
					expect.EQ(t, index.Refs[refId].Bins[binNum].BinNum, uint32(toInt(t, binInfo[0])))
					for i := 1; i < len(binInfo)-1; i += 2 {
						beginFile := int64(toInt(t, binInfo[i])) >> 16
						beginBlock := uint16(toInt(t, binInfo[i]))
						endFile := int64(toInt(t, binInfo[i+1])) >> 16
						endBlock := uint16(toInt(t, binInfo[i+1]))

						expect.EQ(t, index.Refs[refId].Bins[binNum].Chunks[(i-1)/2].Begin.File, beginFile)
						expect.EQ(t, index.Refs[refId].Bins[binNum].Chunks[(i-1)/2].Begin.Block, beginBlock)
						expect.EQ(t, index.Refs[refId].Bins[binNum].Chunks[(i-1)/2].End.File, endFile)
						expect.EQ(t, index.Refs[refId].Bins[binNum].Chunks[(i-1)/2].End.Block, endBlock)

						expectedOffsets[refId] = append(expectedOffsets[refId], bgzf.Offset{beginFile, beginBlock})
					}
					binCount++
				}
			}
			expect.EQ(t, len(index.Refs[refId].Bins), binCount)

			// Check the intervals
			intervals := strings.Split(test.intervals[refId], ",")
			expect.EQ(t, len(index.Refs[refId].Intervals), len(intervals))
			for i, intervalStr := range intervals {
				intervalFile := int64(toInt(t, intervalStr)) >> 16
				intervalBlock := uint16(toInt(t, intervalStr))

				expect.EQ(t, index.Refs[refId].Intervals[i].File, intervalFile)
				expect.EQ(t, index.Refs[refId].Intervals[i].Block, intervalBlock)

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
			expect.EQ(t, *index.UnmappedCount, uint64(test.unmapped))
		} else {
			expect.Nil(t, index.UnmappedCount)
		}

		// Verify AllOffsets output.
		actualOffsets := index.AllOffsets()
		expect.EQ(t, len(actualOffsets), len(expectedOffsets))
		for refId, expected := range expectedOffsets {
			actual, found := actualOffsets[refId]
			assert.True(t, found)
			expect.True(t, reflect.DeepEqual(expected, actual), "%v %v", expected, actual)
		}
	}
}
