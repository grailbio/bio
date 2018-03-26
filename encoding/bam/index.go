package bam

import (
	"encoding/binary"
	"fmt"
	"io"
	"sort"

	"github.com/biogo/hts/bgzf"
)

// Index represents the content of a .bai index file (for use with a .bam file).
type Index struct {
	Magic         [4]byte
	Refs          []Reference
	UnmappedCount *uint64
}

// Reference represents the reference data within a .bai file.
type Reference struct {
	Bins      []Bin
	Intervals []bgzf.Offset
	Meta      Metadata
}

// Bin represents the bin data within a .bai file.
type Bin struct {
	BinNum uint32
	Chunks []Chunk
}

// Chunk represents the Chunk data within a .bai file.
type Chunk struct {
	Begin bgzf.Offset
	End   bgzf.Offset
}

// Metadata represents the Metadata data within a .bai file.
type Metadata struct {
	UnmappedBegin uint64
	UnmappedEnd   uint64
	MappedCount   uint64
	UnmappedCount uint64
}

// ReadIndex parses the content of r and returns an Index or nil and an error.
func ReadIndex(r io.Reader) (*Index, error) {
	i := &Index{}

	if _, err := io.ReadFull(r, i.Magic[0:]); err != nil {
		return nil, err
	}
	if i.Magic != [4]byte{'B', 'A', 'I', 0x1} {
		return nil, fmt.Errorf("bam index invalid magic: %v", i.Magic)
	}

	var refCount int32
	if err := binary.Read(r, binary.LittleEndian, &refCount); err != nil {
		return nil, err
	}
	i.Refs = make([]Reference, refCount)

	// Read each Reference
	for refId := 0; int32(refId) < refCount; refId++ {
		// Read each Bin
		var binCount int32
		if err := binary.Read(r, binary.LittleEndian, &binCount); err != nil {
			return nil, err
		}
		ref := Reference{
			Bins: make([]Bin, 0, binCount),
		}
		for b := 0; int32(b) < binCount; b++ {
			var binNum uint32
			if err := binary.Read(r, binary.LittleEndian, &binNum); err != nil {
				return nil, err
			}
			var chunkCount int32
			if err := binary.Read(r, binary.LittleEndian, &chunkCount); err != nil {
				return nil, err
			}

			bin := Bin{
				BinNum: binNum,
				Chunks: make([]Chunk, chunkCount),
			}

			// Read each Chunk
			for c := 0; int32(c) < chunkCount; c++ {
				var beginOffset uint64
				if err := binary.Read(r, binary.LittleEndian, &beginOffset); err != nil {
					return nil, err
				}
				var endOffset uint64
				if err := binary.Read(r, binary.LittleEndian, &endOffset); err != nil {
					return nil, err
				}
				bin.Chunks[c] = Chunk{
					Begin: toOffset(beginOffset),
					End:   toOffset(endOffset),
				}
			}

			if binNum == 37450 {
				// If we have a metadata chunk, put it in ref.Meta instead of ref.Bins.
				if len(bin.Chunks) != 2 {
					return nil, fmt.Errorf("Invalid metadata chunk has %d chunks, should have 2", len(bin.Chunks))
				}
				ref.Meta = Metadata{
					UnmappedBegin: fromOffset(bin.Chunks[0].Begin),
					UnmappedEnd:   fromOffset(bin.Chunks[0].End),
					MappedCount:   fromOffset(bin.Chunks[1].Begin),
					UnmappedCount: fromOffset(bin.Chunks[1].End),
				}
			} else {
				ref.Bins = append(ref.Bins, bin)
			}
		}

		// Read each Interval.
		var intervalCount int32
		if err := binary.Read(r, binary.LittleEndian, &intervalCount); err != nil {
			return nil, err
		}
		ref.Intervals = make([]bgzf.Offset, intervalCount)
		for inv := 0; int32(inv) < intervalCount; inv++ {
			var ioffset uint64
			if err := binary.Read(r, binary.LittleEndian, &ioffset); err != nil {
				return nil, err
			}
			ref.Intervals[inv] = toOffset(ioffset)
		}
		i.Refs[refId] = ref
	}

	var unmappedCount uint64
	if err := binary.Read(r, binary.LittleEndian, &unmappedCount); err == nil {
		i.UnmappedCount = &unmappedCount
	} else if err != nil && err != io.EOF {
		return nil, err
	}
	return i, nil
}

// AllOffsets returns a map of chunk offsets in the index file, it
// includes chunk begin locations, and interval locations.  The Key of
// the map is the Reference ID, and the value is a slice of
// bgzf.Offsets.  The return map will have an entry for every
// reference ID, even if the list of offsets is empty.
func (i *Index) AllOffsets() map[int][]bgzf.Offset {
	m := make(map[int][]bgzf.Offset)
	for refId, ref := range i.Refs {
		m[refId] = make([]bgzf.Offset, 0)

		// Get the offsets for this ref.
		for _, bin := range ref.Bins {
			for _, chunk := range bin.Chunks {
				if chunk.Begin.File != 0 || chunk.Begin.Block != 0 {
					m[refId] = append(m[refId], chunk.Begin)
				}
			}
		}
		for _, interval := range ref.Intervals {
			if interval.File != 0 || interval.Block != 0 {
				m[refId] = append(m[refId], interval)
			}
		}

		// Sort the offsets
		sort.SliceStable(m[refId], func(i, j int) bool {
			c0 := m[refId][i]
			c1 := m[refId][j]
			if c0.File != c1.File {
				return c0.File < c1.File
			}
			return c0.Block < c1.Block
		})

		// Keep only unique offsets
		uniq := make([]bgzf.Offset, 0)
		previous := bgzf.Offset{-1, 0}
		for _, offset := range m[refId] {
			if offset != previous {
				uniq = append(uniq, offset)
				previous = offset
			}
		}
		m[refId] = uniq
	}
	return m
}

func toOffset(voffset uint64) bgzf.Offset {
	return bgzf.Offset{
		File:  int64(voffset >> 16),
		Block: uint16(voffset),
	}
}

func fromOffset(offset bgzf.Offset) uint64 {
	return uint64(offset.File<<16) | uint64(offset.Block)
}
