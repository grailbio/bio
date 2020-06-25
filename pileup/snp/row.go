// Copyright 2020 Grail Inc.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
package snp

import (
	"encoding/binary"

	"github.com/grailbio/bio/pileup"
)

type PerReadFeatures struct {
	// dist5p is the 0-based distance of the current base from its 5' end.  (Note
	// that Dist3p := fraglen - 1 - dist5p, so we don't need to store it
	// separately.)
	Dist5p uint16

	Fraglen uint16
	Qual    byte
	Strand  byte
}

const (
	FieldCounts = 1 << iota
	FieldPerReadA
	FieldPerReadC
	FieldPerReadG
	FieldPerReadT
	FieldPerReadAny = FieldPerReadA | FieldPerReadC | FieldPerReadG | FieldPerReadT
)

// PileupPayload is a container for all types of pileup data which may be
// associated with a single position.  It does not store the position itself,
// or a tag indicating which parts of the container are used.
//
// Depth and count values are of type uint32 instead of int to reduce cache
// footprint.
type PileupPayload struct {
	Depth   uint32
	Counts  [pileup.NBaseEnum][2]uint32
	PerRead [pileup.NBase][]PerReadFeatures
}

// PileupRow contains all pileup data associated with a single position, along
// with the position itself and the set of PileupPayload fields used.
//
// The main loop splits the genome into shards, and generates lightly
// compressed (zstd level 1) per-shard PileupRow recordio files.  Then, the
// per-shard files are read in sequence and converted to the final requested
// output format.  This is a bit inefficient, but we can easily afford it.
type PileupRow struct {
	FieldsPresent uint32 // field... flags
	RefID         uint32
	Pos           uint32
	Payload       PileupPayload
}

// cutAndAdvance() returns s[offset:offset+pieceLen], and increments offset by
// pieceLen.
//
// This is the lowest-overhead idiom I've found so far for filling a
// preallocated []byte.  In particular, despite recent advances in
// bounds-check-elimination, the Go 1.12 compiler still does not recognize that
// s[offset:offset+k] has length k when k is a hardcoded integer and offset is
// known to be smaller than (MaxInt - k); it really is necessary to write x :=
// s[offset:] followed by x = x[:k] if you want to avoid spurious
// bounds-checking.
//
// Other things I tried:
// - "splitBefore(s *[]byte, x int) []byte { ... }".  This is a slightly
//   simpler interface, but unfortunately it proved to be higher-overhead due
//   to length and capacity requiring separate updates/checks.
// - Making offset and pieceLen into unsigned integers.  That makes no
//   difference to the compiled code; the compiler is able to prove that offset
//   is never negative.
func cutAndAdvance(offset *int, s []byte, pieceLen int) []byte {
	tmpSlice := s[(*offset):]
	*offset += pieceLen
	return tmpSlice[:pieceLen]
}

// Serialized format:
//   [0..4): fieldsPresent
//   [4..8): refID
//   [8..12): pos
//   [12..16): depth
//   if counts present, stored in next 40 bytes
//   if perRead[pileup.baseA] present, length stored in next 4 bytes, then
//     values stored in next 6*n bytes
//   if perRead[pileup.baseC] present... etc.
// This is essentially the simplest format that can support the variable-length
// per-read feature arrays that are needed.  It is not difficult to decrease
// the nominal size of these records by (i) using varints instead of uint32s,
// and (ii) making fieldsPresent indicate which counts[][] values are nonzero
// and only storing those; but I wouldn't expect that to be worth the
// additional complexity since all uses of this marshal function are bundled
// with the "zstd 1" transformer anyway.  (Instead, all the 'extra' complexity
// in this function concerns (i) avoiding extra allocations and (ii) avoiding a
// ridiculous number of spurious bounds-checks, in ways that make sense for a
// wide variety of other serialization functions.)
//
// In the future, we may need to add indel support.
func MarshalPileupRow(scratch []byte, p interface{}) ([]byte, error) {
	pr := p.(*PileupRow)
	fieldsPresent := pr.FieldsPresent
	// Compute length up-front so that, if we need to allocate, we only do so
	// once.
	bytesReq := 16
	if fieldsPresent&FieldCounts != 0 {
		bytesReq += 40
	}
	if fieldsPresent&FieldPerReadAny != 0 {
		for b := range pr.Payload.PerRead {
			// For b in {0,1,2,3}, (FieldPerReadA << b) is the bit indicating that
			// perRead[b] must be stored.
			if fieldsPresent&(FieldPerReadA<<uint(b)) != 0 {
				bytesReq += 4 + 6*len(pr.Payload.PerRead[b])
			}
		}
	}
	t := scratch
	if len(t) < bytesReq {
		t = make([]byte, bytesReq)
	}

	offset := 0
	tStart := cutAndAdvance(&offset, t, 16)
	binary.LittleEndian.PutUint32(tStart[0:4], pr.FieldsPresent)
	binary.LittleEndian.PutUint32(tStart[4:8], pr.RefID)
	binary.LittleEndian.PutUint32(tStart[8:12], pr.Pos)
	binary.LittleEndian.PutUint32(tStart[12:16], pr.Payload.Depth)
	if fieldsPresent&FieldCounts != 0 {
		tCounts := cutAndAdvance(&offset, t, 40)
		// Unfortunately, while the obvious double-loop works fine for reading
		// values from pr.Payload.Counts[], I don't see any way to express the
		// writes to tCounts[] that the Go 1.12 bounds-check-eliminator
		// understands.
		binary.LittleEndian.PutUint32(tCounts[:4], pr.Payload.Counts[pileup.BaseA][0])
		binary.LittleEndian.PutUint32(tCounts[4:8], pr.Payload.Counts[pileup.BaseA][1])
		binary.LittleEndian.PutUint32(tCounts[8:12], pr.Payload.Counts[pileup.BaseC][0])
		binary.LittleEndian.PutUint32(tCounts[12:16], pr.Payload.Counts[pileup.BaseC][1])
		binary.LittleEndian.PutUint32(tCounts[16:20], pr.Payload.Counts[pileup.BaseG][0])
		binary.LittleEndian.PutUint32(tCounts[20:24], pr.Payload.Counts[pileup.BaseG][1])
		binary.LittleEndian.PutUint32(tCounts[24:28], pr.Payload.Counts[pileup.BaseT][0])
		binary.LittleEndian.PutUint32(tCounts[28:32], pr.Payload.Counts[pileup.BaseT][1])
		binary.LittleEndian.PutUint32(tCounts[32:36], pr.Payload.Counts[pileup.BaseX][0])
		binary.LittleEndian.PutUint32(tCounts[36:40], pr.Payload.Counts[pileup.BaseX][1])
	}
	if fieldsPresent&FieldPerReadAny != 0 {
		for b := range pr.Payload.PerRead {
			if fieldsPresent&(FieldPerReadA<<uint(b)) != 0 {
				lenSlice := cutAndAdvance(&offset, t, 4)
				binary.LittleEndian.PutUint32(lenSlice, uint32(len(pr.Payload.PerRead[b])))
				for _, src := range pr.Payload.PerRead[b] {
					dst := cutAndAdvance(&offset, t, 6)
					binary.LittleEndian.PutUint16(dst[:2], src.Dist5p)
					binary.LittleEndian.PutUint16(dst[2:4], src.Fraglen)
					dst[4] = src.Qual
					dst[5] = src.Strand
				}
			}
		}
	}
	return t, nil
}

// tried the block-unmarshal strategy in grail.com/bio/variants, it actually
// seemed to have worse performance for this use case
func unmarshalPileupRow(in []byte) (out interface{}, err error) {
	offset := 0
	inStart := cutAndAdvance(&offset, in, 16)
	pr := &PileupRow{
		FieldsPresent: binary.LittleEndian.Uint32(inStart[:4]),
		RefID:         binary.LittleEndian.Uint32(inStart[4:8]),
		Pos:           binary.LittleEndian.Uint32(inStart[8:12]),
	}
	pr.Payload.Depth = binary.LittleEndian.Uint32(inStart[12:16])
	if pr.FieldsPresent&FieldCounts != 0 {
		inCounts := cutAndAdvance(&offset, in, 40)
		pr.Payload.Counts[pileup.BaseA][0] = binary.LittleEndian.Uint32(inCounts[0:4])
		pr.Payload.Counts[pileup.BaseA][1] = binary.LittleEndian.Uint32(inCounts[4:8])
		pr.Payload.Counts[pileup.BaseC][0] = binary.LittleEndian.Uint32(inCounts[8:12])
		pr.Payload.Counts[pileup.BaseC][1] = binary.LittleEndian.Uint32(inCounts[12:16])
		pr.Payload.Counts[pileup.BaseG][0] = binary.LittleEndian.Uint32(inCounts[16:20])
		pr.Payload.Counts[pileup.BaseG][1] = binary.LittleEndian.Uint32(inCounts[20:24])
		pr.Payload.Counts[pileup.BaseT][0] = binary.LittleEndian.Uint32(inCounts[24:28])
		pr.Payload.Counts[pileup.BaseT][1] = binary.LittleEndian.Uint32(inCounts[28:32])
		pr.Payload.Counts[pileup.BaseX][0] = binary.LittleEndian.Uint32(inCounts[32:36])
		pr.Payload.Counts[pileup.BaseX][1] = binary.LittleEndian.Uint32(inCounts[36:40])
	}
	if pr.FieldsPresent&FieldPerReadAny != 0 {
		for b := range pr.Payload.PerRead {
			if pr.FieldsPresent&(FieldPerReadA<<uint(b)) != 0 {
				lenSlice := cutAndAdvance(&offset, in, 4)
				curLen := binary.LittleEndian.Uint32(lenSlice)

				// If we wanted to further reduce the number of small allocations, we
				// could allocate a single []PerReadFeatures slice outside this loop,
				// and then make the per-base slices point to subslices of the single
				// allocation.  (I don't bother since, most of the time, only one
				// newFeatures slice corresponding to the REF base is allocated
				// anyway.)
				newFeatures := make([]PerReadFeatures, curLen)

				pr.Payload.PerRead[b] = newFeatures
				for i := range newFeatures {
					src := cutAndAdvance(&offset, in, 6)
					newFeatures[i].Dist5p = binary.LittleEndian.Uint16(src[:2])
					newFeatures[i].Fraglen = binary.LittleEndian.Uint16(src[2:4])
					newFeatures[i].Qual = src[4]
					newFeatures[i].Strand = src[5]
				}
			}
		}
	}
	return pr, nil
}
