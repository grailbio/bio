package markduplicates

import (
	"github.com/grailbio/base/log"
	"github.com/grailbio/base/simd"
	"github.com/grailbio/bio/encoding/bam"
	"github.com/grailbio/hts/sam"
)

type strand int8

var (
	rgTag = sam.Tag{'R', 'G'}
	diTag = sam.Tag{'D', 'I'}
	dlTag = sam.Tag{'D', 'L'}
	dsTag = sam.Tag{'D', 'S'}
	dtTag = sam.Tag{'D', 'T'}
	duTag = sam.Tag{'D', 'U'}
)

func mateInPaddedShard(shard *bam.Shard, r *sam.Record) bool {
	return shard.CoordInShard(shard.Padding, bam.NewCoord(r.MateRef, r.MatePos, 0))
}

func abs(x int) int {
	if x < 0 {
		return -x
	}
	return x
}

func min(x, y int) int {
	if x < y {
		return x
	}
	return y
}

func baseQScore(r *sam.Record) int {
	s := simd.Accumulate8Greater(r.Qual, 14)
	s = min(s, 32767/2) // use the same clamping as picard
	if bam.IsQCFailed(r) {
		s -= (32768 / 2)
	}
	return s
}

func getReadGroup(r *sam.Record) (string, bool) {
	aux := r.AuxFields.Get(rgTag)
	if aux == nil {
		return "", false
	}
	return aux.Value().(string), true
}

// GetLibrary returns the library for the given record's read group.
// If the library is not defined in readGroupLibrary, returns "Unknown
// Library".
func GetLibrary(readGroupLibrary map[string]string, record *sam.Record) string {
	const unknown = "Unknown Library"

	readGroup, found := getReadGroup(record)
	if !found {
		return unknown
	}

	library := readGroupLibrary[readGroup]
	if library == "" {
		return unknown
	}
	return library
}

func clearDupFlagTags(r *sam.Record) {
	r.Flags &^= sam.Duplicate

	tagsToRemove := []sam.Tag{diTag, dlTag, dsTag, dtTag, duTag}
	bam.ClearAuxTags(r, tagsToRemove)
}

// GetR1R2Orientation returns an orientation byte containing
// orientations for both R1 and R2.
func GetR1R2Orientation(p *IndexedPair) Orientation {
	if p.Left.R.Flags&sam.Read1 == p.Right.R.Flags&sam.Read1 {
		log.Fatalf("Both reads are first or second for pair: %v %d %d", p.Left.R.Name, p.Left.R.Flags, p.Right.R.Flags)
	}

	if p.Left.R.Flags&sam.Read1 != 0 {
		return orientationBytePair(p.Left.R.Flags&sam.Reverse != 0, p.Right.R.Flags&sam.Reverse != 0)
	} else if p.Right.R.Flags&sam.Read1 != 0 {
		return orientationBytePair(p.Right.R.Flags&sam.Reverse != 0, p.Left.R.Flags&sam.Reverse != 0)
	} else {
		log.Fatalf("Could not find first read in pair: %v", p.Left.R.Name)
	}
	return 0
}

// r1Strand returns +1 or -1 depending on the strand if the reads
// point in opposite directions. If the two reads point in the same
// direction, return 0. For singletons, return the strand for just the
// singleton, ignoring the mate's direction.
func r1Strand(r *sam.Record) strand {
	if r.Flags&sam.MateUnmapped != sam.MateUnmapped && r.Flags&sam.Reverse == r.Flags&sam.MateReverse {
		return 0
	}
	if bam.IsRead1(r) {
		return strand(r.Strand())
	}
	return strand(-r.Strand())
}
