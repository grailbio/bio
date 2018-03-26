package biopb

// This file adds convenience comparison methods to Coord and CoordRange.
//
// TODO(saito) We may want to move them to bion/encoding/bam and make them
// ordinary functions as opposed to Coord / CoordRange methods.

import (
	"math"
)

const (
	// InfinityPos is 1+ the largest possible alignment position.
	InfinityPos = math.MaxInt32

	// LimitValidRefID is a pseudo referenceID for a max possible valid
	// reference.  For example, passing RecRange{{0,0},{LimitValidRefID,
	// InfinityPos}} to ReadOpts.Range will read all mapped sequences.
	LimitValidRefID = math.MaxInt32 - 1

	// InfinityRefID is a pseudo referenceID for unmapped reads. For
	// example, passing RecRange{{UnmappedRefID,0},{UnmappedRefID,
	// InfinityPos}} to ReadOpts.Range will read all mapped sequences.
	InfinityRefID = int32(-1)
	// UnmappedRefID is a synonym of InfinityRefID.
	UnmappedRefID = InfinityRefID

	// InvalidRefID is used as a sentinel. We use -2 because -1 is is taken
	// by UnmappedRefID.
	InvalidRefID = int32(-2)
	// InvalidPos is a sentinel position value. We use -2 instead of -1
	// because -1 is sometimes used by the position of unmapped reads.
	InvalidPos = int32(-2)
)

// For sorting biopb.Coords.
func sortableRefID(id int32) int32 {
	if id == InfinityRefID {
		// Unmapped reads are sorted the last, so use a large value.
		return math.MaxInt32
	}
	return id
}

// Compare returns (negative int, 0, positive int) if (r<r1, r=r1, r>r1)
// respectively.
func (r Coord) Compare(r1 Coord) int {
	refid0 := sortableRefID(r.RefId)
	refid1 := sortableRefID(r1.RefId)
	if refid0 != refid1 {
		return int(refid0 - refid1)
	}
	if r.Pos != r1.Pos {
		return int(r.Pos - r1.Pos)
	}
	return int(r.Seq - r1.Seq)
}

// LT returns true iff r < r1.
func (r Coord) LT(r1 Coord) bool {
	return r.Compare(r1) < 0
}

// LE returns true iff r <= r1
func (r Coord) LE(r1 Coord) bool {
	return r.Compare(r1) <= 0
}

// GE returns true iff r >= r1
func (r Coord) GE(r1 Coord) bool {
	return r.Compare(r1) >= 0
}

// GT return true iff r > r1
func (r Coord) GT(r1 Coord) bool {
	return r.Compare(r1) > 0
}

// EQ returns true iff r = r1.
func (r Coord) EQ(r1 Coord) bool {
	return r.RefId == r1.RefId && r.Pos == r1.Pos && r.Seq == r1.Seq
}

// Min returns the smaller of r and r1.
func (r Coord) Min(r1 Coord) Coord {
	if r.LT(r1) {
		return r
	}
	return r1
}

// EQ returns true iff. r=r1.
func (r CoordRange) EQ(r1 CoordRange) bool {
	return r.Start.EQ(r1.Start) && r.Limit.EQ(r1.Limit)
}

// Intersects returns true iff (a ∩ r) != ∅
func (r CoordRange) Intersects(r1 CoordRange) bool {
	return r.Start.LT(r1.Limit) && r1.Start.LT(r.Limit)
}

// Contains checks if "a" is inside the "r"
func (r CoordRange) Contains(a Coord) bool {
	return r.Start.LE(a) && a.LT(r.Limit)
}

// ContainsRange returns true iff (a ∩ r) = a.
func (r CoordRange) ContainsRange(a CoordRange) bool {
	return r.Start.LE(a.Start) && r.Limit.GE(a.Limit)
}
