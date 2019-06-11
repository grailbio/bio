package interval

import (
	"math"
	"sort"
)

// This file includes support datatypes and functions for representing an
// interval-union as an []int32 containing a sorted sequence of
// interval-endpoints, and iterating over the intervals.
//
// For example, given the intervals
//   [5, 15)
//   [7, 17)
//   [20, 25)
// the interval-union would be
//   [5, 17) U [20, 25)
// so the sorted sequence of endpoints would be
//   {5, 17, 20, 25}.
//
// UnionScanner can be used to iterate over these positions as follows:
//   endpoints := []PosType{5, 17, 20, 25}
//   us := NewUnionScanner(endpoints)
//   var start PosType
//   var end PosType
//   for us.Scan(&start, &end, 22) {
//     for pos := start; pos < end; pos++ {
//       fmt.Printf("%d ", pos)
//     }
//   }
//   fmt.Printf("\n")
// This prints "5 6 7 8 9 10 11 12 13 14 15 16 20 21 ".
// We can follow up with
//   for us.Scan(&start, &end, 30) {
//     for pos := start; pos < end; pos++ {
//       fmt.Printf("%d ", pos)
//     }
//   }
//   fmt.Printf("\n")
// to pick up where we left off and print "22 23 24 ".
//
// If your use case is a bit different, EndpointIndex provides a few
// lower-level functions which may come in handy.

// PosType is the type used to represent interval coordinates.  int32 should be
// wide enough for some time to come, since that's what BAM is limited to.
//
// (This, and PosTypeMax, should move to a more central package once an
// appropriate one exists.  And then, when generics finally become part of the
// language *crosses fingers*, we can allow some applications to redefine this
// as uint32 or a 64-bit type as appropriate.)
type PosType int32

// PosTypeMax is the maximum value that can be represented by a PosType.
const PosTypeMax = math.MaxInt32

// SearchPosTypes returns the index of x in a[], or the position where x would
// be inserted if x isn't in a (this could be len(a)).  It's exactly the same
// as sort.SearchInts(), except for PosType.
func SearchPosTypes(a []PosType, x PosType) EndpointIndex {
	return EndpointIndex(sort.Search(len(a), func(i int) bool { return a[i] >= x }))
}

// ExpsearchPosType performs "exponential search"
// (https://en.wikipedia.org/wiki/Exponential_search ), checking a[idx], then
// a[idx + 1], then a[idx + 3], then a[idx + 7], etc., and finishing with
// binary search once it's either found an element larger than the target or
// has hit the end of the slice.  It's usually a better choice than
// SearchPosTypes when iterating.
// (However, an inlined simple linear search may be better in practice.  Can
// benchmark later if it matters.)
func ExpsearchPosType(a []PosType, x PosType, idx EndpointIndex) EndpointIndex {
	nextIncr := EndpointIndex(1)
	startIdx := idx
	endIdx := EndpointIndex(len(a))
	for idx < endIdx {
		if a[idx] >= x {
			endIdx = idx
			break
		}
		startIdx = idx + 1
		idx += nextIncr
		nextIncr *= 2
	}
	// This is really just an inlined sort.Search call.  We spell it out since
	// startIdx is usually equal to endIdx, and the compiler doesn't inline
	// anything with a loop for now.
	for startIdx < endIdx {
		midIdx := EndpointIndex((uint(startIdx) + uint(endIdx)) >> 1)
		if a[midIdx] >= x {
			endIdx = midIdx
		} else {
			startIdx = midIdx + 1
		}
	}
	return startIdx
}

// EndpointIndex is intended to represent the result of
// SearchPosTypes(endpoints, pos+1).
// NOTE THE "+1"!  This is necessary to get SearchPosTypes to line up with our
// usual left-closed right-open intervals.
type EndpointIndex uint32

// NewEndpointIndex returns an EndpointIndex initialized to
// SearchPosTypes(endpoints, pos+1).
func NewEndpointIndex(pos PosType, endpoints []PosType) EndpointIndex {
	return SearchPosTypes(endpoints, pos+1)
}

// Contained returns whether we're inside an interval.
func (ei EndpointIndex) Contained() bool {
	return ei&1 != 0
}

// Finished returns whether we're past all the intervals.
func (ei EndpointIndex) Finished(endpoints []PosType) bool {
	return ei >= EndpointIndex(len(endpoints))
}

// Begin returns:
// - the index for the beginning of the current interval, if we're inside an
//   interval
// - otherwise, the index for the beginning of the next interval
func (ei EndpointIndex) Begin() EndpointIndex {
	return ei & (^EndpointIndex(1))
}

// Update updates the EndpointIndex to refer to newPos, which cannot be smaller
// than the previous position referred to by this EndpointIndex.  It is
// substantially faster than NewEndpointIndex when the position is increasing
// slowly.
func (ei *EndpointIndex) Update(newPos PosType, endpoints []PosType) {
	*ei = ExpsearchPosType(endpoints, newPos+1, *ei)
}

// UnionScanner supports iteration over an interval-union.
// Invariants:
//   endpointIdx == SearchPosTypes(endpoints, pos+1)
//   Pos is either contained in an interval, or is PosTypeMax
type UnionScanner struct {
	endpoints   []PosType
	pos         PosType
	endpointIdx EndpointIndex
}

// NewUnionScanner returns a UnionScanner initialized to the first interval.
func NewUnionScanner(endpoints []PosType) UnionScanner {
	startPos := PosType(PosTypeMax)
	startEndpointIdx := EndpointIndex(0)
	// May as well make this not crash when there are no intervals.
	if len(endpoints) >= 1 {
		startPos = endpoints[0]
		startEndpointIdx = 1
	}
	return UnionScanner{
		endpoints:   endpoints,
		pos:         startPos,
		endpointIdx: startEndpointIdx,
	}
}

// Pos returns the next position to be iterated over, or PosTypeMax if there
// aren't any.
func (us *UnionScanner) Pos() PosType {
	return us.pos
}

// Scan is written so that the following loop can be used to iterate over all
// within-interval positions up to (and not including) limit:
//   for us.Scan(&start, &end, limit) {
//     for pos := start; pos < end; pos++ {
//       // ...do stuff with pos...
//     }
//   }
func (us *UnionScanner) Scan(start *PosType, end *PosType, limit PosType) bool {
	if us.pos >= limit {
		return false
	}
	*start = us.pos
	intervalEnd := us.endpoints[us.endpointIdx]
	if intervalEnd > limit {
		us.pos = limit
		*end = limit
		return true
	}
	*end = intervalEnd
	us.endpointIdx++
	if us.endpointIdx.Finished(us.endpoints) {
		us.pos = PosTypeMax
	} else {
		us.pos = us.endpoints[us.endpointIdx]
		us.endpointIdx++
	}
	return true
}
