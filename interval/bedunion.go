package interval

import (
	"bufio"
	"fmt"
	"io"
	"strconv"
	"strings"

	"github.com/grailbio/base/file"
	"github.com/grailbio/base/fileio"
	"github.com/grailbio/base/log"
	gunsafe "github.com/grailbio/base/unsafe"
	"github.com/grailbio/base/vcontext"
	"github.com/grailbio/hts/sam"
	"github.com/klauspost/compress/gzip"
)

// NewBEDOpts defines behavior of this package's BED-loading function(s).
type NewBEDOpts struct {
	// SAMHeader enables ID-based lookup.  (This is more convenient than
	// string-based lookup when using gbam.Shard.)
	SAMHeader *sam.Header
	// Invert causes the complement of the interval-union to be returned.  The
	// complement extends down to position -1 at the beginning of each
	// reference, and currently 2^31 - 2 inclusive at the end.  If SAMHeader is
	// provided, any reference mentioned in the SAMHeader but entirely absent
	// from the BED will be fully included.  Otherwise, only the references
	// mentioned in the BED file are included.  (A single empty interval
	// qualifies as a "mention" for the latter purpose.)
	Invert bool
	// OneBasedInput interprets the BED interval boundaries as one-based [start,
	// end] instead of the usual zero-based [start, end).
	OneBasedInput bool
}

// intervalUnion is the internal representation of an interval-union on a
// single chromosome.
//
// As of this writing, it's a sorted slice of interval-endpoints.  More
// precisely, if there are N disjoint intervals on a chromosome, the
// corresponding slice has length 2N, the (0-based) start position of interval
// #k (numbering from zero) is element [2k] in the slice, and the end position
// of that interval is element [2k+1].
//
// Advantages of this representation over a length-N sequence of {start, end}
// structs include simpler inversion code, and reuse of standard []int32 binary
// and similar search algorithms (which the compiler is more likely to optimize
// well).  Since it's a broadly useful data structure, we expose it with the
// EndpointIndex type and the BEDUnion.OverlapByID(), .IntersectionByID(), and
// .EndpointsByID() methods below.
//
// However, we define this type to leave open the possibility of changing the
// internal representation in the future.  In the rest of this file, []PosType
// is used in code that's specific to the current interval-endpoints
// representation, while intervalUnion is used at a higher level of
// abstraction.
type intervalUnion []PosType

// newIntervalUnionFromEndpoints initializes an intervalUnion from a sorted
// slice of interval-endpoints.  (This function doesn't actually do anything
// for now, but it's defined so that tests will continue to work if we change
// intervalUnion.)
func newIntervalUnionFromEndpoints(endpoints []PosType) intervalUnion {
	return endpoints
}

// BEDUnion represents the set of positions covered by the intervals in a BED
// file, in a manner that supports efficient iteration and intersection
// queries.  (It does not keep track of interval names, etc.; just the final
// position set.)
type BEDUnion struct {
	// nameMap is a reference-name-keyed map with interval-union values.  Always
	// initialized.
	nameMap map[string](intervalUnion)
	// idMap is an optional slice of interval-unions, indexed by biogo sam.Header
	// reference ID.  It is only initialized if NewBEDUnion{FromPath} was called
	// with SAMHeader initialized.
	idMap []intervalUnion
	// RefNames maps reference IDs to names, when idMap is initialized.
	RefNames []string
	// lastRefIntervals points to the interval-union for the most recently
	// queried reference.  This is a minor performance optimization.
	lastRefIntervals intervalUnion
	// lastRefName is the name of the last queried-by-name reference.  If it's
	// nonempty, it must be in sync with lastRefIntervals.
	lastRefName string
	// lastRefID is the name of the last queried-by-ID reference.  If it's
	// nonnegative, it must be in sync with lastRefIntervals.
	lastRefID int
	// lastPosPlus1 is 1 plus the last spot-queried position.
	lastPosPlus1 PosType
	// lastEndpointIndex is SearchPosTypes(lastRefIntervals, lastPosPlus1).
	// Cached to accelerate sequential queries.
	lastEndpointIdx EndpointIndex
	// isSequential is true if all queries since the last reference change have
	// been in order of nondecreasing position.
	isSequential bool
}

// RefNameSet returns the set of reference names.  This works even when the
// BEDUnion was constructed without a *sam.Header.
func (u *BEDUnion) RefNameSet() map[string]bool {
	names := make(map[string]bool)
	for k := range u.nameMap {
		names[k] = true
	}
	return names
}

// ContainsByID checks whether the (0-based) interval [pos, pos+1) is contained
// within the BEDUnion, where reference is specified by sam.Header ID.
func (u *BEDUnion) ContainsByID(refID int, pos PosType) bool {
	posPlus1 := pos + 1
	if refID != u.lastRefID {
		u.lastRefID = refID
		// bugfix (27 Jul 2018): need to set lastRefName to either empty, or the
		// name of this reference.  Otherwise lastRefIntervals is out of sync if
		// the next query is by name.
		u.lastRefName = ""

		// just let this error out the usual way if the BEDUnion was not
		// initialized with ID info.
		u.lastRefIntervals = u.idMap[refID]
		// Force use of SearchPosTypes() on the first query for a contig.
		if u.lastRefIntervals == nil {
			return false
		}
		u.lastEndpointIdx = SearchPosTypes(u.lastRefIntervals, posPlus1)
		u.lastPosPlus1 = posPlus1
		u.isSequential = true
		return u.lastEndpointIdx.Contained()
	}
	if u.lastRefIntervals == nil {
		return false
	}
	if u.isSequential {
		if posPlus1 >= u.lastPosPlus1 {
			u.lastEndpointIdx = ExpsearchPosType(u.lastRefIntervals, posPlus1, u.lastEndpointIdx)
			u.lastPosPlus1 = posPlus1
			return u.lastEndpointIdx.Contained()
		}
		u.isSequential = false
	}
	return SearchPosTypes(u.lastRefIntervals, posPlus1).Contained()
}

// ContainsByName checks whether the (0-based) interval [pos, pos+1) is
// contained within the BEDUnion, where reference is specified by name.
func (u *BEDUnion) ContainsByName(refName string, pos PosType) bool {
	posPlus1 := pos + 1
	if refName != u.lastRefName {
		u.lastRefName = refName
		u.lastRefID = -1
		u.lastRefIntervals = u.nameMap[refName]
		// Force use of SearchPosTypes() on the first query for a contig.
		if u.lastRefIntervals == nil {
			return false
		}
		u.lastEndpointIdx = SearchPosTypes(u.lastRefIntervals, posPlus1)
		u.lastPosPlus1 = posPlus1
		u.isSequential = true
		return u.lastEndpointIdx.Contained()
	}
	if u.lastRefIntervals == nil {
		return false
	}
	if u.isSequential {
		if posPlus1 >= u.lastPosPlus1 {
			u.lastEndpointIdx = ExpsearchPosType(u.lastRefIntervals, posPlus1, u.lastEndpointIdx)
			u.lastPosPlus1 = posPlus1
			return u.lastEndpointIdx.Contained()
		}
		u.isSequential = false
	}
	return SearchPosTypes(u.lastRefIntervals, posPlus1).Contained()
}

// Intersects checks whether the given contiguous possibly-multi-reference
// region intersects the interval set.  References must be specified by ID.
// It panics if limitRefID:limitPos isn't after startRefID:startPos.
func (u *BEDUnion) Intersects(startRefID int, startPos PosType, limitRefID int, limitPos PosType) bool {
	if startRefID > limitRefID {
		panic("BEDUnion.Intersects: startRefID <= limitRefID required")
	}
	if startRefIntervals := u.idMap[startRefID]; startRefIntervals != nil {
		idxStart := SearchPosTypes(startRefIntervals, startPos+1)
		if startRefID < limitRefID {
			if !idxStart.Finished(startRefIntervals) {
				return true
			}
		} else {
			if limitPos <= startPos {
				panic("BEDUnion.Intersects: limitPos > startPos required when startRefID == limitRefID")
			}
			if idxStart.Contained() {
				return true
			}
			return (!idxStart.Finished(startRefIntervals)) && (limitPos > startRefIntervals[idxStart])
		}
	}
	if startRefID == limitRefID {
		return false
	}
	for refID := startRefID + 1; refID < limitRefID; refID++ {
		if u.idMap[refID] != nil {
			return true
		}
	}
	if limitRefIntervals := u.idMap[limitRefID]; limitRefIntervals != nil {
		return limitRefIntervals[0] < limitPos
	}
	return false
}

// IntersectsByID checks whether the given nonempty single-reference [startPos,
// limitPos) interval intersects the BEDUnion, where the reference is specified
// by ID.
func (u *BEDUnion) IntersectsByID(refID int, startPos PosType, limitPos PosType) bool {
	// Implementation has more in common with ContainsByID(); this is optimized
	// for a sequence of queries with nondecreasing startPos.
	// (add ByName version of this function when necessary)
	posPlus1 := startPos + 1
	if refID != u.lastRefID {
		u.lastRefID = refID
		u.lastRefName = ""
		u.lastRefIntervals = u.idMap[refID]
		if u.lastRefIntervals == nil {
			return false
		}
		u.lastEndpointIdx = SearchPosTypes(u.lastRefIntervals, posPlus1)
		u.isSequential = true
	} else {
		if u.lastRefIntervals == nil {
			return false
		}
		if posPlus1 < u.lastPosPlus1 {
			u.isSequential = false
		}
		if !u.isSequential {
			startIdx := SearchPosTypes(u.lastRefIntervals, posPlus1)
			if startIdx.Contained() {
				return true
			}
			return (!startIdx.Finished(u.lastRefIntervals)) && (u.lastRefIntervals[startIdx] < limitPos)
		}
		u.lastEndpointIdx = ExpsearchPosType(u.lastRefIntervals, posPlus1, u.lastEndpointIdx)
	}
	u.lastPosPlus1 = posPlus1
	if u.lastEndpointIdx.Contained() {
		return true
	}
	return (!u.lastEndpointIdx.Finished(u.lastRefIntervals)) && (u.lastRefIntervals[u.lastEndpointIdx] < limitPos)
}

// OverlapByID is a low-level function which returns a []PosType describing all
// the BED intervals overlapping [startPos, limitPos) on the given reference,
// where the reference is specified by ID.  For example, assuming positive
// startPos, the BED interval [0, startPos) would not be returned, while [0,
// startPos + 1) would.
// IMPORTANT: the return value must be treated as read-only; it may alias an
// internal data structure.
// The return slice is arranged such that the value at index 2n is the start of
// overlapping interval n, and the value at index (2n+1) is the end of that
// overlapping interval.
func (u *BEDUnion) OverlapByID(refID int, startPos, limitPos PosType) []PosType {
	posPlus1 := startPos + 1
	if refID != u.lastRefID {
		u.lastRefID = refID
		u.lastRefName = ""
		u.lastRefIntervals = u.idMap[refID]
		if u.lastRefIntervals == nil {
			return nil
		}
		u.lastEndpointIdx = SearchPosTypes(u.lastRefIntervals, posPlus1)
		u.isSequential = true
	} else {
		if u.lastRefIntervals == nil {
			return nil
		}
		if posPlus1 < u.lastPosPlus1 {
			u.isSequential = false
		}
		if !u.isSequential {
			startIdx := SearchPosTypes(u.lastRefIntervals, posPlus1)
			// If start of query is inside an interval, get the beginning of that
			// interval by rounding down to the next even index.
			resultStart := startIdx.Begin()
			// If end of query is inside an interval, get the end of that interval by
			// rounding up to the next even index.
			resultEnd := (SearchPosTypes(u.lastRefIntervals, limitPos) + 1).Begin()
			return u.lastRefIntervals[resultStart:resultEnd]
		}
		u.lastEndpointIdx = ExpsearchPosType(u.lastRefIntervals, posPlus1, u.lastEndpointIdx)
	}
	u.lastPosPlus1 = posPlus1
	resultStart := u.lastEndpointIdx.Begin()
	resultEnd := (ExpsearchPosType(u.lastRefIntervals, limitPos, u.lastEndpointIdx) + 1).Begin()
	return u.lastRefIntervals[resultStart:resultEnd]
}

// IntersectionByID is similar to OverlapByID, except that the intersection
// between [startPos, limitPos) and the BEDUnion is returned (so the first
// and/or the last interval may be truncated).  The return value must be
// treated as read-only.
func (u *BEDUnion) IntersectionByID(refID int, startPos, limitPos PosType) []PosType {
	result := u.OverlapByID(refID, startPos, limitPos)
	if (len(result) == 0) || ((result[0] >= startPos) && (result[len(result)-1] <= limitPos)) {
		return result
	}
	// Can't just alias an internal slice, since we must truncate one of the
	// intervals.
	truncatedResult := append([]PosType(nil), result...)
	if truncatedResult[0] < startPos {
		truncatedResult[0] = startPos
	}
	if truncatedResult[len(truncatedResult)-1] > limitPos {
		truncatedResult[len(truncatedResult)-1] = limitPos
	}
	return truncatedResult
}

// EndpointsByID returns the sorted interval-endpoint slice for the reference
// with the given ID.  This must be treated as read-only; it may alias an
// internal data structure.
func (u *BEDUnion) EndpointsByID(refID int) []PosType {
	return u.idMap[refID]
}

// EndpointsByName has the same functionality as EndpointsByID, without
// requiring a *sam.Header at BEDUnion construction time.
func (u *BEDUnion) EndpointsByName(refName string) []PosType {
	return u.nameMap[refName]
}

func initBEDUnion() (bedUnion BEDUnion) {
	bedUnion.nameMap = make(map[string](intervalUnion))
	bedUnion.lastRefName = ""
	bedUnion.lastRefID = -1
	return
}

func (u *BEDUnion) nameToIDData(header *sam.Header, invert bool) {
	samRefs := header.Refs()
	nRef := len(samRefs)
	u.idMap = make([]intervalUnion, nRef)
	u.RefNames = make([]string, nRef)
	for refID, ref := range samRefs {
		// Validate ID property.  (Replace this with a comment if this is
		// guaranteed; I wasn't able to quickly find code in biogo/hts/sam which
		// made this clear one way or the other.)
		if refID != ref.ID() {
			panic("BEDUnion.nameToIDData: sam.header ref.ID != array position")
		}
		refName := ref.Name()
		refIntervals := u.nameMap[refName]
		u.RefNames[refID] = refName
		if refIntervals != nil {
			u.idMap[refID] = refIntervals
		} else if invert {
			u.idMap[refID] = newIntervalUnionFromEndpoints([]PosType{-1, PosTypeMax})
		}
	}
}

// getTokens identifies up to the first len(tokens) tokens from curLine,
// returning the number of tokens saved.  Any (group of) characters <= ' ' is
// treated as a delimiter.
//
// (This is too low-level to make sense alongside the base/tsv reader
// interface, but it comes up often enough that I think we should move it to a
// more central location at some point.)
func getTokens(tokens [][]byte, curLine []byte) int {
	// Make this unsigned to help the bounds-check-eliminator.
	// (As of Go 1.12, if this is signed, the compiler fails to prove that pos
	// and posEnd are never negative, so the signed pos >= lineLen check is
	// followed by an unsigned pos >= len(curLine) bounds-check when accessing
	// curLine[pos].  But if it's unsigned, the checks are identical and the
	// second one is removed.)
	posEnd := uint(0)
	lineLen := uint(len(curLine))
	for tokenIdx := range tokens {
		// These simple loops are better than simd.FirstGreater(src, ' ', startPos)
		// and simd.FirstLeq(src, ' ', startPos) when length <20 tokens are
		// expected.  They are also better than any of the standard library
		// string-split functions.
		// Unfortunately, the compiler currently does not inline any function with
		// a loop no matter how trivial, so we can't justify making these 5-line
		// for loops functions of their own.
		pos := posEnd
		for ; ; pos++ {
			if pos >= lineLen {
				return tokenIdx
			}
			if curLine[pos] > ' ' {
				break
			}
		}
		for posEnd = pos; posEnd < lineLen; posEnd++ {
			if curLine[posEnd] <= ' ' {
				break
			}
		}
		// 3 spurious bounds-checks remain here, but that's a much smaller deal
		// than bounds-checks in the two loops above.
		tokens[tokenIdx] = curLine[pos:posEnd]
	}
	return len(tokens)
}

func scanBEDUnion(scanner *bufio.Scanner, opts NewBEDOpts) (bedUnion BEDUnion, err error) {
	bedUnion = initBEDUnion()

	var startSubtract int
	if opts.OneBasedInput {
		startSubtract++
	}

	// This could also be inside the for loop; minor tradeoff between extra
	// zero-reinitialization and positive side effects of better locality.
	var tokens [3][]byte

	lineIdx := 0
	prevRef := ""
	totBases := 0
	var prevStart, prevEnd PosType
	var refIntervals intervalUnion
	for scanner.Scan() {
		lineIdx++
		// Originally had a scanner.Text() call, since I'll take immutability
		// enforcement where I can get it... but turns out Text() allocates and
		// Bytes() does not?!  Sigh.
		// (Update: gunsafe.BytesToString should only be used in
		// very-limited-scope/lifetime scenarios; otherwise you end up fighting
		// against the language re: string copies and the like.  In particular,
		// making curLine an array of strings proved to be error-prone; better to
		// have e.g. a separate instance of gunsafe.BytesToString for each
		// strconv.Atoi() call despite the extra verbosity.)
		curLine := scanner.Bytes()
		nToken := getTokens(tokens[:], curLine)
		if nToken != 3 {
			if nToken == 0 {
				continue
			}
			err = fmt.Errorf("interval.scanBEDUnion: line %d has fewer tokens than expected", lineIdx)
			return
		}

		curRef := tokens[0]
		var parsedStart int
		if parsedStart, err = strconv.Atoi(gunsafe.BytesToString(tokens[1])); err != nil {
			return
		}
		parsedStart -= startSubtract
		if parsedStart < 0 {
			err = fmt.Errorf("interval.scanBEDUnion: negative start coordinate %v on line %d", tokens[1], lineIdx)
			return
		}
		start := PosType(parsedStart)

		var parsedEnd int
		if parsedEnd, err = strconv.Atoi(gunsafe.BytesToString(tokens[2])); err != nil {
			return
		}
		if (parsedEnd < parsedStart) || (parsedEnd >= PosTypeMax) {
			err = fmt.Errorf("interval.scanBEDUnion: invalid coordinate pair on line %d", lineIdx)
			return
		}
		end := PosType(parsedEnd)
		if prevRef != gunsafe.BytesToString(curRef) {
			if prevRef != "" {
				// Save last interval, add to map.
				if prevEnd != -1 {
					refIntervals = append(refIntervals, prevStart, prevEnd)
				}
				if opts.Invert {
					refIntervals = append(refIntervals, PosTypeMax)
				}
				bedUnion.nameMap[prevRef] = refIntervals
			}
			// bugfix (12 Jul 2018): Must create a copy of curRef contents, since it
			// refers to bytes on curLine that will be overwritten soon.
			// Make a full heap copy instead of reusing a prevRefBytes []byte buffer,
			// since this needs to persist as a map key.
			prevRef = string(curRef)
			if _, found := bedUnion.nameMap[prevRef]; found {
				err = fmt.Errorf("interval.scanBEDUnion: unsorted input (split reference %v)", curRef)
				return
			}
			refIntervals = intervalUnion{}
			if opts.Invert {
				refIntervals = append(refIntervals, -1)
			}
			if end == start {
				// Distinguish between 'mentioned' references without any overlapping
				// bases and unmentioned references.
				prevStart = -1
				prevEnd = -1
			} else {
				prevStart = start
				prevEnd = end
			}
			totBases += int(end - start)
			continue
		}
		if end == start {
			continue
		}
		if start > prevEnd {
			// New interval doesn't overlap previous one, so we can save the previous
			// one.
			refIntervals = append(refIntervals, prevStart, prevEnd)
			prevStart = start
			prevEnd = end
			totBases += int(end - start)
		} else {
			if start < prevStart {
				err = fmt.Errorf("interval.scanBEDUnion: unsorted input")
				return
			}
			// Intervals overlap, merge them.
			if end > prevEnd {
				totBases += int(end - prevEnd)
				prevEnd = end
			}
		}
	}
	if err = scanner.Err(); err != nil {
		return
	}
	log.Printf("BED loaded, %d base(s) covered.\n", totBases)
	if prevRef != "" {
		refIntervals = append(refIntervals, prevStart, prevEnd)
		if opts.Invert {
			refIntervals = append(refIntervals, PosTypeMax)
		}
		bedUnion.nameMap[prevRef] = refIntervals
	}
	return
}

// NewBEDUnion loads just the intervals from a sorted (by first coordinate)
// interval-BED, merging touching/overlapping intervals and eliminating empty
// ones in the process.  A BEDUnion is returned.
func NewBEDUnion(reader io.Reader, opts NewBEDOpts) (bedUnion BEDUnion, err error) {
	// Note that Scanner does not handle very long lines unless we specify an
	// adequate buffer size in advance; it does not auto-resize.
	// Shouldn't matter for BED files, though.
	scanner := bufio.NewScanner(reader)

	if bedUnion, err = scanBEDUnion(scanner, opts); err != nil {
		return
	}

	if opts.SAMHeader != nil {
		bedUnion.nameToIDData(opts.SAMHeader, opts.Invert)
	}
	return
}

// NewBEDUnionFromPath is a wrapper for NewBEDUnion that takes a path instead
// of an io.Reader.
func NewBEDUnionFromPath(path string, opts NewBEDOpts) (bedUnion BEDUnion, err error) {
	ctx := vcontext.Background()
	var infile file.File
	if infile, err = file.Open(ctx, path); err != nil {
		return
	}
	defer file.CloseAndReport(ctx, infile, &err)
	reader := io.Reader(infile.Reader(ctx))
	switch fileio.DetermineType(path) {
	case fileio.Gzip:
		if reader, err = gzip.NewReader(reader); err != nil {
			return
		}
	}
	return NewBEDUnion(reader, opts)
}

// Entry represents a single interval, with 0-based coordinates.
type Entry struct {
	RefName string
	Start0  PosType
	End     PosType
}

// ParseRegionString parses a region string of one of the forms
//   [contig ID]:[1-based first pos]-[last pos]
//   [contig ID]:[1-based pos]
//   [contig ID]
// returning a contig ID and 0-based interval boundaries.  The interval
// [0, PosTypeMax - 1] is returned if there is no positional restriction.
func ParseRegionString(region string) (result Entry, err error) {
	if len(region) == 0 {
		err = fmt.Errorf("interval.ParseRegionString: empty region string")
		return
	}
	colonPos := strings.IndexByte(region, ':')
	if colonPos == -1 {
		result.RefName = region
		result.Start0 = 0
		result.End = PosTypeMax - 1
		return
	}
	if colonPos == 0 {
		err = fmt.Errorf("interval.ParseRegionString: empty contig ID")
		return
	}
	result.RefName = region[0:colonPos]
	rangeStr := region[colonPos+1:]
	dashPos := strings.IndexByte(rangeStr, '-')
	if dashPos == -1 {
		var pos1 int64
		// Specify base for now, but could change to 0 as long as all the other
		// strconv.Atoi calls are replaced.
		if pos1, err = strconv.ParseInt(rangeStr, 10, 32); err != nil {
			return
		}
		if pos1 <= 0 {
			err = fmt.Errorf("interval.ParseRegionString: position %v in region string out of range", rangeStr)
			return
		}
		result.Start0 = PosType(pos1 - 1)
		result.End = PosType(pos1)
		return
	}
	start1Str := rangeStr[:dashPos]
	endStr := rangeStr[dashPos+1:]
	var start1 int
	if start1, err = strconv.Atoi(start1Str); err != nil {
		return
	}
	if start1 <= 0 {
		err = fmt.Errorf("interval.ParseRegionString: position %v in region string out of range", start1Str)
		return
	}
	var end0 int
	if end0, err = strconv.Atoi(endStr); err != nil {
		return
	}
	// We may as well prohibit end0 == PosTypeMax so that the interval-array
	// is guaranteed to contain no repeats.  This means ParseInt(., 10, 32)
	// doesn't quite do the right thing, so Atoi is used above.
	if end0 <= start1 || end0 >= PosTypeMax {
		err = fmt.Errorf("interval.ParseRegionString: invalid range string %v", rangeStr)
		return
	}
	result.Start0 = PosType(start1 - 1)
	result.End = PosType(end0)
	return
}

// NewBEDUnionFromEntries initializes a BEDUnion from a sorted []Entry.
// This ignores opts.OneBasedInput, since start0 is defined to be zero-based.
func NewBEDUnionFromEntries(entries []Entry, opts NewBEDOpts) (bedUnion BEDUnion, err error) {
	bedUnion = initBEDUnion()
	prevRef := ""
	var prevStart, prevEnd PosType
	var refIntervals intervalUnion
	for _, entry := range entries {
		curRef := entry.RefName
		if entry.Start0 < 0 {
			err = fmt.Errorf("interval.NewBEDUnionFromEntries: negative start coordinate")
			return
		}

		if (entry.End < entry.Start0) || (entry.End >= PosTypeMax) {
			err = fmt.Errorf("interval.NewBEDUnionFromEntry: invalid coordinate pair [%d, %d)", entry.Start0, entry.End)
			return
		}
		if prevRef != curRef {
			if prevRef != "" {
				// Save last interval, add to map.
				if prevEnd != -1 {
					refIntervals = append(refIntervals, prevStart, prevEnd)
				}
				if opts.Invert {
					refIntervals = append(refIntervals, PosTypeMax)
				}
				bedUnion.nameMap[prevRef] = refIntervals
			}
			prevRef = curRef
			if _, found := bedUnion.nameMap[prevRef]; found {
				err = fmt.Errorf("interval.NewBEDUnionFromEntry: unsorted input (split reference %v)", curRef)
				return
			}
			refIntervals = intervalUnion{}
			if opts.Invert {
				refIntervals = append(refIntervals, -1)
			}
			if entry.End == entry.Start0 {
				prevStart = -1
				prevEnd = -1
				continue
			}
			prevStart = entry.Start0
			prevEnd = entry.End
			continue
		}
		if entry.End == entry.Start0 {
			continue
		}
		if entry.Start0 > prevEnd {
			// New interval doesn't overlap previous one, so we can save the previous
			// one.
			if prevEnd != -1 {
				refIntervals = append(refIntervals, prevStart, prevEnd)
			}
			prevStart = entry.Start0
			prevEnd = entry.End
		} else {
			if entry.Start0 < prevStart {
				err = fmt.Errorf("interval.NewBEDUnionFromEntries: unsorted input")
				return
			}
			// Intervals overlap, merge them.
			if entry.End > prevEnd {
				prevEnd = entry.End
			}
		}
	}
	if prevRef != "" {
		if prevEnd != -1 {
			refIntervals = append(refIntervals, prevStart, prevEnd)
		}
		if opts.Invert {
			refIntervals = append(refIntervals, PosTypeMax)
		}
		bedUnion.nameMap[prevRef] = refIntervals
	}
	if opts.SAMHeader != nil {
		bedUnion.nameToIDData(opts.SAMHeader, opts.Invert)
	}
	return
}

// Clone returns a new BEDUnion which shares the interval set, but has its own
// search state.
func (u *BEDUnion) Clone() (bedUnion BEDUnion) {
	bedUnion.nameMap = u.nameMap
	bedUnion.idMap = u.idMap
	bedUnion.RefNames = u.RefNames
	bedUnion.lastRefIntervals = nil
	bedUnion.lastRefName = ""
	bedUnion.lastRefID = -1
	return
}

// Subset returns a new BEDUnion which describes the intersection between the
// original interval set and the given (possibly multi-reference) interval.
// The original BEDUnion must support ID-based lookup.  Parameters are handled
// in the same way as Intersects().
func (u *BEDUnion) Subset(startRefID int, startPos PosType, limitRefID int, limitPos PosType) (bedUnion BEDUnion) {
	bedUnion.nameMap = make(map[string](intervalUnion))
	nRef := len(u.idMap)
	bedUnion.idMap = make([]intervalUnion, nRef)
	bedUnion.RefNames = u.RefNames
	bedUnion.lastRefIntervals = nil
	bedUnion.lastRefName = ""
	bedUnion.lastRefID = -1

	// This code has a lot in common with Intersects().
	if startRefID > limitRefID {
		panic("BEDUnion.Subset: startRefID <= limitRefID required")
	}
	if startRefIntervals := u.idMap[startRefID]; startRefIntervals != nil {
		idxStart := NewEndpointIndex(startPos, startRefIntervals)
		// Two cases:
		// 1. startPos is strictly contained in an interval.  Then, we currently
		//    need to make a copy of part of startRefIntervals and change the first
		//    element.
		// 2. Otherwise, we can point the new idMap/nameMap entries at a subslice
		//    of startRefIntervals (unless startRefID == limitRefID and limitPos
		//    is in the middle of an interval).
		newAllocNeeded := idxStart.Contained() && (startRefIntervals[idxStart-1] != startPos)
		var refIntervalsSubset []PosType
		if startRefID < limitRefID {
			if !idxStart.Finished(startRefIntervals) {
				if newAllocNeeded {
					nElem := EndpointIndex(len(startRefIntervals)) - idxStart + 1
					refIntervalsSubset = make([]PosType, nElem)
					refIntervalsSubset[0] = startPos
					copy(refIntervalsSubset[1:], startRefIntervals[idxStart:])
				} else {
					// idxStart must be rounded down to next even integer in edge case.
					refIntervalsSubset = startRefIntervals[idxStart.Begin():]
				}
			}
		} else {
			if limitPos <= startPos {
				panic("BEDUnion.Subset: limitPos > startPos required when startRefID == limitRefID")
			}
			idxEnd := SearchPosTypes(startRefIntervals, limitPos+1)
			if idxEnd.Contained() && (startRefIntervals[idxEnd-1] != limitPos) {
				newAllocNeeded = true
			}
			if newAllocNeeded {
				idxStartCopy := idxStart | 1
				idxEndCopy := (idxEnd - 1) | 1
				nElem := idxEndCopy - idxStartCopy + 2
				refIntervalsSubset = make([]PosType, nElem)
				if idxStart.Contained() {
					refIntervalsSubset[0] = startPos
				} else {
					refIntervalsSubset[0] = startRefIntervals[idxStart]
				}
				copy(refIntervalsSubset[1:nElem-1], startRefIntervals[idxStartCopy:idxEndCopy])
				if idxEnd.Contained() {
					refIntervalsSubset[nElem-1] = limitPos
				} else {
					refIntervalsSubset[nElem-1] = startRefIntervals[idxEndCopy]
				}
			} else {
				idxStartFinal := idxStart.Begin()
				idxEndFinal := idxEnd.Begin()
				if idxStartFinal == idxEndFinal {
					return
				}
				refIntervalsSubset = startRefIntervals[idxStartFinal:idxEndFinal]
			}
		}
		if refIntervalsSubset != nil {
			bedUnion.idMap[startRefID] = refIntervalsSubset
			bedUnion.nameMap[u.RefNames[startRefID]] = refIntervalsSubset
		}
	}
	if startRefID == limitRefID {
		return
	}
	for refID := startRefID + 1; refID < limitRefID; refID++ {
		refIntervals := u.idMap[refID]
		bedUnion.idMap[refID] = refIntervals
		bedUnion.nameMap[u.RefNames[refID]] = refIntervals
	}
	if limitRefIntervals := u.idMap[limitRefID]; limitRefIntervals != nil {
		idxEnd := NewEndpointIndex(limitPos, limitRefIntervals)
		if idxEnd == 0 {
			return
		}
		var refIntervalsSubset []PosType
		if idxEnd.Contained() && (limitRefIntervals[idxEnd-1] != limitPos) {
			refIntervalsSubset = make([]PosType, idxEnd+1)
			copy(refIntervalsSubset[:idxEnd], limitRefIntervals[:idxEnd])
			refIntervalsSubset[idxEnd] = limitPos
		} else {
			refIntervalsSubset = limitRefIntervals[:idxEnd.Begin()]
		}
		bedUnion.idMap[limitRefID] = refIntervalsSubset
		bedUnion.nameMap[u.RefNames[limitRefID]] = refIntervalsSubset
	}
	return
}
