package interval

import (
	"bufio"
	"fmt"
	"io"
	"math"
	"sort"
	"strconv"
	"strings"

	"github.com/biogo/hts/sam"
	"github.com/grailbio/base/file"
	"github.com/grailbio/base/fileio"
	"github.com/grailbio/base/log"
	gunsafe "github.com/grailbio/base/unsafe"
	"github.com/grailbio/base/vcontext"
	"github.com/klauspost/compress/gzip"
)

// getTokens identifies up to the first len(tokens) tokens from curLine,
// returning the number of tokens saved.  Any (group of) characters <= ' ' is
// treated as a delimiter.
//
// A variant of this function which scrapes an arbitrary subset of the columns
// will probably be added to base/simd; that's useful for processing VCF-like
// files (but too heavyweight for the first three columns of a BED).
func getTokens(tokens [][]byte, curLine []byte) int {
	posEnd := 0
	lineLen := len(curLine)
	for tokenIdx := range tokens {
		// These simple loops are better than simd.FirstGreater(src, ' ', startPos)
		// and simd.FirstLeq(src, ' ', startPos) when length <20 tokens are
		// expected.  They are also better than any of the standard library
		// string-split functions.
		// Unfortunately, the compiler currently does not inline any function with
		// a loop no matter how trivial, so we can't justify making these 5-line
		// for loops functions of their own.
		//
		// We may want to tweak this a bit to minimize the number of unnecessary
		// bounds-checks, but wait for Go 1.11 since that contains its own BCE
		// optimizations.
		pos := posEnd
		for ; pos != lineLen; pos++ {
			if curLine[pos] > ' ' {
				break
			}
		}
		if pos == lineLen {
			return tokenIdx
		}
		posEnd = pos
		for ; posEnd != lineLen; posEnd++ {
			if curLine[posEnd] <= ' ' {
				break
			}
		}
		tokens[tokenIdx] = curLine[pos:posEnd]
	}
	return len(tokens)
}

// NewBEDOpts defines behavior of this package's BED-loading function(s).
type NewBEDOpts struct {
	// SAMHeader enables ID-based lookup.  (This is more convenient than
	// string-based lookup when using gbam.Shard.)
	SAMHeader *sam.Header
	// Invert causes the complement of the interval-union to be returned.  The
	// complement extends down to position -1 at the beginning of each
	// chromosome, and currently 2^31 - 2 inclusive at the end.  If SAMHeader is
	// provided, any chromosome mentioned in the SAMHeader but entirely absent
	// from the BED will be fully included.  Otherwise, only the chromosomes
	// mentioned in the BED file are included.  (A single empty interval
	// qualifies as a "mention" for the latter purpose.)
	Invert bool
	// OneBasedInput interprets the BED interval boundaries as one-based [start,
	// end] instead of the usual zero-based [start, end).
	OneBasedInput bool
}

// PosType is BEDUnion's coordinate type.
type PosType int32

const posTypeMax = math.MaxInt32

// searchPosType returns the index of x in a[], or the position where x would
// be inserted if x isn't in a (this could be len(a)).  It's exactly the same
// as sort.SearchInt(), except for PosType.
func searchPosType(a []PosType, x PosType) int {
	return sort.Search(len(a), func(i int) bool { return a[i] >= x })
}

// fwdsearchPosType checks a[idx], then a[idx + 1], then a[idx + 3], then
// a[idx + 7], etc., and then uses binary search to finish the job.  It's
// usually a better choice than searchPosType when iterating.
// (However, an inlined simple linear search may be better in practice.  Can
// benchmark later if it matters.)
func fwdsearchPosType(a []PosType, x PosType, idx int) int {
	nextIncr := 1
	startIdx := idx
	endIdx := len(a)
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
		midIdx := int(uint(startIdx+endIdx) >> 1)
		if a[midIdx] >= x {
			endIdx = midIdx
		} else {
			startIdx = midIdx + 1
		}
	}
	return startIdx
}

// BEDUnion is currently implemented as a collection of length-2N sequences,
// where N is the number of intervals, the (0-based) start position of the
// interval #k (numbering from zero) is in element [2k] and the end position is
// in element [2k+1], and the intervals are stored in increasing order.
// Advantages of this representation over a length-N sequence of {start, end}
// structs include simpler inversion code, and reuse of standard []int32 binary
// and similar search algorithms (which the compiler is more likely to optimize
// well).
type BEDUnion struct {
	// nameMap is a chromosome-keyed map with disjoint-interval-set values.
	// Always initialized.
	nameMap map[string]([]PosType)
	// idMap is an optional slice of disjoint-interval-sets, indexed by biogo
	// sam.Header reference ID.  It is only initialized if NewBEDUnion{FromPath}
	// was called with SAMHeader initialized.
	idMap [][]PosType
	// lastChrIntervals points to the disjoint-interval-set for the most recently
	// queried chromosome.  This is a minor performance optimization.
	lastChrIntervals []PosType
	// lastChrName is the name of the last queried-by-name chromosome.  If it's
	// nonempty, it must be in sync with lastChrIntervals.
	lastChrName string
	// lastChrID is the name of the last queried-by-ID chromosome.  If it's
	// nonnegative, it must be in sync with lastChrIntervals.
	lastChrID int
	// lastPosPlus1 is 1 plus the last spot-queried position.
	lastPosPlus1 PosType
	// lastIdx is searchPosType(lastChrIntervals, lastPosPlus1).  Cached to
	// accelerate sequential queries.
	lastIdx int
	// isSequential is true if all queries since the last chromosome change have
	// been in order of nondecreasing position.
	isSequential bool
}

// ContainsByID checks whether the (0-based) interval [pos, pos+1) is contained
// within the BEDUnion, where chromosome is specified by sam.Header ID.
func (u *BEDUnion) ContainsByID(chrID int, pos PosType) bool {
	posPlus1 := pos + 1
	if chrID != u.lastChrID {
		u.lastChrID = chrID
		// bugfix (27 Jul 2018): need to set lastChrName to either empty, or the
		// name of this chromosome.  Otherwise lastChrIntervals is out of sync if
		// the next query is by name.
		u.lastChrName = ""

		// just let this error out the usual way if the BEDUnion was not
		// initialized with ID info.
		u.lastChrIntervals = u.idMap[chrID]
		// Force use of searchPosType() on the first query for a contig.
		if u.lastChrIntervals == nil {
			return false
		}
		u.lastIdx = searchPosType(u.lastChrIntervals, posPlus1)
		u.lastPosPlus1 = posPlus1
		u.isSequential = true
		return u.lastIdx&1 == 1
	}
	if u.lastChrIntervals == nil {
		return false
	}
	if u.isSequential {
		if posPlus1 >= u.lastPosPlus1 {
			u.lastIdx = fwdsearchPosType(u.lastChrIntervals, posPlus1, u.lastIdx)
			u.lastPosPlus1 = posPlus1
			return u.lastIdx&1 == 1
		}
		u.isSequential = false
	}
	return searchPosType(u.lastChrIntervals, posPlus1)&1 == 1
}

// ContainsByName checks whether the (0-based) interval [pos, pos+1) is
// contained within the BEDUnion, where chromosome is specified by name.
func (u *BEDUnion) ContainsByName(chrName string, pos PosType) bool {
	posPlus1 := pos + 1
	if chrName != u.lastChrName {
		u.lastChrName = chrName
		u.lastChrID = -1
		u.lastChrIntervals = u.nameMap[chrName]
		// Force use of searchPosType() on the first query for a contig.
		if u.lastChrIntervals == nil {
			return false
		}
		u.lastIdx = searchPosType(u.lastChrIntervals, posPlus1)
		u.lastPosPlus1 = posPlus1
		u.isSequential = true
		return u.lastIdx&1 == 1
	}
	if u.lastChrIntervals == nil {
		return false
	}
	if u.isSequential {
		if posPlus1 >= u.lastPosPlus1 {
			u.lastIdx = fwdsearchPosType(u.lastChrIntervals, posPlus1, u.lastIdx)
			u.lastPosPlus1 = posPlus1
			return u.lastIdx&1 == 1
		}
		u.isSequential = false
	}
	return searchPosType(u.lastChrIntervals, posPlus1)&1 == 1
}

// Intersects checks whether the given contiguous possibly-multi-chromosome
// region intersects the interval set.  Chromosomes must be specified by ID.
// It panics if limitRefID:limitPos isn't after startRefID:startPos.
func (u *BEDUnion) Intersects(startRefID int, startPos PosType, limitRefID int, limitPos PosType) bool {
	// May want a variant of this which takes a single chromosome name.
	if startRefID > limitRefID {
		panic("internal error: BEDUnion.Intersects requires startRefID <= limitRefID")
	}
	if startChrIntervals := u.idMap[startRefID]; startChrIntervals != nil {
		idxStart := searchPosType(startChrIntervals, startPos+1)
		if startRefID < limitRefID {
			if idxStart < len(startChrIntervals) {
				return true
			}
		} else {
			if limitPos <= startPos {
				panic("internal error: BEDUnion.Intersects requires limitPos > startPos when startRefID == limitRefID")
			}
			if idxStart&1 == 1 {
				return true
			}
			return (idxStart != len(startChrIntervals)) && (limitPos > startChrIntervals[idxStart])
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
	if limitChrIntervals := u.idMap[limitRefID]; limitChrIntervals != nil {
		return limitChrIntervals[0] < limitPos
	}
	return false
}

func initBEDUnion() (bedUnion BEDUnion) {
	bedUnion.nameMap = make(map[string]([]PosType))
	bedUnion.lastChrName = ""
	bedUnion.lastChrID = -1
	return
}

func (u *BEDUnion) nameToIDData(header *sam.Header, invert bool) {
	samRefs := header.Refs()
	nRef := len(samRefs)
	u.idMap = make([][]PosType, nRef)
	for refID, ref := range samRefs {
		// Validate ID property.  (Replace this with a comment if this is
		// guaranteed; I wasn't able to quickly find code in biogo/hts/sam which
		// made this clear one way or the other.)
		if refID != ref.ID() {
			panic("internal error: sam.header ref.ID != array position")
		}
		refName := ref.Name()
		chrIntervals := u.nameMap[refName]
		if chrIntervals != nil {
			u.idMap[refID] = chrIntervals
		} else if invert {
			u.idMap[refID] = []PosType{-1, posTypeMax}
		}
	}
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
	prevChr := ""
	totBases := 0
	var prevStart, prevEnd PosType
	var chrIntervals []PosType
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

		curChr := tokens[0]
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
		if (parsedEnd < parsedStart) || (parsedEnd >= posTypeMax) {
			err = fmt.Errorf("interval.scanBEDUnion: invalid coordinate pair on line %d", lineIdx)
			return
		}
		end := PosType(parsedEnd)
		if prevChr != gunsafe.BytesToString(curChr) {
			if prevChr != "" {
				// Save last interval, add to map.
				if prevEnd != -1 {
					chrIntervals = append(chrIntervals, prevStart, prevEnd)
				}
				if opts.Invert {
					chrIntervals = append(chrIntervals, posTypeMax)
				}
				bedUnion.nameMap[prevChr] = chrIntervals
			}
			// bugfix (12 Jul 2018): Must create a copy of curChr contents, since it
			// refers to bytes on curLine that will be overwritten soon.
			// Make a full heap copy instead of reusing a prevChrBytes []byte buffer,
			// since this needs to persist as a map key.
			prevChr = string(curChr)
			if _, found := bedUnion.nameMap[prevChr]; found {
				err = fmt.Errorf("interval.scanBEDUnion: unsorted input (split chromosome %v)", curChr)
				return
			}
			chrIntervals = []PosType{}
			if opts.Invert {
				chrIntervals = append(chrIntervals, -1)
			}
			if end == start {
				// Distinguish between 'mentioned' chromosomes without any overlapping
				// bases and unmentioned chromosomes.
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
			chrIntervals = append(chrIntervals, prevStart, prevEnd)
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
	if prevChr != "" {
		chrIntervals = append(chrIntervals, prevStart, prevEnd)
		if opts.Invert {
			chrIntervals = append(chrIntervals, posTypeMax)
		}
		bedUnion.nameMap[prevChr] = chrIntervals
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
	defer func() {
		if cerr := infile.Close(ctx); cerr != nil && err == nil {
			err = cerr
		}
	}()
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
	ChrName string
	Start0  PosType
	End     PosType
}

// ParseRegionString parses a region string of one of the forms
//   [contig ID]:[1-based first pos]-[last pos]
//   [contig ID]:[1-based pos]
//   [contig ID]
// returning a contig ID and 0-based interval boundaries.  The interval
// [0, posTypeMax - 1] is returned if there is no positional restriction.
func ParseRegionString(region string) (result Entry, err error) {
	if len(region) == 0 {
		err = fmt.Errorf("interval.ParseRegionString: empty region string")
		return
	}
	colonPos := strings.IndexByte(region, ':')
	if colonPos == -1 {
		result.ChrName = region
		result.Start0 = 0
		result.End = posTypeMax - 1
		return
	}
	if colonPos == 0 {
		err = fmt.Errorf("interval.ParseRegionString: empty contig ID")
		return
	}
	result.ChrName = region[0:colonPos]
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
	// We may as well prohibit end0 == posTypeMax so that the interval-array
	// is guaranteed to contain no repeats.  This means ParseInt(., 10, 32)
	// doesn't quite do the right thing, so Atoi is used above.
	if end0 <= start1 || end0 >= posTypeMax {
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
	prevChr := ""
	var prevStart, prevEnd PosType
	var chrIntervals []PosType
	for _, entry := range entries {
		curChr := entry.ChrName
		if entry.Start0 < 0 {
			err = fmt.Errorf("interval.NewBEDUnionFromEntries: negative start coordinate")
			return
		}

		if (entry.End < entry.Start0) || (entry.End >= posTypeMax) {
			err = fmt.Errorf("interval.NewBEDUnionFromEntry: invalid coordinate pair [%d, %d)", entry.Start0, entry.End)
			return
		}
		if prevChr != curChr {
			if prevChr != "" {
				// Save last interval, add to map.
				if prevEnd != -1 {
					chrIntervals = append(chrIntervals, prevStart, prevEnd)
				}
				if opts.Invert {
					chrIntervals = append(chrIntervals, posTypeMax)
				}
				bedUnion.nameMap[prevChr] = chrIntervals
			}
			prevChr = curChr
			if _, found := bedUnion.nameMap[prevChr]; found {
				err = fmt.Errorf("interval.NewBEDUnionFromEntry: unsorted input (split chromosome %v)", curChr)
				return
			}
			chrIntervals = []PosType{}
			if opts.Invert {
				chrIntervals = append(chrIntervals, -1)
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
				chrIntervals = append(chrIntervals, prevStart, prevEnd)
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
	if prevChr != "" {
		if prevEnd != -1 {
			chrIntervals = append(chrIntervals, prevStart, prevEnd)
		}
		if opts.Invert {
			chrIntervals = append(chrIntervals, posTypeMax)
		}
		bedUnion.nameMap[prevChr] = chrIntervals
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
	bedUnion.lastChrIntervals = nil
	bedUnion.lastChrName = ""
	bedUnion.lastChrID = -1
	return
}
