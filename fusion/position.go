package fusion

// Pos is the position in a read or a fragment (i.e., paired reads).
//
// When Pos refers to the position in a read, it is simply the zero-based index
// within the read sequence.
//
// When Pos refers to the position in a fragment, it can be either position in
// R1 or in R2.  A position in R2 has a value >= r2PosOffset. To get the actual
// index within the R2 sequence, subtract r2PosOffset from the value.
type Pos int64

// PosRange is a half-open range [start, end).
//
// INVARIANT: Start and End never cross a R1/R2 boundary.
type PosRange struct{ Start, End Pos }

// Equal checks if the two ranges are identical.
func (r PosRange) Equal(other PosRange) bool {
	return r.Start == other.Start && r.End == other.End
}

func (r PosRange) span() int {
	return posSpan(r.End, r.Start)
}

// CrossReadPosRange is the same as PosRange, except the range may cross a R1/R2
// boundary.
type CrossReadPosRange struct{ Start, End Pos }

// newCrossReadPosRange creates a new CrossReadPosRange
//
// REQUIRES: start <= end
func newCrossReadPosRange(start, end Pos) CrossReadPosRange {
	if end < start {
		panic("inverted range")
	}
	return CrossReadPosRange{start, end}
}

// Equal checks if the two ranges are identical.
func (r CrossReadPosRange) Equal(other CrossReadPosRange) bool {
	return r.Start == other.Start && r.End == other.End
}

// ReadType defines the read type (R1 or R2) for a paired fragment.
type ReadType uint8

const (
	// R1 means the read is either raw read from R1 fastq file, or it is a result
	// of stitching R1 and R2.
	R1 ReadType = iota
	// R2 means the read is from R2 fastq file. Note: when the read pair could be
	// stitched, the combined result will be stored in Fragment.R1Seq and
	// R2Seq. will be empty.
	R2
)

const r2PosOffset = Pos(1000000000)

// newR2Pos creates a Pos that refers to the idx'th base in R2.
func newR2Pos(idx Pos) Pos {
	if idx.ReadType() != R1 {
		panic(idx)
	}
	return idx + r2PosOffset
}

// R2Off returns the offset of this position from the start of R2.
//
// REQUIRES: pos.Type()==R2.
func (pos Pos) R2Off() int {
	if pos < r2PosOffset {
		panic(pos)
	}
	return int(pos - r2PosOffset)
}

// ReadType returns the type of read.
func (pos Pos) ReadType() ReadType {
	if pos >= r2PosOffset {
		return R2
	}
	return R1
}

func (r PosRange) readType() ReadType {
	typ := r.Start.ReadType()
	if typ != r.End.ReadType() {
		panic("cross-read range")
	}
	return typ
}

// newPosRange creates a new PosRange.
//
// REQUIRES: [start,end) does not cross a R1/R2 boundary.
func newPosRange(start, end Pos) PosRange {
	if start.ReadType() != end.ReadType() {
		panic("cross-read range")
	}
	if end < start {
		panic("inverted range")
	}
	return PosRange{start, end}
}

func posSpan(end, start Pos) int {
	if start.ReadType() != end.ReadType() {
		panic("cross-read range")
	}
	if end < start {
		panic("inverted range")
	}
	return int(end - start)
}

func maxPos(p1, p2 Pos) Pos {
	if p1 > p2 {
		return p1
	}
	return p2
}

func minPos(p1, p2 Pos) Pos {
	if p1 < p2 {
		return p1
	}
	return p2
}
