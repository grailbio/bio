package fusion

import (
	"strings"

	"github.com/grailbio/base/log"
)

// Fragment is a union of two (unpaired) reads (R1 & R2). Created by the
// Stitcher.
type Fragment struct {
	// Fragment name. It's a copy of the R1 name from the fastq.
	//
	// Example: E00469:245:HHK5TCCXY:1:1101:19634:35080:GTATCT+AGCAAT
	Name string
	// R1Seq is the sequence from the R1 fastq.  When the R1 and R2 sequences are
	// found to have overlapping regions, they are stitched together, and R1Seq
	// will store the combined sequence and R2Seq will be empty.
	R1Seq string
	// R2Seq is the sequence from the R2 fastq. It is nonempty only when the
	// stitcher fails to stitch the R1 and R2 sequences.
	R2Seq string

	// kmers stores the set of kmers generated from R1Seq and R2Seq.
	kmers []kmersAtPos
}

// UMI extracts the UMI component from the Name field. Returns false if the name
// doesn't contain an UMI, or on any other error.
//
// Caution: This methods is very slow. Don't use it in a perf-critical path.
func (f *Fragment) UMI() string {
	name := f.Name
	sp := strings.IndexByte(f.Name, ' ')
	if sp >= 0 {
		name = name[0:sp]
	}
	parts := strings.Split(name, ":")
	if len(parts) == 0 {
		log.Panicf("Invalid name: %+v", f)
	}
	umi := parts[len(parts)-1]
	// Make sure that the string looks like an UMI.
	//
	// TODO(saito) Do a better check - e.g., the two sides of "+" must be of the
	// same length.
	for _, ch := range []byte(umi) {
		if p := strings.IndexByte("ACGTN+", ch); p < 0 {
			log.Panicf("Invalid umi: %+v", f)
		}
	}
	return umi
}

// SubSeq extracts part of the RNA sequence. The arg may cross the R1/R2
// boundary, in which case this function returns a suffix of R1 plus a prefix of
// R2.
func (r *Fragment) SubSeq(p CrossReadPosRange) string {
	if p.End.ReadType() == R1 {
		return r.R1Seq[p.Start:p.End]
	}
	if p.Start.ReadType() == R2 {
		return r.R2Seq[p.Start-r2PosOffset : p.End-r2PosOffset]
	}
	return r.R1Seq[p.Start:] + r.R2Seq[:p.End-r2PosOffset]
}

const infiniteHammingDistance = 1000000

// HammingDistance computes the hamming distance of sequences. If the sequences
// aren't of the same length, it returns a infiniteHammingDistance.
func (r *Fragment) HammingDistance(other Fragment) int {
	if len(r.R1Seq) != len(other.R1Seq) || len(r.R2Seq) != len(other.R2Seq) {
		return infiniteHammingDistance
	}
	return hammingDistance(r.R1Seq, other.R1Seq) + hammingDistance(r.R2Seq, other.R2Seq)
}
