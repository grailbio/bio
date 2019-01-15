package fusion

import (
	"strings"

	"github.com/grailbio/base/log"
)

func addUMIToName(name, r1UMI, r2UMI string) string {
	b := strings.Builder{}
	var i int
	for i < len(name) {
		if name[i] == ' ' {
			break
		}
		b.WriteByte(name[i])
		i++
	}
	b.WriteByte(':')
	b.WriteString(r1UMI)
	b.WriteByte('+')
	b.WriteString(r2UMI)
	b.WriteString(name[i:])
	return b.String()
}

// MaybeRemoveUMI removes an UMI from the sequences and add add it to the name
// part, if the options prescribe such operations. It returns <new name, new r1
// seq, new r2seq>.
func MaybeRemoveUMI(name, r1Seq, r2Seq string, opts Opts) (string, string, string) {
	if opts.UMIInRead {
		// If one of the read lengths < 7, cannot obtain UMI sequences.
		if len(r1Seq) < 7 || len(r2Seq) < 7 {
			log.Error.Printf("UMI not found in %v %v", r1Seq, r2Seq)
			return name, "N", ""
		}
		if !opts.UMIInName {
			name = addUMIToName(name, r1Seq[:6], r2Seq[:6])
		}
		r1Seq = r1Seq[7:]
		r2Seq = r2Seq[7:]
	}
	return name, r1Seq, r2Seq
}

// RemoveLowComplexityReads check if r1Seq or r2Seq is low complexity, i.e., if
// the two most frequent nucleotide types dominate the sequence.  If so, it
// converts them to an empty string.
func RemoveLowComplexityReads(r1Seq, r2Seq string, stats *Stats, opts Opts) (newR1Seq, newR2Seq string) {
	// Check complexity of R1 and R2, respectively.
	isR1LC := IsLowComplexity(r1Seq, opts.LowComplexityFraction)
	isR2LC := IsLowComplexity(r2Seq, opts.LowComplexityFraction)
	if isR1LC {
		if isR2LC { // Both R1 and R2 are LC
			stats.LowComplexityReads2++
			return "N", ""
		}
		// R1 is LC, R2 is not
		stats.LowComplexityReads1++
		return reverseComplement(r2Seq), ""
	}
	if isR2LC { // R1 is not LC, R2 is
		stats.LowComplexityReads1++
		return r1Seq, ""
	}
	// Both R2 and R2 are not LC
	return r1Seq, reverseComplement(r2Seq)
}
