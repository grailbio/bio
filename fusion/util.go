package fusion

import (
	gunsafe "github.com/grailbio/base/unsafe"
	"github.com/grailbio/bio/biosimd"
)

func abs(x int) int {
	if x < 0 {
		return -x
	}
	return x
}

func max(x, y int) int {
	if x > y {
		return x
	}
	return y
}

// reverseComplement computes a reverse complement of the given DNA string.
//
// TODO(saito) Reverse in place, instead of creating a new string.
func reverseComplement(seq string) string {
	buf := make([]byte, len(seq))
	biosimd.ReverseComp8NoValidate(buf, gunsafe.StringToBytes(seq))
	return gunsafe.BytesToString(buf)
}

// acgtnIndex maps A, C, G, T to {0,1,2,3}. It maps other letters to 4.
var acgtnIndex [256]uint8

func init() {
	for i := range acgtnIndex {
		acgtnIndex[i] = 4
	}
	acgtnIndex['a'] = 0
	acgtnIndex['A'] = 0
	acgtnIndex['c'] = 1
	acgtnIndex['C'] = 1
	acgtnIndex['g'] = 2
	acgtnIndex['G'] = 2
	acgtnIndex['t'] = 3
	acgtnIndex['T'] = 3
}

func countACGTN(seq string) [5]int {
	var acgtnCounts [5]int
	for _, ch := range []byte(seq) {
		acgtnCounts[acgtnIndex[ch]]++
	}
	return acgtnCounts
}

// numUnknownBases returns the number of bases are not one of ACGT.
func numUnknownBases(seq string) int { return countACGTN(seq)[4] }

// IsLowComplexity returns true if input DNA is low complexity sequence,
// i.e. any two bases present at over lowComplexityFrac of the total sequence
// length.
func IsLowComplexity(seq string, lowComplexityFrac float64) bool {
	if len(seq) == 0 {
		return true
	}
	acgtnCounts := countACGTN(seq)
	max, max2 := -1, -1
	for _, c := range acgtnCounts {
		if max < c {
			max, max2 = c, max
		} else if max2 < c {
			max2 = c
		}
	}
	return float64(max+max2)/float64(len(seq)) > lowComplexityFrac
}
