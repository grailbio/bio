package circular

import "math/bits"

// NextExp2 returns the next power of 2 strictly greater than x.  (Useful when
// setting circular buffer size.)
func NextExp2(x int) int {
	log2 := 63 - bits.LeadingZeros64(uint64(x))
	return 2 << uint32(log2)
}
