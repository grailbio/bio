// Copyright 2019 GRAIL, Inc.  All rights reserved.
// Use of this source code is governed by the Apache-2.0
// license that can be found in the LICENSE file.

// +build !amd64 appengine

package biosimd

// FillFastqRecordBodyFromNibbles fills the body (defined as the last three
// lines) of a 4-line FASTQ record, given a packed 4-bit representation of the
// base+qual information and the decoding tables.
// - dst must be a byte slice of exactly the right size to fit the three lines;
//   the raw read-length is inferred as (len(dst)-4) >> 1.  (Windows
//   line-breaks are not supported.)
// - Length of src is validated
// - This is designed for read-length >= 32.  It still produces the correct
//   result for smaller lengths, but there is a fairly simple faster algorithm
//   (using a pair of 256-element uint16 lookup tables and encoding/binary's
//   binary.LittleEndian.PutUint16() function) for that case, which is being
//   omitted for now due to irrelevance for our current use cases.
func FillFastqRecordBodyFromNibbles(dst, src []byte, baseTablePtr, qualTablePtr *NibbleLookupTable) {
	readLen := (len(dst) - 4) >> 1
	nSrcFullByte := readLen >> 1
	srcOdd := readLen & 1
	if len(src) != nSrcFullByte+srcOdd {
		// 2x is due to each base appearing on both the base and qual lines.
		// +4 is due to three line-feeds and one '+' character.
		panic("FillFastqRecordBodyFromNibbles() requires len(src) = (<# of bases> + 1) / 2, and len(dst) = 2 * <# of bases> + 4.")
	}
	quals := dst[readLen+3 : 2*readLen+3]
	for srcPos := 0; srcPos != nSrcFullByte; srcPos++ {
		srcByte := src[srcPos]
		srcLowBits := srcByte & 15
		dst[2*srcPos] = baseTablePtr.Get(srcLowBits)
		quals[2*srcPos] = qualTablePtr.Get(srcLowBits)
		srcHighBits := srcByte >> 4
		dst[2*srcPos+1] = baseTablePtr.Get(srcHighBits)
		quals[2*srcPos+1] = qualTablePtr.Get(srcHighBits)
	}
	if srcOdd == 1 {
		srcLowBits := src[nSrcFullByte] & 15
		dst[2*nSrcFullByte] = baseTablePtr.Get(srcLowBits)
		quals[2*nSrcFullByte] = qualTablePtr.Get(srcLowBits)
	}
	copy(dst[readLen:readLen+3], "\n+\n")
	dst[2*readLen+3] = '\n'
}
