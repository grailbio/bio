// Copyright 2019 GRAIL, Inc.  All rights reserved.
// Use of this source code is governed by the Apache-2.0
// license that can be found in the LICENSE file.

// +build amd64,!appengine

package biosimd

import (
	"unsafe"
)

//go:noescape
func fillFastqRecordBodyFromNibblesSSSE3Asm(dst, src, baseTablePtr, qualTablePtr unsafe.Pointer, nBase int)

// FillFastqRecordBodyFromNibbles fills the body (defined as the last three
// lines) of a 4-line FASTQ record, given a packed 4-bit representation of the
// base+qual information and the decoding tables.  (Windows line-breaks are not
// supported.)
// - len(dst) must be at least 2 * nBase + 4, but it's allowed to be larger.
// - len(src) must be at least (nBase + 1) >> 1, but it's allowed to be larger.
// - This is designed for read-length >= 32.  It still produces the correct
//   result for smaller lengths, but there is a fairly simple faster algorithm
//   (using a pair of 256-element uint16 lookup tables and encoding/binary's
//   binary.LittleEndian.PutUint16() function) for that case, which is being
//   omitted for now due to irrelevance for our current use cases.
func FillFastqRecordBodyFromNibbles(dst, src []byte, nBase int, baseTablePtr, qualTablePtr *NibbleLookupTable) {
	if len(dst) < 2*nBase+4 {
		// 2x is due to each base appearing on both the base and qual lines.
		// +4 is due to three line-feeds and one '+' character.
		panic("FillFastqRecordBodyFromNibbles() requires len(dst) >= 2 * nBase + 4.")
	}
	nSrcFullByte := nBase >> 1
	srcOdd := nBase & 1
	if len(src) < nSrcFullByte+srcOdd {
		panic("FillFastqRecordBodyFromNibbles() requires len(src) >= (nBase + 1) / 2.")
	}
	if nBase >= 32 {
		// Note that unsafe.Pointer(&dst[0]) breaks when len(dst) == 0.
		fillFastqRecordBodyFromNibblesSSSE3Asm(unsafe.Pointer(&dst[0]), unsafe.Pointer(&src[0]), unsafe.Pointer(baseTablePtr), unsafe.Pointer(qualTablePtr), nBase)
		return
	}
	quals := dst[nBase+3 : 2*nBase+3]
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
	copy(dst[nBase:nBase+3], "\n+\n")
	dst[2*nBase+3] = '\n'
}
