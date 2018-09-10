// Copyright 2018 GRAIL, Inc.  All rights reserved.
// Use of this source code is governed by the Apache-2.0
// license that can be found in the LICENSE file.

// +build amd64,!appengine

package biosimd

import (
	"fmt"
	"reflect"
	"unsafe"

	"github.com/grailbio/base/simd"
)

// amd64 compile-time constants.  Private base/simd constants are recalculated
// here; probably want to change that.

// BytesPerWord is the number of bytes in a machine word.
const BytesPerWord = simd.BytesPerWord

// Log2BytesPerWord is log2(BytesPerWord).  This is relevant for manual
// bit-shifting when we know that's a safe way to divide and the compiler does
// not (e.g. dividend is of signed int type).
const Log2BytesPerWord = simd.Log2BytesPerWord

// const minPageSize = 4096  may be relevant for safe functions soon.

// These could be compile-time constants for now, but not after AVX2
// autodetection is added.

// bytesPerVec is the size of the maximum-width vector that may be used.  It is
// currently always 16, but it will be set to larger values at runtime in the
// future when AVX2/AVX-512/etc. is detected.
// (Probably use exported version of this from base/simd in the future.)
var bytesPerVec int

// log2BytesPerVec supports efficient division by bytesPerVec.
var log2BytesPerVec uint

//go:linkname hasSSE42Asm github.com/grailbio/base/simd.hasSSE42Asm
func hasSSE42Asm() bool

// *** the following functions are defined in biosimd_amd64.s

//go:noescape
func unpackSeqSSE2Asm(dst, src unsafe.Pointer, nSrcByte int)

//go:noescape
func unpackSeqOddSSE2Asm(dst, src unsafe.Pointer, nSrcFullByte int)

//go:noescape
func packSeqSSE41Asm(dst, src unsafe.Pointer, nSrcByte int)

//go:noescape
func packSeqOddSSSE3Asm(dst, src unsafe.Pointer, nDstFullByte int)

//go:noescape
func unpackAndReplaceSeqSSSE3Asm(dst, src, tablePtr unsafe.Pointer, nSrcByte int)

//go:noescape
func unpackAndReplaceSeqOddSSSE3Asm(dst, src, tablePtr unsafe.Pointer, nSrcFullByte int)

//go:noescape
func cleanASCIISeqInplaceSSSE3Asm(ascii8 unsafe.Pointer, nByte int)

//go:noescape
func cleanASCIISeqNoCapitalizeInplaceSSSE3Asm(ascii8 unsafe.Pointer, nByte int)

//go:noescape
func isNonACGTPresentSSE41Asm(ascii8 unsafe.Pointer, nonTTablePtr *[16]byte, nByte int) bool

//go:noescape
func asciiToSeq8SSSE3Asm(dst, src unsafe.Pointer, nByte int)

// *** end assembly function signatures

func init() {
	if !hasSSE42Asm() {
		panic("SSE4.2 required.")
	}
	bytesPerVec = 16
	log2BytesPerVec = 4
}

// UnpackSeqUnsafe sets the bytes in dst[] as follows:
//   if pos is even, dst[pos] := src[pos / 2] >> 4
//   if pos is odd, dst[pos] := src[pos / 2] & 15
//
// WARNING: This is a function designed to be used in inner loops, which makes
// assumptions about length and capacity which aren't checked at runtime.  Use
// the safe version of this function when that's a problem.
// Assumptions #2-3 are always satisfied when the last
// potentially-size-increasing operation on src[] is simd.{Re}makeUnsafe(),
// ResizeUnsafe(), or XcapUnsafe(), and the same is true for dst[].
//
// 1. len(src) = (len(dst) + 1) / 2.
//
// 2. Capacity of src is at least RoundUpPow2(len(src) + 1, bytesPerVec), and
// the same is true for dst.
//
// 3. The caller does not care if a few bytes past the end of dst[] are
// changed.
func UnpackSeqUnsafe(dst, src []byte) {
	// Based on simd.PackedNibbleLookupUnsafe().  Differences are (i) even/odd is
	// swapped, and (ii) no table lookup is necessary.
	srcHeader := (*reflect.SliceHeader)(unsafe.Pointer(&src))
	dstHeader := (*reflect.SliceHeader)(unsafe.Pointer(&dst))
	unpackSeqSSE2Asm(unsafe.Pointer(dstHeader.Data), unsafe.Pointer(srcHeader.Data), srcHeader.Len)
}

// UnpackSeq sets the bytes in dst[] as follows:
//   if pos is even, dst[pos] := src[pos / 2] >> 4
//   if pos is odd, dst[pos] := src[pos / 2] & 15
// It panics if len(src) != (len(dst) + 1) / 2.
//
// Nothing bad happens if len(dst) is odd and some low bits in the last src[]
// byte are set, though it's generally good practice to ensure that case
// doesn't come up.
func UnpackSeq(dst, src []byte) {
	// This takes ~20-35% longer than the unsafe function on the short-array
	// benchmark.
	dstLen := len(dst)
	nSrcFullByte := dstLen >> 1
	srcOdd := dstLen & 1
	if len(src) != nSrcFullByte+srcOdd {
		panic("UnpackSeq() requires len(src) == (len(dst) + 1) / 2.")
	}
	if nSrcFullByte < 16 {
		for srcPos := 0; srcPos != nSrcFullByte; srcPos++ {
			srcByte := src[srcPos]
			dst[2*srcPos] = srcByte >> 4
			dst[2*srcPos+1] = srcByte & 15
		}
	} else {
		srcHeader := (*reflect.SliceHeader)(unsafe.Pointer(&src))
		dstHeader := (*reflect.SliceHeader)(unsafe.Pointer(&dst))
		unpackSeqOddSSE2Asm(unsafe.Pointer(dstHeader.Data), unsafe.Pointer(srcHeader.Data), nSrcFullByte)
	}
	if srcOdd == 1 {
		srcByte := src[nSrcFullByte]
		dst[2*nSrcFullByte] = srcByte >> 4
	}
}

// PackSeqUnsafe sets the bytes in dst[] as follows:
//   if pos is even, high 4 bits of dst[pos / 2] := src[pos]
//   if pos is odd, low 4 bits of dst[pos / 2] := src[pos]
//   if len(src) is odd, the low 4 bits of dst[len(src) / 2] are zero
// This is the inverse of UnpackSeqUnsafe().
//
// WARNING: This is a function designed to be used in inner loops, which makes
// assumptions about length and capacity which aren't checked at runtime.  Use
// the safe version of this function when that's a problem.
// Assumptions #3-4 are always satisfied when the last
// potentially-size-increasing operation on src[] is simd.{Re}makeUnsafe(),
// ResizeUnsafe(), or XcapUnsafe(), and the same is true for dst[].
//
// 1. len(dst) = (len(src) + 1) / 2.
//
// 2. All elements of src[] are less than 16.
//
// 3. Capacity of src is at least RoundUpPow2(len(src) + 1, bytesPerVec), and
// the same is true for dst.
//
// 4. The caller does not care if a few bytes past the end of dst[] are
// changed.
func PackSeqUnsafe(dst, src []byte) {
	srcLen := len(src)
	srcHeader := (*reflect.SliceHeader)(unsafe.Pointer(&src))
	dstHeader := (*reflect.SliceHeader)(unsafe.Pointer(&dst))
	packSeqSSE41Asm(unsafe.Pointer(dstHeader.Data), unsafe.Pointer(srcHeader.Data), srcLen)
	if srcLen&1 == 1 {
		// Force low bits of last dst[] byte to zero.
		dst[srcLen>>1] = src[srcLen-1] << 4
	}
}

// PackSeq sets the bytes in dst[] as follows:
//   if pos is even, high 4 bits of dst[pos / 2] := src[pos]
//   if pos is odd, low 4 bits of dst[pos / 2] := src[pos]
//   if len(src) is odd, the low 4 bits of dst[len(src) / 2] are zero
// It panics if len(dst) != (len(src) + 1) / 2.
//
// This is the inverse of UnpackSeq().
//
// WARNING: Actual values in dst[] bytes may be garbage if any src[] bytes are
// greater than 15; this function only guarantees that no buffer overflow will
// occur.
func PackSeq(dst, src []byte) {
	// This takes ~4-7% longer than the unsafe function on the short-array
	// benchmark.
	srcLen := len(src)
	nDstFullByte := srcLen >> 1
	dstOdd := srcLen & 1
	if len(dst) != nDstFullByte+dstOdd {
		panic("PackSeq() requires len(dst) == (len(src) + 1) / 2.")
	}
	if nDstFullByte < 16 {
		for dstPos := 0; dstPos < nDstFullByte; dstPos++ {
			dst[dstPos] = (src[2*dstPos] << 4) | src[2*dstPos+1]
		}
	} else {
		srcHeader := (*reflect.SliceHeader)(unsafe.Pointer(&src))
		dstHeader := (*reflect.SliceHeader)(unsafe.Pointer(&dst))
		packSeqOddSSSE3Asm(unsafe.Pointer(dstHeader.Data), unsafe.Pointer(srcHeader.Data), nDstFullByte)
	}
	if dstOdd == 1 {
		dst[nDstFullByte] = src[nDstFullByte*2] << 4
	}
}

// UnpackAndReplaceSeqUnsafe sets the bytes in dst[] as follows:
//   if pos is even, dst[pos] := table[src[pos / 2] >> 4]
//   if pos is odd, dst[pos] := table[src[pos / 2] & 15]
// It panics if len(src) != (len(dst) + 1) / 2.
//
// WARNING: This is a function designed to be used in inner loops, which makes
// assumptions about length and capacity which aren't checked at runtime.  Use
// the safe version of this function when that's a problem.
// Assumptions #2-#3 are always satisfied when the last
// potentially-size-increasing operation on src[] is {Re}makeUnsafe(),
// ResizeUnsafe(), or XcapUnsafe(), and the same is true for dst[].
//
// 1. len(src) == (len(dst) + 1) / 2.
//
// 2. Capacity of src is at least RoundUpPow2(len(src) + 1, bytesPerVec), and
// the same is true for dst.
//
// 3. The caller does not care if a few bytes past the end of dst[] are
// changed.
func UnpackAndReplaceSeqUnsafe(dst, src []byte, tablePtr *[16]byte) {
	// Minor variant of simd.PackedNibbleLookupUnsafe().
	srcHeader := (*reflect.SliceHeader)(unsafe.Pointer(&src))
	dstHeader := (*reflect.SliceHeader)(unsafe.Pointer(&dst))
	unpackAndReplaceSeqSSSE3Asm(unsafe.Pointer(dstHeader.Data), unsafe.Pointer(srcHeader.Data), unsafe.Pointer(tablePtr), srcHeader.Len)
}

// UnpackAndReplaceSeq sets the bytes in dst[] as follows:
//   if pos is even, dst[pos] := table[src[pos / 2] >> 4]
//   if pos is odd, dst[pos] := table[src[pos / 2] & 15]
// It panics if len(src) != (len(dst) + 1) / 2.
//
// Nothing bad happens if len(dst) is odd and some low bits in the last src[]
// byte are set, though it's generally good practice to ensure that case
// doesn't come up.
func UnpackAndReplaceSeq(dst, src []byte, tablePtr *[16]byte) {
	// Minor variant of simd.PackedNibbleLookup().
	dstLen := len(dst)
	nSrcFullByte := dstLen >> 1
	srcOdd := dstLen & 1
	if len(src) != nSrcFullByte+srcOdd {
		panic("UnpackAndReplaceSeq() requires len(src) == (len(dst) + 1) / 2.")
	}
	if nSrcFullByte < 16 {
		for srcPos := 0; srcPos != nSrcFullByte; srcPos++ {
			srcByte := src[srcPos]
			dst[2*srcPos] = tablePtr[srcByte>>4]
			dst[2*srcPos+1] = tablePtr[srcByte&15]
		}
	} else {
		srcHeader := (*reflect.SliceHeader)(unsafe.Pointer(&src))
		dstHeader := (*reflect.SliceHeader)(unsafe.Pointer(&dst))
		unpackAndReplaceSeqOddSSSE3Asm(unsafe.Pointer(dstHeader.Data), unsafe.Pointer(srcHeader.Data), unsafe.Pointer(tablePtr), nSrcFullByte)
	}
	if srcOdd == 1 {
		srcByte := src[nSrcFullByte]
		dst[2*nSrcFullByte] = tablePtr[srcByte>>4]
	}
}

// UnpackAndReplaceSeqSubset sets the bytes in dst[] as follows:
//   if srcPos is even, dst[srcPos-startPos] := table[src[srcPos / 2] >> 4]
//   if srcPos is odd, dst[srcPos-startPos] := table[src[srcPos / 2] & 15]
// It panics if len(dst) != endPos - startPos, startPos < 0, or
// len(src) * 2 < endPos.
func UnpackAndReplaceSeqSubset(dst, src []byte, tablePtr *[16]byte, startPos, endPos int) {
	if (startPos < 0) || (len(src)*2 < endPos) {
		errstr := fmt.Sprintf("UnpackAndReplaceSeqSubset() requires 0 <= startPos <= endPos <= 2 * len(src).\n  len(dst) = %d\n  src = %v\n  startPos = %d\n  endPos = %d\n", len(dst), src, startPos, endPos)
		panic(errstr)
	}
	dstLen := len(dst)
	if dstLen != endPos-startPos {
		errstr := fmt.Sprintf("UnpackAndReplaceSeqSubset() requires len(dst) == endPos - startPos.\n  len(dst) = %d\n  startPos = %d\n  endPos = %d\n", dstLen, startPos, endPos)
		panic(errstr)
	}
	if dstLen == 0 {
		return
	}
	startOffset := startPos >> 1
	startPosOdd := startPos & 1
	if startPosOdd == 1 {
		dst[0] = tablePtr[src[startOffset]&15]
		startOffset++
	}
	nSrcFullByte := (dstLen - startPosOdd) >> 1
	if nSrcFullByte < 16 {
		for srcPos := 0; srcPos != nSrcFullByte; srcPos++ {
			srcByte := src[srcPos+startOffset]
			dst[2*srcPos+startPosOdd] = tablePtr[srcByte>>4]
			dst[2*srcPos+1+startPosOdd] = tablePtr[srcByte&15]
		}
	} else {
		srcHeader := (*reflect.SliceHeader)(unsafe.Pointer(&src))
		dstHeader := (*reflect.SliceHeader)(unsafe.Pointer(&dst))
		unpackAndReplaceSeqOddSSSE3Asm(unsafe.Pointer(dstHeader.Data+uintptr(startPosOdd)), unsafe.Pointer(srcHeader.Data+uintptr(startOffset)), unsafe.Pointer(tablePtr), nSrcFullByte)
	}
	if endPos&1 == 1 {
		srcByte := src[nSrcFullByte+startOffset]
		dst[dstLen-1] = tablePtr[srcByte>>4]
	}
}

var cleanASCIISeqTable = [...]byte{
	'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
	'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
	'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
	'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
	'N', 'A', 'N', 'C', 'N', 'N', 'N', 'G', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
	'N', 'N', 'N', 'N', 'T', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
	'N', 'A', 'N', 'C', 'N', 'N', 'N', 'G', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
	'N', 'N', 'N', 'N', 'T', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
	'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
	'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
	'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
	'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
	'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
	'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
	'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
	'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N'}

// CleanASCIISeqInplace capitalizes 'a'/'c'/'g'/'t', and replaces everything
// non-ACGT with 'N'.
func CleanASCIISeqInplace(ascii8 []byte) {
	nByte := len(ascii8)
	if nByte < 16 {
		for pos, ascii8Byte := range ascii8 {
			ascii8[pos] = cleanASCIISeqTable[ascii8Byte]
		}
		return
	}
	ascii8Header := (*reflect.SliceHeader)(unsafe.Pointer(&ascii8))
	cleanASCIISeqInplaceSSSE3Asm(unsafe.Pointer(ascii8Header.Data), nByte)
}

var cleanASCIISeqNoCapitalizeTable = [...]byte{
	'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
	'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
	'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
	'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
	'N', 'A', 'N', 'C', 'N', 'N', 'N', 'G', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
	'N', 'N', 'N', 'N', 'T', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
	'N', 'a', 'N', 'c', 'N', 'N', 'N', 'g', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
	'N', 'N', 'N', 'N', 't', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
	'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
	'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
	'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
	'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
	'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
	'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
	'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
	'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N'}

// CleanASCIISeqNoCapitalizeInplace replaces everything non-ACGTacgt with 'N'.
func CleanASCIISeqNoCapitalizeInplace(ascii8 []byte) {
	nByte := len(ascii8)
	if nByte < 16 {
		for pos, ascii8Byte := range ascii8 {
			ascii8[pos] = cleanASCIISeqNoCapitalizeTable[ascii8Byte]
		}
		return
	}
	ascii8Header := (*reflect.SliceHeader)(unsafe.Pointer(&ascii8))
	cleanASCIISeqNoCapitalizeInplaceSSSE3Asm(unsafe.Pointer(ascii8Header.Data), nByte)
}

var isNotCapitalACGTTable = [...]bool{
	true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true,
	true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true,
	true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true,
	true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true,
	true, false, true, false, true, true, true, false, true, true, true, true, true, true, true, true,
	true, true, true, true, false, true, true, true, true, true, true, true, true, true, true, true,
	true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true,
	true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true,
	true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true,
	true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true,
	true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true,
	true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true,
	true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true,
	true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true,
	true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true,
	true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true}

var acgTable16 = [...]byte{
	0xff, 0, 0xff, 0, 0xff, 0xff, 0xff, 0, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff}

// IsNonACGTPresent returns true iff there is a non-capital-ACGT character in
// the slice.
func IsNonACGTPresent(ascii8 []byte) bool {
	nByte := len(ascii8)
	if nByte < 16 {
		for _, ascii8Byte := range ascii8 {
			if isNotCapitalACGTTable[ascii8Byte] {
				return true
			}
		}
		return false
	}
	ascii8Header := (*reflect.SliceHeader)(unsafe.Pointer(&ascii8))
	return isNonACGTPresentSSE41Asm(unsafe.Pointer(ascii8Header.Data), &acgTable16, nByte)
}

var isNotCapitalACGTNTable = [...]bool{
	true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true,
	true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true,
	true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true,
	true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true,
	true, false, true, false, true, true, true, false, true, true, true, true, true, true, false, true,
	true, true, true, true, false, true, true, true, true, true, true, true, true, true, true, true,
	true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true,
	true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true,
	true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true,
	true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true,
	true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true,
	true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true,
	true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true,
	true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true,
	true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true,
	true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true}

var acgnTable16 = [...]byte{
	0xff, 0, 0xff, 0, 0xff, 0xff, 0xff, 0, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0, 0xff}

// IsNonACGTNPresent returns true iff there is a non-capital-ACGTN character in
// the slice.
func IsNonACGTNPresent(ascii8 []byte) bool {
	nByte := len(ascii8)
	if nByte < 16 {
		for _, ascii8Byte := range ascii8 {
			if isNotCapitalACGTNTable[ascii8Byte] {
				return true
			}
		}
		return false
	}
	ascii8Header := (*reflect.SliceHeader)(unsafe.Pointer(&ascii8))
	return isNonACGTPresentSSE41Asm(unsafe.Pointer(ascii8Header.Data), &acgnTable16, nByte)
}

var asciiToSeq8Table = [...]byte{
	15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
	15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
	15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
	15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
	15, 1, 15, 2, 15, 15, 15, 4, 15, 15, 15, 15, 15, 15, 15, 15,
	15, 15, 15, 15, 8, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
	15, 1, 15, 2, 15, 15, 15, 4, 15, 15, 15, 15, 15, 15, 15, 15,
	15, 15, 15, 15, 8, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
	15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
	15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
	15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
	15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
	15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
	15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
	15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
	15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15}

// ASCIIToSeq8 sets dst[pos] as follows:
//   src[pos] == 'A'/'a': dst[pos] == 1
//   src[pos] == 'C'/'c': dst[pos] == 2
//   src[pos] == 'G'/'g': dst[pos] == 4
//   src[pos] == 'T'/'t': dst[pos] == 8
//   src[pos] == anything else: dst[pos] == 15
// It panics if len(dst) != len(src).
func ASCIIToSeq8(dst, src []byte) {
	// This is good for unvalidated .fa loading when you're fine with treating
	// all non-ACGT characters as N.
	nByte := len(src)
	if len(dst) != nByte {
		panic("ASCIIToSeq8() requires len(src) == len(dst).")
	}
	if nByte < 16 {
		for pos, srcByte := range src {
			dst[pos] = asciiToSeq8Table[srcByte]
		}
		return
	}
	srcHeader := (*reflect.SliceHeader)(unsafe.Pointer(&src))
	dstHeader := (*reflect.SliceHeader)(unsafe.Pointer(&dst))
	asciiToSeq8SSSE3Asm(unsafe.Pointer(dstHeader.Data), unsafe.Pointer(srcHeader.Data), nByte)
}
