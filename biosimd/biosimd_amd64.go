// Copyright 2018 GRAIL, Inc.  All rights reserved.
// Use of this source code is governed by the Apache-2.0
// license that can be found in the LICENSE file.

// +build amd64,!appengine

package biosimd

import (
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
var bytesPerVec int

// log2BytesPerVec supports efficient division by bytesPerVec.
var log2BytesPerVec uint

// *** the following functions are defined in biosimd_amd64.s

//go:noescape
func hasSSE42Asm() bool

//go:noescape
func unpackSeqSSE2Asm(dst, src unsafe.Pointer, nSrcByte int)

//go:noescape
func unpackSeqOddSSE2Asm(dst, src unsafe.Pointer, nSrcFullByte int)

//go:noescape
func packSeqSSE41Asm(dst, src unsafe.Pointer, nSrcByte int)

//go:noescape
func packSeqOddSSSE3Asm(dst, src unsafe.Pointer, nDstFullByte int)

//go:noescape
func reverseComp4InplaceTinySSSE3Asm(seq8 unsafe.Pointer, nByte int)

//go:noescape
func reverseComp4InplaceSSSE3Asm(seq8 unsafe.Pointer, nByte int)

//go:noescape
func reverseComp4TinySSSE3Asm(dst, src unsafe.Pointer, nByte int)

//go:noescape
func reverseComp4SSSE3Asm(dst, src unsafe.Pointer, nByte int)

//go:noescape
func reverseComp2InplaceTinySSSE3Asm(acgt8 unsafe.Pointer, nByte int)

//go:noescape
func reverseComp2InplaceSSSE3Asm(acgt8 unsafe.Pointer, nByte int)

//go:noescape
func reverseComp2TinySSSE3Asm(dst, src unsafe.Pointer, nByte int)

//go:noescape
func reverseComp2SSSE3Asm(dst, src unsafe.Pointer, nByte int)

//go:noescape
func unpackAndReplaceSeqSSSE3Asm(dst, src, tablePtr unsafe.Pointer, nSrcByte int)

//go:noescape
func unpackAndReplaceSeqOddSSSE3Asm(dst, src, tablePtr unsafe.Pointer, nSrcFullByte int)

//		unpackAndReplaceSeqOddSSSE3Asm(unsafe.Pointer(dstHeader.Data), unsafe.Pointer(srcHeader.Data), unsafe.Pointer(tablePtr), nSrcFullByte)

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
// WARNING: This is a function designed to be used in inner loops, which makes
// assumptions about length and capacity which aren't checked at runtime.  Use
// the safe version of this function when that's a problem.
// Assumptions #2-3 are always satisfied when the last
// potentially-size-increasing operation on src[] is simd.{Re}makeUnsafe(),
// ResizeUnsafe(), or XcapUnsafe(), and the same is true for dst[].
// 1. len(src) = (len(dst) + 1) / 2.
// 2. Capacity of src is at least RoundUpPow2(len(src) + 1, bytesPerVec), and
//    the same is true for dst.
// 3. The caller does not care if a few bytes past the end of dst[] are
//    changed.
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
	if nSrcFullByte < bytesPerVec {
		for srcPos := 0; srcPos < nSrcFullByte; srcPos++ {
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
// WARNING: This is a function designed to be used in inner loops, which makes
// assumptions about length and capacity which aren't checked at runtime.  Use
// the safe version of this function when that's a problem.
// Assumptions #3-4 are always satisfied when the last
// potentially-size-increasing operation on src[] is simd.{Re}makeUnsafe(),
// ResizeUnsafe(), or XcapUnsafe(), and the same is true for dst[].
// 1. len(dst) = (len(src) + 1) / 2.
// 2. All elements of src[] are less than 16.
// 3. Capacity of src is at least RoundUpPow2(len(src) + 1, bytesPerVec), and
//    the same is true for dst.
// 4. The caller does not care if a few bytes past the end of dst[] are
//    changed.
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
// This is the inverse of UnpackSeq().
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
	if nDstFullByte < bytesPerVec {
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

// ReverseComp4UnsafeInplace reverse-complements seq8[], assuming that it's
// using .bam seq-field encoding with one 4-bit byte per base.
// WARNING: This is a function designed to be used in inner loops, which makes
// assumptions about length and capacity which aren't checked at runtime.  Use
// the safe version of this function when that's a problem.
// Assumptions #2-3 are always satisfied when the last
// potentially-size-increasing operation on seq8[] is simd.{Re}makeUnsafe(),
// ResizeUnsafe(), or XcapUnsafe().
// 1. All elements of seq8[] are less than 16.
// 2. Capacity of seq8 is at least RoundUpPow2(len(seq8) + 1, bytesPerVec).
// 3. The caller does not care if a few bytes past the end of seq8[] are
//    changed.
func ReverseComp4UnsafeInplace(seq8 []byte) {
	nByte := len(seq8)
	seq8Header := (*reflect.SliceHeader)(unsafe.Pointer(&seq8))
	if nByte <= bytesPerVec {
		reverseComp4InplaceTinySSSE3Asm(unsafe.Pointer(seq8Header.Data), nByte)
		return
	}
	reverseComp4InplaceSSSE3Asm(unsafe.Pointer(seq8Header.Data), nByte)
}

var revComp4Table = [...]byte{0, 8, 4, 12, 2, 10, 6, 14, 1, 9, 5, 13, 3, 11, 7, 15}

// ReverseComp4Inplace reverse-complements seq8[], assuming that it's using
// .bam seq-field encoding with one 4-bit byte per base.
// WARNING: If a seq8[] value is larger than 15, it's possible for this to
// immediately crash, and it's also possible for this to return and fill seq8[]
// with garbage.  Only promise is that we don't scribble over arbitrary memory.
func ReverseComp4Inplace(seq8 []byte) {
	// This has no performance penalty relative to the unsafe function when nByte
	// > 16, but the penalty can be horrific (in relative terms) below that.
	nByte := len(seq8)
	if nByte <= bytesPerVec {
		nByteDiv2 := nByte >> 1
		for idx, invIdx := 0, nByte-1; idx != nByteDiv2; idx, invIdx = idx+1, invIdx-1 {
			seq8[idx], seq8[invIdx] = revComp4Table[seq8[invIdx]], revComp4Table[seq8[idx]]
		}
		if nByte&1 == 1 {
			seq8[nByteDiv2] = revComp4Table[seq8[nByteDiv2]]
		}
		return
	}
	seq8Header := (*reflect.SliceHeader)(unsafe.Pointer(&seq8))
	reverseComp4InplaceSSSE3Asm(unsafe.Pointer(seq8Header.Data), nByte)
}

// ReverseComp4Unsafe saves the reverse-complement of src[] to dst[], assuming
// .bam seq-field encoding with one 4-bit byte per base.
// WARNING: This is a function designed to be used in inner loops, which makes
// assumptions about length and capacity which aren't checked at runtime.  Use
// the safe version of this function when that's a problem.
// Assumptions #3-4 are always satisfied when the last
// potentially-size-increasing operation on src[] is simd.{Re}makeUnsafe(),
// ResizeUnsafe(), or XcapUnsafe(), and the same is true of dst[].
// 1. len(src) == len(dst).
// 2. All elements of src[] are less than 16.
// 3. Capacity of src is at least RoundUpPow2(len(src) + 1, bytesPerVec), and
//    the same is true of dst.
// 4. The caller does not care if a few bytes past the end of dst[] are
//    changed.
func ReverseComp4Unsafe(dst, src []byte) {
	nByte := len(src)
	srcHeader := (*reflect.SliceHeader)(unsafe.Pointer(&src))
	dstHeader := (*reflect.SliceHeader)(unsafe.Pointer(&dst))
	if nByte <= bytesPerVec {
		reverseComp4TinySSSE3Asm(unsafe.Pointer(dstHeader.Data), unsafe.Pointer(srcHeader.Data), nByte)
		return
	}
	reverseComp4SSSE3Asm(unsafe.Pointer(dstHeader.Data), unsafe.Pointer(srcHeader.Data), nByte)
}

// ReverseComp4 saves the reverse-complement of src[] to dst[], assuming .bam
// seq-field encoding with one 4-bit byte per base.
// It panics if len(dst) != len(src).
// WARNING: If a src[] value is larger than 15, it's possible for this to
// immediately crash, and it's also possible for this to return and fill src[]
// with garbage.  Only promise is that we don't scribble over arbitrary memory.
func ReverseComp4(dst, src []byte) {
	nByte := len(src)
	if len(dst) != len(src) {
		panic("ReverseComp4() requires len(dst) == len(src).")
	}
	if nByte < bytesPerVec {
		for idx, invIdx := 0, nByte-1; idx != nByte; idx, invIdx = idx+1, invIdx-1 {
			dst[idx] = revComp4Table[src[invIdx]]
		}
		return
	}
	srcHeader := (*reflect.SliceHeader)(unsafe.Pointer(&src))
	dstHeader := (*reflect.SliceHeader)(unsafe.Pointer(&dst))
	reverseComp4SSSE3Asm(unsafe.Pointer(dstHeader.Data), unsafe.Pointer(srcHeader.Data), nByte)
}

// ReverseComp2UnsafeInplace reverse-complements acgt8[], assuming that it's
// encoded with one byte per base, ACGT=0123.
// WARNING: This is a function designed to be used in inner loops, which makes
// assumptions about length and capacity which aren't checked at runtime.  Use
// the safe version of this function when that's a problem.
// These assumptions are always satisfied when the last
// potentially-size-increasing operation on acgt8[] is simd.{Re}makeUnsafe(),
// ResizeUnsafe(), or XcapUnsafe().
// 1. Capacity of acgt8[] is at least RoundUpPow2(len(acgt8) + 1, bytesPerVec).
// 2. The caller does not care if a few bytes past the end of acgt8[] are
//    changed.
func ReverseComp2UnsafeInplace(acgt8 []byte) {
	nByte := len(acgt8)
	acgt8Header := (*reflect.SliceHeader)(unsafe.Pointer(&acgt8))
	if nByte <= bytesPerVec {
		reverseComp2InplaceTinySSSE3Asm(unsafe.Pointer(acgt8Header.Data), nByte)
		return
	}
	reverseComp2InplaceSSSE3Asm(unsafe.Pointer(acgt8Header.Data), nByte)
}

// ReverseComp2Inplace reverse-complements acgt8[], assuming that it's encoded
// with one byte per base, ACGT=0123.
func ReverseComp2Inplace(acgt8 []byte) {
	// This has no performance penalty relative to the unsafe function when nByte
	// > 16, but the penalty can be horrific (in relative terms) below that.
	nByte := len(acgt8)
	if nByte <= bytesPerVec {
		simd.Reverse8Inplace(acgt8)
		simd.XorConst8Inplace(acgt8, 3)
		return
	}
	acgt8Header := (*reflect.SliceHeader)(unsafe.Pointer(&acgt8))
	reverseComp2InplaceSSSE3Asm(unsafe.Pointer(acgt8Header.Data), nByte)
}

// ReverseComp2Unsafe saves the reverse-complement of src[] to dst[], assuming
// that they're encoded with one byte per base, ACGT=0123.
// WARNING: This is a function designed to be used in inner loops, which makes
// assumptions about length and capacity which aren't checked at runtime.  Use
// the safe version of this function when that's a problem.
// Assumptions #2-3 are always satisfied when the last
// potentially-size-increasing operation on src[] is simd.{Re}makeUnsafe(),
// ResizeUnsafe(), or XcapUnsafe(), and the same is true of dst[].
// 1. len(src) == len(dst).
// 2. Capacity of src is at least RoundUpPow2(len(src) + 1, bytesPerVec), and
//    the same is true of dst.
// 3. The caller does not care if a few bytes past the end of dst[] are
//    changed.
func ReverseComp2Unsafe(dst, src []byte) {
	nByte := len(src)
	srcHeader := (*reflect.SliceHeader)(unsafe.Pointer(&src))
	dstHeader := (*reflect.SliceHeader)(unsafe.Pointer(&dst))
	if nByte <= bytesPerVec {
		reverseComp2TinySSSE3Asm(unsafe.Pointer(dstHeader.Data), unsafe.Pointer(srcHeader.Data), nByte)
		return
	}
	reverseComp2SSSE3Asm(unsafe.Pointer(dstHeader.Data), unsafe.Pointer(srcHeader.Data), nByte)
}

// ReverseComp2 saves the reverse-complement of src[] to dst[], assuming that
// they're encoded with one byte per base, ACGT=0123.
// It panics if len(dst) != len(src).
func ReverseComp2(dst, src []byte) {
	nByte := len(src)
	if len(dst) != len(src) {
		panic("ReverseComp2() requires len(dst) == len(src).")
	}
	if nByte < bytesPerVec {
		simd.Reverse8(dst, src)
		simd.XorConst8Inplace(dst, 3)
		return
	}
	srcHeader := (*reflect.SliceHeader)(unsafe.Pointer(&src))
	dstHeader := (*reflect.SliceHeader)(unsafe.Pointer(&dst))
	reverseComp2SSSE3Asm(unsafe.Pointer(dstHeader.Data), unsafe.Pointer(srcHeader.Data), nByte)
}

// UnpackAndReplaceSeqUnsafe sets the bytes in dst[] as follows:
//   if pos is even, dst[pos] := table[src[pos / 2] >> 4]
//   if pos is odd, dst[pos] := table[src[pos / 2] & 15]
// It panics if len(src) != (len(dst) + 1) / 2.
// WARNING: This is a function designed to be used in inner loops, which makes
// assumptions about length and capacity which aren't checked at runtime.  Use
// the safe version of this function when that's a problem.
// Assumptions #2-#3 are always satisfied when the last
// potentially-size-increasing operation on src[] is {Re}makeUnsafe(),
// ResizeUnsafe(), or XcapUnsafe(), and the same is true for dst[].
// 1. len(src) == (len(dst) + 1) / 2.
// 2. Capacity of src is at least RoundUpPow2(len(src), bytesPerVec), and the
//    same is true for dst.
// 3. The caller does not care if a few bytes past the end of dst[] are
//    changed.
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
	if nSrcFullByte < bytesPerVec {
		for srcPos := 0; srcPos < nSrcFullByte; srcPos++ {
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
