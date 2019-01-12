// Copyright 2018 GRAIL, Inc.  All rights reserved.
// Use of this source code is governed by the Apache-2.0
// license that can be found in the LICENSE file.

// +build !amd64 appengine

package biosimd

import (
	"github.com/grailbio/base/simd"
)

var revComp8Table = [...]byte{
	'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
	'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
	'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
	'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
	'N', 'T', 'N', 'G', 'N', 'N', 'N', 'C', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
	'N', 'N', 'N', 'N', 'A', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
	'N', 'T', 'N', 'G', 'N', 'N', 'N', 'C', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
	'N', 'N', 'N', 'N', 'A', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
	'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
	'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
	'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
	'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
	'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
	'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
	'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
	'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N'}

var revComp8Table16 = [16]byte{
	'N', 'T', 'N', 'G', 'A', 'N', 'N', 'C', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N'}

// ReverseComp8InplaceNoValidate reverse-complements ascii8[], assuming that
// it's using ASCII encoding, and all values are in {0, '0', 'A', 'C', 'G',
// 'T', 'N', 'a', 'c', 'g', 't', 'n'}.
//
// If the input assumption is satisfied, output is restricted to
// 'A'/'C'/'G'/'T'/'N'.  Other bytes may be written if the input assumption is
// not satisfied.
//
// This usually takes ~35% less time than the validating function.
func ReverseComp8InplaceNoValidate(ascii8 []byte) {
	nByte := len(ascii8)
	nByteDiv2 := nByte >> 1
	for idx, invIdx := 0, nByte-1; idx != nByteDiv2; idx, invIdx = idx+1, invIdx-1 {
		ascii8[idx], ascii8[invIdx] = revComp8Table[ascii8[invIdx]], revComp8Table[ascii8[idx]]
	}
	if nByte&1 == 1 {
		ascii8[nByteDiv2] = revComp8Table[ascii8[nByteDiv2]]
	}
}

// ReverseComp8Inplace reverse-complements ascii8[], assuming that it's using
// ASCII encoding.  More precisely, it maps 'A'/'a' to 'T', 'C'/'c' to 'G',
// 'G'/'g' to 'C', 'T'/'t' to 'A', and everything else to 'N'.
func ReverseComp8Inplace(ascii8 []byte) {
	nByte := len(ascii8)
	nByteDiv2 := nByte >> 1
	for idx, invIdx := 0, nByte-1; idx != nByteDiv2; idx, invIdx = idx+1, invIdx-1 {
		ascii8[idx], ascii8[invIdx] = revComp8Table[ascii8[invIdx]], revComp8Table[ascii8[idx]]
	}
	if nByte&1 == 1 {
		ascii8[nByteDiv2] = revComp8Table[ascii8[nByteDiv2]]
	}
}

// ReverseComp8NoValidate writes the reverse-complement of src[] to dst[],
// assuming src is using ASCII encoding, and all values are in {0, '0', 'A',
// 'C', 'G', 'T', 'N', 'a', 'c', 'g', 't', 'n'}.
//
// If the input assumption is satisfied, output is restricted to
// 'A'/'C'/'G'/'T'/'N'.  Other bytes may be written if the input assumption is
// not satisfied.
//
// It panics if len(dst) != len(src).
func ReverseComp8NoValidate(dst, src []byte) {
	nByte := len(src)
	if len(dst) != nByte {
		panic("ReverseComp8NoValidate requires len(dst) == len(src).")
	}
	for idx, invIdx := 0, nByte-1; idx != nByte; idx, invIdx = idx+1, invIdx-1 {
		dst[idx] = revComp8Table[src[invIdx]]
	}
}

var revComp4Table = [...]byte{0, 8, 4, 12, 2, 10, 6, 14, 1, 9, 5, 13, 3, 11, 7, 15}

// ReverseComp4UnsafeInplace reverse-complements seq8[], assuming that it's
// using .bam seq-field encoding with one 4-bit byte per base.
//
// WARNING: This is a function designed to be used in inner loops, which makes
// assumptions about length and capacity which aren't checked at runtime.  Use
// the safe version of this function when that's a problem.
// Assumptions #2-3 are always satisfied when the last
// potentially-size-increasing operation on seq8[] is simd.{Re}makeUnsafe(),
// ResizeUnsafe(), or XcapUnsafe().
//
// 1. All elements of seq8[] are less than 16.
//
// 2. Capacity of seq8 is at least RoundUpPow2(len(seq8) + 1, bytesPerVec).
//
// 3. The caller does not care if a few bytes past the end of seq8[] are
// changed.
func ReverseComp4UnsafeInplace(seq8 []byte) {
	nByte := len(seq8)
	nByteDiv2 := nByte >> 1
	for idx, invIdx := 0, nByte-1; idx != nByteDiv2; idx, invIdx = idx+1, invIdx-1 {
		seq8[idx], seq8[invIdx] = revComp4Table[seq8[invIdx]], revComp4Table[seq8[idx]]
	}
	if nByte&1 == 1 {
		seq8[nByteDiv2] = revComp4Table[seq8[nByteDiv2]]
	}
}

// ReverseComp4Inplace reverse-complements seq8[], assuming that it's using
// .bam seq-field encoding with one 4-bit byte per base.
//
// WARNING: If a seq8[] value is larger than 15, it's possible for this to
// immediately crash, and it's also possible for this to return and fill seq8[]
// with garbage.  Only promise is that we don't scribble over arbitrary memory.
func ReverseComp4Inplace(seq8 []byte) {
	nByte := len(seq8)
	nByteDiv2 := nByte >> 1
	for idx, invIdx := 0, nByte-1; idx != nByteDiv2; idx, invIdx = idx+1, invIdx-1 {
		seq8[idx], seq8[invIdx] = revComp4Table[seq8[invIdx]], revComp4Table[seq8[idx]]
	}
	if nByte&1 == 1 {
		seq8[nByteDiv2] = revComp4Table[seq8[nByteDiv2]]
	}
}

// ReverseComp4Unsafe saves the reverse-complement of src[] to dst[], assuming
// .bam seq-field encoding with one 4-bit byte per base.
//
// WARNING: This is a function designed to be used in inner loops, which makes
// assumptions about length and capacity which aren't checked at runtime.  Use
// the safe version of this function when that's a problem.
// Assumptions #3-4 are always satisfied when the last
// potentially-size-increasing operation on src[] is simd.{Re}makeUnsafe(),
// ResizeUnsafe(), or XcapUnsafe(), and the same is true of dst[].
//
// 1. len(src) == len(dst).
//
// 2. All elements of src[] are less than 16.
//
// 3. Capacity of src is at least RoundUpPow2(len(src) + 1, bytesPerVec), and
// the same is true of dst.
//
// 4. The caller does not care if a few bytes past the end of dst[] are
// changed.
func ReverseComp4Unsafe(dst, src []byte) {
	nByte := len(src)
	for idx, invIdx := 0, nByte-1; idx != nByte; idx, invIdx = idx+1, invIdx-1 {
		dst[idx] = revComp4Table[src[invIdx]]
	}
}

// ReverseComp4 saves the reverse-complement of src[] to dst[], assuming .bam
// seq-field encoding with one 4-bit byte per base.
// It panics if len(dst) != len(src).
//
// WARNING: If a src[] value is larger than 15, it's possible for this to
// immediately crash, and it's also possible for this to return and fill src[]
// with garbage.  Only promise is that we don't scribble over arbitrary memory.
func ReverseComp4(dst, src []byte) {
	nByte := len(src)
	if len(dst) != len(src) {
		panic("ReverseComp4() requires len(dst) == len(src).")
	}
	for idx, invIdx := 0, nByte-1; idx != nByte; idx, invIdx = idx+1, invIdx-1 {
		dst[idx] = revComp4Table[src[invIdx]]
	}
}

var revComp2Table = [...]byte{3, 2, 1, 0, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff}

// ReverseComp2UnsafeInplace reverse-complements acgt8[], assuming that it's
// encoded with one byte per base, ACGT=0123.
//
// WARNING: This is a function designed to be used in inner loops, which makes
// assumptions about length and capacity which aren't checked at runtime.  Use
// the safe version of this function when that's a problem.
// These assumptions are always satisfied when the last
// potentially-size-increasing operation on acgt8[] is simd.{Re}makeUnsafe(),
// ResizeUnsafe(), or XcapUnsafe().
//
// 1. Capacity of acgt8[] is at least RoundUpPow2(len(acgt8) + 1, bytesPerVec).
//
// 2. The caller does not care if a few bytes past the end of acgt8[] are
// changed.
func ReverseComp2UnsafeInplace(acgt8 []byte) {
	simd.Reverse8Inplace(acgt8)
	simd.XorConst8Inplace(acgt8, 3)
}

// ReverseComp2Inplace reverse-complements acgt8[], assuming that it's encoded
// with one byte per base, ACGT=0123.
func ReverseComp2Inplace(acgt8 []byte) {
	simd.Reverse8Inplace(acgt8)
	simd.XorConst8Inplace(acgt8, 3)
}

// ReverseComp2Unsafe saves the reverse-complement of src[] to dst[], assuming
// that they're encoded with one byte per base, ACGT=0123.
//
// WARNING: This is a function designed to be used in inner loops, which makes
// assumptions about length and capacity which aren't checked at runtime.  Use
// the safe version of this function when that's a problem.
// Assumptions #2-3 are always satisfied when the last
// potentially-size-increasing operation on src[] is simd.{Re}makeUnsafe(),
// ResizeUnsafe(), or XcapUnsafe(), and the same is true of dst[].
//
// 1. len(src) == len(dst).
//
// 2. Capacity of src is at least RoundUpPow2(len(src) + 1, bytesPerVec), and
// the same is true of dst.
//
// 3. The caller does not care if a few bytes past the end of dst[] are
// changed.
func ReverseComp2Unsafe(dst, src []byte) {
	simd.Reverse8(dst, src)
	simd.XorConst8Inplace(dst, 3)
}

// ReverseComp2 saves the reverse-complement of src[] to dst[], assuming that
// they're encoded with one byte per base, ACGT=0123.
// It panics if len(dst) != len(src).
func ReverseComp2(dst, src []byte) {
	if len(dst) != len(src) {
		panic("ReverseComp2() requires len(dst) == len(src).")
	}
	simd.Reverse8(dst, src)
	simd.XorConst8Inplace(dst, 3)
}
