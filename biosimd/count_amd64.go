// Copyright 2018 GRAIL, Inc.  All rights reserved.
// Use of this source code is governed by the Apache-2.0
// license that can be found in the LICENSE file.

// +build amd64,!appengine

package biosimd

import (
	"reflect"
	"unsafe"

	_ "github.com/grailbio/base/simd" // linkname dependency
)

//go:linkname countNibblesInSetSSE41Asm github.com/grailbio/base/simd.countNibblesInSetSSE41Asm
func countNibblesInSetSSE41Asm(src unsafe.Pointer, tablePtr *[16]byte, nByte int) int

//go:linkname countNibblesInTwoSetsSSE41Asm github.com/grailbio/base/simd.countNibblesInTwoSetsSSE41Asm
func countNibblesInTwoSetsSSE41Asm(cnt2Ptr *int, src unsafe.Pointer, table1Ptr, table2Ptr *[16]byte, nByte int) int

// It's generally cleaner to unpack first, and then use base/simd's byte
// counting functions.  However, we have some preexisting code which doesn't
// work that way, and occasionally this is the *only* operation you need to
// perform on the seq field so you can call these functions directly as an
// optimization.  I don't really want to encourage the latter, though, so I
// don't plan to write unsafe versions of these function(s).

// PackedSeqCount counts the number of .bam base codes in positions
// startPos..(endPos - 1) of seq4 in the given set, where seq4 is in .bam
// packed 4-bit big-endian format.
// The set must be represented as table[x] == 1 when code x is in the set, and
// table[x] == 0 when code x isn't.
// WARNING: This function does not validate the table, startPos, or endPos.  It
// may crash or return a garbage result on invalid input.  (However, it won't
// corrupt memory.)
func PackedSeqCount(seq4 []byte, tablePtr *[16]byte, startPos, endPos int) int {
	if endPos <= startPos {
		return 0
	}
	startOffset := startPos >> 1
	cnt := 0
	if startPos&1 == 1 {
		cnt = int(tablePtr[seq4[startOffset]&15])
		startOffset++
	}
	endOffset := endPos >> 1
	nByte := endOffset - startOffset
	if nByte < 16 {
		for pos := startOffset; pos < endOffset; pos++ {
			seq4Byte := seq4[pos]
			cnt += int(tablePtr[seq4Byte&15] + tablePtr[seq4Byte>>4])
		}
	} else {
		seq4Header := (*reflect.SliceHeader)(unsafe.Pointer(&seq4))
		cnt += countNibblesInSetSSE41Asm(unsafe.Pointer(seq4Header.Data+uintptr(startOffset)), tablePtr, nByte)
	}
	if endPos&1 == 1 {
		cnt += int(tablePtr[seq4[endOffset]>>4])
	}
	return cnt
}

// PackedSeqCountTwo counts the number of .bam base codes in positions
// startPos..(endPos - 1) of seq4 in the given two set, where seq4 is in .bam
// packed 4-bit big-endian format.
// The sets must be represented as table[x] == 1 when code x is in the set, and
// table[x] == 0 when code x isn't.
// WARNING: This function does not validate the tables, startPos, or endPos.
// It may crash or return garbage results on invalid input.  (However, it won't
// corrupt memory.)
func PackedSeqCountTwo(seq4 []byte, table1Ptr, table2Ptr *[16]byte, startPos, endPos int) (int, int) {
	if endPos <= startPos {
		return 0, 0
	}
	startOffset := startPos >> 1
	cnt1, cnt2 := 0, 0
	if startPos&1 == 1 {
		lowBits := seq4[startOffset] & 15
		cnt1 = int(table1Ptr[lowBits])
		cnt2 = int(table2Ptr[lowBits])
		startOffset++
	}
	endOffset := endPos >> 1
	nByte := endOffset - startOffset
	if nByte < 16 {
		for pos := startOffset; pos < endOffset; pos++ {
			seq4Byte := seq4[pos]
			lowBits := seq4Byte & 15
			highBits := seq4Byte >> 4
			cnt1 += int(table1Ptr[lowBits] + table1Ptr[highBits])
			cnt2 += int(table2Ptr[lowBits] + table2Ptr[highBits])
		}
	} else {
		seq4Header := (*reflect.SliceHeader)(unsafe.Pointer(&seq4))
		cnt1 += countNibblesInTwoSetsSSE41Asm(&cnt2, unsafe.Pointer(seq4Header.Data+uintptr(startOffset)), table1Ptr, table2Ptr, nByte)
	}
	if endPos&1 == 1 {
		highBits := seq4[endOffset] >> 4
		cnt1 += int(table1Ptr[highBits])
		cnt2 += int(table2Ptr[highBits])
	}
	return cnt1, cnt2
}
