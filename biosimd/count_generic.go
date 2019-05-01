// Copyright 2018 GRAIL, Inc.  All rights reserved.
// Use of this source code is governed by the Apache-2.0
// license that can be found in the LICENSE file.

// +build !amd64 appengine

package biosimd

// It's generally cleaner to unpack first, and then use base/simd's byte
// counting functions.  However, we have some preexisting code which doesn't
// work that way, and occasionally this is the *only* operation you need to
// perform on the seq field so you can call these functions directly as an
// optimization.  I don't really want to encourage the latter, though, so I
// don't plan to write unsafe versions of these function(s).

// PackedSeqCount counts the number of .bam base codes in positions
// startPos..(endPos - 1) of seq4 in the given set, where seq4 is in .bam
// packed 4-bit big-endian format.
//
// The set must be represented as table[x] == 1 when code x is in the set, and
// table[x] == 0 when code x isn't.
//
// WARNING: This function does not validate the table, startPos, or endPos.  It
// may crash or return a garbage result on invalid input.  (However, it won't
// corrupt memory.)
func PackedSeqCount(seq4 []byte, tablePtr *NibbleLookupTable, startPos, endPos int) int {
	if endPos <= startPos {
		return 0
	}
	startOffset := startPos >> 1
	cnt := 0
	if startPos&1 == 1 {
		cnt = int(tablePtr.Get(seq4[startOffset] & 15))
		startOffset++
	}
	endOffset := endPos >> 1
	for pos := startOffset; pos < endOffset; pos++ {
		seq4Byte := seq4[pos]
		cnt += int(tablePtr.Get(seq4Byte&15) + tablePtr.Get(seq4Byte>>4))
	}
	if endPos&1 == 1 {
		cnt += int(tablePtr.Get(seq4[endOffset] >> 4))
	}
	return cnt
}

// PackedSeqCountTwo counts the number of .bam base codes in positions
// startPos..(endPos - 1) of seq4 in the given two sets, where seq4 is in .bam
// packed 4-bit big-endian format.
//
// The sets must be represented as table[x] == 1 when code x is in the set, and
// table[x] == 0 when code x isn't.
//
// WARNING: This function does not validate the tables, startPos, or endPos.
// It may crash or return garbage results on invalid input.  (However, it won't
// corrupt memory.)
func PackedSeqCountTwo(seq4 []byte, table1Ptr, table2Ptr *NibbleLookupTable, startPos, endPos int) (int, int) {
	if endPos <= startPos {
		return 0, 0
	}
	startOffset := startPos >> 1
	cnt1, cnt2 := 0, 0
	if startPos&1 == 1 {
		lowBits := seq4[startOffset] & 15
		cnt1 = int(table1Ptr.Get(lowBits))
		cnt2 = int(table2Ptr.Get(lowBits))
		startOffset++
	}
	endOffset := endPos >> 1
	for pos := startOffset; pos < endOffset; pos++ {
		seq4Byte := seq4[pos]
		lowBits := seq4Byte & 15
		highBits := seq4Byte >> 4
		cnt1 += int(table1Ptr.Get(lowBits) + table1Ptr.Get(highBits))
		cnt2 += int(table2Ptr.Get(lowBits) + table2Ptr.Get(highBits))
	}
	if endPos&1 == 1 {
		highBits := seq4[endOffset] >> 4
		cnt1 += int(table1Ptr.Get(highBits))
		cnt2 += int(table2Ptr.Get(highBits))
	}
	return cnt1, cnt2
}
