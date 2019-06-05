// Copyright 2018 GRAIL, Inc.  All rights reserved.
// Use of this source code is governed by the Apache-2.0
// license that can be found in the LICENSE file.

package biosimd_test

import (
	"math/rand"
	"testing"
	"unsafe"

	"github.com/grailbio/base/simd"
	"github.com/grailbio/bio/biosimd"
)

func init() {
	if unsafe.Sizeof(uintptr(0)) != 8 {
		// popcnt_amd64.go shouldn't compile at all in this case, but just in
		// case...
		panic("8-byte words required.")
	}
}

var countCGTable = simd.MakeNibbleLookupTable([16]byte{
	0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0})

func packedSeqCountTwoSlow(seq4 []byte, startPos, endPos int, baseCode1, baseCode2 byte) int {
	cnt := 0
	for idx := startPos; idx != endPos; idx++ {
		seqByte := seq4[idx>>1]
		if idx&1 == 0 {
			highBits := seqByte >> 4
			if (highBits == baseCode1) || (highBits == baseCode2) {
				cnt++
			}
		} else {
			lowBits := seqByte & 15
			if (lowBits == baseCode1) || (lowBits == baseCode2) {
				cnt++
			}
		}
	}
	return cnt
}

func TestCountTwo(t *testing.T) {
	maxSize := 10000
	nIter := 200
	srcArr := simd.MakeUnsafe(maxSize)
	var table [16]byte
	for iter := 0; iter < nIter; iter++ {
		sliceStart := rand.Intn(maxSize - 1)
		// guarantee nonempty
		sliceEnd := sliceStart + 1 + rand.Intn(maxSize-1-sliceStart)
		srcSlice := srcArr[sliceStart:sliceEnd]
		for ii := range srcSlice {
			srcSlice[ii] = byte(rand.Intn(256))
		}
		sliceBaseCt := 2 * (sliceEnd - sliceStart)
		startPos := rand.Intn(sliceBaseCt)
		endPos := startPos + rand.Intn(sliceBaseCt-startPos)
		baseCode1 := byte(rand.Intn(15))
		baseCode2 := baseCode1 + 1 + byte(rand.Intn(int(15-baseCode1)))
		table[baseCode1] = 1
		table[baseCode2] = 1
		nlt := simd.MakeNibbleLookupTable(table)

		result1 := packedSeqCountTwoSlow(srcSlice, startPos, endPos, baseCode1, baseCode2)
		result2 := biosimd.PackedSeqCount(srcSlice, &nlt, startPos, endPos)
		if result1 != result2 {
			t.Fatal("Mismatched PackedSeqCount result.")
		}
		table[baseCode1] = 0
		table[baseCode2] = 0
	}
}

/*
Benchmark results:
  MacBook Pro (15-inch, 2016)
  2.7 GHz Intel Core i7, 16 GB 2133 MHz LPDDR3

Benchmark_CountCG/SIMDShort1Cpu-8                     10         133177024 ns/op
Benchmark_CountCG/SIMDShortHalfCpu-8                  50          38857226 ns/op
Benchmark_CountCG/SIMDShortAllCpu-8                   50          33514436 ns/op
Benchmark_CountCG/SIMDLong1Cpu-8                       2         520326888 ns/op
Benchmark_CountCG/SIMDLongHalfCpu-8                    5         258162328 ns/op
Benchmark_CountCG/SIMDLongAllCpu-8                     5         244977674 ns/op
Benchmark_CountCG/SlowShort1Cpu-8                      1        3068588677 ns/op
Benchmark_CountCG/SlowShortHalfCpu-8                   2         811044363 ns/op
Benchmark_CountCG/SlowShortAllCpu-8                    2         619750295 ns/op
Benchmark_CountCG/SlowLong1Cpu-8                       1        25609464123 ns/op
Benchmark_CountCG/SlowLongHalfCpu-8                    1        6957196300 ns/op
Benchmark_CountCG/SlowLongAllCpu-8                     1        5489030920 ns/op
*/

func countCGSimdSubtask(dst, src []byte, nIter int) int {
	tot := 0
	baseCt := len(src) * 2
	for iter := 0; iter < nIter; iter++ {
		tot += biosimd.PackedSeqCount(src, &countCGTable, 0, baseCt)
	}
	return tot
}

func countCGSlowSubtask(dst, src []byte, nIter int) int {
	tot := 0
	baseCt := len(src) * 2
	for iter := 0; iter < nIter; iter++ {
		tot += packedSeqCountTwoSlow(src, 0, baseCt, 2, 4)
	}
	return tot
}

func Benchmark_CountCG(b *testing.B) {
	funcs := []taggedMultiBenchFunc{
		{
			f:   countCGSimdSubtask,
			tag: "SIMD",
		},
		{
			f:   countCGSlowSubtask,
			tag: "Slow",
		},
	}
	for _, f := range funcs {
		multiBenchmark(f.f, f.tag+"Short", 0, 75, 9999999, b)
		multiBenchmark(f.f, f.tag+"Long", 0, 249250622/2, 50, b)
	}
}

func packedSeqCountTwoSets(seq4 []byte, table1Ptr, table2Ptr *simd.NibbleLookupTable, startPos, endPos int) (int, int) {
	cnt1 := 0
	cnt2 := 0
	for idx := startPos; idx != endPos; idx++ {
		seqByte := seq4[idx>>1]
		var curBits byte
		if idx&1 == 0 {
			curBits = seqByte >> 4
		} else {
			curBits = seqByte & 15
		}
		cnt1 += int(table1Ptr.Get(curBits))
		cnt2 += int(table2Ptr.Get(curBits))
	}
	return cnt1, cnt2
}

func TestCountTwoSets(t *testing.T) {
	maxSize := 10000
	nIter := 200
	srcArr := simd.MakeUnsafe(maxSize)
	var table1, table2 [16]byte
	for iter := 0; iter < nIter; iter++ {
		sliceStart := rand.Intn(maxSize - 1)
		// guarantee nonempty
		sliceEnd := sliceStart + 1 + rand.Intn(maxSize-1-sliceStart)
		srcSlice := srcArr[sliceStart:sliceEnd]
		for ii := range srcSlice {
			srcSlice[ii] = byte(rand.Intn(256))
		}
		sliceBaseCt := 2 * (sliceEnd - sliceStart)
		startPos := rand.Intn(sliceBaseCt)
		endPos := startPos + rand.Intn(sliceBaseCt-startPos)
		baseCode1 := byte(rand.Intn(15))
		baseCode2 := baseCode1 + 1 + byte(rand.Intn(int(15-baseCode1)))
		table1[baseCode1] = 1
		table1[baseCode2] = 1

		for ii := 0; ii != 5; ii++ {
			table2[rand.Intn(16)] = 1
		}
		nlt1 := simd.MakeNibbleLookupTable(table1)
		nlt2 := simd.MakeNibbleLookupTable(table2)

		result1a, result1b := packedSeqCountTwoSets(srcSlice, &nlt1, &nlt2, startPos, endPos)
		result2a, result2b := biosimd.PackedSeqCountTwo(srcSlice, &nlt1, &nlt2, startPos, endPos)
		if (result1a != result2a) || (result1b != result2b) {
			t.Fatal("Mismatched PackedSeqCountTwo result.")
		}
		table1[baseCode1] = 0
		table1[baseCode2] = 0
		for pos := range table2 {
			table2[pos] = 0
		}
	}
}
