// Copyright 2018 GRAIL, Inc.  All rights reserved.
// Use of this source code is governed by the Apache-2.0
// license that can be found in the LICENSE file.

package biosimd_test

import (
	"math/rand"
	"runtime"
	"testing"
	"unsafe"

	"github.com/grailbio/base/simd"
	"github.com/grailbio/bio/biosimd"
)

/*
Initial benchmark results:
  MacBook Pro (15-inch, 2016)
  2.7 GHz Intel Core i7, 16 GB 2133 MHz LPDDR3

Benchmark_CountCGShort1-8             10         126673237 ns/op
Benchmark_CountCGShort4-8             50          33829967 ns/op
Benchmark_CountCGShortMax-8           50          32196349 ns/op
Benchmark_CountCGLong1-8               1        1158824068 ns/op
Benchmark_CountCGLong4-8               5         274960020 ns/op
Benchmark_CountCGLongMax-8             5         376000031 ns/op

For comparison, packedSeqCountTwoSlow:
Benchmark_CountCGShort1-8              1        2939036568 ns/op
Benchmark_CountCGShort4-8              2         790582580 ns/op
Benchmark_CountCGShortMax-8            2         527043097 ns/op
Benchmark_CountCGLong1-8               1        46681459494 ns/op
Benchmark_CountCGLong4-8               1        24380671980 ns/op
Benchmark_CountCGLongMax-8             1        13357966952 ns/op
*/

func init() {
	if unsafe.Sizeof(uintptr(0)) != 8 {
		// popcnt_amd64.go shouldn't compile at all in this case, but just in
		// case...
		panic("8-byte words required.")
	}
}

var countCGTable = [16]byte{
	0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}

func countCGSubtask(src []byte, nIter int) int {
	tot := 0
	baseCt := len(src) * 2
	for iter := 0; iter < nIter; iter++ {
		tot += biosimd.PackedSeqCount(src, &countCGTable, 0, baseCt)
		// Leave this here to make it easier to switch benchmarks.
		// packedSeqCountTwoSlow(src, 0, baseCt, 2, 4)
	}
	return tot
}

func countCGSubtaskFuture(src []byte, nIter int) chan int {
	future := make(chan int)
	go func() { future <- countCGSubtask(src, nIter) }()
	return future
}

func multiCountCG(srcs [][]byte, cpus int, nJob int) {
	sumFutures := make([]chan int, cpus)
	shardSizeBase := nJob / cpus
	shardRemainder := nJob - shardSizeBase*cpus
	shardSizeP1 := shardSizeBase + 1
	var taskIdx int
	for ; taskIdx < shardRemainder; taskIdx++ {
		sumFutures[taskIdx] = countCGSubtaskFuture(srcs[0], shardSizeP1)
	}
	for ; taskIdx < cpus; taskIdx++ {
		sumFutures[taskIdx] = countCGSubtaskFuture(srcs[0], shardSizeBase)
	}
	var sum int
	for taskIdx = 0; taskIdx < cpus; taskIdx++ {
		sum += <-sumFutures[taskIdx]
	}
}

func benchmarkCountCG(cpus int, nByte int, nJob int, b *testing.B) {
	if cpus > runtime.NumCPU() {
		b.Skipf("only have %v cpus", runtime.NumCPU())
	}

	mainSlices := make([][]byte, 1)
	for ii := range mainSlices {
		// Add 63 to prevent false sharing.
		newArr := simd.MakeUnsafe(nByte + 63)
		mainSlices[ii] = newArr[:nByte]
	}
	for i := 0; i < b.N; i++ {
		multiCountCG(mainSlices, cpus, nJob)
	}
}

// Base sequence in length-150 .bam read occupies 75 bytes, so 75 is a good
// size for the short-array benchmark.
func Benchmark_CountCGShort1(b *testing.B) {
	benchmarkCountCG(1, 75, 9999999, b)
}

func Benchmark_CountCGShort4(b *testing.B) {
	benchmarkCountCG(4, 75, 9999999, b)
}

func Benchmark_CountCGShortMax(b *testing.B) {
	benchmarkCountCG(runtime.NumCPU(), 75, 9999999, b)
}

// GRCh37 chromosome 1 length is 249250621, so that's a plausible long-array
// use case.
func Benchmark_CountCGLong1(b *testing.B) {
	benchmarkCountCG(1, 249250621, 50, b)
}

func Benchmark_CountCGLong4(b *testing.B) {
	benchmarkCountCG(4, 249250621, 50, b)
}

func Benchmark_CountCGLongMax(b *testing.B) {
	benchmarkCountCG(runtime.NumCPU(), 249250621, 50, b)
}

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

		result1 := packedSeqCountTwoSlow(srcSlice, startPos, endPos, baseCode1, baseCode2)
		result2 := biosimd.PackedSeqCount(srcSlice, &table, startPos, endPos)
		if result1 != result2 {
			t.Fatal("Mismatched PackedSeqCount result.")
		}
		table[baseCode1] = 0
		table[baseCode2] = 0
	}
}

func packedSeqCountTwoSets(seq4 []byte, table1Ptr, table2Ptr *[16]byte, startPos, endPos int) (int, int) {
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
		cnt1 += int(table1Ptr[curBits])
		cnt2 += int(table2Ptr[curBits])
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

		result1a, result1b := packedSeqCountTwoSets(srcSlice, &table1, &table2, startPos, endPos)
		result2a, result2b := biosimd.PackedSeqCountTwo(srcSlice, &table1, &table2, startPos, endPos)
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
