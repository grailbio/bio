// Copyright 2018 GRAIL, Inc.  All rights reserved.
// Use of this source code is governed by the Apache-2.0
// license that can be found in the LICENSE file.

package biosimd_test

import (
	"bytes"
	"math/rand"
	"runtime"
	"testing"

	"github.com/grailbio/base/simd"
	"github.com/grailbio/bio/biosimd"
)

/*
Initial benchmark results:
  MacBook Pro (15-inch, 2016)
  2.7 GHz Intel Core i7, 16 GB 2133 MHz LPDDR3

Benchmark_UnpackSeqShort1-8                   20          76568435 ns/op
Benchmark_UnpackSeqShort4-8                  100          22946637 ns/op
Benchmark_UnpackSeqShortMax-8                100          20510380 ns/op
Benchmark_UnpackSeqLong1-8                     1        1560652622 ns/op
Benchmark_UnpackSeqLong4-8                     1        2169053364 ns/op
Benchmark_UnpackSeqLongMax-8                   1        2882049313 ns/op

Benchmark_PackSeqShort1-8             20          94062516 ns/op
Benchmark_PackSeqShort4-8             50          26069558 ns/op
Benchmark_PackSeqShortMax-8          100          24528628 ns/op
Benchmark_PackSeqLong1-8               1        1510392344 ns/op
Benchmark_PackSeqLong4-8               1        2200776260 ns/op
Benchmark_PackSeqLongMax-8             1        3029837297 ns/op

For comparison, unpackSeqSlow:
Benchmark_UnpackSeqShort1-8                    3         473023326 ns/op
Benchmark_UnpackSeqShort4-8                   10         129047060 ns/op
Benchmark_UnpackSeqShortMax-8                 10         125980303 ns/op
Benchmark_UnpackSeqLong1-8                     1        7138005653 ns/op
Benchmark_UnpackSeqLong4-8                     1        2893149098 ns/op
Benchmark_UnpackSeqLongMax-8                   1        3700028341 ns/op

packSeqSlow:
Benchmark_PackSeqShort1-8              3         480596640 ns/op
Benchmark_PackSeqShort4-8             10         129111468 ns/op
Benchmark_PackSeqShortMax-8           10         118149764 ns/op
Benchmark_PackSeqLong1-8               1        6663558987 ns/op
Benchmark_PackSeqLong4-8               1        2954068774 ns/op
Benchmark_PackSeqLongMax-8             1        4180531216 ns/op
*/

func unpackSeqSubtask(dst, src []byte, nIter int) int {
	for iter := 0; iter < nIter; iter++ {
		biosimd.UnpackSeqUnsafe(dst, src)
	}
	return int(dst[0])
}

func unpackSeqSubtaskFuture(dst, src []byte, nIter int) chan int {
	future := make(chan int)
	go func() { future <- unpackSeqSubtask(dst, src, nIter) }()
	return future
}

func multiUnpackSeq(dsts, srcs [][]byte, cpus int, nJob int) {
	sumFutures := make([]chan int, cpus)
	shardSizeBase := nJob / cpus
	shardRemainder := nJob - shardSizeBase*cpus
	shardSizeP1 := shardSizeBase + 1
	var taskIdx int
	for ; taskIdx < shardRemainder; taskIdx++ {
		sumFutures[taskIdx] = unpackSeqSubtaskFuture(dsts[taskIdx], srcs[taskIdx], shardSizeP1)
	}
	for ; taskIdx < cpus; taskIdx++ {
		sumFutures[taskIdx] = unpackSeqSubtaskFuture(dsts[taskIdx], srcs[taskIdx], shardSizeBase)
	}
	var sum int
	for taskIdx = 0; taskIdx < cpus; taskIdx++ {
		sum += <-sumFutures[taskIdx]
	}
}

func benchmarkUnpackSeq(cpus int, nDstByte int, nJob int, b *testing.B) {
	if cpus > runtime.NumCPU() {
		b.Skipf("only have %v cpus", runtime.NumCPU())
	}

	srcSlices := make([][]byte, cpus)
	dstSlices := make([][]byte, cpus)
	nSrcByte := (nDstByte + 1) >> 1
	for ii := range srcSlices {
		// Add 63 to prevent false sharing.
		newArr := simd.MakeUnsafe(nSrcByte + 63)
		for jj := 0; jj < nSrcByte; jj++ {
			newArr[jj] = byte(jj * 3)
		}
		srcSlices[ii] = newArr[:nSrcByte]
		newArr = simd.MakeUnsafe(nDstByte + 63)
		dstSlices[ii] = newArr[:nDstByte]
	}
	for i := 0; i < b.N; i++ {
		multiUnpackSeq(dstSlices, srcSlices, cpus, nJob)
	}
}

func Benchmark_UnpackSeqShort1(b *testing.B) {
	benchmarkUnpackSeq(1, 75, 9999999, b)
}

func Benchmark_UnpackSeqShort4(b *testing.B) {
	benchmarkUnpackSeq(4, 75, 9999999, b)
}

func Benchmark_UnpackSeqShortMax(b *testing.B) {
	benchmarkUnpackSeq(runtime.NumCPU(), 75, 9999999, b)
}

func Benchmark_UnpackSeqLong1(b *testing.B) {
	benchmarkUnpackSeq(1, 249250621, 50, b)
}

func Benchmark_UnpackSeqLong4(b *testing.B) {
	benchmarkUnpackSeq(4, 249250621, 50, b)
}

func Benchmark_UnpackSeqLongMax(b *testing.B) {
	benchmarkUnpackSeq(runtime.NumCPU(), 249250621, 50, b)
}

func unpackSeqSlow(dst, src []byte) {
	dstLen := len(dst)
	nSrcFullByte := dstLen >> 1
	srcOdd := dstLen & 1
	for srcPos := 0; srcPos < nSrcFullByte; srcPos++ {
		srcByte := src[srcPos]
		dst[2*srcPos] = srcByte >> 4
		dst[2*srcPos+1] = srcByte & 15
	}
	if srcOdd == 1 {
		srcByte := src[nSrcFullByte]
		dst[2*nSrcFullByte] = srcByte >> 4
	}
}

func TestUnpackSeq(t *testing.T) {
	maxDstSize := 500
	maxSrcSize := (maxDstSize + 1) >> 1
	nIter := 200
	srcArr := simd.MakeUnsafe(maxSrcSize)
	dst1Arr := simd.MakeUnsafe(maxDstSize)
	dst2Arr := simd.MakeUnsafe(maxDstSize)
	for iter := 0; iter < nIter; iter++ {
		srcSliceStart := rand.Intn(maxSrcSize)
		dstSliceStart := srcSliceStart * 2
		dstSliceEnd := dstSliceStart + rand.Intn(maxDstSize-dstSliceStart)
		srcSliceEnd := (dstSliceEnd + 1) >> 1
		srcSlice := srcArr[srcSliceStart:srcSliceEnd]
		for ii := range srcSlice {
			srcSlice[ii] = byte(rand.Intn(256))
		}
		dst1Slice := dst1Arr[dstSliceStart:dstSliceEnd]
		dst2Slice := dst2Arr[dstSliceStart:dstSliceEnd]
		unpackSeqSlow(dst1Slice, srcSlice)
		// if bytesPerVec is exported, we should verify that Unsafe functions don't
		// clobber bytes more than that many positions past the slice end.
		biosimd.UnpackSeqUnsafe(dst2Slice, srcSlice)
		if !bytes.Equal(dst1Slice, dst2Slice) {
			t.Fatal("Mismatched UnpackSeqUnsafe result.")
		}
		sentinel := byte(rand.Intn(256))
		dst2Arr[dstSliceEnd] = sentinel
		biosimd.UnpackSeq(dst2Slice, srcSlice)
		if !bytes.Equal(dst1Slice, dst2Slice) {
			t.Fatal("Mismatched UnpackSeq result.")
		}
		if dst2Arr[dstSliceEnd] != sentinel {
			t.Fatal("UnpackSeq clobbered an extra byte.")
		}
	}
}

func packSeqSubtask(dst, src []byte, nIter int) int {
	for iter := 0; iter < nIter; iter++ {
		biosimd.PackSeqUnsafe(dst, src)
	}
	return int(dst[0])
}

func packSeqSubtaskFuture(dst, src []byte, nIter int) chan int {
	future := make(chan int)
	go func() { future <- packSeqSubtask(dst, src, nIter) }()
	return future
}

func multiPackSeq(dsts, srcs [][]byte, cpus int, nJob int) {
	sumFutures := make([]chan int, cpus)
	shardSizeBase := nJob / cpus
	shardRemainder := nJob - shardSizeBase*cpus
	shardSizeP1 := shardSizeBase + 1
	var taskIdx int
	for ; taskIdx < shardRemainder; taskIdx++ {
		sumFutures[taskIdx] = packSeqSubtaskFuture(dsts[taskIdx], srcs[taskIdx], shardSizeP1)
	}
	for ; taskIdx < cpus; taskIdx++ {
		sumFutures[taskIdx] = packSeqSubtaskFuture(dsts[taskIdx], srcs[taskIdx], shardSizeBase)
	}
	var sum int
	for taskIdx = 0; taskIdx < cpus; taskIdx++ {
		sum += <-sumFutures[taskIdx]
	}
}

func benchmarkPackSeq(cpus int, nSrcByte int, nJob int, b *testing.B) {
	if cpus > runtime.NumCPU() {
		b.Skipf("only have %v cpus", runtime.NumCPU())
	}

	srcSlices := make([][]byte, cpus)
	dstSlices := make([][]byte, cpus)
	nDstByte := (nSrcByte + 1) >> 1
	for ii := range srcSlices {
		// Add 63 to prevent false sharing.
		newArr := simd.MakeUnsafe(nSrcByte + 63)
		for jj := 0; jj < nSrcByte; jj++ {
			newArr[jj] = byte(jj*3) & 15
		}
		srcSlices[ii] = newArr[:nSrcByte]
		newArr = simd.MakeUnsafe(nDstByte + 63)
		dstSlices[ii] = newArr[:nDstByte]
	}
	for i := 0; i < b.N; i++ {
		multiPackSeq(dstSlices, srcSlices, cpus, nJob)
	}
}

func Benchmark_PackSeqShort1(b *testing.B) {
	benchmarkPackSeq(1, 75, 9999999, b)
}

func Benchmark_PackSeqShort4(b *testing.B) {
	benchmarkPackSeq(4, 75, 9999999, b)
}

func Benchmark_PackSeqShortMax(b *testing.B) {
	benchmarkPackSeq(runtime.NumCPU(), 75, 9999999, b)
}

func Benchmark_PackSeqLong1(b *testing.B) {
	benchmarkPackSeq(1, 249250621, 50, b)
}

func Benchmark_PackSeqLong4(b *testing.B) {
	benchmarkPackSeq(4, 249250621, 50, b)
}

func Benchmark_PackSeqLongMax(b *testing.B) {
	benchmarkPackSeq(runtime.NumCPU(), 249250621, 50, b)
}

func packSeqSlow(dst, src []byte) {
	srcLen := len(src)
	nDstFullByte := srcLen >> 1
	dstOdd := srcLen & 1
	for dstPos := 0; dstPos < nDstFullByte; dstPos++ {
		dst[dstPos] = (src[2*dstPos] << 4) | src[2*dstPos+1]
	}
	if dstOdd == 1 {
		dst[nDstFullByte] = src[2*nDstFullByte] << 4
	}
}

func TestPackSeq(t *testing.T) {
	maxSrcSize := 500
	maxDstSize := (maxSrcSize + 1) >> 1
	nIter := 200
	srcArr := simd.MakeUnsafe(maxSrcSize)
	dst1Arr := simd.MakeUnsafe(maxDstSize)
	// +1 so we can always append sentinel
	dst2Arr := simd.MakeUnsafe(maxDstSize + 1)
	src2Arr := simd.MakeUnsafe(maxSrcSize)
	for iter := 0; iter < nIter; iter++ {
		dstSliceStart := rand.Intn(maxDstSize)
		srcSliceStart := dstSliceStart * 2
		srcSliceEnd := srcSliceStart + rand.Intn(maxSrcSize-srcSliceStart)
		dstSliceEnd := (srcSliceEnd + 1) >> 1
		srcSlice := srcArr[srcSliceStart:srcSliceEnd]
		for ii := range srcSlice {
			srcSlice[ii] = byte(rand.Intn(16))
		}
		dst1Slice := dst1Arr[dstSliceStart:dstSliceEnd]
		dst2Slice := dst2Arr[dstSliceStart:dstSliceEnd]
		src2Slice := src2Arr[srcSliceStart:srcSliceEnd]
		packSeqSlow(dst1Slice, srcSlice)
		biosimd.PackSeqUnsafe(dst2Slice, srcSlice)
		if !bytes.Equal(dst1Slice, dst2Slice) {
			t.Fatal("Mismatched PackSeqUnsafe result.")
		}
		sentinel := byte(rand.Intn(256))
		dst2Arr[dstSliceEnd] = sentinel
		biosimd.PackSeq(dst2Slice, srcSlice)
		if !bytes.Equal(dst1Slice, dst2Slice) {
			t.Fatal("Mismatched PackSeq result.")
		}
		if dst2Arr[dstSliceEnd] != sentinel {
			t.Fatal("PackSeq clobbered an extra byte.")
		}
		// Verify inverse property.
		biosimd.UnpackSeq(src2Slice, dst1Slice)
		if !bytes.Equal(srcSlice, src2Slice) {
			t.Fatal("UnpackSeq didn't invert PackSeq.")
		}
	}
}
