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

Benchmark_ReverseCompShort1-8                 20          76847905 ns/op
Benchmark_ReverseCompShort4-8                 50          29293746 ns/op
Benchmark_ReverseCompShortMax-8              100          19428703 ns/op
Benchmark_ReverseCompLong1-8                   1        1458896247 ns/op
Benchmark_ReverseCompLong4-8                   1        2012557816 ns/op
Benchmark_ReverseCompLongMax-8                 1        2835188762 ns/op

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

reverseComp4Slow:
Benchmark_ReverseCompShort1-8                  3         487496815 ns/op
Benchmark_ReverseCompShort4-8                  5         220447034 ns/op
Benchmark_ReverseCompShortMax-8                5         283437486 ns/op
Benchmark_ReverseCompLong1-8                   1        7214422123 ns/op
Benchmark_ReverseCompLong4-8                   1        4453820099 ns/op
Benchmark_ReverseCompLongMax-8                 1        3593169766 ns/op
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
		simd.Memset8Unsafe(dst2Slice, 0)
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
		simd.Memset8Unsafe(dst2Slice, 0)
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

func reverseCompSubtask(seq8 []byte, nIter int) int {
	for iter := 0; iter < nIter; iter++ {
		biosimd.ReverseComp4UnsafeInplace(seq8)
	}
	return int(seq8[0])
}

func reverseCompSubtaskFuture(main []byte, nIter int) chan int {
	future := make(chan int)
	go func() { future <- reverseCompSubtask(main, nIter) }()
	return future
}

func multiReverseComp(mains [][]byte, cpus int, nJob int) {
	sumFutures := make([]chan int, cpus)
	shardSizeBase := nJob / cpus
	shardRemainder := nJob - shardSizeBase*cpus
	shardSizeP1 := shardSizeBase + 1
	var taskIdx int
	for ; taskIdx < shardRemainder; taskIdx++ {
		sumFutures[taskIdx] = reverseCompSubtaskFuture(mains[taskIdx], shardSizeP1)
	}
	for ; taskIdx < cpus; taskIdx++ {
		sumFutures[taskIdx] = reverseCompSubtaskFuture(mains[taskIdx], shardSizeBase)
	}
	var sum int
	for taskIdx = 0; taskIdx < cpus; taskIdx++ {
		sum += <-sumFutures[taskIdx]
	}
}

func benchmarkReverseComp(cpus int, nByte int, nJob int, b *testing.B) {
	if cpus > runtime.NumCPU() {
		b.Skipf("only have %v cpus", runtime.NumCPU())
	}

	mainSlices := make([][]byte, cpus)
	for ii := range mainSlices {
		// Add 63 to prevent false sharing.
		newArr := simd.MakeUnsafe(nByte + 63)
		for jj := 0; jj < nByte; jj++ {
			newArr[jj] = byte(jj*3) & 15
		}
		mainSlices[ii] = newArr[:nByte]
	}
	for i := 0; i < b.N; i++ {
		multiReverseComp(mainSlices, cpus, nJob)
	}
}

func Benchmark_ReverseCompShort1(b *testing.B) {
	benchmarkReverseComp(1, 75, 9999999, b)
}

func Benchmark_ReverseCompShort4(b *testing.B) {
	benchmarkReverseComp(4, 75, 9999999, b)
}

func Benchmark_ReverseCompShortMax(b *testing.B) {
	benchmarkReverseComp(runtime.NumCPU(), 75, 9999999, b)
}

func Benchmark_ReverseCompLong1(b *testing.B) {
	benchmarkReverseComp(1, 249250621, 50, b)
}

func Benchmark_ReverseCompLong4(b *testing.B) {
	benchmarkReverseComp(4, 249250621, 50, b)
}

func Benchmark_ReverseCompLongMax(b *testing.B) {
	benchmarkReverseComp(runtime.NumCPU(), 249250621, 50, b)
}

var revComp4Table = [...]byte{0, 8, 4, 12, 2, 10, 6, 14, 1, 9, 5, 13, 3, 11, 7, 15}

func reverseComp4Slow(seq8 []byte) {
	nByte := len(seq8)
	nByteDiv2 := nByte >> 1
	for idx, invIdx := 0, nByte-1; idx != nByteDiv2; idx, invIdx = idx+1, invIdx-1 {
		seq8[idx], seq8[invIdx] = revComp4Table[seq8[invIdx]], revComp4Table[seq8[idx]]
	}
	if nByte&1 == 1 {
		seq8[nByteDiv2] = revComp4Table[seq8[nByteDiv2]]
	}
}

func TestReverseComp4(t *testing.T) {
	maxSize := 500
	nIter := 200
	main1Arr := simd.MakeUnsafe(maxSize)
	main2Arr := simd.MakeUnsafe(maxSize)
	main3Arr := simd.MakeUnsafe(maxSize)
	main4Arr := simd.MakeUnsafe(maxSize)
	main5Arr := simd.MakeUnsafe(maxSize)
	for iter := 0; iter < nIter; iter++ {
		sliceStart := rand.Intn(maxSize)
		sliceEnd := sliceStart + rand.Intn(maxSize-sliceStart)
		main1Slice := main1Arr[sliceStart:sliceEnd]
		main2Slice := main2Arr[sliceStart:sliceEnd]
		main3Slice := main3Arr[sliceStart:sliceEnd]
		main4Slice := main4Arr[sliceStart:sliceEnd]
		main5Slice := main5Arr[sliceStart:sliceEnd]
		for ii := range main1Slice {
			main1Slice[ii] = byte(rand.Intn(16))
		}
		copy(main3Slice, main1Slice)
		sentinel := byte(rand.Intn(256))
		main3Arr[sliceEnd] = sentinel
		main5Arr[sliceEnd] = sentinel
		biosimd.ReverseComp4Unsafe(main4Slice, main1Slice)
		biosimd.ReverseComp4Unsafe(main2Slice, main4Slice)
		if !bytes.Equal(main1Slice, main2Slice) {
			t.Fatal("ReverseComp4Unsafe isn't its own inverse.")
		}
		copy(main2Slice, main1Slice)
		biosimd.ReverseComp4(main5Slice, main1Slice)
		reverseComp4Slow(main1Slice)
		biosimd.ReverseComp4UnsafeInplace(main2Slice)
		biosimd.ReverseComp4Inplace(main3Slice)
		if !bytes.Equal(main1Slice, main2Slice) {
			t.Fatal("Mismatched ReverseComp4UnsafeInplace result.")
		}
		if !bytes.Equal(main1Slice, main3Slice) {
			t.Fatal("Mismatched ReverseComp4Inplace result.")
		}
		if main3Arr[sliceEnd] != sentinel {
			t.Fatal("ReverseComp4Inplace clobbered an extra byte.")
		}
		if !bytes.Equal(main1Slice, main4Slice) {
			t.Fatal("Mismatched ReverseComp4Unsafe result.")
		}
		if !bytes.Equal(main1Slice, main5Slice) {
			t.Fatal("Mismatched ReverseComp4 result.")
		}
		if main5Arr[sliceEnd] != sentinel {
			t.Fatal("ReverseComp4 clobbered an extra byte.")
		}
	}
}

func reverseComp2Slow(main []byte) {
	nByte := len(main)
	nByteDiv2 := nByte >> 1
	for idx, invIdx := 0, nByte-1; idx != nByteDiv2; idx, invIdx = idx+1, invIdx-1 {
		main[idx], main[invIdx] = 3-main[invIdx], 3-main[idx]
	}
	if nByte&1 == 1 {
		main[nByteDiv2] = 3 - main[nByteDiv2]
	}
}

func TestReverseComp2(t *testing.T) {
	maxSize := 500
	nIter := 200
	main1Arr := simd.MakeUnsafe(maxSize)
	main2Arr := simd.MakeUnsafe(maxSize)
	main3Arr := simd.MakeUnsafe(maxSize)
	main4Arr := simd.MakeUnsafe(maxSize)
	main5Arr := simd.MakeUnsafe(maxSize)
	for iter := 0; iter < nIter; iter++ {
		sliceStart := rand.Intn(maxSize)
		sliceEnd := sliceStart + rand.Intn(maxSize-sliceStart)
		main1Slice := main1Arr[sliceStart:sliceEnd]
		main2Slice := main2Arr[sliceStart:sliceEnd]
		main3Slice := main3Arr[sliceStart:sliceEnd]
		main4Slice := main4Arr[sliceStart:sliceEnd]
		main5Slice := main5Arr[sliceStart:sliceEnd]
		for ii := range main1Slice {
			main1Slice[ii] = byte(rand.Intn(4))
		}
		copy(main3Slice, main1Slice)
		sentinel := byte(rand.Intn(256))
		main3Arr[sliceEnd] = sentinel
		main5Arr[sliceEnd] = sentinel
		biosimd.ReverseComp2Unsafe(main4Slice, main1Slice)
		biosimd.ReverseComp2Unsafe(main2Slice, main4Slice)
		if !bytes.Equal(main1Slice, main2Slice) {
			t.Fatal("ReverseComp2Unsafe isn't its own inverse.")
		}
		copy(main2Slice, main1Slice)
		biosimd.ReverseComp2(main5Slice, main1Slice)
		reverseComp2Slow(main1Slice)
		biosimd.ReverseComp2UnsafeInplace(main2Slice)
		biosimd.ReverseComp2Inplace(main3Slice)
		if !bytes.Equal(main1Slice, main2Slice) {
			t.Fatal("Mismatched ReverseComp2UnsafeInplace result.")
		}
		if !bytes.Equal(main1Slice, main3Slice) {
			t.Fatal("Mismatched ReverseComp2Inplace result.")
		}
		if main3Arr[sliceEnd] != sentinel {
			t.Fatal("ReverseComp2Inplace clobbered an extra byte.")
		}
		if !bytes.Equal(main1Slice, main4Slice) {
			t.Fatal("Mismatched ReverseComp2Unsafe result.")
		}
		if !bytes.Equal(main1Slice, main5Slice) {
			t.Fatal("Mismatched ReverseComp2 result.")
		}
		if main5Arr[sliceEnd] != sentinel {
			t.Fatal("ReverseComp2 clobbered an extra byte.")
		}
	}
}

// No need to benchmark this separately since it's isomorphic to
// simd.PackedNibbleLookup.
func unpackAndReplaceSeqSlow(dst, src []byte, tablePtr *[16]byte) {
	dstLen := len(dst)
	nSrcFullByte := dstLen / 2
	srcOdd := dstLen & 1
	for srcPos := 0; srcPos < nSrcFullByte; srcPos++ {
		srcByte := src[srcPos]
		dst[2*srcPos] = tablePtr[srcByte>>4]
		dst[2*srcPos+1] = tablePtr[srcByte&15]
	}
	if srcOdd == 1 {
		srcByte := src[nSrcFullByte]
		dst[2*nSrcFullByte] = tablePtr[srcByte>>4]
	}
}

func TestUnpackAndReplaceSeq(t *testing.T) {
	maxDstSize := 500
	maxSrcSize := (maxDstSize + 1) / 2
	nIter := 200
	srcArr := simd.MakeUnsafe(maxSrcSize)
	dst1Arr := simd.MakeUnsafe(maxDstSize)
	dst2Arr := simd.MakeUnsafe(maxDstSize)
	table := [...]byte{'=', 'A', 'C', 'M', 'G', 'R', 'S', 'V', 'T', 'W', 'Y', 'H', 'K', 'D', 'B', 'N'}
	for iter := 0; iter < nIter; iter++ {
		srcSliceStart := rand.Intn(maxSrcSize)
		dstSliceStart := srcSliceStart * 2
		dstSliceEnd := dstSliceStart + rand.Intn(maxDstSize-dstSliceStart)
		srcSliceEnd := (dstSliceEnd + 1) / 2
		srcSlice := srcArr[srcSliceStart:srcSliceEnd]
		for ii := range srcSlice {
			srcSlice[ii] = byte(rand.Intn(256))
		}
		dst1Slice := dst1Arr[dstSliceStart:dstSliceEnd]
		dst2Slice := dst2Arr[dstSliceStart:dstSliceEnd]
		unpackAndReplaceSeqSlow(dst1Slice, srcSlice, &table)
		biosimd.UnpackAndReplaceSeqUnsafe(dst2Slice, srcSlice, &table)
		if !bytes.Equal(dst1Slice, dst2Slice) {
			t.Fatal("Mismatched UnpackAndReplaceSeqUnsafe result.")
		}
		simd.Memset8Unsafe(dst2Arr, 0)
		sentinel := byte(rand.Intn(256))
		dst2Arr[dstSliceEnd] = sentinel
		biosimd.UnpackAndReplaceSeq(dst2Slice, srcSlice, &table)
		if !bytes.Equal(dst1Slice, dst2Slice) {
			t.Fatal("Mismatched UnpackAndReplaceSeq result.")
		}
		if dst2Arr[dstSliceEnd] != sentinel {
			t.Fatal("UnpackAndReplaceSeq clobbered an extra byte.")
		}
	}
}
