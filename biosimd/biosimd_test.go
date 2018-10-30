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

Benchmark_UnpackSeqShort1-8                   20          70897086 ns/op
Benchmark_UnpackSeqShort4-8                  100          21312704 ns/op
Benchmark_UnpackSeqShortMax-8                100          18395262 ns/op
Benchmark_UnpackSeqLong1-8                     1        1538266286 ns/op
Benchmark_UnpackSeqLong4-8                     1        2140915576 ns/op
Benchmark_UnpackSeqLongMax-8                   1        2730406285 ns/op

Benchmark_PackSeqShort1-8             20          87414175 ns/op
Benchmark_PackSeqShort4-8             50          24514465 ns/op
Benchmark_PackSeqShortMax-8          100          23399695 ns/op
Benchmark_PackSeqLong1-8               1        1471399081 ns/op
Benchmark_PackSeqLong4-8               1        2160393376 ns/op
Benchmark_PackSeqLongMax-8             1        2973043492 ns/op

Benchmark_CleanASCIISeqShort1-8               20          95413137 ns/op
Benchmark_CleanASCIISeqShort4-8               50          26567655 ns/op
Benchmark_CleanASCIISeqShortMax-8            100          24327826 ns/op
Benchmark_CleanASCIISeqLong1-8                 1        1533053583 ns/op
Benchmark_CleanASCIISeqLong4-8                 1        1982245778 ns/op
Benchmark_CleanASCIISeqLongMax-8               1        2781139905 ns/op

Benchmark_ASCIIToSeq8Short1-8                 10         108897260 ns/op
Benchmark_ASCIIToSeq8Short4-8                 50          30240106 ns/op
Benchmark_ASCIIToSeq8ShortMax-8               50          28269450 ns/op
Benchmark_ASCIIToSeq8Long1-8                   1        2042849647 ns/op
Benchmark_ASCIIToSeq8Long4-8                   1        2866563421 ns/op
Benchmark_ASCIIToSeq8LongMax-8                 1        4069479778 ns/op

Benchmark_IsNonACGTSeqShort1-8                20          68965449 ns/op
Benchmark_IsNonACGTSeqShort4-8               100          19292183 ns/op
Benchmark_IsNonACGTSeqShortMax-8             100          19445680 ns/op
Benchmark_IsNonACGTSeqLong1-8                  2         570726956 ns/op
Benchmark_IsNonACGTSeqLong4-8                  1        1011456304 ns/op
Benchmark_IsNonACGTSeqLongMax-8                1        1498684970 ns/op

Benchmark_ASCIITo2bitShort1-8                 10         141109698 ns/op
Benchmark_ASCIITo2bitShort4-8                 30          44586065 ns/op
Benchmark_ASCIITo2bitShortMax-8               50          34226516 ns/op
Benchmark_ASCIITo2bitLong1-8                   1        1412872064 ns/op
Benchmark_ASCIITo2bitLong4-8                   1        1857122215 ns/op
Benchmark_ASCIITo2bitLongMax-8                 1        2684606937 ns/op

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

cleanASCIISeqSlow:
Benchmark_CleanASCIISeqShort1-8                3         450481328 ns/op
Benchmark_CleanASCIISeqShort4-8               10         122691751 ns/op
Benchmark_CleanASCIISeqShortMax-8             10         158868958 ns/op
Benchmark_CleanASCIISeqLong1-8                 1        6094399462 ns/op
Benchmark_CleanASCIISeqLong4-8                 1        4005568728 ns/op
Benchmark_CleanASCIISeqLongMax-8               1        3286359547 ns/op

asciiToSeq8Slow:
Benchmark_ASCIIToSeq8Short1-8                  2         534821999 ns/op
Benchmark_ASCIIToSeq8Short4-8                 10         145672279 ns/op
Benchmark_ASCIIToSeq8ShortMax-8               10         133403902 ns/op
Benchmark_ASCIIToSeq8Long1-8                   1        8159363086 ns/op
Benchmark_ASCIIToSeq8Long4-8                   1        3625222422 ns/op
Benchmark_ASCIIToSeq8LongMax-8                 1        4613796268 ns/op

isNonACGTPresentSlow:
Benchmark_IsNonACGTSeqShort1-8                 5         311237808 ns/op
Benchmark_IsNonACGTSeqShort4-8                20          87487932 ns/op
Benchmark_IsNonACGTSeqShortMax-8              20          68635003 ns/op
Benchmark_IsNonACGTSeqLong1-8                  1        3158281885 ns/op
Benchmark_IsNonACGTSeqLong4-8                  1        2215643228 ns/op
Benchmark_IsNonACGTSeqLongMax-8                1        2045172556 ns/op

asciiTo2bitSlow:
Benchmark_ASCIITo2bitShort1-8                  3         445481375 ns/op
Benchmark_ASCIITo2bitShort4-8                 10         115023132 ns/op
Benchmark_ASCIITo2bitShortMax-8               10         114890810 ns/op
Benchmark_ASCIITo2bitLong1-8                   1        7284632010 ns/op
Benchmark_ASCIITo2bitLong4-8                   1        3001575126 ns/op
Benchmark_ASCIITo2bitLongMax-8                 1        4445145537 ns/op
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
		unpackAndReplaceSeqSlow(dst1Slice, srcSlice, &biosimd.SeqASCIITable)
		biosimd.UnpackAndReplaceSeqUnsafe(dst2Slice, srcSlice, &biosimd.SeqASCIITable)
		if !bytes.Equal(dst1Slice, dst2Slice) {
			t.Fatal("Mismatched UnpackAndReplaceSeqUnsafe result.")
		}
		simd.Memset8Unsafe(dst2Arr, 0)
		sentinel := byte(rand.Intn(256))
		dst2Arr[dstSliceEnd] = sentinel
		biosimd.UnpackAndReplaceSeq(dst2Slice, srcSlice, &biosimd.SeqASCIITable)
		if !bytes.Equal(dst1Slice, dst2Slice) {
			t.Fatal("Mismatched UnpackAndReplaceSeq result.")
		}
		if dst2Arr[dstSliceEnd] != sentinel {
			t.Fatal("UnpackAndReplaceSeq clobbered an extra byte.")
		}
	}
}

func unpackAndReplaceSeqSubsetSlow(dst, src []byte, tablePtr *[16]byte, startPos, endPos int) {
	for srcPos := startPos; srcPos != endPos; srcPos++ {
		srcByte := src[srcPos>>1]
		if srcPos&1 == 0 {
			srcByte = srcByte >> 4
		} else {
			srcByte = srcByte & 15
		}
		dst[srcPos-startPos] = tablePtr[srcByte]
	}
}

func TestUnpackAndReplaceSeqSubset(t *testing.T) {
	maxDstSize := 500
	maxSrcSize := (maxDstSize + 1) / 2
	nIter := 200
	srcArr := simd.MakeUnsafe(maxSrcSize)
	dst1Arr := simd.MakeUnsafe(maxDstSize)
	dst2Arr := simd.MakeUnsafe(maxDstSize)
	for iter := 0; iter < nIter; iter++ {
		srcSliceStart := rand.Intn(maxSrcSize - 1)
		// Force nonempty.
		srcSliceEnd := srcSliceStart + 1 + rand.Intn(maxSrcSize-1-srcSliceStart)
		srcSlice := srcArr[srcSliceStart:srcSliceEnd]
		for ii := range srcSlice {
			srcSlice[ii] = byte(rand.Intn(256))
		}
		srcSliceLenX2 := 2 * (srcSliceEnd - srcSliceStart)
		startPos := rand.Intn(srcSliceLenX2)
		endPos := startPos + rand.Intn(srcSliceLenX2-startPos)
		dst1Slice := dst1Arr[:endPos-startPos]
		dst2Slice := dst2Arr[:endPos-startPos]
		sentinel := byte(rand.Intn(256))
		dst2Arr[endPos-startPos] = sentinel
		unpackAndReplaceSeqSubsetSlow(dst1Slice, srcSlice, &biosimd.SeqASCIITable, startPos, endPos)
		biosimd.UnpackAndReplaceSeqSubset(dst2Slice, srcSlice, &biosimd.SeqASCIITable, startPos, endPos)
		if !bytes.Equal(dst1Slice, dst2Slice) {
			t.Fatal("Mismatched UnpackAndReplaceSeqSubset result.")
		}
		if dst2Arr[endPos-startPos] != sentinel {
			t.Fatal("UnpackAndReplaceSeqSubset clobbered an extra byte.")
		}
	}
}

func cleanASCIISeqSubtask(ascii8 []byte, nIter int) int {
	for iter := 0; iter < nIter; iter++ {
		biosimd.CleanASCIISeqInplace(ascii8)
	}
	return int(ascii8[0])
}

func cleanASCIISeqSubtaskFuture(ascii8 []byte, nIter int) chan int {
	future := make(chan int)
	go func() { future <- cleanASCIISeqSubtask(ascii8, nIter) }()
	return future
}

func multiCleanASCIISeq(ascii8s [][]byte, cpus int, nJob int) {
	sumFutures := make([]chan int, cpus)
	shardSizeBase := nJob / cpus
	shardRemainder := nJob - shardSizeBase*cpus
	shardSizeP1 := shardSizeBase + 1
	var taskIdx int
	for ; taskIdx < shardRemainder; taskIdx++ {
		sumFutures[taskIdx] = cleanASCIISeqSubtaskFuture(ascii8s[taskIdx], shardSizeP1)
	}
	for ; taskIdx < cpus; taskIdx++ {
		sumFutures[taskIdx] = cleanASCIISeqSubtaskFuture(ascii8s[taskIdx], shardSizeBase)
	}
	var sum int
	for taskIdx = 0; taskIdx < cpus; taskIdx++ {
		sum += <-sumFutures[taskIdx]
	}
}

func benchmarkCleanASCIISeq(cpus int, nByte int, nJob int, b *testing.B) {
	if cpus > runtime.NumCPU() {
		b.Skipf("only have %v cpus", runtime.NumCPU())
	}

	ascii8Slices := make([][]byte, cpus)
	for ii := range ascii8Slices {
		// Add 63 to prevent false sharing.
		newArr := simd.MakeUnsafe(nByte + 63)
		for jj := 0; jj < nByte; jj++ {
			newArr[jj] = byte(jj * 3)
		}
		ascii8Slices[ii] = newArr[:nByte]
	}
	for i := 0; i < b.N; i++ {
		multiCleanASCIISeq(ascii8Slices, cpus, nJob)
	}
}

func Benchmark_CleanASCIISeqShort1(b *testing.B) {
	benchmarkCleanASCIISeq(1, 75, 9999999, b)
}

func Benchmark_CleanASCIISeqShort4(b *testing.B) {
	benchmarkCleanASCIISeq(4, 75, 9999999, b)
}

func Benchmark_CleanASCIISeqShortMax(b *testing.B) {
	benchmarkCleanASCIISeq(runtime.NumCPU(), 75, 9999999, b)
}

func Benchmark_CleanASCIISeqLong1(b *testing.B) {
	benchmarkCleanASCIISeq(1, 249250621, 50, b)
}

func Benchmark_CleanASCIISeqLong4(b *testing.B) {
	benchmarkCleanASCIISeq(4, 249250621, 50, b)
}

func Benchmark_CleanASCIISeqLongMax(b *testing.B) {
	benchmarkCleanASCIISeq(runtime.NumCPU(), 249250621, 50, b)
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

func cleanASCIISeqSlow(ascii8 []byte) {
	for pos, ascii8Byte := range ascii8 {
		ascii8[pos] = cleanASCIISeqTable[ascii8Byte]
	}
}

func TestCleanASCIISeq(t *testing.T) {
	maxSize := 500
	nIter := 200
	main1Arr := simd.MakeUnsafe(maxSize)
	main2Arr := simd.MakeUnsafe(maxSize)
	for iter := 0; iter < nIter; iter++ {
		sliceStart := rand.Intn(maxSize)
		sliceEnd := sliceStart + rand.Intn(maxSize-sliceStart)
		main1Slice := main1Arr[sliceStart:sliceEnd]
		main2Slice := main2Arr[sliceStart:sliceEnd]
		for ii := range main1Slice {
			main1Slice[ii] = byte(rand.Intn(256))
		}
		copy(main2Slice, main1Slice)
		sentinel := byte(rand.Intn(256))
		main2Arr[sliceEnd] = sentinel
		biosimd.CleanASCIISeqInplace(main2Slice)
		cleanASCIISeqSlow(main1Slice)
		if !bytes.Equal(main1Slice, main2Slice) {
			t.Fatal("Mismatched CleanASCIISeqInplace result.")
		}
		if main2Arr[sliceEnd] != sentinel {
			t.Fatal("CleanASCIISeqInplace clobbered an extra byte.")
		}
	}
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

func cleanASCIISeqNoCapitalizeSlow(ascii8 []byte) {
	for pos, ascii8Byte := range ascii8 {
		ascii8[pos] = cleanASCIISeqNoCapitalizeTable[ascii8Byte]
	}
}

func TestCleanASCIISeqNoCapitalize(t *testing.T) {
	maxSize := 500
	nIter := 200
	main1Arr := simd.MakeUnsafe(maxSize)
	main2Arr := simd.MakeUnsafe(maxSize)
	for iter := 0; iter < nIter; iter++ {
		sliceStart := rand.Intn(maxSize)
		sliceEnd := sliceStart + rand.Intn(maxSize-sliceStart)
		main1Slice := main1Arr[sliceStart:sliceEnd]
		main2Slice := main2Arr[sliceStart:sliceEnd]
		for ii := range main1Slice {
			main1Slice[ii] = byte(rand.Intn(256))
		}
		copy(main2Slice, main1Slice)
		sentinel := byte(rand.Intn(256))
		main2Arr[sliceEnd] = sentinel
		biosimd.CleanASCIISeqNoCapitalizeInplace(main2Slice)
		cleanASCIISeqNoCapitalizeSlow(main1Slice)
		if !bytes.Equal(main1Slice, main2Slice) {
			t.Fatal("Mismatched CleanASCIISeqNoCapitalizeInplace result.")
		}
		if main2Arr[sliceEnd] != sentinel {
			t.Fatal("CleanASCIISeqNoCapitalizeInplace clobbered an extra byte.")
		}
	}
}

func asciiToSeq8Subtask(dst, src []byte, nIter int) int {
	for iter := 0; iter < nIter; iter++ {
		biosimd.ASCIIToSeq8(dst, src)
	}
	return int(dst[0])
}

func asciiToSeq8SubtaskFuture(dst, src []byte, nIter int) chan int {
	future := make(chan int)
	go func() { future <- asciiToSeq8Subtask(dst, src, nIter) }()
	return future
}

func multiASCIIToSeq8(dsts, srcs [][]byte, cpus int, nJob int) {
	sumFutures := make([]chan int, cpus)
	shardSizeBase := nJob / cpus
	shardRemainder := nJob - shardSizeBase*cpus
	shardSizeP1 := shardSizeBase + 1
	var taskIdx int
	for ; taskIdx < shardRemainder; taskIdx++ {
		sumFutures[taskIdx] = asciiToSeq8SubtaskFuture(dsts[taskIdx], srcs[taskIdx], shardSizeP1)
	}
	for ; taskIdx < cpus; taskIdx++ {
		sumFutures[taskIdx] = asciiToSeq8SubtaskFuture(dsts[taskIdx], srcs[taskIdx], shardSizeBase)
	}
	var sum int
	for taskIdx = 0; taskIdx < cpus; taskIdx++ {
		sum += <-sumFutures[taskIdx]
	}
}

func benchmarkASCIIToSeq8(cpus int, nByte int, nJob int, b *testing.B) {
	if cpus > runtime.NumCPU() {
		b.Skipf("only have %v cpus", runtime.NumCPU())
	}

	srcSlices := make([][]byte, cpus)
	dstSlices := make([][]byte, cpus)
	for ii := range srcSlices {
		// Add 63 to prevent false sharing.
		newArr := simd.MakeUnsafe(nByte + 63)
		for jj := 0; jj < nByte; jj++ {
			newArr[jj] = byte(jj * 3)
		}
		srcSlices[ii] = newArr[:nByte]
		newArr = simd.MakeUnsafe(nByte + 63)
		dstSlices[ii] = newArr[:nByte]
	}
	for i := 0; i < b.N; i++ {
		multiASCIIToSeq8(dstSlices, srcSlices, cpus, nJob)
	}
}

func Benchmark_ASCIIToSeq8Short1(b *testing.B) {
	benchmarkASCIIToSeq8(1, 75, 9999999, b)
}

func Benchmark_ASCIIToSeq8Short4(b *testing.B) {
	benchmarkASCIIToSeq8(4, 75, 9999999, b)
}

func Benchmark_ASCIIToSeq8ShortMax(b *testing.B) {
	benchmarkASCIIToSeq8(runtime.NumCPU(), 75, 9999999, b)
}

func Benchmark_ASCIIToSeq8Long1(b *testing.B) {
	benchmarkASCIIToSeq8(1, 249250621, 50, b)
}

func Benchmark_ASCIIToSeq8Long4(b *testing.B) {
	benchmarkASCIIToSeq8(4, 249250621, 50, b)
}

func Benchmark_ASCIIToSeq8LongMax(b *testing.B) {
	benchmarkASCIIToSeq8(runtime.NumCPU(), 249250621, 50, b)
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

func asciiToSeq8Slow(dst, src []byte) {
	for pos, srcByte := range src {
		dst[pos] = asciiToSeq8Table[srcByte]
	}
}

func TestASCIIToSeq8(t *testing.T) {
	maxSize := 500
	nIter := 200
	srcArr := simd.MakeUnsafe(maxSize)
	main1Arr := simd.MakeUnsafe(maxSize)
	main2Arr := simd.MakeUnsafe(maxSize)
	for iter := 0; iter < nIter; iter++ {
		sliceStart := rand.Intn(maxSize)
		sliceEnd := sliceStart + rand.Intn(maxSize-sliceStart)
		srcSlice := srcArr[sliceStart:sliceEnd]
		main1Slice := main1Arr[sliceStart:sliceEnd]
		main2Slice := main2Arr[sliceStart:sliceEnd]
		for ii := range srcSlice {
			srcSlice[ii] = byte(rand.Intn(256))
		}
		sentinel := byte(rand.Intn(256))
		main2Arr[sliceEnd] = sentinel
		biosimd.ASCIIToSeq8(main2Slice, srcSlice)
		asciiToSeq8Slow(main1Slice, srcSlice)
		if !bytes.Equal(main1Slice, main2Slice) {
			t.Fatal("Mismatched ASCIIToSeq8 result.")
		}
		if main2Arr[sliceEnd] != sentinel {
			t.Fatal("ASCIIToSeq8 clobbered an extra byte.")
		}
	}
}

func isNonACGTSubtask(ascii8 []byte, nIter int) int {
	result := true
	for iter := 0; iter < nIter; iter++ {
		result = result && biosimd.IsNonACGTPresent(ascii8)
	}
	if result {
		return int(ascii8[0])
	}
	return int(ascii8[1])
}

func isNonACGTSubtaskFuture(ascii8 []byte, nIter int) chan int {
	future := make(chan int)
	go func() { future <- isNonACGTSubtask(ascii8, nIter) }()
	return future
}

func multiIsNonACGTSeq(ascii8s [][]byte, cpus int, nJob int) {
	sumFutures := make([]chan int, cpus)
	shardSizeBase := nJob / cpus
	shardRemainder := nJob - shardSizeBase*cpus
	shardSizeP1 := shardSizeBase + 1
	var taskIdx int
	for ; taskIdx < shardRemainder; taskIdx++ {
		sumFutures[taskIdx] = isNonACGTSubtaskFuture(ascii8s[taskIdx], shardSizeP1)
	}
	for ; taskIdx < cpus; taskIdx++ {
		sumFutures[taskIdx] = isNonACGTSubtaskFuture(ascii8s[taskIdx], shardSizeBase)
	}
	var sum int
	for taskIdx = 0; taskIdx < cpus; taskIdx++ {
		sum += <-sumFutures[taskIdx]
	}
}

func benchmarkIsNonACGTSeq(cpus int, nByte int, nJob int, b *testing.B) {
	if cpus > runtime.NumCPU() {
		b.Skipf("only have %v cpus", runtime.NumCPU())
	}

	ascii8Slices := make([][]byte, cpus)
	for ii := range ascii8Slices {
		// Add 63 to prevent false sharing.
		newArr := simd.MakeUnsafe(nByte + 63)
		for jj := 0; jj < nByte; jj++ {
			newArr[jj] = 'T'
		}
		newArr[nByte/2] = 'N'
		ascii8Slices[ii] = newArr[:nByte]
	}
	for i := 0; i < b.N; i++ {
		multiIsNonACGTSeq(ascii8Slices, cpus, nJob)
	}
}

func Benchmark_IsNonACGTSeqShort1(b *testing.B) {
	benchmarkIsNonACGTSeq(1, 75, 9999999, b)
}

func Benchmark_IsNonACGTSeqShort4(b *testing.B) {
	benchmarkIsNonACGTSeq(4, 75, 9999999, b)
}

func Benchmark_IsNonACGTSeqShortMax(b *testing.B) {
	benchmarkIsNonACGTSeq(runtime.NumCPU(), 75, 9999999, b)
}

func Benchmark_IsNonACGTSeqLong1(b *testing.B) {
	benchmarkIsNonACGTSeq(1, 249250621, 50, b)
}

func Benchmark_IsNonACGTSeqLong4(b *testing.B) {
	benchmarkIsNonACGTSeq(4, 249250621, 50, b)
}

func Benchmark_IsNonACGTSeqLongMax(b *testing.B) {
	benchmarkIsNonACGTSeq(runtime.NumCPU(), 249250621, 50, b)
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

func isNonACGTPresentSlow(ascii8 []byte) bool {
	for _, ascii8Byte := range ascii8 {
		if isNotCapitalACGTTable[ascii8Byte] {
			return true
		}
		// explicit boolean expression is a bit slower
		/*
			if (ascii8Byte != 'A') && (ascii8Byte != 'T') && ((ascii8Byte & 0xfb) != 'C') {
				return true
			}
		*/
	}
	return false
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

func isNonACGTNPresentSlow(ascii8 []byte) bool {
	for _, ascii8Byte := range ascii8 {
		if isNotCapitalACGTNTable[ascii8Byte] {
			return true
		}
	}
	return false
}

var randACGTN0Table = [...]byte{
	'A', 'A', 'A', 'A', 'C', 'C', 'C', 'C', 'G', 'G', 'G', 'G', 'T', 'T', 'T', 'T', 'N', '0'}

func TestIsNonACGTPresent(t *testing.T) {
	maxSize := 500
	nIter := 200
	srcArr := simd.MakeUnsafe(maxSize)
	for iter := 0; iter < nIter; iter++ {
		sliceStart := rand.Intn(maxSize)
		sliceEnd := sliceStart + rand.Intn(maxSize-sliceStart)
		srcSlice := srcArr[sliceStart:sliceEnd]
		for ii := range srcSlice {
			srcSlice[ii] = randACGTN0Table[rand.Intn(18)]
		}
		resultACGT := isNonACGTPresentSlow(srcSlice)
		resultACGT2 := biosimd.IsNonACGTPresent(srcSlice)
		if resultACGT != resultACGT2 {
			t.Fatal("Mismatched IsNonACGTPresent result.")
		}
		resultACGT = isNonACGTNPresentSlow(srcSlice)
		resultACGT2 = biosimd.IsNonACGTNPresent(srcSlice)
		if resultACGT != resultACGT2 {
			t.Fatal("Mismatched IsNonACGTNPresent result.")
		}
	}
}

func asciiTo2bitSubtask(dst, src []byte, nIter int) int {
	for iter := 0; iter < nIter; iter++ {
		biosimd.ASCIITo2bit(dst, src)
	}
	return int(dst[0])
}

func asciiTo2bitSubtaskFuture(dst, src []byte, nIter int) chan int {
	future := make(chan int)
	go func() { future <- asciiTo2bitSubtask(dst, src, nIter) }()
	return future
}

func multiASCIITo2bit(dsts, srcs [][]byte, cpus int, nJob int) {
	sumFutures := make([]chan int, cpus)
	shardSizeBase := nJob / cpus
	shardRemainder := nJob - shardSizeBase*cpus
	shardSizeP1 := shardSizeBase + 1
	var taskIdx int
	for ; taskIdx < shardRemainder; taskIdx++ {
		sumFutures[taskIdx] = asciiTo2bitSubtaskFuture(dsts[taskIdx], srcs[taskIdx], shardSizeP1)
	}
	for ; taskIdx < cpus; taskIdx++ {
		sumFutures[taskIdx] = asciiTo2bitSubtaskFuture(dsts[taskIdx], srcs[taskIdx], shardSizeBase)
	}
	var sum int
	for taskIdx = 0; taskIdx < cpus; taskIdx++ {
		sum += <-sumFutures[taskIdx]
	}
}

func benchmarkASCIITo2bit(cpus int, nSrcByte int, nJob int, b *testing.B) {
	if cpus > runtime.NumCPU() {
		b.Skipf("only have %v cpus", runtime.NumCPU())
	}

	srcSlices := make([][]byte, cpus)
	dstSlices := make([][]byte, cpus)
	nDstByte := (nSrcByte + 3) >> 2
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
		multiASCIITo2bit(dstSlices, srcSlices, cpus, nJob)
	}
}

func Benchmark_ASCIITo2bitShort1(b *testing.B) {
	benchmarkASCIITo2bit(1, 75, 9999999, b)
}

func Benchmark_ASCIITo2bitShort4(b *testing.B) {
	benchmarkASCIITo2bit(4, 75, 9999999, b)
}

func Benchmark_ASCIITo2bitShortMax(b *testing.B) {
	benchmarkASCIITo2bit(runtime.NumCPU(), 75, 9999999, b)
}

func Benchmark_ASCIITo2bitLong1(b *testing.B) {
	benchmarkASCIITo2bit(1, 249250621, 50, b)
}

func Benchmark_ASCIITo2bitLong4(b *testing.B) {
	benchmarkASCIITo2bit(4, 249250621, 50, b)
}

func Benchmark_ASCIITo2bitLongMax(b *testing.B) {
	benchmarkASCIITo2bit(runtime.NumCPU(), 249250621, 50, b)
}

var asciiTo2bitTable = [...]byte{
	0, 0, 0, 1, 3, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 1, 3, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 1, 3, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 1, 3, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 1, 3, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 1, 3, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 1, 3, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 1, 3, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}

func asciiTo2bitSlow(dst, src []byte) {
	srcLen := len(src)
	nDstFullByte := srcLen >> 2
	dstRem := srcLen & 3
	for dstPos := 0; dstPos < nDstFullByte; dstPos++ {
		dst[dstPos] = asciiTo2bitTable[src[4*dstPos]] |
			(asciiTo2bitTable[src[4*dstPos+1]] << 2) |
			(asciiTo2bitTable[src[4*dstPos+2]] << 4) |
			(asciiTo2bitTable[src[4*dstPos+3]] << 6)
	}
	if dstRem != 0 {
		lastByte := asciiTo2bitTable[src[nDstFullByte*4]]
		if dstRem != 1 {
			lastByte |= asciiTo2bitTable[src[nDstFullByte*4+1]] << 2
			if dstRem != 2 {
				lastByte |= asciiTo2bitTable[src[nDstFullByte*4+2]] << 4
			}
		}
		dst[nDstFullByte] = lastByte
	}
}

var twoBitToASCIITable = [...]byte{'A', 'C', 'G', 'T', 'a', 'c', 'g', 't'}

func TestASCIITo2bit(t *testing.T) {
	maxSrcSize := 500
	maxDstSize := (maxSrcSize + 3) >> 2
	nIter := 200
	srcArr := simd.MakeUnsafe(maxSrcSize)
	dst1Arr := simd.MakeUnsafe(maxDstSize)
	// +1 so we can always append sentinel
	dst2Arr := simd.MakeUnsafe(maxDstSize + 1)
	for iter := 0; iter < nIter; iter++ {
		dstSliceStart := rand.Intn(maxDstSize)
		srcSliceStart := dstSliceStart * 4
		srcSliceEnd := srcSliceStart + rand.Intn(maxSrcSize-srcSliceStart)
		dstSliceEnd := (srcSliceEnd + 3) >> 2
		srcSlice := srcArr[srcSliceStart:srcSliceEnd]
		for ii := range srcSlice {
			srcSlice[ii] = twoBitToASCIITable[rand.Intn(8)]
		}
		dst1Slice := dst1Arr[dstSliceStart:dstSliceEnd]
		dst2Slice := dst2Arr[dstSliceStart:dstSliceEnd]
		asciiTo2bitSlow(dst1Slice, srcSlice)
		simd.Memset8Unsafe(dst2Slice, 0)
		sentinel := byte(rand.Intn(256))
		dst2Arr[dstSliceEnd] = sentinel
		biosimd.ASCIITo2bit(dst2Slice, srcSlice)
		if !bytes.Equal(dst1Slice, dst2Slice) {
			t.Fatal("Mismatched ASCIITo2bit result.")
		}
		if dst2Arr[dstSliceEnd] != sentinel {
			t.Fatal("ASCIITo2bit clobbered an extra byte.")
		}
	}
}
