// Copyright 2018 GRAIL, Inc.  All rights reserved.
// Use of this source code is governed by the Apache-2.0
// license that can be found in the LICENSE file.

package biosimd_test

import (
	"bytes"
	"math/rand"
	"testing"

	"github.com/grailbio/base/simd"
	"github.com/grailbio/bio/biosimd"
)

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

/*
Benchmark results:
  MacBook Pro (15-inch, 2016)
  2.7 GHz Intel Core i7, 16 GB 2133 MHz LPDDR3

Benchmark_UnpackSeq/SIMDShort1Cpu-8                   10         114488272 ns/op
Benchmark_UnpackSeq/SIMDShortHalfCpu-8                50          33101140 ns/op
Benchmark_UnpackSeq/SIMDShortAllCpu-8                 50          29595972 ns/op
Benchmark_UnpackSeq/SIMDLong1Cpu-8                     1        1422335485 ns/op
Benchmark_UnpackSeq/SIMDLongHalfCpu-8                  1        1259798679 ns/op
Benchmark_UnpackSeq/SIMDLongAllCpu-8                   1        1225190627 ns/op
Benchmark_UnpackSeq/SlowShort1Cpu-8                    1        1015049353 ns/op
Benchmark_UnpackSeq/SlowShortHalfCpu-8                 5         261288559 ns/op
Benchmark_UnpackSeq/SlowShortAllCpu-8                  5         264394570 ns/op
Benchmark_UnpackSeq/SlowLong1Cpu-8                     1        7125495274 ns/op
Benchmark_UnpackSeq/SlowLongHalfCpu-8                  1        2041117525 ns/op
Benchmark_UnpackSeq/SlowLongAllCpu-8                   1        2022058870 ns/op
*/

func unpackSeqSimdSubtask(dst, src []byte, nIter int) int {
	for iter := 0; iter < nIter; iter++ {
		biosimd.UnpackSeq(dst, src)
	}
	return int(dst[0])
}

func unpackSeqSlowSubtask(dst, src []byte, nIter int) int {
	for iter := 0; iter < nIter; iter++ {
		unpackSeqSlow(dst, src)
	}
	return int(dst[0])
}

func Benchmark_UnpackSeq(b *testing.B) {
	funcs := []taggedMultiBenchFunc{
		{
			f:   unpackSeqSimdSubtask,
			tag: "SIMD",
		},
		{
			f:   unpackSeqSlowSubtask,
			tag: "Slow",
		},
	}
	for _, f := range funcs {
		multiBenchmark(f.f, f.tag+"Short", 150, 75, 9999999, b)
		multiBenchmark(f.f, f.tag+"Long", 249250621, 249250622/2, 50, b)
	}
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

/*
Benchmark results:
  MacBook Pro (15-inch, 2016)
  2.7 GHz Intel Core i7, 16 GB 2133 MHz LPDDR3

Benchmark_PackSeq/SIMDShort1Cpu-8                     10         129231849 ns/op
Benchmark_PackSeq/SIMDShortHalfCpu-8                  30          35308881 ns/op
Benchmark_PackSeq/SIMDShortAllCpu-8                   50          29446951 ns/op
Benchmark_PackSeq/SIMDLong1Cpu-8                       1        1269761380 ns/op
Benchmark_PackSeq/SIMDLongHalfCpu-8                    2         943575752 ns/op
Benchmark_PackSeq/SIMDLongAllCpu-8                     2         978903937 ns/op
Benchmark_PackSeq/SlowShort1Cpu-8                      2         884202854 ns/op
Benchmark_PackSeq/SlowShortHalfCpu-8                   5         238293524 ns/op
Benchmark_PackSeq/SlowShortAllCpu-8                    5         249022709 ns/op
Benchmark_PackSeq/SlowLong1Cpu-8                       1        6614932340 ns/op
Benchmark_PackSeq/SlowLongHalfCpu-8                    1        1959242488 ns/op
Benchmark_PackSeq/SlowLongAllCpu-8                     1        1906831458 ns/op
*/

func packSeqSimdSubtask(dst, src []byte, nIter int) int {
	for iter := 0; iter < nIter; iter++ {
		biosimd.PackSeq(dst, src)
	}
	return int(dst[0])
}

func packSeqSlowSubtask(dst, src []byte, nIter int) int {
	for iter := 0; iter < nIter; iter++ {
		packSeqSlow(dst, src)
	}
	return int(dst[0])
}

func Benchmark_PackSeq(b *testing.B) {
	funcs := []taggedMultiBenchFunc{
		{
			f:   packSeqSimdSubtask,
			tag: "SIMD",
		},
		{
			f:   packSeqSlowSubtask,
			tag: "Slow",
		},
	}
	opts := multiBenchmarkOpts{
		srcInit: bytesInitMax15,
	}
	for _, f := range funcs {
		multiBenchmark(f.f, f.tag+"Short", 75, 150, 9999999, b, opts)
		multiBenchmark(f.f, f.tag+"Long", 249250622/2, 249250621, 50, b, opts)
	}
}

// No need to benchmark this separately since it's isomorphic to
// simd.PackedNibbleLookup.
func unpackAndReplaceSeqSlow(dst, src []byte, tablePtr *simd.NibbleLookupTable) {
	dstLen := len(dst)
	nSrcFullByte := dstLen / 2
	srcOdd := dstLen & 1
	for srcPos := 0; srcPos < nSrcFullByte; srcPos++ {
		srcByte := src[srcPos]
		dst[2*srcPos] = tablePtr.Get(srcByte >> 4)
		dst[2*srcPos+1] = tablePtr.Get(srcByte & 15)
	}
	if srcOdd == 1 {
		srcByte := src[nSrcFullByte]
		dst[2*nSrcFullByte] = tablePtr.Get(srcByte >> 4)
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

func unpackAndReplaceSeqSubsetSlow(dst, src []byte, tablePtr *simd.NibbleLookupTable, startPos, endPos int) {
	for srcPos := startPos; srcPos != endPos; srcPos++ {
		srcByte := src[srcPos>>1]
		if srcPos&1 == 0 {
			srcByte = srcByte >> 4
		} else {
			srcByte = srcByte & 15
		}
		dst[srcPos-startPos] = tablePtr.Get(srcByte)
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

/*
Benchmark results:
  MacBook Pro (15-inch, 2016)
  2.7 GHz Intel Core i7, 16 GB 2133 MHz LPDDR3

Benchmark_CleanASCIISeq/SIMDShort1Cpu-8                       10         160858434 ns/op
Benchmark_CleanASCIISeq/SIMDShortHalfCpu-8                    30          46693005 ns/op
Benchmark_CleanASCIISeq/SIMDShortAllCpu-8                     30          54947652 ns/op
Benchmark_CleanASCIISeq/SIMDLong1Cpu-8                         1        1253789279 ns/op
Benchmark_CleanASCIISeq/SIMDLongHalfCpu-8                      2         956054465 ns/op
Benchmark_CleanASCIISeq/SIMDLongAllCpu-8                       2         966463331 ns/op
Benchmark_CleanASCIISeq/SlowShort1Cpu-8                        2         786219782 ns/op
Benchmark_CleanASCIISeq/SlowShortHalfCpu-8                     5         216001911 ns/op
Benchmark_CleanASCIISeq/SlowShortAllCpu-8                      5         284746643 ns/op
Benchmark_CleanASCIISeq/SlowLong1Cpu-8                         1        5772484661 ns/op
Benchmark_CleanASCIISeq/SlowLongHalfCpu-8                      1        1661768283 ns/op
Benchmark_CleanASCIISeq/SlowLongAllCpu-8                       1        1655490487 ns/op
*/

func cleanASCIISeqSimdSubtask(dst, src []byte, nIter int) int {
	for iter := 0; iter < nIter; iter++ {
		biosimd.CleanASCIISeqInplace(src)
	}
	return int(src[0])
}

func cleanASCIISeqSlowSubtask(dst, src []byte, nIter int) int {
	for iter := 0; iter < nIter; iter++ {
		cleanASCIISeqSlow(src)
	}
	return int(src[0])
}

func Benchmark_CleanASCIISeq(b *testing.B) {
	funcs := []taggedMultiBenchFunc{
		{
			f:   cleanASCIISeqSimdSubtask,
			tag: "SIMD",
		},
		{
			f:   cleanASCIISeqSlowSubtask,
			tag: "Slow",
		},
	}
	for _, f := range funcs {
		multiBenchmark(f.f, f.tag+"Short", 0, 150, 9999999, b)
		multiBenchmark(f.f, f.tag+"Long", 0, 249250621, 50, b)
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

/*
Benchmark results:
  MacBook Pro (15-inch, 2016)
  2.7 GHz Intel Core i7, 16 GB 2133 MHz LPDDR3

Benchmark_ASCIIToSeq8/SIMDShort1Cpu-8                 10         174509334 ns/op
Benchmark_ASCIIToSeq8/SIMDShortHalfCpu-8              30          50186046 ns/op
Benchmark_ASCIIToSeq8/SIMDShortAllCpu-8               30          43507043 ns/op
Benchmark_ASCIIToSeq8/SIMDLong1Cpu-8                   1        1776279803 ns/op
Benchmark_ASCIIToSeq8/SIMDLongHalfCpu-8                1        1460196089 ns/op
Benchmark_ASCIIToSeq8/SIMDLongAllCpu-8                 1        1558179066 ns/op
Benchmark_ASCIIToSeq8/SlowShort1Cpu-8                  1        1045087004 ns/op
Benchmark_ASCIIToSeq8/SlowShortHalfCpu-8               5         290375515 ns/op
Benchmark_ASCIIToSeq8/SlowShortAllCpu-8                5         264632378 ns/op
Benchmark_ASCIIToSeq8/SlowLong1Cpu-8                   1        7548472579 ns/op
Benchmark_ASCIIToSeq8/SlowLongHalfCpu-8                1        2205978878 ns/op
Benchmark_ASCIIToSeq8/SlowLongAllCpu-8                 1        2205048275 ns/op
*/

func asciiToSeq8SimdSubtask(dst, src []byte, nIter int) int {
	for iter := 0; iter < nIter; iter++ {
		biosimd.ASCIIToSeq8(dst, src)
	}
	return int(dst[0])
}

func asciiToSeq8SlowSubtask(dst, src []byte, nIter int) int {
	for iter := 0; iter < nIter; iter++ {
		asciiToSeq8Slow(dst, src)
	}
	return int(dst[0])
}

func Benchmark_ASCIIToSeq8(b *testing.B) {
	funcs := []taggedMultiBenchFunc{
		{
			f:   asciiToSeq8SimdSubtask,
			tag: "SIMD",
		},
		{
			f:   asciiToSeq8SlowSubtask,
			tag: "Slow",
		},
	}
	for _, f := range funcs {
		multiBenchmark(f.f, f.tag+"Short", 150, 150, 9999999, b)
		multiBenchmark(f.f, f.tag+"Long", 249250621, 249250621, 50, b)
	}
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

/*
Benchmark results:
  MacBook Pro (15-inch, 2016)
  2.7 GHz Intel Core i7, 16 GB 2133 MHz LPDDR3

Benchmark_IsNonACGT/SIMDShort1Cpu-8                   10         119109843 ns/op
Benchmark_IsNonACGT/SIMDShortHalfCpu-8                50          32582132 ns/op
Benchmark_IsNonACGT/SIMDShortAllCpu-8                 50          31520017 ns/op
Benchmark_IsNonACGT/SIMDLong1Cpu-8                     2         523407296 ns/op
Benchmark_IsNonACGT/SIMDLongHalfCpu-8                  5         246869590 ns/op
Benchmark_IsNonACGT/SIMDLongAllCpu-8                   5         242370257 ns/op
Benchmark_IsNonACGT/SlowShort1Cpu-8                    2         859943960 ns/op
Benchmark_IsNonACGT/SlowShortHalfCpu-8                 5         231780298 ns/op
Benchmark_IsNonACGT/SlowShortAllCpu-8                 10         192068918 ns/op
Benchmark_IsNonACGT/SlowLong1Cpu-8                     1        4354582419 ns/op
Benchmark_IsNonACGT/SlowLongHalfCpu-8                  1        1213058951 ns/op
Benchmark_IsNonACGT/SlowLongAllCpu-8                   1        1068068278 ns/op
*/

func isNonACGTSimdSubtask(dst, src []byte, nIter int) int {
	result := true
	for iter := 0; iter < nIter; iter++ {
		result = result && biosimd.IsNonACGTPresent(src)
	}
	if result {
		return int(src[0])
	}
	return int(src[1])
}

func isNonACGTSlowSubtask(dst, src []byte, nIter int) int {
	result := true
	for iter := 0; iter < nIter; iter++ {
		result = result && isNonACGTPresentSlow(src)
	}
	if result {
		return int(src[0])
	}
	return int(src[1])
}

func bytesInitIsNonACGT(src []byte) {
	for i := 0; i < len(src); i++ {
		src[i] = 'T'
	}
	src[len(src)/2] = 'N'
}

func Benchmark_IsNonACGT(b *testing.B) {
	funcs := []taggedMultiBenchFunc{
		{
			f:   isNonACGTSimdSubtask,
			tag: "SIMD",
		},
		{
			f:   isNonACGTSlowSubtask,
			tag: "Slow",
		},
	}
	opts := multiBenchmarkOpts{
		srcInit: bytesInitIsNonACGT,
	}
	for _, f := range funcs {
		multiBenchmark(f.f, f.tag+"Short", 0, 150, 9999999, b, opts)
		multiBenchmark(f.f, f.tag+"Long", 0, 249250621, 50, b, opts)
	}
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

/*
Benchmark results:
  MacBook Pro (15-inch, 2016)
  2.7 GHz Intel Core i7, 16 GB 2133 MHz LPDDR3

Benchmark_ASCIITo2bit/SIMDShort1Cpu-8                 10         167372585 ns/op
Benchmark_ASCIITo2bit/SIMDShortHalfCpu-8              30          47711367 ns/op
Benchmark_ASCIITo2bit/SIMDShortAllCpu-8               50          40750863 ns/op
Benchmark_ASCIITo2bit/SIMDLong1Cpu-8                   1        1100332706 ns/op
Benchmark_ASCIITo2bit/SIMDLongHalfCpu-8                2         745580110 ns/op
Benchmark_ASCIITo2bit/SIMDLongAllCpu-8                 2         733930687 ns/op
Benchmark_ASCIITo2bit/SlowShort1Cpu-8                  2         861736085 ns/op
Benchmark_ASCIITo2bit/SlowShortHalfCpu-8               5         236399713 ns/op
Benchmark_ASCIITo2bit/SlowShortAllCpu-8                5         237820232 ns/op
Benchmark_ASCIITo2bit/SlowLong1Cpu-8                   1        6900630152 ns/op
Benchmark_ASCIITo2bit/SlowLongHalfCpu-8                1        1946458627 ns/op
Benchmark_ASCIITo2bit/SlowLongAllCpu-8                 1        1953776199 ns/op
*/

func asciiTo2bitSimdSubtask(dst, src []byte, nIter int) int {
	for iter := 0; iter < nIter; iter++ {
		biosimd.ASCIITo2bit(dst, src)
	}
	return int(dst[0])
}

func asciiTo2bitSlowSubtask(dst, src []byte, nIter int) int {
	for iter := 0; iter < nIter; iter++ {
		asciiTo2bitSlow(dst, src)
	}
	return int(dst[0])
}

func Benchmark_ASCIITo2bit(b *testing.B) {
	funcs := []taggedMultiBenchFunc{
		{
			f:   asciiTo2bitSimdSubtask,
			tag: "SIMD",
		},
		{
			f:   asciiTo2bitSlowSubtask,
			tag: "Slow",
		},
	}
	for _, f := range funcs {
		// Most of the src characters here aren't in {A,C,G,T}, but that doesn't
		// affect the benchmark results.
		multiBenchmark(f.f, f.tag+"Short", 38, 150, 9999999, b)
		multiBenchmark(f.f, f.tag+"Long", 249250624/4, 249250621, 50, b)
	}
}
