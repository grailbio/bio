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

func reverseComp8Slow(ascii8 []byte) {
	nByte := len(ascii8)
	nByteDiv2 := nByte >> 1
	for idx, invIdx := 0, nByte-1; idx != nByteDiv2; idx, invIdx = idx+1, invIdx-1 {
		ascii8[idx], ascii8[invIdx] = revComp8Table[ascii8[invIdx]], revComp8Table[ascii8[idx]]
	}
	if nByte&1 == 1 {
		ascii8[nByteDiv2] = revComp8Table[ascii8[nByteDiv2]]
	}
}

var revComp8RandTable = [...]byte{
	'A', 'C', 'G', 'T', 'N', 'a', 'c', 'g', 't', 'n', '0', 0}

func TestReverseComp8(t *testing.T) {
	maxSize := 500
	nIter := 200
	main1Arr := simd.MakeUnsafe(maxSize)
	main2Arr := simd.MakeUnsafe(maxSize)
	main3Arr := simd.MakeUnsafe(maxSize)
	main4Arr := simd.MakeUnsafe(maxSize)
	for iter := 0; iter < nIter; iter++ {
		sliceStart := rand.Intn(maxSize)
		sliceEnd := sliceStart + rand.Intn(maxSize-sliceStart)
		main1Slice := main1Arr[sliceStart:sliceEnd]
		main2Slice := main2Arr[sliceStart:sliceEnd]
		main3Slice := main3Arr[sliceStart:sliceEnd]
		main4Slice := main4Arr[sliceStart:sliceEnd]
		for ii := range main1Slice {
			main1Slice[ii] = revComp8RandTable[rand.Intn(12)]
		}
		copy(main2Slice, main1Slice)
		copy(main3Slice, main1Slice)
		sentinel := byte(rand.Intn(256))
		main3Arr[sliceEnd] = sentinel
		main4Arr[sliceEnd] = sentinel
		biosimd.ReverseComp8NoValidate(main4Slice, main1Slice)
		biosimd.ReverseComp8Inplace(main3Slice)
		reverseComp8Slow(main1Slice)
		biosimd.ReverseComp8InplaceNoValidate(main2Slice)
		if !bytes.Equal(main1Slice, main4Slice) {
			t.Fatal("Mismatched ReverseComp8NoValidate result.")
		}
		if main4Arr[sliceEnd] != sentinel {
			t.Fatal("ReverseComp8NoValidate clobbered an extra byte.")
		}
		if !bytes.Equal(main1Slice, main2Slice) {
			t.Fatal("Mismatched ReverseComp8InplaceNoValidate result.")
		}
		if !bytes.Equal(main1Slice, main3Slice) {
			t.Fatal("Mismatched ReverseComp8Inplace result.")
		}
		if main3Arr[sliceEnd] != sentinel {
			t.Fatal("ReverseComp8Inplace clobbered an extra byte.")
		}
		// Also check ReverseComp8Inplace's validation.
		for ii := range main1Slice {
			main1Slice[ii] = byte(rand.Intn(256))
		}
		copy(main3Slice, main1Slice)
		biosimd.ReverseComp8Inplace(main3Slice)
		reverseComp8Slow(main1Slice)
		if !bytes.Equal(main1Slice, main3Slice) {
			t.Fatal("Mismatched ReverseComp8Inplace result.")
		}
		if main3Arr[sliceEnd] != sentinel {
			t.Fatal("ReverseComp8Inplace clobbered an extra byte.")
		}
	}
}

/*
Benchmark results:
  MacBook Pro (15-inch, 2016)
  2.7 GHz Intel Core i7, 16 GB 2133 MHz LPDDR3

Benchmark_ReverseComp8/SIMDShort1Cpu-8                10         181445752 ns/op
Benchmark_ReverseComp8/SIMDShortHalfCpu-8             30          56466662 ns/op
Benchmark_ReverseComp8/SIMDShortAllCpu-8              20          65764240 ns/op
Benchmark_ReverseComp8/SIMDLong1Cpu-8                  1        1471211448 ns/op
Benchmark_ReverseComp8/SIMDLongHalfCpu-8               2         981752259 ns/op
Benchmark_ReverseComp8/SIMDLongAllCpu-8                1        1012784124 ns/op
Benchmark_ReverseComp8/SlowShort1Cpu-8                 2         764930526 ns/op
Benchmark_ReverseComp8/SlowShortHalfCpu-8              5         212006209 ns/op
Benchmark_ReverseComp8/SlowShortAllCpu-8               5         251716898 ns/op
Benchmark_ReverseComp8/SlowLong1Cpu-8                  1        5852695613 ns/op
Benchmark_ReverseComp8/SlowLongHalfCpu-8               1        1717930116 ns/op
Benchmark_ReverseComp8/SlowLongAllCpu-8                1        1693872482 ns/op
*/

func reverseComp8SimdSubtask(dst, src []byte, nIter int) int {
	for iter := 0; iter < nIter; iter++ {
		biosimd.ReverseComp8Inplace(src)
	}
	return int(src[0])
}

func reverseComp8SlowSubtask(dst, src []byte, nIter int) int {
	for iter := 0; iter < nIter; iter++ {
		reverseComp8Slow(src)
	}
	return int(src[0])
}

func Benchmark_ReverseComp8(b *testing.B) {
	funcs := []taggedMultiBenchFunc{
		{
			f:   reverseComp8SimdSubtask,
			tag: "SIMD",
		},
		{
			f:   reverseComp8SlowSubtask,
			tag: "Slow",
		},
	}
	for _, f := range funcs {
		multiBenchmark(f.f, f.tag+"Short", 0, 150, 9999999, b)
		multiBenchmark(f.f, f.tag+"Long", 0, 249250621, 50, b)
	}
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

/*
Benchmark results:
  MacBook Pro (15-inch, 2016)
  2.7 GHz Intel Core i7, 16 GB 2133 MHz LPDDR3

Benchmark_ReverseComp4/SIMDShort1Cpu-8                20         109340337 ns/op
Benchmark_ReverseComp4/SIMDShortHalfCpu-8             50          38704560 ns/op
Benchmark_ReverseComp4/SIMDShortAllCpu-8              30          45579667 ns/op
Benchmark_ReverseComp4/SIMDLong1Cpu-8                  1        1140730896 ns/op
Benchmark_ReverseComp4/SIMDLongHalfCpu-8               2         974729035 ns/op
Benchmark_ReverseComp4/SIMDLongAllCpu-8                1        1010230155 ns/op
Benchmark_ReverseComp4/SlowShort1Cpu-8                 2         877184805 ns/op
Benchmark_ReverseComp4/SlowShortHalfCpu-8              5         246037663 ns/op
Benchmark_ReverseComp4/SlowShortAllCpu-8               3         370323771 ns/op
Benchmark_ReverseComp4/SlowLong1Cpu-8                  1        6525568174 ns/op
Benchmark_ReverseComp4/SlowLongHalfCpu-8               1        1847395360 ns/op
Benchmark_ReverseComp4/SlowLongAllCpu-8                1        1893590854 ns/op
*/

func reverseComp4SimdSubtask(dst, src []byte, nIter int) int {
	for iter := 0; iter < nIter; iter++ {
		biosimd.ReverseComp4Inplace(src)
	}
	return int(src[0])
}

func reverseComp4SlowSubtask(dst, src []byte, nIter int) int {
	for iter := 0; iter < nIter; iter++ {
		reverseComp4Slow(src)
	}
	return int(src[0])
}

func Benchmark_ReverseComp4(b *testing.B) {
	funcs := []taggedMultiBenchFunc{
		{
			f:   reverseComp4SimdSubtask,
			tag: "SIMD",
		},
		{
			f:   reverseComp4SlowSubtask,
			tag: "Slow",
		},
	}
	opts := multiBenchmarkOpts{
		srcInit: bytesInitMax15,
	}
	for _, f := range funcs {
		multiBenchmark(f.f, f.tag+"Short", 0, 150, 9999999, b, opts)
		multiBenchmark(f.f, f.tag+"Long", 0, 249250621, 50, b, opts)
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
