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

Benchmark_ReverseComp8Short1-8                10         102511556 ns/op
Benchmark_ReverseComp8Short4-8                50          28367465 ns/op
Benchmark_ReverseComp8ShortMax-8              50          26605330 ns/op
Benchmark_ReverseComp8Long1-8                  1        1652044508 ns/op
Benchmark_ReverseComp8Long4-8                  1        2036898677 ns/op
Benchmark_ReverseComp8LongMax-8                1        2861959800 ns/op

Benchmark_ReverseComp4Short1-8                20          67100763 ns/op
Benchmark_ReverseComp4Short4-8               100          18280864 ns/op
Benchmark_ReverseComp4ShortMax-8             100          17809492 ns/op
Benchmark_ReverseComp4Long1-8                  1        1347194884 ns/op
Benchmark_ReverseComp4Long4-8                  1        1977003084 ns/op
Benchmark_ReverseComp4LongMax-8                1        2772038087 ns/op

For comparison, reverseComp8Slow:
Benchmark_ReverseComp8Short1-8                 3         454187379 ns/op
Benchmark_ReverseComp8Short4-8                10         133079347 ns/op
Benchmark_ReverseComp8ShortMax-8               5         226770101 ns/op
Benchmark_ReverseComp8Long1-8                  1        6430295372 ns/op
Benchmark_ReverseComp8Long4-8                  1        2680017758 ns/op
Benchmark_ReverseComp8LongMax-8                1        3464161375 ns/op

reverseComp4Slow:
Benchmark_ReverseComp4Short1-8                 3         487496815 ns/op
Benchmark_ReverseComp4Short4-8                 5         220447034 ns/op
Benchmark_ReverseComp4ShortMax-8               5         283437486 ns/op
Benchmark_ReverseComp4Long1-8                  1        7214422123 ns/op
Benchmark_ReverseComp4Long4-8                  1        4453820099 ns/op
Benchmark_ReverseComp4LongMax-8                1        3593169766 ns/op
*/

func reverseComp8Subtask(ascii8 []byte, nIter int) int {
	for iter := 0; iter < nIter; iter++ {
		biosimd.ReverseComp8Inplace(ascii8)
	}
	return int(ascii8[0])
}

func reverseComp8SubtaskFuture(main []byte, nIter int) chan int {
	future := make(chan int)
	go func() { future <- reverseComp8Subtask(main, nIter) }()
	return future
}

func multiReverseComp8(mains [][]byte, cpus int, nJob int) {
	sumFutures := make([]chan int, cpus)
	shardSizeBase := nJob / cpus
	shardRemainder := nJob - shardSizeBase*cpus
	shardSizeP1 := shardSizeBase + 1
	var taskIdx int
	for ; taskIdx < shardRemainder; taskIdx++ {
		sumFutures[taskIdx] = reverseComp8SubtaskFuture(mains[taskIdx], shardSizeP1)
	}
	for ; taskIdx < cpus; taskIdx++ {
		sumFutures[taskIdx] = reverseComp8SubtaskFuture(mains[taskIdx], shardSizeBase)
	}
	var sum int
	for taskIdx = 0; taskIdx < cpus; taskIdx++ {
		sum += <-sumFutures[taskIdx]
	}
}

func benchmarkReverseComp8(cpus int, nByte int, nJob int, b *testing.B) {
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
		multiReverseComp8(mainSlices, cpus, nJob)
	}
}

func Benchmark_ReverseComp8Short1(b *testing.B) {
	benchmarkReverseComp8(1, 75, 9999999, b)
}

func Benchmark_ReverseComp8Short4(b *testing.B) {
	benchmarkReverseComp8(4, 75, 9999999, b)
}

func Benchmark_ReverseComp8ShortMax(b *testing.B) {
	benchmarkReverseComp8(runtime.NumCPU(), 75, 9999999, b)
}

func Benchmark_ReverseComp8Long1(b *testing.B) {
	benchmarkReverseComp8(1, 249250621, 50, b)
}

func Benchmark_ReverseComp8Long4(b *testing.B) {
	benchmarkReverseComp8(4, 249250621, 50, b)
}

func Benchmark_ReverseComp8LongMax(b *testing.B) {
	benchmarkReverseComp8(runtime.NumCPU(), 249250621, 50, b)
}

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

func reverseComp4Subtask(seq8 []byte, nIter int) int {
	for iter := 0; iter < nIter; iter++ {
		biosimd.ReverseComp4UnsafeInplace(seq8)
	}
	return int(seq8[0])
}

func reverseComp4SubtaskFuture(main []byte, nIter int) chan int {
	future := make(chan int)
	go func() { future <- reverseComp4Subtask(main, nIter) }()
	return future
}

func multiReverseComp4(mains [][]byte, cpus int, nJob int) {
	sumFutures := make([]chan int, cpus)
	shardSizeBase := nJob / cpus
	shardRemainder := nJob - shardSizeBase*cpus
	shardSizeP1 := shardSizeBase + 1
	var taskIdx int
	for ; taskIdx < shardRemainder; taskIdx++ {
		sumFutures[taskIdx] = reverseComp4SubtaskFuture(mains[taskIdx], shardSizeP1)
	}
	for ; taskIdx < cpus; taskIdx++ {
		sumFutures[taskIdx] = reverseComp4SubtaskFuture(mains[taskIdx], shardSizeBase)
	}
	var sum int
	for taskIdx = 0; taskIdx < cpus; taskIdx++ {
		sum += <-sumFutures[taskIdx]
	}
}

func benchmarkReverseComp4(cpus int, nByte int, nJob int, b *testing.B) {
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
		multiReverseComp4(mainSlices, cpus, nJob)
	}
}

func Benchmark_ReverseComp4Short1(b *testing.B) {
	benchmarkReverseComp4(1, 75, 9999999, b)
}

func Benchmark_ReverseComp4Short4(b *testing.B) {
	benchmarkReverseComp4(4, 75, 9999999, b)
}

func Benchmark_ReverseComp4ShortMax(b *testing.B) {
	benchmarkReverseComp4(runtime.NumCPU(), 75, 9999999, b)
}

func Benchmark_ReverseComp4Long1(b *testing.B) {
	benchmarkReverseComp4(1, 249250621, 50, b)
}

func Benchmark_ReverseComp4Long4(b *testing.B) {
	benchmarkReverseComp4(4, 249250621, 50, b)
}

func Benchmark_ReverseComp4LongMax(b *testing.B) {
	benchmarkReverseComp4(runtime.NumCPU(), 249250621, 50, b)
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
