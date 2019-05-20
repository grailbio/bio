// Copyright 2019 GRAIL, Inc.  All rights reserved.
// Use of this source code is governed by the Apache-2.0
// license that can be found in the LICENSE file.

package biosimd_test

import (
	"fmt"
	"math/rand"
	"runtime"
	"testing"

	"github.com/grailbio/base/simd"
	"github.com/grailbio/bio/biosimd"
	"github.com/grailbio/testutil/assert"
)

/*
Initial benchmark results:
  MacBook Pro (15-inch, 2016)
  2.7 GHz Intel Core i7, 16 GB 2133 MHz LPDDR3

Benchmark_FastqRender1-8              10         172173636 ns/op
Benchmark_FastqRender4-8              20          66978727 ns/op
Benchmark_FastqRenderMax-8            30          49230571 ns/op

With two simd.PackedNibbleLookup calls instead of the specialized function:
Benchmark_FastqRender1-8               5         306911631 ns/op
Benchmark_FastqRender4-8              20          99318550 ns/op
Benchmark_FastqRenderMax-8            20          78737149 ns/op

*/

var baseTable = simd.MakeNibbleLookupTable([16]byte{
	'N', '?', '?', '?', 'A', 'C', 'G', 'T', 'A', 'C', 'G', 'T', 'A', 'C', 'G', 'T'})
var qualTable = simd.MakeNibbleLookupTable([16]byte{
	'#', '#', '#', '#', ',', ',', ',', ',', ':', ':', ':', ':', 'F', 'F', 'F', 'F'})

var asciiToBaseBitsTable = [...]byte{
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}

// Not optimized.  Might want to add an optimized version of this to the main
// package later.
func packAsciiBasesAndQuals(bases, quals, qualEncoding []byte) (dst []byte, err error) {
	nBase := len(bases)
	if nBase != len(quals) {
		err = fmt.Errorf("packAsciiBasesAndQuals: inconsistent len(bases) and len(quals)")
		return
	}
	if len(qualEncoding) != 4 {
		err = fmt.Errorf("packAsciiBasesAndQuals: unexpected len(qualEncoding)")
		return
	}

	var asciiToBase16Qual [256]byte
	for i := range asciiToBase16Qual {
		// Use this value to indicate "unexpected qual".
		asciiToBase16Qual[i] = 255
	}
	for i, q := range qualEncoding {
		// left-shift 2 since qual is stored in the high 2 bits of each nibble
		asciiToBase16Qual[q+33] = byte(i << 2)
	}

	dst = make([]byte, (nBase+1)/2)
	for i := 0; i != nBase; i++ {
		packedBaseAndQual := asciiToBase16Qual[quals[i]]
		if packedBaseAndQual == 255 {
			err = fmt.Errorf("packAsciiBasesAndQuals: unexpected qual value")
			return
		}
		packedBaseAndQual |= asciiToBaseBitsTable[bases[i]]
		if i&1 == 0 {
			dst[i/2] = packedBaseAndQual
		} else {
			dst[i/2] |= packedBaseAndQual << 4
		}
	}
	return
}

var baseChars = []byte{'A', 'C', 'G', 'T', 'N'}

func TestFillFastqRecordBody(t *testing.T) {
	// A few fixed test cases.
	qualEncoding := []byte{2, 11, 25, 37}
	tests := []struct {
		baseAscii string
		qualAscii string
	}{
		{
			baseAscii: "",
			qualAscii: "",
		},
		{
			baseAscii: "G",
			qualAscii: ",",
		},
		{
			baseAscii: "N",
			qualAscii: "#",
		},
		{
			baseAscii: "ACACNGGAGAGCTTTTTTACA",
			qualAscii: ",,,F#FFFFFFFF::FFFFF,",
		},
		{
			baseAscii: "CAAAAAAAAAAAAAAAAAAAAAAAAAAAATTTAGGATTATAGGCCCCC",
			qualAscii: "FFFFFFFFFFFFFFFFFFFFF,FFFFFF:FFFFFFFFFFFFFFF,,FF",
		},
		{
			baseAscii: "AGAGTGACTCGACTACTCACCGCGAACAAAAAAAAAAAATACATAGCATCGAGCAGCTACGACACGATCGATCGATCGACTAGTCAGTCGACTCGACTGGGGTGCTAAAAAAAAAAAAAAAAAGATGTCAGCATCGATCGCGGGGGAACCC",
			qualAscii: ",,,,,,,,,,,,,,,,,,,FFFFFFFFFFFFFFFFFFFFF::::::::::::::::::F:FFF:F:FFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFF::FFFFFFFFFFF:FFFFFFFFFFF:FFFFFFFFFF",
		},
	}
	for _, tc := range tests {
		b16, err := packAsciiBasesAndQuals([]byte(tc.baseAscii), []byte(tc.qualAscii), qualEncoding)
		if err != nil {
			t.Fatalf("packAsciiBasesAndQuals error: %v", err)
		}
		readLen := len(tc.baseAscii)
		dst := make([]byte, 2*readLen+4)
		biosimd.FillFastqRecordBodyFromNibbles(dst, b16, readLen, &baseTable, &qualTable)
		assert.EQ(t, dst[:readLen], []byte(tc.baseAscii))
		assert.EQ(t, dst[readLen+3:2*readLen+3], []byte(tc.qualAscii))
	}
	// Random test cases.
	maxReadLen := 151
	nIter := 2000
	rand.Seed(1)
	for iter := 0; iter < nIter; iter++ {
		readLen := 1 + rand.Intn(maxReadLen)
		baseAscii := make([]byte, readLen)
		qualAscii := make([]byte, readLen)
		for i := range baseAscii {
			baseType := rand.Intn(5)
			baseAscii[i] = baseChars[baseType]
			if baseType == 4 {
				qualAscii[i] = '#'
			} else {
				qualAscii[i] = 33 + qualEncoding[1+rand.Intn(3)]
			}
		}
		b16, err := packAsciiBasesAndQuals(baseAscii, qualAscii, qualEncoding)
		if err != nil {
			t.Fatalf("packAsciiBasesAndQuals error: %v", err)
		}
		dst := make([]byte, 2*readLen+4)
		biosimd.FillFastqRecordBodyFromNibbles(dst, b16, readLen, &baseTable, &qualTable)
		assert.EQ(t, dst[:readLen], baseAscii)
		assert.EQ(t, dst[readLen+3:2*readLen+3], qualAscii)
	}
}

func fastqRenderSubtask(dst, b16 []byte, nIter int) int {
	readLen := (len(dst) - 4) >> 1
	for iter := 0; iter < nIter; iter++ {
		biosimd.FillFastqRecordBodyFromNibbles(dst, b16, readLen, &baseTable, &qualTable)
	}
	return int(dst[0])
}

/*
// This represents the code that would be used without a specialized
// FillFastqRecordBodyFromNibbles function.
func fastqRenderSubtask(dst, b16 []byte, nIter int) int {
	readLen := (len(dst) - 4) >> 1
	for iter := 0; iter < nIter; iter++ {
		simd.PackedNibbleLookup(dst[:readLen], b16, &baseTable)
		copy(dst[readLen:readLen+3], "\n+\n")
		simd.PackedNibbleLookup(dst[readLen+3:2*readLen+3], b16, &qualTable)
		dst[2*readLen+3] = '\n'
	}
	return int(dst[0])
}
*/

func fastqRenderSubtaskFuture(dst, src []byte, nIter int) chan int {
	future := make(chan int)
	go func() { future <- fastqRenderSubtask(dst, src, nIter) }()
	return future
}

func multiFastqRender(dsts, srcs [][]byte, cpus int, nJob int) {
	sumFutures := make([]chan int, cpus)
	shardSizeBase := nJob / cpus
	shardRemainder := nJob - shardSizeBase*cpus
	shardSizeP1 := shardSizeBase + 1
	var taskIdx int
	for ; taskIdx < shardRemainder; taskIdx++ {
		sumFutures[taskIdx] = fastqRenderSubtaskFuture(dsts[taskIdx], srcs[taskIdx], shardSizeP1)
	}
	for ; taskIdx < cpus; taskIdx++ {
		sumFutures[taskIdx] = fastqRenderSubtaskFuture(dsts[taskIdx], srcs[taskIdx], shardSizeBase)
	}
	var sum int
	for taskIdx = 0; taskIdx < cpus; taskIdx++ {
		sum += <-sumFutures[taskIdx]
	}
}

func benchmarkFastqRender(cpus int, readLen int, nJob int, b *testing.B) {
	if cpus > runtime.NumCPU() {
		b.Skipf("only have %v cpus", runtime.NumCPU())
	}

	srcSlices := make([][]byte, cpus)
	dstSlices := make([][]byte, cpus)
	nSrcByte := (readLen + 1) / 2
	nDstByte := 2*readLen + 4
	for ii := range srcSlices {
		// Add 63 to prevent false sharing.
		newArr := simd.MakeUnsafe(nSrcByte + 63)
		for jj := 0; jj < nSrcByte; jj++ {
			b := byte(jj * 3)
			// Doesn't make a difference in the benchmark, but may as well enforce
			// nibbles 1-3 being invalid.
			if ((b & 15) >= 1) && ((b & 15) < 4) {
				b = b & 252
			}
			if ((b & 240) >= 16) && ((b & 240) < 64) {
				b = b & 207
			}
			newArr[jj] = b
		}
		srcSlices[ii] = newArr[:nSrcByte]
		newArr = simd.MakeUnsafe(nDstByte + 63)
		dstSlices[ii] = newArr[:nDstByte]
	}
	for i := 0; i < b.N; i++ {
		multiFastqRender(dstSlices, srcSlices, cpus, nJob)
	}
}

func Benchmark_FastqRender1(b *testing.B) {
	benchmarkFastqRender(1, 151, 9999999, b)
}

func Benchmark_FastqRender4(b *testing.B) {
	benchmarkFastqRender(4, 151, 9999999, b)
}

func Benchmark_FastqRenderMax(b *testing.B) {
	benchmarkFastqRender(runtime.NumCPU(), 151, 9999999, b)
}
