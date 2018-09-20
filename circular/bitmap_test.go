// Copyright 2018 GRAIL, Inc.  All rights reserved.
// Use of this source code is governed by the Apache-2.0
// license that can be found in the LICENSE file.

package circular_test

import (
	"math/rand"
	"testing"

	"github.com/grailbio/base/bitset"
	"github.com/grailbio/bio/circular"
	bi "github.com/grailbio/bio/interval"
)

func setSomeBits(cb *circular.Bitmap, bitCounts []int, nRowBit, start, end, n int) {
	diff := end - start
	mask := int(cb.NCirc()) - 1
	for i := 0; i < n; {
		col := rand.Intn(nRowBit)
		pos := rand.Intn(diff) + start
		circPos := pos & mask
		if !bitset.Test(cb.Row(bi.PosType(circPos)), col) {
			cb.Set(bi.PosType(pos), bi.PosType(circPos), uint32(col))
			bitCounts[pos]++
			i++
		}
	}
}

func TestBitmap(t *testing.T) {
	maxSize := 10000
	nIter := 200
	for iter := 0; iter < nIter; iter++ {
		// guarantee size >= 4 so test plan makes sense
		size := circular.NextExp2(rand.Intn(maxSize) + 2)
		rowWidth := rand.Intn(255) + 1
		cb := circular.NewBitmap(bi.PosType(size), bi.PosType(rowWidth))
		nRowBit := rowWidth * bitset.BitsPerWord
		bitCounts := make([]int, 3*size/2)

		// guaranteed to be low enough to avoid saturation
		bitsToSet := rand.Intn(size) * rand.Intn(bitset.BitsPerWord/8)

		// 1. Set some bits in [0, size/2)
		// 2. Iterate up to size/4
		// 3. Set some more bits in [size/4, size)
		// 4. Iterate up to 3*size / 4
		// 5. Set some more bits in [3*size/4, 3*size/2)
		// 6. Iterate up to 5*size / 4; this should test wraparound
		setSomeBits(&cb, bitCounts, nRowBit, 0, size/2, bitsToSet)
		for pos := cb.FirstPos(); pos < bi.PosType(size/4); pos = cb.FirstPos() {
			i := 0
			for s, col := cb.NewRowScanner(); col != -1; col = s.Next() {
				i++
			}
			if i != bitCounts[pos] {
				t.Fatal("Mismatched set-bit counts (part 1).")
			}
		}
		setSomeBits(&cb, bitCounts, nRowBit, size/4, size, bitsToSet)
		for pos := cb.FirstPos(); pos < bi.PosType(3*size/4); pos = cb.FirstPos() {
			i := 0
			for s, col := cb.NewRowScanner(); col != -1; col = s.Next() {
				i++
			}
			if i != bitCounts[pos] {
				t.Fatal("Mismatched set-bit counts (part 2).")
			}
		}
		setSomeBits(&cb, bitCounts, nRowBit, 3*size/4, 3*size/2, bitsToSet)
		for pos := cb.FirstPos(); pos < bi.PosType(5*size/4); pos = cb.FirstPos() {
			i := 0
			for s, col := cb.NewRowScanner(); col != -1; col = s.Next() {
				i++
			}
			if i != bitCounts[pos] {
				t.Fatal("Mismatched set-bit counts (part 3).")
			}
		}
	}
}
