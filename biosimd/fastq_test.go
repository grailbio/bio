// Copyright 2018 GRAIL, Inc.  All rights reserved.
// Use of this source code is governed by the Apache-2.0
// license that can be found in the LICENSE file.

package biosimd_test

import (
	"fmt"
	"testing"

	"github.com/grailbio/base/simd"
	"github.com/grailbio/bio/biosimd"
	"github.com/grailbio/testutil/assert"
)

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

func TestFillFastqRecordBody(t *testing.T) {
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
		biosimd.FillFastqRecordBodyFromNibbles(dst, b16, &baseTable, &qualTable)
		assert.EQ(t, dst[:readLen], []byte(tc.baseAscii))
		assert.EQ(t, dst[readLen+3:2*readLen+3], []byte(tc.qualAscii))
	}
}
