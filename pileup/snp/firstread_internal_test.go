// Copyright 2020 Grail Inc.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
package snp

import (
	"math/rand"
	"strconv"
	"strings"
	"testing"

	"github.com/grailbio/bio/pileup"
	"github.com/grailbio/hts/sam"
)

func TestFirstreadSNPTable(t *testing.T) {
	nReadPair := 10000
	nSingleton := 500
	nFarPair := 5000
	nPos := PosType(2048)
	halfNPos := nPos / 2
	maxReadLen := 160
	maxReadSpan := 255
	maxReadSpanMinus1 := maxReadSpan - 1
	nCirc := PosType(512)
	nIter := 5

	// Must mock a non-nil sam.Reference, so that singletons can be represented
	// by MateRef == nil.
	ref, _ := sam.NewReference("ref1", "", "", int(nPos), nil, nil)

	for iter := 0; iter < nIter; iter++ {
		samrs := make([][]*sam.Record, nPos)
		// Generate fake read pairs separated by less than maxReadSpan.
		for i := 0; i < nReadPair; i++ {
			pos1 := PosType(rand.Intn(int(nPos)))
			pos2 := pos1 + PosType(rand.Intn(2*maxReadSpanMinus1+1)-maxReadSpanMinus1)
			if pos2 < 0 {
				pos2 = 0
			} else if pos2 >= nPos {
				pos2 = nPos - 1
			}
			readname := []byte("pair")
			readname = strconv.AppendInt(readname, int64(i), 10)

			// Best to mock CIGARs, due to the span computation.
			// We don't need to mock seq/qual to match the cigar here.
			samrs[pos1] = append(samrs[pos1], &sam.Record{
				Name:    string(readname),
				Ref:     ref,
				Pos:     int(pos1),
				Cigar:   []sam.CigarOp{sam.NewCigarOp(sam.CigarMatch, maxReadLen), sam.NewCigarOp(sam.CigarDeletion, maxReadSpan-maxReadLen)},
				Flags:   0x61,
				MateRef: ref,
				MatePos: int(pos2),
			})
			samrs[pos2] = append(samrs[pos2], &sam.Record{
				Name:    string(readname),
				Ref:     ref,
				Pos:     int(pos2),
				Cigar:   []sam.CigarOp{sam.NewCigarOp(sam.CigarMatch, maxReadLen), sam.NewCigarOp(sam.CigarDeletion, maxReadSpan-maxReadLen)},
				Flags:   0x91,
				MateRef: ref,
				MatePos: int(pos1),
			})
		}
		// Also add some singletons, and some read pairs separated by large spans.
		for i := 0; i < nSingleton; i++ {
			pos := PosType(rand.Intn(int(nPos)))
			readname := []byte("singleton")
			readname = strconv.AppendInt(readname, int64(i), 10)
			samrs[pos] = append(samrs[pos], &sam.Record{
				Name:    string(readname),
				Ref:     ref,
				Pos:     int(pos),
				Cigar:   []sam.CigarOp{sam.NewCigarOp(sam.CigarMatch, maxReadLen), sam.NewCigarOp(sam.CigarDeletion, maxReadSpan-maxReadLen)},
				Flags:   0x60,
				MatePos: -1,
			})
		}
		for i := 0; i < nFarPair; i++ {
			pos1 := PosType(rand.Intn(int(nPos)))
			mateGap := PosType(maxReadSpan + rand.Intn(int(halfNPos)-maxReadSpan))
			var pos2 PosType
			if pos1 >= halfNPos {
				pos2 = pos1 - mateGap
			} else {
				pos2 = pos1 + mateGap
			}
			readname := []byte("farpair")
			readname = strconv.AppendInt(readname, int64(i), 10)

			samrs[pos1] = append(samrs[pos1], &sam.Record{
				Name:    string(readname),
				Ref:     ref,
				Pos:     int(pos1),
				Cigar:   []sam.CigarOp{sam.NewCigarOp(sam.CigarMatch, maxReadLen), sam.NewCigarOp(sam.CigarDeletion, maxReadSpan-maxReadLen)},
				Flags:   0x61,
				MateRef: ref,
				MatePos: int(pos2),
			})
			samrs[pos2] = append(samrs[pos2], &sam.Record{
				Name:    string(readname),
				Ref:     ref,
				Pos:     int(pos2),
				Cigar:   []sam.CigarOp{sam.NewCigarOp(sam.CigarMatch, maxReadLen), sam.NewCigarOp(sam.CigarDeletion, maxReadSpan-maxReadLen)},
				Flags:   0x91,
				MateRef: ref,
				MatePos: int(pos1),
			})
		}

		firstReads := newFirstreadSNPTable(nCirc)

		nFoundRead := 0
		var readPair [2]readSNP
		readPair[0].seq8 = make([]byte, 0, maxReadLen)
		readPair[1].seq8 = make([]byte, 0, maxReadLen)
		for pos := PosType(0); pos < nPos; pos++ {
			for _, samr := range samrs[pos] {
				span, _ := samr.Cigar.Lengths()
				readPair[0].mapEnd = PosType(samr.Pos + span)
				nRead := firstReads.addOrRemove(&readPair, samr, pileup.StrandFwd, maxReadSpan)
				nFoundRead += nRead
				if nRead == 1 {
					// Name must start with 'singleton' or 'farpair'.
					if (!strings.HasPrefix(samr.Name, "singleton")) && (!strings.HasPrefix(samr.Name, "farpair")) {
						t.Fatal("firstreadSNPAddOrRemove lookup failed when it shouldn't have.")
					}
				} else if nRead == 2 {
					if !strings.HasPrefix(samr.Name, "pair") {
						t.Fatal("firstreadSNPAddOrRemove lookup succeeded when it shouldn't have.")
					}
					// Sanity check: read names should match, Pos/MatePos should be
					// mirror-images.
					if readPair[0].samr.Name != readPair[1].samr.Name {
						t.Fatal("mismatching read names from firstreadSNPAddOrRemove")
					}
					if (readPair[0].samr.Pos != readPair[1].samr.MatePos) ||
						(readPair[1].samr.Pos != readPair[0].samr.MatePos) {
						t.Fatal("mismatching positions from firstreadSNPAddOrRemove")
					}
				}
			}
		}
		if nFoundRead != ((nReadPair+nFarPair)*2 + nSingleton) {
			t.Fatal("TestFirstreadSNPTable: nFoundRead != (nReadPair + nFarPair) * 2 + nSingleton.")
		}
	}
}
