package snp

import (
	"math/rand"
	"reflect"
	"strconv"
	"testing"

	"github.com/grailbio/bio/interval"
	"github.com/grailbio/bio/pileup"
	"github.com/grailbio/hts/sam"
	"github.com/grailbio/testutil/assert"
)

func TestAlignRelevantBases(t *testing.T) {
	ref1, _ := sam.NewReference("chr1", "", "", 249250621, nil, nil)
	samHeader, _ := sam.NewHeader(nil, []*sam.Reference{ref1})
	opts := interval.NewBEDOpts{
		SAMHeader: samHeader,
	}
	bedPart, err := interval.NewBEDUnionFromEntries([]interval.Entry{
		{
			RefName: "chr1",
			Start0:  2488104,
			End:     2488172,
		},
		{
			RefName: "chr1",
			Start0:  2489165,
			End:     2489273,
		},
		{
			RefName: "chr1",
			Start0:  2489782,
			End:     2489907,
		},
		{
			RefName: "chr1",
			Start0:  2490320,
			End:     2490438,
		},
		{
			RefName: "chr1",
			Start0:  2491262,
			End:     2491417,
		},
		{
			RefName: "chr1",
			Start0:  2492063,
			End:     2492157,
		},
		{
			RefName: "chr1",
			Start0:  2493112,
			End:     2493254,
		},
		{
			RefName: "chr1",
			Start0:  2494304,
			End:     2494335,
		},
		{
			RefName: "chr1",
			Start0:  2494587,
			End:     2494712,
		},
	}, opts)
	assert.NoError(t, err)
	tests := []struct {
		pos   int
		cigar sam.Cigar
		want  []alignedPos
	}{
		{
			pos: 2488103,
			cigar: []sam.CigarOp{
				sam.NewCigarOp(sam.CigarMatch, 3),
			},
			want: []alignedPos{
				alignedPos{
					posInRef:  2488104,
					posInRead: 1,
				},
				alignedPos{
					posInRef:  2488105,
					posInRead: 2,
				},
			},
		},
		{
			pos: 2490400,
			cigar: []sam.CigarOp{
				sam.NewCigarOp(sam.CigarHardClipped, 1),
				sam.NewCigarOp(sam.CigarSoftClipped, 2),
				sam.NewCigarOp(sam.CigarMatch, 3),
				sam.NewCigarOp(sam.CigarSoftClipped, 4),
				sam.NewCigarOp(sam.CigarHardClipped, 5),
			},
			want: []alignedPos{
				alignedPos{
					posInRef:  2490400,
					posInRead: 2,
				},
				alignedPos{
					posInRef:  2490401,
					posInRead: 3,
				},
				alignedPos{
					posInRef:  2490402,
					posInRead: 4,
				},
			},
		},
		{
			pos: 2490420,
			cigar: []sam.CigarOp{
				sam.NewCigarOp(sam.CigarMatch, 1),
				sam.NewCigarOp(sam.CigarDeletion, 15),
				sam.NewCigarOp(sam.CigarMatch, 1),
				sam.NewCigarOp(sam.CigarInsertion, 15),
				sam.NewCigarOp(sam.CigarMatch, 826),
			},
			want: []alignedPos{
				alignedPos{
					posInRef:  2490420,
					posInRead: 0,
				},
				alignedPos{
					posInRef:  2490436,
					posInRead: 1,
				},
				alignedPos{
					posInRef:  2490437,
					posInRead: 17,
				},
				alignedPos{
					posInRef:  2491262,
					posInRead: 842,
				},
			},
		},
	}
	var result []alignedPos
	for _, tt := range tests {
		var mockRead readSNP
		// Not currently necessary to mock mockRead.seq8.
		mockRead.samr = &sam.Record{
			Ref:   ref1,
			Pos:   tt.pos,
			Cigar: tt.cigar,
		}
		span, _ := tt.cigar.Lengths()
		mockRead.mapEnd = PosType(tt.pos + span)
		err = alignRelevantBases(&result, mockRead, &bedPart)
		assert.NoError(t, err)
		if !reflect.DeepEqual(result, tt.want) {
			t.Errorf("Wanted: %v  Got: %v", tt.want, result)
		}
	}
}

// This only verifies that AddOrphanReads() removes the expected number of
// orphan reads from the firstread-table; it does not check the accuracy of the
// resulting pileup.  See TestPileupSNP() below for the latter.
func TestAddOrphanReads(t *testing.T) {
	rand.Seed(1)
	nReadPair := 20000
	nPos := PosType(4096)
	halfNPos := nPos / 2
	orphanStopPos := halfNPos - 50
	maxReadLen := 500
	maxReadSpan := 511
	maxReadSpanMinus1 := maxReadSpan - 1
	nCirc := PosType(1024)
	nIter := 5

	ref, _ := sam.NewReference("chr1", "", "", int(nPos), nil, nil)
	samHeader, _ := sam.NewHeader(nil, []*sam.Reference{ref})
	opts := interval.NewBEDOpts{
		SAMHeader: samHeader,
	}
	minBaseQual := byte(43)
	qpt, _ := newQualPassTable(minBaseQual)
	pCtx := pileupContext{
		qpt:    &qpt,
		stitch: true,
	}
	var err error
	pCtx.bedPart, err = interval.NewBEDUnionFromEntries([]interval.Entry{
		{
			RefName: "chr1",
			Start0:  0,
			End:     2048,
		},
	}, opts)
	assert.NoError(t, err)

	dummySeq := make([]byte, maxReadLen)
	for i := range dummySeq {
		dummySeq[i] = byte(sam.BaseN)
	}
	dummyQual := make([]byte, maxReadLen)
	for i := range dummyQual {
		dummyQual[i] = minQual
	}
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

			// Now necessary to mock CIGAR/Seq/Qual.  (These are degenerate since we
			// are only testing orphan-read detection here, not CIGAR-iteration.)
			samrs[pos1] = append(samrs[pos1], &sam.Record{
				Name:    string(readname),
				Ref:     ref,
				Pos:     int(pos1),
				Cigar:   []sam.CigarOp{sam.NewCigarOp(sam.CigarMatch, maxReadLen), sam.NewCigarOp(sam.CigarDeletion, maxReadSpan-maxReadLen)},
				Flags:   0x61,
				MateRef: ref,
				MatePos: int(pos2),
				Seq:     sam.NewSeq(dummySeq),
				Qual:    append([]byte{}, dummyQual...),
			})
			samrs[pos2] = append(samrs[pos2], &sam.Record{
				Name:    string(readname),
				Ref:     ref,
				Pos:     int(pos2),
				Cigar:   []sam.CigarOp{sam.NewCigarOp(sam.CigarMatch, maxReadLen), sam.NewCigarOp(sam.CigarDeletion, maxReadSpan-maxReadLen)},
				Flags:   0x91,
				MateRef: ref,
				MatePos: int(pos1),
				Seq:     sam.NewSeq(dummySeq),
				Qual:    append([]byte{}, dummyQual...),
			})
		}

		results := newPileupMutable(nCirc, maxReadLen, true, nil)

		nFoundRead := 0
		nOrphanRead := 0
		var readPair [2]readSNP
		readPair[0].seq8 = make([]byte, 0, maxReadLen)
		readPair[1].seq8 = make([]byte, 0, maxReadLen)
		for pos := PosType(0); pos < nPos; pos++ {
			if pos == halfNPos {
				if err = results.addOrphanReads(&pCtx, orphanStopPos); err != nil {
					t.Fatalf("pileupMutable.addOrphanReads error: %v", err)
				}
			}
			for _, samr := range samrs[pos] {
				span, _ := samr.Cigar.Lengths()
				readPair[0].mapEnd = PosType(samr.Pos + span)
				nRead := results.firstReads.addOrRemove(&readPair, samr, pileup.StrandFwd, maxReadSpan)
				nFoundRead += nRead
				if nRead == 1 {
					// Lookup should fail iff current position >= halfNPos and mate
					// position < orphanStopPos.
					if samr.Pos < int(halfNPos) {
						t.Fatal("firstreadSNPTable.addOrRemove lookup failed when it shouldn't have.")
					} else if samr.MatePos >= int(orphanStopPos) {
						t.Fatal("pileupMutable.addOrphanReads processed a read it shouldn't have.")
					}
					// addOrphanReads doesn't directly return the number of orphan reads
					// it handled, but with our test setup, every nRead==1 return value
					// from addOrRemove must correspond to a handled orphan read (and the
					// converse is also true).
					nOrphanRead++
				} else if nRead == 2 {
					if (samr.Pos >= int(halfNPos)) && (samr.MatePos < int(orphanStopPos)) {
						t.Fatal("pileupMutable.addOrphanReads did not process a read it should have.")
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
		if nFoundRead+nOrphanRead != (nReadPair * 2) {
			t.Fatal("TestAddOrphanReads: nFoundRead + nOrphanRead != nReadPair * 2.")
		}
	}
}
