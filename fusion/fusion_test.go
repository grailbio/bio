package fusion

import (
	"fmt"
	"io/ioutil"
	"strings"
	"testing"

	"github.com/grailbio/base/vcontext"
	"github.com/grailbio/testutil"
	"github.com/grailbio/testutil/assert"
	"github.com/grailbio/testutil/expect"
	"github.com/grailbio/testutil/h"
)

func TestPrefixSuffixLength(t *testing.T) {
	pos := []PosRange{{1, 10}, {12, 21}, {22, 25}, {27, 36}}
	expect.EQ(t, prefixLength(pos), []int{9, 18, 21, 30})
	expect.EQ(t, suffixLength(pos), []int{30, 21, 12, 9})

	pos = []PosRange{{0, 10}, {5, 6}, {10, 15}}
	expect.EQ(t, prefixLength(pos), []int{10, 11, 16})
	expect.EQ(t, suffixLength(pos), []int{16, 6, 5})
}

func testListKmers(seq string, opts Opts) []kmersAtPos {
	k := newKmerizer(opts.KmerLength)
	k.Reset(seq)
	var pos []kmersAtPos
	for k.Scan() {
		pos = append(pos, k.Get())
	}
	return pos
}

func TestKmerizer(t *testing.T) {
	opts := Opts{KmerLength: 5}
	expect.That(t, testListKmers("AAAGTTCAGGT", opts),
		h.ElementsAre(
			kmersAtPos{0, 11, 127},
			kmersAtPos{1, 47, 31},
			kmersAtPos{2, 189, 519},
			kmersAtPos{3, 756, 897},
			kmersAtPos{4, 978, 480},
			kmersAtPos{5, 842, 376},
			kmersAtPos{6, 299, 94},
		))
}

func testNewFragment(name, seq string, opts Opts) Fragment {
	frag := Fragment{
		Name:  name,
		R1Seq: seq,
		kmers: testListKmers(seq, opts),
	}
	return frag
}

func testWriteFile(dir, data string) string {
	f, err := ioutil.TempFile(dir, "")
	if err != nil {
		panic(err)
	}
	if _, err := f.Write([]byte(data)); err != nil {
		panic(err)
	}
	if err := f.Close(); err != nil {
		panic(err)
	}
	return f.Name()
}

func TestInferGeneRangeInfo(t *testing.T) {
	ctx := vcontext.Background()
	tempDir, cleanup := testutil.TempDir(t, "", "")
	defer cleanup()

	const kmerLength = 5
	opts := DefaultOpts
	opts.KmerLength = kmerLength
	frag := testNewFragment("f1", "AAAGTTCAGGT", opts)
	geneDB := NewGeneDB(opts)

	transcriptomePath := testWriteFile(tempDir, `>E1|YWHAE|chr1:1-2:3|first kmer in f1
AAAGT
>E2|YWHAE|chr1:1-2:3|reverse-complement of AAGTT
AACTT
>E2|FAM22A|chr1:2-3:4|reverse-complement of TCAGG
CCTGA
`)
	geneDB.ReadTranscriptome(ctx, transcriptomePath, false /*denovo*/)
	expect.EQ(t, inferGeneRangeInfo(frag, geneDB, opts.KmerLength),
		[]geneRangeInfo{
			geneRangeInfo{geneID: geneDB.geneID("YWHAE"), r1Span: 6, ranges: []PosRange{{0, 6}}},
			geneRangeInfo{geneID: geneDB.geneID("FAM22A"), r1Span: 5, ranges: []PosRange{{5, 10}}}})
}

// newR2PosRange creates a new PosRange for a range in R2.
//
// REQUIRES: start <= r2PosOffset, end <= r2PosOffset, start <= end
func newR2PosRange(start, end Pos) PosRange {
	if end < start {
		panic("inverted range")
	}
	return PosRange{newR2Pos(start), newR2Pos(end)}
}

func TestInferCandidatePair(t *testing.T) {
	ctx := vcontext.Background()
	tempDir, cleanup := testutil.TempDir(t, "", "")
	defer cleanup()

	cosmicPath := testWriteFile(tempDir,
		"Genes\tSamples\tMutations\tPapers\n"+
			"NUTM2B/YWHAE\t736\t13\t8\n"+
			"YWHAE/FAM22A\t685\t2\t8\n")

	opts := Opts{
		KmerLength:     13,
		MaxGap:         6,
		UnstrandedPrep: false,
		MaxHomology:    13,
		MinSpan:        13,
	}
	geneDB := NewGeneDB(opts)
	geneDB.ReadFusionEvents(ctx, cosmicPath)

	// Prepare GeneRangeInfo Vector for fusion YWHAE/NUTM2B
	// Assuming:
	//  GeneName  r1_span, r2_span [r1_ranges] | [r2_ranges]
	//   YWHAE    142,      98,    [0, 141] | [160, 201] [203, 258]
	//   NUTM2B   0  ,      30,             | [250, 283]
	gri := []geneRangeInfo{
		geneRangeInfo{geneID: geneDB.geneID("YWHAE"), r1Span: 142, r2Span: 98,
			ranges: []PosRange{{0, 142}, newR2PosRange(160, 202), newR2PosRange(203, 259)}},
		geneRangeInfo{geneID: geneDB.geneID("NUTM2B"), r1Span: 0, r2Span: 40,
			ranges: []PosRange{newR2PosRange(250, 284)}},
	}

	got := inferCandidatePair(geneDB, "xx", 145, 145, gri, opts)
	expect.EQ(t, len(got), 1)
	expect.EQ(t, got[0].Name(geneDB, opts), "NUTM2B/YWHAE/-")
	expect.EQ(t, got[0].G1Span, 240)
	expect.EQ(t, got[0].G2Span, 25)
	expect.EQ(t, got[0].JointSpan, 265)
	expect.EQ(t, got[0].FusionOrder, true)
	expect.EQ(t, got[0].G1Range, newCrossReadPosRange(0, newR2Pos(259)))
	expect.EQ(t, got[0].G2Range, newCrossReadPosRange(newR2Pos(259), newR2Pos(284)))
}

func TestInferLongestCombinedSpan(t *testing.T) {
	totalSpan := func(ranges []PosRange, typ ReadType) int {
		span := 0
		for _, r := range ranges {
			if r.readType() == typ {
				span += posSpan(r.End, r.Start)
			}
		}
		return span
	}
	runTest := func(rangesL, rangesR []PosRange, maxGap int) FusionInfo {
		gL := geneRangeInfo{
			geneID: 100,
			r1Span: totalSpan(rangesL, R1),
			r2Span: totalSpan(rangesL, R2),
			ranges: rangesL,
		}
		gR := geneRangeInfo{
			geneID: 101,
			r1Span: totalSpan(rangesR, R1),
			r2Span: totalSpan(rangesR, R2),
			ranges: rangesR,
		}
		opts := DefaultOpts
		opts.MaxGap = maxGap
		opts.MaxHomology = 20
		fi := inferLongestCombinedSpan(gL, gR, opts)
		// Clear fields that have constant values, so that individual assertions
		// don't need to fill the expectations.
		assert.EQ(t, fi.G1ID, GeneID(100))
		assert.EQ(t, fi.G2ID, GeneID(101))
		assert.EQ(t, fi.RefOrder, true)
		fi.G1ID = 0
		fi.G2ID = 0
		fi.RefOrder = false
		return fi
	}
	// LHS    ^ --------------  ---------------
	// RHS    ^---------  --------- ---
	// Taking first range of rhs and entire lhs.
	expect.EQ(t, runTest(
		[]PosRange{{2, 16}, {18, 33}},
		[]PosRange{{1, 10}, {12, 21}, {22, 25}},
		0),
		FusionInfo{
			G1Span:      29,
			G2Span:      1,
			JointSpan:   30,
			FusionOrder: false,
			G1Range:     newCrossReadPosRange(2, 33),
			G2Range:     newCrossReadPosRange(1, 2)})
	// LHS    ^ --------------  ---------------
	// RHS    ^---------  ---------              ---------
	// Taking lhs + last range of rhs.
	expect.EQ(t, runTest(
		[]PosRange{{2, 16}, {18, 33}},
		[]PosRange{{1, 10}, {12, 21}, newR2PosRange(34, 36)},
		0),
		FusionInfo{
			G1Span:      29,
			G2Span:      2,
			JointSpan:   31,
			FusionOrder: true,
			G1Range:     newCrossReadPosRange(2, 33),
			G2Range:     newCrossReadPosRange(newR2Pos(34), newR2Pos(36))})
	// LHS    ^ --------------
	// RHS    ^---------  ---------      ---------
	// Taking first range of lhs and second + ranges of rhs.
	expect.EQ(t, runTest(
		[]PosRange{{2, 16}},
		[]PosRange{{1, 10}, {12, 21}, newR2PosRange(34, 36)},
		0),
		FusionInfo{
			G1Span:      14,
			G2Span:      7,
			JointSpan:   21,
			FusionOrder: true,
			G1Range:     newCrossReadPosRange(2, 16),
			G2Range:     newCrossReadPosRange(16, newR2Pos(36))})
	// Only one range is non-empty
	// LHS    ^ --------------
	// RHS    ^
	expect.EQ(t, runTest(
		[]PosRange{{2, 16}},
		nil,
		0),
		FusionInfo{
			G1Span:      14,
			G2Span:      0,
			JointSpan:   14,
			FusionOrder: true,
			G1Range:     newCrossReadPosRange(2, 16),
			G2Range:     newCrossReadPosRange(0, 0)})
	// LHS completely overlaps RHS with overlap < max_homology
	// Only one range is non-empty
	// LHS    ^-----------------      ---------------------------
	// RHS    ^           -----                           ------
	expect.EQ(t, runTest(
		[]PosRange{{0, 129}, {147, 273}},
		[]PosRange{{114, 127}, {265, 271}},
		0),
		FusionInfo{
			G1Span:      255,
			G2Span:      0,
			JointSpan:   255,
			FusionOrder: true,
			G1Range:     newCrossReadPosRange(0, 273),
			G2Range:     newCrossReadPosRange(0, 0)})
	// LHS completely overlaps RHS with overlap > max_homology (changes nothing
	// since max span is still LHS)
	expect.EQ(t, runTest(
		[]PosRange{{0, 129}, {147, 273}},
		[]PosRange{{114, 127}, {258, 271}},
		0),
		FusionInfo{
			G1Span:      255,
			G2Span:      0,
			JointSpan:   255,
			FusionOrder: true,
			G1Range:     newCrossReadPosRange(0, 273),
			G2Range:     newCrossReadPosRange(0, 0)})
	// Taking into consideration of junction_pos and max_gap
	// case(1) span1 and span2 are far away but located on different side of
	// junction_pos
	// --------------|------------
	// -------- span1
	//                   ------- span2
	expect.EQ(t, runTest(
		[]PosRange{{0, 129}},
		[]PosRange{newR2PosRange(160, 251)},
		10),
		FusionInfo{
			G1Span:      129,
			G2Span:      91,
			JointSpan:   220,
			FusionOrder: true,
			G1Range:     newCrossReadPosRange(0, 129),
			G2Range:     newCrossReadPosRange(newR2Pos(160), newR2Pos(251))})
	// case(2) span1 and span2 are far away and located on same side of
	// junction_pos, max_span = first range of g1 + first range of g2 as it is
	// larger than span of g1
	// --------------|------------
	// -------------  ---
	//                        -----
	expect.EQ(t, runTest(
		[]PosRange{{0, 129}, newR2PosRange(142, 151)},
		[]PosRange{newR2PosRange(170, 251)},
		10),
		FusionInfo{
			G1Span:      129,
			G2Span:      81,
			JointSpan:   210,
			FusionOrder: true,
			G1Range:     newCrossReadPosRange(0, 129),
			G2Range:     newCrossReadPosRange(newR2Pos(170), newR2Pos(251))})
	// case(3) gene2 is a subset of gene1
	// --------------|------------
	// -------------  ---
	//     -----
	expect.EQ(t, runTest(
		[]PosRange{{0, 129}, newR2PosRange(142, 151)},
		[]PosRange{{50, 90}},
		10),
		FusionInfo{
			G1Span:      138,
			G2Span:      0,
			JointSpan:   138,
			FusionOrder: true,
			G1Range:     newCrossReadPosRange(0, newR2Pos(151)),
			G2Range:     newCrossReadPosRange(0, 0)})
	// case(4) gene2 overlaps with gene1 < homology
	// --------------|------------
	// -------------  ---
	//                 ------
	expect.EQ(t, runTest(
		[]PosRange{{0, 129}, newR2PosRange(142, 151)},
		[]PosRange{{newR2Pos(145), newR2Pos(161)}},
		10),
		FusionInfo{
			G1Span:      138,
			G2Span:      10,
			JointSpan:   148,
			FusionOrder: true,
			G1Range:     newCrossReadPosRange(0, newR2Pos(151)),
			G2Range:     newCrossReadPosRange(newR2Pos(151), newR2Pos(161))})

	// case(4) gene2 overlaps with gene1 > homology
	// --------------|------------
	// -------------  ---
	//                 ------
	expect.EQ(t, runTest(
		[]PosRange{{0, 129}, newR2PosRange(142, 171)},
		[]PosRange{{newR2Pos(145), newR2Pos(181)}},
		10),
		FusionInfo{FusionOrder: true})
}

func testInternGene(geneDB *GeneDB, name, ref string, start, end, index int) GeneID {
	geneID := geneDB.internGene(name)
	geneDB.genes[geneID] = &GeneInfo{
		ID:    geneID,
		Gene:  name,
		Chrom: ref,
		Start: start,
		End:   end,
		Index: index,
	}
	return geneID
}

func TestLinkedByLowComplexSubstring(t *testing.T) {
	geneDB := NewGeneDB(DefaultOpts)
	expect.True(t,
		LinkedByLowComplexSubstring(Fragment{
			R1Seq: "aaaaaaaaaa",
			R2Seq: "aaaaaaaaaa"},
			FusionInfo{
				G1ID:    testInternGene(geneDB, "a0", "chr1", 0, 1, 0),
				G2ID:    testInternGene(geneDB, "a1", "chr1", 0, 1, 0),
				G1Span:  10,
				G2Span:  10,
				G1Range: newCrossReadPosRange(0, 9),
				G2Range: newCrossReadPosRange(newR2Pos(1), newR2Pos(9)),
			}, 0.9))

	expect.False(t,
		LinkedByLowComplexSubstring(Fragment{
			R1Seq: "acgtacgt",
			R2Seq: "tcgatcga"},
			FusionInfo{
				G1ID:    testInternGene(geneDB, "b0", "chr1", 0, 1, 0),
				G2ID:    testInternGene(geneDB, "b1", "chr1", 0, 1, 0),
				G1Span:  0,
				G2Span:  7,
				G1Range: newCrossReadPosRange(0, 7),
				G2Range: newCrossReadPosRange(newR2Pos(2), newR2Pos(8)),
			}, 0.9))
}

func TestFilterByMinSpan(t *testing.T) {
	geneDB := NewGeneDB(DefaultOpts)
	cid := 0
	newCandidate := func(g1Name, g2Name, umi string, g1Span, g2Span, g1Start, g1End, g2Start, g2End int) Candidate {
		id := cid
		cid++
		return Candidate{
			Frag: Fragment{
				Name: fmt.Sprintf("x:%s seq%d", umi, id),
			},
			Fusions: []FusionInfo{FusionInfo{
				G1ID:    testInternGene(geneDB, g1Name, "chr1", 0, 1, 0),
				G2ID:    testInternGene(geneDB, g2Name, "chr1", 0, 1, 0),
				G1Span:  g1Span,
				G2Span:  g2Span,
				G1Range: newCrossReadPosRange(Pos(g1Start), Pos(g2End)),
				G2Range: newCrossReadPosRange(newR2Pos(Pos(g2Start)), newR2Pos(Pos(g2End))),
			}},
		}
	}

	// only b will be retained; a is out because lack of collapsed support;
	// c is out because lack of number of support on both span.
	orgCandidates := []Candidate{
		newCandidate("a", "a", "AA+CC", 30, 25, 0, 29, 0, 24),
		newCandidate("a", "a", "AA+CC", 30, 25, 0, 29, 0, 24),
		newCandidate("b", "b", "AA+GG", 30, 25, 0, 29, 0, 24),
		newCandidate("b", "b", "AA+GG", 30, 25, 0, 29, 0, 24),
		newCandidate("b", "b", "AA+CC", 30, 30, 0, 29, 0, 24),
		newCandidate("b", "b", "AA+CC", 30, 30, 0, 29, 0, 24),
		newCandidate("b", "b", "AA+TC", 30, 12, 0, 29, 0, 24),
		newCandidate("c", "c", "AA+TT", 30, 25, 0, 29, 0, 24),
		newCandidate("c", "c", "AA+CC", 30, 13, 0, 29, 0, 24),
	}
	expected := []Candidate{orgCandidates[2], orgCandidates[4], orgCandidates[6]}
	candidates := orgCandidates
	FilterByMinSpan(true, 25, &candidates, 2)
	expect.EQ(t, candidates, expected)

	orgCandidates = []Candidate{
		newCandidate("a", "a", "", 30, 25, 0, 29, 0, 24),
		newCandidate("a", "a", "", 30, 25, 0, 29, 0, 24),
		newCandidate("b", "b", "", 30, 25, 0, 29, 0, 24),
		newCandidate("b", "b", "", 30, 25, 0, 29, 0, 24),
		newCandidate("b", "b", "", 30, 30, 0, 29, 0, 24),
		newCandidate("b", "b", "", 30, 30, 0, 29, 0, 24),
		newCandidate("b", "b", "", 30, 12, 0, 29, 0, 24),
		newCandidate("c", "c", "", 30, 25, 0, 29, 0, 24),
		newCandidate("c", "c", "", 30, 13, 0, 29, 0, 24),
	}
	expected = []Candidate{
		orgCandidates[0],
		orgCandidates[1],
		orgCandidates[2],
		orgCandidates[3],
		orgCandidates[4],
		orgCandidates[5],
		orgCandidates[6],
	}
	candidates = orgCandidates
	FilterByMinSpan(false, 25, &candidates, 2)
	expect.EQ(t, candidates, expected)
}

func TestCloseProximity(t *testing.T) {
	geneDB := NewGeneDB(DefaultOpts)

	type genePair struct {
		name1                string
		ref1                 string
		start1, end1, index1 int

		name2                string
		ref2                 string
		start2, end2, index2 int
	}
	runTest := func(gp genePair,
		//g1Name, g2Name string,
		//g1Ref string, g1RefStart, g1RefEnd, g1Index int,
		//g2Ref string, g2RefStart, g2RefEnd, g2Index int
		g1Span, g2Span, g1Start, g1End, g2Start, g2End int) bool {
		fi := FusionInfo{
			G1ID:    testInternGene(geneDB, gp.name1, gp.ref1, gp.start1, gp.end1, gp.index1),
			G2ID:    testInternGene(geneDB, gp.name2, gp.ref2, gp.start2, gp.end2, gp.index2),
			G1Span:  g1Span,
			G2Span:  g2Span,
			G1Range: newCrossReadPosRange(Pos(g1Start), Pos(g2End)),
			G2Range: newCrossReadPosRange(newR2Pos(Pos(g2Start)), newR2Pos(Pos(g2End))),
		}
		return CloseProximity(geneDB, fi, 500, 5)
	}

	// Fail distance < 500
	expect.True(t, runTest(
		genePair{"a0", "chr1", 2000, 2500, 4, "a1", "chr1", 2750, 3500, 15},
		30, 25, 0, 29, 0, 24))

	// Fail num < 5
	expect.True(t, runTest(
		genePair{"b0", "chr1", 2000, 2500, 4, "b1", "chr1", 3750, 4500, 7},
		30, 25, 0, 29, 0, 24))

	// Fail one of candidates < 50000 (the second one)
	expect.True(t, runTest(
		genePair{"c0", "chr1", 1000, 1500, 4, "c1", "chr1", 2000, 2500, 4},
		30, 25, 0, 29, 0, 24))

	expect.True(t, runTest(
		genePair{"d0", "chr1", 3500, 4500, 15, "d1", "chr1", 2750, 3500, 15},
		30, 25, 0, 29, 0, 24))

	// Fail overlap
	expect.True(t, runTest(
		genePair{"e0", "chr1", 2000, 3500, 4, "e1", "chr1", 2750, 3000, 10},
		30, 25, 0, 29, 0, 243))

	// pass same chrom
	expect.False(t, runTest(
		genePair{"f0", "chr1", 2000, 2500, 4, "f1", "chr1", 3750, 4500, 10},
		30, 25, 0, 29, 0, 24))

	// pass diff chrom
	expect.False(t, runTest(
		genePair{"a0", "chr1", 2000, 2500, 4, "a1", "chr2", 3750, 4500, 7},
		30, 25, 0, 29, 0, 24))
}

func TestFilterDuplicates(t *testing.T) {
	geneDB := NewGeneDB(DefaultOpts)
	cid := 0
	newCandidate := func(g1Name, g2Name, umi string, g1Span, g2Span, g1Start, g1End, g2Start, g2End int, seq string) Candidate {
		id := cid
		cid++
		return Candidate{
			Frag: Fragment{
				Name:  fmt.Sprintf("x:%s seq%d", umi, id),
				R1Seq: seq,
			},
			Fusions: []FusionInfo{FusionInfo{
				G1ID:    testInternGene(geneDB, g1Name, "chr1", 0, 1, 0),
				G2ID:    testInternGene(geneDB, g2Name, "chr1", 0, 1, 0),
				G1Span:  g1Span,
				G2Span:  g2Span,
				G1Range: newCrossReadPosRange(Pos(g1Start), Pos(g2End)),
				G2Range: newCrossReadPosRange(newR2Pos(Pos(g2Start)), newR2Pos(Pos(g2End))),
			}},
		}
	}

	candidates := []Candidate{
		// Same seq, same fusion, different UMIs -- Both unique
		// 0
		newCandidate("a", "a", "AATT+CCGG", 30, 25, 0, 29, 0, 24, "ACTGACTGACTGACTGACTG"),
		// 1
		newCandidate("a", "a", "AAGG+GGCC", 30, 25, 0, 29, 30, 54, "ACTGACTGACTGACTGACTG"),
		// 2 with same UMI (and same seq), others with diff seq -- 3 unique
		// 2
		newCandidate("b", "b", "AATT+CCGG", 30, 25, 0, 29, 0, 24, "ACTGACTGACTGACTGACTG"),
		// 3 -- remove because of UMI satisfy hd < 3 with 2
		newCandidate("b", "b", "AATT+CCTN", 30, 25, 0, 29, 0, 24, "ACTGACTGACTGACTGACTG"),
		// 4 -- keep
		newCandidate("b", "b", "AAGG+GGCC", 30, 25, 0, 29, 0, 24, "ACTGACTGACTGACTGA"),
		// 5 -- remove because of UMI satisfy hd < 3 with 2
		newCandidate("b", "b", "AATT+CCGG", 30, 30, 0, 29, 0, 24, "ACTGACTG"),
		// 6 -- remove because of UMI satisfy hd < 3 with 4
		newCandidate("b", "b", "AAGG+GGCC", 30, 30, 0, 29, 0, 24, "ACTGACTG"),
		// 7
		// Single seq
		newCandidate("c", "c", "AA+CC", 30, 13, 0, 29, 0, 24, "ACTGACTGACTGACTGACTG"),
		// 8
		newCandidate("d", "d", "GCGATTAA+GATCTGCT", 30, 13, 0, 29, 0, 24, "ACTGACTGACTGACTGACTG"),
		// 9 -- remove
		newCandidate("d", "d", "NCGATTAA+NATCTGCT", 30, 13, 0, 29, 0, 24, "ACTGACTGACTGACTGACTG"),
	}
	want := []Candidate{
		candidates[0],
		candidates[1],
		candidates[2],
		candidates[4],
		candidates[7],
		candidates[8],
	}
	FilterDuplicates(&candidates, true /*hasUMI*/)
	expect.EQ(t, candidates, want)

	// No UMI case
	candidates = []Candidate{
		// Same seq, same fusion -- 1 unique
		newCandidate("a", "a", "", 30, 25, 0, 29, 0, 24, "ACTGACTGACTGACTGACTG"),
		newCandidate("a", "a", "", 30, 25, 0, 29, 0, 24, "ACTGACTGACTGACTGACTG"),
		// 2 with same seq, others with diff seq -- 3 unique
		newCandidate("b", "b", "", 30, 25, 0, 29, 0, 24, "ACTGACTGACTGACTGACTG"),
		newCandidate("b", "b", "", 30, 25, 0, 29, 0, 24, "ACTGACTGACTGACTGACTG"),
		newCandidate("b", "b", "", 30, 25, 0, 29, 0, 24, "ACTGACTGACTGACTGA"),
		// the case that tolerates sequencing error
		newCandidate("b", "b", "", 30, 30, 0, 29, 0, 24, "NCTGACTGACTGACTGACTG"),
		newCandidate("b", "b", "", 30, 30, 0, 29, 0, 24, "ACTGACTG"),
		// Single seq
		newCandidate("c", "c", "", 30, 13, 0, 29, 0, 24, "ACTGACTGACTGACTGACTG"),
	}
	want = []Candidate{
		candidates[0],
		candidates[2],
		candidates[4],
		candidates[6],
		candidates[7],
	}
	FilterDuplicates(&candidates, false /*!hasUMI*/)
	expect.EQ(t, candidates, want)
}

func TestDiscardAbundantPartners(t *testing.T) {
	type genePair struct{ g1, g2 string }
	geneDB := NewGeneDB(DefaultOpts)
	newCandidate := func(pairs []genePair, g1Span, g2Span, g1Start, g1End, g2Start, g2End int, seqName, seq string) Candidate {
		c := Candidate{
			Frag: Fragment{
				Name:  seqName,
				R1Seq: seq,
			},
		}
		for _, gp := range pairs {
			c.Fusions = append(c.Fusions, FusionInfo{
				G1ID:    testInternGene(geneDB, gp.g1, "", 0, 1, 1),
				G2ID:    testInternGene(geneDB, gp.g2, "", 0, 1, 1),
				G1Span:  g1Span,
				G2Span:  g2Span,
				G1Range: newCrossReadPosRange(Pos(g1Start), Pos(g2End)),
				G2Range: newCrossReadPosRange(newR2Pos(Pos(g2Start)), newR2Pos(Pos(g2End))),
			})
		}
		return c
	}

	candidates := []Candidate{
		// a - 1 partner : retain
		newCandidate([]genePair{{"a", "a"}}, 30, 25, 0, 29, 30, 54, "seq1", "ACTGACTGACTGACTGACTG"),
		// b, c 2 partners : retain
		newCandidate([]genePair{{"b", "b"}}, 30, 25, 0, 29, 30, 54, "seq3", "ACTGACTGACTGACTGACTG"),
		newCandidate([]genePair{{"b", "c"}}, 30, 30, 0, 29, 30, 54, "seq6", "ACTGACTG"),
		newCandidate([]genePair{{"c", "c"}}, 30, 13, 0, 29, 30, 54, "seq7", "ACTGACTGACTGACTGACTG"),
		// d,e,f,g - e has 3 partners : reject
		newCandidate([]genePair{{"d", "e"}}, 30, 25, 0, 29, 30, 54, "seq3", "ACTGACTGACTGACTGACTG"),
		newCandidate([]genePair{{"e", "f"}}, 30, 25, 0, 29, 30, 54, "seq4", "ACTGACTGACTGACTGACTG"),
		newCandidate([]genePair{{"e", "g"}}, 30, 25, 0, 29, 30, 54, "seq4", "ACTGACTGACTGACTGACTG"),
		// Muti with 2 occurrences : retain
		newCandidate([]genePair{{"h", "i"}, {"j", "k"}}, 30, 25, 0, 29, 30, 54, "seq4", "ACTGACTGACTGACTGACTG"),
		// Muti with 3 occurrences : reject
		newCandidate([]genePair{{"l", "m"}, {"n", "m"}, {"m", "o"}}, 30, 25, 0, 29, 30, 54, "seq4", "ACTGACTGACTGACTGACTG"),
		// Muti with 4 occurrences : reject
		newCandidate([]genePair{{"p", "q"}, {"p", "r"}, {"p", "s"}}, 30, 25, 0, 29, 30, 54, "seq4", "ACTGACTGACTGACTGACTG"),
		newCandidate([]genePair{{"p", "t"}}, 30, 25, 0, 29, 30, 54, "seq4", "ACTGACTGACTGACTGACTG"),
	}

	want := []Candidate{
		candidates[0],
		candidates[1],
		candidates[2],
		candidates[3],
		candidates[7],
	}
	DiscardAbundantPartners(&candidates, 2)
	expect.EQ(t, candidates, want)
}

// More end-to-end style tests to check that behavior of this code matches the
// C++'s.
func TestDetectFusion(t *testing.T) {
	t.Skip("not enabled by default")

	ctx := vcontext.Background()
	opts := DefaultOpts
	geneDB := NewGeneDB(opts)
	geneDB.ReadFusionEvents(ctx, "/scratch-nvme/fusion/small_pairs.txt")
	geneDB.ReadTranscriptome(ctx, "/scratch-nvme/fusion/transcriptome.fa", true)

	run := func(seq string, genePair string, g1Span, g2Span, jointSpan int, order bool, g1Start, g1End, g2Start, g2End int) {
		genes := strings.Split(genePair, "/")
		frag := testNewFragment("f0", seq, opts)
		geneRangeVec := inferGeneRangeInfo(frag, geneDB, opts.KmerLength)
		if len(geneRangeVec) > opts.MaxGeneCandidatesPerFragment {
			geneRangeVec = geneRangeVec[:opts.MaxGeneCandidatesPerFragment]
		}
		fi := inferCandidatePair(geneDB, "xx", len(frag.R1Seq), len(frag.R2Seq), geneRangeVec, opts)
		expect.EQ(t, len(fi), 1)
		expect.EQ(t, fi[0], FusionInfo{
			G1ID:        geneDB.internGene(genes[0]),
			G2ID:        geneDB.internGene(genes[1]),
			G1Span:      g1Span,
			G2Span:      g2Span,
			JointSpan:   jointSpan,
			RefOrder:    true,
			FusionOrder: order,
			G1Range:     newCrossReadPosRange(Pos(g1Start), Pos(g1End)),
			G2Range:     newCrossReadPosRange(Pos(g2Start), Pos(g2End)),
		})
	}

	run("CCCAATACTCCGGCCCCTCCTGCTCTATCCACGGCGCCCGCGGCTCCATCCTCTGGCTCGCGGCGTCGCTGTCGAACCGCACGAACTGCGTGTCGTCCACGTAGCCCACGGCGATGAAGCGGGGCTCCCCGCGGCCGGGCCGGGACACGGAGGTGTAGAAAT", "HLA-C/HLA-E", 95, 39, 134, true, 0, 107, 112, 151)

	run(
		"TGGATCCCTAGCTCTTATCATGGCACTCTGTTGAGTTTGTGAAATGCATCTTCAAAGAGGTTGTCAACTGTTGCTGGAGACAACGGCTCTTCACAGACCACCTCCTTTTCTAAGGAAAATGGCTGGTATGACG",
		"IFNGR1/IL6ST",
		64, 68, 132,
		false,
		1, 65,
		65, 133)

	run("GTCCATAGCTGCTCGGTTGCCCATAGGTGTTCTGCTGAGAGTAACTGCTCTGATCATAACTAGTCGGCTGTGTAGAGGAATAGCTGGTAGGAGGGTAGGATGGAGGTGCAGTGACGGGCTATCCCCACCATCCCAATCGCAGGCTGAATTATT", "EWSR1/GNPTAB", 115, 33, 148, true, 0, 115, 120, 153)
}
