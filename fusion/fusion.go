package fusion

import (
	"sort"
	"strings"

	"github.com/grailbio/base/log"
)

// ReproduceBug introduces extra logic to reproduce suspicious behaviours of
// //bio/rna/fusion/ code.
const ReproduceBug = true

// geneRangeInfo represents coverage of a fragment by a gene.
type geneRangeInfo struct {
	// Gene info. It can be looked up by calling GeneDB.GeneInfo()
	geneID GeneID
	// r1Span is the sum of span for []ranges for R1. Ditto for r2Span.
	r1Span, r2Span int
	// Fragment ranges covered by this gene.
	//
	// The ranges are disjoint, and they are sorted by ascending order of start,
	// which implies that all the ranges for R1 come before those for R2.
	ranges []PosRange
}

// totalSpan returns the gr.r1Span+gr.r2Span
func (gr *geneRangeInfo) totalSpan() int { return gr.r1Span + gr.r2Span }

// FusionInfo represents a fusion event between two genes.
type FusionInfo struct {
	// G1ID is the ID of the gene 1. The gene info can be looked up by calling GeneDB.GeneInfo().
	G1ID GeneID
	// G2ID is the ID of the gene 2. The gene info can be looked up by calling GeneDB.GeneInfo().
	G2ID GeneID
	// G1Span is the total length of the gene1 that intersects with the fragment.
	G1Span int
	// G2Span is the total length of the gene2 that intersects with the fragment.
	G2Span int
	// JointSpan is the total length covered by either G1 or G2.
	//
	// REQUIRES: JointSpan >= max(G1Span, G2Span)
	JointSpan int
	// FusionOrder=true iff g1 is aligned before g2.
	FusionOrder bool
	// G1Range is the [min,max) range of the fragment covered by G1.
	G1Range CrossReadPosRange
	// G2Range is the [min,max) range of the fragment covered by G2.
	G2Range CrossReadPosRange
}

// Name returns a human-readable string that encodes key attributes of this
// fusion.
func (fi FusionInfo) Name(geneDB *GeneDB, opts Opts) string {
	buf := strings.Builder{}
	order := CosmicOrder
	if opts.Denovo {
		order = AlphabeticalOrder
	}
	g1, g2 := SortGenePair(geneDB, fi.G1ID, fi.G2ID, order)
	buf.WriteString(geneDB.GeneInfo(g1).Gene)
	buf.WriteByte('/')
	buf.WriteString(geneDB.GeneInfo(g2).Gene)
	return buf.String()
}

// Convert <kmer, gene> map and <kmer, fragment position> map, to
// <gene, fragment position> map, where the latter is then joined to generate
// gene to ranges on fragment. Taking the following form, where each range
// consists of a start and end position.
//   Gi -> [pos_start, pos_end]_i1 ... [pos_start, pos_end]_il ...
//   Gj -> [pos_start, pos_end]_j1 ...
func inferGeneRangeInfo(frag Fragment, geneDB *GeneDB, kmerLength int) []geneRangeInfo {
	geneRanges := map[GeneID][]PosRange{}

	// Generate gene_ranges: <gene_name, ranges on fragment>
	for _, km := range frag.kmers {
		k := km.minKmer()
		giter := geneDB.findByKmer(k)
		// Examine each unique gene in genes, ignoring its position.
		for idx := 0; ; idx++ {
			geneID := giter.get(idx)
			if geneID == invalidGeneID {
				break
			}
			end := km.pos + Pos(kmerLength)
			if ranges, ok := geneRanges[geneID]; ok && km.pos <= ranges[len(ranges)-1].End {
				geneRanges[geneID][len(ranges)-1].End = end
			} else {
				geneRanges[geneID] = append(geneRanges[geneID], newPosRange(km.pos, end))
			}
		}
	}

	// Convert geneRanges to geneRangeVec, sorted by reverse order in terms
	// of total span length: GeneRangeInfo.(r1_span + r2_span).
	geneRangeVec := make([]geneRangeInfo, 0, len(geneRanges))
	for geneID, ranges := range geneRanges {
		ri := geneRangeInfo{geneID: geneID}
		for _, rr := range ranges {
			switch rr.readType() {
			case R1:
				ri.r1Span += rr.span()
			case R2:
				ri.r2Span += rr.span()
			default:
				panic(rr)
			}
		}
		ri.ranges = ranges
		geneRangeVec = append(geneRangeVec, ri)
	}
	// Sort in descending order of span.
	sort.Slice(geneRangeVec, func(i, j int) bool {
		gri := geneRangeVec[i]
		grj := geneRangeVec[j]
		giSpan := gri.r1Span + gri.r2Span
		gjSpan := grj.r1Span + grj.r2Span
		if giSpan != gjSpan {
			return giSpan > gjSpan
		}
		return gri.geneID < grj.geneID
	})
	return geneRangeVec
}

// Assign one candidate pair to a given fragment.
//
// TODO(saito,xyang) consider a larger number of candidates and prioritize them
// later.
func inferCandidatePair(
	geneDB *GeneDB,
	fragName string,
	r1SeqLen, r2SeqLen int,
	geneRangeVec []geneRangeInfo,
	opts Opts) []FusionInfo {
	var fusionSpanInfoVec []FusionInfo

	// At least one gene must span minReadSpanPerc of the fragment.
	const minReadSpanPerc = 50.0
	// Stringent threshold for a fragment to be covered by a fusion gene pair.
	const fragSpanPercLarge = 90.0
	// Loose threshold for a fragment to be covered by a fusion gene pair.
	const fragSpanPercSmall = 80.0

	// avoid reads smaller than kmerlength.
	r1MinSpan := max(opts.KmerLength, r1SeqLen*minReadSpanPerc/100.0)
	r2MinSpan := max(opts.KmerLength, r2SeqLen*minReadSpanPerc/100.0)

	totalLen := r1SeqLen + r2SeqLen
	fragSpanLarge := totalLen * fragSpanPercLarge / 100
	fragSpanSmall := totalLen * fragSpanPercSmall / 100

	// When goodSpanEvidence is met, loose threshold (fragSpanPercSmall)
	// is used. This is to allow for more sequencing errors.
	//
	// TODO(saito,xyang) Revisit this constant. Can we just use 80% always?
	const goodSpanEvidence = 30

	// Track the longest span of current fragment.
	maxSpan := 0

	for gi, gr := range geneRangeVec {
		if gr.totalSpan() > totalLen-opts.KmerLength {
			// Stop searching when a single gene is determined to be solely spanning
			// the fragment.
			if !ReproduceBug {
				if len(geneRangeVec) > 1 {
					break
				}
			}
		}
		if gr.r1Span < r1MinSpan && gr.r2Span < r2MinSpan {
			break
		}

		// Pair with any remaining gene.
		for _, gr2 := range geneRangeVec[gi+1:] {
			// Determine the minimum combined span of a fragment to be considered
			// further. If both of genes have good evidence supporting the frag
			// then the threshold would be reduced.
			minFragSpan := fragSpanLarge
			if gr.totalSpan() >= goodSpanEvidence && gr2.totalSpan() >= goodSpanEvidence {
				minFragSpan = fragSpanSmall
			}
			// If the current pair has no chance of beating the current max pair,
			// skip it to reduce compute cycles.
			minFragSpan = max(maxSpan, minFragSpan)

			// No need to further search in the remaining list as geneRangeVec is
			// sorted in decreasing order by length.
			if gr.totalSpan()+gr2.totalSpan() < minFragSpan {
				break
			}

			if !opts.Denovo && !geneDB.IsFusionPair(gr.geneID, gr2.geneID) && !geneDB.IsFusionPair(gr2.geneID, gr.geneID) {
				break
			}
			// Algorithm to identify combined range by the gene pair.
			fusionSpan := inferLongestCombinedSpan(gr, gr2, opts)

			// Next further check positional information: whether
			// gr2 extends on either side of gr for at least
			// kmerlength bases.
			if fusionSpan.JointSpan-gr.totalSpan() >= opts.MinSpan {
				// This is a legitimate fusion. Ensure that the order is proper when
				// we store the information.
				// If it is a stranded prep, we need to retain the strand
				// information
				if !fusionSpan.FusionOrder {
					// In stranded prep, if the order in the span is not gr/gRV[i],
					// we need to swap the ranges (The user knows how to read the
					// results since the /- has been suffixed).
					fusionSpan.G2Span, fusionSpan.G1Span = fusionSpan.G1Span, fusionSpan.G2Span
					fusionSpan.G2Range, fusionSpan.G1Range = fusionSpan.G1Range, fusionSpan.G2Range
				}
			} else if fusionSpan.G2Span*fusionSpan.G1Span != 0 {
				// Call where one partner contributes too much, but not
				// completely. Assume that partner dominates.
				fusionSpan.G1Span = fusionSpan.JointSpan
				fusionSpan.G2Span = 0
			}

			if fusionSpan.JointSpan >= minFragSpan {
				fusionSpanInfoVec = append(fusionSpanInfoVec, fusionSpan)
				// This pair is now the largest
				if maxSpan > fusionSpan.JointSpan {
					log.Error.Printf("maxspan: %v %v %v", maxSpan, minFragSpan, fusionSpan.JointSpan)
				}
				maxSpan = fusionSpan.JointSpan
			}
		}
	}

	sort.SliceStable(fusionSpanInfoVec, func(i, j int) bool {
		return fusionSpanInfoVec[i].JointSpan > fusionSpanInfoVec[j].JointSpan
	})

	// Check if two fusion events have the same structure, modulo the gene names.
	sameFusionSpan := func(f1, f2 FusionInfo) bool {
		return f1.G1Span == f2.G1Span &&
			f1.G2Span == f2.G2Span &&
			f1.JointSpan == f2.JointSpan && (
		// Same orientation
		(f1.G1Range.Equal(f2.G1Range) && f1.G2Range.Equal(f2.G2Range)) ||
			// Different orientation
			(f1.G1Range.Equal(f2.G2Range) && f1.G2Range.Equal(f2.G1Range)))
	}

	// Now that we have a sorted vector of fusion calls, need to decide the
	// appropriate event(s) to return.
	if len(fusionSpanInfoVec) > 0 {
		// There can be many candidates with jointSpan == maxSpan and if any of
		// them are predominantly comprised of only one gene, discard this fragment.
		for _, fs := range fusionSpanInfoVec {
			if fs.JointSpan < maxSpan {
				break
			}
			if fs.G1Span == 0 || fs.G2Span == 0 {
				return nil
			}
		}
		if len(fusionSpanInfoVec) == 1 {
			// 1 fusion detected
			return fusionSpanInfoVec
		}
		fi0 := fusionSpanInfoVec[0]
		// Many potential fusions detected
		var i int
		for i = 1; i < len(fusionSpanInfoVec); i++ {
			fi := fusionSpanInfoVec[i]
			if fi.JointSpan < fi0.JointSpan {
				break
			}
			if !sameFusionSpan(fi0, fi) {
				// There were >=2 events that had the same span but had different
				// ranges.  We can't handle if it isn't one predominant case yet
				// so abort.
				//
				// TODO(arao): Fix this case
				// There are 2 possibilities
				// 1) (more probable) The fusions could potentially be A1/C, A2/C but
				// are the ranges are shifted by a few base pairs.
				// 2) (less probable) The fusions are A/B and C/D and they genuinely
				// are different events.
				//
				// TODO(saito) Handle this case by marking the fusion event in some
				// field.
				log.Error.Printf("Dropping candidates with equal length for fragment %s", fragName)
				return nil
			}
		}
		return fusionSpanInfoVec[:i]
	}
	// No fusion detected
	return nil
}

// Algorithm to infer the longest combined span on a fragment by a given gene
// pair. More specifically, gL and gR specify the ranges of two
// genes under consideration.
// The function
// returns inferred FusionSpan. If lhs -- rhs span the same order wrt r1 -- r2
// FusionSpan.order = true.
func inferLongestCombinedSpan(gL, gR geneRangeInfo, opts Opts) (fi FusionInfo) {
	fi.G1ID = gL.geneID
	fi.G2ID = gR.geneID
	fi.FusionOrder = true

	// Generate suffix and prefix lengths for lhs and rhs.
	prefixLenL := prefixLength(gL.ranges)
	prefixLenR := prefixLength(gR.ranges)
	suffixLenL := suffixLength(gL.ranges)
	suffixLenR := suffixLength(gR.ranges)

	// A flag and variable used to identify candidates that are rejected due to a
	// large overlap. Such an event is likely the suggests strong homology between
	// both candidates.
	overlap := 0
	// Identify any two pivotal ranges b/t lhs and rhs such that the two range
	// extends to the right. For each of such pair, calculate the span.
	maxSpan := 0
	var iL, iR int
	for iL < len(gL.ranges) && iR < len(gR.ranges) {
		if overlap >= opts.MaxHomology {
			// If the overlap exceeds the max allowable homology, reject this call.
			break
		}
		rangeL := gL.ranges[iL]
		rangeR := gR.ranges[iR]
		if rangeL.End <= rangeR.End {
			// L |-----|
			// R          |----|
			if rangeL.End < rangeR.Start {
				// The first condiition is for the case where two reads are not
				// stitched. We can't compute the gap, so we assume conservatively that
				// the gap between the two ranges is small.
				if (rangeL.readType() == R1 && rangeR.readType() == R2) ||
					posSpan(rangeR.Start, rangeL.End)+1 <= opts.MaxGap {
					jointSpan := prefixLenL[iL] + suffixLenR[iR]
					if jointSpan > maxSpan {
						maxSpan = jointSpan
						fi.G1Span = prefixLenL[iL]
						fi.G2Span = suffixLenR[iR]
						fi.JointSpan = maxSpan
						fi.FusionOrder = true
						fi.G1Range = newCrossReadPosRange(gL.ranges[0].Start, rangeL.End)
						fi.G2Range = newCrossReadPosRange(rangeR.Start, gR.ranges[len(gR.ranges)-1].End)
					}
				}
			} else {
				if rangeL.Start <= rangeR.Start {
					// L  |------------|
					// R       |-------------|
					//    | <jointspan>      |
					//         |overlap|
					jointSpan := posSpan(rangeR.End, rangeL.Start)
					overlap += posSpan(rangeL.End, rangeR.Start)
					if iL > 0 {
						jointSpan += prefixLenL[iL-1]
					}
					if iR < len(gR.ranges)-1 {
						jointSpan += suffixLenR[iR+1]
					}
					if jointSpan > maxSpan {
						maxSpan = jointSpan

						// Preference is giving to longer span among lhs and rhs.
						if prefixLenL[iL] >= suffixLenR[iR] {
							fi.G1Span = prefixLenL[iL]
							fi.G2Span = maxSpan - prefixLenL[iL]
							fi.JointSpan = maxSpan
							fi.FusionOrder = true
							fi.G1Range = newCrossReadPosRange(gL.ranges[0].Start, rangeL.End)
							fi.G2Range = newCrossReadPosRange(rangeL.End, gR.ranges[len(gR.ranges)-1].End)
						} else {
							fi.G1Span = maxSpan - suffixLenR[iR]
							fi.G2Span = suffixLenR[iR]
							fi.JointSpan = maxSpan
							fi.FusionOrder = true
							fi.G1Range = newCrossReadPosRange(gL.ranges[0].Start, rangeR.Start)
							fi.G2Range = newCrossReadPosRange(rangeR.Start, gR.ranges[len(gR.ranges)-1].End)
						}
					}
				} else {
					// L    |-------|
					// R   |------------|
					//
					// No need to consider a fusion in this case.  A fusion junction at
					// the rangeR.End would be always a better candidate.
					overlap += posSpan(rangeL.End, rangeL.Start)
				}
			}
			iL++
		} else {
			// L        |------|
			// R |----|
			if rangeL.Start > rangeR.End {
				// The first condiition is for the case where two reads are not
				// stitched. We can't compute the gap, so we assume conservatively that
				// the gap between the two ranges is small.
				if (rangeL.readType() == R2 && rangeR.readType() == R1) ||
					posSpan(rangeL.Start, rangeR.End)+1 <= opts.MaxGap {
					jointSpan := prefixLenR[iR] + suffixLenL[iL]
					if jointSpan > maxSpan {
						maxSpan = jointSpan
						fi.G1Span = suffixLenL[iL]
						fi.G2Span = prefixLenR[iR]
						fi.JointSpan = maxSpan
						fi.FusionOrder = false
						fi.G1Range = newCrossReadPosRange(rangeL.Start, gL.ranges[len(gL.ranges)-1].End)
						fi.G2Range = newCrossReadPosRange(gR.ranges[0].Start, rangeR.End)
					}
				}
			} else {
				if rangeL.Start >= rangeR.Start {
					// L  |-------|
					// R |------|
					jointSpan := posSpan(rangeL.End, rangeR.Start)
					overlap += posSpan(rangeR.End, rangeL.Start)
					if iL < len(gL.ranges)-1 {
						jointSpan += suffixLenL[iL+1]
					}
					if iR > 0 {
						jointSpan += prefixLenR[iR-1]
					}
					if jointSpan > maxSpan {
						maxSpan = jointSpan

						if suffixLenL[iL] >= prefixLenR[iR] {
							fi.G1Span = suffixLenL[iL]
							fi.G2Span = maxSpan - suffixLenL[iL]
							fi.JointSpan = maxSpan
							fi.FusionOrder = false
							fi.G1Range = newCrossReadPosRange(rangeL.Start, gL.ranges[len(gL.ranges)-1].End)
							fi.G2Range = newCrossReadPosRange(gR.ranges[0].Start, rangeL.Start)
						} else {
							fi.G1Span = maxSpan - prefixLenR[iR]
							fi.G2Span = prefixLenR[iR]
							fi.JointSpan = maxSpan
							fi.FusionOrder = false
							fi.G1Range = newCrossReadPosRange(rangeR.End, gL.ranges[len(gL.ranges)-1].End)
							fi.G2Range = newCrossReadPosRange(gR.ranges[0].Start, rangeR.End)
						}
					}
				} else {
					// L    |-----------|
					// R       |-------|
					//
					// No need to consider a fusion in this case.  A fusion junction at
					// the rangeL.End would be always a better candidate.
					overlap += posSpan(rangeR.End, rangeR.Start)
				}
			}
			iR++
		}
	}

	if len(gL.ranges) > 0 && prefixLenL[len(prefixLenL)-1] > maxSpan {
		maxSpan = prefixLenL[len(prefixLenL)-1]
		fi.G1Span = maxSpan
		fi.G2Span = 0
		fi.JointSpan = maxSpan
		fi.FusionOrder = true
		fi.G1Range = newCrossReadPosRange(gL.ranges[0].Start, gL.ranges[len(gL.ranges)-1].End)
		fi.G2Range = newCrossReadPosRange(0, 0)
	} else if len(gR.ranges) > 0 && prefixLenR[len(prefixLenR)-1] > maxSpan {
		maxSpan = prefixLenR[len(prefixLenR)-1]
		fi.G1Span = maxSpan
		fi.G2Span = 0
		fi.JointSpan = maxSpan
		fi.FusionOrder = true
		fi.G1Range = newCrossReadPosRange(0, 0)
		fi.G2Range = newCrossReadPosRange(gR.ranges[0].Start, gR.ranges[len(gR.ranges)-1].End)
	} else if overlap >= opts.MaxHomology {
		// Parsimonious to prefer one gene over fusion.
		//
		// Return a Null fusion
		fi = FusionInfo{G1ID: gL.geneID, G2ID: gR.geneID, FusionOrder: true}
	}
	return
}

func prefixLength(ranges []PosRange) []int {
	if len(ranges) == 0 {
		return nil
	}
	ll := make([]int, len(ranges))
	ll[0] = ranges[0].span()
	for i := 1; i < len(ranges); i++ {
		ll[i] = ll[i-1] + ranges[i].span()
	}
	return ll
}

func suffixLength(ranges []PosRange) []int {
	n := len(ranges)
	if n == 0 {
		return nil
	}
	ll := make([]int, len(ranges))
	ll[n-1] = ranges[n-1].span()
	for i := n - 2; i >= 0; i-- {
		ll[i] = ll[i+1] + ranges[i].span()
	}
	return ll
}

// DetectFusion is the toplevel entry point. It determines whether the given
// fragment is a fusion of two genes. It returns the list of candidate fusion
// events. If no event is found, it returns an empty slice.
func DetectFusion(geneDB *GeneDB, frag Fragment, stats *Stats, opts Opts) []FusionInfo {
	geneRangeVec := inferGeneRangeInfo(frag, geneDB, opts.KmerLength)
	stats.RawGenes += len(geneRangeVec)
	for _, g := range geneRangeVec {
		stats.RawRanges += len(g.ranges)
	}
	if len(geneRangeVec) > opts.MaxGeneCandidatesPerFragment {
		geneRangeVec = geneRangeVec[:opts.MaxGeneCandidatesPerFragment]
	}
	stats.Genes += len(geneRangeVec)
	for _, g := range geneRangeVec {
		stats.Ranges += len(g.ranges)
	}
	nn := len(geneRangeVec)
	if nn >= len(stats.FragmentsWithMatchingGenes) {
		nn = len(stats.FragmentsWithMatchingGenes) - 1
	}
	stats.FragmentsWithMatchingGenes[nn]++
	fusions := inferCandidatePair(geneDB, frag.Name, len(frag.R1Seq), len(frag.R2Seq), geneRangeVec, opts)
	if len(fusions) > 0 && ReproduceBug {
		// The C++ code sorts the fusions by lexicographic ordering of constituent
		// gene names.
		sort.SliceStable(fusions, func(i, j int) bool {
			f1 := &fusions[i]
			gp1 := int64(f1.G1ID<<32) | int64(f1.G2ID)
			f2 := &fusions[j]
			gp2 := int64(f2.G1ID<<32) | int64(f2.G2ID)
			return gp1 < gp2
		})
	}
	return fusions
}
