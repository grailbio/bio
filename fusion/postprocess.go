package fusion

import (
	"encoding/binary"
	"sort"

	"github.com/grailbio/base/log"
	"github.com/minio/highwayhash"
)

// Candidate is a combination of a fragment and possible fusions detected for
// the fragment.
type Candidate struct {
	Frag    Fragment
	Fusions []FusionInfo
}

// Return true if the substring used to infer candidate gene pair is of low
// complexity.
func LinkedByLowComplexSubstring(frag Fragment, fi FusionInfo, lowComplexityFraction float64) bool {
	g1Seq := frag.SubSeq(fi.G1Range)
	g2Seq := frag.SubSeq(fi.G2Range)
	return IsLowComplexity(g1Seq, lowComplexityFraction) || IsLowComplexity(g2Seq, lowComplexityFraction)
}

// Return true if two genes in a candidate are deemed to be within close
// proximity of each other (Either they are within prox_dist bases of each
// other, or they are within prox_num genes of each other, on the same
// chromosome)
func CloseProximity(geneDB *GeneDB, fi FusionInfo,
	maxProximityDistance, maxProximityGenes int) bool {

	gi1 := geneDB.GeneInfo(fi.G1ID)
	gi2 := geneDB.GeneInfo(fi.G2ID)

	// If any of the possible candidates are read-throughs, it explains the
	// fragment.
	if gi1.Chrom != gi2.Chrom {
		// Can't be reads through if genes are on separate chromosomes.
		//
		// TODO(saito) This should be checked a lot earlier.
		return false
	}

	if maxProximityGenes > 0 {
		numDiff := abs(gi1.Index - gi2.Index)
		if numDiff <= maxProximityGenes {
			return true
		}
	}

	if maxProximityDistance > 0 {
		if gi1.Start <= gi2.Start {
			if gi2.End < gi1.End || gi2.Start-gi1.End <= maxProximityDistance {
				return true
			}
		} else { // rhs.start < lhs.start
			if gi1.End < gi2.End || gi1.Start-gi2.End <= maxProximityDistance {
				return true
			}
		}
	}
	return false
}

type uniqueUMIs []string

type hashKey = [highwayhash.Size]uint8

// groupCandidatesByGenePair groups the candidatse by the list of genes involved in in the fusion event.
func groupCandidatesByGenePair(candidates []Candidate) map[hashKey][]int {
	var zeroSeed = hashKey{}

	hashGeneIDs := func(fusions []FusionInfo, hashBuf *[]uint8) hashKey {
		hashInt := func(v int) {
			l := len(*hashBuf)
			if cap(*hashBuf)-l < 4 {
				tmp := make([]uint8, l, 4+2*cap(*hashBuf))
				copy(tmp, *hashBuf)
				*hashBuf = tmp
			}
			buf := (*hashBuf)[l : l+4]
			(*hashBuf) = (*hashBuf)[:l+4]
			binary.LittleEndian.PutUint32(buf, uint32(v))
		}
		*hashBuf = (*hashBuf)[:0]
		for _, fi := range fusions {
			g1, g2 := fi.G1ID, fi.G2ID
			if g1 > g2 {
				g1, g2 = g2, g1
			}
			hashInt(int(g1))
			hashInt(int(g2))
		}
		return highwayhash.Sum(*hashBuf, zeroSeed[:])
	}

	var hashBuf []uint8
	candidatesMap := map[hashKey][]int{}
	for i, c := range candidates {
		h := hashGeneIDs(c.Fusions, &hashBuf)
		candidatesMap[h] = append(candidatesMap[h], i)
	}
	return candidatesMap
}

func subsetCandidates(candidates []Candidate, indices []int) []Candidate {
	sort.Slice(indices, func(i, j int) bool {
		vi := indices[i]
		vj := indices[j]
		if vi == vj {
			panic(indices)
		}
		return vi < vj
	})
	for i, ui := range indices {
		candidates[i] = candidates[ui]
	}
	return candidates[:len(indices)]
}

// Filter duplicate reads that call the same event
func FilterDuplicates(candidatesPtr *[]Candidate, hasUMI bool) {
	candidates := *candidatesPtr
	candidatesMap := groupCandidatesByGenePair(candidates)
	var (
		uniqIndices  []int
		uniqUMIs     uniqueUMIs
		validIndices []int
	)
	for _, indices := range candidatesMap {
		if hasUMI {
			uniqUMIs = uniqUMIs[:0]
			uniqUMIs = append(uniqUMIs, candidates[indices[0]].Frag.UMI())
		} else {
			uniqIndices = uniqIndices[:0]
			uniqIndices = append(uniqIndices, indices[0])
		}
		validIndices = append(validIndices, indices[0])
		for _, ci := range indices[1:] {
			cB := candidates[ci]
			if hasUMI {
				umi := cB.Frag.UMI()
				found := false
				// Check if umi is considered to be found in uniqUMIs and will replace
				// the found one in uniqUMIs if query_umi has less 'N's.
				for ui := range uniqUMIs {
					if hammingDistance(umi, uniqUMIs[ui]) <= 2 {
						found = true
						if numUnknownBases(umi) < numUnknownBases(uniqUMIs[ui]) {
							uniqUMIs[ui] = umi
						}
						if ReproduceBug {
							// TODO(saito) need to replace other UMIs too.
							break
						}
					}
				}
				if !found {
					validIndices = append(validIndices, ci)
					uniqUMIs = append(uniqUMIs, umi)
				}
			} else {
				// Arbitrarily set 5% distance threshold.
				totalLen := len(cB.Frag.R1Seq) + len(cB.Frag.R2Seq)
				maxDist := max(1, int(float64(totalLen)*0.05))
				found := false
				for _, ui := range uniqIndices {
					if d := cB.Frag.HammingDistance(candidates[ui].Frag); d <= maxDist {
						found = true
						break
					}
				}
				if !found {
					uniqIndices = append(uniqIndices, ci)
					validIndices = append(validIndices, ci)
				}
			}
		}
	}
	*candidatesPtr = subsetCandidates(*candidatesPtr, validIndices)
}

// Discard calls where one of the partners is involved in numerous events
func DiscardAbundantPartners(
	candidatesPtr *[]Candidate,
	maxGenePartners int) {
	candidates := *candidatesPtr

	// First pass, record every gene and its partners
	partners := map[GeneID]map[GeneID]struct{}{}
	addGenePair := func(g1, g2 GeneID) {
		if partners[g1] == nil {
			partners[g1] = map[GeneID]struct{}{}
		}
		partners[g1][g2] = struct{}{}
	}

	for _, c := range candidates {
		for _, fusion := range c.Fusions {
			addGenePair(fusion.G1ID, fusion.G2ID)
			addGenePair(fusion.G2ID, fusion.G1ID)
		}
	}

	// Now that we have an idea of the number of unique partners per gene, we can
	// discard those with > cap partners.
	isGeneWithAbundantPartners := func(g GeneID) bool {
		return len(partners[g]) > maxGenePartners
	}
	var (
		j                int
		nDiscardedFusion int
		nTotalFusion     int
	)
	for _, c := range candidates {
		var k int
		for _, fusion := range c.Fusions {
			nTotalFusion++
			if isGeneWithAbundantPartners(fusion.G1ID) || isGeneWithAbundantPartners(fusion.G2ID) {
				nDiscardedFusion++
				continue
			}
			c.Fusions[k] = fusion
			k++
		}
		c.Fusions = c.Fusions[:k]
		if len(c.Fusions) > 0 {
			candidates[j] = c
			j++
		}
	}
	log.Printf(
		"Discarding %d of %d candidates, %d of %d fusions for having too many fusion partners",
		len(candidates)-j, len(candidates), nDiscardedFusion, nTotalFusion)
	*candidatesPtr = candidates[:j]
}

// FilterByMinSpan filters candidates that aren't covered minSpan bases either
// G1 or G2.  It also performs UMI collapsing and returns all valid fragment
// indices. More specifically, this function scans through all candidate gene
// pairs, and considers it as a valid pair if there are at least a prescribed
// minimum number of unique fragment supporting this fusion.
func FilterByMinSpan(hasUMI bool, minSpan int, candidatesPtr *[]Candidate, minReadSupport int) {
	candidates := *candidatesPtr
	candidatesMap := groupCandidatesByGenePair(*candidatesPtr)
	var validIndices []int

	hasGoodSpan := func(fusions []FusionInfo) bool {
		for _, fi := range fusions {
			if fi.G1Span < minSpan || fi.G2Span < minSpan {
				return false
			}
		}
		return true
	}

	if hasUMI {
		for _, indices := range candidatesMap {
			var (
				numSpanSupport int
				umisSeen       = map[string]bool{}
				curIndices     []int
			)
			for _, ci := range indices {
				c := candidates[ci]
				umi := c.Frag.UMI()
				if !umisSeen[umi] {
					umisSeen[umi] = true
					curIndices = append(curIndices, ci)
					if hasGoodSpan(c.Fusions) {
						numSpanSupport++
					}
				}
			}
			if numSpanSupport >= minReadSupport {
				validIndices = append(validIndices, curIndices...)
			}
		}
	} else {
		for _, indices := range candidatesMap {
			var numSpanSupport int
			for _, ci := range indices {
				c := candidates[ci]
				if hasGoodSpan(c.Fusions) {
					numSpanSupport++
				}
			}
			if numSpanSupport >= minReadSupport {
				validIndices = append(validIndices, indices...)
			}
		}
	}
	*candidatesPtr = subsetCandidates(*candidatesPtr, validIndices)
}
