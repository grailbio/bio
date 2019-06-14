package umi

import (
	"bufio"
	"bytes"
	"fmt"
	"strings"

	"github.com/grailbio/base/log"
	"github.com/grailbio/bio/util"
)

var (
	alphabetMap = map[byte]bool{
		'A': true,
		'C': true,
		'G': true,
		'T': true,
	}

	alphabetWithN    = []byte{'A', 'C', 'G', 'T', 'N'}
	alphabetWithNMap = map[byte]bool{
		'A': true,
		'C': true,
		'G': true,
		'T': true,
		'N': true,
	}
)

func levenshteinCostFn(s1, s2 string) int {
	return util.Levenshtein(s1, s2, "", "")
}

type snapCorrectorEntry struct {
	knownUMI string
	edits    int
}

// SnapCorrector implements "snap" correction of UMIs.  A umi U is
// snappable if there is a known non-random umi U1 that is closer to U
// than all other known umis, in terms of Levenshtein edit distance.
type SnapCorrector struct {
	knownUMIs []string
	k         int

	// correctionTable contains a mapping from all snappable k-mers (k
	// is the length of the umi) to the known UMI they should snap to.
	correctionTable map[string]snapCorrectorEntry
}

// NewSnapCorrector creates a new snap corrector.  The knownUMIs are a
// \n separated list of UMIs (identical to the file content of a list
// of UMIs, where each line contains a UMI).  Each UMI should consist
// of characters ACGTN.
func NewSnapCorrector(knownUMIs []byte) *SnapCorrector {
	log.Debug.Printf("Building snappable UMI correction table")
	reader := bytes.NewBuffer(knownUMIs)
	scanner := bufio.NewScanner(reader)
	known := []string{}
	k := -1
	for scanner.Scan() {
		umi := strings.ToUpper(scanner.Text())
		if k < 0 {
			k = len(umi)
		}
		if len(umi) != k {
			panic(fmt.Sprintf("umi %s has length %d, other umis have length %d", umi, len(umi), k))
		}
		validateUMI(umi, false)

		known = append(known, umi)
	}
	if k < 0 {
		panic("no umis in input")
	}

	// Initialize the cost table.
	costTable := map[string][][]string{}
	all := allKmers(k, alphabetWithN)
	for _, s := range all {
		costTable[s] = make([][]string, k+1)
	}

	// Populate the cost table
	for _, umi := range all {
		for _, knownUMI := range known {
			cost := levenshteinCostFn(umi, knownUMI)
			if costTable[umi][cost] == nil {
				costTable[umi][cost] = make([]string, 0)
			}
			costTable[umi][cost] = append(costTable[umi][cost], knownUMI)
		}
	}

	// Find umis that can be snapped to a known umi, and save them to correctionTable.
	correctionTable := map[string]snapCorrectorEntry{}
	for umi, costList := range costTable {
		for cost, knownList := range costList {
			if len(knownList) == 1 {
				log.Debug.Printf("%s snaps to %s with cost %d", umi, knownList[0], cost)
				correctionTable[umi] = snapCorrectorEntry{knownList[0], cost}
			}
			if len(knownList) > 0 {
				break
			}
		}
	}
	log.Debug.Printf("Done building snappable UMI correction table")

	return &SnapCorrector{
		knownUMIs:       known,
		k:               k,
		correctionTable: correctionTable,
	}
}

// CorrectUMI returns a corrected umi, number of edits to the
// corrected umi, and true if there is exactly one known UMI that is
// closest to the original umi with respect to Levenshtein edit
// distance.  Otherwise, return the original umi, -1, and false.
func (c *SnapCorrector) CorrectUMI(umi string) (correctedUMI string, edits int, corrected bool) {
	umi = strings.ToUpper(umi)
	validateUMI(umi, true)
	entry, corrected := c.correctionTable[umi]
	if corrected {
		return entry.knownUMI, entry.edits, entry.knownUMI != umi
	}
	return umi, -1, false
}

func validateUMI(umi string, allowN bool) {
	for _, c := range umi {
		if (allowN && !alphabetWithNMap[byte(c)]) || (!allowN && !alphabetMap[byte(c)]) {
			panic(fmt.Sprintf("invalid base %c in umi %v", c, umi))
		}
	}
}

// returns a slice of all possible kmers with the given alphabet.
func allKmers(k int, alphabet []byte) []string {
	var fn func(partial string, length int) []string
	fn = func(partial string, length int) []string {
		if len(partial) == length {
			return []string{partial}
		}

		kmers := []string{}
		for _, c := range alphabet {
			newPartial := append([]byte(partial), c)
			kmers = append(kmers, fn(string(newPartial), length)...)
		}
		return kmers
	}

	return fn("", k)
}
