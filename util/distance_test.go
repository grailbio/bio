package util

import (
	"reflect"
	"testing"

	"github.com/antzucaro/matchr"
)

func TestOperationsContains(t *testing.T) {
	tests := []struct {
		o     operations
		given operations
		want  bool
	}{
		{operations{diagonal, right, down}, operations{diagonal}, true},
		{operations{right, down}, operations{diagonal}, false},
		{operations{diagonal, right}, operations{diagonal, right}, true},
	}

	for _, test := range tests {
		got := test.o.contains(test.given)
		if !reflect.DeepEqual(got, test.want) {
			t.Errorf("incorrect operations contains result: got %v, want %v", got, test.want)
		}
	}
}

// TestLevenshtein tests our implementation of the Levenshtein distance, which
// accounts for the fact that in the presence of deletions outnumbering
// insertions, additional primer bases downstream from the barcode will be
// read. We also test the standard calculation of the Levenshtein distance.
func TestLevenshtein(t *testing.T) {
	tests := []struct {
		barcode1    string
		barcode2    string
		downstream1 string
		downstream2 string
		want        int
	}{
		// Test 1: We test two barcodes whose optimal levenshtein distance contains a deletion
		// of the second base in barcode 1.
		// ATCGGTX (where XYZ is the downstream sequence)
		// | ||||
		// A-CGGTX (where X is read from 'downstream1' sequence)
		{"ATCGGT", "ACGGTX", "XYZ", "", 1},
		// Test 2: Same as Test1 except barcodes and accompanying downstream sequences are switched.
		{"ACGGTX", "ATCGGT", "", "XYZ", 1},
		// Test 3: Test a standard case with no deletions.
		{"ACAATTGG", "AXAAXTGX", "", "", 3},
		// Test 4: We test a case with many deletions.
		// ATATACGGT (where HIJKLMN is the downstream sequence)
		// |    ||||
		// A----CGGTHIJK (where HIJK is read from 'downstream1' sequence)
		{"ATATACGGT", "ACGGTHIJK", "HIJKLMN", "", 4},
		// Test 5: Case that has many deletions toward the end.
		{"CTCAGCGGCT", "AGCCTAACTC", "ACACTCTTTCCCTACACGACGCTCTTCCGATCT", "GTGACTGGAGTTCAGACGTGTGCTCTTCCGATC", 8},
	}

	for _, test := range tests {
		gotGrail := Levenshtein(test.barcode1, test.barcode2, test.downstream1, test.downstream2)
		grailStandard := Levenshtein(test.barcode1, test.barcode2, "", "")
		if !reflect.DeepEqual(gotGrail, test.want) {
			t.Errorf("incorrect Grail levenshtein result: got %v, want %v", gotGrail, test.want)
		}
		gotStandard := matchr.Levenshtein(test.barcode1, test.barcode2)
		if !reflect.DeepEqual(gotStandard, grailStandard) {
			t.Errorf("discrepancy between standard levenshtein and grail standard: standard levenshtein %v, grail standard %v", gotStandard, grailStandard)
		}
	}
}
