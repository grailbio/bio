package parsegencode

import (
	"context"
	"testing"

	"github.com/grailbio/testutil"
	"github.com/grailbio/testutil/assert"
	"github.com/grailbio/testutil/expect"
)

func findGene(genes []*GencodeGene, name string) *GencodeGene {
	for _, gene := range genes {
		if gene.geneID == name {
			return gene
		}
	}
	panic(name)
}

// TestReadGTFWithoutPadding will test gtf reading functionality without padded exons
func TestReadGTFWithoutPadding(t *testing.T) {
	gtfRecords := ReadGTF(context.Background(), testutil.GetFilePath(
		"//go/src/github.com/grailbio/bio/fusion/parsegencode/testdata/annotation.gtf"), false, 0, false,
		0)
	assert.EQ(t, len(gtfRecords), 4)
	// Test some random things in there
	expect.EQ(t, findGene(gtfRecords, "ENSG1.1").transcripts["ENST1.1"].havanaTranscript, "OTTHUMT1.1")
	expect.EQ(t, findGene(gtfRecords, "ENSG4.1").transcripts["ENST6.1"].exons[1].stop, 300)
	expect.EQ(t, findGene(gtfRecords, "ENSG3.1").strand, "-")
}

// TestReadGTFWithoutPadding will test gtf reading functionality with 20bp padded exons
func TestReadGTFWithPadding(t *testing.T) {
	gtfRecords := ReadGTF(context.Background(), testutil.GetFilePath(
		"//go/src/github.com/grailbio/bio/fusion/parsegencode/testdata/annotation.gtf"), false, 20, false,
		0)
	assert.EQ(t, len(gtfRecords), 4)
	// Test gene on positive strand worked correctly
	assert.EQ(t, findGene(gtfRecords, "ENSG1.1").transcripts["ENST1.1"].exons[0].start, 80)
	assert.EQ(t, findGene(gtfRecords, "ENSG1.1").transcripts["ENST1.1"].exons[0].stop, 170)
	// Test gene on negative strand worked correctly (incidentally this also checks overlaps worked)
	assert.EQ(t, findGene(gtfRecords, "ENSG2.1").transcripts["ENST3.1"].exons[0].start, 130)
	assert.EQ(t, findGene(gtfRecords, "ENSG2.1").transcripts["ENST3.1"].exons[0].stop, 320)
}

// TestReverseComplement will test ReverseComplement is correctly handling sequences
func TestReverseComplement(t *testing.T) {
	assert.EQ(t, reverseComplement("ACTG"), "CAGT")
}

// TestOverlappingGenomicRanges will test overlappingGenomicRanges
func TestGenomicRangesAppends(t *testing.T) {
	testGR := genomicRanges{}

	testGR = append(testGR, genomicRange{0, 10})  // Append to empty genomicRanges with merge=false
	testGR = append(testGR, genomicRange{15, 20}) // Append without overlap, merge=false
	testGR.merge(genomicRange{25, 35})            // Append without overlap, merge=true
	testGR = append(testGR, genomicRange{30, 45}) // Append with overlap, merge=false
	testGR.merge(genomicRange{40, 55})            // Append with overlap, merge=true

	want := genomicRanges{
		genomicRange{0, 10},
		genomicRange{15, 20},
		genomicRange{25, 35},
		genomicRange{30, 55},
	}

	for idx := range testGR {
		if testGR[idx].start != want[idx].start || testGR[idx].stop != want[idx].stop {
			t.Errorf("Error in overlapping ranges:\nGot {%d, %d}\nWant {%d, %d}", testGR[idx].start,
				testGR[idx].stop, want[idx].start, want[idx].stop)
		}
	}
	testGR = nil
	testGR = append(testGR, genomicRange{0, 10}) // Append to empty genomicRanges with merge=true
	if testGR[0].start != want[0].start || testGR[0].stop != want[0].stop {
		t.Errorf("Error in overlapping ranges:\nGot {%d, %d}\nWant {%d, %d}", testGR[0].start,
			testGR[0].stop, want[0].start, want[0].stop)
	}
}
