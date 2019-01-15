package fusion

import (
	"testing"

	"github.com/grailbio/testutil/expect"
)

func TestPreprocess(t *testing.T) {
	type result struct{ name, r1Seq, r2Seq string }
	doTest := func(name, r1Seq, r2Seq string, umiInRead, umiInName bool) result {
		opts := Opts{
			UMIInRead:             umiInRead,
			UMIInName:             umiInName,
			LowComplexityFraction: 0.9,
		}
		stats := Stats{}
		name, r1Seq, r2Seq = MaybeRemoveUMI(name, r1Seq, r2Seq, opts)
		r1Seq, r2Seq = RemoveLowComplexityReads(r1Seq, r2Seq, &stats, opts)
		return result{name, r1Seq, r2Seq}
	}

	// R2 is lc, umi_in_read
	expect.EQ(t, doTest("f7", "111111nGCGTGCTCTTCAG", "222222nCCCCTCCCCCCGCCCCCCCCCCCCCCCGG", true, false),
		result{"f7:111111+222222", "GCGTGCTCTTCAG", ""})

	// R2 is lc, no UMI
	expect.EQ(t, doTest("f4", "GCGTGCTCTTCAGGCGGGGA", "GCTCCGCCCCCTCCCCCCGCCCCCCCCCCCCCCCGG", false, false),
		result{"f4", "GCGTGCTCTTCAGGCGGGGA", ""})

	// Both R1 and R2 are not lc, umi_in_read
	expect.EQ(t, doTest("f1 1:N:0:CTGAAGCT+ACGTCCTG", "111111nAAAGTTCAG", "222222nTTCAGGTA", true, false),
		result{"f1:111111+222222 1:N:0:CTGAAGCT+ACGTCCTG", "AAAGTTCAG", "TACCTGAA"})

	// Both R1 and R2 are not lc, umi_in_name
	expect.EQ(t, doTest("f1 1:N:0:CTGAAGCT+ACGTCCTG", "AGCTAGCAAAGTTCAG", "AGCTAGCTTCAGGTA", false, true),
		result{"f1 1:N:0:CTGAAGCT+ACGTCCTG", "AGCTAGCAAAGTTCAG", "TACCTGAAGCTAGCT"})

	// <= 7bp, replace with n directly, has UMI
	expect.EQ(t, doTest("f2", "AAAGTTCAG", "TTCAGG", true, false), result{"f2", "N", ""})

	// no UMI case
	expect.EQ(t, doTest("f3", "AAAGTTCAG", "TTCAGG", false, false),
		result{"f3", "AAAGTTCAG", "CCTGAA"})

	// R1 is lc, no UMI
	expect.EQ(t, doTest("f5", "GCTCCGCCCCCTCCCCCCGCCCCCCCCCCCCCCCGG", "GCGTGCTCTTCAGGCGGGGA", false, false),
		result{"f5", "TCCCCGCCTGAAGAGCACGC", ""})

	// Both R1 and R2 are lc, no UMI
	expect.EQ(t, doTest("f6", "GGGGGGGGGG", "GGGGGGGGGG", false, false), result{"f6", "N", ""})

	// R1 is lc, umi_in_read
	expect.EQ(t, doTest("f8", "111111nCCCCTCCCCCCGCCCCCCCCCCCCCCCGG", "222222nGCGTGCTCTTCAG", true, false),
		result{"f8:111111+222222", "CTGAAGAGCACGC", ""})

	// Both R1 and R2 are lc, umi_in_read
	expect.EQ(t, doTest("f9", "111111nGGGGGGGGGG", "222222nGGGGGGGGGG", true, false),
		result{"f9:111111+222222", "N", ""})
}
