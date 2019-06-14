package bam

import (
	"reflect"
	"regexp"
	"runtime"
	"strconv"
	"testing"

	"github.com/grailbio/hts/sam"
	"github.com/stretchr/testify/assert"
)

// getFunctionName returns the runtime function name.
func getFunctionName(i interface{}) string {
	return runtime.FuncForPC(reflect.ValueOf(i).Pointer()).Name()
}

func TestFlagParser(t *testing.T) {
	// Define tests.
	tests := []struct {
		flag sam.Flags
		f    func(record *sam.Record) bool
		want bool
	}{
		// Test true behavior.
		{sam.Paired, IsPaired, true},
		{sam.ProperPair, IsProperPair, true},
		{sam.Unmapped, IsUnmapped, true},
		{sam.MateUnmapped, IsMateUnmapped, true},
		{sam.Reverse, IsReverse, true},
		{sam.MateReverse, IsMateReverse, true},
		{sam.Read1, IsRead1, true},
		{sam.Read2, IsRead2, true},
		{sam.Secondary, IsSecondary, true},
		{sam.QCFail, IsQCFail, true},
		{sam.Duplicate, IsDuplicate, true},
		{sam.Supplementary, IsSupplementary, true},
		{sam.Paired, IsPrimary, true},
		// Test false behavior.
		{sam.Supplementary, IsPaired, false},
		{sam.Duplicate, IsProperPair, false},
		{sam.QCFail, IsUnmapped, false},
		{sam.Secondary, IsMateUnmapped, false},
		{sam.Read2, IsReverse, false},
		{sam.Read1, IsMateReverse, false},
		{sam.MateReverse, IsRead1, false},
		{sam.Reverse, IsRead2, false},
		{sam.MateUnmapped, IsSecondary, false},
		{sam.Unmapped, IsQCFail, false},
		{sam.ProperPair, IsDuplicate, false},
		{sam.Paired, IsSupplementary, false},
		{sam.Secondary | sam.Supplementary, IsPrimary, false},
	}

	ref, err := sam.NewReference("chrTest", "", "", 1000, nil, nil)
	if err != nil {
		t.Fatal(err)
	}

	for _, test := range tests {
		// Make sam.Record.
		myRecord := sam.Record{
			Name: "TestRead",
			Ref:  ref,
			Pos:  0,
			MapQ: 0,
			Cigar: []sam.CigarOp{
				sam.NewCigarOp(sam.CigarMatch, 5),
			},
			Flags:   sam.Flags(test.flag),
			MateRef: ref,
			MatePos: 0,
			TempLen: 0,
			Seq:     sam.NewSeq([]byte{}),
			Qual:    []byte{},
		}

		got := test.f(&myRecord)

		if !reflect.DeepEqual(got, test.want) {
			t.Errorf("for flag %v and test %v: got %v, want %v", test.flag, getFunctionName(test.f), got, test.want)
		}
	}
}

func TestClippingDistance(t *testing.T) {
	var (
		chr1, _ = sam.NewReference("chr1", "", "", 1000, nil, nil)

		R1F = sam.Paired | sam.Read1
		R1R = sam.Paired | sam.Read1 | sam.Reverse

		c10M = []sam.CigarOp{
			sam.NewCigarOp(sam.CigarMatch, 10),
		}
		c1S8M1S = []sam.CigarOp{
			sam.NewCigarOp(sam.CigarSoftClipped, 1),
			sam.NewCigarOp(sam.CigarMatch, 8),
			sam.NewCigarOp(sam.CigarSoftClipped, 1),
		}
		c1H8M1H = []sam.CigarOp{
			sam.NewCigarOp(sam.CigarHardClipped, 1),
			sam.NewCigarOp(sam.CigarMatch, 8),
			sam.NewCigarOp(sam.CigarHardClipped, 1),
		}
		c1H1S6M1S1H = []sam.CigarOp{
			sam.NewCigarOp(sam.CigarHardClipped, 1),
			sam.NewCigarOp(sam.CigarSoftClipped, 1),
			sam.NewCigarOp(sam.CigarMatch, 6),
			sam.NewCigarOp(sam.CigarSoftClipped, 1),
			sam.NewCigarOp(sam.CigarHardClipped, 1),
		}
		c1S1H6M1H1S = []sam.CigarOp{
			sam.NewCigarOp(sam.CigarSoftClipped, 1),
			sam.NewCigarOp(sam.CigarHardClipped, 1),
			sam.NewCigarOp(sam.CigarMatch, 6),
			sam.NewCigarOp(sam.CigarHardClipped, 1),
			sam.NewCigarOp(sam.CigarSoftClipped, 1),
		}
		c1H1S1H4M1H1S1H = []sam.CigarOp{
			sam.NewCigarOp(sam.CigarHardClipped, 1),
			sam.NewCigarOp(sam.CigarSoftClipped, 1),
			sam.NewCigarOp(sam.CigarHardClipped, 1),
			sam.NewCigarOp(sam.CigarMatch, 4),
			sam.NewCigarOp(sam.CigarHardClipped, 1),
			sam.NewCigarOp(sam.CigarSoftClipped, 1),
			sam.NewCigarOp(sam.CigarHardClipped, 1),
		}
		c1S1H1S4M1S1H1S = []sam.CigarOp{
			sam.NewCigarOp(sam.CigarSoftClipped, 1),
			sam.NewCigarOp(sam.CigarHardClipped, 1),
			sam.NewCigarOp(sam.CigarSoftClipped, 1),
			sam.NewCigarOp(sam.CigarMatch, 4),
			sam.NewCigarOp(sam.CigarSoftClipped, 1),
			sam.NewCigarOp(sam.CigarHardClipped, 1),
			sam.NewCigarOp(sam.CigarSoftClipped, 1),
		}
		c2S7M1S = []sam.CigarOp{
			sam.NewCigarOp(sam.CigarSoftClipped, 2),
			sam.NewCigarOp(sam.CigarMatch, 7),
			sam.NewCigarOp(sam.CigarSoftClipped, 1),
		}
	)

	tests := []struct {
		record                     *sam.Record
		unclippedFivePrimePosition int
		unclippedStart             int
		unclippedEnd               int
		leftClipDistance           int
		rightClipDistance          int
		fivePrimeClipDistance      int
	}{
		// Forward reads
		{&sam.Record{Name: "A", Ref: chr1, Pos: 0, Flags: R1F, Cigar: c10M}, 0, 0, 9, 0, 0, 0},
		{&sam.Record{Name: "A", Ref: chr1, Pos: 0, Flags: R1F, Cigar: c1S8M1S}, -1, -1, 8, 1, 1, 1},
		{&sam.Record{Name: "A", Ref: chr1, Pos: 0, Flags: R1F, Cigar: c1H8M1H}, -1, -1, 8, 1, 1, 1},
		{&sam.Record{Name: "A", Ref: chr1, Pos: 0, Flags: R1F, Cigar: c1H1S6M1S1H}, -2, -2, 7, 2, 2, 2},
		{&sam.Record{Name: "A", Ref: chr1, Pos: 0, Flags: R1F, Cigar: c1S1H6M1H1S}, -2, -2, 7, 2, 2, 2},
		{&sam.Record{Name: "A", Ref: chr1, Pos: 0, Flags: R1F, Cigar: c1H1S1H4M1H1S1H}, -3, -3, 6, 3, 3, 3},
		{&sam.Record{Name: "A", Ref: chr1, Pos: 0, Flags: R1F, Cigar: c1S1H1S4M1S1H1S}, -3, -3, 6, 3, 3, 3},
		{&sam.Record{Name: "A", Ref: chr1, Pos: 0, Flags: R1F, Cigar: c2S7M1S}, -2, -2, 7, 2, 1, 2},

		// Reverse reads
		{&sam.Record{Name: "A", Ref: chr1, Pos: 0, Flags: R1R, Cigar: c10M}, 9, 0, 9, 0, 0, 0},
		{&sam.Record{Name: "A", Ref: chr1, Pos: 0, Flags: R1R, Cigar: c1S8M1S}, 8, -1, 8, 1, 1, 1},
		{&sam.Record{Name: "A", Ref: chr1, Pos: 0, Flags: R1R, Cigar: c1H8M1H}, 8, -1, 8, 1, 1, 1},
		{&sam.Record{Name: "A", Ref: chr1, Pos: 0, Flags: R1R, Cigar: c1H1S6M1S1H}, 7, -2, 7, 2, 2, 2},
		{&sam.Record{Name: "A", Ref: chr1, Pos: 0, Flags: R1R, Cigar: c1S1H6M1H1S}, 7, -2, 7, 2, 2, 2},
		{&sam.Record{Name: "A", Ref: chr1, Pos: 0, Flags: R1R, Cigar: c1H1S1H4M1H1S1H}, 6, -3, 6, 3, 3, 3},
		{&sam.Record{Name: "A", Ref: chr1, Pos: 0, Flags: R1R, Cigar: c1S1H1S4M1S1H1S}, 6, -3, 6, 3, 3, 3},
		{&sam.Record{Name: "A", Ref: chr1, Pos: 0, Flags: R1R, Cigar: c2S7M1S}, 7, -2, 7, 2, 1, 1},
	}

	for testIdx, test := range tests {
		t.Logf("---- starting tests[%d] ----", testIdx)
		assert.Equal(t, test.unclippedFivePrimePosition, UnclippedFivePrimePosition(test.record))
		assert.Equal(t, test.unclippedStart, UnclippedStart(test.record))
		assert.Equal(t, test.unclippedEnd, UnclippedEnd(test.record))

		assert.Equal(t, test.leftClipDistance, LeftClipDistance(test.record))
		assert.Equal(t, test.rightClipDistance, RightClipDistance(test.record))
		assert.Equal(t, test.fivePrimeClipDistance, FivePrimeClipDistance(test.record))
	}
}

func TestBaseAtPos(t *testing.T) {
	chr1, err := sam.NewReference("chr1", "", "", 1000, nil, nil)
	assert.Nil(t, err)
	sam.NewHeader(nil, []*sam.Reference{chr1})

	newRecord := func(name string, ref *sam.Reference, pos int, co []sam.CigarOp, seq string) *sam.Record {
		r, err := sam.NewRecord(name, ref, nil, pos, -1, 100, 'k', co, []byte(seq), nil, nil)
		assert.Nil(t, err)
		return r
	}

	re := regexp.MustCompile(`^(\d+)([MIDNSHP=X])(.*)`)
	makeCigar := func(cigar string) []sam.CigarOp {
		ops := []sam.CigarOp{}
		for {
			a := re.FindStringSubmatch(cigar)
			if len(a) == 0 {
				break
			}
			typ := map[string]sam.CigarOpType{
				"M": sam.CigarMatch,
				"I": sam.CigarInsertion,
				"D": sam.CigarDeletion,
				"N": sam.CigarSkipped,
				"S": sam.CigarSoftClipped,
				"H": sam.CigarHardClipped,
				"P": sam.CigarPadded,
				"=": sam.CigarEqual,
				"X": sam.CigarMismatch,
			}[a[2]]

			l, err := strconv.Atoi(a[1])
			assert.Nil(t, err)
			op := sam.NewCigarOp(typ, l)
			ops = append(ops, op)
			cigar = a[3]
		}
		return ops
	}

	tests := []struct {
		record *sam.Record
		refPos int
		base   byte
		found  bool
	}{
		{newRecord("R", chr1, 0, makeCigar("4M"), "AAAA"), -1, 0, false},
		{newRecord("R", chr1, 0, makeCigar("4M"), "CAAA"), 0, 'C', true},
		{newRecord("R", chr1, 0, makeCigar("4M"), "ACAA"), 1, 'C', true},
		{newRecord("R", chr1, 0, makeCigar("4M"), "AAAC"), 3, 'C', true},
		{newRecord("R", chr1, 0, makeCigar("4M"), "AAAA"), 5, 0, false},

		{newRecord("R", chr1, 1, makeCigar("4M"), "AAAA"), 0, 0, false},
		{newRecord("R", chr1, 1, makeCigar("4M"), "AAAA"), 5, 0, false},

		// Insertion to the reference
		{newRecord("R", chr1, 1, makeCigar("1M1I1M"), "AAA"), 0, 0, false},
		{newRecord("R", chr1, 1, makeCigar("1M1I1M"), "CAA"), 1, 'C', true},
		{newRecord("R", chr1, 1, makeCigar("1M1I1M"), "AAC"), 2, 'C', true},
		{newRecord("R", chr1, 1, makeCigar("1M1I1M"), "AAA"), 3, 0, false},

		// Deletion from the reference
		{newRecord("R", chr1, 1, makeCigar("1M1D1M"), "AA"), 0, 0, false},
		{newRecord("R", chr1, 1, makeCigar("1M1D1M"), "CA"), 1, 'C', true},
		{newRecord("R", chr1, 1, makeCigar("1M1D1M"), "AA"), 2, 0, true},
		{newRecord("R", chr1, 1, makeCigar("1M1D1M"), "AC"), 3, 'C', true},
		{newRecord("R", chr1, 1, makeCigar("1M1D1M"), "AA"), 4, 0, false},

		// Skipped from the reference
		{newRecord("R", chr1, 1, makeCigar("1M1N1M"), "AA"), 0, 0, false},
		{newRecord("R", chr1, 1, makeCigar("1M1N1M"), "CA"), 1, 'C', true},
		{newRecord("R", chr1, 1, makeCigar("1M1N1M"), "AA"), 2, 0, true},
		{newRecord("R", chr1, 1, makeCigar("1M1N1M"), "AC"), 3, 'C', true},
		{newRecord("R", chr1, 1, makeCigar("1M1N1M"), "AA"), 4, 0, false},

		// Soft clipping at the start.
		// Ref Pos 0 overlaps the read, but not in the mapped portion of the read.  In those cases, return false.
		{newRecord("R", chr1, 1, makeCigar("1S2M"), "AAA"), 0, 0, false},
		{newRecord("R", chr1, 1, makeCigar("1S2M"), "ACA"), 1, 'C', true},
		{newRecord("R", chr1, 1, makeCigar("1S2M"), "AAC"), 2, 'C', true},
		{newRecord("R", chr1, 1, makeCigar("1S2M"), "AAA"), 3, 0, false},

		// Hard clipping at the start.
		{newRecord("R", chr1, 1, makeCigar("1H2M"), "AA"), 0, 0, false},
		{newRecord("R", chr1, 1, makeCigar("1H2M"), "CA"), 1, 'C', true},
		{newRecord("R", chr1, 1, makeCigar("1H2M"), "AC"), 2, 'C', true},
		{newRecord("R", chr1, 1, makeCigar("1H2M"), "AA"), 3, 0, false},

		// Soft clipping at the end.
		{newRecord("R", chr1, 1, makeCigar("2M1S"), "AAA"), 0, 0, false},
		{newRecord("R", chr1, 1, makeCigar("2M1S"), "CAA"), 1, 'C', true},
		{newRecord("R", chr1, 1, makeCigar("2M1S"), "ACA"), 2, 'C', true},
		{newRecord("R", chr1, 1, makeCigar("2M1S"), "AAA"), 3, 0, false},
		{newRecord("R", chr1, 1, makeCigar("2M1S"), "AAA"), 4, 0, false},

		// Hard clipping at the end.
		{newRecord("R", chr1, 1, makeCigar("2M1H"), "AA"), 0, 0, false},
		{newRecord("R", chr1, 1, makeCigar("2M1H"), "CA"), 1, 'C', true},
		{newRecord("R", chr1, 1, makeCigar("2M1H"), "AC"), 2, 'C', true},
		{newRecord("R", chr1, 1, makeCigar("2M1H"), "AA"), 3, 0, false},
		{newRecord("R", chr1, 1, makeCigar("2M1H"), "AA"), 4, 0, false},

		// Padding
		{newRecord("R", chr1, 1, makeCigar("1M1P1M"), "AA"), 0, 0, false},
		{newRecord("R", chr1, 1, makeCigar("1M1P1M"), "CA"), 1, 'C', true},
		{newRecord("R", chr1, 1, makeCigar("1M1P1M"), "AC"), 2, 'C', true},
		{newRecord("R", chr1, 1, makeCigar("1M1P1M"), "AA"), 3, 0, false},
	}

	for _, test := range tests {
		base, found := BaseAtPos(test.record, test.refPos)
		assert.Equal(t, test.base, base, "Base mismatch for %s, %d", test.record, test.refPos)
		assert.Equal(t, test.found, found, "Bool mismatch for %s, %d", test.record, test.refPos)

	}
}
