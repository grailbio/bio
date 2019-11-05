package markduplicates

import (
	"bytes"
	"fmt"
	"strings"
	"sync"
	"testing"

	gbam "github.com/grailbio/bio/encoding/bam"
	"github.com/grailbio/bio/encoding/bamprovider"
	"github.com/grailbio/hts/sam"
	"github.com/grailbio/testutil"
	"github.com/stretchr/testify/assert"
)

var (
	defaultOpts = Opts{
		ShardSize:            100,
		Padding:              10,
		Parallelism:          1,
		QueueLength:          10,
		ClearExisting:        false,
		RemoveDups:           false,
		TagDups:              true,
		UseUmis:              false,
		EmitUnmodifiedFields: true,
		OpticalDetector: &TileOpticalDetector{
			OpticalDistance: 2500,
		},
	}

	chr1, _   = sam.NewReference("chr1", "", "", 1000, nil, nil)
	chr2, _   = sam.NewReference("chr2", "", "", 2000, nil, nil)
	header, _ = sam.NewHeader(nil, []*sam.Reference{chr1, chr2})

	r1F = sam.Paired | sam.Read1
	r1R = sam.Paired | sam.Read1 | sam.Reverse
	r2F = sam.Paired | sam.Read2
	r2R = sam.Paired | sam.Read2 | sam.Reverse
	s1F = sam.Paired | sam.Read1 | sam.MateUnmapped
	s2R = sam.Paired | sam.Read2 | sam.Reverse | sam.MateUnmapped
	u1  = sam.Paired | sam.Read1 | sam.Unmapped
	u2  = sam.Paired | sam.Read2 | sam.Unmapped
	sec = sam.Paired | sam.Read1 | sam.Secondary
	up1 = sam.Paired | sam.Read1 | sam.Unmapped | sam.MateUnmapped
	up2 = sam.Paired | sam.Read2 | sam.Unmapped | sam.MateUnmapped

	cigar0 = []sam.CigarOp{
		sam.NewCigarOp(sam.CigarMatch, 10),
	}
	cigar100M = []sam.CigarOp{
		sam.NewCigarOp(sam.CigarMatch, 100),
	}
	cigarSoft1 = []sam.CigarOp{
		sam.NewCigarOp(sam.CigarSoftClipped, 1),
		sam.NewCigarOp(sam.CigarMatch, 8),
		sam.NewCigarOp(sam.CigarSoftClipped, 1),
	}
	cigarHard1 = []sam.CigarOp{
		sam.NewCigarOp(sam.CigarHardClipped, 1),
		sam.NewCigarOp(sam.CigarMatch, 8),
		sam.NewCigarOp(sam.CigarHardClipped, 1),
	}
	cigarSoft2 = []sam.CigarOp{
		sam.NewCigarOp(sam.CigarSoftClipped, 2),
		sam.NewCigarOp(sam.CigarMatch, 6),
		sam.NewCigarOp(sam.CigarSoftClipped, 2),
	}
	cigar1M = []sam.CigarOp{
		sam.NewCigarOp(sam.CigarMatch, 1),
	}

	// Reads for testing duplicate marking.
	// The following duplicate group (basic) is entirely in the same shard.
	basicA1 = NewRecord("A:::1:10:1:1", chr1, 0, r1F, 10, chr1, cigar0)
	basicA2 = NewRecord("A:::1:10:1:1", chr1, 10, r2F, 0, chr1, cigar0)
	// Exact duplicate, no clipping
	basicB1 = NewRecord("B:::1:10:2:2", chr1, 0, r1F, 10, chr1, cigar0)
	basicB2 = NewRecord("B:::1:10:2:2", chr1, 10, r2F, 0, chr1, cigar0)
	// Duplicate with soft clipping
	basicC1 = NewRecord("C:::1:10:3:3", chr1, 1, r1F, 11, chr1, cigarSoft1)
	basicC2 = NewRecord("C:::1:10:3:3", chr1, 11, r2F, 1, chr1, cigarSoft1)
	// Duplicate with hard clipping
	basicD1 = NewRecord("D:::1:10:4:4", chr1, 1, r1F, 11, chr1, cigarHard1)
	basicD2 = NewRecord("D:::1:10:4:4", chr1, 11, r2F, 1, chr1, cigarHard1)

	nonDup1 = NewRecord("X:::1:10:4:4", chr1, 12, r1F, 67, chr1, cigar0)
	nonDup2 = NewRecord("X:::1:10:4:4", chr1, 67, r2R, 12, chr1, cigar0)

	single = NewRecord("sgl:::1:10:4:4", chr1, 0, s1F, 10, chr1, cigar0)

	// A1 is in shard0, A2 is outside of shard0.
	padA1 = NewRecord("A:::1:11:2:2", chr1, 50, r1F, 115, chr1, cigar0)
	padA2 = NewRecord("A:::1:11:2:2", chr1, 115, r2F, 50, chr1, cigar0)
	padB1 = NewRecord("B:::1:11:2:2", chr1, 50, r1F, 115, chr1, cigar0)
	padB2 = NewRecord("B:::1:11:2:2", chr1, 115, r2F, 50, chr1, cigar0)

	// E and F are duplicates, but they both have a negative unclipped position of -1.
	negE1 = NewRecord("E:::1:10:5:5", chr1, 0, r1F, 10, chr1, cigarSoft1)
	negE2 = NewRecord("E:::1:10:5:5", chr1, 10, r2F, 0, chr1, cigar0)
	negF1 = NewRecord("F:::1:10:6:6", chr1, 1, r1F, 10, chr1, cigarSoft2)
	negF2 = NewRecord("F:::1:10:6:6", chr1, 10, r2F, 1, chr1, cigar0)

	//  G has alignment position in shard, and H has alignment position in padding, but they are still duplicates.
	crossclipG1 = NewRecord("G:::1:10:5:5", chr1, 99, r1F, 115, chr1, cigar0)
	crossclipG2 = NewRecord("G:::1:10:5:5", chr1, 115, r2F, 99, chr1, cigar0)
	crossclipH1 = NewRecord("H:::1:10:6:6", chr1, 101, r1F, 115, chr1, cigarSoft2)
	crossclipH2 = NewRecord("H:::1:10:6:6", chr1, 115, r2F, 101, chr1, cigar0)

	//  Read1 is in shard, but Read2 is distant but same reference.
	distantDupI1 = NewRecord("I:::1:10:6:6", chr1, 50, r1F, 150, chr1, cigar0)
	distantDupI2 = NewRecord("I:::1:10:6:6", chr1, 150, r2F, 50, chr1, cigar0)
	distantDupJ1 = NewRecord("J:::1:10:6:6", chr1, 50, r1F, 150, chr1, cigar0)
	distantDupJ2 = NewRecord("J:::1:10:6:6", chr1, 150, r2F, 50, chr1, cigar0)

	//  Read1 is in shard, but Read2 is distant in a different reference.
	distantDupK1 = NewRecord("K:::1:10:6:6", chr1, 50, r1F, 150, chr2, cigar0)
	distantDupK2 = NewRecord("K:::1:10:6:6", chr2, 150, r2F, 50, chr1, cigar0)
	distantDupL1 = NewRecord("L:::1:10:6:6", chr1, 50, r1F, 150, chr2, cigar0)
	distantDupL2 = NewRecord("L:::1:10:6:6", chr2, 150, r2F, 50, chr1, cigar0)

	//  Read1 is in shard, but Read2 is distant in a different reference (but in the same numerical shard range as Read1).
	//    Explanation:
	//           Ref  Pos Shard
	//    Read1  chr1 50  chr1:0-100
	//    Read2  chr2 55  chr2:0-100 -- Note that Read2 is between 0 and 100,
	//                                  and could be mistaken to be inside shard chr1:0-100
	//                                  if we forget to compare Read2's reference which is
	//                                  chr2 when checking if Read2 is in shard chr1:0-100.
	distantDupM1 = NewRecord("M:::1:10:6:6", chr1, 50, r1F, 55, chr2, cigar0)
	distantDupM2 = NewRecord("M:::1:10:6:6", chr2, 55, r2F, 50, chr1, cigar0)
	distantDupN1 = NewRecord("N:::1:10:6:6", chr1, 50, r1F, 55, chr2, cigar0)
	distantDupN2 = NewRecord("N:::1:10:6:6", chr2, 55, r2F, 50, chr1, cigar0)

	// Read pair that needs to be swapped for left and right according to unclipped 5' pos.
	reverseA1 = NewRecord("revA:::1:10:5:5", chr1, 0, r1R, 5, chr1, cigar0)
	reverseA2 = NewRecord("revA:::1:10:5:5", chr1, 5, r2F, 0, chr1, cigar0)
	reverseB1 = NewRecord("revB:::2:10:5:5", chr1, 0, r1R, 5, chr1, cigar0)
	reverseB2 = NewRecord("revB:::2:10:5:5", chr1, 5, r2F, 0, chr1, cigar0)

	// Test case for overlapping 5' start position, with opposing directions.
	overA1 = NewRecord("overA:::1:10:5:5", chr1, 50, r1R, 59, chr1, cigar0)
	overA2 = NewRecord("overA:::1:10:5:5", chr1, 59, r2F, 50, chr1, cigar0)
	overB1 = NewRecord("overB:::2:10:5:5", chr1, 50, r1R, 59, chr1, cigar0)
	overB2 = NewRecord("overB:::2:10:5:5", chr1, 59, r2F, 50, chr1, cigar0)

	// Test case for overlapping 5' start position, with opposing directions, when left mate is in prior shard.
	overC1 = NewRecord("overC:::1:10:5:5", chr1, 50, r1R, 149, chr1, cigar100M)
	overC2 = NewRecord("overC:::1:10:5:5", chr1, 149, r2F, 50, chr1, cigar0)
	overD1 = NewRecord("overD:::2:10:5:5", chr1, 50, r1R, 149, chr1, cigar100M)
	overD2 = NewRecord("overD:::2:10:5:5", chr1, 149, r2F, 50, chr1, cigar0)
)

func TestBasicDuplicates(t *testing.T) {
	cases := []TestCase{
		{
			[]TestRecord{
				{R: basicA1, DupFlag: false, ExpectedAuxs: []sam.Aux{NewAux("DI", "0"), NewAux("DS", 4)}},
				{R: basicB1, DupFlag: true, ExpectedAuxs: []sam.Aux{NewAux("DI", "0"), NewAux("DS", 4)}},
				{R: basicC1, DupFlag: true, ExpectedAuxs: []sam.Aux{NewAux("DI", "0"), NewAux("DS", 4)}},
				{R: basicD1, DupFlag: true, ExpectedAuxs: []sam.Aux{NewAux("DI", "0"), NewAux("DS", 4)}},
				{R: basicA2, DupFlag: false, ExpectedAuxs: []sam.Aux{NewAux("DI", "0"), NewAux("DS", 4)}},
				{R: basicB2, DupFlag: true, ExpectedAuxs: []sam.Aux{NewAux("DI", "0"), NewAux("DS", 4)}},
				{R: basicC2, DupFlag: true, ExpectedAuxs: []sam.Aux{NewAux("DI", "0"), NewAux("DS", 4)}},
				{R: basicD2, DupFlag: true, ExpectedAuxs: []sam.Aux{NewAux("DI", "0"), NewAux("DS", 4)}},
			},
			defaultOpts,
		},
		{
			[]TestRecord{
				{R: padA1, DupFlag: false, ExpectedAuxs: []sam.Aux{NewAux("DI", "0"), NewAux("DS", 2)}},
				{R: padB1, DupFlag: true, ExpectedAuxs: []sam.Aux{NewAux("DI", "0"), NewAux("DS", 2)}},
				{R: padA2, DupFlag: false, ExpectedAuxs: []sam.Aux{NewAux("DI", "0"), NewAux("DS", 2)}},
				{R: padB2, DupFlag: true, ExpectedAuxs: []sam.Aux{NewAux("DI", "0"), NewAux("DS", 2)}},
			},
			defaultOpts,
		},
		{
			[]TestRecord{
				{R: negE1, DupFlag: false, ExpectedAuxs: []sam.Aux{NewAux("DI", "0"), NewAux("DS", 2)}},
				{R: negF1, DupFlag: true, ExpectedAuxs: []sam.Aux{NewAux("DI", "0"), NewAux("DS", 2)}},
				{R: negE2, DupFlag: false, ExpectedAuxs: []sam.Aux{NewAux("DI", "0"), NewAux("DS", 2)}},
				{R: negF2, DupFlag: true, ExpectedAuxs: []sam.Aux{NewAux("DI", "0"), NewAux("DS", 2)}},
			},
			defaultOpts,
		},
		{
			[]TestRecord{
				{R: crossclipG1, DupFlag: false, ExpectedAuxs: []sam.Aux{NewAux("DI", "0"), NewAux("DS", 2)}},
				{R: crossclipH1, DupFlag: true, ExpectedAuxs: []sam.Aux{NewAux("DI", "0"), NewAux("DS", 2)}},
				{R: crossclipG2, DupFlag: false, ExpectedAuxs: []sam.Aux{NewAux("DI", "0"), NewAux("DS", 2)}},
				{R: crossclipH2, DupFlag: true, ExpectedAuxs: []sam.Aux{NewAux("DI", "0"), NewAux("DS", 2)}},
			},
			defaultOpts,
		},
		{
			[]TestRecord{
				{R: distantDupI1, DupFlag: false, ExpectedAuxs: []sam.Aux{NewAux("DI", "0"), NewAux("DS", 2)}},
				{R: distantDupJ1, DupFlag: true, ExpectedAuxs: []sam.Aux{NewAux("DI", "0"), NewAux("DS", 2)}},
				{R: distantDupI2, DupFlag: false, ExpectedAuxs: []sam.Aux{NewAux("DI", "0"), NewAux("DS", 2)}},
				{R: distantDupJ2, DupFlag: true, ExpectedAuxs: []sam.Aux{NewAux("DI", "0"), NewAux("DS", 2)}},
			},
			defaultOpts,
		},
		{
			[]TestRecord{
				{R: distantDupK1, DupFlag: false, ExpectedAuxs: []sam.Aux{NewAux("DI", "0"), NewAux("DS", 2)}},
				{R: distantDupL1, DupFlag: true, ExpectedAuxs: []sam.Aux{NewAux("DI", "0"), NewAux("DS", 2)}},
				{R: distantDupK2, DupFlag: false, ExpectedAuxs: []sam.Aux{NewAux("DI", "0"), NewAux("DS", 2)}},
				{R: distantDupL2, DupFlag: true, ExpectedAuxs: []sam.Aux{NewAux("DI", "0"), NewAux("DS", 2)}},
			},
			defaultOpts,
		},
		{
			[]TestRecord{
				{R: distantDupM1, DupFlag: false, ExpectedAuxs: []sam.Aux{NewAux("DI", "0"), NewAux("DS", 2)}},
				{R: distantDupN1, DupFlag: true, ExpectedAuxs: []sam.Aux{NewAux("DI", "0"), NewAux("DS", 2)}},
				{R: distantDupM2, DupFlag: false, ExpectedAuxs: []sam.Aux{NewAux("DI", "0"), NewAux("DS", 2)}},
				{R: distantDupN2, DupFlag: true, ExpectedAuxs: []sam.Aux{NewAux("DI", "0"), NewAux("DS", 2)}},
			},
			defaultOpts,
		},
		{
			[]TestRecord{
				{R: basicA1, DupFlag: false, ExpectedAuxs: []sam.Aux{NewAux("DI", "0"), NewAux("DS", 1)},
					UnexpectedTags: []sam.Tag{sam.NewTag("DT")}},
				{R: basicA2, DupFlag: false, ExpectedAuxs: []sam.Aux{NewAux("DI", "0"), NewAux("DS", 1)},
					UnexpectedTags: []sam.Tag{sam.NewTag("DT")}},
				{R: nonDup1, DupFlag: false, ExpectedAuxs: []sam.Aux{NewAux("DI", "2"), NewAux("DS", 1)},
					UnexpectedTags: []sam.Tag{sam.NewTag("DT")}},
				{R: nonDup2, DupFlag: false, ExpectedAuxs: []sam.Aux{NewAux("DI", "2"), NewAux("DS", 1)},
					UnexpectedTags: []sam.Tag{sam.NewTag("DT")}},
			},
			defaultOpts,
		},
		{
			// If there is just one pair, and single in the dupSet, the single gets just DT.
			[]TestRecord{
				{R: basicA1, DupFlag: false, ExpectedAuxs: []sam.Aux{NewAux("DI", "0"), NewAux("DS", 1), NewAux("DL", 0)},
					UnexpectedTags: []sam.Tag{sam.NewTag("DT")}},
				{R: single, DupFlag: true, ExpectedAuxs: []sam.Aux{NewAux("DT", "LB")}, UnexpectedTags: []sam.Tag{sam.NewTag("DI"), sam.NewTag("DL")}},
				{R: basicA2, DupFlag: false, ExpectedAuxs: []sam.Aux{NewAux("DI", "0"), NewAux("DS", 1), NewAux("DL", 0)},
					UnexpectedTags: []sam.Tag{sam.NewTag("DT")}},
			},
			defaultOpts,
		},
		{
			// DI should be set equal to the leftmost 5' read position of the primary pair.  Here, reverseB is the primary
			// because reverseB2's unclipped 5' start is left most, and reverseB2 has a smaller file index than reverseA2.
			// DI should be set to 3, the file index of reverseB2, because reverseB2's 5' position is left of reverseB1's
			// 5' position.
			[]TestRecord{
				{R: single, DupFlag: false, UnexpectedTags: []sam.Tag{sam.NewTag("DI"), sam.NewTag("DS"), sam.NewTag("DT")}},
				{R: reverseA1, DupFlag: true, ExpectedAuxs: []sam.Aux{NewAux("DI", "3"), NewAux("DS", 2), NewAux("DT", "LB")}},
				{R: reverseB1, DupFlag: false, ExpectedAuxs: []sam.Aux{NewAux("DI", "3"), NewAux("DS", 2)},
					UnexpectedTags: []sam.Tag{sam.NewTag("DT")}},
				{R: reverseB2, DupFlag: false, ExpectedAuxs: []sam.Aux{NewAux("DI", "3"), NewAux("DS", 2)},
					UnexpectedTags: []sam.Tag{sam.NewTag("DT")}},
				{R: reverseA2, DupFlag: true, ExpectedAuxs: []sam.Aux{NewAux("DI", "3"), NewAux("DS", 2), NewAux("DT", "LB")}},
			},
			defaultOpts,
		},
		{
			[]TestRecord{
				{R: overA1, DupFlag: false, ExpectedAuxs: []sam.Aux{NewAux("DI", "0"), NewAux("DS", 2)},
					UnexpectedTags: []sam.Tag{sam.NewTag("DT")}},
				{R: overB1, DupFlag: true, ExpectedAuxs: []sam.Aux{NewAux("DI", "0"), NewAux("DS", 2), NewAux("DT", "LB")}},
				{R: overA2, DupFlag: false, ExpectedAuxs: []sam.Aux{NewAux("DI", "0"), NewAux("DS", 2)},
					UnexpectedTags: []sam.Tag{sam.NewTag("DT")}},
				{R: overB2, DupFlag: true, ExpectedAuxs: []sam.Aux{NewAux("DI", "0"), NewAux("DS", 2), NewAux("DT", "LB")}},
			},
			defaultOpts,
		},
		{
			// Make sure that readPair swap works when the unclipped positions in a pair are equal, but pairs straddle two
			// shards.  The readPair swap is supposed to break ties with the read's fileIndex. DI should be 0.  It's not
			// clear that this can happen in real life since the read length would need to be longer than the total padding
			// distance.
			[]TestRecord{
				{R: overC1, DupFlag: false, ExpectedAuxs: []sam.Aux{NewAux("DI", "0"), NewAux("DS", 2)},
					UnexpectedTags: []sam.Tag{sam.NewTag("DT")}},
				{R: overD1, DupFlag: true, ExpectedAuxs: []sam.Aux{NewAux("DI", "0"), NewAux("DS", 2), NewAux("DT", "LB")}},
				{R: overC2, DupFlag: false, ExpectedAuxs: []sam.Aux{NewAux("DI", "0"), NewAux("DS", 2)},
					UnexpectedTags: []sam.Tag{sam.NewTag("DT")}},
				{R: overD2, DupFlag: true, ExpectedAuxs: []sam.Aux{NewAux("DI", "0"), NewAux("DS", 2), NewAux("DT", "LB")}},
			},
			Opts{
				ShardSize:            100,
				Padding:              100,
				Parallelism:          1,
				QueueLength:          10,
				ClearExisting:        false,
				RemoveDups:           false,
				TagDups:              true,
				UseUmis:              false,
				EmitUnmodifiedFields: true,
				OpticalDetector: &TileOpticalDetector{
					OpticalDistance: 2500,
				},
			},
		},
		{
			// if all reads in two pairs all have the same ref and position, and orientation, but r1 and r2 are swapped, the
			// two pairs are still duplicates.
			[]TestRecord{
				{R: NewRecord("A:1:1:1:1:1:1", chr1, 10, r1F, 10, chr1, cigar1M), DupFlag: false},
				{R: NewRecord("B:1:1:1:1:1:1", chr1, 10, r1R, 10, chr1, cigar1M), DupFlag: true},
				{R: NewRecord("A:1:1:1:1:1:1", chr1, 10, r2R, 10, chr1, cigar1M), DupFlag: false},
				{R: NewRecord("B:1:1:1:1:1:1", chr1, 10, r2F, 10, chr1, cigar1M), DupFlag: true},
			},
			defaultOpts,
		},
		{
			// A single mapped read that is not a duplicate of anything should have no tags.
			[]TestRecord{
				{R: single, DupFlag: false,
					UnexpectedTags: []sam.Tag{sam.NewTag("DI"), sam.NewTag("DL"), sam.NewTag("DT"), sam.NewTag("DS")}},
			},
			defaultOpts,
		},
		{
			// A single read with no dups should have DS=1.
			[]TestRecord{
				{R: basicA1, DupFlag: false, ExpectedAuxs: []sam.Aux{NewAux("DI", "0"), NewAux("DS", 1)}},
				{R: basicA2, DupFlag: false, ExpectedAuxs: []sam.Aux{NewAux("DI", "0"), NewAux("DS", 1)}},
			},
			defaultOpts,
		},
		{
			// bagsize=2, with a PCR duplicate should have DL=1.
			[]TestRecord{
				{R: NewRecord("A:::1:1000:6:6", chr1, 50, r1F, 55, chr2, cigar0), DupFlag: false,
					ExpectedAuxs: []sam.Aux{NewAux("DI", "0"), NewAux("DS", 2), NewAux("DL", 1)}},
				{R: NewRecord("B:::1:1000:6:5000", chr1, 50, r1F, 55, chr2, cigar0), DupFlag: true,
					ExpectedAuxs: []sam.Aux{NewAux("DI", "0"), NewAux("DS", 2), NewAux("DL", 1)}},
				{R: NewRecord("A:::1:1000:6:6", chr2, 55, r2F, 50, chr1, cigar0), DupFlag: false,
					ExpectedAuxs: []sam.Aux{NewAux("DI", "0"), NewAux("DS", 2), NewAux("DL", 1)}},
				{R: NewRecord("B:::1:1000:6:5000", chr2, 55, r2F, 50, chr1, cigar0), DupFlag: true,
					ExpectedAuxs: []sam.Aux{NewAux("DI", "0"), NewAux("DS", 2), NewAux("DL", 1)}},
			},
			defaultOpts,
		},
		{
			// bagsize=2, with an optical duplicate should have DL=0.
			[]TestRecord{
				{R: NewRecord("A:::1:1000:6:6", chr1, 50, r1F, 55, chr2, cigar0), DupFlag: false,
					ExpectedAuxs: []sam.Aux{NewAux("DI", "0"), NewAux("DS", 2), NewAux("DL", 0)}},
				{R: NewRecord("B:::1:1000:6:7", chr1, 50, r1F, 55, chr2, cigar0), DupFlag: true,
					ExpectedAuxs: []sam.Aux{NewAux("DI", "0"), NewAux("DS", 2), NewAux("DL", 0)}},
				{R: NewRecord("A:::1:1000:6:6", chr2, 55, r2F, 50, chr1, cigar0), DupFlag: false,
					ExpectedAuxs: []sam.Aux{NewAux("DI", "0"), NewAux("DS", 2), NewAux("DL", 0)}},
				{R: NewRecord("B:::1:1000:6:7", chr2, 55, r2F, 50, chr1, cigar0), DupFlag: true,
					ExpectedAuxs: []sam.Aux{NewAux("DI", "0"), NewAux("DS", 2), NewAux("DL", 0)}},
			},
			defaultOpts,
		},
		{
			// bagsize=3, with one optical duplicate should have DL=1.
			[]TestRecord{
				{R: NewRecord("A:::1:1000:6:6", chr1, 50, r1F, 55, chr2, cigar0), DupFlag: false,
					ExpectedAuxs: []sam.Aux{NewAux("DI", "0"), NewAux("DS", 3), NewAux("DL", 1)}},
				{R: NewRecord("B:::1:1000:6:7", chr1, 50, r1F, 55, chr2, cigar0), DupFlag: true,
					ExpectedAuxs: []sam.Aux{NewAux("DI", "0"), NewAux("DS", 3), NewAux("DL", 1)}},
				{R: NewRecord("C:::1:1000:6:5000", chr1, 50, r1F, 55, chr2, cigar0), DupFlag: true,
					ExpectedAuxs: []sam.Aux{NewAux("DI", "0"), NewAux("DS", 3), NewAux("DL", 1)}},
				{R: NewRecord("A:::1:1000:6:6", chr2, 55, r2F, 50, chr1, cigar0), DupFlag: false,
					ExpectedAuxs: []sam.Aux{NewAux("DI", "0"), NewAux("DS", 3), NewAux("DL", 1)}},
				{R: NewRecord("B:::1:1000:6:7", chr2, 55, r2F, 50, chr1, cigar0), DupFlag: true,
					ExpectedAuxs: []sam.Aux{NewAux("DI", "0"), NewAux("DS", 3), NewAux("DL", 1)}},
				{R: NewRecord("C:::1:1000:6:5000", chr2, 55, r2F, 50, chr1, cigar0), DupFlag: true,
					ExpectedAuxs: []sam.Aux{NewAux("DI", "0"), NewAux("DS", 3), NewAux("DL", 1)}},
			},
			defaultOpts,
		},
	}
	RunTestCases(t, header, cases)
}

// Test that tags are not present when tagDups is false.
func TestTagDups(t *testing.T) {
	noTags := defaultOpts
	noTags.TagDups = false

	cases := []TestCase{
		{
			[]TestRecord{
				{R: basicA1, DupFlag: false, UnexpectedTags: []sam.Tag{sam.NewTag("DI"), sam.NewTag("DS"), sam.NewTag("DT"), sam.NewTag("DU")}},
				{R: basicB1, DupFlag: true, UnexpectedTags: []sam.Tag{sam.NewTag("DI"), sam.NewTag("DS"), sam.NewTag("DT"), sam.NewTag("DU")}},
				{R: basicA2, DupFlag: false, UnexpectedTags: []sam.Tag{sam.NewTag("DI"), sam.NewTag("DS"), sam.NewTag("DT"), sam.NewTag("DU")}},
				{R: basicB2, DupFlag: true, UnexpectedTags: []sam.Tag{sam.NewTag("DI"), sam.NewTag("DS"), sam.NewTag("DT"), sam.NewTag("DU")}},
			},
			noTags,
		},
		{
			[]TestRecord{
				{R: basicA1, DupFlag: false, ExpectedAuxs: []sam.Aux{NewAux("DI", "0"), NewAux("DS", 2)}},
				{R: basicB1, DupFlag: true, ExpectedAuxs: []sam.Aux{NewAux("DI", "0"), NewAux("DS", 2), NewAux("DT", "SQ")}},
				{R: basicA2, DupFlag: false, ExpectedAuxs: []sam.Aux{NewAux("DI", "0"), NewAux("DS", 2)}},
				{R: basicB2, DupFlag: true, ExpectedAuxs: []sam.Aux{NewAux("DI", "0"), NewAux("DS", 2), NewAux("DT", "SQ")}},
			},
			defaultOpts,
		},
	}
	RunTestCases(t, header, cases)
}

// Test that tags are not present when clear-existing is true.
func TestClearExisting(t *testing.T) {
	opts := defaultOpts
	opts.ClearExisting = true
	opts.TagDups = false

	// B is marked as a duplicate on the input, but A and B are not
	// duplicates.  This test checks that B's duplicate flag and aux
	// tags are not set in the output.  A and B are distant mates to
	// exercise that the flag clearing works on distant mates.
	a1 := NewRecord("A:::1:10:6:6", chr1, 50, r1F, 150, chr1, cigar0)
	a2 := NewRecord("A:::1:10:6:6", chr1, 150, r2F, 50, chr1, cigar0)
	b1 := NewRecord("B:::1:10:6:6", chr1, 50, r1F, 150, chr1, cigar0)
	b2 := NewRecord("B:::1:10:6:6", chr1, 151, r2F, 50, chr1, cigar0)

	b1.Flags |= sam.Duplicate
	aux, err := sam.NewAux(sam.NewTag("DI"), 123)
	assert.Nil(t, err)
	b1.AuxFields = append(b1.AuxFields, aux)
	aux, err = sam.NewAux(sam.NewTag("DL"), 4)
	assert.Nil(t, err)
	b1.AuxFields = append(b1.AuxFields, aux)

	b2.Flags |= sam.Duplicate
	aux, err = sam.NewAux(sam.NewTag("DI"), 123)
	assert.Nil(t, err)
	b2.AuxFields = append(b2.AuxFields, aux)
	aux, err = sam.NewAux(sam.NewTag("DL"), 4)
	assert.Nil(t, err)
	b2.AuxFields = append(b2.AuxFields, aux)

	cases := []TestCase{
		{
			[]TestRecord{
				{R: a1, DupFlag: false, UnexpectedTags: []sam.Tag{sam.NewTag("DI"), sam.NewTag("DL")}},
				{R: b1, DupFlag: false, UnexpectedTags: []sam.Tag{sam.NewTag("DI"), sam.NewTag("DL")}},
				{R: a2, DupFlag: false, UnexpectedTags: []sam.Tag{sam.NewTag("DI"), sam.NewTag("DL")}},
				{R: b2, DupFlag: false, UnexpectedTags: []sam.Tag{sam.NewTag("DI"), sam.NewTag("DL")}},
			},
			opts,
		},
	}
	RunTestCases(t, header, cases)
}

func TestExactUmis(t *testing.T) {
	useUmis := defaultOpts
	useUmis.UseUmis = true

	cases := []TestCase{
		{
			// Use Umis, basic exact match Umis
			[]TestRecord{
				{R: NewRecord("A:1:1:1:1:1:1:AAC+CCG", chr1, 0, r1F, 10, chr1, cigar0),
					DupFlag: false, ExpectedAuxs: []sam.Aux{NewAux("DI", "0"), NewAux("DS", 2)}},
				{R: NewRecord("B:1:1:1:1:1:1:AAC+CCG", chr1, 0, r1F, 10, chr1, cigar0),
					DupFlag: true, ExpectedAuxs: []sam.Aux{NewAux("DI", "0"), NewAux("DS", 2)}},
				{R: NewRecord("A:1:1:1:1:1:1:AAC+CCG", chr1, 10, r2R, 0, chr1, cigar0),
					DupFlag: false, ExpectedAuxs: []sam.Aux{NewAux("DI", "0"), NewAux("DS", 2)}},
				{R: NewRecord("B:1:1:1:1:1:1:AAC+CCG", chr1, 10, r2R, 0, chr1, cigar0),
					DupFlag: true, ExpectedAuxs: []sam.Aux{NewAux("DI", "0"), NewAux("DS", 2)}},
			},
			useUmis,
		},
		{
			// equal umis should not make up for mismatched positions
			[]TestRecord{
				{R: NewRecord("A:1:1:1:1:1:1:AAC+CCG", chr1, 0, r1F, 10, chr1, cigar0), DupFlag: false},
				{R: NewRecord("B:1:1:1:1:1:1:AAC+CCG", chr1, 1, r1F, 10, chr1, cigar0), DupFlag: false},
				{R: NewRecord("A:1:1:1:1:1:1:AAC+CCG", chr1, 10, r2R, 0, chr1, cigar0), DupFlag: false},
				{R: NewRecord("B:1:1:1:1:1:1:AAC+CCG", chr1, 10, r2R, 0, chr1, cigar0), DupFlag: false},
			},
			useUmis,
		},
		{
			// When useUmis is false and the umis differ, reads equal by position should still be duplicates
			[]TestRecord{
				{R: NewRecord("A:1:1:1:1:1:1:AAC+CCG", chr1, 0, r1F, 10, chr1, cigar0), DupFlag: false},
				{R: NewRecord("B:1:1:1:1:1:1:ACC+CCG", chr1, 0, r1F, 10, chr1, cigar0), DupFlag: true},
				{R: NewRecord("A:1:1:1:1:1:1:AAC+CCG", chr1, 10, r2R, 0, chr1, cigar0), DupFlag: false},
				{R: NewRecord("B:1:1:1:1:1:1:ACC+CCG", chr1, 10, r2R, 0, chr1, cigar0), DupFlag: true},
			},
			defaultOpts,
		},
		{
			// Mismatched r1 umis should cause non-duplicates
			[]TestRecord{
				{R: NewRecord("A:1:1:1:1:1:1:AAC+CCG", chr1, 0, r1F, 10, chr1, cigar0), DupFlag: false},
				{R: NewRecord("B:1:1:1:1:1:1:ACC+CCG", chr1, 0, r1F, 10, chr1, cigar0), DupFlag: false},
				{R: NewRecord("A:1:1:1:1:1:1:AAC+CCG", chr1, 10, r2R, 0, chr1, cigar0), DupFlag: false},
				{R: NewRecord("B:1:1:1:1:1:1:ACC+CCG", chr1, 10, r2R, 0, chr1, cigar0), DupFlag: false},
			},
			useUmis,
		},
		{
			// Mismatched r2 umis should cause non-duplicates
			[]TestRecord{
				{R: NewRecord("A:1:1:1:1:1:1:AAC+CCG", chr1, 0, r1F, 10, chr1, cigar0), DupFlag: false},
				{R: NewRecord("B:1:1:1:1:1:1:AAC+ACG", chr1, 0, r1F, 10, chr1, cigar0), DupFlag: false},
				{R: NewRecord("A:1:1:1:1:1:1:AAC+CCG", chr1, 10, r2R, 0, chr1, cigar0), DupFlag: false},
				{R: NewRecord("B:1:1:1:1:1:1:AAC+ACG", chr1, 10, r2R, 0, chr1, cigar0), DupFlag: false},
			},
			useUmis,
		},
		{
			// Equal everything with matching Ns in umi should not be duplicates
			[]TestRecord{
				{R: NewRecord("A:1:1:1:1:1:1:AAC+CCN", chr1, 0, r1F, 10, chr1, cigar0), DupFlag: false},
				{R: NewRecord("B:1:1:1:1:1:1:AAC+CCN", chr1, 0, r1F, 10, chr1, cigar0), DupFlag: false},
				{R: NewRecord("A:1:1:1:1:1:1:AAC+CCN", chr1, 10, r2R, 0, chr1, cigar0), DupFlag: false},
				{R: NewRecord("B:1:1:1:1:1:1:AAC+CCN", chr1, 10, r2R, 0, chr1, cigar0), DupFlag: false},
			},
			useUmis,
		},
		{
			// Equal everything with one Ns in r1 umi should not be duplicates
			[]TestRecord{
				{R: NewRecord("A:1:1:1:1:1:1:NAC+CCG", chr1, 0, r1F, 10, chr1, cigar0), DupFlag: false},
				{R: NewRecord("B:1:1:1:1:1:1:AAC+CCG", chr1, 0, r1F, 10, chr1, cigar0), DupFlag: false},
				{R: NewRecord("A:1:1:1:1:1:1:NAC+CCG", chr1, 10, r2R, 0, chr1, cigar0), DupFlag: false},
				{R: NewRecord("B:1:1:1:1:1:1:AAC+CCG", chr1, 10, r2R, 0, chr1, cigar0), DupFlag: false},
			},
			useUmis,
		},
		{
			// Equal everything with one Ns in r2 umi should not be duplicates
			[]TestRecord{
				{R: NewRecord("A:1:1:1:1:1:1:AAC+CCG", chr1, 0, r1F, 10, chr1, cigar0), DupFlag: false},
				{R: NewRecord("B:1:1:1:1:1:1:AAC+CCN", chr1, 0, r1F, 10, chr1, cigar0), DupFlag: false},
				{R: NewRecord("A:1:1:1:1:1:1:AAC+CCG", chr1, 10, r2R, 0, chr1, cigar0), DupFlag: false},
				{R: NewRecord("B:1:1:1:1:1:1:AAC+CCN", chr1, 10, r2R, 0, chr1, cigar0), DupFlag: false},
			},
			useUmis,
		},
		{
			// duplicate read pair with R1 vs R2 swapped, and umi1 vs umi2 swapped, should be equal
			[]TestRecord{
				{R: NewRecord("A:1:1:1:1:1:1:AAC+CCG", chr1, 0, r1F, 10, chr1, cigar0), DupFlag: false},
				{R: NewRecord("B:1:1:1:1:1:1:CCG+AAC", chr1, 0, r2F, 10, chr1, cigar0), DupFlag: true},
				{R: NewRecord("A:1:1:1:1:1:1:AAC+CCG", chr1, 10, r2R, 0, chr1, cigar0), DupFlag: false},
				{R: NewRecord("B:1:1:1:1:1:1:CCG+AAC", chr1, 10, r1R, 0, chr1, cigar0), DupFlag: true},
			},
			useUmis,
		},
		{
			// duplicate read pair with R1 vs R2 swapped, and umi1 vs umi2 swapped, should be equal
			[]TestRecord{
				{R: NewRecord("A:1:1:1:1:1:1:AAC+CCG", chr1, 0, r1R, 9, chr1, cigar0), DupFlag: false},
				{R: NewRecord("B:1:1:1:1:1:1:AAC+CCG", chr1, 0, r1R, 9, chr1, cigar0), DupFlag: true},
				{R: NewRecord("A:1:1:1:1:1:1:AAC+CCG", chr1, 9, r2F, 0, chr1, cigar0), DupFlag: false},
				{R: NewRecord("B:1:1:1:1:1:1:AAC+CCG", chr1, 9, r2F, 0, chr1, cigar0), DupFlag: true},
			},
			useUmis,
		},
		{
			// duplicate read pair with overlapping 5' positions, with match length 1, with matching FR FR orientations should be equal
			[]TestRecord{
				{R: NewRecord("A:1:1:1:1:1:1:AAC+CCG", chr1, 10, r1F, 10, chr1, cigar1M), DupFlag: false},
				{R: NewRecord("B:1:1:1:1:1:1:AAC+CCG", chr1, 10, r1F, 10, chr1, cigar1M), DupFlag: true},
				{R: NewRecord("A:1:1:1:1:1:1:AAC+CCG", chr1, 10, r2R, 10, chr1, cigar1M), DupFlag: false},
				{R: NewRecord("B:1:1:1:1:1:1:AAC+CCG", chr1, 10, r2R, 10, chr1, cigar1M), DupFlag: true},
			},
			useUmis,
		},
		{
			// duplicate read pair with overlapping 5' positions, with match length 1, with matching RF RF orientations should be equal
			[]TestRecord{
				{R: NewRecord("A:1:1:1:1:1:1:AAC+CCG", chr1, 10, r1R, 10, chr1, cigar1M), DupFlag: false},
				{R: NewRecord("B:1:1:1:1:1:1:AAC+CCG", chr1, 10, r1R, 10, chr1, cigar1M), DupFlag: true},
				{R: NewRecord("A:1:1:1:1:1:1:AAC+CCG", chr1, 10, r2F, 10, chr1, cigar1M), DupFlag: false},
				{R: NewRecord("B:1:1:1:1:1:1:AAC+CCG", chr1, 10, r2F, 10, chr1, cigar1M), DupFlag: true},
			},
			useUmis,
		},
		{
			// duplicate read pair with overlapping 5' positions, with match length 1, with matching FF FF orientations should be equal
			[]TestRecord{
				{R: NewRecord("A:1:1:1:1:1:1:AAC+CCG", chr1, 10, r1F, 10, chr1, cigar1M), DupFlag: false},
				{R: NewRecord("B:1:1:1:1:1:1:AAC+CCG", chr1, 10, r1F, 10, chr1, cigar1M), DupFlag: true},
				{R: NewRecord("A:1:1:1:1:1:1:AAC+CCG", chr1, 10, r2F, 10, chr1, cigar1M), DupFlag: false},
				{R: NewRecord("B:1:1:1:1:1:1:AAC+CCG", chr1, 10, r2F, 10, chr1, cigar1M), DupFlag: true},
			},
			useUmis,
		},
		{
			// duplicate read pair with overlapping 5' positions, with match length 1, with matching RR RR orientations should be equal
			[]TestRecord{
				{R: NewRecord("A:1:1:1:1:1:1:AAC+CCG", chr1, 10, r1R, 10, chr1, cigar1M), DupFlag: false},
				{R: NewRecord("B:1:1:1:1:1:1:AAC+CCG", chr1, 10, r1R, 10, chr1, cigar1M), DupFlag: true},
				{R: NewRecord("A:1:1:1:1:1:1:AAC+CCG", chr1, 10, r2R, 10, chr1, cigar1M), DupFlag: false},
				{R: NewRecord("B:1:1:1:1:1:1:AAC+CCG", chr1, 10, r2R, 10, chr1, cigar1M), DupFlag: true},
			},
			useUmis,
		},
		{
			// duplicate read pair with overlapping 5' positions, with matching FF FF orientations, and swapped umis should be equal
			[]TestRecord{
				{R: NewRecord("A:1:1:1:1:1:1:AAC+CCG", chr1, 10, r1F, 10, chr1, cigar1M), DupFlag: false},
				{R: NewRecord("B:1:1:1:1:1:1:CCG+AAC", chr1, 10, r1F, 10, chr1, cigar1M), DupFlag: true},
				{R: NewRecord("A:1:1:1:1:1:1:AAC+CCG", chr1, 10, r2F, 10, chr1, cigar1M), DupFlag: false},
				{R: NewRecord("B:1:1:1:1:1:1:CCG+AAC", chr1, 10, r2F, 10, chr1, cigar1M), DupFlag: true},
			},
			useUmis,
		},
		{
			// duplicate read pair with overlapping 5' positions, with matching RR RR orientations, and swapped umis should be equal
			[]TestRecord{
				{R: NewRecord("A:1:1:1:1:1:1:AAC+CCG", chr1, 10, r1R, 10, chr1, cigar1M), DupFlag: false},
				{R: NewRecord("B:1:1:1:1:1:1:CCG+AAC", chr1, 10, r1R, 10, chr1, cigar1M), DupFlag: true},
				{R: NewRecord("A:1:1:1:1:1:1:AAC+CCG", chr1, 10, r2R, 10, chr1, cigar1M), DupFlag: false},
				{R: NewRecord("B:1:1:1:1:1:1:CCG+AAC", chr1, 10, r2R, 10, chr1, cigar1M), DupFlag: true},
			},
			useUmis,
		},
		{
			// duplicate read pair with overlapping 5' positions, with match length 1, with swapped FR RF orientations and swapped
			// umis should be equal.
			[]TestRecord{
				{R: NewRecord("A:1:1:1:1:1:1:AAC+CCG", chr1, 10, r1F, 10, chr1, cigar1M), DupFlag: false},
				{R: NewRecord("B:1:1:1:1:1:1:CCG+AAC", chr1, 10, r1R, 10, chr1, cigar1M), DupFlag: true},
				{R: NewRecord("A:1:1:1:1:1:1:AAC+CCG", chr1, 10, r2R, 10, chr1, cigar1M), DupFlag: false},
				{R: NewRecord("B:1:1:1:1:1:1:CCG+AAC", chr1, 10, r2F, 10, chr1, cigar1M), DupFlag: true},
			},
			useUmis,
		},
		{
			// One read has an N in its name, should not affect umi comparison
			[]TestRecord{
				{R: NewRecord("A:1:1:1:1:1:1:AAC+CCG", chr1, 0, r1F, 10, chr1, cigar0), DupFlag: false},
				{R: NewRecord("N:1:1:1:1:1:1:AAC+CCG", chr1, 0, r1F, 10, chr1, cigar0), DupFlag: true},
				{R: NewRecord("A:1:1:1:1:1:1:AAC+CCG", chr1, 10, r2R, 0, chr1, cigar0), DupFlag: false},
				{R: NewRecord("N:1:1:1:1:1:1:AAC+CCG", chr1, 10, r2R, 0, chr1, cigar0), DupFlag: true},
			},
			useUmis,
		},
		{
			// Singleton should match against r1 from a pair.
			[]TestRecord{
				{R: NewRecord("A:1:1:1:1:1:1:AAC+CCG", chr1, 0, r1F, 10, chr1, cigar0), DupFlag: false},
				{R: NewRecord("S:1:1:1:1:1:1:AAC+GGG", chr1, 0, s1F, 10, chr1, cigar0), DupFlag: true,
					ExpectedAuxs: []sam.Aux{NewAux("DT", "LB")}},
				{R: NewRecord("A:1:1:1:1:1:1:AAC+CCG", chr1, 10, r2R, 0, chr1, cigar0), DupFlag: false},
			},
			useUmis,
		},
		{
			// Singleton should match against r2 from a pair.
			[]TestRecord{
				{R: NewRecord("A:1:1:1:1:1:1:AAC+CCG", chr1, 0, r1F, 10, chr1, cigar0), DupFlag: false},
				{R: NewRecord("S:1:1:1:1:1:1:CCC+CCG", chr1, 10, s2R, 0, chr1, cigar0), DupFlag: true,
					ExpectedAuxs: []sam.Aux{NewAux("DT", "LB")}},
				{R: NewRecord("A:1:1:1:1:1:1:AAC+CCG", chr1, 10, r2R, 0, chr1, cigar0), DupFlag: false},
			},
			useUmis,
		},
		{
			// singletons with matching N in umi should not be dupes.
			[]TestRecord{
				{R: NewRecord("A:::1:10:1:1:NAA+CCC", chr1, 0, s1F, 0, nil, cigar0), DupFlag: false},
				{R: NewRecord("A:::1:10:1:1:NAA+CCC", chr1, 0, u2, 0, nil, cigar0), DupFlag: false},
				{R: NewRecord("B:::1:10:1:1:NAA+CCC", chr1, 0, s1F, 0, nil, cigar0), DupFlag: false},
				{R: NewRecord("B:::1:10:1:1:NAA+CCC", chr1, 0, u2, 0, nil, cigar0), DupFlag: false},
			},
			useUmis,
		},
		{
			// singleton and pair with matching N in umi should not be dupes.
			[]TestRecord{
				{R: NewRecord("P:::1:10:1:1:NAA+NCC", chr1, 0, r1F, 10, chr1, cigar0), DupFlag: false},
				{R: NewRecord("A:::1:10:1:1:NAA+CCC", chr1, 0, s1F, 0, nil, cigar0), DupFlag: false},
				{R: NewRecord("A:::1:10:1:1:NAA+CCC", chr1, 0, u2, 0, nil, cigar0), DupFlag: false},
				{R: NewRecord("B:::1:10:1:1:AAA+NCC", chr1, 10, s2R, 0, nil, cigar0), DupFlag: false},
				{R: NewRecord("B:::1:10:1:1:AAA+NCC", chr1, 10, u1, 0, nil, cigar0), DupFlag: false},
				{R: NewRecord("P:::1:10:1:1:NAA+NCC", chr1, 10, r2R, 0, chr1, cigar0), DupFlag: false},
			},
			useUmis,
		},
	}
	RunTestCases(t, header, cases)
}

func TestUmiSnapCorrection(t *testing.T) {
	useUmis := defaultOpts
	useUmis.UseUmis = true

	snapCorrection := defaultOpts
	snapCorrection.UseUmis = true
	snapCorrection.KnownUmis = []byte("AAA\nCCC\nGGG\nTTT")

	cases := []TestCase{
		{
			// Snappable umis should not match against each other if umi correction is off.
			[]TestRecord{
				{R: NewRecord("A:1:1:1:1:1:1:AAA+CCC", chr1, 0, r1F, 10, chr1, cigar0), DupFlag: false},
				{R: NewRecord("B:1:1:1:1:1:1:AAC+CCG", chr1, 0, r1F, 10, chr1, cigar0), DupFlag: false},
				{R: NewRecord("A:1:1:1:1:1:1:AAA+CCC", chr1, 10, r2R, 0, chr1, cigar0), DupFlag: false},
				{R: NewRecord("B:1:1:1:1:1:1:AAC+CCG", chr1, 10, r2R, 0, chr1, cigar0), DupFlag: false},
			},
			useUmis,
		},
		{
			// Snappable umis should match against each other if umi correction is on.
			[]TestRecord{
				{R: NewRecord("A:1:1:1:1:1:1:AAA+CCC", chr1, 0, r1F, 10, chr1, cigar0),
					DupFlag: false, ExpectedAuxs: []sam.Aux{NewAux("DI", "0"), NewAux("DS", 2)}},
				{R: NewRecord("B:1:1:1:1:1:1:AAC+CCG", chr1, 0, r1F, 10, chr1, cigar0),
					DupFlag: true, ExpectedAuxs: []sam.Aux{NewAux("DI", "0"), NewAux("DS", 2), NewAux("DU", "AAA+CCC")}},
				{R: NewRecord("A:1:1:1:1:1:1:AAA+CCC", chr1, 10, r2R, 0, chr1, cigar0),
					DupFlag: false, ExpectedAuxs: []sam.Aux{NewAux("DI", "0"), NewAux("DS", 2)}},
				{R: NewRecord("B:1:1:1:1:1:1:AAC+CCG", chr1, 10, r2R, 0, chr1, cigar0),
					DupFlag: true, ExpectedAuxs: []sam.Aux{NewAux("DI", "0"), NewAux("DS", 2), NewAux("DU", "AAA+CCC")}},
			},
			snapCorrection,
		},
		{
			// N in a umi is snappable.
			[]TestRecord{
				{R: NewRecord("A:1:1:1:1:1:1:AAA+CCC", chr1, 0, r1F, 10, chr1, cigar0),
					DupFlag: false, ExpectedAuxs: []sam.Aux{NewAux("DI", "0"), NewAux("DS", 2)}},
				{R: NewRecord("B:1:1:1:1:1:1:ANA+CCC", chr1, 0, r1F, 10, chr1, cigar0),
					DupFlag: true, ExpectedAuxs: []sam.Aux{NewAux("DI", "0"), NewAux("DS", 2), NewAux("DU", "AAA+CCC")}},
				{R: NewRecord("A:1:1:1:1:1:1:AAA+CCC", chr1, 10, r2R, 0, chr1, cigar0),
					DupFlag: false, ExpectedAuxs: []sam.Aux{NewAux("DI", "0"), NewAux("DS", 2)}},
				{R: NewRecord("B:1:1:1:1:1:1:ANA+CCC", chr1, 10, r2R, 0, chr1, cigar0),
					DupFlag: true, ExpectedAuxs: []sam.Aux{NewAux("DI", "0"), NewAux("DS", 2), NewAux("DU", "AAA+CCC")}},
			},
			snapCorrection,
		},
		{
			// Some UMIs are not snappable and should not be corrected (TAG) is not snappable.
			[]TestRecord{
				{R: NewRecord("A:1:1:1:1:1:1:AAA+CCC", chr1, 0, r1F, 10, chr1, cigar0), DupFlag: false},
				{R: NewRecord("B:1:1:1:1:1:1:TAG+CCC", chr1, 0, r1F, 10, chr1, cigar0), DupFlag: false},
				{R: NewRecord("A:1:1:1:1:1:1:AAA+CCC", chr1, 10, r2R, 0, chr1, cigar0), DupFlag: false},
				{R: NewRecord("B:1:1:1:1:1:1:TAG+CCC", chr1, 10, r2R, 0, chr1, cigar0), DupFlag: false},
			},
			snapCorrection,
		},
		{
			// singleton and pair with snappable N in umi should be dupes.
			[]TestRecord{
				{R: NewRecord("P:::1:10:1:1:NAA+NCC", chr1, 0, r1F, 10, chr1, cigar0), DupFlag: false},
				{R: NewRecord("A:::1:10:1:1:NAA+CCC", chr1, 0, s1F, 0, nil, cigar0), DupFlag: true},
				{R: NewRecord("A:::1:10:1:1:NAA+CCC", chr1, 0, u2, 0, nil, cigar0), DupFlag: false},
				{R: NewRecord("B:::1:10:1:1:AAA+NCC", chr1, 10, s2R, 0, nil, cigar0), DupFlag: true},
				{R: NewRecord("B:::1:10:1:1:AAA+NCC", chr1, 10, u1, 0, nil, cigar0), DupFlag: false},
				{R: NewRecord("P:::1:10:1:1:NAA+NCC", chr1, 10, r2R, 0, chr1, cigar0), DupFlag: false},
			},
			snapCorrection,
		},
		{
			// Pair P has one snappable umi, and one non-snappable.  Singleton B should be marked as duplicate of P on snappable umi.
			[]TestRecord{
				{R: NewRecord("P:::1:10:1:1:NAC+NCC", chr1, 0, r1F, 10, chr1, cigar0), DupFlag: false},
				{R: NewRecord("P:::1:10:1:1:NAC+NCC", chr1, 10, r2R, 0, chr1, cigar0), DupFlag: false},
				{R: NewRecord("B:::1:10:1:1:AAA+CCC", chr1, 10, s2R, 0, nil, cigar0), DupFlag: true},
				{R: NewRecord("B:::1:10:1:1:AAA+CCC", chr1, 10, u1, 0, nil, cigar0), DupFlag: false},
			},
			snapCorrection,
		},
	}
	RunTestCases(t, header, cases)
}

func TestUmiScavengeCorrection(t *testing.T) {
	noScavenge := defaultOpts
	noScavenge.UseUmis = true
	noScavenge.KnownUmis = []byte("AAA\nCCC\nGGG\nTTT")
	noScavenge.ScavengeUmis = -1

	scavenge := defaultOpts
	scavenge.UseUmis = true
	scavenge.KnownUmis = []byte("AAA\nCCC\nGGG\nTTT")
	scavenge.ScavengeUmis = 2

	cases := []TestCase{
		{
			// B is not snappable to A, but should be scavengable to A.  If scavenging is off, then make sure B is not a dup.
			[]TestRecord{
				{R: NewRecord("A:1:1:1:1:1:1:AAA+CCC", chr1, 0, r1F, 10, chr1, cigar0), DupFlag: false},
				{R: NewRecord("B:1:1:1:1:1:1:TAC+CCG", chr1, 0, r1F, 10, chr1, cigar0), DupFlag: false},
				{R: NewRecord("A:1:1:1:1:1:1:AAA+CCC", chr1, 10, r2R, 0, chr1, cigar0), DupFlag: false},
				{R: NewRecord("B:1:1:1:1:1:1:TAC+CCG", chr1, 10, r2R, 0, chr1, cigar0), DupFlag: false},
			},
			noScavenge,
		},
		{
			// B is not snappable to A, but should be scavengable to A.
			[]TestRecord{
				{R: NewRecord("A:1:1:1:1:1:1:AAA+CCC", chr1, 0, r1F, 10, chr1, cigar0), DupFlag: false},
				{R: NewRecord("B:1:1:1:1:1:1:TAC+CCG", chr1, 0, r1F, 10, chr1, cigar0), DupFlag: true,
					ExpectedAuxs: []sam.Aux{NewAux("DU", "AAA+CCC")}},
				{R: NewRecord("A:1:1:1:1:1:1:AAA+CCC", chr1, 10, r2R, 0, chr1, cigar0), DupFlag: false},
				{R: NewRecord("B:1:1:1:1:1:1:TAC+CCG", chr1, 10, r2R, 0, chr1, cigar0), DupFlag: true,
					ExpectedAuxs: []sam.Aux{NewAux("DU", "AAA+CCC")}},
			},
			scavenge,
		},
		{
			// B is not snappable to A, but should be scavengable to A if umis are swapped.
			[]TestRecord{
				{R: NewRecord("A:1:1:1:1:1:1:AAA+CCC", chr1, 10, r1F, 10, chr1, cigar1M), DupFlag: false},
				{R: NewRecord("B:1:1:1:1:1:1:CCG+TAC", chr1, 10, r1R, 10, chr1, cigar1M), DupFlag: true,
					ExpectedAuxs: []sam.Aux{NewAux("DU", "CCC+AAA")}},
				{R: NewRecord("A:1:1:1:1:1:1:AAA+CCC", chr1, 10, r2R, 10, chr1, cigar1M), DupFlag: false},
				{R: NewRecord("B:1:1:1:1:1:1:CCG+TAC", chr1, 10, r2F, 10, chr1, cigar1M), DupFlag: true,
					ExpectedAuxs: []sam.Aux{NewAux("DU", "CCC+AAA")}},
			},
			scavenge,
		},
		{
			// A and B are singletons.  B is not snappable to A, but B is scavengable to A.
			[]TestRecord{
				{R: NewRecord("A:::1:10:1:1:AAC+GGG", chr1, 0, s1F, 0, nil, cigar0), DupFlag: false},
				{R: NewRecord("A:::1:10:1:1:AAC+GGG", chr1, 0, u2, 0, nil, cigar0), DupFlag: false},
				{R: NewRecord("B:::1:10:1:1:GAC+GGG", chr1, 0, s1F, 0, nil, cigar0), DupFlag: true},
				{R: NewRecord("B:::1:10:1:1:GAC+GGG ", chr1, 0, u2, 0, nil, cigar0), DupFlag: false},
			},
			scavenge,
		},
	}
	RunTestCases(t, header, cases)
}

func TestSeparateSingletons(t *testing.T) {
	separateSingletons := defaultOpts
	separateSingletons.SeparateSingletons = true

	cases := []TestCase{
		{
			// Singleton should match against r1 from a pair.
			[]TestRecord{
				{R: NewRecord("A:1:1:1:1:1:1", chr1, 0, r1F, 10, chr1, cigar0), DupFlag: false},
				{R: NewRecord("S:1:1:1:1:1:1", chr1, 0, s1F, 10, chr1, cigar0), DupFlag: true,
					ExpectedAuxs: []sam.Aux{NewAux("DT", "LB")}},
				{R: NewRecord("A:1:1:1:1:1:1", chr1, 10, r2R, 0, chr1, cigar0), DupFlag: false},
			},
			defaultOpts,
		},
		{
			// Singleton should match against r2 from a pair.
			[]TestRecord{
				{R: NewRecord("A:1:1:1:1:1:1", chr1, 0, r1F, 10, chr1, cigar0), DupFlag: false},
				{R: NewRecord("S:1:1:1:1:1:1", chr1, 10, s2R, 0, chr1, cigar0), DupFlag: true,
					ExpectedAuxs: []sam.Aux{NewAux("DT", "LB")}},
				{R: NewRecord("A:1:1:1:1:1:1", chr1, 10, r2R, 0, chr1, cigar0), DupFlag: false},
			},
			defaultOpts,
		},
		{
			// Singleton should not match against r1 from a pair if separateSingletons=true.
			[]TestRecord{
				{R: NewRecord("A:1:1:1:1:1:1", chr1, 0, r1F, 10, chr1, cigar0), DupFlag: false},
				{R: NewRecord("S:1:1:1:1:1:1", chr1, 0, s1F, 10, chr1, cigar0), DupFlag: false},
				{R: NewRecord("A:1:1:1:1:1:1", chr1, 10, r2R, 0, chr1, cigar0), DupFlag: false},
			},
			separateSingletons,
		},
		{
			// Singleton should not match against r2 from a pair if separateSingletons=true.
			[]TestRecord{
				{R: NewRecord("A:1:1:1:1:1:1", chr1, 0, r1F, 10, chr1, cigar0), DupFlag: false},
				{R: NewRecord("S:1:1:1:1:1:1", chr1, 10, s2R, 0, chr1, cigar0), DupFlag: false},
				{R: NewRecord("A:1:1:1:1:1:1", chr1, 10, r2R, 0, chr1, cigar0), DupFlag: false},
			},
			separateSingletons,
		},
	}
	RunTestCases(t, header, cases)
}

// Ensure that int-di mode correctly formats DI aux tag as 'i' integer.
func TestIntDI(t *testing.T) {
	tempDir, cleanup := testutil.TempDir(t, "", "")
	defer cleanup()

	testrecords := []*sam.Record{
		basicA1, basicB1, basicA2, basicB2,
	}

	for _, format := range []string{"bam", "pam"} {
		provider := bamprovider.NewFakeProvider(header, testrecords)
		outputPath := NewTestOutput(tempDir, 0, format)
		opts := Opts{
			ShardSize:            100,
			Padding:              10,
			Parallelism:          1,
			QueueLength:          10,
			ClearExisting:        false,
			RemoveDups:           false,
			TagDups:              true,
			IntDI:                true,
			EmitUnmodifiedFields: true,
			OutputPath:           outputPath,
			Format:               format,
			OpticalDetector: &TileOpticalDetector{
				OpticalDistance: 2500,
			},
		}
		markDuplicates := &MarkDuplicates{
			Provider: provider,
			Opts:     &opts,
		}
		_, err := markDuplicates.Mark(nil)
		assert.NoError(t, err)

		actualRecords := ReadRecords(t, outputPath)
		assert.Equal(t, len(testrecords), len(actualRecords))
		for i, r := range actualRecords {
			t.Logf("output[%v]: %v", i, r)

			// Verify that DI tag exist, and have the right formatting and value.
			expectedAux := NewAux("DI", 0)
			actual, ok := r.Tag([]byte{expectedAux.Tag()[0], expectedAux.Tag()[1]})
			assert.Equal(t, true, ok, "Expected tag %s to exist, but it does not", expectedAux)
			if ok {
				assert.Equal(t, expectedAux, actual)
			}
		}
	}
}

// Verify that optical distance histogram is correct.
func TestOpticalHistogram(t *testing.T) {
	tests := []struct {
		records      []*sam.Record
		expectedHist []map[int]int // expected counts in opticalDistance bucket
	}{
		{
			[]*sam.Record{
				NewRecord("oA:::1:10:1:1", chr1, 0, r1F, 100, chr1, cigar0),
				NewRecord("oB:::1:10:1:5", chr1, 0, r1F, 100, chr1, cigar0),
				NewRecord("oA:::1:10:1:1", chr1, 100, r2R, 0, chr1, cigar0),
				NewRecord("oB:::1:10:1:5", chr1, 100, r2R, 0, chr1, cigar0),
			},
			[]map[int]int{
				map[int]int{4: 1},
				map[int]int{},
				map[int]int{},
				map[int]int{},
			},
		},
		{
			[]*sam.Record{
				// 3-4-5 triangle.
				NewRecord("oA:::1:10:1:1", chr1, 0, r1F, 100, chr1, cigar0),
				NewRecord("oB:::1:10:1:4", chr1, 0, r1F, 100, chr1, cigar0),
				NewRecord("oC:::1:10:5:1", chr1, 0, r1F, 100, chr1, cigar0),
				NewRecord("oA:::1:10:1:1", chr1, 100, r2R, 0, chr1, cigar0),
				NewRecord("oB:::1:10:1:4", chr1, 100, r2R, 0, chr1, cigar0),
				NewRecord("oC:::1:10:5:1", chr1, 100, r2R, 0, chr1, cigar0),
			},
			[]map[int]int{
				map[int]int{},
				map[int]int{
					3: 1,
					4: 1,
					5: 1,
				},
				map[int]int{},
				map[int]int{},
			},
		},
		{
			[]*sam.Record{
				NewRecord("oA:::1:10:1:1", chr1, 0, r1F, 100, chr1, cigar0),
				NewRecord("oB:::1:10:1:2", chr1, 0, r1F, 100, chr1, cigar0),
				NewRecord("oC:::1:10:1:3", chr1, 0, r1F, 100, chr1, cigar0),
				NewRecord("oD:::1:10:1:4", chr1, 0, r1F, 100, chr1, cigar0),
				NewRecord("oE:::1:10:1:5", chr1, 0, r1F, 100, chr1, cigar0),

				NewRecord("oA:::1:10:1:1", chr1, 100, r2R, 0, chr1, cigar0),
				NewRecord("oB:::1:10:1:2", chr1, 100, r2R, 0, chr1, cigar0),
				NewRecord("oC:::1:10:1:3", chr1, 100, r2R, 0, chr1, cigar0),
				NewRecord("oD:::1:10:1:4", chr1, 100, r2R, 0, chr1, cigar0),
				NewRecord("oE:::1:10:1:5", chr1, 100, r2R, 0, chr1, cigar0),
			},
			[]map[int]int{
				map[int]int{},
				map[int]int{},
				map[int]int{
					1: 4,
					2: 3,
					3: 2,
					4: 1,
				},
				map[int]int{},
			},
		},
		{
			[]*sam.Record{
				NewRecord("oA:::1:10:1:1", chr1, 0, r1F, 100, chr1, cigar0),
				NewRecord("oB:::1:10:1:2", chr1, 0, r1F, 100, chr1, cigar0),
				NewRecord("oC:::1:10:1:3", chr1, 0, r1F, 100, chr1, cigar0),
				NewRecord("oD:::1:10:1:4", chr1, 0, r1F, 100, chr1, cigar0),
				NewRecord("oE:::1:10:1:5", chr1, 0, r1F, 100, chr1, cigar0),
				NewRecord("oF:::1:10:1:6", chr1, 0, r1F, 100, chr1, cigar0),
				NewRecord("oG:::1:10:1:7", chr1, 0, r1F, 100, chr1, cigar0),
				NewRecord("oH:::1:10:1:8", chr1, 0, r1F, 100, chr1, cigar0),

				NewRecord("oA:::1:10:1:1", chr1, 100, r2R, 0, chr1, cigar0),
				NewRecord("oB:::1:10:1:2", chr1, 100, r2R, 0, chr1, cigar0),
				NewRecord("oC:::1:10:1:3", chr1, 100, r2R, 0, chr1, cigar0),
				NewRecord("oD:::1:10:1:4", chr1, 100, r2R, 0, chr1, cigar0),
				NewRecord("oE:::1:10:1:5", chr1, 100, r2R, 0, chr1, cigar0),
				NewRecord("oF:::1:10:1:6", chr1, 100, r2R, 0, chr1, cigar0),
				NewRecord("oG:::1:10:1:7", chr1, 100, r2R, 0, chr1, cigar0),
				NewRecord("oH:::1:10:1:8", chr1, 100, r2R, 0, chr1, cigar0),
			},
			[]map[int]int{
				map[int]int{},
				map[int]int{},
				map[int]int{},
				map[int]int{
					1: 7,
					2: 6,
					3: 5,
					4: 4,
					5: 3,
					6: 2,
					7: 1,
				},
			},
		},
	}

	tempDir, cleanup := testutil.TempDir(t, "", "")
	defer cleanup()
	for testIdx, test := range tests {
		for _, format := range []string{"bam", "pam"} {
			t.Logf("---- starting tests[%d] ----", testIdx)
			provider := bamprovider.NewFakeProvider(header, test.records)
			outputPath := NewTestOutput(tempDir, testIdx, format)
			opts := defaultOpts
			opts.OutputPath = outputPath
			opts.Format = format
			opts.OpticalHistogram = "optical-histogram.txt"
			opts.OpticalHistogramMax = -1

			markDuplicates := &MarkDuplicates{
				Provider: provider,
				Opts:     &opts,
			}
			actualMetrics, err := markDuplicates.Mark(nil)
			assert.NoError(t, err)

			t.Logf("distances: %v", actualMetrics.OpticalDistance)
			assert.Equal(t, len(test.expectedHist), len(actualMetrics.OpticalDistance))
			for i := range test.expectedHist {
				for j := range actualMetrics.OpticalDistance[i] {
					assert.Equal(t, int64(test.expectedHist[i][j]), int64(actualMetrics.OpticalDistance[i][j]),
						"i %d j %d %v", i, j, actualMetrics.OpticalDistance)
				}
			}
		}
	}
}

func TestOpticalHistogramMax(t *testing.T) {
	const max = 1000
	var records []*sam.Record

	// All records are in a 3,4,5 triangle.
	for i := 0; i < max; i++ {
		records = append(records, []*sam.Record{
			NewRecord(fmt.Sprintf("A%d:::1:10:1:1", i), chr1, 0, r1F, 100, chr1, cigar0),
			NewRecord(fmt.Sprintf("B%d:::1:10:1:4", i), chr1, 0, r1F, 100, chr1, cigar0),
			NewRecord(fmt.Sprintf("C%d:::1:10:5:1", i), chr1, 0, r1F, 100, chr1, cigar0),
		}...)
	}
	for i := 0; i < max; i++ {
		records = append(records, []*sam.Record{
			NewRecord(fmt.Sprintf("A%d:::1:10:1:1", i), chr1, 100, r2R, 0, chr1, cigar0),
			NewRecord(fmt.Sprintf("B%d:::1:10:1:4", i), chr1, 100, r2R, 0, chr1, cigar0),
			NewRecord(fmt.Sprintf("C%d:::1:10:5:1", i), chr1, 100, r2R, 0, chr1, cigar0),
		}...)
	}

	tempDir, cleanup := testutil.TempDir(t, "", "")
	defer cleanup()
	provider := bamprovider.NewFakeProvider(header, records)
	outputPath := NewTestOutput(tempDir, 0, "bam")
	opts := defaultOpts
	opts.OutputPath = outputPath
	opts.Format = "bam"
	opts.OpticalHistogram = "optical-histogram.txt"
	opts.OpticalHistogramMax = 300

	markDuplicates := &MarkDuplicates{
		Provider: provider,
		Opts:     &opts,
	}
	actualMetrics, err := markDuplicates.Mark(nil)
	assert.NoError(t, err)

	// opts.OpticalHistogramMax is set to 300, that means 100 each of
	// A, B, and C.  That implies that each pair A-B, B-C, and C-A
	// should have 100*100 histogram entries.
	assert.Equal(t, 4, len(actualMetrics.OpticalDistance))
	assert.True(t, int64(10000*.9) < actualMetrics.OpticalDistance[3][3] && actualMetrics.OpticalDistance[3][3] < int64(10000*1.1),
		fmt.Sprintf("%d is out of expected range (%d, %d)", actualMetrics.OpticalDistance[3][3], int64(10000*.9), int64(10000*1.1)))
	assert.True(t, int64(10000*.9) < actualMetrics.OpticalDistance[3][4] && actualMetrics.OpticalDistance[3][4] < int64(10000*1.1),
		fmt.Sprintf("%d is out of expected range (%d, %d)", actualMetrics.OpticalDistance[3][4], int64(10000*.9), int64(10000*1.1)))
	assert.True(t, int64(10000*.9) < actualMetrics.OpticalDistance[3][5] && actualMetrics.OpticalDistance[3][5] < int64(10000*1.1),
		fmt.Sprintf("%d is out of expected range (%d, %d)", actualMetrics.OpticalDistance[3][5], int64(10000*.9), int64(10000*1.1)))
}

func TestStrandSpecific(t *testing.T) {
	notStrandSpecific := defaultOpts
	strandSpecific := defaultOpts
	strandSpecific.StrandSpecific = true

	cases := []TestCase{
		{
			// A and B are from different strands, if strandSpecific =
			// false, they should be duplicates.
			[]TestRecord{
				{R: NewRecord("A:1:1:1:1:1:1:AAC+CCG", chr1, 0, r1F, 10, chr1, cigar0), DupFlag: false},
				{R: NewRecord("B:1:1:1:1:1:2:CCG+AAC", chr1, 0, r2F, 10, chr1, cigar0), DupFlag: true},
				{R: NewRecord("A:1:1:1:1:1:1:AAC+CCG", chr1, 10, r2R, 0, chr1, cigar0), DupFlag: false},
				{R: NewRecord("B:1:1:1:1:1:2:CCG+AAC", chr1, 10, r1R, 0, chr1, cigar0), DupFlag: true},
			},
			notStrandSpecific,
		},
		{
			// A and B are from different strand, if strandSpecific =
			// true, they should not be duplicates.
			[]TestRecord{
				{R: NewRecord("A:1:1:1:1:1:1:AAC+CCG", chr1, 0, r1F, 10, chr1, cigar0), DupFlag: false},
				{R: NewRecord("B:1:1:1:1:1:2:CCG+AAC", chr1, 0, r2F, 10, chr1, cigar0), DupFlag: false},
				{R: NewRecord("A:1:1:1:1:1:1:AAC+CCG", chr1, 10, r2R, 0, chr1, cigar0), DupFlag: false},
				{R: NewRecord("B:1:1:1:1:1:2:CCG+AAC", chr1, 10, r1F, 0, chr1, cigar0), DupFlag: false},
			},
			strandSpecific,
		},
		{
			// A and B are from the same strand, if strandSpecific =
			// true, they should be duplicates.
			[]TestRecord{
				{R: NewRecord("A:1:1:1:1:1:1:AAC+CCG", chr1, 0, r1F, 10, chr1, cigar0), DupFlag: false},
				{R: NewRecord("B:1:1:1:1:1:2:AAC+CCG", chr1, 0, r1F, 10, chr1, cigar0), DupFlag: true},
				{R: NewRecord("A:1:1:1:1:1:1:AAC+CCG", chr1, 10, r2R, 0, chr1, cigar0), DupFlag: false},
				{R: NewRecord("B:1:1:1:1:1:2:AAC+CCG", chr1, 10, r2R, 0, chr1, cigar0), DupFlag: true},
			},
			strandSpecific,
		},
		{
			// Singleton S has a position match and a matching
			// strand. If strandSpecific = true, they should be
			// duplicates.
			[]TestRecord{
				{R: NewRecord("A:1:1:1:1:1:1:AAC+CCG", chr1, 0, r1F|sam.MateReverse, 10, chr1, cigar0), DupFlag: false},
				{R: NewRecord("A:1:1:1:1:1:1:AAC+CCG", chr1, 10, r2R, 0, chr1, cigar0), DupFlag: false},
				{R: NewRecord("S:1:1:1:1:1:1:AAC+CCG", chr1, 10, s2R, 10, chr1, cigar0), DupFlag: true},
			},
			strandSpecific,
		},
		{
			// Singleton S has a position match but not a matching
			// strand. If strandSpecific = true, they should not be
			// duplicates.
			[]TestRecord{
				{R: NewRecord("A:1:1:1:1:1:1:AAC+CCG", chr1, 0, r1R|sam.MateReverse, 10, chr1, cigar0), DupFlag: false},
				{R: NewRecord("A:1:1:1:1:1:1:AAC+CCG", chr1, 10, r2R|sam.MateReverse, 0, chr1, cigar0), DupFlag: false},
				{R: NewRecord("S:1:1:1:1:1:1:AAC+CCG", chr1, 10, s2R, 10, chr1, cigar0), DupFlag: false},
			},
			strandSpecific,
		},
	}
	for _, useUmis := range []bool{false, true} {
		t.Logf("useUmis: %v", useUmis)
		for i := range cases {
			cases[i].Opts.UseUmis = useUmis
		}
		RunTestCases(t, header, cases)
	}
}

// Test that BagIDs match when 1 read is in a shard that crosses
// reference boundary, and there are records with a alignment less
// than the shard start's alignment position in the second reference
// in the shard.
func TestBagID(t *testing.T) {
	tempDir, cleanup := testutil.TempDir(t, "", "")
	defer cleanup()

	shards := []gbam.Shard{
		gbam.Shard{
			StartRef: chr1,
			EndRef:   chr1,
			Start:    0,
			End:      100,
			StartSeq: 0,
			EndSeq:   0,
			Padding:  10,
			ShardIdx: 0,
		},
		gbam.Shard{
			StartRef: chr1,
			EndRef:   chr2, // Extends past chr1 into chr2.
			Start:    100,
			End:      100,
			StartSeq: 0,
			EndSeq:   0,
			Padding:  10,
			ShardIdx: 1,
		},
		gbam.Shard{
			StartRef: chr2,
			EndRef:   chr2, // Extends past chr1 into chr2.
			Start:    100,
			End:      2000,
			StartSeq: 0,
			EndSeq:   0,
			Padding:  10,
			ShardIdx: 2,
		},
		gbam.Shard{
			StartRef: nil,
			EndRef:   nil,
			Start:    0,
			End:      0,
			StartSeq: 0,
			EndSeq:   0,
			Padding:  10,
			ShardIdx: 3,
		},
	}

	testrecords := []*sam.Record{
		NewRecord("A:1:1:1:1:1:1:AAC+CCG", chr1, 200, r1F|sam.MateReverse, 200, chr2, cigar0),
		NewRecord("B:1:1:1:1:1:1:AAC+CCG", chr1, 200, r1F|sam.MateReverse, 200, chr2, cigar0),

		// We need Q to reside in shard1 but in the beginning of the second reference.
		NewRecord("Q:1:1:1:1:1:1:AAC+CCG", chr2, 50, r1F, 52, chr2, cigar0),
		NewRecord("Q:1:1:1:1:1:1:AAC+CCG", chr2, 52, r2R, 50, chr2, cigar0),

		NewRecord("A:1:1:1:1:1:1:AAC+CCG", chr2, 200, r2R, 200, chr1, cigar0),
		NewRecord("B:1:1:1:1:1:1:AAC+CCG", chr2, 200, r2R, 200, chr1, cigar0),
	}

	for _, format := range []string{"bam", "pam"} {
		provider := bamprovider.NewFakeProvider(header, testrecords)
		outputPath := NewTestOutput(tempDir, 0, format)
		opts := Opts{
			Padding:              10,
			Parallelism:          1,
			QueueLength:          10,
			ClearExisting:        false,
			RemoveDups:           false,
			TagDups:              true,
			IntDI:                true,
			EmitUnmodifiedFields: true,
			OutputPath:           outputPath,
			Format:               format,
			OpticalDetector: &TileOpticalDetector{
				OpticalDistance: 2500,
			},
		}
		markDuplicates := &MarkDuplicates{
			Provider: provider,
			Opts:     &opts,
		}
		_, err := markDuplicates.Mark(shards)
		assert.NoError(t, err)

		actualRecords := ReadRecords(t, outputPath)
		assert.Equal(t, len(testrecords), len(actualRecords))
		var commonDI []byte
		for i, r := range actualRecords {
			t.Logf("output[%v]: %v", i, r)
			if strings.HasPrefix(r.Name, "Q") {
				continue
			}

			// Verify that DI tag exist, and have the right value.
			expectedAux := NewAux("DI", 0)
			actual, ok := r.Tag([]byte{expectedAux.Tag()[0], expectedAux.Tag()[1]})
			assert.True(t, ok)
			assert.NotNil(t, actual)
			if commonDI == nil {
				commonDI = actual
			} else {
				assert.True(t, bytes.Equal(commonDI, actual), "bytes %v %v", commonDI, actual)
			}
		}
	}
}

func TestOpticalDetector(t *testing.T) {
	tests := []struct {
		records         []*sam.Record
		opticalDistance int
		metrics         *MetricsCollection
		dupType         []string
	}{
		{
			//  Read1 is in shard, but Read2 is distant in a different reference (but in the same numerical shard range as Read1).
			//    Explanation:
			//           Ref  Pos Shard
			//    Read1  chr1 50  chr1:0-100
			//    Read2  chr2 55  chr2:0-100 -- Note that Read2 is between 0 and 100,
			//                                  and could be mistaken to be inside shard chr1:0-100
			//                                  if we forget to compare Read2's reference which is
			//                                  chr2 when checking if Read2 is in shard chr1:0-100.
			[]*sam.Record{
				NewRecord("M:::1:10:6:6", chr1, 50, r1F, 55, chr2, cigar0),
				NewRecord("N:::1:10:6:6", chr1, 50, r1F, 55, chr2, cigar0),
				NewRecord("M:::1:10:6:6", chr2, 55, r2F, 50, chr1, cigar0),
				NewRecord("N:::1:10:6:6", chr2, 55, r2F, 50, chr1, cigar0),
			},
			2500,
			&MetricsCollection{
				// Notes that ReadPairsExamined, ReadPairDups, and
				// ReadPairOpticalDups are doubled here because they
				// are halved when written to the metrics file.
				LibraryMetrics: map[string]*Metrics{
					"Unknown Library": &Metrics{
						ReadPairsExamined:   4,
						ReadPairDups:        2,
						ReadPairOpticalDups: 2,
					},
				},
			},
			[]string{"", "SQ", "", "SQ"},
		},
		{
			//  A and B are optical duplicates with threshold distance 2500.
			[]*sam.Record{
				NewRecord("oA:::1:10:1:1", chr1, 0, r1F, 100, chr1, cigar0),
				NewRecord("oB:::1:10:5:5", chr1, 0, r1F, 100, chr1, cigar0),
				NewRecord("oA:::1:10:1:1", chr1, 100, r2R, 0, chr1, cigar0),
				NewRecord("oB:::1:10:5:5", chr1, 100, r2R, 0, chr1, cigar0),
			},
			2500,
			&MetricsCollection{
				LibraryMetrics: map[string]*Metrics{
					"Unknown Library": &Metrics{
						ReadPairsExamined:   4,
						ReadPairDups:        2,
						ReadPairOpticalDups: 2,
					},
				},
			},
			[]string{"", "SQ", "", "SQ"},
		},
		{
			// C is not an optical duplicate of A because 3000 is too far from 10.
			[]*sam.Record{
				NewRecord("oA:::1:10:1:1", chr1, 0, r1F, 100, chr1, cigar0),
				NewRecord("oC:::1:10:3000:5", chr1, 0, r1F, 100, chr1, cigar0),
				NewRecord("oA:::1:10:1:1", chr1, 100, r2R, 0, chr1, cigar0),
				NewRecord("oC:::1:10:3000:5", chr1, 100, r2R, 0, chr1, cigar0),
			},
			2500,
			&MetricsCollection{
				LibraryMetrics: map[string]*Metrics{
					"Unknown Library": &Metrics{
						ReadPairsExamined:   4,
						ReadPairDups:        2,
						ReadPairOpticalDups: 0,
					},
				},
			},
			[]string{"", "LB", "", "LB"},
		},
		{
			// D is a duplicate, but not an *optical* duplicate of A because the Read1/Read2 orientations do not match do not match.
			[]*sam.Record{
				NewRecord("oA:::1:10:1:1", chr1, 0, r1F, 100, chr1, cigar0),
				NewRecord("oD:::1:10:5:5", chr1, 0, r2F, 100, chr1, cigar0),
				NewRecord("oA:::1:10:1:1", chr1, 100, r2R, 0, chr1, cigar0),
				NewRecord("oD:::1:10:5:5", chr1, 100, r1R, 0, chr1, cigar0),
			},
			2500,
			&MetricsCollection{
				LibraryMetrics: map[string]*Metrics{
					"Unknown Library": &Metrics{
						ReadPairsExamined:   4,
						ReadPairDups:        2,
						ReadPairOpticalDups: 0,
					},
				},
			},
			[]string{"", "LB", "", "LB"},
		},
		{
			//  G has alignment position in shard, and H has alignment position in padding, but they are still duplicates.
			[]*sam.Record{
				NewRecord("G:::1:10:5:5", chr1, 99, r1F, 115, chr1, cigar0),
				NewRecord("H:::1:10:6:6", chr1, 101, r1F, 115, chr1, cigarSoft2),
				NewRecord("G:::1:10:5:5", chr1, 115, r2F, 99, chr1, cigar0),
				NewRecord("H:::1:10:6:6", chr1, 115, r2F, 101, chr1, cigar0),
			},
			2500,
			&MetricsCollection{
				LibraryMetrics: map[string]*Metrics{
					"Unknown Library": &Metrics{
						ReadPairsExamined:   4,
						ReadPairDups:        2,
						ReadPairOpticalDups: 2,
					},
				},
			},
			[]string{"", "SQ", "", "SQ"},
		},
		{
			// E and F are optical duplicates, and both ends are in a padding region.  Make sure each optical duplicate is
			// counted just once.
			[]*sam.Record{
				NewRecord("oE:::1:10:1:1", chr1, 103, r1F, 203, chr1, cigar0),
				NewRecord("oF:::1:10:5:5", chr1, 103, r1F, 203, chr1, cigar0),
				NewRecord("oE:::1:10:1:1", chr1, 203, r2R, 103, chr1, cigar0),
				NewRecord("oF:::1:10:5:5", chr1, 203, r2R, 103, chr1, cigar0),
			},
			2500,
			&MetricsCollection{
				LibraryMetrics: map[string]*Metrics{
					"Unknown Library": &Metrics{
						ReadPairsExamined:   4,
						ReadPairDups:        2,
						ReadPairOpticalDups: 2,
					},
				},
			},
			[]string{"", "SQ", "", "SQ"},
		},
		{
			// G is best, only I should be optical duplicate, because it shows up first in the file.  MarkDuplicates must sort the
			// duplicate set before determining optical duplicates,
			[]*sam.Record{
				NewRecord("oG:::1:10:1:1", chr1, 50, r1F, 150, chr1, cigar0),
				NewRecord("oH:::1:10:3000:5", chr1, 50, r1F, 151, chr1, cigar0),
				NewRecord("oI:::1:10:4000:5", chr1, 51, r1F, 151, chr1, cigarSoft1),
				NewRecord("oI:::1:10:4000:5", chr1, 150, r2R, 51, chr1, cigar0),
				NewRecord("oG:::1:10:1:1", chr1, 150, r2R, 50, chr1, cigar0),
				NewRecord("oH:::1:10:3000:5", chr1, 151, r2R, 50, chr1, cigarSoft1),
			},
			2500,
			&MetricsCollection{
				LibraryMetrics: map[string]*Metrics{
					"Unknown Library": &Metrics{
						ReadPairsExamined:   6,
						ReadPairDups:        4,
						ReadPairOpticalDups: 2,
					},
				},
			},
			[]string{"", "LB", "SQ", "SQ", "", "LB"},
		},
		{
			// Decoys drive up the file index of G1 within shard0.  MarkDuplicates must sort the duplicate set using global fileIdxs
			// rather than per-shard fileIdxs.
			[]*sam.Record{
				NewRecord("od1:::1:10:1:1", chr1, 20, r1F, 30, chr1, cigar0),
				NewRecord("od1:::1:10:1:1", chr1, 30, r2F, 20, chr1, cigar0),
				NewRecord("oJ:::1:10:1:1", chr1, 84, r1F, 150, chr1, cigar0),
				NewRecord("oK:::1:10:5:5", chr1, 86, r1F, 150, chr1, cigarSoft2),
				NewRecord("oJ:::1:10:1:1", chr1, 150, r2R, 84, chr1, cigar0),
				NewRecord("oK:::1:10:5:5", chr1, 150, r2R, 86, chr1, cigar0),
			},
			2500,
			&MetricsCollection{
				LibraryMetrics: map[string]*Metrics{
					"Unknown Library": &Metrics{
						ReadPairsExamined:   6,
						ReadPairDups:        2,
						ReadPairOpticalDups: 2,
					},
				},
			},
			[]string{"", "", "", "SQ", "", "SQ"},
		},
		{
			// Test disabling opticals.  Output should have no DT tags, even though the inputs have LB and SQ dups.
			// G is best, only I should be optical duplicate, because it shows up first in the file.  MarkDuplicates must sort the
			// duplicate set before determining optical duplicates,
			[]*sam.Record{
				NewRecord("oG:::1:10:1:1", chr1, 50, r1F, 150, chr1, cigar0),
				NewRecord("oH:::1:10:3000:5", chr1, 50, r1F, 151, chr1, cigar0),
				NewRecord("oI:::1:10:4000:5", chr1, 51, r1F, 151, chr1, cigarSoft1),
				NewRecord("oI:::1:10:4000:5", chr1, 150, r2R, 51, chr1, cigar0),
				NewRecord("oG:::1:10:1:1", chr1, 150, r2R, 50, chr1, cigar0),
				NewRecord("oH:::1:10:3000:5", chr1, 151, r2R, 50, chr1, cigarSoft1),
			},
			-1,
			&MetricsCollection{
				LibraryMetrics: map[string]*Metrics{
					"Unknown Library": &Metrics{
						ReadPairsExamined: 6,
						ReadPairDups:      4,
					},
				},
			},
			[]string{"", "", "", "", "", ""},
		},
		{ // Try with M before N.
			// L is primary.  Both M and N should be marked as optical dups.
			[]*sam.Record{
				NewRecord("oL:::1:10:1000:1", chr1, 50, r1F, 150, chr1, cigar0),
				NewRecord("oM:::1:10:3000:5", chr1, 50, r1F, 150, chr1, cigar0),
				NewRecord("oN:::1:10:4000:5", chr1, 50, r1F, 150, chr1, cigar0),
				NewRecord("oL:::1:10:1000:1", chr1, 150, r2R, 50, chr1, cigar0),
				NewRecord("oM:::1:10:3000:5", chr1, 150, r2R, 50, chr1, cigar0),
				NewRecord("oN:::1:10:4000:5", chr1, 150, r2R, 50, chr1, cigar0),
			},
			2500,
			&MetricsCollection{
				LibraryMetrics: map[string]*Metrics{
					"Unknown Library": &Metrics{
						ReadPairsExamined:   6,
						ReadPairDups:        4,
						ReadPairOpticalDups: 4,
					},
				},
			},
			[]string{"", "SQ", "SQ", "", "SQ", "SQ"},
		},
		{ // Try with N before M.
			// L is primary.  Both M and N should be marked as optical dups.
			[]*sam.Record{
				NewRecord("oL:::1:10:1000:1", chr1, 50, r1F, 150, chr1, cigar0),
				NewRecord("oN:::1:10:4000:5", chr1, 50, r1F, 150, chr1, cigar0),
				NewRecord("oM:::1:10:3000:5", chr1, 50, r1F, 150, chr1, cigar0),
				NewRecord("oL:::1:10:1000:1", chr1, 150, r2R, 50, chr1, cigar0),
				NewRecord("oN:::1:10:4000:5", chr1, 150, r2R, 50, chr1, cigar0),
				NewRecord("oM:::1:10:3000:5", chr1, 150, r2R, 50, chr1, cigar0),
			},
			2500,
			&MetricsCollection{
				LibraryMetrics: map[string]*Metrics{
					"Unknown Library": &Metrics{
						ReadPairsExamined:   6,
						ReadPairDups:        4,
						ReadPairOpticalDups: 4,
					},
				},
			},
			[]string{"", "SQ", "SQ", "", "SQ", "SQ"},
		},
		{
			// O is not an optical duplicate of A because their tiles don't match.
			[]*sam.Record{
				NewRecord("oA:::1:10:1:1", chr1, 0, r1F, 100, chr1, cigar0),
				NewRecord("oO:::1:11:5:5", chr1, 0, r1F, 100, chr1, cigar0),
				NewRecord("oA:::1:10:1:1", chr1, 100, r2R, 0, chr1, cigar0),
				NewRecord("oO:::1:11:5:5", chr1, 100, r2R, 0, chr1, cigar0),
			},
			2500,
			&MetricsCollection{
				LibraryMetrics: map[string]*Metrics{
					"Unknown Library": &Metrics{
						ReadPairsExamined:   4,
						ReadPairDups:        2,
						ReadPairOpticalDups: 0,
					},
				},
			},
			[]string{"", "LB", "", "LB"},
		},
		{
			// P is not an optical duplicate of A because their lanes don't match.
			[]*sam.Record{
				NewRecord("oA:::1:10:1:1", chr1, 0, r1F, 100, chr1, cigar0),
				NewRecord("oP:::2:10:5:5", chr1, 0, r1F, 100, chr1, cigar0),
				NewRecord("oA:::1:10:1:1", chr1, 100, r2R, 0, chr1, cigar0),
				NewRecord("oP:::2:10:5:5", chr1, 100, r2R, 0, chr1, cigar0),
			},
			2500,
			&MetricsCollection{
				LibraryMetrics: map[string]*Metrics{
					"Unknown Library": &Metrics{
						ReadPairsExamined:   4,
						ReadPairDups:        2,
						ReadPairOpticalDups: 0,
					},
				},
			},
			[]string{"", "LB", "", "LB"},
		},
	}

	tempDir, cleanup := testutil.TempDir(t, "", "")
	defer cleanup()
	for testIdx, test := range tests {
		for _, format := range []string{"bam", "pam"} {
			t.Logf("---- starting tests[%d] ----", testIdx)
			provider := bamprovider.NewFakeProvider(header, test.records)
			outputPath := NewTestOutput(tempDir, testIdx, format)
			opts := Opts{
				ShardSize:            100,
				Padding:              10,
				Parallelism:          1,
				QueueLength:          10,
				ClearExisting:        false,
				RemoveDups:           false,
				TagDups:              true,
				EmitUnmodifiedFields: true,
				OutputPath:           outputPath,
				Format:               format,
			}
			if test.opticalDistance >= 0 {
				opts.OpticalDetector = &TileOpticalDetector{
					OpticalDistance: test.opticalDistance,
				}
			}

			markDuplicates := &MarkDuplicates{
				Provider: provider,
				Opts:     &opts,
			}
			actualMetrics, err := markDuplicates.Mark(nil)
			assert.NoError(t, err)

			assert.Equal(t, len(test.metrics.LibraryMetrics), len(actualMetrics.LibraryMetrics))
			for k, m := range actualMetrics.LibraryMetrics {
				assert.Equal(t, *(test.metrics.LibraryMetrics[k]), *m)
			}
			actualRecords := ReadRecords(t, outputPath)
			for i, r := range actualRecords {
				t.Logf("output[%v]: %v", i, r)
				aux := r.AuxFields.Get(sam.Tag{'D', 'T'})
				if aux == nil {
					assert.Equal(t, test.dupType[i], "")
				} else {
					assert.Equal(t, test.dupType[i], aux.Value().(string))
				}
			}
		}
	}
}

// Test the Metrics that markDuplicates() returns.
func TestMetrics(t *testing.T) {
	// Notes that ReadPairsExamined, ReadPairDups, and
	// ReadPairOpticalDups below are doubled because they are halved
	// when written to the metrics file.
	tests := []struct {
		records []*sam.Record
		metrics *MetricsCollection
	}{
		{
			[]*sam.Record{
				NewRecord("A:::1:10:1:1", chr1, 0, r1F, 10, chr1, cigar0),
				NewRecord("A:::1:10:1:1", chr1, 10, r2R, 0, chr1, cigar0),
			},
			&MetricsCollection{
				LibraryMetrics: map[string]*Metrics{
					"Unknown Library": &Metrics{
						UnpairedReads:          0,
						ReadPairsExamined:      2,
						SecondarySupplementary: 0,
						UnmappedReads:          0,
						UnpairedDups:           0,
						ReadPairDups:           0,
						ReadPairOpticalDups:    0,
					},
				},
			},
		},
		{
			[]*sam.Record{
				NewRecord("A:::1:10:1:1", chr1, 0, r1F, 10, chr1, cigar0),
				NewRecord("B:::2:10:1:1", chr1, 0, r1F, 10, chr1, cigar0),
				NewRecord("A:::1:10:1:1", chr1, 10, r2R, 0, chr1, cigar0),
				NewRecord("B:::2:10:1:1", chr1, 10, r2R, 0, chr1, cigar0),
			},
			&MetricsCollection{
				LibraryMetrics: map[string]*Metrics{
					"Unknown Library": &Metrics{
						UnpairedReads:          0,
						ReadPairsExamined:      4,
						SecondarySupplementary: 0,
						UnmappedReads:          0,
						UnpairedDups:           0,
						ReadPairDups:           2,
						ReadPairOpticalDups:    0,
					},
				},
			},
		},
		{
			[]*sam.Record{
				NewRecord("A:::1:10:1:1", chr1, 0, r1F, 10, chr1, cigar0),
				NewRecord("B:::2:11:1:1", chr1, 0, r1F, 10, chr1, cigar0),
				NewRecord("C:::2:11:1:1", chr1, 0, r1F, 10, chr1, cigar0),
				NewRecord("A:::1:10:1:1", chr1, 10, r2R, 0, chr1, cigar0),
				NewRecord("B:::2:11:1:1", chr1, 10, r2R, 0, chr1, cigar0),
				NewRecord("C:::2:11:1:1", chr1, 10, r2R, 0, chr1, cigar0),
			},
			&MetricsCollection{
				LibraryMetrics: map[string]*Metrics{
					"Unknown Library": &Metrics{
						UnpairedReads:          0,
						ReadPairsExamined:      6,
						SecondarySupplementary: 0,
						UnmappedReads:          0,
						UnpairedDups:           0,
						ReadPairDups:           4,
						ReadPairOpticalDups:    2,
					},
				},
			},
		},
		{
			[]*sam.Record{
				NewRecord("A:::1:10:1:1", chr1, 0, s1F, 0, nil, cigar0),
				NewRecord("A:::1:10:1:1", chr1, 0, u2, 0, nil, cigar0),
				NewRecord("B:::2:10:1:1", chr1, 0, s1F, 0, nil, cigar0),
				NewRecord("B:::2:10:1:1", chr1, 0, u2, 0, nil, cigar0),
			},
			&MetricsCollection{
				LibraryMetrics: map[string]*Metrics{
					"Unknown Library": &Metrics{
						UnpairedReads:          2,
						ReadPairsExamined:      0,
						SecondarySupplementary: 0,
						UnmappedReads:          2,
						UnpairedDups:           1,
						ReadPairDups:           0,
						ReadPairOpticalDups:    0,
					},
				},
			},
		},
		{
			// Cross-shard pairs
			[]*sam.Record{
				NewRecord("A:::1:10:1:1", chr1, 0, r1F, 105, chr1, cigar0),
				NewRecord("B:::1:10:1:1", chr1, 0, r1F, 105, chr1, cigar0),
				NewRecord("A:::1:10:1:1", chr1, 105, r2R, 0, chr1, cigar0),
				NewRecord("B:::1:10:1:1", chr1, 105, r2R, 0, chr1, cigar0),
			},
			&MetricsCollection{
				LibraryMetrics: map[string]*Metrics{
					"Unknown Library": &Metrics{
						UnpairedReads:          0,
						ReadPairsExamined:      4,
						SecondarySupplementary: 0,
						UnmappedReads:          0,
						UnpairedDups:           0,
						ReadPairDups:           2,
						ReadPairOpticalDups:    2,
					},
				},
			},
		},
		{
			// Secondary alignments.
			[]*sam.Record{
				NewRecord("A:::1:10:1:1", chr1, 0, r1F, 105, chr1, cigar0),
				NewRecord("A:::1:10:1:1", chr1, 105, r2R, 0, chr1, cigar0),
				NewRecord("A:::1:10:1:1", chr2, 12, sec, 105, chr1, cigar0),
			},
			&MetricsCollection{
				LibraryMetrics: map[string]*Metrics{
					"Unknown Library": &Metrics{
						UnpairedReads:          0,
						ReadPairsExamined:      2,
						SecondarySupplementary: 1,
						UnmappedReads:          0,
						UnpairedDups:           0,
						ReadPairDups:           0,
						ReadPairOpticalDups:    0,
					},
				},
			},
		},
		{
			// Verify that we're counting dual-unmapped pairs
			[]*sam.Record{
				NewRecord("A:::1:10:1:1", chr1, 0, r1F, 105, chr1, cigar0),
				NewRecord("A:::1:10:1:1", chr1, 105, r2R, 0, chr1, cigar0),
				NewRecord("U:::2:11:1:1", nil, -1, up1, -1, nil, cigar0),
				NewRecord("U:::2:11:1:1", nil, -1, up2, -1, nil, cigar0),
			},
			&MetricsCollection{
				LibraryMetrics: map[string]*Metrics{
					"Unknown Library": &Metrics{
						UnpairedReads:          0,
						ReadPairsExamined:      2,
						SecondarySupplementary: 0,
						UnmappedReads:          2,
						UnpairedDups:           0,
						ReadPairDups:           0,
						ReadPairOpticalDups:    0,
					},
				},
			},
		},
	}
	tempDir, cleanup := testutil.TempDir(t, "", "")
	defer cleanup()
	for testIdx, test := range tests {
		for _, tagDups := range []bool{true, false} {
			for _, format := range []string{"bam", "pam"} {
				t.Logf("---- starting tests[%d, %v] ----", testIdx, tagDups)
				provider := bamprovider.NewFakeProvider(header, test.records)
				outputPath := NewTestOutput(tempDir, testIdx, format)
				opts := Opts{
					ShardSize:            100,
					Padding:              10,
					Parallelism:          1,
					QueueLength:          10,
					ClearExisting:        false,
					RemoveDups:           false,
					TagDups:              tagDups,
					EmitUnmodifiedFields: true,
					OutputPath:           outputPath,
					Format:               format,
					OpticalDetector: &TileOpticalDetector{
						OpticalDistance: 2500,
					},
				}
				markDuplicates := &MarkDuplicates{
					Provider: provider,
					Opts:     &opts,
				}
				actualMetrics, err := markDuplicates.Mark(nil)
				assert.NoError(t, err)

				for i, r := range ReadRecords(t, outputPath) {
					t.Logf("output[%v]: %v", i, r)
				}
				assert.Equal(t, len(test.metrics.LibraryMetrics), len(actualMetrics.LibraryMetrics))
				for k, m := range actualMetrics.LibraryMetrics {
					assert.Equal(t, *(test.metrics.LibraryMetrics[k]), *m)
				}
			}
		}
	}
}

func TestMetricsString(t *testing.T) {
	m := Metrics{
		UnpairedReads:          2,
		ReadPairsExamined:      8,
		SecondarySupplementary: 2,
		UnmappedReads:          1,
		UnpairedDups:           2,
		ReadPairDups:           4,
		ReadPairOpticalDups:    2,
	}

	assert.Equal(t, "2\t4\t2\t1\t2\t2\t1\t60.000000\t3", m.String())
}

func TestAlignDistCheck(t *testing.T) {
	var (
		max int
		m   sync.Mutex
	)
	c := maxAlignDistCheck{
		padding:            10,
		globalMaxAlignDist: &max,
		mutex:              &m,
	}
	assert.NoError(t, c.Process(NewRecord("A", chr1, 0, r1F, 100, chr1,
		[]sam.CigarOp{
			sam.NewCigarOp(sam.CigarMatch, 10),
		})))
	assert.NoError(t, c.Process(NewRecord("A", chr1, 0, r1F, 100, chr1,
		[]sam.CigarOp{
			sam.NewCigarOp(sam.CigarSoftClipped, 10),
			sam.NewCigarOp(sam.CigarMatch, 10),
		})))
	assert.NoError(t, c.Process(NewRecord("A", chr1, 0, r1F, 100, chr1,
		[]sam.CigarOp{
			sam.NewCigarOp(sam.CigarHardClipped, 10),
			sam.NewCigarOp(sam.CigarMatch, 10),
		})))
	assert.Error(t, c.Process(NewRecord("A", chr1, 0, r1F, 100, chr1,
		[]sam.CigarOp{
			sam.NewCigarOp(sam.CigarSoftClipped, 11),
			sam.NewCigarOp(sam.CigarMatch, 10),
		})),
		"alignment distance(%d) exceeds padding(%d) on read: %v", 11, 10, "A")
	assert.Error(t, c.Process(NewRecord("A", chr1, 0, r1F, 100, chr1,
		[]sam.CigarOp{
			sam.NewCigarOp(sam.CigarHardClipped, 11),
			sam.NewCigarOp(sam.CigarMatch, 10),
		})),
		"alignment distance(%d) exceeds padding(%d) on read: %v", 11, 10, "A")
	assert.Error(t, c.Process(NewRecord("A", chr1, 10, r1R, 100, chr1,
		[]sam.CigarOp{
			sam.NewCigarOp(sam.CigarMatch, 5),
			sam.NewCigarOp(sam.CigarDeletion, 2),
			sam.NewCigarOp(sam.CigarMatch, 5),
		})),
		"alignment distance(%d) exceeds padding(%d) on read: %v", 11, 10, "A")
	assert.Error(t, c.Process(NewRecord("A", chr1, 10, r1R, 100, chr1,
		[]sam.CigarOp{
			sam.NewCigarOp(sam.CigarMatch, 10),
			sam.NewCigarOp(sam.CigarSoftClipped, 2),
		})),
		"alignment distance(%d) exceeds padding(%d) on read: %v", 12, 10, "A")
	assert.Error(t, c.Process(NewRecord("A", chr1, 10, r1R, 100, chr1,
		[]sam.CigarOp{
			sam.NewCigarOp(sam.CigarMatch, 10),
			sam.NewCigarOp(sam.CigarHardClipped, 3),
		})),
		"alignment distance(%d) exceeds padding(%d) on read: %v", 13, 10, "A")
	assert.Error(t, c.Process(NewRecord("A", chr1, 0, r1R, 100, chr1,
		[]sam.CigarOp{
			sam.NewCigarOp(sam.CigarMatch, 12),
		})),
		"alignment distance(%d) exceeds padding(%d) on read: %v", 11, 10, "A")
}

func TestAlignDistCheckIntegration(t *testing.T) {
	tempDir, cleanup := testutil.TempDir(t, "", "")
	defer cleanup()

	testrecords := []*sam.Record{
		NewRecord("A", chr1, 10, r1R, 100, chr1,
			[]sam.CigarOp{
				sam.NewCigarOp(sam.CigarMatch, 5),
				sam.NewCigarOp(sam.CigarDeletion, 2),
				sam.NewCigarOp(sam.CigarMatch, 5),
			}),
	}
	provider := bamprovider.NewFakeProvider(header, testrecords)
	outputPath := NewTestOutput(tempDir, 0, "bam")

	opts := defaultOpts
	opts.OutputPath = outputPath
	opts.Format = "bam"
	markDuplicates := &MarkDuplicates{
		Provider: provider,
		Opts:     &opts,
	}

	_, err := markDuplicates.Mark(nil)
	assert.Error(t, err, "alignment distance(%d) exceeds padding(%d) on read: %v", 13, 10, "A")
}

func TestMetricsCollection(t *testing.T) {
	m := MetricsCollection{
		OpticalDistance: make([][]int64, 1),
	}
	m.OpticalDistance[0] = make([]int64, 10)
	m.AddDistance(2, 10)
}
