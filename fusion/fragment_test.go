package fusion

import (
	"testing"

	"github.com/grailbio/testutil/expect"
)

func TestFragment(t *testing.T) {
	frag := Fragment{
		Name:  "f0",
		R1Seq: "AACC",
		R2Seq: "GGTT",
		kmers: nil,
	}
	expect.EQ(t, frag.SubSeq(newCrossReadPosRange(0, 3)), "AAC")
	expect.EQ(t, frag.SubSeq(newCrossReadPosRange(0, newR2Pos(1))), "AACCG")
	expect.EQ(t, frag.SubSeq(newCrossReadPosRange(newR2Pos(0), newR2Pos(3))), "GGT")
}
