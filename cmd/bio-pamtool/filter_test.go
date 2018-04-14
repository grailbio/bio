package main

import (
	"testing"

	"github.com/biogo/hts/sam"
	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

func eval(t *testing.T, str string, r *sam.Record) bool {
	node, err := parseFilterExpr(str)
	require.NoError(t, err)
	return evaluateFilterExpr(node, r)
}

func TestParser(t *testing.T) {
	chr1, err := sam.NewReference("chr1", "", "", 10000, nil, nil)
	require.NoError(t, err)
	chr2, err := sam.NewReference("chr2", "", "", 10000, nil, nil)
	require.NoError(t, err)
	_, err = sam.NewHeader(nil, []*sam.Reference{chr1, chr2})
	require.NoError(t, err)

	rec, err := sam.NewRecord("read1",
		chr1,                             /*ref*/
		chr2,                             /*materef*/
		122 /*pos*/, 455 /*matepos*/, 20, /*templen*/
		60, /*mapq*/
		[]sam.CigarOp{sam.NewCigarOp(sam.CigarMatch, 4)}, /*cigar*/
		[]byte("ACGT"),
		[]byte{30, 31, 32, 33}, /*qual*/
		nil /*aux*/)
	require.NoError(t, err)
	assert.True(t, eval(t, "mapping_quality == 0x3c", rec))
	assert.True(t, eval(t, "(mapping_quality >= 60) && rec_name == \"read1\"", rec))
	assert.False(t, eval(t, "(mapping_quality >= 61) && rec_name == \"read1\"", rec))
	assert.True(t, eval(t, "(mapping_quality >= 61) || rec_name == \"read1\"", rec))
	assert.True(t, eval(t, "ref_name == \"chr1\"", rec))
	assert.True(t, eval(t, "re(ref_name, \"^ch.*1$\")", rec))
	assert.True(t, eval(t, "re(ref_name, \"c\")", rec))
	assert.False(t, eval(t, "re(ref_name, \"^chr$\")", rec))
	assert.True(t, eval(t, "!re(ref_name, \"^chr$\")", rec))
	assert.False(t, eval(t, "ref_name == \"chr2\"", rec))
	assert.True(t, eval(t, "position == 122", rec))
	assert.False(t, eval(t, "position != 122", rec))
	assert.True(t, eval(t, "sequence_length == 4", rec))
	assert.False(t, eval(t, "sequence_length < 4", rec))

	assert.True(t, eval(t, "mate_ref_name >= \"chr1\"", rec))
	assert.False(t, eval(t, "mate_ref_name < \"chr1\"", rec))
	assert.True(t, eval(t, "mate_position >= 455", rec))
	assert.False(t, eval(t, "mate_position < 455", rec))
	assert.True(t, eval(t, "template_length == 20", rec))
	assert.False(t, eval(t, "template_length != 20", rec))

	rec.Flags = sam.Paired
	assert.True(t, eval(t, "paired", rec))
	assert.False(t, eval(t, "!paired", rec))

	rec.Flags = sam.ProperPair | sam.Secondary | sam.Read1
	assert.True(t, eval(t, "proper_pair && secondary_alignment", rec))
	assert.True(t, eval(t, "proper_pair && first_of_pair", rec))
	assert.False(t, eval(t, "proper_pair && second_of_pair", rec))
}
