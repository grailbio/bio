package main

import (
	"context"
	"io/ioutil"
	"os"
	"sort"
	"testing"

	"github.com/grailbio/bio/encoding/fasta"
	"github.com/grailbio/testutil"
	"github.com/grailbio/testutil/assert"
	"github.com/grailbio/testutil/expect"
)

func TestParseGencodeWithoutPadding(t *testing.T) {
	testGenerateTranscriptome(t,
		gencodeFlags{
			exonPadding:            0,
			codingOnly:             false,
			separateJns:            false,
			retainedExonBases:      0,
			keepMitochondrialGenes: true,
			wholeGenes:             false,
			collapseTranscripts:    false},
		"file_1.fa")
}

func TestParseGencodeWithPadding(t *testing.T) {
	testGenerateTranscriptome(t,
		gencodeFlags{
			exonPadding:            20,
			codingOnly:             false,
			separateJns:            false,
			retainedExonBases:      0,
			keepMitochondrialGenes: true,
			wholeGenes:             false,
			collapseTranscripts:    false},
		"file_2.fa")
}

func TestParseGencodeOnlyCoding(t *testing.T) {
	testGenerateTranscriptome(t,
		gencodeFlags{
			exonPadding:            0,
			codingOnly:             true,
			separateJns:            false,
			retainedExonBases:      0,
			keepMitochondrialGenes: true,
			wholeGenes:             false,
			collapseTranscripts:    false},
		"file_3.fa")
}

func TestParseGencodeSeparateJunctions(t *testing.T) {
	testGenerateTranscriptome(t,
		gencodeFlags{
			exonPadding:            20,
			codingOnly:             false,
			separateJns:            true,
			retainedExonBases:      5,
			keepMitochondrialGenes: true,
			wholeGenes:             false,
			collapseTranscripts:    false},
		"file_4.fa")
}

func TestParseGencodeWholeGenesWithoutPadding(t *testing.T) {
	testGenerateTranscriptome(t,
		gencodeFlags{
			exonPadding:            0,
			codingOnly:             false,
			separateJns:            false,
			retainedExonBases:      0,
			keepMitochondrialGenes: true,
			wholeGenes:             true,
			collapseTranscripts:    false},
		"file_5.fa") // expectedFileName
}

func TestParseGencodeWholeGenesWithPadding(t *testing.T) {
	testGenerateTranscriptome(t,
		gencodeFlags{
			exonPadding:            20,
			codingOnly:             false,
			separateJns:            false,
			retainedExonBases:      5,
			keepMitochondrialGenes: true,
			wholeGenes:             true,
			collapseTranscripts:    false},
		"file_6.fa") // expectedFileName
}

func TestParseGencodeWholeGenesWithPaddingCodingOnly(t *testing.T) {
	testGenerateTranscriptome(t,
		gencodeFlags{
			exonPadding:         20,
			codingOnly:          true,
			separateJns:         false,
			retainedExonBases:   0,
			wholeGenes:          true,
			collapseTranscripts: false},
		"file_7.fa") // expectedFileName
}

func TestParseGencodeCollapseTranscriptsWithoutPadding(t *testing.T) {
	testGenerateTranscriptome(t,
		gencodeFlags{
			exonPadding:            0,
			codingOnly:             false,
			separateJns:            false,
			retainedExonBases:      0,
			keepMitochondrialGenes: true,
			wholeGenes:             false,
			collapseTranscripts:    true},
		"file_8.fa") // expectedFileName
}

func TestParseGencodeCollapseTranscriptsWithPadding(t *testing.T) {
	testGenerateTranscriptome(t,
		gencodeFlags{
			exonPadding:            20,
			codingOnly:             false,
			separateJns:            false,
			retainedExonBases:      0,
			keepMitochondrialGenes: true,
			wholeGenes:             false,
			collapseTranscripts:    true},
		"file_9.fa")
}

func TestParseGencodeCollapseTranscriptsWithPaddingCodingOnly(t *testing.T) {
	testGenerateTranscriptome(t,
		gencodeFlags{
			exonPadding:            20,
			codingOnly:             true,
			separateJns:            false,
			retainedExonBases:      0,
			keepMitochondrialGenes: true,
			wholeGenes:             false,
			collapseTranscripts:    true},
		"file_10.fa")
}

type fastaLine struct{ key, seq string }

func readFASTA(t *testing.T, path string) []fastaLine {
	in, err := os.Open(path)
	assert.NoError(t, err)
	defer in.Close()
	f, err := fasta.New(in)
	assert.NoError(t, err)
	var r []fastaLine
	for _, key := range f.SeqNames() {
		len, err := f.Len(key)
		assert.NoError(t, err)
		seq, err := f.Get(key, 0, len)
		assert.NoError(t, err)
		r = append(r, fastaLine{key, seq})
	}
	sort.SliceStable(r, func(i, j int) bool {
		return r[i].key < r[j].key
	})
	return r
}

// Pluggable test for all cases
func testGenerateTranscriptome(t *testing.T, flags gencodeFlags, want string) {
	tempdir := testutil.GetTmpDir()
	flags.output = tempdir + "/" + want
	GenerateTranscriptome(context.Background(),
		testutil.GetFilePath("//go/src/github.com/grailbio/bio/fusion/parsegencode/testdata/annotation.gtf"),
		testutil.GetFilePath("//go/src/github.com/grailbio/bio/fusion/parsegencode/testdata/genome.fa"),
		flags)
	if *updateGoldenFlag {
		data, err := ioutil.ReadFile(flags.output)
		assert.NoError(t, err)
		assert.NoError(t, ioutil.WriteFile(testutil.GetFilePath("//go/src/github.com/grailbio/bio/fusion/parsegencode/testdata/"+want), data, 0644))
	} else {
		expect.EQ(t,
			readFASTA(t, testutil.GetFilePath("//go/src/github.com/grailbio/bio/fusion/parsegencode/testdata/"+want)),
			readFASTA(t, flags.output))
	}
}
