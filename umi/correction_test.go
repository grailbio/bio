package umi

import (
	"os"
	"strings"
	"testing"

	"github.com/grailbio/base/grail"
	"github.com/stretchr/testify/assert"
)

func TestAllKmers(t *testing.T) {
	assertValidKmer := func(kmer string) {
		for _, c := range strings.ToUpper(kmer) {
			assert.True(t, c == 'A' || c == 'C' || c == 'G' || c == 'T' || c == 'N',
				"%s is not a valid kmer", kmer)
		}
	}

	kmers := allKmers(3, alphabetWithN)
	uniq := map[string]bool{}
	for _, kmer := range kmers {
		assertValidKmer(kmer)
		uniq[kmer] = true
	}
	assert.Equal(t, 125, len(uniq)) // 5^3 possible kmers including ACGTN.
}

func TestSnapCorrector(t *testing.T) {
	known3 := "AAA\nCCC\nGGG\nTTT"
	known4 := "AAAA\nCCCC\nGGGG\nTTTT"

	tests := []struct {
		knownUMIs   string
		umi         string
		expected    string
		edits       int
		correctable bool
	}{
		{known3, "AAA", "AAA", 0, false},
		{known3, "TAA", "AAA", 1, true},
		{known3, "ATA", "AAA", 1, true},
		{known3, "AAT", "AAA", 1, true},
		{known3, "NAA", "AAA", 1, true},

		{known4, "AACC", "AACC", -1, false}, // Could be AAAA or CCCC
		{known4, "AANN", "AAAA", 2, true},
		{known4, "ANNN", "AAAA", 3, true},
		{known4, "NNNN", "NNNN", -1, false},
	}

	for _, test := range tests {
		c := NewSnapCorrector([]byte(test.knownUMIs))
		correctedUMI, edits, corrected := c.CorrectUMI(test.umi)
		assert.Equal(t, test.expected, correctedUMI, "'%s' should have corrected to '%s'", test.umi, test.expected)
		assert.Equal(t, test.edits, edits, "'%s' should have corrected to '%s' with %d edits", test.umi, test.expected, test.edits)
		assert.Equal(t, test.correctable, corrected, "'%s' should have corrected %v", test.umi, test.correctable)
	}
}

func TestMain(m *testing.M) {
	shutdown := grail.Init()
	defer shutdown()
	os.Exit(m.Run())
}
