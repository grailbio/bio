package fusion

import (
	"testing"

	"github.com/grailbio/testutil/expect"
)

func TestParseTranscriptomeReferenceName(t *testing.T) {
	ensemblID, gene, chrom, start, end, index, err := ParseTranscriptomeKey("ENST00000279783.3|OR8K1|chr11:56346039-56346998:1051|960")
	expect.NoError(t, err)
	expect.EQ(t, ensemblID, "ENST00000279783.3")
	expect.EQ(t, gene, "OR8K1")
	expect.EQ(t, chrom, "chr11")
	expect.EQ(t, start, 56346039)
	expect.EQ(t, end, 56346998)
	expect.EQ(t, index, 1051)
}
