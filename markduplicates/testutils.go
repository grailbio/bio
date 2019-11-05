package markduplicates

import (
	"fmt"
	"io"
	"os"
	"path/filepath"
	"strings"
	"testing"

	gbam "github.com/grailbio/bio/encoding/bam"
	"github.com/grailbio/bio/encoding/bamprovider"
	"github.com/grailbio/hts/bam"
	"github.com/grailbio/hts/sam"
	"github.com/grailbio/testutil"
	"github.com/stretchr/testify/assert"
)

type TestRecord struct {
	R              *sam.Record
	DupFlag        bool
	ExpectedAuxs   []sam.Aux
	UnexpectedTags []sam.Tag
}

type TestCase struct {
	TRecords []TestRecord
	Opts     Opts
}

func NewRecord(name string, ref *sam.Reference, pos int, flags sam.Flags, matePos int, mateRef *sam.Reference, cigar sam.Cigar) *sam.Record {
	r := sam.GetFromFreePool()
	r.Name = name
	r.Ref = ref
	r.Pos = pos
	r.MatePos = matePos
	r.MateRef = mateRef
	r.Flags = flags
	r.Cigar = cigar
	return r
}

func NewRecordSeq(name string, ref *sam.Reference, pos int, flags sam.Flags, matePos int, mateRef *sam.Reference,
	cigar sam.Cigar, seq, qual string) *sam.Record {
	if len(seq) != len(qual) {
		panic("seq and qual must be equal length")
	}
	r := sam.GetFromFreePool()
	r.Name = name
	r.Ref = ref
	r.Pos = pos
	r.MatePos = matePos
	r.MateRef = mateRef
	r.Flags = flags
	r.Cigar = cigar
	r.Seq = sam.NewSeq([]byte(seq))
	r.Qual = []byte(qual)
	return r
}

func NewRecordAux(name string, ref *sam.Reference, pos int, flags sam.Flags, matePos int, mateRef *sam.Reference,
	cigar sam.Cigar, aux sam.Aux) *sam.Record {
	r := NewRecord(name, ref, pos, flags, matePos, mateRef, cigar)
	r.AuxFields = append(r.AuxFields, aux)
	return r
}

func NewAux(name string, val interface{}) sam.Aux {
	aux, err := sam.NewAux(sam.NewTag(name), val)
	if err != nil {
		panic(fmt.Sprintf("error creating %s %v tag: %v", name, val, err))
	}
	return aux
}

func RunTestCases(t *testing.T, header *sam.Header, cases []TestCase) {
	tempDir, cleanup := testutil.TempDir(t, "", "")
	defer cleanup()

	for testIdx, test := range cases {
		for _, format := range []string{"bam", "pam"} {
			t.Logf("---- starting TestCase[%d] ----", testIdx)
			testrecords := make([]*sam.Record, 0, len(test.TRecords))
			for _, tr := range test.TRecords {
				testrecords = append(testrecords, tr.R)
			}
			provider := bamprovider.NewFakeProvider(header, testrecords)

			outputPath := NewTestOutput(tempDir, testIdx, format)

			markDuplicates := &MarkDuplicates{
				Provider: provider,
				Opts:     &test.Opts,
			}
			markDuplicates.Opts.OutputPath = outputPath
			markDuplicates.Opts.Format = format

			_, err := markDuplicates.Mark(nil)
			assert.NoError(t, err)
			for i, r := range testrecords {
				t.Logf("input[%v]: %v begin %d end %d", i, r, r.Start(), r.End())
			}

			actualRecords := ReadRecords(t, outputPath)
			assert.Equal(t, len(test.TRecords), len(actualRecords))
			for i, r := range actualRecords {
				t.Logf("output[%v]: %v", i, r)

				assert.Equal(t, test.TRecords[i].DupFlag, r.Flags&sam.Duplicate != 0, "duplicate flag is wrong")

				// Verify that exactly one of each expected tag exists, and has the right value.
				for _, expectedAux := range test.TRecords[i].ExpectedAuxs {
					found := 0
					for _, aux := range r.AuxFields {
						if aux[0] == expectedAux.Tag()[0] && aux[1] == expectedAux.Tag()[1] {
							assert.Equal(t, expectedAux, aux)
							found++
						}
					}
					assert.Equal(t, 1, found, "Incorrect number of %s tags, expected 1, got %d",
						expectedAux.Tag(), found)
				}
				// Verify that these tags do not exist.
				for _, negTag := range test.TRecords[i].UnexpectedTags {
					actual, ok := r.Tag([]byte{negTag[0], negTag[1]})
					assert.Equal(t, false, ok, "Expected tag to be absent, but it exists: %v", actual)
				}
			}
		}
	}
}

// NewTestOutput returns different string filename for the different output formats.
func NewTestOutput(dir string, index int, format string) string {
	switch format {
	case "bam":
		return filepath.Join(dir, fmt.Sprintf("%d.bam", index))
	case "pam":
		return filepath.Join(dir, fmt.Sprintf("%d.pam", index))
	}
	panic(format)
}

// ReadRecords reads the records from path and returns them as a slice, in order.
func ReadRecords(t *testing.T, path string) []*sam.Record {
	records := make([]*sam.Record, 0)
	if strings.HasSuffix(path, ".bam") {
		// BAM files produced by this test don't have indexes, so read them using
		// the raw reader.
		in, err := os.Open(path)
		assert.NoError(t, err)
		defer func() {
			assert.NoError(t, in.Close())
		}()
		reader, err := bam.NewReader(in, 1)
		assert.NoError(t, err)
		for {
			r, err := reader.Read()
			if err == io.EOF {
				break
			}
			assert.NoError(t, err)
			records = append(records, r)
		}
	} else {
		p := bamprovider.NewProvider(path)
		header, err := p.GetHeader()
		assert.NoError(t, err)
		iter := p.NewIterator(gbam.UniversalShard(header))
		for iter.Scan() {
			records = append(records, iter.Record())
		}
		assert.NoError(t, iter.Close())
		assert.NoError(t, p.Close())
	}
	return records
}
