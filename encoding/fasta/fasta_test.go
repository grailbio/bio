package fasta_test

import (
	"bytes"
	"flag"
	"fmt"
	"io"
	"math/rand"
	"reflect"
	"sort"
	"strings"
	"testing"

	"github.com/grailbio/base/file"
	"github.com/grailbio/base/vcontext"
	"github.com/grailbio/bio/encoding/fasta"
	"github.com/grailbio/testutil/assert"
)

var fastaData string
var fastaIndex string

func init() {
	fastaData = ">seq1\n" + "ACGTA\nCGTAC\nGT\n" + ">seq2 A viral sequence\n" + "ACGT\n" + "ACGT\n"
	fastaIndex = "seq1\t12\t6\t5\t6\n" + "seq2\t8\t44\t4\t5\n"
}

func TestGet(t *testing.T) {
	tests := []struct {
		seq   string
		start uint64
		end   uint64
		want  string
		err   error
	}{
		{"seq1", 1, 2, "C", nil},
		{"seq1", 1, 6, "CGTAC", nil},
		{"seq1", 0, 12, "ACGTACGTACGT", nil},
		{"seq1", 10, 12, "GT", nil},
		{"seq2", 0, 8, "ACGTACGT", nil},
		{"seq2", 2, 5, "GTA", nil},
		{"seq0", 0, 1, "", fmt.Errorf("sequence not found in index: seq0")},
		{"seq1", 10, 13, "", fmt.Errorf("end is past end of sequence seq1: 12")},
		{"seq1", 4, 3, "", fmt.Errorf("start must be less than end")},
	}
	unindexed, err := fasta.New(strings.NewReader(fastaData))
	if err != nil {
		t.Errorf("couldn't create Fasta: %v", err)
	}
	indexed, err := fasta.NewIndexed(strings.NewReader(fastaData), strings.NewReader(fastaIndex))
	if err != nil {
		t.Errorf("couldn't read index: %v", err)
	}
	for _, tt := range tests {
		got, err := unindexed.Get(tt.seq, tt.start, tt.end)
		if (err == nil && tt.err != nil) || (err != nil && tt.err == nil) {
			t.Errorf("unexpected error: want %v, got %v", tt.err, err)
		}
		if got != tt.want {
			t.Errorf("unexpected sequence: want %s, got %s", tt.want, got)
		}

		got, err = indexed.Get(tt.seq, tt.start, tt.end)
		if (err == nil && tt.err != nil) || (err != nil && tt.err == nil) {
			t.Errorf("unexpected error: want %v, got %v", tt.err, err)
		}
		if got != tt.want {
			t.Errorf("unexpected sequence: want %s, got %s", tt.want, got)
		}
	}
}

func TestLength(t *testing.T) {
	tests := []struct {
		seq  string
		want uint64
		err  error
	}{
		{"seq1", 12, nil},
		{"seq2", 8, nil},
		{"seq0", 0, fmt.Errorf("sequence not found in index: seq0")},
	}
	unindexed, err := fasta.New(strings.NewReader(fastaData))
	if err != nil {
		t.Errorf("couldn't create Fasta: %v", err)
	}
	indexed, err := fasta.NewIndexed(strings.NewReader(fastaData), strings.NewReader(fastaIndex))
	if err != nil {
		t.Errorf("couldn't read index: %v", err)
	}
	for _, tt := range tests {
		got, err := unindexed.Len(tt.seq)
		if (err == nil && tt.err != nil) || (err != nil && tt.err == nil) {
			t.Errorf("unexpected error: want %v, got %v", tt.err, err)
		}
		if got != tt.want {
			t.Errorf("unexpected length: want %v, got %v", tt.want, got)
		}

		got, err = indexed.Len(tt.seq)
		if (err == nil && tt.err != nil) || (err != nil && tt.err == nil) {
			t.Errorf("unexpected error: want %v, got %v", tt.err, err)
		}
		if got != tt.want {
			t.Errorf("unexpected length: want %v, got %v", tt.want, got)
		}
	}
}

func TestSeqNames(t *testing.T) {
	unindexed, err := fasta.New(strings.NewReader(fastaData))
	if err != nil {
		t.Errorf("couldn't create Fasta: %v", err)
	}
	indexed, err := fasta.NewIndexed(strings.NewReader(fastaData), strings.NewReader(fastaIndex))
	if err != nil {
		t.Errorf("couldn't read index: %v", err)
	}
	want := sort.StringSlice([]string{"seq1", "seq2"})
	want.Sort()
	got := sort.StringSlice(unindexed.SeqNames())
	got.Sort()
	if !reflect.DeepEqual(got, want) {
		t.Errorf("got %v, want %v", got, want)
	}

	got = sort.StringSlice(indexed.SeqNames())
	got.Sort()
	if !reflect.DeepEqual(got, want) {
		t.Errorf("got %v, want %v", got, want)
	}
}

func TestFastaFaiToReferenceLengths(t *testing.T) {
	type ref struct {
		chrom  string
		length uint64
	}

	var testFai bytes.Buffer
	testFai.Write([]byte("chr1\t250000000\t6\t60\t61\n"))
	testFai.Write([]byte("chr2\t199000000\t6\t60\t61\n"))
	testFaiReader := bytes.NewReader(testFai.Bytes())

	tests := []struct {
		reader     io.Reader
		references []ref
	}{
		{testFaiReader,
			[]ref{ref{chrom: "chr1", length: uint64(250000000)},
				ref{chrom: "chr2", length: uint64(199000000)}},
		},
	}

	for _, test := range tests {
		faiReader := test.reader
		var result map[string]uint64
		result, err := fasta.FaiToReferenceLengths(faiReader)
		if err != nil {
			t.Errorf("error generating reference lengths: %v", err)
		}

		for _, testData := range test.references {
			reference := testData.chrom
			length := testData.length
			if val, ok := result[reference]; ok {
				if val != length {
					t.Errorf("error reading fasta index: got %d, want %d", val, length)
				} else {
					fmt.Printf("read fasta index: got %d, want %d\n", val, length)
				}
			}
		}
	}
}

func TestGenerateIndex(t *testing.T) {
	generateIndex := func(fa string) (faidx string) {
		idx := bytes.Buffer{}
		assert.NoError(t, fasta.GenerateIndex(&idx, strings.NewReader(fa)))
		return idx.String()
	}

	fa := `>E0
GGTGAAATC
CCTGAAATC
AAAATTGCT
>E1
GTCCCTCCCCAGACATGGCCCTGGGAGGC
>E2
CCGCGCCCGCGCCCCCGCCGCC
>E3
GTCAAGGTTGCACAG
>E4
ATGAATCATGTGGTAAAA
`
	fai := generateIndex(fa)
	assert.EQ(t, fai, `E0	27	4	9	10
E1	29	38	29	30
E2	22	72	22	23
E3	15	99	15	16
E4	18	119	18	19
`)
	// Read using the generated index
	indexed, err := fasta.NewIndexed(strings.NewReader(fa), strings.NewReader(fai))
	assert.NoError(t, err)
	l, err := indexed.Len("E3")
	assert.NoError(t, err)
	assert.EQ(t, l, uint64(15))
	seq, err := indexed.Get("E3", 0, l)
	assert.NoError(t, err)
	assert.EQ(t, seq, "GTCAAGGTTGCACAG")

	// MO-DOS newline encodinng.
	assert.EQ(t, generateIndex(">E0\r\nGGGG\r\n>E1\r\nAAAAA\r\n"),
		`E0	4	5	4	6
E1	5	16	5	7
`)

	// No newline at the end.
	assert.EQ(t, generateIndex(">E0\nGGGG\n>E1\nCCCCC\nAAAAA"),
		`E0	4	4	4	5
E1	10	13	5	6
`)
	// Note: samtool faidx emits "5 13 5 6" for E1, but "5 13 5 5" is correct
	// according to the spec.
	assert.EQ(t, generateIndex(">E0\nGGGG\n>E1\nAAAAA"),
		`E0	4	4	4	5
E1	5	13	5	5
`)

	idx := bytes.Buffer{}
	assert.Regexp(t, fasta.GenerateIndex(&idx, strings.NewReader("")), "empty FASTA")
}

var (
	pathFlag    = flag.String("path", "", "FASTA file used by benchmarks")
	idxPathFlag = flag.String("index-path", "", "FASTA index file used by benchmarks")
	shuffleFlag = flag.Bool("shuffle", false, "Read sequences in random order")
)

func BenchmarkRead(b *testing.B) {
	if *pathFlag == "" {
		b.Skip("--path not set")
	}
	for i := 0; i < b.N; i++ {
		ctx := vcontext.Background()
		in, err := file.Open(ctx, *pathFlag)
		assert.NoError(b, err)

		var (
			fin   fasta.Fasta
			idxIn file.File
		)
		if *idxPathFlag != "" {
			idxIn, err = file.Open(ctx, *idxPathFlag)
			assert.NoError(b, err)
			fin, err = fasta.NewIndexed(in.Reader(ctx), idxIn.Reader(ctx))
			assert.NoError(b, err)
		} else {
			fin, err = fasta.New(in.Reader(ctx))
			assert.NoError(b, err)
		}
		seqNames := append([]string{}, fin.SeqNames()...)
		if *shuffleFlag {
			rand.Shuffle(len(seqNames), func(i, j int) {
				seqNames[i], seqNames[j] = seqNames[j], seqNames[i]
			})
		}
		for _, seq := range seqNames {
			n, err := fin.Len(seq)
			assert.NoError(b, err)
			_, err = fin.Get(seq, 0, n)
			assert.NoError(b, err)
		}
		if idxIn != nil {
			assert.NoError(b, idxIn.Close(ctx))
		}
		assert.NoError(b, in.Close(ctx))
	}
}
