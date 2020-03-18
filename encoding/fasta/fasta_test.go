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
	fastaData = `>seq1
AcGTA
CGTAC
GT
>seq2 A viral sequence
ACGT
ACGT
`
	fastaIndex = `seq1	12	6	5	6
seq2	8	44	4	5
`
}

func TestOps(t *testing.T) {
	type impl struct {
		name string
		fa   func() fasta.Fasta
	}
	impls := []impl{
		{"unidx", func() fasta.Fasta {
			fa, err := fasta.New(strings.NewReader(fastaData), fasta.OptClean)
			assert.NoError(t, err)
			return fa
		}},
		{"idx", func() fasta.Fasta {
			fa, err := fasta.NewIndexed(strings.NewReader(fastaData), strings.NewReader(fastaIndex), fasta.OptClean)
			assert.NoError(t, err)
			return fa
		}},
	}
	for _, impl := range impls {
		fa := impl.fa()
		t.Run(impl.name, func(t *testing.T) {
			t.Run("Get", func(t *testing.T) {
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
				for _, tt := range tests {
					got, err := fa.Get(tt.seq, tt.start, tt.end)
					if (err == nil && tt.err != nil) || (err != nil && tt.err == nil) {
						t.Errorf("unexpected error: want %v, got %v", tt.err, err)
					}
					if got != tt.want {
						t.Errorf("unexpected sequence: want %s, got %s", tt.want, got)
					}
				}
			})
			t.Run("Length", func(t *testing.T) {
				tests := []struct {
					seq  string
					want uint64
					err  error
				}{
					{"seq1", 12, nil},
					{"seq2", 8, nil},
					{"seq0", 0, fmt.Errorf("sequence not found in index: seq0")},
				}
				for _, tt := range tests {
					got, err := fa.Len(tt.seq)
					if (err == nil && tt.err != nil) || (err != nil && tt.err == nil) {
						t.Errorf("unexpected error: want %v, got %v", tt.err, err)
					}
					if got != tt.want {
						t.Errorf("unexpected length: want %v, got %v", tt.want, got)
					}
				}
			})
			t.Run("SeqNames", func(t *testing.T) {
				want := sort.StringSlice([]string{"seq1", "seq2"})
				want.Sort()
				got := sort.StringSlice(fa.SeqNames())
				got.Sort()
				if !reflect.DeepEqual(got, want) {
					t.Errorf("got %v, want %v", got, want)
				}
			})
		})
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
)

// On an AWS EC2 m5d.4xlarge where /tmp resides on local NVME SSD:
//
//   $ go test github.com/grailbio/bio/encoding/fasta -bench BenchmarkRead -benchmem -path /tmp/hg19.fa -index-path /tmp/hg19.fa.fai
//   goos: linux
//   goarch: amd64
//   pkg: github.com/grailbio/bio/encoding/fasta
//   BenchmarkRead/unidx/init-16            1        9285597629 ns/op        21848127712 B/op        32723733 allocs/op
//   BenchmarkRead/unidx/read_all-16                      763           1560333 ns/op          723851 B/op      20064 allocs/op
//   BenchmarkRead/unidx/read_rand-16                 2925188               408 ns/op             128 B/op          4 allocs/op
//   BenchmarkRead/idx/init-16                            121           9656116 ns/op         2575314 B/op      15212 allocs/op
//   BenchmarkRead/idx/read_all-16                          1        8822692112 ns/op        3771464752 B/op    25703 allocs/op
//   BenchmarkRead/idx/read_rand-16                    228588              5114 ns/op             970 B/op          4 allocs/op
//   PASS
//   ok      github.com/grailbio/bio/encoding/fasta  29.222s
func BenchmarkRead(b *testing.B) {
	if *pathFlag == "" {
		b.Fatal("-path not set")
	}
	if *idxPathFlag == "" {
		b.Fatal("-index-path not set")
	}
	ctx := vcontext.Background()
	type impl struct {
		name string
		fa   func() fasta.Fasta
	}
	impls := []impl{
		impl{"unidx", func() fasta.Fasta {
			in, err := file.Open(ctx, *pathFlag)
			assert.NoError(b, err)
			fa, err := fasta.New(in.Reader(ctx), fasta.OptClean)
			assert.NoError(b, err)
			return fa
		}},
		impl{"idx", func() fasta.Fasta {
			in, err := file.Open(ctx, *pathFlag)
			assert.NoError(b, err)
			idxIn, err := file.Open(ctx, *idxPathFlag)
			assert.NoError(b, err)
			fa, err := fasta.NewIndexed(in.Reader(ctx), idxIn.Reader(ctx), fasta.OptClean)
			assert.NoError(b, err)
			return fa
		}},
	}
	for _, impl := range impls {
		impl := impl
		var fa fasta.Fasta
		b.Run(impl.name, func(b *testing.B) {
			b.Run("init", func(b *testing.B) {
				for i := 0; i < b.N; i++ {
					fa = impl.fa()
				}
			})
			b.Run("read_all", func(b *testing.B) {
				for i := 0; i < b.N; i++ {
					var totalSeqLength int64
					seqNames := append([]string{}, fa.SeqNames()...)
					for _, seqName := range seqNames {
						n, err := fa.Len(seqName)
						assert.NoError(b, err)
						seq, err := fa.Get(seqName, 0, n)
						assert.NoError(b, err)
						totalSeqLength += int64(len(seq))
					}
					assert.GT(b, totalSeqLength, int64(0))
				}
			})
			b.Run("read_rand", func(b *testing.B) {
				seqNames := append([]string{}, fa.SeqNames()...)
				var totalSeqLength int64
				for i := 0; i < b.N; i++ {
					seqName := seqNames[rand.Intn(len(seqNames))]
					seqLen, err := fa.Len(seqName)
					assert.NoError(b, err)
					readLen := rand.Intn(int(seqLen)-1) + 1
					if readLen > 1000 {
						readLen = 1000
					}
					readStart := rand.Intn(int(seqLen) - readLen)
					seq, err := fa.Get(seqName, uint64(readStart), uint64(readStart+readLen))
					assert.NoError(b, err)
					totalSeqLength += int64(len(seq))
				}
				assert.GT(b, totalSeqLength, int64(0))
			})
		})
	}
}
