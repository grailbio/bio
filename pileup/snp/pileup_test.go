package snp_test

import (
	"path/filepath"
	"testing"

	"github.com/grailbio/base/file"
	"github.com/grailbio/base/recordio"
	"github.com/grailbio/base/vcontext"
	gbam "github.com/grailbio/bio/encoding/bam"
	"github.com/grailbio/bio/pileup/snp"
	"github.com/grailbio/hts/bam"
	"github.com/grailbio/hts/sam"
	"github.com/grailbio/testutil"
	"github.com/grailbio/testutil/assert"
)

func TestPileup(t *testing.T) {
	// Write a temporary BED file.
	tmpdir, cleanup := testutil.TempDir(t, "", "")
	defer testutil.NoCleanupOnError(t, cleanup, tmpdir)

	ctx := vcontext.Background()
	bedpath := filepath.Join(tmpdir, "tmp.bed")
	out, err := file.Create(ctx, bedpath)
	assert.NoError(t, err)
	_, err = out.Writer(ctx).Write([]byte("chr2_subset\t100000\t100005\nchr2_subset\t100006\t100007\n"))
	assert.NoError(t, err)
	err = out.Close(ctx)
	assert.NoError(t, err)

	// Must define these for bamWriter.Write() calls to work.
	ref, _ := sam.NewReference("chr2_subset", "", "", 124994, nil, nil)
	samHeader, _ := sam.NewHeader(nil, []*sam.Reference{ref})

	// Will make this vary by test in the future.
	gbaipath := filepath.Join(tmpdir, "tmp.bam.gbai")
	outPrefix := filepath.Join(tmpdir, "bio-pileup")
	opts := snp.Opts{
		BedPath:      bedpath,
		BamIndexPath: gbaipath,
		Clip:         0,
		FlagExclude:  0xf00,
		Mapq:         60,
		MaxReadLen:   500,
		MaxReadSpan:  511,
		MinBagDepth:  0,
		MinBaseQual:  43,
		Parallelism:  1,
		PerStrand:    false,
		RemoveSq:     false,
		Stitch:       true,
		TempDir:      "",
	}

	tests := []struct {
		name  string
		reads []sam.Record
		want  []snp.BaseStrandPile
	}{
		// Basic read-pair with no overlap, and high enough base-quals to pass the
		// filter.
		{
			name: "no_overlap",
			reads: []sam.Record{
				{
					Name:  "read1",
					Ref:   ref,
					Pos:   100000,
					MapQ:  60,
					Cigar: []sam.CigarOp{sam.NewCigarOp(sam.CigarMatch, 2)},
					// + strand
					Flags:   sam.Paired | sam.ProperPair | sam.MateReverse | sam.Read1,
					MateRef: ref,
					MatePos: 100003,
					Seq:     sam.NewSeq([]byte("AC")),
					Qual:    []byte{43, 43},
				},
				{
					Name:    "read1",
					Ref:     ref,
					Pos:     100003,
					MapQ:    60,
					Cigar:   []sam.CigarOp{sam.NewCigarOp(sam.CigarMatch, 2)},
					Flags:   sam.Paired | sam.ProperPair | sam.Reverse | sam.Read2,
					MateRef: ref,
					MatePos: 100000,
					Seq:     sam.NewSeq([]byte("GG")),
					Qual:    []byte{43, 43},
				},
			},
			want: []snp.BaseStrandPile{
				{
					Pos: 100000,
					Counts: [4][2]uint32{
						{1, 0},
						{0, 0},
						{0, 0},
						{0, 0},
					},
				},
				{
					Pos: 100001,
					Counts: [4][2]uint32{
						{0, 0},
						{1, 0},
						{0, 0},
						{0, 0},
					},
				},
				{
					Pos: 100002,
				},
				{
					Pos: 100003,
					Counts: [4][2]uint32{
						{0, 0},
						{0, 0},
						{1, 0},
						{0, 0},
					},
				},
				{
					Pos: 100004,
					Counts: [4][2]uint32{
						{0, 0},
						{0, 0},
						{1, 0},
						{0, 0},
					},
				},
				{
					Pos: 100006,
				},
			},
		},
		// Read-pair with some overlap, and only stitched bases pass the overlap.
		// Overlapped and non-overlapped N base included.  One stitching-conflict
		// included.
		{
			name: "stitch_overlap",
			reads: []sam.Record{
				{
					Name:  "read2",
					Ref:   ref,
					Pos:   100000,
					MapQ:  60,
					Cigar: []sam.CigarOp{sam.NewCigarOp(sam.CigarMatch, 6)},
					// - strand
					Flags:   sam.Paired | sam.ProperPair | sam.Reverse | sam.Read1,
					MateRef: ref,
					MatePos: 100002,
					Seq:     sam.NewSeq([]byte("ACNGGT")),
					Qual:    []byte{37, 37, 2, 25, 37, 37},
				},
				{
					Name:    "read2",
					Ref:     ref,
					Pos:     100002,
					MapQ:    60,
					Cigar:   []sam.CigarOp{sam.NewCigarOp(sam.CigarMatch, 6)},
					Flags:   sam.Paired | sam.ProperPair | sam.MateReverse | sam.Read2,
					MateRef: ref,
					MatePos: 100000,
					Seq:     sam.NewSeq([]byte("NGCTTA")),
					Qual:    []byte{2, 25, 25, 37, 37, 37},
				},
			},
			want: []snp.BaseStrandPile{
				{
					Pos: 100000,
				},
				{
					Pos: 100001,
				},
				{
					Pos: 100002,
					// No nonzero counts here since BaseStrandPile ignores Ns.
				},
				{
					Pos: 100003,
					Counts: [4][2]uint32{
						{0, 0},
						{0, 0},
						{0, 1},
						{0, 0},
					},
				},
				{
					Pos: 100004,
					// No nonzero counts here due to G/C stitching mismatch.
				},
				// No pos=100005 since BED file excludes that position.
				{
					Pos: 100006,
				},
			},
		},
		// Four reads with mates which were filtered out of the BAM.  One of them
		// doesn't pass the mapq threshold.  One of them contains an insertion and
		// a deletion; these don't currently appear in the pileup, but we'll
		// probably want to track them later.
		{
			name: "mates_removed",
			reads: []sam.Record{
				{
					Name:  "read3",
					Ref:   ref,
					Pos:   99999,
					MapQ:  60,
					Cigar: []sam.CigarOp{sam.NewCigarOp(sam.CigarMatch, 6)},
					// - strand
					Flags:   sam.Paired | sam.ProperPair | sam.Reverse | sam.Read1,
					MateRef: ref,
					MatePos: 100002,
					Seq:     sam.NewSeq([]byte("TACCGG")),
					Qual:    []byte{43, 43, 43, 43, 43, 37}, // last base fails base-qual
				},
				{
					Name:  "read4",
					Ref:   ref,
					Pos:   100001,
					MapQ:  60,
					Cigar: []sam.CigarOp{sam.NewCigarOp(sam.CigarMatch, 1), sam.NewCigarOp(sam.CigarInsertion, 1), sam.NewCigarOp(sam.CigarMatch, 2), sam.NewCigarOp(sam.CigarDeletion, 1), sam.NewCigarOp(sam.CigarMatch, 1)},
					// + strand
					Flags:   sam.Paired | sam.ProperPair | sam.MateReverse | sam.Read1,
					MateRef: ref,
					MatePos: 100002,
					Seq:     sam.NewSeq([]byte("CTCGT")),
					Qual:    []byte{43, 43, 43, 43, 43},
				},
				{
					Name:  "read5",
					Ref:   ref,
					Pos:   100003,
					MapQ:  60,
					Cigar: []sam.CigarOp{sam.NewCigarOp(sam.CigarMatch, 6)},
					// + strand
					Flags:   sam.Paired | sam.ProperPair | sam.Reverse | sam.Read2,
					MateRef: ref,
					MatePos: 100000,
					Seq:     sam.NewSeq([]byte("GGTCTT")),
					Qual:    []byte{43, 43, 43, 43, 43, 43},
				},
				{
					Name:    "read6",
					Ref:     ref,
					Pos:     100003,
					MapQ:    40, // fail
					Cigar:   []sam.CigarOp{sam.NewCigarOp(sam.CigarMatch, 6)},
					Flags:   sam.Paired | sam.ProperPair | sam.MateReverse | sam.Read2,
					MateRef: ref,
					MatePos: 100002,
					Seq:     sam.NewSeq([]byte("GGTCTT")),
					Qual:    []byte{43, 43, 43, 43, 43, 43},
				},
			},
			want: []snp.BaseStrandPile{
				{
					Pos: 100000,
					Counts: [4][2]uint32{
						{0, 1},
						{0, 0},
						{0, 0},
						{0, 0},
					},
				},
				{
					Pos: 100001,
					Counts: [4][2]uint32{
						{0, 0},
						{1, 1},
						{0, 0},
						{0, 0},
					},
				},
				{
					Pos: 100002,
					Counts: [4][2]uint32{
						{0, 0},
						{1, 1},
						{0, 0},
						{0, 0},
					},
				},
				{
					Pos: 100003,
					Counts: [4][2]uint32{
						{0, 0},
						{0, 0},
						{2, 1},
						{0, 0},
					},
				},
				{
					Pos: 100004,
					Counts: [4][2]uint32{
						{0, 0},
						{0, 0},
						{1, 0},
						{0, 0},
					},
				},
				{
					Pos: 100006,
					Counts: [4][2]uint32{
						{0, 0},
						{1, 0},
						{0, 0},
						{0, 0},
					},
				},
			},
		},
	}
	bampath := filepath.Join(tmpdir, "tmp.bam")
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			// Write a temporary .bam file containing just the given reads.
			out, err = file.Create(ctx, bampath)
			assert.NoError(t, err)
			// No "defer file.CloseAndReport" since we have to close this file before
			// we start writing the .gbai, and when we close a file twice, the second
			// call may produce an error.

			var bamWriter *bam.Writer
			bamWriter, err = bam.NewWriter(out.Writer(ctx), samHeader, 1)
			assert.NoError(t, err)
			for _, r := range tt.reads {
				err = bamWriter.Write(&r)
				assert.NoError(t, err)
			}
			err = bamWriter.Close()
			assert.NoError(t, err)
			err = out.Close(ctx)
			assert.NoError(t, err)

			// Write corresponding .gbai file.
			var inBam file.File
			inBam, err = file.Open(ctx, bampath)
			assert.NoError(t, err)
			defer file.CloseAndReport(ctx, inBam, &err)

			var gbai file.File
			gbai, err = file.Create(ctx, gbaipath)
			assert.NoError(t, err)
			// Similarly, no "defer file.CloseAndReport" here since we have to close
			// the gbai before calling Pileup().

			err = gbam.WriteGIndex(gbai.Writer(ctx), inBam.Reader(ctx), 1024, 1)
			assert.NoError(t, err)
			err = gbai.Close(ctx)
			assert.NoError(t, err)

			// Write pileup to temporary file in basestrand-rio format.
			err = snp.Pileup(ctx, bampath, filepath.Join("testdata", "chr2_subset.fa"), "basestrand-rio", outPrefix, &opts, nil)
			assert.NoError(t, err)

			// Verify output is as expected.
			var rio file.File
			rio, err = file.Open(ctx, filepath.Join(tmpdir, "bio-pileup.basestrand.rio"))
			assert.NoError(t, err)
			defer file.CloseAndReport(ctx, rio, &err)

			var unmarshaller snp.BaseStrandUnmarshaller
			scanner := recordio.NewScanner(rio.Reader(ctx), recordio.ScannerOpts{
				Unmarshal: unmarshaller.UnmarshalBaseStrand,
			})
			for _, want := range tt.want {
				assert.True(t, scanner.Scan())
				pileupRow := scanner.Get().(*snp.BaseStrandPile)
				assert.EQ(t, *pileupRow, want)
			}
			assert.False(t, scanner.Scan())
			assert.NoError(t, scanner.Err())
		})
		assert.NoError(t, err)
	}
}
