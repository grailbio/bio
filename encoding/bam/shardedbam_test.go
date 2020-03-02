package bam_test

import (
	"bytes"
	"io"
	"os"
	"testing"

	gbam "github.com/grailbio/bio/encoding/bam"
	"github.com/grailbio/hts/bam"
	"github.com/grailbio/hts/sam"
	"github.com/grailbio/testutil"
	"github.com/grailbio/testutil/assert"
	"github.com/grailbio/testutil/expect"
	"github.com/klauspost/compress/gzip"
)

func verifyBAM(t *testing.T, records []*sam.Record, bamBuffer *bytes.Buffer) {
	reader, err := bam.NewReader(bamBuffer, 1)
	expect.Nil(t, err)
	i := 0
	for {
		r, err := reader.Read()
		if err == io.EOF {
			expect.EQ(t, i, len(records), "not enough records in bam output %d vs %d", len(records), i)
			break
		}
		expected, err := records[i].MarshalText()
		expect.Nil(t, err)
		actual, err := r.MarshalText()
		expect.Nil(t, err)

		expect.EQ(t, actual, expected, "record[%d] does not match %v vs %v", i, records[i], r)
		i++
	}
}

func writeAndVerify(t *testing.T, header *sam.Header, records []*sam.Record, compressors, shards int, forward bool) {
	var bamBuffer bytes.Buffer
	w, err := gbam.NewShardedBAMWriter(&bamBuffer, gzip.DefaultCompression, 10, header)
	if err != nil {
		t.Errorf("error creating ShardedBAMWriter: %v", err)
	}

	shardCompressors := make([]*gbam.ShardedBAMCompressor, compressors)
	for i := 0; i < compressors; i++ {
		shardCompressors[i] = w.GetCompressor()
	}

	for i := 0; i < shards; i++ {
		var shardNum int
		if forward {
			shardNum = i
		} else {
			shardNum = shards - 1 - i
		}

		c := shardNum % compressors
		err := shardCompressors[c].StartShard(shardNum)
		assert.Nil(t, err)

		for _, r := range records[shardNum*(len(records)/shards) : (shardNum+1)*(len(records)/shards)] {
			err := shardCompressors[c].AddRecord(r)
			expect.Nil(t, err)
		}
		// If there are remainders from uneven division, then add them to the last shard.
		if shardNum == shards-1 {
			for _, r := range records[(shardNum+1)*(len(records)/shards):] {
				err := shardCompressors[c].AddRecord(r)
				expect.Nil(t, err)
			}
		}

		err = shardCompressors[c].CloseShard()
		expect.Nil(t, err)
	}
	err = w.Close()
	expect.Nil(t, err)
	verifyBAM(t, records, &bamBuffer)
}

func TestShardedBAMSmall(t *testing.T) {
	chr1, err := sam.NewReference("chr1", "", "", 1000, nil, nil)
	expect.Nil(t, err)
	chr2, err := sam.NewReference("chr2", "", "", 2000, nil, nil)
	expect.Nil(t, err)
	header, err := sam.NewHeader(nil, []*sam.Reference{chr1, chr2})
	expect.Nil(t, err)
	cigar := []sam.CigarOp{
		sam.NewCigarOp(sam.CigarSoftClipped, 1),
		sam.NewCigarOp(sam.CigarMatch, 8),
		sam.NewCigarOp(sam.CigarSoftClipped, 1),
	}

	r1F := sam.Paired | sam.Read1
	r2R := sam.Paired | sam.Read2 | sam.Reverse

	records := []*sam.Record{
		&sam.Record{Name: "A::::10:1:1", Ref: chr1, Pos: 0, Flags: r1F, MatePos: 10, MateRef: chr1, Cigar: cigar},
		&sam.Record{Name: "A::::10:1:1", Ref: chr1, Pos: 10, Flags: r2R, MatePos: 0, MateRef: chr1, Cigar: cigar},
		&sam.Record{Name: "B::::10:1:1", Ref: chr1, Pos: 20, Flags: r1F, MatePos: 30, MateRef: chr1, Cigar: cigar},
		&sam.Record{Name: "B::::10:1:1", Ref: chr1, Pos: 30, Flags: r2R, MatePos: 20, MateRef: chr1, Cigar: cigar},
		&sam.Record{Name: "C::::10:1:1", Ref: chr1, Pos: 40, Flags: r1F, MatePos: 50, MateRef: chr1, Cigar: cigar},
		&sam.Record{Name: "C::::10:1:1", Ref: chr1, Pos: 50, Flags: r2R, MatePos: 40, MateRef: chr1, Cigar: cigar},
		&sam.Record{Name: "D::::10:1:1", Ref: chr2, Pos: 60, Flags: r1F, MatePos: 70, MateRef: chr2, Cigar: cigar},
		&sam.Record{Name: "D::::10:1:1", Ref: chr2, Pos: 70, Flags: r2R, MatePos: 60, MateRef: chr2, Cigar: cigar},
		&sam.Record{Name: "E::::10:1:1", Ref: chr2, Pos: 80, Flags: r1F, MatePos: 90, MateRef: chr2, Cigar: cigar},
		&sam.Record{Name: "E::::10:1:1", Ref: chr2, Pos: 90, Flags: r2R, MatePos: 80, MateRef: chr2, Cigar: cigar},
		&sam.Record{Name: "F::::10:1:1", Ref: chr2, Pos: 100, Flags: r1F, MatePos: 110, MateRef: chr2, Cigar: cigar},
		&sam.Record{Name: "F::::10:1:1", Ref: chr2, Pos: 110, Flags: r2R, MatePos: 100, MateRef: chr2, Cigar: cigar},
	}
	writeAndVerify(t, header, records, 2, 3, true)
	writeAndVerify(t, header, records, 2, 3, false)
}

func TestShardedBAMLarge(t *testing.T) {
	filename := testutil.GetFilePath("//go/src/grail.com/bio/encoding/bam/testdata/170614_WGS_LOD_Pre_Library_B3_27961B_05.merged.10000.bam")
	f, err := os.Open(filename)
	expect.Nil(t, err)
	reader, err := bam.NewReader(f, 1)
	expect.Nil(t, err)

	records := []*sam.Record{}
	for {
		record, err := reader.Read()
		if err == io.EOF {
			break
		}
		records = append(records, record)
	}

	writeAndVerify(t, reader.Header(), records, 3, 6, true)
	writeAndVerify(t, reader.Header(), records, 3, 6, false)
}
