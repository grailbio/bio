package bam_test

import (
	"bytes"
	"compress/gzip"
	"context"
	"fmt"
	"io"
	"testing"

	"github.com/grailbio/base/traverse"
	gbam "github.com/grailbio/bio/encoding/bam"
	"github.com/grailbio/hts/bam"
	"github.com/grailbio/hts/sam"
	"github.com/grailbio/testutil/expect"
)

var (
	chr8, _    = sam.NewReference("chr8", "", "", 2000000, nil, nil)
	chr9, _    = sam.NewReference("chr9", "", "", 3000000, nil, nil)
	header, _  = sam.NewHeader(nil, []*sam.Reference{chr8, chr9})
	headerb, _ = header.MarshalBinary()
	read1      = newRecord("ABCDEFG", chr8, 123, chr8, 456, sam.Read2)
	read1_R1   = newRecord("ABCDEFG", chr8, 123, chr8, 456, sam.Read1)
	read2      = newRecord("ABCDEFG", chr8, 456, chr8, 123, sam.Read1)
	read2_R2   = newRecord("ABCDEFG", chr8, 456, chr8, 123, sam.Read2)
	read3      = newRecord("foo", chr9, 777, chr9, 1000001, sam.Read2)
	read4      = newRecord("foo", chr9, 1000001, chr9, 777, sam.Read1)
	read5      = newRecord("dist999", chr9, 1000, chr9, 1999, sam.Read1)
	read6      = newRecord("dist999", chr9, 1999, chr9, 1000, sam.Read2)
	read7      = newRecord("dist1000", chr9, 2000, chr9, 3000, sam.Read1)
	read8      = newRecord("dist1000", chr9, 3000, chr9, 2000, sam.Read2)
	read9      = newRecord("dist1001", chr9, 3000, chr9, 4001, sam.Read1)
	read10     = newRecord("dist1001", chr9, 4001, chr9, 3000, sam.Read2)
	read11     = newRecord("HIJKLMNOP", chr8, 123, chr8, 457, sam.Read2)
	read12     = newRecord("HIJKLMNOP", chr8, 456, chr8, 123, sam.Read1)
	read13     = newRecord("A", chr8, 124, chr8, 457, sam.Read2)
	read14     = newRecord("A", chr8, 457, chr8, 124, sam.Read1)
	read15     = newRecord("B", chr8, 125, chr8, 458, sam.Read2)
	read16     = newRecord("B", chr8, 458, chr8, 125, sam.Read1)
	read17     = newRecord("C", chr8, 126, chr8, 459, sam.Read2)
	read18     = newRecord("C", chr8, 459, chr8, 126, sam.Read1)
	read19     = newRecord("C", chr8, 127, chr8, 460, sam.Read2)
	read20     = newRecord("C", chr8, 460, chr8, 127, sam.Read1)
	unmapped00 = newRecord("unmapped0", nil, -1, nil, -1, sam.Read1|sam.Unmapped|sam.MateUnmapped)
	unmapped01 = newRecord("unmapped0", nil, -1, nil, -1, sam.Read2|sam.Unmapped|sam.MateUnmapped)
)

func TestAdjacentShardedBAMReader(t *testing.T) {
	for _, tt := range []struct {
		name string
		recs []*sam.Record
		want []gbam.Pair
	}{
		{
			name: "mapped_pairs",
			recs: []*sam.Record{read1, read2, read3, read4, read7, read8},
			want: []gbam.Pair{{R1: read2, R2: read1}, {R1: read4, R2: read3}, {R1: read7, R2: read8}},
		},
		{
			name: "unmapped_pairs",
			recs: []*sam.Record{read1, read2, unmapped00, unmapped01},
			want: []gbam.Pair{{R1: read2, R2: read1}, {R1: unmapped00, R2: unmapped01}},
		},
	} {
		runSimpleTestCase(t, tt.name, tt.recs, tt.want)
	}
}

func TestConcurrentAdjacentShardedBAMReader(t *testing.T) {
	recs := []*sam.Record{
		read2, read1, // Pair 1
		read4, read3, // Pair 2
		read5, read6, // Pair 3
		read7, read8, // Pair 4
		read9, read10, // Pair 5
		read14, read13, // Pair 6
		unmapped00, unmapped01, // Pair 7
		read16, read15, // Pair 8
		read18, read17, // Pair 9
		read20, read19, // Pair 10
	}
	// 1 pair/shard
	runConcurrentTestCase(t, "concurrent_reads_1_pair_per_shard", 2, recs)
	// 3 pairs/shard (last shard will have 1 pair)
	runConcurrentTestCase(t, "concurrent_reads_3_pairs_per_shard", 6, recs)
	// 1 shard with all pairs (recordsPerShard >> 20)
	runConcurrentTestCase(t, "concurrent_reads_500_pairs_per_shard", 1000, recs)
}

func TestAdjacentShardedBAMReaderErrors(t *testing.T) {
	for _, tt := range []struct {
		name string
		recs []*sam.Record
	}{
		{
			name: "unordered_reads",
			recs: []*sam.Record{read1, read4, read3, read2, read7, read8},
		},
		{
			name: "odd_number_of_reads",
			recs: []*sam.Record{read1, read2, read3, read7, read8},
		},
		{
			name: "two_R2",
			recs: []*sam.Record{read1, read2_R2, read3, read4, read7, read8},
		},
		{
			name: "two_R1",
			recs: []*sam.Record{read1_R1, read2, read3, read4, read7, read8},
		},
		{
			name: "non_matching_positions",
			recs: []*sam.Record{read1, read2, read3, read4, read11, read12},
		},
	} {
		runErrorTestCase(t, tt.name, tt.recs)
	}
}

type marshalledPair struct {
	r1 []byte
	r2 []byte
}

func marshalPair(t *testing.T, testName string, pair gbam.Pair) marshalledPair {
	if pair.Err != nil {
		t.Fatalf("test %s: pair error: %s", testName, pair.Err)
	}
	var bufR1, bufR2 bytes.Buffer
	if err := bam.Marshal(pair.R1, &bufR1); err != nil {
		t.Fatalf("test %s: error marshalling record %s: %s", testName, pair.R1, err)
	}
	if err := bam.Marshal(pair.R2, &bufR2); err != nil {
		t.Fatalf("test %s: error marshalling record %s: %s", testName, pair.R2, err)
	}
	return marshalledPair{
		r1: bufR1.Bytes(),
		r2: bufR2.Bytes(),
	}
}

func marshalPairs(t *testing.T, testName string, pairs []gbam.Pair) []marshalledPair {
	marshalledPairs := make([]marshalledPair, len(pairs))
	for _, pair := range pairs {
		marshalledPairs = append(marshalledPairs, marshalPair(t, testName, pair))
	}
	return marshalledPairs
}

func newRecord(name string, ref *sam.Reference, pos int, mateRef *sam.Reference, matePos int, flags sam.Flags) *sam.Record {
	r := sam.GetFromFreePool()
	r.Name = name
	r.Ref = ref
	r.Pos = pos
	r.MateRef = mateRef
	r.MatePos = matePos
	r.Flags = flags
	return r
}

func runSimpleTestCase(t *testing.T, testName string, recs []*sam.Record, want []gbam.Pair) {
	ctx := context.Background()
	br := getShardedReader(t, ctx, testName, recs, len(recs), 1)
	shard := br.GetShard()

	got := make([]gbam.Pair, 0, len(want))
	for shard.Scan() {
		rec := shard.Record()
		got = append(got, rec)
		if rec.Err != nil {
			t.Fatalf("test %s: error while reading record: %s", testName, rec.Err)
		}
	}
	expect.EQ(t, marshalPairs(t, testName, got), marshalPairs(t, testName, want), "test %s", testName)

	// Ensure that we are not waiting on more shards.
	if br.GetShard() != nil {
		t.Errorf("test %s: still waiting on shards", testName)
	}

}

func runConcurrentTestCase(t *testing.T, testName string, recordsPerShard int, recs []*sam.Record) {
	ctx := context.Background()
	sbr := getShardedReader(t, ctx, testName, recs, recordsPerShard, 10)
	var wbuf bytes.Buffer
	sbw, err := gbam.NewShardedBAMWriter(&wbuf, gzip.DefaultCompression, 10, sbr.Header())
	if err != nil {
		t.Fatalf("test %s: error creating sharded bam writer", testName)
	}

	err = traverse.CPU(func() error {
		var (
			c        *gbam.ShardedBAMCompressor
			shardErr error
		)
		for {
			rshard := sbr.GetShard()
			if rshard == nil {
				break
			}
			c = sbw.GetCompressor()
			if shardErr = c.StartShard(rshard.ShardIdx); shardErr != nil {
				return fmt.Errorf("test %s: error while starting shard %d: %s", testName, rshard.ShardIdx, shardErr)
			}
			for rshard.Scan() {
				rec := rshard.Record()
				if shardErr = rec.Err; shardErr != nil {
					break
				}
				if shardErr = c.AddRecord(rec.R1); shardErr != nil {
					shardErr = fmt.Errorf("test %s: error while writing record %s to shard %d: %s", testName, rec.R1, rshard.ShardIdx, shardErr)
					break
				}
				if shardErr = c.AddRecord(rec.R2); shardErr != nil {
					shardErr = fmt.Errorf("test %s: error while writing record %s to shard %d: %s", testName, rec.R1, rshard.ShardIdx, shardErr)
					break
				}
			}
			if closeErr := c.CloseShard(); closeErr != nil {
				return fmt.Errorf("test %s: error while closing shard %d: %s", testName, rshard.ShardIdx, closeErr)
			} else if shardErr != nil {
				return shardErr
			}
		}
		return nil
	})
	if err != nil {
		t.Fatal(err)
	}

	if err = sbw.Close(); err != nil {
		t.Fatalf("test %s: error while closing sharded bam writer: %s", testName, err)
	}

	var gotReader *bam.Reader
	if gotReader, err = bam.NewReader(bytes.NewReader(wbuf.Bytes()), 1); err != nil {
		t.Fatalf("test %s: error while creating new bam reader: %s", testName, err)
	}
	for _, wantRec := range recs {
		var gotRec *sam.Record
		if gotRec, err = gotReader.Read(); err != nil {
			t.Errorf("test %s: error while reading read: %s", testName, err)
		}
		var gotBuf, wantBuf bytes.Buffer
		if err = bam.Marshal(gotRec, &gotBuf); err != nil {
			t.Fatalf("test %s: error marshalling record %s: %s", testName, gotRec, err)
		}
		if err = bam.Marshal(wantRec, &wantBuf); err != nil {
			t.Fatalf("test %s: error marshalling record %s: %s", testName, wantRec, err)
		}
		expect.EQ(t, wantBuf, gotBuf, "test %s", testName)
	}

	// Ensure that there are no more records.
	if _, err = gotReader.Read(); err != io.EOF {
		t.Errorf("test %s: extra records read", testName)
	}

	// Ensure that we are not waiting on more shards.
	if sbr.GetShard() != nil {
		t.Errorf("test %s: still waiting on shards", testName)
	}
}

func runErrorTestCase(t *testing.T, testName string, recs []*sam.Record) {
	ctx := context.Background()
	br := getShardedReader(t, ctx, testName, recs, 10, 1)
	shard := br.GetShard()

	for shard.Scan() {
		if pair := shard.Record(); pair.Err != nil {
			break
		}
	}

	if shard.Record().Err == nil {
		t.Errorf("test %s: expected error, but none found", testName)
	}

	// Ensure that we are not waiting on more shards.
	if br.GetShard() != nil {
		t.Errorf("test %s: still waiting on shards", testName)
	}
}

func getShardedReader(t *testing.T, ctx context.Context, testName string, recs []*sam.Record, recordsPerShard, queueSize int) *gbam.AdjacentShardedBAMReader {
	var buf bytes.Buffer
	bw, err := bam.NewWriter(&buf, header, 1)
	if err != nil {
		t.Fatalf("test %s: error creating bam writer: %s", testName, err)
	}
	for _, rec := range recs {
		if err = bw.Write(rec); err != nil {
			t.Fatalf("test %s: error writing rec %s: %s", testName, rec, err)
		}
	}
	if err = bw.Close(); err != nil {
		t.Fatalf("test %s: error closing bam writer: %s", testName, err)
	}
	br, err := gbam.NewAdjacentShardedBAMReader(ctx, &buf, recordsPerShard, queueSize)
	if err != nil {
		t.Fatalf("test %s: error creating bam reader: %s", testName, err)
	}
	gotHeaderBinary, err := br.Header().MarshalBinary()
	if err != nil {
		t.Fatalf("test %s: error marshalling BAM reader header: %s", testName, err)
	}
	expect.EQ(t, gotHeaderBinary, headerb, "test %s", testName)
	return br
}
