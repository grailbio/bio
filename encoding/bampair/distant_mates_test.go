package bampair

import (
	"fmt"
	"os"
	"testing"

	"github.com/biogo/hts/sam"
	"github.com/grailbio/base/grail"
	gbam "github.com/grailbio/bio/encoding/bam"
	"github.com/grailbio/bio/encoding/bamprovider"
	"github.com/grailbio/testutil"
	"github.com/stretchr/testify/assert"
	grailbam "grail.com/bio/encoding/bam"
)

var (
	chr1, _   = sam.NewReference("chr1", "", "", 1000, nil, nil)
	chr2, _   = sam.NewReference("chr2", "", "", 1000, nil, nil)
	header, _ = sam.NewHeader(nil, []*sam.Reference{chr1, chr2})

	cigar1M = []sam.CigarOp{
		sam.NewCigarOp(sam.CigarMatch, 1),
	}
	cigar0 = []sam.CigarOp{
		sam.NewCigarOp(sam.CigarMatch, 10),
	}
	cigarSoft1 = []sam.CigarOp{
		sam.NewCigarOp(sam.CigarSoftClipped, 1),
		sam.NewCigarOp(sam.CigarMatch, 8),
		sam.NewCigarOp(sam.CigarSoftClipped, 1),
	}
	cigarHard1 = []sam.CigarOp{
		sam.NewCigarOp(sam.CigarHardClipped, 1),
		sam.NewCigarOp(sam.CigarMatch, 8),
		sam.NewCigarOp(sam.CigarHardClipped, 1),
	}
	r1F = sam.Paired | sam.Read1
	r2F = sam.Paired | sam.Read2

	// Reads for testing duplicate marking.
	// The following duplicate group (basic) is entirely in the same shard.
	basicA1 = newRecord("A:::1:10:1:1", chr1, 0, r1F, 10, chr1, cigar0)
	basicA2 = newRecord("A:::1:10:1:1", chr1, 10, r2F, 0, chr1, cigar0)
	// Exact duplicate, no clipping
	basicB1 = newRecord("B:::1:10:2:2", chr1, 0, r1F, 10, chr1, cigar0)
	basicB2 = newRecord("B:::1:10:2:2", chr1, 10, r2F, 0, chr1, cigar0)
	// Duplicate with soft clipping
	basicC1 = newRecord("C:::1:10:3:3", chr1, 1, r1F, 11, chr1, cigarSoft1)
	basicC2 = newRecord("C:::1:10:3:3", chr1, 11, r2F, 1, chr1, cigarSoft1)
	// Duplicate with hard clipping
	basicD1 = newRecord("D:::1:10:4:4", chr1, 1, r1F, 11, chr1, cigarHard1)
	basicD2 = newRecord("D:::1:10:4:4", chr1, 11, r2F, 1, chr1, cigarHard1)

	// distantA2 is in shard1, and distantA1 is outside the padding of shard1.
	//   -> distantA1 should be in distantMates.
	distantA1 = newRecord("A", chr1, 0, sam.Read1, 100, chr1, cigar0)
	distantA2 = newRecord("A", chr1, 100, sam.Read1, 0, chr1, cigar0)

	// Reads for testing distant mates:

	// distantB1 is inside the padding of shard0, and distantB2
	// is outside the padding of shard0.
	//   -> distantB2 should be in distantMates.
	distantB1 = newRecord("B", chr1, 100, sam.Read1, 115, chr1, cigar0)
	distantB2 = newRecord("B", chr1, 115, sam.Read1, 100, chr1, cigar0)

	// distantC1 is in shard0, and distantC2 is outside the padding of shard0.
	//   -> distantC2 should be in distantMates.
	// distantC2 is in shard1, and distantC1 is outside the padding of shard1.
	//   -> distantC1 should be in distantMates.
	distantC1 = newRecord("C", chr1, 50, sam.Read1, 150, chr1, cigar0)
	distantC2 = newRecord("C", chr1, 150, sam.Read1, 50, chr1, cigar0)

	// Neither distantD1 nor distantD2 should both be in distantMates.
	// They are both in the padding of shard0.  They are both in shard1.
	distantD1 = newRecord("D", chr1, 103, sam.Read1, 104, chr1, cigar0)
	distantD2 = newRecord("D", chr1, 104, sam.Read1, 103, chr1, cigar0)

	// Test when padding is larger than the shard size.  In this case,
	// distant mates should place R2 as a distant mate for 5 different
	// shards.  R1 will also be a distant mate for 5 different shards.
	distantF1 = newRecord("F", chr1, 12, sam.Read1, 37, chr1, cigar1M)
	distantF2 = newRecord("F", chr1, 37, sam.Read2, 12, chr1, cigar1M)
)

func newRecord(name string, ref *sam.Reference, pos int, flags sam.Flags, matePos int, mateRef *sam.Reference, cigar sam.Cigar) *sam.Record {
	r := gbam.CastUp(gbam.GetFromFreePool())
	r.Name = name
	r.Ref = ref
	r.Pos = pos
	r.MatePos = matePos
	r.MateRef = mateRef
	r.Flags = flags
	r.Cigar = cigar
	return r
}

func TestGetDistantMates(t *testing.T) {
	// This test sets diskMateShards to 1000 to ensure that each mate
	// shard is kept separate from the other mate shards, so any bugs
	// that mix up mate shards will cause a test failure.

	type distantEntry struct {
		shardIdx int
		mate     *sam.Record
		r        *sam.Record
		fileIdx  uint64
	}

	tests := []struct {
		records         []*sam.Record
		expectedDistant []*distantEntry
		shardSize       int
		padding         int
		perShardReads   map[string]uint64
	}{
		{
			[]*sam.Record{basicA1, basicB1, basicC1, basicD1, basicA2, basicB2, basicC2, basicD2},
			[]*distantEntry{},
			100, 5,
			map[string]uint64{
				"chr1:0:100": 8,
			},
		},
		{
			[]*sam.Record{distantA1, distantA2},
			[]*distantEntry{{1, distantA1, distantA2, 0}},
			100, 5,
			map[string]uint64{
				"chr1:0:100":   1,
				"chr1:100:200": 1,
			},
		},
		{
			[]*sam.Record{distantB1, distantB2},
			[]*distantEntry{{0, distantB2, distantB1, 1}},
			100, 5,
			map[string]uint64{
				"chr1:100:200": 2,
			},
		},
		{
			[]*sam.Record{distantC1, distantC2},
			[]*distantEntry{{1, distantC1, distantC2, 0}, {0, distantC2, distantC1, 1}},
			100, 5,
			map[string]uint64{
				"chr1:0:100":   1,
				"chr1:100:200": 1,
			},
		},
		{
			[]*sam.Record{distantD1, distantD2},
			[]*distantEntry{},
			100, 5,
			map[string]uint64{
				"chr1:100:200": 2,
			},
		},
		{
			[]*sam.Record{distantF1, distantF2},
			[]*distantEntry{
				{0, distantF2, distantF1, 1},
				{1, distantF2, distantF1, 1},
				{2, distantF2, distantF1, 1},
				{3, distantF2, distantF1, 1},
				{4, distantF2, distantF1, 1},
				{5, distantF1, distantF2, 0},
				{6, distantF1, distantF2, 0},
				{7, distantF1, distantF2, 0},
				{8, distantF1, distantF2, 0},
				{9, distantF1, distantF2, 0},
			},
			5, 10,
			map[string]uint64{
				"chr1:10:15": 1,
				"chr1:35:40": 1,
			},
		},
	}

	for i, test := range tests {
		t.Logf("--- starting tests[%d] ---", i)
		tempDir, cleanup := testutil.TempDir(t, "", "")
		defer cleanup()
		provider := bamprovider.NewFakeProvider(header, test.records)

		opts := Opts{
			Parallelism: 1,
			DiskShards:  1000,
			ScratchDir:  tempDir,
		}

		providerHeader, err := provider.GetHeader()
		assert.NoError(t, err)
		shardList, err := gbam.GetPositionBasedShards(providerHeader, test.shardSize, test.padding, true)
		assert.NoError(t, err)
		distantMates, allShards, err := GetDistantMates(provider, shardList, &opts, nil)
		if err != nil {
			t.Errorf("error from GetDistantMates: %v", err)
		}
		for shardIdx := 0; shardIdx < len(shardList); shardIdx++ {
			distantMates.OpenShard(shardIdx)
		}
		assert.Equal(t, len(test.expectedDistant), int(distantMates.len()))
		for _, expected := range test.expectedDistant {
			actualMate, actualFileIdx := distantMates.GetMate(expected.shardIdx, expected.r)
			if actualMate != nil {
				assert.Equal(t, expected.mate.String(), actualMate.String())
				assert.Equal(t, expected.fileIdx, actualFileIdx)
			} else {
				t.Errorf("could not find mate %v in distantMates", expected.mate)
			}
		}

		// Check allShards output.
		for i := 0; i < allShards.Len(); i++ {
			info := allShards.GetInfoByIdx(i)

			// Shards must arrive in order, and have the right number of reads.
			assert.Equal(t, i, info.Shard.ShardIdx)

			k := fmt.Sprintf("%s:%d:%d", info.Shard.StartRef.Name(), info.Shard.Start, info.Shard.End)
			expectedReadCount, ok := test.perShardReads[k]
			if ok {
				assert.Equal(t, expectedReadCount, info.NumReads)
			} else {
				assert.Equal(t, uint64(0), info.NumReads)
			}
		}
		for shardIdx := 0; shardIdx < len(shardList); shardIdx++ {
			distantMates.CloseShard(shardIdx)
		}
	}
}

func ExampleResolvePairs() {
	testRecords := []*sam.Record{
		newRecord("A:::1:10:1:1", chr1, 0, r1F, 10, chr1, cigar0), // near mates
		newRecord("A:::1:10:1:1", chr1, 10, r2F, 0, chr1, cigar0),
		newRecord("B:::1:10:2:2", chr1, 120, r1F, 210, chr1, cigar0), // distant mates
		newRecord("B:::1:10:2:2", chr1, 210, r2F, 120, chr1, cigar0),
	}
	provider := bamprovider.NewFakeProvider(header, testRecords)
	providerHeader, err := provider.GetHeader()
	if err != nil {
		panic(err)
	}
	shards, err := gbam.GetPositionBasedShards(providerHeader, 100, 0, true)
	if err != nil {
		panic(err)
	}

	// Scan the BAM/PAM file to identify and save the distant mates.
	opts := &Opts{Parallelism: 3}
	distantMates, _, err := GetDistantMates(provider, shards, opts, nil)
	if err != nil {
		panic(err)
	}

	nearMates := 0
	farMates := 0
	for shardIdx, shard := range shards {
		// Open mates for reading.
		distantMates.OpenShard(shardIdx)

		// Read all the records
		records := []*sam.Record{}
		byName := map[string][]*sam.Record{}
		iter := provider.NewIterator(shard)
		for iter.Scan() {
			r := iter.Record()

			// Ignore secondary, supplementary, unmapped records, and
			// records with unmapped mates.
			if (r.Flags&sam.Secondary) != 0 || (r.Flags&sam.Supplementary) != 0 {
				continue
			} else if (r.Flags&sam.Unmapped) != 0 || grailbam.HasNoMappedMate(r) {
				continue
			}

			records = append(records, r)
			byName[r.Name] = append(byName[r.Name], r)
		}
		if err := iter.Close(); err != nil {
			panic(err)
		}

		// Scan through all mapped pairs in original order and resolve
		// mates.  Use IsLeftMost() to process each pair only once.
		for _, r := range records {
			if !IsLeftMost(r) {
				continue
			}

			pair := byName[r.Name]
			if len(pair) == 1 {
				mate, _ := distantMates.GetMate(shardIdx, r)
				// Process the read pair here: r, mate
				if mate != nil {
					farMates++
				}
			} else if len(pair) == 2 {
				// Process the read pair here: pair[0], pair[1]
				nearMates++
			} else {
				panic(fmt.Sprintf("unexpected number of records for %s: %d", r.Name, len(pair)))
			}
		}

		// Close mates.
		distantMates.CloseShard(shardIdx)
	}
	distantMates.Close()
	fmt.Printf("nearMates: %d\n", nearMates)
	fmt.Printf("farMates: %d\n", farMates)

	// Output:
	// nearMates: 1
	// farMates: 1
}

func TestMain(m *testing.M) {
	shutdown := grail.Init()
	defer shutdown()
	os.Exit(m.Run())
}
