package bamprovider

import (
	"math"

	"github.com/grailbio/bio/biopb"
	gbam "github.com/grailbio/bio/encoding/bam"
	"github.com/grailbio/hts/sam"
)

// fakeProvider is only for unittests. It yields the given records.
type fakeProvider struct {
	header *sam.Header
	recs   []*sam.Record
}

type fakeIterator struct {
	recs []*sam.Record
	rec  *sam.Record

	shardRange    biopb.CoordRange
	addrGenerator gbam.CoordGenerator
}

// NewFakeProvider creates a provider that returns "header" in response to a
// GetHeader() call, and recs by GenerateShards+NewIterator calls
func NewFakeProvider(header *sam.Header, recs []*sam.Record) Provider {
	return &fakeProvider{header, recs}
}

// FileInfo implements the Provider interface.
func (b *fakeProvider) FileInfo() (FileInfo, error) {
	return FileInfo{}, nil
}

// GetHeader implements the Provider interface. It returns the header passed to
// the constructor.
func (b *fakeProvider) GetHeader() (*sam.Header, error) {
	return b.header, nil
}

// Close implements the Provider interface.
func (b *fakeProvider) Close() error {
	return nil
}

func (b *fakeProvider) GetFileShards() ([]gbam.Shard, error) {
	return []gbam.Shard{gbam.UniversalShard(b.header)}, nil
}

// GenerateShards implements the Provider interface.
func (b *fakeProvider) GenerateShards(opts GenerateShardsOpts) ([]gbam.Shard, error) {
	shards := []gbam.Shard{gbam.Shard{
		StartRef: b.header.Refs()[0],
		Start:    0,
		EndRef:   nil,
		End:      0,
		Padding:  opts.Padding,
	}}

	if !opts.IncludeUnmapped {
		// TODO(saito): create more complex shards for mapped reads.
		return shards, nil
	}

	if opts.SplitUnmappedCoords {
		// Create two shards for unmapped reads. The first shard converts up to the
		// first two unmapped reads. The second shard covers the rest.
		shards = append(shards,
			gbam.Shard{
				StartRef: nil,
				Start:    0,
				EndRef:   nil,
				End:      0,
				EndSeq:   2,
				ShardIdx: 1,
			},
			gbam.Shard{
				StartRef: nil,
				Start:    0,
				StartSeq: 2,
				EndRef:   nil,
				End:      math.MaxInt32,
				ShardIdx: 2,
			})
	} else {
		shards = append(shards,
			gbam.Shard{
				StartRef: nil,
				Start:    0,
				EndRef:   nil,
				End:      math.MaxInt32,
				ShardIdx: 1,
			})
	}
	return shards, nil
}

// NewIterator implements the Provider interface.
//
// REQUIRES: shard must be the one created by GenerateShards.
func (b *fakeProvider) NewIterator(shard gbam.Shard) Iterator {
	return &fakeIterator{recs: b.recs, rec: nil,
		shardRange: biopb.CoordRange{
			biopb.Coord{int32(shard.StartRef.ID()), int32(shard.PaddedStart()), int32(shard.StartSeq)},
			biopb.Coord{int32(shard.EndRef.ID()), int32(shard.PaddedEnd()), int32(shard.EndSeq)},
		}}
}

// Err implements the Iterator interface.
func (i *fakeIterator) Err() error {
	return nil
}

// Close implements the Iterator interface.
func (i *fakeIterator) Close() error {
	return nil
}

func (i *fakeIterator) Scan() bool {
	for {
		if len(i.recs) == 0 {
			return false
		}
		i.rec = i.recs[0]
		i.recs = i.recs[1:]
		addr := i.addrGenerator.GenerateFromRecord(i.rec)
		if i.shardRange.Contains(addr) {
			return true
		}
	}
}

func (i *fakeIterator) Record() *sam.Record {
	// Return a copy so that the code under test cannot alter the
	// original test input data.
	copy := sam.GetFromFreePool()
	*copy = *i.rec
	return copy
}
