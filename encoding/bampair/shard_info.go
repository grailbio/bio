package bampair

import (
	"fmt"

	"github.com/biogo/store/llrb"
	gbam "github.com/grailbio/bio/encoding/bam"
	"github.com/grailbio/hts/sam"
	"grail.com/bio/encoding/bam"
)

type key struct {
	refID int
	start int
	info  *ShardInfoEntry
}

// Compare compares two key objects for use in llrb.
func (k key) Compare(c2 llrb.Comparable) int {
	k2 := c2.(key)
	if diff := k.refID - k2.refID; diff != 0 {
		return diff
	}
	return k.start - k2.start
}

// ShardInfoEntry contains handy information about a particular shard.
type ShardInfoEntry struct {
	Shard               bam.Shard // Shard is the bam.Shard object.
	NumStartPadding     uint64    // NumStartPadding is the number of reads in the start padding.
	NumReads            uint64    // NumReads is the number of reads in the actual shard.
	PaddingStartFileIdx uint64    // PaddingStartFileIdx is the FileIdx of the first read in the start padding.
	ShardStartFileIdx   uint64    // ShardStartFileIdx is the FileIdx of the first read in the shard (excluding the padding).
}

// ShardInfo contains handy information about all shards, and is
// indexed by both key object and shardIdx.
type ShardInfo struct {
	byKey   llrb.Tree               // byKey maps a key object to a ShardInfoEntry.
	byIndex map[int]*ShardInfoEntry // byIndex maps a shardIdx to a ShardInfoEntry.
}

func newShardInfo() *ShardInfo {
	i := ShardInfo{}
	i.byKey = llrb.Tree{}
	i.byIndex = make(map[int]*ShardInfoEntry)
	return &i
}

func (i *ShardInfo) add(shard *bam.Shard) {
	info := ShardInfoEntry{*shard, 0, 0, 0, 0}
	e := key{shard.StartRef.ID(), shard.Start, &info}
	i.byKey.Insert(e)
	i.byIndex[shard.ShardIdx] = &info
}

func (i *ShardInfo) getInfoByRecord(r *sam.Record) *ShardInfoEntry {
	k := key{r.Ref.ID(), r.Pos, nil}
	c := i.byKey.Floor(k)
	info := c.(key).info
	if !info.Shard.CoordInShard(0, gbam.CoordFromSAMRecord(r, 0)) {
		panic(fmt.Sprintf("Could not find shard for %v; found %+v instead", r, info.Shard))
	}
	return info
}

// GetInfoByShard returns the info for the given shard.
func (i *ShardInfo) GetInfoByShard(shard *bam.Shard) *ShardInfoEntry {
	k := key{shard.StartRef.ID(), shard.Start, nil}
	c := i.byKey.Get(k)
	if c != nil {
		return c.(key).info
	}
	return nil
}

// GetInfoByIdx returns the info for the given shard index..
func (i *ShardInfo) GetInfoByIdx(shardIdx int) *ShardInfoEntry {
	return i.byIndex[shardIdx]
}

func (i *ShardInfo) getMateShard(r *sam.Record) bam.Shard {
	k := key{r.MateRef.ID(), r.MatePos, nil}
	c := i.byKey.Floor(k)
	info := c.(key).info
	if !info.Shard.CoordInShard(0, gbam.NewCoord(r.MateRef, r.MatePos, 0)) {
		panic(fmt.Sprintf("Could not find shard for mate of %v; found %+v instead", r, info.Shard))
	}
	return info.Shard
}

func (i *ShardInfo) updateInfoByShard(shard *bam.Shard, numStartPadding, numReads uint64) {
	info := i.GetInfoByShard(shard)
	if info == nil {
		panic(fmt.Sprintf("Could not find info for shard: %v", shard))
	}
	info.NumStartPadding = numStartPadding
	info.NumReads = numReads
}

func (i *ShardInfo) getShardStartFileIdx(r *sam.Record) uint64 {
	info := i.getInfoByRecord(r)
	if info == nil {
		panic(fmt.Sprintf("Could not find shard info for record: %v", r))
	}
	return info.ShardStartFileIdx
}

// Len returns the number of shards in i.
func (i *ShardInfo) Len() int {
	return len(i.byIndex)
}

func (i *ShardInfo) computeFileIndexes() {
	readCount := uint64(0)
	for x := 0; x < len(i.byIndex); x++ {
		info := i.byIndex[x]
		info.ShardStartFileIdx = readCount
		info.PaddingStartFileIdx = info.ShardStartFileIdx - info.NumStartPadding
		readCount += info.NumReads
	}
}
