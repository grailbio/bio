package bampair

import (
	"fmt"
	"io/ioutil"
	"os"
	"sync/atomic"

	"github.com/biogo/hts/sam"
)

type indexedRecord struct {
	r       *sam.Record
	fileIdx uint64
}

type mateShard interface {
	add(mate *sam.Record, fileIdx uint64) error
	closeWriter() error
	openReader() error
	getMate(shardInfo *ShardInfo, r *sam.Record) (*sam.Record, uint64)
	closeReader()
}

// DistantMateTable provides access to sam records.  It is indexed
// by shardIdx and read name.  The interface is designed so that the
// table can store mates either in memory or on disk.  It's intended
// use is to store read pair mates.
//
// Calls to addDistantMate() can occur concurrently, but the call to
// finishedAdding() should occur after all calls to addDistantMate()
// complete.  After calling finishedAdding(), any number of threads
// can call len(), openShard(), getMate(), and closeShard().
type DistantMateTable struct {
	header        *sam.Header
	tempDir       string
	numMateShards int
	inputShards   int
	entries       []mateShard
	shardInfo     *ShardInfo
	total         uint64
}

func newDistantMateTable(header *sam.Header, scratchDir string, numMateShards, numShards int) (*DistantMateTable, error) {
	d := &DistantMateTable{
		header:        header,
		numMateShards: numMateShards,
		inputShards:   numShards,
		entries:       make([]mateShard, numMateShards),
	}
	if numMateShards > 0 {
		tempDir, err := ioutil.TempDir(scratchDir, "markdups")
		if err != nil {
			return nil, fmt.Errorf("Could not create tmp dir in %s: %v", scratchDir, err)
		}
		d.tempDir = tempDir

		for mateShardIdx := 0; mateShardIdx < numMateShards; mateShardIdx++ {
			var err error
			d.entries[mateShardIdx], err = newDiskMateShard(header, d.tempDir, mateShardIdx, numMateShards)
			if err != nil {
				return nil, err
			}
		}
	} else {
		d.numMateShards = 1
		d.entries = []mateShard{newMemMateShard()}
	}
	return d, nil
}

// addDistantMate adds a record to the DistantMateTable.  Please be
// careful with parameters to this call; shardIdx is the shard index
// of read R, and mate is read R's mate.  fileIdx is the fileIdx of
// mate.
func (d *DistantMateTable) addDistantMate(shardIdx int, mate *sam.Record, fileIdx uint64) error {
	if d.shardInfo != nil {
		panic("Finished writing, no new mates allowed")
	}
	shardEntry := d.getShardEntry(shardIdx)
	if err := shardEntry.add(mate, fileIdx); err != nil {
		return err
	}
	atomic.AddUint64(&d.total, 1)
	return nil
}

func (d *DistantMateTable) finish(shardInfo *ShardInfo) error {
	d.shardInfo = shardInfo
	for _, shardEntry := range d.entries {
		if err := shardEntry.closeWriter(); err != nil {
			return err
		}
	}
	return nil
}

func (d *DistantMateTable) len() uint64 {
	return atomic.LoadUint64(&d.total)
}

// OpenShard prepares the shard, with the given shardIdx, to be
// queried with GetMate().
func (d *DistantMateTable) OpenShard(shardIdx int) error {
	if d.shardInfo == nil {
		panic("Called openShard before finishedAdding")
	}
	shardEntry := d.getShardEntry(shardIdx)
	return shardEntry.openReader()
}

// GetMate returns the mate of r, and also the mate's FileIdx (as
// computed using shardInfo and the mate's shard-relative FileIdx).
// The shardIdx argument is equal to the shardIdx of the shard where r
// resides.
func (d *DistantMateTable) GetMate(shardIdx int, r *sam.Record) (*sam.Record, uint64) {
	shardEntry := d.getShardEntry(shardIdx)
	return shardEntry.getMate(d.shardInfo, r)
}

// CloseShard closes the given shard so that further calls to
// GetMate() with the given shardIdx will fail.  CloseShard() frees
// resources that OpenShard() allocates.
func (d *DistantMateTable) CloseShard(shardIdx int) {
	shardEntry := d.getShardEntry(shardIdx)
	shardEntry.closeReader()
}

// Close frees resources taken by a DistantMateTable.  A user must
// call this after finishing with a DistantMateTable, and all shards
// have been closed with CloseShard().
func (d *DistantMateTable) Close() error {
	if len(d.tempDir) > 0 {
		return os.RemoveAll(d.tempDir)
	}
	return nil
}

func (d *DistantMateTable) getShardEntry(shardIdx int) mateShard {
	mateShardIdx := (shardIdx * d.numMateShards) / d.inputShards
	shardEntry := d.entries[mateShardIdx]
	return shardEntry
}
