package bampair

import (
	"fmt"
	"reflect"
	"sync"
	"time"

	"github.com/biogo/hts/sam"
	"github.com/grailbio/base/log"
	"github.com/grailbio/base/sync/multierror"
	"github.com/grailbio/bio/biopb"
	gbam "github.com/grailbio/bio/encoding/bam"
	"github.com/grailbio/bio/encoding/bamprovider"
	"grail.com/bio/encoding/bam"
)

const (
	collectionChanSize = 10000
)

// Opts contains the opts for GetDistantMates
type Opts struct {
	Parallelism int
	DiskShards  int
	ScratchDir  string
}

// RecordProcessor is a way for GetDistantMates to run Process() on
// every record in the bam file.  After a given shard invokes
// Process() on all the records in the shard, including the padding,
// the shard will invoke Close().
type RecordProcessor interface {
	Process(r *sam.Record)
	Close()
}

type shardInfoUpdate struct {
	shard           bam.Shard
	numStartPadding uint64 // number of reads in the start padding, including padding.
	numReads        uint64 // number of reads in the actual shard.
}

// GetDistantMates scans the BAM/PAM file given by provider, and then
// returns a DistantMateTable.  The caller must call
// DistantMateTable.Close() to release the resources used by the
// DistantMateTable.  GetDistantMates also returns a ShardInfo object
// that includes information like number of records in each shard.
// For an example of how to use GetDistantMates, see
// ExampleResolvePairs() in distant_mates_test.go.
func GetDistantMates(provider bamprovider.Provider, shardList []bam.Shard, opts *Opts,
	createProcessor func() RecordProcessor) (distantMates *DistantMateTable, shardInfo *ShardInfo, returnErr error) {
	log.Debug.Printf("scanning %d shards", len(shardList))
	shardChannel := gbam.NewShardChannel(shardList)

	// Scan all shards, looking for distant mates.
	header, returnErr := provider.GetHeader()
	if returnErr != nil {
		return
	}

	shardInfo = newShardInfo()
	for _, shard := range shardList {
		log.Debug.Printf("shard: %v", shard)
		shardInfo.add(&shard)
	}

	distantMates, returnErr = newDistantMateTable(header, opts.ScratchDir,
		opts.DiskShards, len(shardList))
	if returnErr != nil {
		return nil, nil, returnErr
	}
	var workerGroup sync.WaitGroup
	var distantMateGroup sync.WaitGroup
	collectionChannel := make(chan interface{}, collectionChanSize)
	distantMateGroup.Add(1)
	errs := multierror.NewMultiError(1)
	go func() {
		defer distantMateGroup.Done()
		errs.Add(collectDistantMates(collectionChannel, shardInfo, distantMates))
	}()
	log.Debug.Printf("Creating %d scanners for distant mates", opts.Parallelism)
	t0 := time.Now()

	for i := 0; i < opts.Parallelism; i++ {
		workerGroup.Add(1)
		log.Debug.Printf("Creating scanner %d", i)
		go func(worker int) {
			defer workerGroup.Done()
			var recordProcessor RecordProcessor
			if createProcessor != nil {
				recordProcessor = createProcessor()
			}
			errs.Add(findDistantMates(provider, worker, shardInfo, opts, recordProcessor,
				shardChannel, distantMates, collectionChannel))
		}(i)
	}
	workerGroup.Wait()
	close(collectionChannel)
	distantMateGroup.Wait()
	t1 := time.Now()
	if errs.ErrorOrNil() != nil {
		log.Debug.Printf("scanners failed in %v", t1.Sub(t0))
		return nil, nil, errs
	}

	log.Debug.Printf("scanners all done, found %d distant pairs in %v", distantMates.len(), t1.Sub(t0))
	if err := distantMates.finish(shardInfo); err != nil {
		return nil, nil, err
	}
	return
}

// IsLeftMost returns true for only one read from a pair.  LeftMost is
// defined by the read on the smaller reference id, the smaller
// alignment position, and if both refID and position are the same, R1
// is considered the LeftMost.
func IsLeftMost(r *sam.Record) bool {
	if r.Ref.ID() != r.MateRef.ID() {
		return r.Ref.ID() < r.MateRef.ID()
	}
	if r.Pos != r.MatePos {
		return r.Pos < r.MatePos
	}
	return (r.Flags & sam.Read1) != 0
}

func isDistantMate(shardInfo *ShardInfo, r *sam.Record) []int {
	// Determine if we must add r into distantMates.  Returns true if
	// either of the following conditions is true:
	//  1) Ref != MateRef
	//  2) r.MatePos is in the padding of shard S and r.Pos is
	//     outside of the padding of shard S.
	//
	// isDistantMate calculates this for every padded shard that
	// r might reside in, and adds those shard numbers to the return
	// value.  Therefore, the return value contains every shardIdx for
	// which r is a necessary as a distant mate to resolve a pair.
	//
	// isDistantMate assumes that r has a mapped mate.

	shards := []int{}
	mateShard := shardInfo.getMateShard(r)
	mateCoord := gbam.NewCoord(r.MateRef, r.MatePos, 0)
	potentials := []bam.Shard{}
	for shardIdx := mateShard.ShardIdx; shardIdx >= 0; shardIdx-- {
		info := shardInfo.GetInfoByIdx(shardIdx)
		if info.Shard.CoordInShard(info.Shard.Padding, mateCoord) {
			potentials = append(potentials, info.Shard)
			continue
		}
		break
	}
	for shardIdx := mateShard.ShardIdx + 1; shardIdx < shardInfo.Len(); shardIdx++ {
		info := shardInfo.GetInfoByIdx(shardIdx)
		if info.Shard.CoordInShard(info.Shard.Padding, mateCoord) {
			potentials = append(potentials, info.Shard)
			continue
		}
		break
	}
	for _, potential := range potentials {
		if r.RefID() != r.MateRef.ID() || !potential.RecordInPaddedShard(r) {
			shards = append(shards, potential.ShardIdx)
		}
	}
	return shards
}

func findDistantMates(provider bamprovider.Provider, worker int, shardInfo *ShardInfo,
	opts *Opts, processor RecordProcessor, channel chan bam.Shard,
	distantMates *DistantMateTable, collectionChannel chan interface{}) error {
	for shard := range channel {
		if shard.StartRef == nil {
			continue
		}
		log.Debug.Printf("worker %d scanning shard (%s,%d,%d)",
			worker, shard.StartRef.Name(), shard.Start, shard.End)

		iter := provider.NewIterator(shard)
		totalInStartPadding := uint64(0)
		fileIdx := uint64(0) // file index zeroed at the start of the actual shard.
		total := 0
		distant := 0
		lastCoord := biopb.Coord{RefId: biopb.InvalidRefID, Pos: biopb.InvalidPos}
		uniqueDistant := 0
		t0 := time.Now()
		for iter.Scan() {
			record := iter.Record()
			coord := gbam.CoordFromSAMRecord(record, 0)
			if lastCoord.GT(coord) {
				return fmt.Errorf("record were out of order last position %v, read %v", lastCoord, record)
			}
			lastCoord = coord

			if processor != nil {
				processor.Process(record)
			}

			if record.Pos < shard.Start {
				totalInStartPadding++
			}
			if shard.RecordInShard(record) {
				if (record.Flags & sam.Secondary) != 0 {
					log.Debug.Printf("Ignoring secondary read %s", record.Name)
				} else if (record.Flags & sam.Supplementary) != 0 {
					log.Debug.Printf("Ignoring supplementary read %s", record.Name)
				} else if (record.Flags & sam.Unmapped) != 0 {
					log.Debug.Printf("Ignoring unmapped read %s", record.Name)
				} else if (record.Flags & sam.MateUnmapped) != 0 {
					log.Debug.Printf("Ignoring singleton read %s", record.Name)
				} else if mateShards := isDistantMate(shardInfo, record); len(mateShards) > 0 {
					// The mate is in a different ref, or the mate's
					// padded shard does not contain this read, so save
					// this read to distantMates.
					uniqueDistant++
					for _, mateShardIdx := range mateShards {
						distant++
						log.Debug.Printf("saving distant mate for %s %d, diff ref %v, dist %v to shard %d",
							record.Name, record.Flags, record.Ref.ID() != record.MateRef.ID(),
							abs(record.Pos-record.MatePos), mateShardIdx)
						if err := distantMates.addDistantMate(mateShardIdx, record,
							fileIdx); err != nil {
							return err
						}
					}
					record = nil
				}
				fileIdx++
			}
			if record != nil {
				sam.PutInFreePool(record)
			}
			total++
		}

		if err := iter.Close(); err != nil {
			return fmt.Errorf("Close: %v", err)
		}
		collectionChannel <- shardInfoUpdate{shard, totalInStartPadding, fileIdx}
		if processor != nil {
			processor.Close()
		}

		t1 := time.Now()
		log.Debug.Printf("worker %d scanned shard %s read %d records, %d in shard, %d distant in %v",
			worker, shard.String(), total, fileIdx, distant, t1.Sub(t0))
	}
	log.Debug.Printf("worker %d done", worker)
	return nil
}

func collectDistantMates(collectionChannel chan interface{}, shardInfo *ShardInfo,
	distantMates *DistantMateTable) error {
	// Wait for all distant mates to arrive on collectionChannel, and put them
	// into the distantMates map.
	for item := range collectionChannel {
		switch item.(type) {
		case shardInfoUpdate:
			update := item.(shardInfoUpdate)
			shardInfo.updateInfoByShard(&update.shard, update.numStartPadding, update.numReads)
		default:
			return fmt.Errorf("Unexpected item type: %v", reflect.TypeOf(item))
		}
	}
	shardInfo.computeFileIndexes()
	return nil
}

func abs(x int) int {
	if x < 0 {
		return -x
	}
	return x
}
