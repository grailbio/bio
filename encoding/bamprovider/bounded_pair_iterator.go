package bamprovider

import (
	"fmt"
	"runtime"

	gbam "github.com/grailbio/bio/encoding/bam"
	"github.com/grailbio/hts/sam"
)

const DefaultMaxPairSpan = 1000

// BoundedPairIterator is a deterministic alternative to PairIterator, usable
// in settings where distant mates (and read-pairs with unmapped read(s)) can
// be ignored.  ("Distant mates" are defined in the next comment.)
type BoundedPairIterator struct {
	rec      Pair
	iter     Iterator
	shard    gbam.Shard
	shardIdx int

	duplicateShardCrossers bool
	localNameToRecord      map[string]*sam.Record
}

// BoundedPairIteratorOpts controls the behavior of NewBoundedPairIterators()
// below.
// - There are two modes.
//   - In the default mode, each valid read-pair is returned by exactly one
//     iterator, dependent on min(r1Start, r2Start).  All read-pairs returned
//     by iterator 0 start before all read-pairs returned by iterator 1, etc.
//   - In the DuplicateShardCrossers=true mode, read-pairs which span multiple
//     shards are returned by all of those shard-iterators.  This facilitates
//     computation of position-based stats (e.g. bio-pileup).
// - For the purpose of these iterators, a read-pair has distant mates if
//   either (i) the reads are mapped to different chromosomes, or (ii)
//   max(r1End, r2End) - min(r1Start, r2Start) is greater than MaxPairSpan.  If
//   MaxPairSpan is zero, DefaultMaxPairSpan is used.  No read-pairs with
//   distant mates are returned.
// - The function attempts to return an iterator-slice of length
//   TargetParallelism.  If TargetParallelism is zero, runtime.NumCPU() is
//   used.  Occasionally, the returned slice will have a different length: e.g.
//   if we're working with a BAM with only one read-pair, there's no point in
//   subdividing it.
type BoundedPairIteratorOpts struct {
	DuplicateShardCrossers bool
	MaxPairSpan            int
	TargetParallelism      int
}

// NewBoundedPairIterators returns a slice of BoundedPairIterators covering the
// provider's BAM/PAM.
func NewBoundedPairIterators(provider Provider, opts BoundedPairIteratorOpts) (iters []*BoundedPairIterator, err error) {
	maxPairSpan := opts.MaxPairSpan
	if maxPairSpan == 0 {
		maxPairSpan = DefaultMaxPairSpan
	}
	targetParallelism := opts.TargetParallelism
	if targetParallelism == 0 {
		targetParallelism = runtime.NumCPU()
	}

	// GenerateShards() can return either fewer or more shards than NumShards;
	// hence "TargetParallelism" rather than "Parallelism".
	// Unlike most other GenerateShards() calls, we do not set NumShards larger
	// than targetParallelism: we want deterministic behavior here, so we
	// intentionally do not set up a shard channel with extra elements that can
	// be picked up by the first goroutines to finish.  Instead, we settle for
	// creating the most even byte-based shards we can.
	var shards []gbam.Shard
	if shards, err = provider.GenerateShards(GenerateShardsOpts{
		Strategy:  ByteBased,
		NumShards: targetParallelism,
		Padding:   maxPairSpan,
	}); err != nil {
		return
	}
	parallelism := len(shards)
	iters = make([]*BoundedPairIterator, parallelism)
	for i := 0; i < parallelism; i++ {
		iter := provider.NewIterator(shards[i])
		iters[i] = &BoundedPairIterator{
			iter:                   iter,
			shard:                  shards[i],
			shardIdx:               i,
			duplicateShardCrossers: opts.DuplicateShardCrossers,
			localNameToRecord:      make(map[string]*sam.Record),
		}
	}
	return
}

// Record returns the current pair, or an error.
//
// REQUIRES: Scan() has been called and its last call returned true.
func (bpi *BoundedPairIterator) Record() Pair { return bpi.rec }

func maxInt(a, b int) int {
	if a > b {
		return a
	}
	return b
}

// Scan reads the next record.  It returns true if a record has been read, and
// false on end of data stream.
func (bpi *BoundedPairIterator) Scan() bool {
	shard := bpi.shard
	// We must not include StartSeq/EndSeq in positional comparisons.  (Might
	// want to switch to a different type which doesn't have that field.)
	startAddr := gbam.NewCoord(shard.StartRef, shard.Start, 0)
	endAddr := gbam.NewCoord(shard.EndRef, shard.End, 0)
	duplicateShardCrossers := bpi.duplicateShardCrossers
	maxPairSpan := bpi.shard.Padding
	for bpi.iter.Scan() {
		record := bpi.iter.Record()
		if !isPrimary(record) {
			sam.PutInFreePool(record)
			continue
		}
		mate, ok := bpi.localNameToRecord[record.Name]
		if ok {
			delete(bpi.localNameToRecord, record.Name)

			// Sanity check.
			if record.MateRef != mate.Ref {
				bpi.rec = Pair{Err: fmt.Errorf("MateRef != mate's Ref for read-pair %s", record.Name)}
				return true
			}
			if record.MatePos != mate.Pos {
				bpi.rec = Pair{Err: fmt.Errorf("MatePos != mate's Pos for read-pair %s", record.Name)}
				return true
			}

			recordSpan, _ := record.Cigar.Lengths()
			mateSpan, _ := mate.Cigar.Lengths()
			endPos := maxInt(record.Pos+recordSpan, mate.Pos+mateSpan)
			if endPos-mate.Pos > maxPairSpan {
				sam.PutInFreePool(record)
				sam.PutInFreePool(mate)
				continue
			}

			if duplicateShardCrossers {
				// Start-position is in [shardStart - padding, shardEnd).  So the only
				// way this read-pair doesn't overlap the shard is if the end-position
				// <= shardStart.
				pairEndAddr := gbam.NewCoord(record.Ref, endPos, 0)
				if pairEndAddr.LE(startAddr) {
					sam.PutInFreePool(record)
					sam.PutInFreePool(mate)
					continue
				}
			}

			if record.Flags&sam.Read1 != 0 {
				bpi.rec = Pair{R1: record, R2: mate}
			} else {
				bpi.rec = Pair{R1: mate, R2: record}
			}
			return true
		}
		// Ignore read-pairs with distant mates.
		if (record.Ref != record.MateRef) || (record.MatePos > record.Pos+maxPairSpan) {
			sam.PutInFreePool(record)
			continue
		}
		// Ignore the read-pair if it's out-of-bounds, or where we already ignored
		// the first half.
		if record.MatePos < record.Pos {
			// If MateRef/MatePos are such that we *shouldn't* have ignored the first
			// half, this must be an invalid (for this use case) BAM/PAM where the
			// first half was filtered out.  Much better to error out immediately
			// than to risk running out of memory from localNameToRecord becoming
			// huge.
			mateAddr := gbam.NewCoord(record.MateRef, record.MatePos, 0)
			if mateAddr.GE(startAddr) && mateAddr.LT(endAddr) && (record.MatePos+maxPairSpan >= record.Pos) {
				bpi.rec = Pair{Err: fmt.Errorf("Missing mate for read-pair %s", record.Name)}
				sam.PutInFreePool(record)
				return true
			}
			sam.PutInFreePool(record)
			continue
		}
		recordAddr := gbam.NewCoord(record.Ref, record.Pos, 0)
		if recordAddr.GE(endAddr) || ((!duplicateShardCrossers) && recordAddr.LT(startAddr)) {
			sam.PutInFreePool(record)
			continue
		}
		bpi.localNameToRecord[record.Name] = record
	}
	err := bpi.iter.Close()
	bpi.iter = nil
	if err != nil {
		bpi.rec = Pair{Err: fmt.Errorf("Failed to close: %v", err)}
		return true
	}
	return false
}

// Shard returns the shard covered by this iterator.
// Note that this shard's StartSeq and EndSeq values are not guaranteed to be
// zero.
func (bpi *BoundedPairIterator) Shard() gbam.Shard {
	return bpi.shard
}

// FinishBoundedPairIterators should be called after reading all pairs, or on
// error-exit.  If any iterators are still open, this assumes error-exit
// occurred and just closes them.  Otherwise, it returns an error iff there are
// some unpaired reads.
func FinishBoundedPairIterators(iters []*BoundedPairIterator) error {
	nUnfinishedIter := 0
	nUnmatched := 0
	for _, bpi := range iters {
		if bpi.iter != nil {
			nUnfinishedIter++
			// Could bubble up this error.
			_ = bpi.iter.Close()
		} else {
			nUnmatched += len(bpi.localNameToRecord)
		}
	}
	if nUnfinishedIter != 0 {
		// Unmatched mate count is inaccurate in this case (there are some reads we
		// didn't iterate over), so don't report that.
		return fmt.Errorf("FinishBoundedPairIterators: %d unfinished iterator(s)", nUnfinishedIter)
	}
	if nUnmatched != 0 {
		return MissingMateError{Message: fmt.Sprintf("found %d unmatched mates", nUnmatched)}
	}
	return nil
}
