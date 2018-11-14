package bamprovider

import (
	"fmt"
	"runtime"
	"strings"

	"github.com/biogo/hts/sam"
	gbam "github.com/grailbio/bio/encoding/bam"
)

type pairIteratorSharedState struct {
	provider  Provider
	shardChan chan gbam.Shard // for receiving shard ranges through NewShardChannel.

	// "distantMates" store records whose mate are in different shards.
	distantMates *concurrentMap
}

// Pair encapsulates a pair of SAM records for a pair of reads, and whether
// any error was encountered in retrieving them.
type Pair struct {
	R1  *sam.Record
	R2  *sam.Record
	Err error
}

// MissingMateError is a specific error that can be used when one or more mates
// are missing.
type MissingMateError struct {
	Message string
}

func (mme MissingMateError) Error() string {
	return mme.Message
}

// PairIterator reads matched pairs of records from a BAM or PAM file. Use
// NewPairIterators to create an iterator.
type PairIterator struct {
	rec   Pair
	iter  Iterator
	shard gbam.Shard // Shard currently read

	shared            *pairIteratorSharedState
	localNameToRecord map[string]*sam.Record
}

// NewPairIterators creates a set of PairIterators.  A PairIterator yields pairs
// of records in the BAM or PAM data corresponding to primary alignments for
// paired reads.  Records will not be included if they represent secondary or
// supplemental alignments (based on SAM flags).  Pairs that have both reads
// unmapped will not be included unless includeUnmapped is true.
//
// The pairs in the BAM file will be randomly sharded across the PairIterators
// created by this function. Pairs are returned in an unspecified order, even
// within one PairIterator.
//
// Each PairIterator is thread-compatible. It is recommended to create one
// goroutine for each iterator.
func NewPairIterators(provider Provider, includeUnmapped bool) ([]*PairIterator, error) {
	parallelism := runtime.NumCPU()
	shards, err := provider.GenerateShards(GenerateShardsOpts{
		Strategy:            ByteBased,
		Padding:             0,
		IncludeUnmapped:     includeUnmapped,
		SplitUnmappedCoords: true})
	if err != nil {
		return nil, err
	}
	shared := &pairIteratorSharedState{
		provider:     provider,
		shardChan:    gbam.NewShardChannel(shards),
		distantMates: newConcurrentMap(),
	}
	iters := make([]*PairIterator, parallelism)
	for i := 0; i < parallelism; i++ {
		iters[i] = &PairIterator{
			shared:            shared,
			localNameToRecord: make(map[string]*sam.Record),
		}
	}
	return iters, nil
}

func isPrimary(record *sam.Record) bool {
	return (record.Flags&sam.Secondary) == 0 && (record.Flags&sam.Supplementary) == 0
}

// Record returns the current pair, or an error.
//
// REQUIRES: Scan() has been called and its last call returned true.
func (l *PairIterator) Record() Pair { return l.rec }

// Scan reads the next record. It returns true if a record has been read, and
// false on end of data stream.
func (l *PairIterator) Scan() bool {
	for {
		if l.iter == nil {
			// Start reading a new shard.
			var ok bool
			if l.shard, ok = <-l.shared.shardChan; !ok {
				break
			}
			l.iter = l.shared.provider.NewIterator(l.shard)
		}
		if l.iter.Scan() {
			record := l.iter.Record()
			if !isPrimary(record) {
				gbam.PutInFreePool(gbam.CastDown(record))
				continue
			}
			mate, ok := l.localNameToRecord[record.Name]
			if ok {
				// We've already seen the mate of this record.
				delete(l.localNameToRecord, record.Name)
				if record.Flags&sam.Read1 != 0 {
					l.rec = Pair{R1: record, R2: mate}
				} else {
					l.rec = Pair{R1: mate, R2: record}
				}
				return true
			}
			if mateInShard(record, &l.shard) {
				// Store the record for later, when we see its mate.
				l.localNameToRecord[record.Name] = record
				continue
			}
			// The reads in this pair are in different shards, so we have to synchronize
			// with other goroutines.
			mate = l.shared.distantMates.lookupAndDelete(record)
			if mate != nil {
				if record.Flags&sam.Read1 != 0 {
					l.rec = Pair{R1: record, R2: mate}
				} else {
					l.rec = Pair{R1: mate, R2: record}
				}
				return true
			}
			continue
		}
		// End of shard. Report records that didn't find a mate locally.
		if err := l.iter.Close(); err != nil {
			l.rec = Pair{Err: err}
			l.iter = nil
			return true
		}
		l.iter = nil
		var orphans []string
		if len(l.localNameToRecord) > 0 {
			for _, rec := range l.localNameToRecord {
				orphans = append(orphans, fmt.Sprintf("%v:[%v:%d,%v:%d]", rec.Name, rec.Ref.ID(), rec.Pos, rec.MateRef.ID(), rec.MatePos))
				if len(orphans) > 100 {
					break
				}
			}
		}
		l.localNameToRecord = make(map[string]*sam.Record, len(l.localNameToRecord))
		if len(orphans) > 0 {
			l.rec = Pair{Err: MissingMateError{fmt.Sprintf("shard %+v: didn't find expected mates for reads: %v", l.shard, strings.Join(orphans, "\n"))}}
			return true
		}
	}
	return false
}

// FinishPairIterators should be called after reading all pairs. It returns an
// error if there are some unpaired reads.
func FinishPairIterators(iters []*PairIterator) error {
	if len(iters) > 0 {
		// All iters have the same "shared" value, so just check the iters[0].
		n := iters[0].shared.distantMates.approxSize()
		if n > 0 {
			return MissingMateError{Message: fmt.Sprintf("found %d unmatched mates in the global hash", n)}
		}
	}
	return nil
}

// MateInShard checks if the mate of "record" is contained in the shard's range.
//
// If shard.StartSeq != 0 and the mate is at coordinate
// <shard.StartRef,shard.Start>, this function returns false, even though the
// real status may be true. The same holds for shard.EndSeq.  This function
// still guarantees that MateInShard(r1, shard1) == MateInShard(r2, shard2) for
// any matching pair of reads r1 (found in shard1) and r2 (found in shard2) .
func mateInShard(record *sam.Record, shard *gbam.Shard) bool {
	mateAddr := gbam.NewCoord(record.MateRef, record.MatePos, 0)
	startAddr := gbam.NewCoord(shard.StartRef, shard.Start, 0)
	endAddr := gbam.NewCoord(shard.EndRef, shard.End, 0)
	if mateAddr.LT(startAddr) {
		return false
	}
	if mateAddr.GE(endAddr) {
		return false
	}
	if shard.StartSeq != 0 && mateAddr.EQ(startAddr) {
		return false
	}
	return true
}
