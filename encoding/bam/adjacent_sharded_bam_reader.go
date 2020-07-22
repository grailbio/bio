package bam

import (
	"context"
	"fmt"
	"io"

	"github.com/grailbio/base/errors"
	"github.com/grailbio/hts/bam"
	"github.com/grailbio/hts/sam"
)

// Pair encapsulates a pair of SAM records for a pair of reads, and whether
// any error was encountered in retrieving them.
type Pair struct {
	R1  *sam.Record
	R2  *sam.Record
	Err error
}

// AdjacentShardedBAMReader provides a deterministic way
// to read a BAM file, provided that all records are grouped
// into adjacent pairs. Each AdjacentBAMShard has sequentially
// increasing shard numbers starting at zero. The records in
// each AdjacentBAMShard are in the same order as they appear
// in the underlying BAM file.
//
// AdjacentBAMShards are returned one at a time by GetShard()
// as soon as a shard is ready. Shards are created
// automatically in a goroutine when
// NewAdjacentShardedBAMReader() is called.
// Each AdjacentBAMShard is thread-safe.
//
// When records are adjacent, there is no BAM index so the
// the number of shards is indeterminate. Because of this,
// it is recommended that the caller limits the number of
// goroutines used for reading a BAM file. In this example,
// the number of goroutines is limited to the number of
// available CPUs. If there are more shards that can be
// concurrently processed than CPUs, there will be multiple
// shards sequentially processed in the same goroutine.
//
// Example Use of AdjacentShardedBAMReader:
//
//   ctx := context.Background()
//   f, _ := os.Create("input.bam")
//   r, _ := NewAdjacentShardedBAMReader(ctx, f, 100000, 2)
//
//   err = traverse.CPU(func() error {
//      for {
//         shard := r.GetShard()
//         if shard == nil { break }
//         for shard.Scan() {
//               pair := shard.Record()
//               if pair.Err != nil { return pair.Err }
//               // Do something with the record and
//               // and use shard.ShardIdx to denote
//               // the order of the shards.
//         }
//      }
//      return nil
//   })
//
type AdjacentShardedBAMReader struct {
	// shardc is the channel to which shards
	// are pushed to once they are created.
	// The caller obtains shards for reading
	// by calling and waiting on GetShard().
	shardc chan *AdjacentBAMShard
	// header is the header of the BAM
	// file being read.
	header *sam.Header
	// errors is used to halt the reading
	// of the underlying BAM file and the
	// scanning of all shards if the BAM
	// reader or any shard encounters an
	// error.
	errors *errors.Once
}

// NewAdjacentShardedBAMReader returns a new AdjacentShardedBAMReader
// that allows for the concurrent reading of a BAM file with adjacent
// paired records.
func NewAdjacentShardedBAMReader(ctx context.Context, r io.Reader, recordsPerShard int, queueSize int) (*AdjacentShardedBAMReader, error) {
	if recordsPerShard%2 != 0 {
		panic("shards must have an even number of records: each record must belong to a pair")
	}
	br, err := bam.NewReader(r, 1)
	if err != nil {
		return nil, err
	}
	sbr := &AdjacentShardedBAMReader{
		shardc: make(chan *AdjacentBAMShard, queueSize),
		errors: new(errors.Once),
		header: br.Header(),
	}

	go func() {
		defer func() {
			close(sbr.shardc)
			_ = br.Close()
		}()

		var shardIdx int
		shard := &AdjacentBAMShard{
			ShardIdx: shardIdx,
			recs:     newRecBuffer(recordsPerShard),
			errors:   sbr.errors,
		}
		// Iterate through each record and put each
		// record into a shard. By iterating through
		// records sequentially and sequentially adding
		// them to a shard, we preserve the order of records
		// within a shard as well as the order of shards
		// relative to each other.
		for {
			select {
			case <-ctx.Done():
				sbr.errors.Set(fmt.Errorf("error while reading BAM file: %s", ctx.Err()))
				return
			default:
			}
			// A failure in reading one shard
			// must prevent the reading of other
			// records.
			if sbr.errors.Err() != nil {
				return
			}

			rec, err := br.Read()
			if rec == nil {
				if err != io.EOF {
					sbr.errors.Set(fmt.Errorf("error while reading BAM file: %s", err))
					return
				}
				// Finish off the last shard, if it contains reads.
				if shard.recs.numRecords() > 0 {
					sbr.shardc <- shard
				}
				return
			}
			shard.recs.addRecord(rec)
			if shard.recs.numRecords() >= recordsPerShard {
				sbr.shardc <- shard
				shardIdx++
				shard = &AdjacentBAMShard{
					ShardIdx: shardIdx,
					recs:     newRecBuffer(recordsPerShard),
					errors:   sbr.errors,
				}
			}
		}
	}()
	return sbr, nil
}

// Header returns the SAM Header held by the Reader.
func (r *AdjacentShardedBAMReader) Header() *sam.Header {
	return r.header
}

// GetShard returns one AdjacentBAMShard from AdjacentShardedBAMReader. GetShard
// will wait until a shard is available, or until AdjacentShardedBAMReader
// has no more shards to return.
func (r *AdjacentShardedBAMReader) GetShard() *AdjacentBAMShard {
	s, ok := <-r.shardc
	if !ok {
		return nil
	}
	return s
}

// AdjacentBAMShard represents an ordered subset of records from an
// AdjacentShardedBAMReader. The order of AdjacentBAMShard is
// determined by the shard's ShardIdx. If all of an
// AdjacentShardedBAMReader's AdjacentBAMShards were read sequentially
// from shard 0 to n, all the records would be read in the same order
// they appear in the underlying BAM file.
type AdjacentBAMShard struct {
	// ShardIdx is the index of the shard.
	// Indexing starts at zero.
	ShardIdx int
	// recs is the ordered set of records
	// belong to this shard.
	recs recBuffer
	// pair is the current record pair.
	pair Pair
	// errors is used to halt the reading
	// of the underlying BAM file and the
	// scanning of all shards if the BAM
	// reader or any shard encounters an
	// error.
	errors *errors.Once
}

// Record returns the current pair, or an error.
//
// REQUIRES: Scan() has been called and its last call returned true.
func (s *AdjacentBAMShard) Record() Pair { return s.pair }

// Scan reads the next record. It returns true if a record has been read
// or if an error is encountered, and false on end of data stream.
func (s *AdjacentBAMShard) Scan() bool {
	// A failure in reading one shard must result
	// in the failure of reading all shards.
	if err := s.errors.Err(); err != nil {
		s.pair = Pair{
			Err: fmt.Errorf("error while reading shard %d: stopped reading shard due to read error in other worker: %s", s.ShardIdx, err),
		}
		return true
	}
	defer func() {
		if s.pair.Err != nil {
			s.errors.Set(s.pair.Err)
		}
	}()

	record := s.recs.getRecord()
	// Stop scanning because there are no
	// more records to read.
	if record == nil {
		return false
	}
	mate := s.recs.getRecord()
	if mate == nil {
		s.pair = Pair{
			Err: fmt.Errorf("error while reading shard %d: could not find mate for record %s", s.ShardIdx, record),
		}
		return true
	}

	// Check to see if a pair is valid. A valid pair has the following properties:
	// 1. The Pos/MatePos from one record corresponds to the Pos/MatePos of the other record.
	// 2. One record has an R1 flag and the other record has an R2 flag.
	if record.Pos != mate.MatePos || record.MatePos != mate.Pos {
		s.pair = Pair{
			Err: fmt.Errorf("shard %d: records %s and %s do not have matching pos/mate pos", s.ShardIdx, record, mate),
		}
		return true
	}
	switch {
	case IsRead1(record) && IsRead2(mate):
		s.pair = Pair{R1: record, R2: mate}
	case IsRead1(mate) && IsRead2(record):
		s.pair = Pair{R1: mate, R2: record}
	default:
		s.pair = Pair{
			Err: fmt.Errorf("shard %d: records %s and %s are not a valid R1/R2 pair", s.ShardIdx, record, mate),
		}
	}
	return true
}

// recBuffer keeps track of the
// next record to be read.
type recBuffer struct {
	recs   []*sam.Record
	recIdx int
}

func newRecBuffer(maxRecords int) recBuffer {
	return recBuffer{
		recs: make([]*sam.Record, 0, maxRecords),
	}
}

func (r *recBuffer) numRecords() int {
	return len(r.recs)
}

func (r *recBuffer) addRecord(rec *sam.Record) {
	r.recs = append(r.recs, rec)
}

func (r *recBuffer) getRecord() *sam.Record {
	defer func() { r.recIdx++ }()
	if r.recIdx >= len(r.recs) {
		return nil
	}
	return r.recs[r.recIdx]
}
