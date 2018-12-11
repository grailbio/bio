package bampair

import (
	"bytes"
	"encoding/binary"
	"fmt"
	"io"
	"os"
	"path"
	"sync"
	"time"

	"github.com/golang/snappy"
	"github.com/grailbio/base/log"
	grailbam "github.com/grailbio/bio/encoding/bam"
	"github.com/grailbio/hts/bam"
	"github.com/grailbio/hts/sam"
)

// diskMateShard is an on-disk shard of the distant mates.  It is
// meant to store a portion of the distant mates in a way that allows
// the caller to load the distant mates into memory only when they are
// useful.
//
// The number of diskMateShards is controlled by a command line
// parameter passed through distantMateTable.  A reasonable number is
// 1000 diskMateShards, which results in the largest shards being
// roughly 200MB for a 130GB input file with 5% distant mates.  The
// diskMateShards are not evenly sized because the input file is
// unevenly sharded, and has uneven distribution of distant mates.
// The input bam file's shards map to these disk shards evenly, so an
// input bam has 30,000 shards, then the first 30 input shards map to
// the first diskMateShard, and the next 30 input shards map to the
// second diskMateShard.
//
// add(), openReader() and closeReader() are threadsafe.  The call to
// closeWriter() should occur after all calls to add() are finished.
// Calls to getMate() should occur only between a call to openReader()
// and a call to closeReader().
type diskMateShard struct {
	header       *sam.Header
	f            *os.File
	writer       io.WriteCloser
	m            map[string][]*indexedRecord
	mateShardIdx int
	mutex        sync.Mutex
	refcount     int
	buf          bytes.Buffer
}

// newDiskMateShard returns a new diskMateShard.
func newDiskMateShard(header *sam.Header, tempDir string, mateShardIdx, numMateShards int) (*diskMateShard, error) {
	filename := path.Join(tempDir, fmt.Sprintf("mates_%04d_of_%04d", mateShardIdx, numMateShards))
	f, err := os.Create(filename)
	if err != nil {
		return nil, fmt.Errorf("error creating file %s: %v", filename, err)
	}

	return &diskMateShard{
		header:       header,
		f:            f,
		writer:       snappy.NewBufferedWriter(f),
		mateShardIdx: mateShardIdx,
	}, nil
}

func (s *diskMateShard) add(mate *sam.Record, fileIdx uint64) error {
	s.mutex.Lock()
	defer s.mutex.Unlock()

	if err := binary.Write(s.writer, binary.LittleEndian, fileIdx); err != nil {
		return fmt.Errorf("error writing fileIdx to mate shard: %v", err)
	}
	s.buf.Reset()
	if err := bam.Marshal(mate, &s.buf); err != nil {
		return fmt.Errorf("error marshalling to mate shard: %v", err)
	}
	if _, err := s.writer.Write(s.buf.Bytes()); err != nil {
		return fmt.Errorf("error writing record to mate shard: %v", err)
	}
	return nil
}

func (s *diskMateShard) closeWriter() error {
	s.mutex.Lock()
	defer s.mutex.Unlock()

	if err := s.writer.Close(); err != nil {
		return fmt.Errorf("failed to close snappy writer for distantMateShard %s: %v", s.f.Name(), err)
	}
	if err := s.f.Close(); err != nil {
		return fmt.Errorf("failed to close f for distantMateShard %s: %v", s.f.Name(), err)
	}
	return nil
}

func (s *diskMateShard) openReader() error {
	s.mutex.Lock()
	defer s.mutex.Unlock()

	s.refcount++
	if s.refcount > 1 {
		return nil
	}

	f, err := os.Open(s.f.Name())
	if err != nil {
		return fmt.Errorf("error opening file %s: %v", s.f.Name(), err)
	}
	reader := snappy.NewReader(f)

	t0 := time.Now()
	s.m = map[string][]*indexedRecord{}
	records := 0
	log.Debug.Printf("reading disk shard %d", s.mateShardIdx)

	buf := make([]byte, 1024)
	for {
		// Read fileIdx and record size.
		n, err := io.ReadFull(reader, buf[0:12])
		if n > 0 && n < 12 {
			return fmt.Errorf("Not enough bytes: %d, file %s", n, s.f.Name())
		}
		if err == io.EOF {
			break
		} else if err != nil {
			return fmt.Errorf("Error reading mate shard: %v", err)
		}

		fileIdx := binary.LittleEndian.Uint64(buf[0:8])
		size := binary.LittleEndian.Uint32(buf[8:12])

		// Grow the slice if necessary
		for uint32(cap(buf)) < size {
			buf = append(buf[:cap(buf)], 0)
		}
		buf = buf[:size]

		n, err = io.ReadFull(reader, buf)
		if uint32(n) < size {
			return fmt.Errorf("Not enough bytes: %d expected %d file %s", n, size, s.f.Name())
		}
		if err != nil && err != io.EOF {
			return fmt.Errorf("Error reading mate shard: %v", err)
		}

		// Unmarshal record
		record, err := grailbam.Unmarshal(buf, s.header)
		if err != nil {
			return fmt.Errorf("Failed to unmarshal disk mate: %v", err)
		}

		records++
		iRecord := &indexedRecord{record, fileIdx}

		e, found := s.m[iRecord.r.Name]
		if found {
			// We might receive duplicate records here because
			// multiple inputs shards map to one mateShardIdx.  Add
			// only unique records.
			duplicate := false
			for _, existing := range e {
				if recordsEqual(existing.r, iRecord.r) {
					duplicate = true
					break
				}
			}
			if duplicate {
				continue
			}
			// There's no particular order for left and right in a distantMateEntry.
			s.m[iRecord.r.Name] = append(s.m[iRecord.r.Name], iRecord)
			if len(s.m[iRecord.r.Name]) > 2 {
				return fmt.Errorf("Got 3rd read, %v %d", iRecord.r, iRecord.fileIdx)
			}
		} else {
			s.m[iRecord.r.Name] = []*indexedRecord{iRecord}
		}
	}

	if err := f.Close(); err != nil {
		return fmt.Errorf("Error closing %s: %v", s.f.Name(), err)
	}
	t1 := time.Now()
	log.Debug.Printf("read disk shard %d with %d records in %v", s.mateShardIdx, records, t1.Sub(t0))
	return nil
}

func (s *diskMateShard) getMate(shardInfo *ShardInfo, r *sam.Record) (*sam.Record, uint64) {
	e := s.m[r.Name]
	for _, iRecord := range e {
		if !recordsEqual(iRecord.r, r) {
			return iRecord.r, iRecord.fileIdx + shardInfo.getShardStartFileIdx(iRecord.r)
		}
	}
	return nil, 0
}

func (s *diskMateShard) closeReader() {
	s.mutex.Lock()
	defer s.mutex.Unlock()
	s.refcount--
	// Zero out the map to free the memory.
	if s.refcount == 0 {
		s.m = nil
	}
}

func recordsEqual(a, b *sam.Record) bool {
	return a.Name == b.Name &&
		a.Ref.ID() == b.Ref.ID() &&
		a.Pos == b.Pos &&
		a.MapQ == b.MapQ &&
		a.Flags == b.Flags &&
		a.MateRef.ID() == b.MateRef.ID() &&
		a.MatePos == b.MatePos &&
		a.TempLen == b.TempLen
}
