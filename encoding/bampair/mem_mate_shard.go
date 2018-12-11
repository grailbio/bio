package bampair

import (
	"fmt"
	"sync"

	"github.com/grailbio/hts/sam"
)

// memMateShard is an in-memory shard of the distant mates.
//
// add(), openReader() and closeReader() are threadsafe.  The call to
// closeWriter() should occur after all calls to add() are finished.
// Calls to getMate() should occur only between a call to openReader()
// and a call to closeReader().
type memMateShard struct {
	m     map[string][]*indexedRecord
	mutex sync.Mutex
}

func newMemMateShard() *memMateShard {
	return &memMateShard{
		m: map[string][]*indexedRecord{},
	}
}

func (s *memMateShard) add(mate *sam.Record, fileIdx uint64) error {
	s.mutex.Lock()
	defer s.mutex.Unlock()

	e, found := s.m[mate.Name]
	if found {
		for _, iRecord := range e {
			if recordsEqual(iRecord.r, mate) {
				return nil
			}
		}
		s.m[mate.Name] = append(s.m[mate.Name], &indexedRecord{mate, fileIdx})
		if len(e) > 2 {
			return fmt.Errorf("Got 3rd read, right is not nil! %v %d", mate, fileIdx)
		}
	} else {
		s.m[mate.Name] = []*indexedRecord{&indexedRecord{mate, fileIdx}}
	}
	return nil
}

func (s *memMateShard) closeWriter() error {
	return nil
}

func (s *memMateShard) openReader() error {
	return nil
}

func (s *memMateShard) getMate(shardInfo *ShardInfo, r *sam.Record) (*sam.Record, uint64) {
	e := s.m[r.Name]
	for _, iRecord := range e {
		if !recordsEqual(iRecord.r, r) {
			return iRecord.r, iRecord.fileIdx + shardInfo.getShardStartFileIdx(iRecord.r)
		}
	}
	return nil, 0
}

func (s *memMateShard) closeReader() {
}
