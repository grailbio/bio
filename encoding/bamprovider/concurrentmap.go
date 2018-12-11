package bamprovider

import (
	"sync"

	"github.com/blainsmith/seahash"
	"github.com/grailbio/base/unsafe"
	"github.com/grailbio/hts/sam"
)

const numConcurrentMapShards = 1024

type mapShard struct {
	mu    sync.Mutex
	mates map[string]*sam.Record
}

// concurrentMap is a sharded, thread-safe map from sequence name to sam.Record.
type concurrentMap struct {
	shards [numConcurrentMapShards]mapShard
}

func newConcurrentMap() *concurrentMap {
	m := &concurrentMap{}
	for i := 0; i < len(m.shards); i++ {
		m.shards[i].mates = make(map[string]*sam.Record)
	}
	return m
}

// lookupAndDelete is Called when r's mate is in a different shard (say
// X). Returns a non-nil record if the mate was already added to d.mates by the
// thread handling shard X. Else, it puts r in d.mates and returns nil.
func (m *concurrentMap) lookupAndDelete(r *sam.Record) *sam.Record {
	h := seahash.Sum64(unsafe.StringToBytes(r.Name))
	shard := &m.shards[int(h%uint64(numConcurrentMapShards))]

	shard.mu.Lock()
	mate, ok := shard.mates[r.Name]
	if ok {
		// The mate for this record has been seen already, in a different shard.
		delete(shard.mates, r.Name)
	} else {
		// The mate for this record hasn't been seen yet, so store it until then.
		shard.mates[r.Name] = r
		mate = nil
	}
	shard.mu.Unlock()
	return mate
}

// approxSize returns the approximate number of entries in the map.  It returns
// a correct number iff it is invoked when no other thread is accessing the map.
func (m *concurrentMap) approxSize() int {
	n := 0
	for i := range m.shards {
		s := &m.shards[i]
		s.mu.Lock()
		n += len(s.mates)
		s.mu.Unlock()
	}
	return n
}
