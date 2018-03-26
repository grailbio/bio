package bam

import (
	"runtime"
	"sync"
	"sync/atomic"
	_ "unsafe" // needed to enable go:linkname

	"github.com/biogo/hts/sam"
	"log"
)

// FreePool is a variation of sync.Pool, specialized for
// sam.Records. Differences are:
//
// - Put() performs power-of-two loadbalancing, and Get() looks only at the
//   local queue.  This improves the performance of Get() on many-core machines,
//   at the cost of slightly more allocations.
//
// - It assumes that GOMAXPROCS is fixed at boot.
//
// - It never frees objects accumulated in the pool. We could add this feature
//   if needed.
type FreePool struct {
	local        []poolLocal
	maxLocalSize int64
}

const maxPrivateElems = 4

type poolLocal struct {
	private     [maxPrivateElems]*Record // Can be used only by the respective P.
	privateSize int

	shared     []*Record  // Can be used by any P.
	sharedSize int64      // ==len(shared), but can be accessed w/o holding mu.
	mu         sync.Mutex // Protects shared.
	pad        [120]byte  // Prevents false sharing.
}

// NewFreePool creates a new free object pool. maxSize bounds the approx max
// number of objects that can be stored in the pool. Beyond this limit, Put()
// call will drop the objects.
func NewFreePool(maxSize int) *FreePool {
	maxProcs := runtime.GOMAXPROCS(0)
	maxLocalSize := -1
	if maxSize > 0 {
		maxLocalSize = maxSize / maxProcs
		if maxLocalSize <= 0 {
			maxLocalSize = 1
		}
	}
	p := &FreePool{
		local:        make([]poolLocal, maxProcs),
		maxLocalSize: int64(maxLocalSize),
	}
	return p
}

func pinPool(p *FreePool) *poolLocal {
	pid := runtime_procPin()
	if int(pid) >= len(p.local) {
		panic(pid)
	}
	return &p.local[pid]
}

// Put adds an object to the freepool. The caller shall not touch the object after the call.
func (p *FreePool) Put(x *Record) {
	l := pinPool(p)
	if l.privateSize < maxPrivateElems {
		l.private[l.privateSize] = x
		l.privateSize++
		x = nil
	}
	runtime_procUnpin()
	if x != nil {
		// Pick another random queue, then add x to the shorter one.
		// This policy ("power of two") reduces load imbalance across
		// queues to log(log(#queues)) .
		//
		// https://www.eecs.harvard.edu/~michaelm/postscripts/mythesis.pdf
		l2 := &p.local[int(fastrand())%len(p.local)]
		lSize := atomic.LoadInt64(&l.sharedSize)
		l2Size := atomic.LoadInt64(&l2.sharedSize)
		if l2Size < lSize {
			l = l2
		}
		l.mu.Lock()
		if p.maxLocalSize >= 0 && l.sharedSize < p.maxLocalSize {
			l.shared = append(l.shared, x)
			atomic.StoreInt64(&l.sharedSize, l.sharedSize+1) // Release store.
		}
		l.mu.Unlock()
	}
}

// Get removes an object from the freepool. If pool is empty, it calls the
// callback passed to NewFreePool.
func (p *FreePool) Get() *Record {
	l := pinPool(p)
	var x *Record
	if l.privateSize > 0 {
		l.privateSize--
		x = l.private[l.privateSize]
		l.private[l.privateSize] = nil
		if x == nil {
			panic(x)
		}
	}
	runtime_procUnpin()
	if x != nil {
		return x
	}
	l.mu.Lock()
	last := len(l.shared) - 1
	if last >= 0 {
		x = l.shared[last]
		l.shared = l.shared[:last]
		atomic.StoreInt64(&l.sharedSize, l.sharedSize-1)
	}
	l.mu.Unlock()
	if x == nil {
		x = &Record{Magic: Magic}
	}
	return x
}

func (p *FreePool) testLen() int {
	n := 0
	for i := range p.local {
		n += p.local[i].privateSize
		n += int(p.local[i].sharedSize)
	}
	return n
}

var recordPool = NewFreePool(1 << 20)

// GetFromFreePool gets a sam.Record object from the singleton freepool, or
// allocate one anew if the pool is empty.
func GetFromFreePool() *Record {
	rec := recordPool.Get()
	rec.Name = ""
	rec.Ref = nil
	rec.MateRef = nil
	rec.Cigar = nil
	rec.Seq = sam.Seq{}
	rec.Qual = nil
	rec.AuxFields = nil
	return rec
}

var nPoolWarnings int32

// PutInFreePool adds "r" to the singleton freepool.  The caller must guarantee
// that there is no outstanding references to "r"; "r" will be overwritten in a
// future.
func PutInFreePool(r *Record) {
	if r == nil {
		panic("r=nil")
	}
	if r.Magic != Magic {
		if atomic.AddInt32(&nPoolWarnings, 1) < 2 {
			log.Printf(`putSamRecord: object must be bam.Record, not sam.Record. magic %x.
If you see this warning in non-test code path, you MUST fix the problem`, r.Magic)
		}
		return
	}
	recordPool.Put(r)
}

// The following functions are defined in go runtime.  To use them, we need to
// import "unsafe", and elsewhere in this package, import "C" to force compiler
// to recognize the "go:linktime" directive. Some of the details are explained
// in the below blog post.
//
// procPin() pins the caller to the current processor, and returns the processor
// id in range [0,GOMAXPROCS). procUnpin() undos the effect of procPin().
//
// http://www.alangpierce.com/blog/2016/03/17/adventures-in-go-accessing-unexported-functions/

//go:linkname runtime_procPin sync.runtime_procPin
//go:nosplit
func runtime_procPin() int

//go:linkname runtime_procUnpin sync.runtime_procUnpin
//go:nosplit
func runtime_procUnpin()

//go:linkname fastrand sync.fastrand
func fastrand() uint32
