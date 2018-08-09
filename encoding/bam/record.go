// Copyright 2018 GRAIL, Inc. All rights reserved.
// Use of this source code is governed by the Apache 2.0
// license that can be found in the LICENSE file.

//go:generate ../../../base/gtl/generate.py --prefix=Record -DELEM=*Record --package=bam --output=record_pool.go ../../../base/gtl/randomized_freepool.go.tpl

package bam

import (
	"sync/atomic"
	"unsafe"

	"github.com/biogo/hts/sam"
	gunsafe "github.com/grailbio/base/unsafe"
	"v.io/x/lib/vlog"
)

// Record extends biogo sam.Record with grail-specific stuff.  The
// bioinformatics code inside grail always uses this type, sometimes casting to
// sam.Record.
type Record struct {
	sam.Record

	// Magic is fixed to bam.Magic to detect when this object is bam.Record as
	// opposed to sam.Record. This check is fundamentally unsafe and production
	// code shouldn't rely on it.
	Magic uint64

	// Scratch is used by the record parser to store internal data structures.
	Scratch []byte
}

// Magic is the value of Record.Magic.
const Magic = uint64(0x93c9838d4d9f4f71)

// ResizeScratch makes *buf exactly n bytes long.
func ResizeScratch(buf *[]byte, n int) {
	if cap(*buf) < n {
		// Allocate slightly more memory than needed to prevent frequent
		// reallocation.
		size := (n/16 + 1) * 16
		*buf = make([]byte, n, size)
	} else {
		gunsafe.ExtendBytes(buf, n)
	}
}

// CastUp casts bam.Record to biogo sam.Record.
func CastUp(rb *Record) *sam.Record {
	return (*sam.Record)(unsafe.Pointer(rb))
}

// CastDown casts sam.Record to bam.Record. It panics if the input's type is not
// bam.Record.
func CastDown(rb *sam.Record) *Record {
	rec := (*Record)(unsafe.Pointer(rb))
	if rec.Magic != Magic {
		panic("CastFromRecord: object must be bam.Record, not sam.Record.")
	}
	return rec
}

var recordPool = NewRecordFreePool(func() *Record { return &Record{Magic: Magic} }, 1<<20)

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
			vlog.Errorf(`putSamRecord: object must be bam.Record, not sam.Record. magic %x.
If you see this warning in non-test code path, you MUST fix the problem`, r.Magic)
		}
		return
	}
	recordPool.Put(r)
}
