package bam

import (
	"unsafe"

	"github.com/biogo/hts/sam"
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
		*buf = (*buf)[:n]
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
