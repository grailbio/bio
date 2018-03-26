package bam

// Functions in this file provides unsafe casting from and from sam.Record
// fields to []byte.

import (
	"reflect"
	"unsafe"

	"github.com/biogo/hts/sam"
)

// CigarOpSize is the size of one sam.CigarOp, in bytes.
const CigarOpSize = int(unsafe.Sizeof(sam.CigarOp(0)))

// UnsafeBytesToCigar casts src to sam.Cigar.  "src" must store an array of
// uint32s (sam.CigarOps) in host byte order.
func UnsafeBytesToCigar(src []byte) (cigar sam.Cigar) {
	sh := (*reflect.SliceHeader)(unsafe.Pointer(&src))
	dh := (*reflect.SliceHeader)(unsafe.Pointer(&cigar))
	dh.Data = sh.Data
	dh.Len = sh.Len / CigarOpSize
	dh.Cap = sh.Cap / CigarOpSize
	return cigar
}

// UnsafeCigarToBytes casts a cigar string to []byte.
func UnsafeCigarToBytes(src sam.Cigar) (d []byte) {
	sh := (*reflect.SliceHeader)(unsafe.Pointer(&src))
	dh := (*reflect.SliceHeader)(unsafe.Pointer(&d))
	dh.Data = sh.Data
	dh.Len = sh.Len * CigarOpSize
	dh.Cap = sh.Cap * CigarOpSize
	return d
}

// UnsafeBytesToDoublets casts []byte to []sam.Doublet.
func UnsafeBytesToDoublets(src []byte) (d []sam.Doublet) {
	sh := (*reflect.SliceHeader)(unsafe.Pointer(&src))
	dh := (*reflect.SliceHeader)(unsafe.Pointer(&d))
	*dh = *sh
	return d
}

// UnsafeDoubletsToBytes casts []sam.Doublet to []byte.
func UnsafeDoubletsToBytes(src []sam.Doublet) (d []byte) {
	sh := (*reflect.SliceHeader)(unsafe.Pointer(&src))
	dh := (*reflect.SliceHeader)(unsafe.Pointer(&d))
	*dh = *sh
	return d
}
