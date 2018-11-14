// Copyright 2018 GRAIL, Inc. All rights reserved.
// Use of this source code is governed by the Apache 2.0
// license that can be found in the LICENSE file.

package fieldio

import (
	"encoding/binary"
	"math"
)

// byteBuffer is a wrapper around standard varint encoder and decoder to allow
// for automatic buffer sizing. It can be used either for reading or writing,
// but not both at the same time.
type byteBuffer []byte

// Reader functions
//
// TODO(saito) Do more strict error checks. Currently if the buffer is corrupt
// the methods return arbitrary values.

// Uint16 reads a fixed16 value.
func (b *byteBuffer) Uint16() uint16 {
	value := binary.LittleEndian.Uint16(*b)
	*b = (*b)[2:]
	return value
}

// Float64 reads a float64 value.
func (b *byteBuffer) Float64() float64 {
	value := binary.LittleEndian.Uint64(*b)
	*b = (*b)[8:]
	return math.Float64frombits(value)
}

// Uint8 reads a fixed8 value.
func (b *byteBuffer) Uint8() uint8 {
	value := (*b)[0]
	*b = (*b)[1:]
	return value
}

// Varint64 reads a signed varint.
func (b *byteBuffer) Varint64() int64 {
	value, n := binary.Varint(*b)
	if n <= 0 {
		panic("byteBuffer.Varint64: underflow")
	}
	*b = (*b)[n:]
	return value
}

// Uvarint32 reads a unsigned 32bit varint.
func (b *byteBuffer) Uvarint32() uint32 {
	value, n := binary.Uvarint(*b)
	if n <= 0 || value > math.MaxUint32 {
		panic("byteBuffer.Uvarint32")
	}
	*b = (*b)[n:]
	return uint32(value)
}

// Uvarint64 reads a unsigned varint.
func (b *byteBuffer) Uvarint64() uint64 {
	value, n := binary.Uvarint(*b)
	if n <= 0 {
		panic("byteBuffer.Uvarint64: underflow")
	}
	*b = (*b)[n:]
	return value
}

// RawBytes extracts the next n bytes.
func (b *byteBuffer) RawBytes(n int) []byte {
	value := (*b)[:n]
	*b = (*b)[n:]
	return value
}

// Writer functions

// Ensure that b.buf can store at least "bytes" more bytes.
func (b *byteBuffer) alloc(bytes int) []byte {
	blen := len(*b)
	newLen := blen + bytes
	if cap(*b) >= newLen {
		(*b) = (*b)[:newLen]
		return (*b)[blen:]
	}
	newCap := (newLen/16 + 1) * 16
	if newCap < cap(*b)*2 {
		newCap = cap(*b) * 2
	}
	newBuf := make([]byte, newLen, newCap)
	copy(newBuf, *b)
	*b = newBuf
	return (*b)[blen:]
}

// PutUint8 adds one byte to the buffer.
func (b *byteBuffer) PutUint8(value uint8) {
	(b.alloc(1))[0] = value
}

// PutBytes add bytes raw, w/o prefixing its length.
func (b *byteBuffer) PutBytes(data []byte) {
	copy(b.alloc(len(data)), data)
}

// PutString adds a string, w/o prefixing its length.
func (b *byteBuffer) PutString(data string) {
	copy(b.alloc(len(data)), data)
}

// PutUint16 adds the value as a fixed16.
func (b *byteBuffer) PutUint16(value uint16) {
	binary.LittleEndian.PutUint16(b.alloc(2), value)
}

// PutFloat64 adds the value as a float64.
func (b *byteBuffer) PutFloat64(value float64) {
	binary.LittleEndian.PutUint64(b.alloc(8), math.Float64bits(value))
}

// PutVarint adds the value as a signed varint
func (b *byteBuffer) PutVarint64(value int64) {
	x := b.alloc(binary.MaxVarintLen64)
	n := binary.PutVarint(x, value)
	if delta := binary.MaxVarintLen64 - n; delta != 0 {
		(*b) = (*b)[:len(*b)-delta]
	}
}

// PutVarint adds the value as a signed uvarint
func (b *byteBuffer) PutUvarint64(value uint64) {
	x := b.alloc(binary.MaxVarintLen64)
	n := binary.PutUvarint(x, value)
	if delta := binary.MaxVarintLen64 - n; delta != 0 {
		(*b) = (*b)[:len(*b)-delta]
	}
}

// resizeBuf resizes "*buf" to exactly "size" bytes. The existing data may be
// destroyed.
func resizeBuf(buf *[]byte, size int) {
	if *buf == nil || cap(*buf) < size {
		*buf = make([]byte, size)
	} else {
		*buf = (*buf)[:size]
	}
}
