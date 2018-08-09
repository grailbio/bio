// Copyright 2018 GRAIL, Inc. All rights reserved.
// Use of this source code is governed by the Apache 2.0
// license that can be found in the LICENSE file.

package fieldio

import (
	"encoding/binary"
	"math"

	"github.com/grailbio/base/log"
)

// byteBuffer is a wrapper around standard varint encoder and decoder to allow
// for automatic buffer sizing. It can be used either for reading or writing,
// but not both at the same time.
type byteBuffer struct {
	n   int
	buf []byte
}

// Reader functions
//
// TODO(saito) Do more strict error checks. Currently if the buffer is corrupt
// the methods return arbitrary values.

// Uint16 reads a fixed16 value.
func (b *byteBuffer) Uint16() uint16 {
	value := binary.LittleEndian.Uint16(b.buf[b.n:])
	b.n += 2
	return value
}

// Float64 reads a float64 value.
func (b *byteBuffer) Float64() float64 {
	value := binary.LittleEndian.Uint64(b.buf[b.n:])
	b.n += 8
	return math.Float64frombits(value)
}

// Uint8 reads a fixed8 value.
func (b *byteBuffer) Uint8() uint8 {
	value := b.buf[b.n]
	b.n++
	return value
}

// Varint64 reads a signed varint.
func (b *byteBuffer) Varint64() int64 {
	value, n := binary.Varint(b.buf[b.n:])
	if n <= 0 {
		log.Panic("byteBuffer.Varint64: underflow")
	}
	b.n += n
	return value
}

// Varint32 reads a signed 32bit varint.
func (b *byteBuffer) Varint32() int {
	value := b.Varint64()
	if value < math.MinInt32 || value > math.MaxInt32 {
		panic(value)
	}
	return int(value)
}

// Uvarint32 reads a unsigned 32bit varint.
func (b *byteBuffer) Uvarint32() uint32 {
	value, n := binary.Uvarint(b.buf[b.n:])
	if n <= 0 {
		log.Panic("byteBuffer.Uvarint32: underflow")
	}
	b.n += n
	if value > math.MaxUint32 {
		panic(value)
	}
	return uint32(value)
}

// Uvarint64 reads a unsigned varint.
func (b *byteBuffer) Uvarint64() uint64 {
	value, n := binary.Uvarint(b.buf[b.n:])
	if n <= 0 {
		log.Panic("byteBuffer.Uvarint64: underflow")
	}
	b.n += n
	return value
}

// RawBytes extracts the next n bytes.
func (b *byteBuffer) RawBytes(n int) []byte {
	value := b.buf[b.n : b.n+n]
	b.n += n
	return value
}

// Writer functions

// Ensure that b.buf can store at least "bytes" more bytes.
func (b *byteBuffer) ensure(bytes int) {
	if cap(b.buf) >= b.n+bytes {
		return
	}
	newCap := ((b.n+bytes)/16 + 1) * 16
	if newCap < cap(b.buf)*2 {
		newCap = cap(b.buf) * 2
	}
	newBuf := make([]byte, newCap)
	copy(newBuf, b.Bytes())
	b.buf = newBuf
}

// PutByte adds one byte to the buffer.
func (b *byteBuffer) PutByte(value uint8) {
	b.ensure(1)
	b.buf[b.n] = value
	b.n++
}

// PutBytes add bytes raw, w/o prefixing its length.
func (b *byteBuffer) PutBytes(data []byte) {
	dataLen := len(data)
	b.ensure(dataLen)
	if dataLen > len(b.buf)-b.n {
		log.Panicf("dataLen: %v, remaining: %v %v", dataLen, len(b.buf)-b.n, b.n)
	}
	copy(b.buf[b.n:], data)
	b.n += dataLen
}

// PutString adds a string, w/o prefixing its length.
func (b *byteBuffer) PutString(data string) {
	dataLen := len(data)
	b.ensure(dataLen)
	if dataLen > len(b.buf)-b.n {
		log.Panicf("dataLen: %v, remaining: %v %v", dataLen, len(b.buf)-b.n, b.n)
	}
	copy(b.buf[b.n:], data)
	b.n += dataLen
}

// PutUint16 adds the value as a fixed16.
func (b *byteBuffer) PutUint16(value uint16) {
	b.ensure(2)
	binary.LittleEndian.PutUint16(b.buf[b.n:], value)
	b.n += 2
}

// PutFloat64 adds the value as a float64.
func (b *byteBuffer) PutFloat64(value float64) {
	b.ensure(8)
	binary.LittleEndian.PutUint64(b.buf[b.n:], math.Float64bits(value))
	b.n += 8
}

// PutVarint adds the value as a signed varint
func (b *byteBuffer) PutVarint(value int64) {
	b.ensure(binary.MaxVarintLen64)
	b.n += binary.PutVarint(b.buf[b.n:], value)
}

// PutVarint adds the value as a signed uvarint
func (b *byteBuffer) PutUvarint(value uint64) {
	b.ensure(binary.MaxVarintLen64)
	b.n += binary.PutUvarint(b.buf[b.n:], value)
}

// Bytes returns the data written so far.
func (b *byteBuffer) Bytes() []byte { return b.buf[:b.n] }

// Len returns len(Bytes()).
func (b *byteBuffer) Len() int { return b.n }

// resizeBuf resizes "*buf" to exactly "size" bytes. The existing data may be
// destroyed.
func resizeBuf(buf *[]byte, size int) {
	if *buf == nil || cap(*buf) < size {
		*buf = make([]byte, size)
	} else {
		*buf = (*buf)[:size]
	}
}
