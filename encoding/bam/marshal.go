// Copyright 2018 GRAIL, Inc. All rights reserved.
// Use of this source code is governed by the Apache 2.0
// license that can be found in the LICENSE file.

//go:generate ../../../base/gtl/generate_randomized_freepool.py --output=bytes_pool --prefix=bytes -DELEM=[]byte --package=bam

package bam

import (
	"bytes"
	"encoding/binary"
	"errors"

	"github.com/biogo/hts/sam"
)

var (
	errNameAbsentOrTooLong           = errors.New("bam: name absent or too long")
	errSequenceQualityLengthMismatch = errors.New("bam: sequence/quality length mismatch")
	bufPool                          = NewbytesFreePool(func() []byte { return nil }, 1024)
)

// buildAux constructs a single byte slice that represents a slice of sam.Aux.
// *buf should be an empty slice on call, and it is filled with the result on
// return.
func buildAux(aa []sam.Aux, buf *[]byte) {
	for _, a := range aa {
		// TODO: validate each 'a'
		*buf = append(*buf, []byte(a)...)
		switch a.Type() {
		case 'Z', 'H':
			*buf = append(*buf, 0)
		}
	}
}

type binaryWriter struct {
	w   *bytes.Buffer
	buf [4]byte
}

func (w *binaryWriter) writeUint8(v uint8) {
	w.buf[0] = v
	w.w.Write(w.buf[:1])
}

func (w *binaryWriter) writeUint16(v uint16) {
	binary.LittleEndian.PutUint16(w.buf[:2], v)
	w.w.Write(w.buf[:2])
}

func (w *binaryWriter) writeInt32(v int32) {
	binary.LittleEndian.PutUint32(w.buf[:4], uint32(v))
	w.w.Write(w.buf[:4])
}

func (w *binaryWriter) writeUint32(v uint32) {
	binary.LittleEndian.PutUint32(w.buf[:4], v)
	w.w.Write(w.buf[:4])
}

// Marshal serializes the record in BAM format.
func Marshal(r *sam.Record, buf *bytes.Buffer) error {
	if len(r.Name) == 0 || len(r.Name) > 254 {
		return errNameAbsentOrTooLong
	}
	if r.Qual != nil && len(r.Qual) != r.Seq.Length {
		return errSequenceQualityLengthMismatch
	}
	scratch := bufPool.Get()
	ResizeScratch(&scratch, 0)
	buildAux(r.AuxFields, &scratch)
	tags := scratch
	bin := binaryWriter{w: buf}
	recLen := bamFixedBytes +
		len(r.Name) + 1 + // Null terminated.
		len(r.Cigar)<<2 + // CigarOps are 4 bytes.
		len(r.Seq.Seq) +
		len(r.Qual) +
		len(tags)

	// Write record header data.
	bin.writeInt32(int32(recLen))
	bin.writeInt32(int32(r.Ref.ID()))
	bin.writeInt32(int32(r.Pos))
	bin.writeUint8(byte(len(r.Name) + 1))
	bin.writeUint8(r.MapQ)
	bin.writeUint16(uint16(r.Bin())) //r.bin
	bin.writeUint16(uint16(len(r.Cigar)))
	bin.writeUint16(uint16(r.Flags))
	bin.writeInt32(int32(r.Seq.Length))
	bin.writeInt32(int32(r.MateRef.ID()))
	bin.writeInt32(int32(r.MatePos))
	bin.writeInt32(int32(r.TempLen))

	// Write variable length data.
	buf.WriteString(r.Name)
	buf.WriteByte(0)
	for _, o := range r.Cigar {
		bin.writeUint32(uint32(o))
	}
	buf.Write(UnsafeDoubletsToBytes(r.Seq.Seq))
	if r.Qual != nil {
		buf.Write(r.Qual)
	} else {
		for i := 0; i < r.Seq.Length; i++ {
			buf.WriteByte(0xff)
		}
	}
	buf.Write(tags)
	bufPool.Put(scratch)
	return nil
}

// MarshalHeader encodes header in BAM binary format.
func MarshalHeader(header *sam.Header) ([]byte, error) {
	bb := bytes.Buffer{}
	if err := header.EncodeBinary(&bb); err != nil {
		return nil, err
	}
	return bb.Bytes(), nil
}
