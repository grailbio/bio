package fieldio

import (
	"testing"

	"github.com/grailbio/testutil/expect"
)

func TestByteBuffer(t *testing.T) {
	b := byteBuffer{}
	b.PutUint8(123)
	b.PutUint16(456)
	b.PutFloat64(23465.5)
	b.PutVarint64(1234)
	b.PutUvarint64(333334)
	b.PutString("hello")

	r := b
	expect.EQ(t, r.Uint8(), uint8(123))
	expect.EQ(t, r.Uint16(), uint16(456))
	expect.EQ(t, r.Float64(), 23465.5)
	expect.EQ(t, r.Varint64(), int64(1234))
	expect.EQ(t, r.Uvarint64(), uint64(333334))
	expect.EQ(t, string(r.RawBytes(5)), "hello")
	expect.EQ(t, len(r), 0)
}
