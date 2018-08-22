package fieldio

import (
	"encoding/binary"
	"fmt"
	"io"
	"sync"

	"github.com/biogo/hts/sam"
	"github.com/grailbio/base/errorreporter"
	"github.com/grailbio/base/file"
	"github.com/grailbio/base/log"
	"github.com/grailbio/base/recordio"
	"github.com/grailbio/base/recordio/recordioiov"
	"github.com/grailbio/base/syncqueue"
	grailunsafe "github.com/grailbio/base/unsafe"
	"github.com/grailbio/base/vcontext"
	"github.com/grailbio/bio/biopb"
	gbam "github.com/grailbio/bio/encoding/bam"
	"github.com/grailbio/bio/encoding/pam/pamutil"
)

const (
	// FieldIndexMagic is the value of FieldIndex.Magic.
	FieldIndexMagic = uint64(0xe360ac9026052aca)
)

// Writer buffers values of one field and writes them to a recordio file.
type Writer struct {
	label string          // for logging
	out   file.File       // output destination.
	wout  io.Writer       // out.Writer
	rio   recordio.Writer // recordio wrapper for out.

	// Value to be assigned to the "seq" field of a new fieldWriteBuf.
	nextBlockSeq int

	// The current buffer. Once buf becomes full, it is asynchorously compressed
	// and written to "rio".
	buf         *fieldWriteBuf
	bufFreePool *WriteBufPool // Points to pam.Writer.bufPool.

	err *errorreporter.T // Any error encountered so far.

	// The following fields are used by async buf flusher.
	mu   *sync.Mutex
	cond *sync.Cond
	// The last block (seq) flushed into the data file.  Guarded by mu.
	lastBlockFlushed int
	// Block indexes generated so far. guarded by mu.
	blockIndexes []biopb.PAMBlockIndexEntry
}

// fieldWriteBuf contains bytes for one field but for sam.Records in a
// compression block.
type fieldWriteBuf struct {
	label      string      // for logging only.
	seq        int         // sequence number, 0, 1, 2, ...
	numRecords int         // # of records written so far.
	startAddr  biopb.Coord // addr of the first record stored in this buf.
	endAddr    biopb.Coord // addr of the last record stored in this buf.

	defaultBuf byteBuffer // for storing numeric values
	blobBuf    byteBuffer // for storing string and bytes.

	// Used to prefix-delta-encode a string.  prevString deep-copies the
	// string value, because strings in *sam.Record may be recycled by a
	// freepool.
	prevString []byte
	// For encoding deltas for fields RefID, Pos and MatePos.
	prevInt64Value0 int64
	prevInt64Value1 int64
}

// totalLen computes the total # of bytes stored in the buffer.
func (wb *fieldWriteBuf) totalLen() int {
	return len(wb.defaultBuf) + len(wb.blobBuf)
}

// Reset the state so that it can be reused for serializing another set of
// records. It avoids unnecessary memory allocation.
func (wb *fieldWriteBuf) reset(seq int, label string) {
	wb.seq = seq
	wb.label = label
	wb.defaultBuf = wb.defaultBuf[:0]
	wb.blobBuf = wb.blobBuf[:0]
	wb.prevString = wb.prevString[:0]
	wb.prevInt64Value0 = 0
	wb.prevInt64Value1 = 0
	wb.startAddr = biopb.Coord{biopb.InvalidRefID, biopb.InvalidPos, 0}
	wb.endAddr = biopb.Coord{biopb.InvalidRefID, biopb.InvalidPos, 0}
	wb.numRecords = 0
}

func (wb *fieldWriteBuf) updateAddrBounds(addr biopb.Coord) {
	if wb.numRecords == 0 {
		wb.startAddr = addr
	} else {
		if !wb.endAddr.LT(addr) {
			log.Panicf("Record addr decreased from %+v to %+v", wb.endAddr, addr)
		}
	}
	wb.endAddr = addr
	wb.numRecords++
}

// PutCoordField adds a coordinate field to the buffer.
//
// TODO(saito) we don't need (refid, pos). They can be derived from coord.
func (fw *Writer) PutCoordField(addr biopb.Coord, refID int, pos int) {
	wb := fw.buf
	wb.updateAddrBounds(addr)
	delta := int64(refID) - wb.prevInt64Value0
	wb.defaultBuf.PutVarint64(delta)
	wb.prevInt64Value0 = int64(refID)

	delta = int64(pos) - wb.prevInt64Value1
	wb.blobBuf.PutVarint64(delta)
	wb.prevInt64Value1 = int64(pos)
}

// PutVarintDeltaField adds a varint-delta-encoded field.
func (fw *Writer) PutVarintDeltaField(addr biopb.Coord, value int64) {
	wb := fw.buf
	wb.updateAddrBounds(addr)
	b := &wb.defaultBuf
	delta := value - wb.prevInt64Value0
	b.PutVarint64(delta)
	wb.prevInt64Value0 = value
}

// PutVarintField adds a varint-encoded field.
func (fw *Writer) PutVarintField(addr biopb.Coord, value int64) {
	wb := fw.buf
	wb.updateAddrBounds(addr)
	b := &wb.defaultBuf
	b.PutVarint64(value)
}

// PutUint16Field adds a uint16 field.
func (fw *Writer) PutUint16Field(addr biopb.Coord, v uint16) {
	wb := fw.buf
	wb.updateAddrBounds(addr)
	wb.defaultBuf.PutUint16(v)
}

// PutUint8Field adds a byte field.
func (fw *Writer) PutUint8Field(addr biopb.Coord, v byte) {
	wb := fw.buf
	wb.updateAddrBounds(addr)
	wb.defaultBuf.PutUint8(v)
}

// PutFloat64Field adds a float64 field.
func (fw *Writer) PutFloat64Field(addr biopb.Coord, v float64) {
	wb := fw.buf
	wb.updateAddrBounds(addr)
	wb.defaultBuf.PutFloat64(v)
}

func (wb *fieldWriteBuf) putLengthPrefixedBytes(data []byte) {
	wb.defaultBuf.PutUvarint64(uint64(len(data)))
	wb.blobBuf.PutBytes(data)
}

func computeDiff(prev, cur string) (int, string) {
	minLen := len(prev)
	if len(cur) < minLen {
		minLen = len(cur)
	}
	var i int
	for i = 0; i < minLen; i++ {
		if prev[i] != cur[i] {
			break
		}
	}
	if i >= len(cur) {
		return i, ""
	}
	return i, cur[i:]
}

// BufLen returns the bytelength of the buffer.
func (fw *Writer) BufLen() int {
	return fw.buf.totalLen()
}

// PutStringDeltaField adds a string-delta-encoded field.
func (fw *Writer) PutStringDeltaField(addr biopb.Coord, data string) {
	wb := fw.buf
	wb.updateAddrBounds(addr)
	prefix, delta := computeDiff(grailunsafe.BytesToString(wb.prevString), data)
	wb.defaultBuf.PutUvarint64(uint64(prefix))
	wb.defaultBuf.PutUvarint64(uint64(len(delta)))
	wb.blobBuf.PutString(delta)

	resizeBuf(&wb.prevString, len(data))
	copy(wb.prevString, data)
}

// PutBytesField adds a field consisting of variable-length byte slice
func (fw *Writer) PutBytesField(addr biopb.Coord, data []byte) {
	wb := fw.buf
	wb.updateAddrBounds(addr)
	wb.putLengthPrefixedBytes(data)
}

// PutVarint32sField adds a variable-length int32 slice.
func (fw *Writer) PutVarint32sField(addr biopb.Coord, data []int32) {
	wb := fw.buf
	wb.updateAddrBounds(addr)
	wb.defaultBuf.PutUvarint64(uint64(len(data)))
	for _, v := range data {
		wb.defaultBuf.PutVarint64(int64(v))
	}
}

// PutAuxField adds the aux field.
func (fw *Writer) PutAuxField(addr biopb.Coord, aa []sam.Aux) {
	wb := fw.buf
	wb.updateAddrBounds(addr)
	wb.defaultBuf.PutUvarint64(uint64(len(aa)))
	for _, a := range aa {
		wb.blobBuf.PutBytes(a[:3])
	}
	for _, a := range aa {
		switch a[2] {
		case 'A', 'c', 'C': // ascii, int8, uint8
			if len(a) != 4 {
				log.Panic(a)
			}
			wb.blobBuf.PutUint8(a[3])
		case 's', 'S': // int16, uint16
			if len(a) != 5 {
				log.Panic(a)
			}
			wb.blobBuf.PutBytes(a[3:5])
		case 'i', 'I', 'f': // int32, uint32, float32
			if len(a) != 7 {
				log.Panic(a)
			}
			wb.blobBuf.PutBytes(a[3:7])
		case 'Z', 'H': // text, hexstr
			wb.putLengthPrefixedBytes(a[3:])
		default:
			log.Panic(a)
		}
	}
}

// PutCigarField adds the cigar field.
func (fw *Writer) PutCigarField(addr biopb.Coord, cigar sam.Cigar) {
	wb := fw.buf
	wb.updateAddrBounds(addr)
	wb.defaultBuf.PutUvarint64(uint64(len(cigar)))
	for _, op := range cigar {
		wb.defaultBuf.PutUvarint64(uint64(op))
	}
}

// PutSeqField adds the seq field.
func (fw *Writer) PutSeqField(addr biopb.Coord, seq sam.Seq) {
	wb := fw.buf
	wb.updateAddrBounds(addr)
	wb.defaultBuf.PutUvarint64(uint64(seq.Length))
	wb.blobBuf.PutBytes(gbam.UnsafeDoubletsToBytes(seq.Seq))
}

// FlushBuf starts flushing the buffer to the underlying recordio file. It
// returns before the data is flushed to the storage.
func (fw *Writer) FlushBuf() {
	fw.rio.Append(fw.buf)
	fw.rio.Flush()
	fw.buf = nil
}

func (fw *Writer) marshalBlock(scratch []byte, v interface{}) ([]byte, error) {
	wb := v.(*fieldWriteBuf)
	defaultData := wb.defaultBuf
	blobData := wb.blobBuf
	defaultOffset := 0
	blobOffset := len(defaultData)

	header := biopb.PAMBlockHeader{
		Offset:     uint32(defaultOffset),
		BlobOffset: uint32(blobOffset),
	}
	// tmpBuf0 stores the byte length of the header.
	var tmpBuf0 [binary.MaxVarintLen64]byte
	// tmpBuf1 stores the serialized BlockHeader.
	var tmpBuf1 [12 + gbam.NumFields*(binary.MaxVarintLen64*2+3)]byte

	n, err := header.MarshalTo(tmpBuf1[:])
	if err != nil {
		panic(err)
	}
	bb := byteBuffer(tmpBuf0[:0])
	bb.PutVarint64(int64(n))
	serialized := recordioiov.Slice(scratch, len(bb)+n+len(defaultData)+len(blobData))
	copy(serialized, bb)
	copy(serialized[len(bb):], tmpBuf1[:n])
	copy(serialized[len(bb)+n:], defaultData)
	copy(serialized[len(bb)+n+len(defaultData):], blobData)
	return serialized, nil
}

func (fw *Writer) indexCallback(loc recordio.ItemLocation, v interface{}) error {
	wb := v.(*fieldWriteBuf)
	if loc.Item != 0 {
		panic(loc)
	}
	index := biopb.PAMBlockIndexEntry{
		NumRecords: uint32(wb.numRecords),
		StartAddr:  wb.startAddr,
		EndAddr:    wb.endAddr,
		FileOffset: loc.Block,
	}
	if index.StartAddr.RefId == biopb.InvalidRefID || index.StartAddr.Pos == biopb.InvalidPos ||
		index.EndAddr.RefId == biopb.InvalidRefID || index.EndAddr.Pos == biopb.InvalidPos {
		log.Panic(index)
	}
	fw.mu.Lock()
	fw.blockIndexes = append(fw.blockIndexes, index)
	fw.mu.Unlock()
	fw.bufFreePool.pool.Put(wb)
	return nil
}

// NewBuf allocates a new buffer and set it in fw.buf. It blocks the caller if there are
// too many flushing already ongoing.
func (fw *Writer) NewBuf() {
	if fw.buf != nil {
		log.Panicf("Overwriting buffer %+v", fw)
	}
	vv, ok := fw.bufFreePool.pool.Get()
	if !ok {
		panic("get")
	}
	wb := vv.(*fieldWriteBuf)
	seq := fw.nextBlockSeq
	fw.nextBlockSeq++
	wb.reset(seq, fmt.Sprintf("%s:%d", fw.label, seq))
	fw.buf = wb
}

// Close the output file and return any error encountered so far.  No method
// shall be called after fw.close().
//
// REQUIRES: All outstanding flushes have completed.
func (fw *Writer) Close() {
	fb := fw.buf
	if fb.numRecords > 0 {
		log.Debug.Printf("%v: Start flush (close)", fb.label)
		fw.FlushBuf()
	} else {
		fw.bufFreePool.pool.Put(fb)
		fw.buf = nil
	}
	if fw.out != nil {
		fw.rio.Wait()
		index := biopb.PAMFieldIndex{
			Magic:   FieldIndexMagic,
			Version: pamutil.DefaultVersion,
			Blocks:  fw.blockIndexes,
		}
		log.Debug.Printf("creating index with %d blocks", len(index.Blocks))
		data, err := index.Marshal()
		if err != nil {
			panic(err)
		}
		fw.rio.SetTrailer(data)
		if err := fw.rio.Finish(); err != nil {
			fw.err.Set(err)
		}
		if err := fw.out.Close(vcontext.Background()); err != nil {
			fw.err.Set(err)
		}
	}
}

// NewWriter creates a new field writer that writes to the given path. Label is
// used for logging. Transformers is set as the recordio transformers.
func NewWriter(path, label string, transformers []string, bufFreePool *WriteBufPool, errReporter *errorreporter.T) *Writer {
	mu := &sync.Mutex{}
	fw := &Writer{
		label:            label,
		bufFreePool:      bufFreePool,
		mu:               mu,
		cond:             sync.NewCond(mu),
		lastBlockFlushed: -1,
		err:              errReporter,
	}
	fw.NewBuf()
	// Create a recordio file
	ctx := vcontext.Background()
	out, err := file.Create(ctx, path)
	if err != nil {
		fw.err.Set(err)
		return fw
	}
	fw.out = out
	fw.wout = out.Writer(ctx)
	fw.rio = recordio.NewWriter(fw.wout, recordio.WriterOpts{
		Transformers:        transformers,
		Marshal:             fw.marshalBlock,
		Index:               fw.indexCallback,
		MaxFlushParallelism: 2,
	})
	fw.rio.AddHeader(recordio.KeyTrailer, true)
	return fw
}

type WriteBufPool struct {
	capacity int
	pool     *syncqueue.LIFO
}

func NewBufPool(capacity int) *WriteBufPool {
	p := &WriteBufPool{
		capacity: capacity,
		pool:     syncqueue.NewLIFO(),
	}
	for i := 0; i < capacity; i++ {
		p.pool.Put(&fieldWriteBuf{})
	}
	return p
}

func (p *WriteBufPool) Finish() {
	for i := 0; i < p.capacity; i++ {
		if _, ok := p.pool.Get(); !ok {
			panic("get")
		}
	}
}
