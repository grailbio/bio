package pam

import (
	"encoding/binary"
	"fmt"
	"io"
	"path/filepath"
	"sync"

	"github.com/biogo/hts/sam"
	"github.com/grailbio/base/errorreporter"
	"github.com/grailbio/base/file"
	"github.com/grailbio/base/recordio"
	"github.com/grailbio/base/recordio/recordioiov"
	"github.com/grailbio/base/recordio/recordiozstd"
	"github.com/grailbio/base/syncqueue"
	grailunsafe "github.com/grailbio/base/unsafe"
	"github.com/grailbio/base/vcontext"
	"github.com/grailbio/bio/biopb"
	gbam "github.com/grailbio/bio/encoding/bam"
	"v.io/x/lib/vlog"
)

const (
	// DefaultMaxBufSize is the default value for WriteOpts.MaxBufSize.
	DefaultMaxBufSize = 8 << 20
	// DefaultWriteParallelism is the default value for WriteOpts.MaxFlushParallelism.
	DefaultWriteParallelism = 4
)

// fieldWriter buffers values of one field and writes them to a recordio file.
type fieldWriter struct {
	field      gbam.FieldType   // const after construction
	label      string           // for logging
	path       string           // The PAM path prefix
	shardRange biopb.CoordRange // The shard record range; passed from Opts.Range
	out        file.File        // output destination.
	wout       io.Writer        // out.Writer
	rio        recordio.Writer  // recordio wrapper for out.

	// Value to be assigned to the "seq" field of a new fieldWriteBuf.
	nextBlockSeq int

	// The current buffer. Once buf becomes full, it is asynchorously compressed
	// and written to "rio".
	buf         *fieldWriteBuf
	bufFreePool *syncqueue.LIFO // Points to Writer.bufFreePool.

	err *errorreporter.T // Any error encountered so far.

	// The following fields are used by async buf flusher.
	mu   *sync.Mutex
	cond *sync.Cond
	// The last block (seq) flushed into the data file.  Guarded by mu.
	lastBlockFlushed int
	// Block indexes generated so far. guarded by mu.
	blockIndexes []biopb.PAMBlockIndexEntry
}

// Serialize "msg" into a single-block recordio file "path".  Existing contents
// of "path" is clobbered.
func writeShardIndex(path string, msg *biopb.PAMShardIndex) error {
	ctx := vcontext.Background()
	data, e := msg.Marshal()
	if e != nil {
		return e
	}
	out, e := file.Create(ctx, path)
	if e != nil {
		return e
	}
	err := errorreporter.T{}
	rio := recordio.NewWriter(out.Writer(ctx), recordio.WriterOpts{
		Transformers: []string{"zstd"},
	})
	rio.Append(data)
	err.Set(rio.Finish())
	err.Set(out.Close(ctx))
	return err.Err()
}

// fieldWriteBuf contains bytes for one field but for sam.Records in a
// compression block.
type fieldWriteBuf struct {
	field      gbam.FieldType
	label      string      // for vlogging only.
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
	return wb.defaultBuf.n + wb.blobBuf.n
}

// Reset the state so that it can be reused for serializing another set of
// records. It avoids unnecessary memory allocation.
func (wb *fieldWriteBuf) reset(seq int, label string, f gbam.FieldType) {
	wb.field = f
	wb.seq = seq
	wb.label = label
	wb.defaultBuf.n = 0
	wb.blobBuf.n = 0
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
			vlog.Fatalf("Record addr decreased from %+v to %+v", wb.endAddr, addr)
		}
	}
	wb.endAddr = addr
	wb.numRecords++
}

func (wb *fieldWriteBuf) putCoordField(addr biopb.Coord, refID int, pos int) {
	wb.updateAddrBounds(addr)
	delta := int64(refID) - wb.prevInt64Value0
	wb.defaultBuf.PutVarint(delta)
	wb.prevInt64Value0 = int64(refID)

	delta = int64(pos) - wb.prevInt64Value1
	wb.blobBuf.PutVarint(delta)
	wb.prevInt64Value1 = int64(pos)
}

func (wb *fieldWriteBuf) putVarintDeltaField(addr biopb.Coord, value int64) {
	wb.updateAddrBounds(addr)
	b := &wb.defaultBuf
	delta := value - wb.prevInt64Value0
	b.PutVarint(delta)
	wb.prevInt64Value0 = value
}

func (wb *fieldWriteBuf) putVarintField(addr biopb.Coord, value int64) {
	wb.updateAddrBounds(addr)
	b := &wb.defaultBuf
	b.PutVarint(value)
}

func (wb *fieldWriteBuf) putUint16Field(addr biopb.Coord, v uint16) {
	wb.updateAddrBounds(addr)
	wb.defaultBuf.PutUint16(v)
}

func (wb *fieldWriteBuf) putByteField(addr biopb.Coord, v byte) {
	wb.updateAddrBounds(addr)
	wb.defaultBuf.PutByte(v)
}

func (wb *fieldWriteBuf) putLengthPrefixedBytes(data []byte) {
	wb.defaultBuf.PutUvarint(uint64(len(data)))
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

func (wb *fieldWriteBuf) putStringDeltaField(addr biopb.Coord, data string) {
	wb.updateAddrBounds(addr)
	prefix, delta := computeDiff(grailunsafe.BytesToString(wb.prevString), data)
	wb.defaultBuf.PutUvarint(uint64(prefix))
	wb.defaultBuf.PutUvarint(uint64(len(delta)))
	wb.blobBuf.PutString(delta)

	resizeBuf(&wb.prevString, len(data))
	copy(wb.prevString, data)
}

func (wb *fieldWriteBuf) putQualField(addr biopb.Coord, qual []byte) {
	wb.updateAddrBounds(addr)
	wb.putLengthPrefixedBytes(qual)
}

func (wb *fieldWriteBuf) putAuxField(addr biopb.Coord, aa []sam.Aux) {
	wb.updateAddrBounds(addr)
	wb.defaultBuf.PutUvarint(uint64(len(aa)))
	for _, a := range aa {
		wb.blobBuf.PutBytes(a[:3])
	}
	for _, a := range aa {
		switch a[2] {
		case 'A', 'c', 'C': // ascii, int8, uint8
			if len(a) != 4 {
				vlog.Fatal(a)
			}
			wb.blobBuf.PutByte(a[3])
		case 's', 'S': // int16, uint16
			if len(a) != 5 {
				vlog.Fatal(a)
			}
			wb.blobBuf.PutBytes(a[3:5])
		case 'i', 'I', 'f': // int32, uint32, float32
			if len(a) != 7 {
				vlog.Fatal(a)
			}
			wb.blobBuf.PutBytes(a[3:7])
		case 'Z', 'H': // text, hexstr
			wb.putLengthPrefixedBytes(a[3:])
		default:
			vlog.Fatal(a)
		}
	}
}

func (wb *fieldWriteBuf) putCigarField(addr biopb.Coord, cigar sam.Cigar) {
	wb.updateAddrBounds(addr)
	wb.defaultBuf.PutUvarint(uint64(len(cigar)))
	for _, op := range cigar {
		wb.defaultBuf.PutUvarint(uint64(op))
	}
}

func (wb *fieldWriteBuf) putSeqField(addr biopb.Coord, seq sam.Seq) {
	wb.updateAddrBounds(addr)
	wb.defaultBuf.PutUvarint(uint64(seq.Length))
	wb.blobBuf.PutBytes(gbam.UnsafeDoubletsToBytes(seq.Seq))
}

func (fw *fieldWriter) flushBuf() {
	fw.rio.Append(fw.buf)
	fw.rio.Flush()
	fw.buf = nil
}

func (fw *fieldWriter) marshalBlock(scratch []byte, v interface{}) ([]byte, error) {
	wb := v.(*fieldWriteBuf)
	defaultData := wb.defaultBuf.Bytes()
	blobData := wb.blobBuf.Bytes()
	defaultOffset := 0
	blobOffset := len(defaultData)

	header := biopb.PAMBlockHeader{
		Field:      int32(wb.field),
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
	bb := byteBuffer{n: 0, buf: tmpBuf0[:]}
	bb.PutVarint(int64(n))
	serialized := recordioiov.Slice(scratch, bb.Len()+n+len(defaultData)+len(blobData))
	copy(serialized, bb.Bytes())
	copy(serialized[bb.Len():], tmpBuf1[:n])
	copy(serialized[bb.Len()+n:], defaultData)
	copy(serialized[bb.Len()+n+len(defaultData):], blobData)
	return serialized, nil
}

func (fw *fieldWriter) indexCallback(loc recordio.ItemLocation, v interface{}) error {
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
	doassert(index.StartAddr.RefId != biopb.InvalidRefID && index.StartAddr.Pos != biopb.InvalidPos)
	doassert(index.EndAddr.RefId != biopb.InvalidRefID && index.EndAddr.Pos != biopb.InvalidPos)
	fw.mu.Lock()
	fw.blockIndexes = append(fw.blockIndexes, index)
	fw.mu.Unlock()
	fw.bufFreePool.Put(wb)
	return nil
}

// WriteOpts defines options for NewWriter.
type WriteOpts struct {
	// MaxBufSize limits the max size of a recordio block, pre compression.
	// If <= 0, DefaultMaxBufSize is used.
	MaxBufSize int

	// WriteParallelism limits the max number of pending recordio flushes
	// allowed. If <= 0, DefaultWriteParallelism is used.
	WriteParallelism int

	// DropFields causes the writer not to write the specified fields to file.
	DropFields []gbam.FieldType

	// Transformers defines the recordio block transformers. It can be used to
	// change the compression algorithm, for example. The value is passed to
	// recordio.WriteOpts.Transformers. If empty, {"zstd"} is used.
	//
	// Currently there is no way to to disable transformation. If you want to minize
	// CPU overheads, pass "zstd 1".
	Transformers []string

	// Range defines the range of records that can be stored in the PAM
	// file.  The range will be encoded in the path name. Also, Write() will
	// cause an error if it sees a record outside the range. An empty range
	// (default) means UniversalRange.
	//
	// The range bound is closed at the start, open at the limit.
	Range biopb.CoordRange
}

// Check that "r" has valid contents, and that its positiion is in range
// [(startRef,startPos), (limitRef, limitPos)).
//
// TODO(saito) The sam writer does more strict checking. Import that.
func validateRecord(r *sam.Record, recRange biopb.CoordRange) error {
	recAddr := gbam.CoordFromSAMRecord(r, 0)
	if recAddr.LT(recRange.Start) {
		return fmt.Errorf("Record (%d,%d) out of start of shard range %+v : record %v",
			r.Ref.ID(), r.Pos, recRange, r)
	}
	if recAddr.GE(recRange.Limit) {
		return fmt.Errorf("Record (%d,%d) out of limit of shard range: %+v : record %v",
			r.Ref.ID(), r.Pos, recRange, r)
	}
	return nil
}

// Validate and fill the option values.
func validateWriteOpts(o *WriteOpts) error {
	if o.MaxBufSize <= 0 {
		o.MaxBufSize = DefaultMaxBufSize
	}
	if o.WriteParallelism <= 0 {
		o.WriteParallelism = DefaultWriteParallelism
	}
	if len(o.Transformers) == 0 {
		o.Transformers = []string{"zstd"}
	}
	return ValidateCoordRange(&o.Range)
}

func newShardIndex(wo WriteOpts, h *sam.Header) biopb.PAMShardIndex {
	index := biopb.PAMShardIndex{}
	index.Magic = ShardIndexMagic
	index.Version = DefaultVersion
	var err error
	if index.EncodedBamHeader, err = gbam.MarshalHeader(h); err != nil {
		// TODO(saito) propagate errors up
		vlog.Fatalf("Encode header: %v, err")
	}
	index.Range = wo.Range
	return index
}

// Writer is a class for generating a PAM rowshard.
type Writer struct {
	label string // For vlogging only.
	opts  WriteOpts
	dir   string // Output destination
	index biopb.PAMShardIndex

	addrGenerator   gbam.CoordGenerator
	bufFreePoolSize int
	bufFreePool     *syncqueue.LIFO
	fieldWriters    [gbam.NumFields]*fieldWriter // Writer for each field

	// Value to be assigned to the "seq" field of a new recBlockWriteBuf.
	nextBlockSeq int
	err          errorreporter.T
}

// Allocate a new buffer and set it in fw.buf. It blocks the caller if there are
// too many flushing already ongoing.
func (fw *fieldWriter) newBuf() {
	if fw.buf != nil {
		vlog.Fatalf("Overwriting buffer %+v", fw)
	}
	vv, ok := fw.bufFreePool.Get()
	if !ok {
		panic("get")
	}
	wb := vv.(*fieldWriteBuf)
	seq := fw.nextBlockSeq
	fw.nextBlockSeq++
	wb.reset(seq, fmt.Sprintf("%s:%d", fw.label, seq), fw.field)
	fw.buf = wb
}

// Close the output file and return any error encountered so far.  No method
// shall be called after fw.close().
//
// REQUIRES: All outstanding flushes have completed.
func (fw *fieldWriter) close() {
	if fw.out != nil {
		fw.rio.Wait()
		index := biopb.PAMFieldIndex{
			Magic:   FieldIndexMagic,
			Version: DefaultVersion,
			Field:   int32(fw.field),
			Blocks:  fw.blockIndexes,
		}
		vlog.VI(1).Infof("creating index with %d blocks", len(index.Blocks))
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

// Write appends a record. It does not flush the record immediately, and the
// record becomes stable only after a successful Close call. "r" can be recycled
// after Write returns. This function is thread compatible.
//
// REQUIRES: records must be added in increasing position order (Cf. RecAddr).
func (w *Writer) Write(r *sam.Record) {
	if w.err.Err() != nil {
		return
	}
	err := validateRecord(r, w.opts.Range)
	if err != nil {
		w.err.Set(err)
		return
	}
	addr := w.addrGenerator.GenerateFromRecord(r)
	if w.fieldWriters[gbam.FieldCoord] != nil {
		w.fieldWriters[gbam.FieldCoord].buf.putCoordField(addr, r.Ref.ID(), r.Pos)
	}
	if w.fieldWriters[gbam.FieldFlags] != nil {
		w.fieldWriters[gbam.FieldFlags].buf.putUint16Field(addr, uint16(r.Flags))
	}
	if w.fieldWriters[gbam.FieldMapq] != nil {
		w.fieldWriters[gbam.FieldMapq].buf.putByteField(addr, r.MapQ)
	}
	if w.fieldWriters[gbam.FieldCigar] != nil {
		w.fieldWriters[gbam.FieldCigar].buf.putCigarField(addr, r.Cigar)
	}
	if w.fieldWriters[gbam.FieldMateRefID] != nil {
		w.fieldWriters[gbam.FieldMateRefID].buf.putVarintDeltaField(addr, int64(r.MateRef.ID()))
	}
	if w.fieldWriters[gbam.FieldMatePos] != nil {
		w.fieldWriters[gbam.FieldMatePos].buf.putVarintDeltaField(addr, int64(r.MatePos))
	}
	if w.fieldWriters[gbam.FieldTempLen] != nil {
		w.fieldWriters[gbam.FieldTempLen].buf.putVarintField(addr, int64(r.TempLen))
	}
	if w.fieldWriters[gbam.FieldName] != nil {
		w.fieldWriters[gbam.FieldName].buf.putStringDeltaField(addr, r.Name)
	}
	if w.fieldWriters[gbam.FieldSeq] != nil {
		w.fieldWriters[gbam.FieldSeq].buf.putSeqField(addr, r.Seq)
	}
	if w.fieldWriters[gbam.FieldQual] != nil {
		w.fieldWriters[gbam.FieldQual].buf.putQualField(addr, r.Qual)
	}
	if w.fieldWriters[gbam.FieldAux] != nil {
		w.fieldWriters[gbam.FieldAux].buf.putAuxField(addr, r.AuxFields)
	}
	for _, fw := range w.fieldWriters {
		if fw != nil && fw.buf.totalLen() >= w.opts.MaxBufSize {
			fw.flushBuf()
			fw.newBuf()
		}
	}
}

// Close must be called exactly once. After close, no operation other than Err()
// may be called.
func (w *Writer) Close() error {
	wg := sync.WaitGroup{}
	for _, fw := range w.fieldWriters {
		if fw == nil {
			continue
		}
		fb := fw.buf
		if fb.numRecords > 0 {
			vlog.VI(1).Infof("%v: Start flush (close)", fb.label)
			fw.flushBuf()
		} else {
			fw.bufFreePool.Put(fb)
			fw.buf = nil
		}
		wg.Add(1)
		go func(fw *fieldWriter) {
			fw.close()
			wg.Done()
		}(fw)
	}
	// Wait for all flushing ops.
	wg.Wait()
	for i := 0; i < w.bufFreePoolSize; i++ {
		if _, ok := w.bufFreePool.Get(); !ok {
			panic("get")
		}
	}
	if w.err.Err() != nil {
		return w.err.Err()
	}
	return writeShardIndex(ShardIndexPath(w.dir, w.opts.Range), &w.index)
}

// NewWriter creates a new PAM writer. Files are created in "dir". If "dir" does
// not exist already it is created. Existing contents of "dir", if any, are
// deleted.
func NewWriter(wo WriteOpts, samHeader *sam.Header, dir string) *Writer {
	w := &Writer{
		opts:          wo,
		dir:           dir,
		nextBlockSeq:  0,
		addrGenerator: gbam.NewCoordGenerator(),
		err:           errorreporter.T{},
	}
	if w.err.Set(validateWriteOpts(&w.opts)); w.err.Err() != nil {
		return w
	}
	dropField := [gbam.NumFields]bool{}
	for _, f := range wo.DropFields {
		dropField[f] = true
	}
	nWrittenFields := 0
	for _, f := range dropField {
		if !f {
			nWrittenFields++
		}
	}

	w.label = fmt.Sprintf("%s:%s", dir, CoordRangePathString(w.opts.Range))
	w.bufFreePoolSize = w.opts.WriteParallelism * nWrittenFields
	w.bufFreePool = syncqueue.NewLIFO()
	for i := 0; i < w.bufFreePoolSize; i++ {
		w.bufFreePool.Put(&fieldWriteBuf{})
	}
	w.index = newShardIndex(w.opts, samHeader)
	for f := range w.fieldWriters {
		if dropField[f] {
			continue
		}
		fw := newFieldWriter(dir, w.opts, gbam.FieldType(f), w.bufFreePool, &w.err)
		w.fieldWriters[f] = fw
		fw.newBuf()
	}
	return w
}

func newFieldWriter(path string, opts WriteOpts, f gbam.FieldType, bufFreePool *syncqueue.LIFO, errReporter *errorreporter.T) *fieldWriter {
	mu := &sync.Mutex{}
	fw := &fieldWriter{
		field:            f,
		label:            fmt.Sprintf("%s:%s:%v", filepath.Base(path), CoordRangePathString(opts.Range), f),
		path:             path,
		shardRange:       opts.Range,
		bufFreePool:      bufFreePool,
		mu:               mu,
		cond:             sync.NewCond(mu),
		lastBlockFlushed: -1,
		err:              errReporter,
	}
	// Create a recordio file
	ctx := vcontext.Background()
	out, err := file.Create(ctx, FieldDataPath(path, opts.Range, f))
	if err != nil {
		fw.err.Set(err)
		return fw
	}
	fw.out = out
	fw.wout = out.Writer(ctx)
	fw.rio = recordio.NewWriter(fw.wout, recordio.WriterOpts{
		Transformers: opts.Transformers,
		Marshal:      fw.marshalBlock,
		Index:        fw.indexCallback,
	})
	fw.rio.AddHeader(recordio.KeyTrailer, true)
	return fw
}

// Err returns any error encountered so far.
func (w *Writer) Err() error {
	return w.err.Err()
}

func init() {
	recordiozstd.Init()
}
