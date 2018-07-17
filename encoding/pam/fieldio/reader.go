package fieldio

import (
	"encoding/binary"
	"fmt"
	"io"
	"path/filepath"
	"reflect"
	"sync"
	"sync/atomic"
	"unsafe"

	"github.com/biogo/hts/sam"
	"github.com/grailbio/base/errorreporter"
	"github.com/grailbio/base/errors"
	"github.com/grailbio/base/file"
	"github.com/grailbio/base/log"
	"github.com/grailbio/base/recordio"
	gunsafe "github.com/grailbio/base/unsafe"
	"github.com/grailbio/base/vcontext"
	"github.com/grailbio/bio/biopb"
	gbam "github.com/grailbio/bio/encoding/bam"
	"github.com/grailbio/bio/encoding/pam/pamutil"
)

// Reader reads a sequence of values for one field type.
type Reader struct {
	field          gbam.FieldType
	requestedRange biopb.CoordRange           // copy of ReadOpts.Range.
	label          string                     // for vlogging only.
	in             file.File                  // For reading the data file.
	rin            io.ReadSeeker              // in.Reader
	rio            recordio.Scanner           // recordio wrapper for "in".
	index          biopb.PAMFieldIndex        // Contents of *<fieldname>.index file.
	blocks         []biopb.PAMBlockIndexEntry // Subset of index.Blocks that intersect requestedRange.
	fb             fieldReadBuf               // Current buffer being parsed.
	err            *errorreporter.T

	// addrGenerator is used only when field==bam.FieldCoord
	addrGenerator gbam.CoordGenerator
}

// NewReader creates a new Reader. (path,shardRange) defines the path prefix and
// the rowshard range of the data/index files. requestedRange is the coordinate
// range to read, as requested by the user.
func NewReader(path string, shardRange biopb.CoordRange, requestedRange biopb.CoordRange, f gbam.FieldType, errp *errorreporter.T) (*Reader, error) {
	fr := &Reader{
		field: f,
		label: fmt.Sprintf("%s:s%s:u%s(%v)",
			filepath.Base(path), pamutil.CoordRangePathString(shardRange), pamutil.CoordRangePathString(requestedRange), f),
		requestedRange: requestedRange,
		err:            errp,
	}
	ctx := vcontext.Background()
	in, err := file.Open(ctx, pamutil.FieldDataPath(path, shardRange, f))
	if err != nil {
		return fr, errors.E(err, fmt.Sprintf("%v: Failed to open field data file for %v", path, f))
	}
	fr.in = in
	fr.rin = fr.in.Reader(ctx)
	fr.rio = recordio.NewScanner(fr.rin, recordio.ScannerOpts{})
	fr.addrGenerator = gbam.NewCoordGenerator()
	trailer := fr.rio.Trailer()
	if len(trailer) == 0 {
		return fr, fmt.Errorf("%v: file does not contain an index: %v", path, fr.rio.Err())
	}
	if err := fr.index.Unmarshal(trailer); err != nil {
		return fr, errors.E(err, fmt.Sprintf("%v: Failed to unmarshal field index for %v", path, f))
	}
	return fr, nil
}

// For parsing values of one field in one recordio block.
type fieldReadBuf struct {
	header              biopb.PAMBlockHeader
	index               biopb.PAMBlockIndexEntry
	buf                 []byte // raw uncompressed bytes
	defaultBuf, blobBuf byteBuffer
	remaining           int   // total # of records remaining in the current recordio block.
	prevInt64Value0     int64 // for decoding delta-encoded int
	prevInt64Value1     int64
	prevString          []byte // for decoding prefix-delta-encoded string.
	tmpAuxMd            AuxMetadata
}

func (rb *fieldReadBuf) reset(index biopb.PAMBlockIndexEntry, buf []byte, blob []byte) {
	rb.index = index
	rb.remaining = int(index.NumRecords)

	rb.defaultBuf.n = 0
	rb.defaultBuf.buf = buf
	rb.blobBuf.n = 0
	rb.blobBuf.buf = blob
	rb.prevInt64Value0 = 0
	rb.prevInt64Value1 = 0
	rb.prevString = rb.prevString[:0]
}

type StringDeltaMetadata struct {
	PrefixLen int // Prefix shared with the prev record
	DeltaLen  int // Length of the suffix that differs from the prev record.
}

// ReadStringDeltaMetadata reads the length information for delta-encoded
// string.  Pass the result to readStringDeltaField() to actually decode the
// string.
func (fr *Reader) ReadStringDeltaMetadata() StringDeltaMetadata {
	rb := &fr.fb
	return StringDeltaMetadata{
		PrefixLen: int(rb.defaultBuf.Uvarint32()),
		DeltaLen:  int(rb.defaultBuf.Uvarint32()),
	}
}

// SkipStringDeltaField skips a delta-encoded string.
func (fr *Reader) SkipStringDeltaField() {
	rb := &fr.fb
	rb.remaining--
	md := fr.ReadStringDeltaMetadata()
	prefix := rb.prevString[:md.PrefixLen]
	resizeBuf(&rb.prevString, md.PrefixLen+md.DeltaLen)
	copy(rb.prevString, prefix)
	copy(rb.prevString[md.PrefixLen:], rb.blobBuf.RawBytes(md.DeltaLen))
}

// ReadStringDeltaField reads a delta-encoded string. The function call must be
// preceded by a call to ReadSTringDeltaMetadata, which returns the value of
// "md".
func (fr *Reader) ReadStringDeltaField(md StringDeltaMetadata, arena *UnsafeArena) string {
	rb := &fr.fb
	rb.remaining--
	if md.PrefixLen < 0 {
		log.Panic(md)
	}
	destBuf := arena.Alloc(md.PrefixLen + md.DeltaLen)
	if md.PrefixLen > 0 {
		copy(destBuf, rb.prevString[:md.PrefixLen])
	}
	copy(destBuf[md.PrefixLen:], rb.blobBuf.RawBytes(md.DeltaLen))

	resizeBuf(&rb.prevString, len(destBuf))
	copy(rb.prevString, destBuf)
	return gunsafe.BytesToString(destBuf)
}

// ReadVarintDeltaField reads a field containing a delta-encoded int.
func (fr *Reader) ReadVarintDeltaField() int64 {
	rb := &fr.fb
	rb.remaining--
	delta := rb.defaultBuf.Varint64()
	value := rb.prevInt64Value0 + delta
	rb.prevInt64Value0 = value
	return value
}

// ReadVarintField reads a field containing a varint.
func (fr *Reader) ReadVarintField() int64 {
	rb := &fr.fb
	rb.remaining--
	return rb.defaultBuf.Varint64()
}

// ReadMapqField reads a mapq value.
func (fr *Reader) ReadMapqField() byte {
	rb := &fr.fb
	rb.remaining--
	return rb.defaultBuf.Uint8()
}

// ReadFlagsField reads a flags value.
func (fr *Reader) ReadFlagsField() sam.Flags {
	rb := &fr.fb
	rb.remaining--
	return sam.Flags(rb.defaultBuf.Uint16())
}

// Read a block from recordio and uncompress it.
func (fr *Reader) readBlock(fileOff int64) error {
	fb := &fr.fb

	fr.rio.Seek(recordio.ItemLocation{uint64(fileOff), 0})
	if !fr.rio.Scan() {
		err := fr.rio.Err()
		if err == nil {
			err = fmt.Errorf("Failed to read a block at offset %d", fileOff)
		}
		return err
	}
	fb.buf = fr.rio.Get().([]byte)
	var err error
	fb.header, err = readBlockHeader(&fb.buf)
	return err
}

// readBlock reads a set of recordio blocks listed in "addr", uncompresses them,
// and generate sam.Records. Returns false on EOF. An error is reported in
// Reader.err.
func (fr *Reader) readNextBlock() bool {
	if len(fr.blocks) == 0 {
		return false
	}
	addr := fr.blocks[0]
	fr.blocks = fr.blocks[1:]

	// Read and uncompress the recordio block.
	if err := fr.readBlock(int64(addr.FileOffset)); err != nil {
		fr.err.Set(err)
		return false
	}
	// Set up the read pointers
	fb := &fr.fb
	limitOffset := uint32(len(fb.buf))
	if fb.header.Offset > fb.header.BlobOffset || fb.header.BlobOffset > limitOffset {
		log.Panic(fb)
	}
	fb.reset(addr,
		fb.buf[fb.header.Offset:fb.header.BlobOffset],
		fb.buf[fb.header.BlobOffset:limitOffset])
	log.Debug.Printf("%v: Read block %+v, %d remaining", fr.label, addr, len(fr.blocks))
	return true
}

// SkipCigarField skips the next cigar field.
func (fr *Reader) SkipCigarField() {
	rb := &fr.fb
	rb.remaining--
	nOps := fr.ReadCigarLen()
	for i := 0; i < nOps; i++ {
		rb.defaultBuf.Uvarint32()
	}
}

// ReadCigarLen reads the the # of cigar ops.
func (fr *Reader) ReadCigarLen() int {
	return int(fr.fb.defaultBuf.Uvarint32())
}

// ReadCigarField reads next the Cigar field. The function call must be preceded
// by a call to ReadCigarLen, which returns the value of nOps.
func (fr *Reader) ReadCigarField(nOps int, arena *UnsafeArena) sam.Cigar {
	rb := &fr.fb
	rb.remaining--
	cigar := gbam.UnsafeBytesToCigar(arena.Alloc(nOps * 4))
	for i := 0; i < nOps; i++ {
		cigar[i] = sam.CigarOp(rb.defaultBuf.Uvarint32())
	}
	return cigar
}

// SkipSeqField skips the next seq field.
func (fr *Reader) SkipSeqField() {
	rb := &fr.fb
	rb.remaining--
	nBases := int(fr.fb.defaultBuf.Uvarint32())
	bytes := SeqBytes(nBases)
	rb.blobBuf.RawBytes(bytes)
}

// SeqBytes computes the size of a sam.Seq.Seq that stores n bases.  It returns
// ⌈nbases/2⌉, since each base consumes 4 bits.
func SeqBytes(n int) int {
	return (n + 1) / 2
}

// ReadSeqMetadata returns the length of the next seq field.
func (fr *Reader) ReadSeqMetadata() int {
	return int(fr.fb.defaultBuf.Uvarint32())
}

// ReadSeqField reads the Seq field. nBases must be obtained by calling
// ReadSeqMetadata.
func (fr *Reader) ReadSeqField(nBases int, arena *UnsafeArena) sam.Seq {
	rb := &fr.fb
	rb.remaining--
	bytes := SeqBytes(nBases)
	destBuf := arena.Alloc(bytes)
	copy(destBuf, rb.blobBuf.RawBytes(bytes))
	return sam.Seq{
		Length: nBases,
		Seq:    gbam.UnsafeBytesToDoublets(destBuf),
	}
}

// SkipQualField skips the next qual field.
func (fr *Reader) SkipQualField() {
	rb := &fr.fb
	rb.remaining--
	nBases := int(rb.defaultBuf.Uvarint32())
	rb.blobBuf.RawBytes(nBases)
}

// ReadQualMetadata returns the size of the qual field.
func (fr *Reader) ReadQualMetadata() int {
	return int(fr.fb.defaultBuf.Uvarint32())
}

// ReadQualField reads the next qual field. The function call must be preceded
// by a call to ReadQualLen.
func (fr *Reader) ReadQualField(nBases int, arena *UnsafeArena) []byte {
	rb := &fr.fb
	rb.remaining--
	qual := arena.Alloc(nBases)
	copy(qual, rb.blobBuf.RawBytes(nBases))
	return qual
}

type AuxTagHeader struct {
	// Two-letter tag name + datatype ('Z', 'H', 'i', etc)
	Name [3]byte
	// Length of the payload part (excluding the first three letters).
	Len int
}

type AuxMetadata struct {
	Tags []AuxTagHeader
}

// ReadAuxMetadata reads the number and the size information of the aux field.
func (fr *Reader) ReadAuxMetadata() AuxMetadata {
	rb := &fr.fb
	nAux := int(rb.defaultBuf.Uvarint32())
	if cap(rb.tmpAuxMd.Tags) < nAux {
		rb.tmpAuxMd.Tags = make([]AuxTagHeader, nAux)
	} else {
		rb.tmpAuxMd.Tags = rb.tmpAuxMd.Tags[:nAux]
	}
	for i := 0; i < nAux; i++ {
		t := &rb.tmpAuxMd.Tags[i]
		copy(t.Name[:], rb.blobBuf.RawBytes(3))
		switch t.Name[2] {
		case 'A', 'c', 'C': // ascii, int8, uint8
			t.Len = 1
		case 's', 'S': // int16, uint16
			t.Len = 2
		case 'i', 'I', 'f': // int32, uint32, float32
			t.Len = 4
		case 'Z', 'H': // text, hex string
			t.Len = int(rb.defaultBuf.Uvarint32())
		default:
			// TODO(saito) Handle unknown tags more gracefully.
			log.Panicf("Unknown aux tag: %+v", t)
		}
	}
	return rb.tmpAuxMd
}

// SkipAuxField skips the next aux field.
func (fr *Reader) SkipAuxField() {
	rb := &fr.fb
	rb.remaining--
	md := fr.ReadAuxMetadata()
	for _, tag := range md.Tags {
		rb.blobBuf.RawBytes(tag.Len)
	}
}

// SizeofSliceHeader is the internal size of a slice. Usually 2*(CPU word size).
const SizeofSliceHeader = int(unsafe.Sizeof(reflect.SliceHeader{}))

// ReadAuxField reads the next aux field. This function call must be preceded by
// a call to ReadAuxMetadata.
func (fr *Reader) ReadAuxField(md AuxMetadata, arena *UnsafeArena) []sam.Aux {
	rb := &fr.fb
	rb.remaining--
	var aux []sam.Aux
	// Allocate the backing space for aux.
	arena.Align()
	auxBuf := arena.Alloc(len(md.Tags) * SizeofSliceHeader)
	// Clear the array before updating rec.AuxFields. GC will be
	// confused otherwise.
	for i := range auxBuf {
		auxBuf[i] = 0
	}
	auxBufHdr := (*reflect.SliceHeader)(unsafe.Pointer(&auxBuf))
	auxHdr := (*reflect.SliceHeader)(unsafe.Pointer(&aux))
	auxHdr.Data = auxBufHdr.Data
	auxHdr.Len = len(md.Tags)
	auxHdr.Cap = auxHdr.Len

	for i, tag := range md.Tags {
		tagBuf := arena.Alloc(len(tag.Name) + tag.Len)
		copy(tagBuf, tag.Name[:])
		copy(tagBuf[3:], rb.blobBuf.RawBytes(tag.Len))
		aux[i] = sam.Aux(tagBuf)
	}
	return aux
}

// ReadCoordField reads the next coordinate value.
func (fr *Reader) ReadCoordField() biopb.Coord {
	rb := &fr.fb
	rb.remaining--
	refID := rb.prevInt64Value0 + rb.defaultBuf.Varint64()
	rb.prevInt64Value0 = refID
	pos := rb.prevInt64Value1 + rb.blobBuf.Varint64()
	rb.prevInt64Value1 = pos
	return fr.addrGenerator.Generate(int32(refID), int32(pos))
}

// PeekCoordField reads the next coordinate value without advancing the read
// pointer.
func (fr *Reader) PeekCoordField() biopb.Coord {
	rb := &fr.fb
	s0 := rb.defaultBuf.n
	s1 := rb.blobBuf.n
	refID := rb.prevInt64Value0 + rb.defaultBuf.Varint64()
	pos := rb.prevInt64Value1 + rb.blobBuf.Varint64()
	rb.defaultBuf.n = s0
	rb.blobBuf.n = s1

	save := fr.addrGenerator
	coord := fr.addrGenerator.Generate(int32(refID), int32(pos))
	fr.addrGenerator = save
	return coord
}

// BlockStartAddr returns the coordinate of the first record in the current
// block.
func (fr *Reader) BlockStartAddr() biopb.Coord {
	return fr.fb.index.StartAddr
}

// MaybeReadNextBlock reads the next recordio block if the current block had
// been fully consumed, or no block has ever been read.
func (fr *Reader) MaybeReadNextBlock() bool {
	if fr.fb.remaining > 0 {
		return true
	}
	if !fr.readNextBlock() {
		return false
	}
	if fr.field == gbam.FieldCoord {
		start := fr.fb.index.StartAddr
		fr.addrGenerator.LastRec = biopb.Coord{start.RefId, start.Pos, start.Seq - 1}
	}
	return true
}

// Field returns the field type of the reader.
func (fr *Reader) Field() gbam.FieldType { return fr.field }

var (
	dummyMu   sync.Mutex
	dummyQual unsafe.Pointer // stores *[]byte
	dummySeq  unsafe.Pointer // stores *[]sam.Doublet
)

// GetDummyQual returns a singleton qual string that contains dummy data.
func GetDummyQual(length int) []byte {
	buf := (*[]byte)(atomic.LoadPointer(&dummyQual))
	if buf != nil && len(*buf) >= length {
		return (*buf)[:length]
	}
	dummyMu.Lock()
	newBuf := make([]byte, length)
	for i := 0; i < length; i++ {
		newBuf[i] = 0xff
	}
	atomic.StorePointer(&dummyQual, unsafe.Pointer(&newBuf))
	dummyMu.Unlock()
	return newBuf[:length]
}

// GetDummySeq returns a singleton seq string that contains dummy data.
func GetDummySeq(length int) sam.Seq {
	n := (length + 1) / 2
	buf := (*[]sam.Doublet)(atomic.LoadPointer(&dummySeq))
	if buf != nil && len(*buf) >= n {
		return sam.Seq{Length: length, Seq: (*buf)[:n]}
	}
	dummyMu.Lock()
	newBuf := make([]sam.Doublet, length)
	for i := 0; i < length; i++ {
		newBuf[i] = 0xf
	}
	atomic.StorePointer(&dummySeq, unsafe.Pointer(&newBuf))
	dummyMu.Unlock()
	return sam.Seq{Length: length, Seq: newBuf[:n]}
}

func readBlockHeader(buf *[]byte) (biopb.PAMBlockHeader, error) {
	headerSize, n := binary.Varint(*buf)
	if n <= 0 {
		err := fmt.Errorf("Failed to read the block header size")
		log.Error.Print(err)
		return biopb.PAMBlockHeader{}, err
	}
	*buf = (*buf)[n:]
	// TODO(saito): range check
	headerBytes := (*buf)[:headerSize]
	*buf = (*buf)[headerSize:]

	bh := biopb.PAMBlockHeader{}
	err := bh.Unmarshal(headerBytes)
	if err != nil {
		panic(err)
	}
	return bh, nil
}

// Index returns the index of this field.
func (fr *Reader) Index() *biopb.PAMFieldIndex { return &fr.index }

// Label returns the diagnostic label of the reader object.
func (fr *Reader) Label() string { return fr.label }

// Close closes the reader.  Errors are reported through fr.err.
func (fr *Reader) Close() {
	ctx := vcontext.Background()
	if fr.rio != nil { // fr.rio =nil on error
		fr.err.Set(fr.rio.Finish())
	}
	if fr.in != nil { // fr.in =nil on error
		fr.err.Set(fr.in.Close(ctx))
	}
}

// Seek sets up the reader to read the requested coordinate range. Since the
// requestedRange may not be exactly aligned with recordio block boundaries,
// this method will typically arrange to read a slightly wider range than
// requested.  It returns the start coordinate of the first recordio block to be
// read.
//
// REQUIRES: MaybeReadNextBlock has never been called.
func (fr *Reader) Seek(requestedRange biopb.CoordRange) (biopb.Coord, bool) {
	fr.blocks = nil
	for _, b := range fr.Index().Blocks {
		if pamutil.BlockIntersectsRange(b.StartAddr, b.EndAddr, requestedRange) {
			fr.blocks = append(fr.blocks, b)
		}
	}
	if len(fr.blocks) == 0 {
		// There's no record to be read in the range.  We'll report EOF when
		// reading later. Usually, if fr.blocks is empty for one field, it's
		// empty for any other field too.
		return biopb.Coord{}, false
	}
	return fr.blocks[0].StartAddr, true
}
