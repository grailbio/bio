// Package fieldio provides a reader and a writer for individual column (field).
package fieldio

import (
	"encoding/binary"
	"fmt"
	"io"
	"reflect"
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
	label  string                     // for vlogging only.
	in     file.File                  // For reading the data file.
	rin    io.ReadSeeker              // in.Reader
	rio    recordio.Scanner           // recordio wrapper for "in".
	index  biopb.PAMFieldIndex        // Contents of *<fieldname>.index file.
	blocks []biopb.PAMBlockIndexEntry // Subset of index.Blocks that intersect requestedRange.
	fb     fieldReadBuf               // Current buffer being parsed.
	err    *errorreporter.T

	coordField    bool                // True if the field is gbam.FieldCoord.
	addrGenerator gbam.CoordGenerator // Computes biopb.Coord.Seq. Used only when coordField=true.
}

// NewReader creates a new Reader that reads from the given path. Label is shown
// in log messages. coordField should be true if the file stores the genomic
// coordinate. Setting setting coordField=true enables the codepath that
// computes biopb.Coord.Seq values.
func NewReader(path, label string, coordField bool, errp *errorreporter.T) (*Reader, error) {
	fr := &Reader{
		coordField: coordField,
		label:      label,
		err:        errp,
	}
	ctx := vcontext.Background()
	in, err := file.Open(ctx, path)
	if err != nil {
		return fr, errors.E(err, fmt.Sprintf("%s: Failed to open field file for %s", label, path))
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
		return fr, errors.E(err, fmt.Sprintf("%s: Failed to unmarshal field index for %s", label, path))
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
func (fr *Reader) ReadStringDeltaMetadata() (StringDeltaMetadata, bool) {
	if fr.fb.remaining <= 0 && !fr.readNextBlock() {
		return StringDeltaMetadata{}, false
	}
	rb := &fr.fb
	return StringDeltaMetadata{
		PrefixLen: int(rb.defaultBuf.Uvarint32()),
		DeltaLen:  int(rb.defaultBuf.Uvarint32()),
	}, true
}

// SkipStringDeltaField skips a delta-encoded string.
// It panics on EOF or any error.
func (fr *Reader) SkipStringDeltaField() {
	rb := &fr.fb
	rb.remaining--
	md, ok := fr.ReadStringDeltaMetadata()
	if !ok {
		panic(fr)
	}
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
// It returns false on EOF or any error.
func (fr *Reader) ReadVarintDeltaField() (int64, bool) {
	if fr.fb.remaining <= 0 && !fr.readNextBlock() {
		return 0, false
	}
	rb := &fr.fb
	rb.remaining--
	delta := rb.defaultBuf.Varint64()
	value := rb.prevInt64Value0 + delta
	rb.prevInt64Value0 = value
	return value, true
}

// SkipVarintDeltaField skips a varint-delta-encoded string.
// It panics on EOF or any error.
func (fr *Reader) SkipVarintDeltaField() {
	if _, ok := fr.ReadVarintDeltaField(); !ok {
		panic(fr)
	}
}

// ReadVarintField reads a field containing a varint.
// It returns false on EOF or any error.
func (fr *Reader) ReadVarintField() (int64, bool) {
	if fr.fb.remaining <= 0 && !fr.readNextBlock() {
		return 0, false
	}
	rb := &fr.fb
	rb.remaining--
	return rb.defaultBuf.Varint64(), true
}

// SkipVarintField skips the next varint-encoded field.
// It panics on EOF or any error.
func (fr *Reader) SkipVarintField() {
	if _, ok := fr.ReadVarintField(); !ok {
		panic(fr)
	}
}

// ReadUint8Field reads a mapq value. It returns false on EOF or any error.
func (fr *Reader) ReadUint8Field() (uint8, bool) {
	if fr.fb.remaining <= 0 && !fr.readNextBlock() {
		return 0, false
	}
	rb := &fr.fb
	rb.remaining--
	return rb.defaultBuf.Uint8(), true
}

// SkipUint8Field skips the next uint8 value.  It panics on EOF or any
// error.
func (fr *Reader) SkipUint8Field() {
	if _, ok := fr.ReadUint8Field(); !ok {
		panic(fr)
	}
}

// ReadUint16Field reads a uint16 value. It returns false on EOF or any error.
func (fr *Reader) ReadUint16Field() (uint16, bool) {
	if fr.fb.remaining <= 0 && !fr.readNextBlock() {
		return 0, false
	}
	rb := &fr.fb
	rb.remaining--
	return rb.defaultBuf.Uint16(), true
}

// SkipUint16Field skips the next uint16 value.  It panics on EOF or any
// error.
func (fr *Reader) SkipUint16Field() {
	if _, ok := fr.ReadUint16Field(); !ok {
		panic(fr)
	}
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

	if fr.coordField {
		start := fr.fb.index.StartAddr
		fr.addrGenerator.LastRec = biopb.Coord{start.RefId, start.Pos, start.Seq - 1}
	}
	log.Debug.Printf("%v: Read block %+v, %d remaining", fr.label, addr, len(fr.blocks))
	return true
}

// SkipCigarField skips the next cigar field.
func (fr *Reader) SkipCigarField() {
	rb := &fr.fb
	rb.remaining--
	nOps, ok := fr.ReadCigarLen()
	if !ok {
		panic(fr)
	}
	for i := 0; i < nOps; i++ {
		rb.defaultBuf.Uvarint32()
	}
}

// ReadCigarLen reads the the # of cigar ops.
func (fr *Reader) ReadCigarLen() (int, bool) {
	if fr.fb.remaining <= 0 && !fr.readNextBlock() {
		return 0, false
	}
	return int(fr.fb.defaultBuf.Uvarint32()), true
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
// It panics on EOF or any error.
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
func (fr *Reader) ReadSeqMetadata() (int, bool) {
	if fr.fb.remaining <= 0 && !fr.readNextBlock() {
		return 0, false
	}
	return int(fr.fb.defaultBuf.Uvarint32()), true
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
// It panics on EOF or any error.
func (fr *Reader) SkipQualField() {
	rb := &fr.fb
	rb.remaining--
	nBases := int(rb.defaultBuf.Uvarint32())
	rb.blobBuf.RawBytes(nBases)
}

// ReadQualMetadata returns the size of the qual field.
func (fr *Reader) ReadQualMetadata() (int, bool) {
	if fr.fb.remaining <= 0 && !fr.readNextBlock() {
		return 0, false
	}
	return int(fr.fb.defaultBuf.Uvarint32()), true
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
func (fr *Reader) ReadAuxMetadata() (AuxMetadata, bool) {
	if fr.fb.remaining <= 0 && !fr.readNextBlock() {
		return AuxMetadata{}, false
	}
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
	return rb.tmpAuxMd, true
}

// SkipAuxField skips the next aux field.
// It panics on EOF or any error.
func (fr *Reader) SkipAuxField() {
	rb := &fr.fb
	rb.remaining--
	md, ok := fr.ReadAuxMetadata()
	if !ok {
		panic(fr)
	}
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

// ReadCoordField reads the next coordinate value.  It returns false on EOF or
// any error.
func (fr *Reader) ReadCoordField() (biopb.Coord, bool) {
	if fr.fb.remaining <= 0 && !fr.readNextBlock() {
		return biopb.Coord{}, false
	}
	rb := &fr.fb
	rb.remaining--
	refID := rb.prevInt64Value0 + rb.defaultBuf.Varint64()
	rb.prevInt64Value0 = refID
	pos := rb.prevInt64Value1 + rb.blobBuf.Varint64()
	rb.prevInt64Value1 = pos
	return fr.addrGenerator.Generate(int32(refID), int32(pos)), true
}

// PeekCoordField reads the next coordinate value without advancing the read
// pointer. It returns false on EOF or any error.
func (fr *Reader) PeekCoordField() (biopb.Coord, bool) {
	if fr.fb.remaining <= 0 && !fr.readNextBlock() {
		return biopb.Coord{}, false
	}
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
	return coord, true
}

// BlockStartAddr returns the coordinate of the first record in the current
// block.
func (fr *Reader) BlockStartAddr() (biopb.Coord, bool) {
	if fr.fb.remaining <= 0 && !fr.readNextBlock() {
		return biopb.Coord{}, false
	}
	return fr.fb.index.StartAddr, true
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
// REQUIRES: maybeReadNextBlock has never been called.
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
