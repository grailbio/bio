//go:generate protoc -I. -I../../../vendor -I../../../vendor/github.com/gogo/protobuf/protobuf --gogofaster_out=. pam.proto

package pam

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
	"github.com/grailbio/base/file"
	"github.com/grailbio/base/recordio"
	"github.com/grailbio/base/recordio/recordiozstd"
	grailunsafe "github.com/grailbio/base/unsafe"
	"github.com/grailbio/base/vcontext"
	"github.com/grailbio/bio/biopb"
	gbam "github.com/grailbio/bio/encoding/bam"
	"github.com/pkg/errors"
	"v.io/x/lib/vlog"
)

// ReadOpts defines configuration parameters for PAM readers.
type ReadOpts struct {
	// DropFields causes the listed fields not to be filled in Read().
	DropFields []gbam.FieldType

	// Optional row shard range. Only records in this range will be returned
	// by Scan() and Read().
	//
	// Examples:
	//   RecRange{{0,0},{InfinityRefID, InfinityPos}} : covers all possible sequences.
	//   RecRange{{0,0},{LimitValidRefID, InfinityPos}} : covers all mapped sequences.
	//   RecRange{{UnmappedRefID,0},{UnmappedRefID, InfinityPos}} : covers all unmapped sequences.
	Range biopb.CoordRange
}

// fieldReader reads a sequence of values for one field type.
type fieldReader struct {
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
}

// Create a new fieldReader. (path,shardRange) defines the path prefix and the
// rowshard range of the data/index files. requestedRange is the coordinate
// range to read, as requested by the user.
func newFieldReader(path string, shardRange biopb.CoordRange, requestedRange biopb.CoordRange, f gbam.FieldType, errp *errorreporter.T) (*fieldReader, error) {
	fr := &fieldReader{
		field: f,
		label: fmt.Sprintf("%s:s%s:u%s(%v)",
			filepath.Base(path), CoordRangePathString(shardRange), CoordRangePathString(requestedRange), f),
		requestedRange: requestedRange,
		err:            errp,
	}
	ctx := vcontext.Background()
	in, err := file.Open(ctx, FieldDataPath(path, shardRange, f))
	if err != nil {
		return fr, errors.Wrapf(err, "%v: Failed to open field data file for %v", path, f)
	}
	fr.in = in
	fr.rin = fr.in.Reader(ctx)
	fr.rio = recordio.NewScanner(fr.rin, recordio.ScannerOpts{})
	trailer := fr.rio.Trailer()
	if len(trailer) == 0 {
		return fr, errors.Errorf("%v: file does not contain an index: %v", path, fr.rio.Err())
	}
	if err := fr.index.Unmarshal(trailer); err != nil {
		return fr, errors.Wrapf(err, "%v: Failed to unmarshal field index for %v", path, f)
	}
	return fr, nil
}

// ShardReader is for reading one PAM rowshard. This class is generally hidden
// behind pam.Reader and is not used directly.
type ShardReader struct {
	label string // for vlogging only.

	// Requested range of records to read. It is a copy of ReadOpts.Range.
	// It may span outside the range of this shard.
	//
	// INVARIANT: intersection of requestedRange shardRange is nonempty.
	requestedRange biopb.CoordRange

	// Fields to read. It is a complement of ReadOpts.DropFields.
	needField [gbam.NumFields]bool

	addrGenerator gbam.CoordGenerator
	// Reader for each field. nil if !needField[f]
	fieldReaders [gbam.NumFields]*fieldReader
	index        biopb.PAMShardIndex
	header       *sam.Header      // Parsed out of index.EmbeddedBamHeader
	path         string           // PAM directory path.
	shardRange   biopb.CoordRange // row range parsed out of the filename.
	nRecords     int              // # records read so far
	err          *errorreporter.T // Points to Reader.err
}

// For parsing values of one field in one recordio block.
type fieldReadBuf struct {
	header              biopb.PAMBlockHeader
	index               biopb.PAMBlockIndexEntry
	buf                 []byte // raw uncompressed bytes
	defaultBuf, blobBuf byteBuffer
	err                 error // any error that happens during read is accumulated here.
	remaining           int   // total # of records remaining in the current recordio block.
	prevInt64Value0     int64 // for decoding delta-encoded int
	prevInt64Value1     int64
	prevString          []byte // for decoding prefix-delta-encoded string.
	tmpAuxMd            auxMetadata
}

func (rb *fieldReadBuf) reset(addr biopb.PAMBlockIndexEntry, buf []byte, blob []byte) {
	rb.index = addr
	rb.remaining = int(addr.NumRecords)

	rb.defaultBuf.n = 0
	rb.defaultBuf.buf = buf
	rb.blobBuf.n = 0
	rb.blobBuf.buf = blob
	rb.err = nil
	rb.prevInt64Value0 = 0
	rb.prevInt64Value1 = 0
	rb.prevString = rb.prevString[:0]
}

type stringDeltaMetadata struct {
	prefixLen int // Prefix shared with the prev record
	deltaLen  int // Length of the suffix that differs from the prev record.
}

// Read the length information for delta-encoded string.  Pass the result to
// readStringDeltaField() to actually decode the string.
func (rb *fieldReadBuf) readStringDeltaMetadata() stringDeltaMetadata {
	return stringDeltaMetadata{
		prefixLen: int(rb.defaultBuf.Uvarint32()),
		deltaLen:  int(rb.defaultBuf.Uvarint32()),
	}
}

func (rb *fieldReadBuf) skipStringDeltaField() {
	rb.remaining--
	md := rb.readStringDeltaMetadata()
	prefix := rb.prevString[:md.prefixLen]
	resizeBuf(&rb.prevString, md.prefixLen+md.deltaLen)
	copy(rb.prevString, prefix)
	copy(rb.prevString[md.prefixLen:], rb.blobBuf.RawBytes(md.deltaLen))
}

// Reconstruct a delta-encoded string. "sd" can be obtained by calling
// readStringDeltaMetadata().
func (rb *fieldReadBuf) readStringDeltaField(md stringDeltaMetadata, arena *unsafeArena) string {
	rb.remaining--
	doassert(md.prefixLen >= 0)
	destBuf := arena.alloc(md.prefixLen + md.deltaLen)
	if md.prefixLen > 0 {
		copy(destBuf, rb.prevString[:md.prefixLen])
	}
	copy(destBuf[md.prefixLen:], rb.blobBuf.RawBytes(md.deltaLen))

	resizeBuf(&rb.prevString, len(destBuf))
	copy(rb.prevString, destBuf)
	return grailunsafe.BytesToString(destBuf)
}

// Consume a field containing a delta-encoded int.
func (rb *fieldReadBuf) readVarintDeltaField() int64 {
	rb.remaining--
	delta := rb.defaultBuf.Varint64()
	value := rb.prevInt64Value0 + delta
	rb.prevInt64Value0 = value
	return value
}

// Consume a field containing a varint.
func (rb *fieldReadBuf) readVarintField() int64 {
	rb.remaining--
	return rb.defaultBuf.Varint64()
}

// Consume the mapq value.
func (rb *fieldReadBuf) readMapqField() byte {
	rb.remaining--
	return rb.defaultBuf.Uint8()
}

// Consume the flags value.
func (rb *fieldReadBuf) readFlagsField() sam.Flags {
	rb.remaining--
	return sam.Flags(rb.defaultBuf.Uint16())
}

// Read a block from recordio and uncompress it.
func (fr *fieldReader) readBlock(fileOff int64) error {
	fb := &fr.fb

	fr.rio.Seek(recordio.ItemLocation{uint64(fileOff), 0})
	if !fr.rio.Scan() {
		err := fr.rio.Err()
		if err == nil {
			err = errors.Errorf("Failed to read a block at offset %d", fileOff)
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
func (fr *fieldReader) readNextBlock() bool {
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
	doassert(fb.header.Offset <= fb.header.BlobOffset)
	doassert(fb.header.BlobOffset <= limitOffset)
	fb.reset(addr,
		fb.buf[fb.header.Offset:fb.header.BlobOffset],
		fb.buf[fb.header.BlobOffset:limitOffset])
	vlog.VI(2).Infof("%v: Read block %+v, %d remaining", fr.label, addr, len(fr.blocks))
	return true
}

func (rb *fieldReadBuf) skipCigarField() {
	rb.remaining--
	nOps := rb.readCigarLen()
	for i := 0; i < nOps; i++ {
		rb.defaultBuf.Uvarint32()
	}
}

// Extract the # of cigar ops.
func (rb *fieldReadBuf) readCigarLen() int {
	return int(rb.defaultBuf.Uvarint32())
}

// Parse the Cigar field.
func (rb *fieldReadBuf) readCigarField(nOps int, arena *unsafeArena) sam.Cigar {
	rb.remaining--
	cigar := gbam.UnsafeBytesToCigar(arena.alloc(nOps * 4))
	for i := 0; i < nOps; i++ {
		cigar[i] = sam.CigarOp(rb.defaultBuf.Uvarint32())
	}
	return cigar
}

func (rb *fieldReadBuf) skipSeqField() {
	rb.remaining--
	nBases := rb.readSeqLen()
	bytes := seqBytes(nBases)
	rb.blobBuf.RawBytes(bytes)
}

// Extract the length of the Seq field.
func (rb *fieldReadBuf) readSeqLen() int {
	return int(rb.defaultBuf.Uvarint32())
}

// Compute the size of a sam.Seq.Seq that stores n bases.  It returns
// ⌈nbases/2⌉, since each base consumes 4 bits.
func seqBytes(n int) int {
	return (n + 1) / 2
}

// Parse the Seq field. nBases must be obtained by calling readSeqLen.
func (rb *fieldReadBuf) readSeqField(nBases int, arena *unsafeArena) sam.Seq {
	rb.remaining--
	bytes := seqBytes(nBases)
	destBuf := arena.alloc(bytes)
	copy(destBuf, rb.blobBuf.RawBytes(bytes))
	return sam.Seq{
		Length: nBases,
		Seq:    gbam.UnsafeBytesToDoublets(destBuf),
	}
}

func (rb *fieldReadBuf) skipQualField() {
	rb.remaining--
	nBases := int(rb.defaultBuf.Uvarint32())
	rb.blobBuf.RawBytes(nBases)
}

func (rb *fieldReadBuf) readQualField(nBases int, arena *unsafeArena) []byte {
	rb.remaining--
	qual := arena.alloc(nBases)
	copy(qual, rb.blobBuf.RawBytes(nBases))
	return qual
}

type auxTagHeader struct {
	// Two-letter tag name + datatype ('Z', 'H', 'i', etc)
	name [3]byte
	// Length of the payload part (excluding the first three letters).
	len int
}

type auxMetadata struct {
	tags []auxTagHeader
}

func (rb *fieldReadBuf) readAuxMetadata() auxMetadata {
	nAux := int(rb.defaultBuf.Uvarint32())
	if cap(rb.tmpAuxMd.tags) < nAux {
		rb.tmpAuxMd.tags = make([]auxTagHeader, nAux)
	} else {
		rb.tmpAuxMd.tags = rb.tmpAuxMd.tags[:nAux]
	}
	for i := 0; i < nAux; i++ {
		t := &rb.tmpAuxMd.tags[i]
		copy(t.name[:], rb.blobBuf.RawBytes(3))
		switch t.name[2] {
		case 'A', 'c', 'C': // ascii, int8, uint8
			t.len = 1
		case 's', 'S': // int16, uint16
			t.len = 2
		case 'i', 'I', 'f': // int32, uint32, float32
			t.len = 4
		case 'Z', 'H': // text, hex string
			t.len = int(rb.defaultBuf.Uvarint32())
		default:
			// TODO(saito) Handle unknown tags more gracefully.
			vlog.Fatalf("Unknown aux tag: %+v", t)
		}
	}
	return rb.tmpAuxMd
}

func (rb *fieldReadBuf) skipAuxField() {
	rb.remaining--
	md := rb.readAuxMetadata()
	for _, tag := range md.tags {
		rb.blobBuf.RawBytes(tag.len)
	}
}

const sizeofSliceHeader = int(unsafe.Sizeof(reflect.SliceHeader{}))

func (rb *fieldReadBuf) readAuxField(md auxMetadata, arena *unsafeArena) []sam.Aux {
	rb.remaining--
	var aux []sam.Aux
	// Allocate the backing space for aux.
	arena.align()
	auxBuf := arena.alloc(len(md.tags) * sizeofSliceHeader)
	// Clear the array before updating rec.AuxFields. GC will be
	// confused otherwise.
	for i := range auxBuf {
		auxBuf[i] = 0
	}
	auxBufHdr := (*reflect.SliceHeader)(unsafe.Pointer(&auxBuf))
	auxHdr := (*reflect.SliceHeader)(unsafe.Pointer(&aux))
	auxHdr.Data = auxBufHdr.Data
	auxHdr.Len = len(md.tags)
	auxHdr.Cap = auxHdr.Len

	for i, tag := range md.tags {
		tagBuf := arena.alloc(len(tag.name) + tag.len)
		copy(tagBuf, tag.name[:])
		copy(tagBuf[3:], rb.blobBuf.RawBytes(tag.len))
		aux[i] = sam.Aux(tagBuf)
	}
	return aux
}

// Read the next coordinate value.
func (rb *fieldReadBuf) readCoordField() (int32, int32) {
	rb.remaining--
	refid := rb.prevInt64Value0 + rb.defaultBuf.Varint64()
	rb.prevInt64Value0 = refid
	pos := rb.prevInt64Value1 + rb.blobBuf.Varint64()
	rb.prevInt64Value1 = pos
	return int32(refid), int32(pos)
}

// Read the next coordinate value without advancing the read pointer.
func (rb *fieldReadBuf) peekCoordField() (int32, int32) {
	s0 := rb.defaultBuf.n
	s1 := rb.blobBuf.n
	refid := rb.prevInt64Value0 + rb.defaultBuf.Varint64()
	pos := rb.prevInt64Value1 + rb.blobBuf.Varint64()
	rb.defaultBuf.n = s0
	rb.blobBuf.n = s1
	return int32(refid), int32(pos)
}

// If the fieldReadBuffer for FieldCoord is empty, read the next recordio block
// and reset r.addr.
func (r *ShardReader) maybeReadNextCoordBlock() bool {
	fr := r.fieldReaders[gbam.FieldCoord]
	if fr.fb.remaining > 0 {
		return true
	}
	if !fr.readNextBlock() {
		return false
	}
	if fr.field == gbam.FieldCoord {
		start := fr.fb.index.StartAddr
		r.addrGenerator.LastRec = biopb.Coord{start.RefId, start.Pos, start.Seq - 1}
	}
	return true
}

// If the fieldReadBuffer for the field is empty, read the next recordio block.
func (fr *fieldReader) maybeReadNextBlock() bool {
	if fr.fb.remaining > 0 {
		return true
	}
	if fr.field == gbam.FieldCoord {
		vlog.Fatal("use maybeReadNextCoordBlock instead")
	}
	return fr.readNextBlock()
}

var (
	dummyMu   sync.Mutex
	dummyQual unsafe.Pointer // stores *[]byte
	dummySeq  unsafe.Pointer // stores *[]sam.Doublet
)

// Return a singleton qual string that contains dummy data.
func getDummyQual(length int) []byte {
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

// Return a singleton seq string that contains dummy data.
func getDummySeq(length int) sam.Seq {
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

// Read one record from the buffer "rb". prevRec is used to delta-decode some
// fields.
func (r *ShardReader) readRecord() *sam.Record {
	refs := r.header.Refs()
	rb := gbam.GetFromFreePool()
	rec := gbam.CastUp(rb)

	if !r.maybeReadNextCoordBlock() {
		return nil
	}
	refID, pos := r.fieldReaders[gbam.FieldCoord].fb.readCoordField()
	rec.Ref = nil
	if refID >= 0 {
		rec.Ref = refs[refID]
	}
	rec.Pos = int(pos)

	if r.needField[gbam.FieldFlags] {
		if !r.fieldReaders[gbam.FieldFlags].maybeReadNextBlock() {
			return nil
		}
		rec.Flags = r.fieldReaders[gbam.FieldFlags].fb.readFlagsField()
	}
	if r.needField[gbam.FieldMapq] {
		if !r.fieldReaders[gbam.FieldMapq].maybeReadNextBlock() {
			return nil
		}
		rec.MapQ = r.fieldReaders[gbam.FieldMapq].fb.readMapqField()
	}
	if r.needField[gbam.FieldMateRefID] {
		if !r.fieldReaders[gbam.FieldMateRefID].maybeReadNextBlock() {
			return nil
		}
		mateRefID := r.fieldReaders[gbam.FieldMateRefID].fb.readVarintDeltaField()
		rec.MateRef = nil
		if mateRefID >= 0 {
			rec.MateRef = refs[mateRefID]
		}
	}
	if r.needField[gbam.FieldMatePos] {
		if !r.fieldReaders[gbam.FieldMatePos].maybeReadNextBlock() {
			return nil
		}
		rec.MatePos = int(r.fieldReaders[gbam.FieldMatePos].fb.readVarintDeltaField())
	}
	if r.needField[gbam.FieldTempLen] {
		if !r.fieldReaders[gbam.FieldTempLen].maybeReadNextBlock() {
			return nil
		}
		rec.TempLen = int(r.fieldReaders[gbam.FieldTempLen].fb.readVarintField())
	}

	// Collect the length info for variable-length fields.
	arenaBytes := 0
	nCigarOps := 0
	if r.needField[gbam.FieldCigar] {
		if !r.fieldReaders[gbam.FieldCigar].maybeReadNextBlock() {
			return nil
		}
		nCigarOps = r.fieldReaders[gbam.FieldCigar].fb.readCigarLen()
		arenaBytes += nCigarOps * gbam.CigarOpSize
	}

	nameMd := stringDeltaMetadata{}
	if r.needField[gbam.FieldName] {
		if !r.fieldReaders[gbam.FieldName].maybeReadNextBlock() {
			return nil
		}
		nameMd = r.fieldReaders[gbam.FieldName].fb.readStringDeltaMetadata()
		arenaBytes += nameMd.prefixLen + nameMd.deltaLen
	}
	rec.Seq.Length = 0
	if r.needField[gbam.FieldSeq] {
		if !r.fieldReaders[gbam.FieldSeq].maybeReadNextBlock() {
			return nil
		}
		rec.Seq.Length = int(r.fieldReaders[gbam.FieldSeq].fb.defaultBuf.Uvarint32())
		arenaBytes += seqBytes(rec.Seq.Length)
	}
	qualLen := 0
	if r.needField[gbam.FieldQual] {
		if !r.fieldReaders[gbam.FieldQual].maybeReadNextBlock() {
			return nil
		}
		qualLen = int(r.fieldReaders[gbam.FieldQual].fb.defaultBuf.Uvarint32())
		arenaBytes += qualLen
	}
	auxMd := auxMetadata{}
	if r.needField[gbam.FieldAux] {
		if !r.fieldReaders[gbam.FieldAux].maybeReadNextBlock() {
			return nil
		}
		const pointerSize = int(unsafe.Sizeof(uintptr(0)))
		auxMd = r.fieldReaders[gbam.FieldAux].fb.readAuxMetadata()
		// Round up to the next CPU word boundary, since we will store
		// pointers in the arena.
		arenaBytes += len(auxMd.tags)*sizeofSliceHeader + pointerSize
		for _, tag := range auxMd.tags {
			arenaBytes += len(tag.name) + tag.len
		}
	}
	gbam.ResizeScratch(&rb.Scratch, arenaBytes)
	arena := unsafeArena{buf: rb.Scratch}
	if r.needField[gbam.FieldCigar] {
		rec.Cigar = r.fieldReaders[gbam.FieldCigar].fb.readCigarField(nCigarOps, &arena)
	}
	if r.needField[gbam.FieldName] {
		rec.Name = r.fieldReaders[gbam.FieldName].fb.readStringDeltaField(nameMd, &arena)
	}

	// Qual and Seq must be of the same length, as per BAM file format spec.  So
	// when one of them is not read, fill it with dummy data.
	switch {
	case r.needField[gbam.FieldSeq] && r.needField[gbam.FieldQual]:
		// Common case
		rec.Seq = r.fieldReaders[gbam.FieldSeq].fb.readSeqField(rec.Seq.Length, &arena)
		rec.Qual = r.fieldReaders[gbam.FieldQual].fb.readQualField(qualLen, &arena)
	case r.needField[gbam.FieldSeq] && !r.needField[gbam.FieldQual]:
		// Fill qual with garbage data w/ the same length as seq
		rec.Seq = r.fieldReaders[gbam.FieldSeq].fb.readSeqField(rec.Seq.Length, &arena)
		rec.Qual = getDummyQual(rec.Seq.Length)
	case !r.needField[gbam.FieldSeq] && r.needField[gbam.FieldQual]:
		// Fill seq with garbage data w/ the same length as qual
		rec.Qual = r.fieldReaders[gbam.FieldQual].fb.readQualField(qualLen, &arena)
		rec.Seq = getDummySeq(len(rec.Qual))
	}

	if r.needField[gbam.FieldAux] {
		rec.AuxFields = r.fieldReaders[gbam.FieldAux].fb.readAuxField(auxMd, &arena)
	}
	addr := r.addrGenerator.GenerateFromRecord(rec)
	r.nRecords++
	if addr.LT(r.requestedRange.Start) {
		// This can't happen; seek() should have moved the read pointer >=
		// requestedRange.Start.
		vlog.Fatalf("Record too small: %v, requested %+v", r, r.requestedRange)
	}
	if addr.GE(r.requestedRange.Limit) {
		return nil
	}
	return rec
}

func readBlockHeader(buf *[]byte) (biopb.PAMBlockHeader, error) {
	headerSize, n := binary.Varint(*buf)
	if n <= 0 {
		err := errors.Errorf("Failed to read the block header size")
		vlog.Error(err)
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

// ReadShardIndex reads the index file, "dir/recRange.index".
func ReadShardIndex(dir string, recRange biopb.CoordRange) (biopb.PAMShardIndex, error) {
	path := ShardIndexPath(dir, recRange)
	var index biopb.PAMShardIndex

	ctx := vcontext.Background()
	in, err := file.Open(ctx, path)
	if err != nil {
		return index, errors.Wrap(err, path)
	}
	defer in.Close(ctx)
	rio := recordio.NewScanner(in.Reader(ctx), recordio.ScannerOpts{})
	defer rio.Finish()
	if !rio.Scan() {
		return index, errors.Errorf("ReadShardIndex %v: Failed to read record: %v", path, rio.Err())
	}
	err = index.Unmarshal(rio.Get().([]byte))
	if err != nil {
		return index, err
	}
	if index.Magic != ShardIndexMagic {
		return index, errors.Errorf("Wrong index version '%v'; expect '%v'", index.Magic, ShardIndexMagic)
	}
	if index.Version != DefaultVersion {
		return index, errors.Errorf("Wrong PAM version '%v'; expect '%v'", index.Version, DefaultVersion)
	}
	return index, rio.Err()
}

func validateFieldIndex(index biopb.PAMFieldIndex) error {
	for _, block := range index.Blocks {
		if block.NumRecords <= 0 {
			return errors.Errorf("Corrupt block index: %+v", block)
		}
	}
	return nil
}

func validateReadOpts(o *ReadOpts) error {
	for _, fi := range o.DropFields {
		if int(fi) < 0 || int(fi) >= gbam.NumFields {
			return errors.Errorf("Invalid DropField %v in %+v", fi, *o)
		}
		if fi == gbam.FieldCoord {
			// Coord field is needed to support range reads.
			return errors.Errorf("Dropping Coord field is not supported in %+v", *o)
		}
	}
	return ValidateCoordRange(&o.Range)
}

// Set up the reader so that next call to readRecord() will read the first
// record at or after requestedRange.Start.
func (r *ShardReader) seek(requestedRange biopb.CoordRange) {
	// For each field (except coord; more on that below), find the subset of index
	// blocks that contain requestedRange.
	coordRange := requestedRange
	for f := int(gbam.FieldCoord + 1); f < gbam.NumFields; f++ {
		fr := r.fieldReaders[f]
		if fr == nil {
			continue
		}
		for _, b := range fr.index.Blocks {
			if blockIntersectsRange(b.StartAddr, b.EndAddr, requestedRange) {
				fr.blocks = append(fr.blocks, b)
			}
		}
		if len(fr.blocks) == 0 {
			// There's no record to be read in the range.  We'll report EOF when
			// reading later. Usually, if fr.blocks is empty for one field, it's
			// empty for any other field too.
			return
		}
		coordRange.Start = coordRange.Start.Min(fr.blocks[0].StartAddr)
	}

	// We need to advance the read pointer of each field to the first record at or
	// after requestedRange.Start. We do the following:
	//
	// 1. Assume that (say) FieldSeq has three recordio blocks {b0, b1, b2}, that
	// intersect with requestedRange.
	//
	// 2. Read the recordio blocks for FieldCoord so that they cover (b0,b1,b2).
	// Then sequentially scan these blocks and find b0.StartAddr.
	//
	// 3. Sequentially scan both FieldCoord and FieldSeq simultaneously, until the
	// the read pointer for FieldCoord is at requestedRange.Start.
	//
	// The below code does this for all the fields in parallel.
	// Read FieldCoord so that it covers all the recordioblocks read by other
	// fields.
	fr := r.fieldReaders[gbam.FieldCoord]
	for _, b := range fr.index.Blocks {
		if blockIntersectsRange(b.StartAddr, b.EndAddr, coordRange) {
			fr.blocks = append(fr.blocks, b)
		}
	}
	if len(fr.blocks) == 0 {
		// This shouldn't happen, unless is the file is corrupt
		err := errors.Errorf("%v: Cannot find blocks for coords in range %+v, index: %+v", fr.label, coordRange, fr.index)
		vlog.Error(err)
		r.err.Set(err)
		return
	}

	// readingField is for eliding calls to addr.GE() below in the fast path.
	var readingField [gbam.NumFields]bool

	// getReader() returns a fieldReader for the given reader, or nil if the field
	// is dropped by the user, or all its data blocks are after "addr"
	getReader := func(f gbam.FieldType, addr biopb.Coord) *fieldReader {
		fr = r.fieldReaders[f]
		if fr != nil {
			if !fr.maybeReadNextBlock() {
				vlog.Fatalf("%v: EOF while reading %+v", fr.label, addr)
			}
			if readingField[f] {
				return fr
			}
			if addr.GE(fr.fb.index.StartAddr) {
				readingField[f] = true
				return fr
			}
		}
		return nil
	}

	// Seek the field pointers to requestedRange.Start
	for {
		if !r.maybeReadNextCoordBlock() {
			// No data to read
			vlog.VI(1).Infof("Reached end of data, %+v", r.addrGenerator)
			return
		}
		fr := r.fieldReaders[gbam.FieldCoord]
		refID, pos := fr.fb.peekCoordField()
		save := r.addrGenerator
		addr := r.addrGenerator.Generate(refID, pos)
		if addr.GE(requestedRange.Start) {
			r.addrGenerator = save
			return
		}
		fr.fb.readCoordField()
		if fr := getReader(gbam.FieldFlags, addr); fr != nil {
			fr.fb.readFlagsField()
		}
		if fr := getReader(gbam.FieldMapq, addr); fr != nil {
			fr.fb.readMapqField()
		}
		if fr := getReader(gbam.FieldMateRefID, addr); fr != nil {
			fr.fb.readVarintDeltaField()
		}
		if fr := getReader(gbam.FieldMatePos, addr); fr != nil {
			fr.fb.readVarintDeltaField()
		}
		if fr := getReader(gbam.FieldTempLen, addr); fr != nil {
			fr.fb.readVarintField()
		}
		if fr := getReader(gbam.FieldCigar, addr); fr != nil {
			fr.fb.skipCigarField()
		}
		if fr := getReader(gbam.FieldName, addr); fr != nil {
			fr.fb.skipStringDeltaField()
		}
		if fr := getReader(gbam.FieldSeq, addr); fr != nil {
			fr.fb.skipSeqField()
		}
		if fr := getReader(gbam.FieldQual, addr); fr != nil {
			fr.fb.skipQualField()
		}
		if fr := getReader(gbam.FieldAux, addr); fr != nil {
			fr.fb.skipAuxField()
		}
	}
}

// NewShardReader creates a reader for a rowshard.  requestedRange and
// dropFields are fields from ReadOpts.  pamIndex is the index file information,
// gleaned from its pathname. muPtr and errPtr are for reporting errors to the
// parent.
//
// REQUIRES: requestedRange ∩ pamIndex.Range != ∅
func NewShardReader(
	requestedRange biopb.CoordRange,
	dropFields []gbam.FieldType,
	pamIndex FileInfo,
	errp *errorreporter.T) *ShardReader {
	r := &ShardReader{
		label:          fmt.Sprintf("%s:s%s:u%s", filepath.Base(pamIndex.Dir), CoordRangePathString(pamIndex.Range), CoordRangePathString(requestedRange)),
		path:           pamIndex.Dir,
		shardRange:     pamIndex.Range,
		requestedRange: requestedRange,
		err:            errp,
		addrGenerator:  gbam.NewCoordGenerator(),
	}
	vlog.VI(1).Infof("%v: NewShardReader", r.label)
	for i := range r.needField {
		r.needField[i] = true
	}
	for _, f := range dropFields {
		r.needField[f] = false
	}
	var err error
	if r.index, err = ReadShardIndex(r.path, r.shardRange); err != nil {
		vlog.Errorf("Failed to read shard index: %v", err)
		r.err.Set(errors.Wrapf(err, "Failed to read shard index for %v", r.path))
		return r
	}

	r.header, err = gbam.UnmarshalHeader(r.index.EncodedBamHeader)
	if err != nil {
		r.err.Set(errors.Wrapf(err, "%v: Failed to decode sam.Header in index", r.path))
		return r
	}
	if !r.requestedRange.Intersects(r.shardRange) {
		vlog.Fatalf("%v: Range doesn't intersect", r.label)
	}

	for f := range r.needField {
		if r.needField[f] {
			r.fieldReaders[f], err = newFieldReader(pamIndex.Dir, pamIndex.Range, r.requestedRange, gbam.FieldType(f), errp)
			if err != nil {
				r.err.Set(err)
				return r
			}
		}
	}
	r.seek(r.requestedRange)
	return r
}

// Close must be called exactly once. After close, no method may be called.
func (r *ShardReader) Close() {
	ctx := vcontext.Background()
	for f := range r.fieldReaders {
		fr := r.fieldReaders[f]
		if fr != nil {
			if fr.rio != nil { // fr.rio =nil on error
				fr.err.Set(fr.rio.Finish())
			}
			if fr.in != nil { // fr.in =nil on error
				r.err.Set(fr.in.Close(ctx))
			}
		}
	}
}

// Reader is the main PAM reader class. It can read across multiple rowshard
// files.
type Reader struct {
	label string // for vlogging only.
	opts  ReadOpts

	// Sorted list of *.index files whose rowrange intersect with
	// opts.Range.
	indexFiles []FileInfo

	// Current rowshard reader.
	r *ShardReader

	// Object returned by Record().
	rec *sam.Record

	numRead int
	err     errorreporter.T
}

// NewReader creates a new Reader.
func NewReader(opts ReadOpts, dir string) *Reader {
	r := &Reader{
		opts: opts,
	}
	r.err.Set(validateReadOpts(&r.opts))
	if r.err.Err() != nil {
		return r
	}
	r.label = fmt.Sprintf("%s:u%s", filepath.Base(dir), CoordRangePathString(r.opts.Range))
	var err error
	if r.indexFiles, err = findIndexFilesInRange(dir, r.opts.Range); err != nil {
		r.err.Set(err)
		return r
	}
	if len(r.indexFiles) == 0 {
		r.err.Set(errors.Errorf("%v: No pam file found for range %+v", dir, r.opts))
		return r
	}
	vlog.VI(1).Infof("Found index files in range %+v: %+v", r.opts.Range, r.indexFiles)
	r.r = NewShardReader(r.opts.Range, r.opts.DropFields, r.indexFiles[0], &r.err)
	r.indexFiles = r.indexFiles[1:]
	return r
}

// Record returns the most recent record read by Scan.
//
// REQUIRES: Scan() returned true.
func (r *Reader) Record() *sam.Record {
	return r.rec
}

// Scan reads the next record. It returns true if a record is found. It returns
// false on EOF or an error. Call Record to get the record, and Err or Close to
// obtain the error code.
//
// Example:
//   r := pam.NewReader(...)
//   for r.Scan() {
//      rec := r.Record()
//      ... use the rec ...
//   }
//   err := r.Close()
func (r *Reader) Scan() bool {
	for {
		if r.Err() != nil {
			// TODO(saito) do lock-free access
			return false
		}
		r.rec = r.r.readRecord()
		if r.rec != nil {
			r.numRead++
			return true
		}
		if r.Err() != nil {
			return false
		}
		// End of shard
		if len(r.indexFiles) == 0 {
			if r.opts.Range.Limit.RefId < 0 {
				vlog.VI(1).Infof("%v: Finished reading %d records for range %v", r.label, r.numRead, r.opts.Range)
			}
			return false
		}
		r.r.Close()
		r.r = NewShardReader(r.opts.Range, r.opts.DropFields, r.indexFiles[0], &r.err)
		r.indexFiles = r.indexFiles[1:]
	}
}

// Err returns any error encountered so far.
//
// Note: Err never returns io.EOF. On EOF, Scan() returns false, and Err()
// returns nil.
func (r *Reader) Err() error {
	return r.err.Err()
}

// Close all the files. Must be called exactly once.
func (r *Reader) Close() error {
	if r.r != nil {
		r.r.Close()
	}
	return r.err.Err()
}

func init() {
	recordiozstd.Init()
}
