// Copyright 2018 GRAIL, Inc. All rights reserved.
// Use of this source code is governed by the Apache 2.0
// license that can be found in the LICENSE file.

//go:generate protoc -I. -I../../../vendor -I../../../vendor/github.com/gogo/protobuf/protobuf --gogofaster_out=. pam.proto

package pam

import (
	"fmt"
	"path/filepath"
	"sync"
	"sync/atomic"
	"unsafe"

	"github.com/biogo/hts/sam"
	"github.com/grailbio/base/errorreporter"
	"github.com/grailbio/base/file"
	"github.com/grailbio/base/recordio"
	"github.com/grailbio/base/recordio/recordiozstd"
	"github.com/grailbio/base/vcontext"
	"github.com/grailbio/bio/biopb"
	gbam "github.com/grailbio/bio/encoding/bam"
	"github.com/grailbio/bio/encoding/pam/fieldio"
	"github.com/grailbio/bio/encoding/pam/pamutil"
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

// ShardReader is for reading one PAM rowshard. This class is generally hidden
// behind pam.Reader and is not used directly.
type ShardReader struct {
	label string // for vlogging only.

	// Requested range of records to read. It is a copy of ReadOpts.Range.
	// It may span outside the range of this shard.
	//
	// INVARIANT: (requestedRange ∩ shardRange) != ∅
	requestedRange biopb.CoordRange

	// Fields to read. It is a complement of ReadOpts.DropFields.
	needField [gbam.NumFields]bool

	// Reader for each field. nil if !needField[f]
	fieldReaders [gbam.NumFields]*fieldio.Reader
	index        biopb.PAMShardIndex
	header       *sam.Header      // Parsed out of index.EmbeddedBamHeader
	path         string           // PAM directory path.
	shardRange   biopb.CoordRange // row range parsed out of the filename.
	nRecords     int              // # records read so far
	err          *errorreporter.T // Points to Reader.err
}

var (
	dummyMu   sync.Mutex
	dummyQual unsafe.Pointer // stores *[]byte
	dummySeq  unsafe.Pointer // stores *[]sam.Doublet
)

// GetDummyQual returns a singleton qual string that contains dummy data.  It is
// returned when the qual field is dropped in the options.
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

// GetDummySeq returns a singleton seq string that contains dummy data. It is
// returned when the seq field is dropped in the options.
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

// Read one record from the buffer "rb". prevRec is used to delta-decode some
// fields.
func (r *ShardReader) readRecord() *sam.Record {
	refs := r.header.Refs()
	rb := gbam.GetFromFreePool()
	rec := gbam.CastUp(rb)

	coord, ok := r.fieldReaders[gbam.FieldCoord].ReadCoordField()
	if !ok {
		return nil
	}
	if coord.RefId >= 0 {
		rec.Ref = refs[coord.RefId]
		rec.Pos = int(coord.Pos)
	} else {
		rec.Ref = nil
		rec.Pos = -1
	}

	if r.needField[gbam.FieldFlags] {
		flags, ok := r.fieldReaders[gbam.FieldFlags].ReadUint16Field()
		if !ok {
			return nil
		}
		rec.Flags = sam.Flags(flags)
	}
	if r.needField[gbam.FieldMapq] {
		rec.MapQ, ok = r.fieldReaders[gbam.FieldMapq].ReadUint8Field()
		if !ok {
			return nil
		}
	}
	if r.needField[gbam.FieldMateRefID] {
		mateRefID, ok := r.fieldReaders[gbam.FieldMateRefID].ReadVarintDeltaField()
		if !ok {
			return nil
		}
		rec.MateRef = nil
		if mateRefID >= 0 {
			rec.MateRef = refs[mateRefID]
		}
	}
	if r.needField[gbam.FieldMatePos] {
		matePos, ok := r.fieldReaders[gbam.FieldMatePos].ReadVarintDeltaField()
		if !ok {
			return nil
		}
		rec.MatePos = int(matePos)
	}
	if r.needField[gbam.FieldTempLen] {
		tempLen, ok := r.fieldReaders[gbam.FieldTempLen].ReadVarintField()
		if !ok {
			return nil
		}
		rec.TempLen = int(tempLen)
	}

	// Collect the length info for variable-length fields.
	arenaBytes := 0
	nCigarOps := 0
	if r.needField[gbam.FieldCigar] {
		if nCigarOps, ok = r.fieldReaders[gbam.FieldCigar].ReadCigarLen(); !ok {
			return nil
		}
		arenaBytes += nCigarOps * gbam.CigarOpSize
	}

	nameMd := fieldio.StringDeltaMetadata{}
	if r.needField[gbam.FieldName] {
		if nameMd, ok = r.fieldReaders[gbam.FieldName].ReadStringDeltaMetadata(); !ok {
			return nil
		}
		arenaBytes += nameMd.PrefixLen + nameMd.DeltaLen
	}
	rec.Seq.Length = 0
	if r.needField[gbam.FieldSeq] {
		if rec.Seq.Length, ok = r.fieldReaders[gbam.FieldSeq].ReadSeqMetadata(); !ok {
			return nil
		}
		arenaBytes += fieldio.SeqBytes(rec.Seq.Length)
	}
	qualLen := 0
	if r.needField[gbam.FieldQual] {
		qualLen, ok = r.fieldReaders[gbam.FieldQual].ReadQualMetadata()
		if !ok {
			return nil
		}
		arenaBytes += qualLen
	}
	auxMd := fieldio.AuxMetadata{}
	if r.needField[gbam.FieldAux] {
		if auxMd, ok = r.fieldReaders[gbam.FieldAux].ReadAuxMetadata(); !ok {
			return nil
		}
		// Round up to the next CPU word boundary, since we will store
		// pointers in the arena.
		const pointerSize = int(unsafe.Sizeof(uintptr(0)))
		arenaBytes += len(auxMd.Tags)*fieldio.SizeofSliceHeader + pointerSize
		for _, tag := range auxMd.Tags {
			arenaBytes += len(tag.Name) + tag.Len
		}
	}
	gbam.ResizeScratch(&rb.Scratch, arenaBytes)
	arena := fieldio.NewUnsafeArena(rb.Scratch)
	if r.needField[gbam.FieldCigar] {
		rec.Cigar = r.fieldReaders[gbam.FieldCigar].ReadCigarField(nCigarOps, &arena)
	}
	if r.needField[gbam.FieldName] {
		rec.Name = r.fieldReaders[gbam.FieldName].ReadStringDeltaField(nameMd, &arena)
	}

	// Qual and Seq must be of the same length, as per BAM file format spec.  So
	// when one of them is not read, fill it with dummy data.
	switch {
	case r.needField[gbam.FieldSeq] && r.needField[gbam.FieldQual]:
		// Common case
		rec.Seq = r.fieldReaders[gbam.FieldSeq].ReadSeqField(rec.Seq.Length, &arena)
		rec.Qual = r.fieldReaders[gbam.FieldQual].ReadQualField(qualLen, &arena)
	case r.needField[gbam.FieldSeq] && !r.needField[gbam.FieldQual]:
		// Fill qual with garbage data w/ the same length as seq
		rec.Seq = r.fieldReaders[gbam.FieldSeq].ReadSeqField(rec.Seq.Length, &arena)
		rec.Qual = GetDummyQual(rec.Seq.Length)
	case !r.needField[gbam.FieldSeq] && r.needField[gbam.FieldQual]:
		// Fill seq with garbage data w/ the same length as qual
		rec.Qual = r.fieldReaders[gbam.FieldQual].ReadQualField(qualLen, &arena)
		rec.Seq = GetDummySeq(len(rec.Qual))
	}

	if r.needField[gbam.FieldAux] {
		rec.AuxFields = r.fieldReaders[gbam.FieldAux].ReadAuxField(auxMd, &arena)
	}
	r.nRecords++
	if coord.LT(r.requestedRange.Start) {
		// This can't happen; seek() should have moved the read pointer >=
		// requestedRange.Start.
		vlog.Fatalf("Record too small: %v, requested %+v", r, r.requestedRange)
	}
	if coord.GE(r.requestedRange.Limit) {
		return nil
	}
	return rec
}

// ReadShardIndex reads the index file, "dir/coordRange.index".
func ReadShardIndex(dir string, recRange biopb.CoordRange) (biopb.PAMShardIndex, error) {
	path := pamutil.ShardIndexPath(dir, recRange)
	var index biopb.PAMShardIndex

	ctx := vcontext.Background()
	in, err := file.Open(ctx, path)
	if err != nil {
		return index, errors.Wrap(err, path)
	}
	defer in.Close(ctx) // nolint: errcheck
	rio := recordio.NewScanner(in.Reader(ctx), recordio.ScannerOpts{})
	defer rio.Finish() // nolint: errcheck
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
	if index.Version != pamutil.DefaultVersion {
		return index, errors.Errorf("Wrong PAM version '%v'; expect '%v'", index.Version, pamutil.DefaultVersion)
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
		blockStart, ok := fr.Seek(requestedRange)
		if !ok {
			// There's no record to be read in the range.  We'll report EOF when
			// reading later. Usually, if fr.blocks is empty for one field, it's
			// empty for any other field too.
			return
		}
		coordRange.Start = coordRange.Start.Min(blockStart)
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
	if _, ok := fr.Seek(coordRange); !ok {
		// This shouldn't happen, unless is the file is corrupt
		err := errors.Errorf("%v: Cannot find blocks for coords in range %+v", r.label, coordRange)
		vlog.Error(err)
		r.err.Set(err)
		return
	}

	// readingField is for eliding calls to addr.GE() below in the fast path.
	var readingField [gbam.NumFields]bool

	// getReader() returns a fieldReader for the given reader, or nil if the field
	// is dropped by the user, or all its data blocks are after "addr"
	getReader := func(f gbam.FieldType, addr biopb.Coord) *fieldio.Reader {
		fr = r.fieldReaders[f]
		if fr != nil {
			if readingField[f] {
				return fr
			}
			blockStartAddr, ok := fr.BlockStartAddr()
			if !ok {
				vlog.Fatalf("%v: EOF while reading %+v", fr.Label(), addr)
			}
			if addr.GE(blockStartAddr) {
				readingField[f] = true
				return fr
			}
		}
		return nil
	}

	// Seek the field pointers to requestedRange.Start
	for {
		fr := r.fieldReaders[gbam.FieldCoord]
		addr, ok := fr.PeekCoordField()
		if !ok {
			// No data to read
			vlog.VI(1).Infof("%v: Reached end of data", fr.Label())
			return
		}
		if addr.GE(requestedRange.Start) {
			return
		}
		fr.ReadCoordField()
		if fr := getReader(gbam.FieldFlags, addr); fr != nil {
			fr.SkipUint16Field()
		}
		if fr := getReader(gbam.FieldMapq, addr); fr != nil {
			fr.SkipUint8Field()
		}
		if fr := getReader(gbam.FieldMateRefID, addr); fr != nil {
			fr.ReadVarintDeltaField()
		}
		if fr := getReader(gbam.FieldMatePos, addr); fr != nil {
			fr.ReadVarintDeltaField()
		}
		if fr := getReader(gbam.FieldTempLen, addr); fr != nil {
			fr.ReadVarintField()
		}
		if fr := getReader(gbam.FieldCigar, addr); fr != nil {
			fr.SkipCigarField()
		}
		if fr := getReader(gbam.FieldName, addr); fr != nil {
			fr.SkipStringDeltaField()
		}
		if fr := getReader(gbam.FieldSeq, addr); fr != nil {
			fr.SkipSeqField()
		}
		if fr := getReader(gbam.FieldQual, addr); fr != nil {
			fr.SkipQualField()
		}
		if fr := getReader(gbam.FieldAux, addr); fr != nil {
			fr.SkipAuxField()
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
		label:          fmt.Sprintf("%s:s%s:u%s", filepath.Base(pamIndex.Dir), pamutil.CoordRangePathString(pamIndex.Range), pamutil.CoordRangePathString(requestedRange)),
		path:           pamIndex.Dir,
		shardRange:     pamIndex.Range,
		requestedRange: requestedRange,
		err:            errp,
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
			path := pamutil.FieldDataPath(pamIndex.Dir, pamIndex.Range, gbam.FieldType(f))
			label := fmt.Sprintf("%s:s%s:u%s(%v)",
				filepath.Base(pamIndex.Dir),
				pamutil.CoordRangePathString(pamIndex.Range),
				pamutil.CoordRangePathString(r.requestedRange),
				gbam.FieldType(f))
			r.fieldReaders[f], err = fieldio.NewReader(path, label, (f == int(gbam.FieldCoord)), errp)
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
	for f := range r.fieldReaders {
		if fr := r.fieldReaders[f]; fr != nil {
			fr.Close()
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
	r.label = fmt.Sprintf("%s:u%s", filepath.Base(dir), pamutil.CoordRangePathString(r.opts.Range))
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
