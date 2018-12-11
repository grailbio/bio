// Copyright 2018 GRAIL, Inc. All rights reserved.
// Use of this source code is governed by the Apache 2.0
// license that can be found in the LICENSE file.

package pam

import (
	"fmt"

	"github.com/grailbio/base/errors"
	"github.com/grailbio/base/file"
	"github.com/grailbio/base/recordio/recordiozstd"
	"github.com/grailbio/base/traverse"
	"github.com/grailbio/base/vcontext"
	"github.com/grailbio/bio/biopb"
	gbam "github.com/grailbio/bio/encoding/bam"
	"github.com/grailbio/bio/encoding/pam/fieldio"
	"github.com/grailbio/bio/encoding/pam/pamutil"
	"github.com/grailbio/hts/sam"
)

const (
	// DefaultMaxBufSize is the default value for WriteOpts.MaxBufSize.
	DefaultMaxBufSize = 8 << 20
	// DefaultWriteParallelism is the default value for WriteOpts.MaxFlushParallelism.
	DefaultWriteParallelism = 2
)

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
	return pamutil.ValidateCoordRange(&o.Range)
}

// Writer is a class for generating a PAM rowshard.
type Writer struct {
	label string // For vlogging only.
	opts  WriteOpts
	dir   string // Output destination
	index biopb.PAMShardIndex

	addrGenerator gbam.CoordGenerator

	bufPool      *fieldio.WriteBufPool
	fieldWriters [gbam.NumFields]*fieldio.Writer // Writer for each field

	// Value to be assigned to the "seq" field of a new recBlockWriteBuf.
	nextBlockSeq int
	err          errors.Once
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
		w.fieldWriters[gbam.FieldCoord].PutCoordField(addr, r.Ref.ID(), r.Pos)
	}
	if w.fieldWriters[gbam.FieldFlags] != nil {
		w.fieldWriters[gbam.FieldFlags].PutUint16Field(addr, uint16(r.Flags))
	}
	if w.fieldWriters[gbam.FieldMapq] != nil {
		w.fieldWriters[gbam.FieldMapq].PutUint8Field(addr, r.MapQ)
	}
	if w.fieldWriters[gbam.FieldCigar] != nil {
		w.fieldWriters[gbam.FieldCigar].PutCigarField(addr, r.Cigar)
	}
	if w.fieldWriters[gbam.FieldMateRefID] != nil {
		w.fieldWriters[gbam.FieldMateRefID].PutVarintDeltaField(addr, int64(r.MateRef.ID()))
	}
	if w.fieldWriters[gbam.FieldMatePos] != nil {
		w.fieldWriters[gbam.FieldMatePos].PutVarintDeltaField(addr, int64(r.MatePos))
	}
	if w.fieldWriters[gbam.FieldTempLen] != nil {
		w.fieldWriters[gbam.FieldTempLen].PutVarintField(addr, int64(r.TempLen))
	}
	if w.fieldWriters[gbam.FieldName] != nil {
		w.fieldWriters[gbam.FieldName].PutStringDeltaField(addr, r.Name)
	}
	if w.fieldWriters[gbam.FieldSeq] != nil {
		w.fieldWriters[gbam.FieldSeq].PutSeqField(addr, r.Seq)
	}
	if w.fieldWriters[gbam.FieldQual] != nil {
		w.fieldWriters[gbam.FieldQual].PutBytesField(addr, r.Qual)
	}
	if w.fieldWriters[gbam.FieldAux] != nil {
		w.fieldWriters[gbam.FieldAux].PutAuxField(addr, r.AuxFields)
	}
	for _, fw := range w.fieldWriters {
		if fw != nil && fw.BufLen() >= w.opts.MaxBufSize {
			fw.FlushBuf()
			fw.NewBuf()
		}
	}
}

// Close must be called exactly once. After close, no operation other than Err()
// may be called.
func (w *Writer) Close() error {
	traverse.Each(len(w.fieldWriters), func(i int) error { // nolint: errcheck
		fw := w.fieldWriters[i]
		if fw != nil {
			fw.Close()
		}
		return nil
	})
	w.bufPool.Finish()
	if w.err.Err() != nil {
		return w.err.Err()
	}
	return pamutil.WriteShardIndex(vcontext.Background(), w.dir, w.opts.Range, &w.index)
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
		err:           errors.Once{},
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

	w.label = fmt.Sprintf("%s:%s", dir, pamutil.CoordRangePathString(w.opts.Range))
	w.bufPool = fieldio.NewBufPool(w.opts.WriteParallelism * nWrittenFields)
	w.index = pamutil.NewShardIndex(w.opts.Range, samHeader)
	for f := range w.fieldWriters {
		if dropField[f] {
			continue
		}

		path := pamutil.FieldDataPath(dir, w.opts.Range, gbam.FieldType(f).String())
		label := fmt.Sprintf("%s:%s:%v", file.Base(dir), pamutil.CoordRangePathString(w.opts.Range), gbam.FieldType(f))
		fw := fieldio.NewWriter(path, label, w.opts.Transformers, w.bufPool, &w.err)
		w.fieldWriters[f] = fw
	}
	return w
}

// Err returns any error encountered so far.
func (w *Writer) Err() error {
	return w.err.Err()
}

func init() {
	recordiozstd.Init()
}
