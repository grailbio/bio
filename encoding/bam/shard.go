// Copyright 2018 GRAIL, Inc. All rights reserved.
// Use of this source code is governed by the Apache 2.0
// license that can be found in the LICENSE file.

package bam

import (
	"fmt"
	"math"
	"strings"

	"github.com/grailbio/base/file"
	"github.com/grailbio/base/vcontext"
	"github.com/grailbio/bio/biopb"
	biogobam "github.com/grailbio/hts/bam"
	"github.com/grailbio/hts/bgzf"
	"github.com/grailbio/hts/sam"
	"v.io/x/lib/vlog"
)

// UniversalRange is a range that covers all possible records.
var UniversalRange = biopb.CoordRange{
	Start: biopb.Coord{0, 0, 0},
	Limit: biopb.Coord{biopb.InfinityRefID, biopb.InfinityPos, 0},
}

// MappedRange is a range that covers all mapped records.
var MappedRange = biopb.CoordRange{
	Start: biopb.Coord{0, 0, 0},
	Limit: biopb.Coord{biopb.LimitValidRefID, biopb.InfinityPos, 0},
}

// Shard represents a genomic interval. The <StartRef,Start,StartSeq> and
// <EndRef,End,EndSeq> coordinates form a half-open, 0-based interval. An
// iterator for such a range will return reads whose start positions fall within
// that range.
//
// The StartSeq, EndSeq fields are used to distinguish a list of reads that
// start at the same coordinate.  The Nth read that start at the coordinate is
// assigned the seq value of N-1 (assuming N is 1-based). For example, Passing
// range [(startref=10,start=100,startseq=15),
// (limitref=10,limit=100,limitseq=20)] will read 16th to 20th read sequences at
// coordinate (10,100)
//
// Uses of non-zero {Start,End}Seq is supported only in PAM files. For BAM
// files, *Seq must be zero.
//
// An unmapped sequence has coordinate (nil,0,seq), and it is stored after any
// mapped sequence. Thus, a shard that contains an unmapped sequence will have
// EndRef=nil, End=1, EndSeq=0> (in theory, End can be any value > 0, but in
// practice we use End=1).
//
// Padding must be >=0. It expands the read range to [PaddedStart, PaddedEnd),
// where PaddedStart=max(0, Start-Padding) and PaddedEnd=min(EndRef.Len(),
// End+Padding)).  The regions [PaddedStart,Start) and [End,PaddedEnd) are not
// part of the shard, since the padding regions will overlap with another
// Shard's [Start, End).
//
// The Shards are ordered according to the order of the bam input file.
// ShardIdx is an index into that ordering.  The first Shard has index 0, and
// the subsequent shards increment the ShardIdx by one each.
type Shard struct {
	StartRef *sam.Reference
	EndRef   *sam.Reference
	Start    int
	End      int
	StartSeq int
	EndSeq   int

	Padding  int
	ShardIdx int
}

// UniversalShard creates a Shard that covers the entire genome, and unmapped
// reads.
func UniversalShard(header *sam.Header) Shard {
	var startRef *sam.Reference
	if len(header.Refs()) > 0 {
		startRef = header.Refs()[0]
	}
	return Shard{
		StartRef: startRef,
		EndRef:   nil,
		Start:    0,
		End:      math.MaxInt32,
	}
}

// PadStart returns max(s.Start-padding, 0).
func (s *Shard) PadStart(padding int) int {
	return max(0, s.Start-padding)
}

// PaddedStart computes the effective start of the range to read, including
// padding.
func (s *Shard) PaddedStart() int {
	return s.PadStart(s.Padding)
}

// PadEnd end returns min(s.End+padding, length of s.EndRef)
func (s *Shard) PadEnd(padding int) int {
	if s.End == 0 && s.EndSeq == 0 {
		// The shard extends to the end of the previous reference. So PadEnd can
		// stay zero.
		return 0
	}
	if s.EndRef == nil {
		// Unmapped reads are all at position 0, so limit can be any positive value.
		return min(math.MaxInt32, s.End+padding)
	}
	return min(s.EndRef.Len(), s.End+padding)
}

// PaddedEnd computes the effective limit of the range to read, including
// padding.
func (s *Shard) PaddedEnd() int {
	return s.PadEnd(s.Padding)
}

// RecordInShard returns true if r is in s.
func (s *Shard) RecordInShard(r *sam.Record) bool {
	return s.CoordInShard(0, CoordFromSAMRecord(r, 0))
}

// RecordInPaddedShard returns true if r is in s+padding.
func (s *Shard) RecordInPaddedShard(r *sam.Record) bool {
	return s.CoordInShard(s.Padding, CoordFromSAMRecord(r, 0))
}

// CoordInShard returns whether coord is within the shard plus the
// supplied padding (this uses the padding parameter in place of
// s.Padding).
func (s *Shard) CoordInShard(padding int, coord biopb.Coord) bool {
	startCoord := NewCoord(s.StartRef, s.PadStart(padding), 0)
	if coord.LT(startCoord) {
		return false
	}
	endCoord := NewCoord(s.EndRef, s.PadEnd(padding), 0)
	return coord.LT(endCoord)
}

// String returns a debug string for s.
func (s *Shard) String() string {
	return fmt.Sprintf("%d:(%s[%d],%d(%d))-(%s[%d],%d(%d))",
		s.ShardIdx, s.StartRef.Name(), s.StartRef.ID(), s.Start, s.PaddedStart(),
		s.EndRef.Name(), s.EndRef.ID(), s.End, s.PaddedEnd())
}

func min(x, y int) int {
	if y < x {
		return y
	}
	return x
}

func max(x, y int) int {
	if y > x {
		return y
	}
	return x
}

// ShardToCoordRange converts bam.Shard to CoordRange.
func ShardToCoordRange(shard Shard) biopb.CoordRange {
	return biopb.CoordRange{
		biopb.Coord{RefId: int32(shard.StartRef.ID()), Pos: int32(shard.Start), Seq: int32(shard.StartSeq)},
		biopb.Coord{RefId: int32(shard.EndRef.ID()), Pos: int32(shard.End), Seq: int32(shard.EndSeq)},
	}
}

// CoordRangeToShard converts RecRange to bam.Shard.
func CoordRangeToShard(header *sam.Header, r biopb.CoordRange, padding, shardIdx int) Shard {
	var startRef *sam.Reference
	if r.Start.RefId >= 0 {
		startRef = header.Refs()[r.Start.RefId]
	}
	var limitRef *sam.Reference
	var limitPos = int(r.Limit.Pos)
	if r.Limit.RefId >= 0 {
		if n := len(header.Refs()); int(r.Limit.RefId) < n {
			limitRef = header.Refs()[r.Limit.RefId]
			limitPos = int(r.Limit.Pos)
		} else {
			limitRef = header.Refs()[n-1]
			limitPos = limitRef.Len()
		}
	}
	return Shard{
		StartRef: startRef,
		Start:    int(r.Start.Pos),
		StartSeq: int(r.Start.Seq),
		EndRef:   limitRef,
		End:      limitPos,
		EndSeq:   int(r.Limit.Seq),
		Padding:  padding,
		ShardIdx: shardIdx,
	}
}

// CoordGenerator is a helper class for computing the Coord.Seq value from a
// sam.Record. This object must be created per pam shard. Generate() must be
// called for every record that is being read or written to the pam file in
// order.
type CoordGenerator struct {
	LastRec biopb.Coord
}

// NewCoordGenerator creates a new CoordGenerator.
func NewCoordGenerator() CoordGenerator {
	return CoordGenerator{biopb.Coord{RefId: 0, Pos: -1, Seq: 0}}
}

// Generate generates the Coord for the given (refid,pos).
//
// REQUIRES: successive calls to this function must supply a non-decreasing
// sequnece of (ref,pos) values.
func (g *CoordGenerator) Generate(refID, pos int32) biopb.Coord {
	if refID == biopb.InfinityRefID {
		// Pos for unmapped reads are meaningless.  The convention in SAM/BAM is to
		// store -1 as Pos, but we don't use negative positions elsewhere, so we
		// just use 0 as a placeholder.
		pos = 0
	}
	if refID < biopb.UnmappedRefID || pos < 0 {
		vlog.Panicf("Illegal addr: %v %v", refID, pos)
	}
	// See if the (refid,pos) has changed. If not, increment the "Seq" part.
	a := biopb.Coord{RefId: refID, Pos: pos, Seq: 0}
	p := biopb.Coord{RefId: g.LastRec.RefId, Pos: g.LastRec.Pos, Seq: 0}
	cmp := a.Compare(p)
	if cmp < 0 {
		vlog.Panicf("Record coordinate decreased from %+v to %v:%v", g.LastRec, refID, pos)
	}
	if cmp == 0 {
		g.LastRec.Seq++
	} else {
		g.LastRec = a
	}
	return g.LastRec
}

// GenerateFromRecord generates the Coord for the given record.
//
// REQUIRES: successive calls to this function must supply record in
// non-decreasing coordinate order.
func (g *CoordGenerator) GenerateFromRecord(rec *sam.Record) biopb.Coord {
	return g.Generate(int32(rec.Ref.ID()), int32(rec.Pos))
}

// CoordFromSAMRecord computes the biopb.Coord for the given record.  It is a
// shorthand for biopb.CoordFromCoord(rec.Ref, rec.Pos, seq).
func CoordFromSAMRecord(rec *sam.Record, seq int32) biopb.Coord {
	return NewCoord(rec.Ref, rec.Pos, seq)
}

// NewCoord generates biopb.Coord from the given parameters.
func NewCoord(ref *sam.Reference, pos int, seq int32) biopb.Coord {
	a := biopb.Coord{RefId: int32(ref.ID()), Pos: int32(pos), Seq: seq}
	if a.RefId == biopb.InfinityRefID && pos < 0 {
		// Pos for unmapped reads are meaningless.  The convention is to
		// store -1 as Pos, but we don't use negative positions
		// elsewhere, so we just use 0 as a placeholder.
		a.Pos = 0
	}
	return a
}

// NewShardChannel returns a closed channel containing the shards.
func NewShardChannel(shards []Shard) chan Shard {
	shardChan := make(chan Shard, len(shards))
	for _, shard := range shards {
		shardChan <- shard
	}
	close(shardChan)
	return shardChan
}

// GetPositionBasedShards returns a list of shards that cover the
// genome using the specified shard size and padding size.  Return a
// shard for the unmapped && mate-unmapped pairs if includeUnmapped is
// true.
//
// The Shards split the BAM data from the given provider into
// contiguous, non-overlapping genomic intervals (Shards). A SAM
// record is associated with a shard if its alignment start position
// is within the given padding distance of the shard. This means reads
// near shard boundaries may be associated with more than one shard.
func GetPositionBasedShards(header *sam.Header, shardSize int, padding int, includeUnmapped bool) ([]Shard, error) {
	var shards []Shard
	shardIdx := 0
	for _, ref := range header.Refs() {
		var start int
		for start < ref.Len() {
			end := min(start+shardSize, ref.Len())
			shards = append(shards,
				Shard{
					StartRef: ref,
					EndRef:   ref,
					Start:    start,
					End:      end,
					Padding:  padding,
					ShardIdx: shardIdx,
				})
			start += shardSize
			shardIdx++
		}
	}
	if includeUnmapped {
		shards = append(shards,
			Shard{
				StartRef: nil,
				EndRef:   nil,
				Start:    0,
				End:      math.MaxInt32,
				ShardIdx: shardIdx,
			})
	}
	ValidateShardList(header, shards, padding)
	return shards, nil
}

// GetByteBasedShards returns a list of shards much like
// GetPositionBasedShards, but the shards are based on a target
// bytesPerShard, and a minimum number of bases pershard (minBases).
// baiPath can point to a traditional style .bai index, or a new style
// .gbai index.
func GetByteBasedShards(bamPath, baiPath string, bytesPerShard int64, minBases, padding int, includeUnmapped bool) (shards []Shard, err error) {
	// TODO(saito) pass the context explicitly.
	ctx := vcontext.Background()
	var bamIn file.File
	bamIn, err = file.Open(ctx, bamPath)
	if err != nil {
		return nil, err
	}
	defer file.CloseAndReport(ctx, bamIn, &err)
	var bamr *biogobam.Reader
	bamr, err = biogobam.NewReader(bamIn.Reader(ctx), 1)
	if err != nil {
		return nil, err
	}
	header := bamr.Header()

	if strings.HasSuffix(baiPath, ".gbai") {
		return gbaiByteBasedShards(header, baiPath, bytesPerShard, minBases, padding, includeUnmapped)
	}
	return baiByteBasedShards(bamr, baiPath, bytesPerShard, minBases, padding, includeUnmapped)
}

// baiByteBasedShards creates shards based on a traditional .bai style index.
func baiByteBasedShards(bamr *biogobam.Reader, baiPath string, bytesPerShard int64,
	minBases, padding int, includeUnmapped bool) (shards []Shard, err error) {
	type boundary struct {
		pos     int32
		filePos int64
	}

	// Get chunks from the .bai file.
	ctx := vcontext.Background()
	var indexIn file.File
	indexIn, err = file.Open(ctx, baiPath)
	if err != nil {
		return nil, err
	}
	defer file.CloseAndReport(ctx, indexIn, &err)
	var index *Index
	index, err = ReadIndex(indexIn.Reader(ctx))
	if err != nil {
		return nil, err
	}
	chunksByRef := index.AllOffsets()
	if len(chunksByRef) <= 0 {
		return nil, fmt.Errorf("%v: no chunks found in the index", baiPath)
	}

	// Compute shards
	shards = []Shard{}
	for refId := 0; refId < len(chunksByRef); refId++ {
		ref := bamr.Header().Refs()[refId]
		refLen := ref.Len()
		offsets := chunksByRef[refId]

		// Pick initial shard boundaries based on bytesPerShard.
		boundaries := []boundary{}
		prevFilePos := int64(0)
		for i, offset := range offsets {
			if i == 0 || (offset.File-prevFilePos) > bytesPerShard {
				var rec biopb.Coord
				var coordErr error
				rec, coordErr = GetCoordAtOffset(bamr, offset)
				if coordErr != nil {
					vlog.Panic(coordErr)
				}
				if int(rec.RefId) != refId {
					vlog.VI(1).Infof("No more reads in refid %d %s, filepos %d", refId, ref.Name(), offset.File)
					break
				}

				if i == 0 {
					boundaries = append(boundaries, boundary{0, offset.File})
				} else {
					boundaries = append(boundaries, boundary{rec.Pos, offset.File})
				}
				prevFilePos = offset.File
			}
		}
		// If there's a next ref following this ref, then add the next
		// ref's start location as the end of this ref.  That way, we
		// will know how many bytes are in the final chunk of this ref,
		// so that we might break up that final chunk.
		if refId < len(chunksByRef)-1 {
			nextOffsets := chunksByRef[refId+1]
			if len(nextOffsets) > 0 {
				vlog.VI(3).Infof("Adding final boundary %v", boundary{int32(refLen), nextOffsets[0].File})
				boundaries = append(boundaries, boundary{int32(refLen), nextOffsets[0].File})
			}
		}

		// Some shards might be too big since the index does not cover
		// the bam file at regular intervals, so further break up
		// large shards based on genomic position.
		boundaries2 := []boundary{boundary{0, -1}}
		for i := 1; i < len(boundaries); i++ {
			if (boundaries[i].filePos - boundaries[i-1].filePos) > 2*bytesPerShard {
				genomeSubdivisions := ((boundaries[i].pos - boundaries[i-1].pos) / int32(minBases)) - 1
				bytesSubdivisions := ((boundaries[i].filePos - boundaries[i-1].filePos) / bytesPerShard) - 1
				var subdivisions int32
				if int64(genomeSubdivisions) < bytesSubdivisions {
					subdivisions = genomeSubdivisions
				} else {
					subdivisions = int32(bytesSubdivisions)
				}
				if subdivisions < 0 {
					subdivisions = 0
				}

				for s := int32(1); s <= subdivisions; s++ {
					subpos := boundaries[i-1].pos + s*((boundaries[i].pos-boundaries[i-1].pos)/(subdivisions+1))
					boundaries2 = append(boundaries2, boundary{subpos, -1})
				}
			}
			boundaries2 = append(boundaries2, boundaries[i])
		}
		if boundaries2[len(boundaries2)-1].pos != int32(refLen) {
			boundaries2 = append(boundaries2, boundary{int32(refLen), -1})
		}

		// Some shards might be smaller than minBases, so break up those shards.
		boundaries3 := []boundary{boundary{0, -1}}
		for i := 1; i < len(boundaries2); i++ {
			if boundaries2[i].pos-boundaries3[len(boundaries3)-1].pos >= int32(minBases) || i == len(boundaries2)-1 {
				boundaries3 = append(boundaries3, boundaries2[i])
			} else {
				vlog.VI(3).Infof("dropping boundary %v (min bases %d)", boundaries2[i], int32(minBases))
			}
		}

		// Convert boundaries to shards.
		for i := 0; i < len(boundaries3)-1; i++ {
			start := int(boundaries3[i].pos)
			end := int(boundaries3[i+1].pos)
			shards = append(shards, Shard{
				StartRef: ref,
				EndRef:   ref,
				Start:    start,
				End:      end,
				Padding:  padding,
				ShardIdx: len(shards),
			})
		}
	}
	if includeUnmapped {
		shards = append(shards, Shard{
			End:      math.MaxInt32,
			ShardIdx: len(shards),
		})
	}
	ValidateShardList(bamr.Header(), shards, padding)
	return shards, nil
}

// gbaiByteBasedShards creates shards based on a newer .gbai style index.
func gbaiByteBasedShards(header *sam.Header, gbaiPath string, bytesPerShard int64,
	minBases, padding int, includeUnmapped bool) (shards []Shard, err error) {
	ctx := vcontext.Background()
	var indexIn file.File
	indexIn, err = file.Open(ctx, gbaiPath)
	if err != nil {
		return nil, err
	}
	defer file.CloseAndReport(ctx, indexIn, &err)
	var index *GIndex
	index, err = ReadGIndex(indexIn.Reader(ctx))
	if err != nil {
		return nil, err
	}

	shards = []Shard{}
	prevRefID := int32(0)
	prevRefPos := int32(0)
	prevFilePos := uint64(0)
	for i, entry := range *index {
		if i == 0 {
			prevRefID = entry.RefID
			prevRefPos = 0
			prevFilePos = entry.VOffset >> 16
		}

		if entry.RefID != prevRefID {
			// Close out the previous reference.
			shards = append(shards, Shard{
				StartRef: header.Refs()[prevRefID],
				EndRef:   header.Refs()[prevRefID],
				Start:    int(prevRefPos),
				End:      header.Refs()[prevRefID].Len(),
				Padding:  padding,
				ShardIdx: len(shards),
			})
			prevRefID = entry.RefID
			prevRefPos = 0
			prevFilePos = entry.VOffset >> 16
		} else if (entry.Pos-prevRefPos) >= int32(minBases) &&
			// Close out the previous shard.
			((entry.VOffset>>16)-prevFilePos) >= uint64(bytesPerShard) {
			shards = append(shards, Shard{
				StartRef: header.Refs()[prevRefID],
				EndRef:   header.Refs()[prevRefID],
				Start:    int(prevRefPos),
				End:      int(entry.Pos),
				Padding:  padding,
				ShardIdx: len(shards),
			})
			prevRefPos = entry.Pos
			prevFilePos = entry.VOffset >> 16
		}

		if i == len(*index)-1 && prevRefID != -1 {
			// Close out the last shard.
			shards = append(shards, Shard{
				StartRef: header.Refs()[prevRefID],
				EndRef:   header.Refs()[prevRefID],
				Start:    int(prevRefPos),
				End:      header.Refs()[prevRefID].Len(),
				Padding:  padding,
				ShardIdx: len(shards),
			})
			prevRefID = entry.RefID
			prevRefPos = 0
			prevFilePos = entry.VOffset >> 16
		}
	}
	if includeUnmapped {
		shards = append(shards, Shard{
			End:      math.MaxInt32,
			ShardIdx: len(shards),
		})
	}
	ValidateShardList(header, shards, padding)
	return shards, nil
}

// ValidateShardList validates that shardList has sensible values. Exposed only for testing.
func ValidateShardList(header *sam.Header, shardList []Shard, padding int) {
	var prevRef *sam.Reference
	for i, shard := range shardList {
		if shard.Start >= shard.End {
			vlog.Panicf("Shard start must precede end for ref %s: %d, %d", shard.StartRef.Name(), shard.Start, shard.End)
		}

		if shard.StartRef == nil {
			if i == len(shardList)-1 {
				continue
			}
			vlog.Panicf("Only the last shard may have nil Ref, not shard %d", i)
		}

		if i == 0 || shard.StartRef != prevRef {
			prevRef = shard.StartRef
			if shard.Start != 0 {
				vlog.Panicf("First shard of ref %s should start at 0, not %d", shard.StartRef.Name(), shard.Start)
			}
		} else {
			if shard.Start != shardList[i-1].End {
				vlog.Panicf("Shard gap for ref %s between %d and %d", shard.StartRef.Name(), shardList[i-1].End, shard.Start)
			}
		}
		if i < len(shardList)-1 && shardList[i+1].StartRef != shard.StartRef && shard.End != shard.StartRef.Len() {
			vlog.Panicf("Last shard of %s should end at reference end: %d, %d", shard.StartRef.Name(), shard.End, shard.StartRef.Len())
		}

		if shard.Padding < 0 {
			vlog.Panicf("Padding must be non-negative: %d", shard.Padding)
		}
	}
}

const (
	infinityRefID = -1
)

// GetCoordAtOffset starts reading BAM from "off", and finds the first place
// where the read position increases. It returns the record
// coordinate. Coord.Seq field is always zero.
func GetCoordAtOffset(bamReader *biogobam.Reader, off bgzf.Offset) (biopb.Coord, error) {
	if off.File == 0 && off.Block == 0 {
		return biopb.Coord{RefId: 0, Pos: 0}, nil
	}
	if err := bamReader.Seek(off); err != nil {
		return biopb.Coord{}, err
	}
	rec, err := bamReader.Read()
	if err != nil {
		return biopb.Coord{}, err
	}
	c := bamReader.LastChunk()
	if c.Begin.File != off.File || c.Begin.Block != off.Block {
		err := fmt.Errorf("Corrupt BAM index %+v, bam reader offset: %+v", c, off)
		vlog.Error(err)
		return biopb.Coord{}, err
	}
	if rec.Ref.ID() > math.MaxInt32 || rec.Pos > math.MaxInt32 {
		return biopb.Coord{}, fmt.Errorf("Read coord does not fit in int32 for %v", rec)
	}
	addr := biopb.Coord{RefId: int32(rec.Ref.ID()), Pos: int32(rec.Pos)}
	if addr.RefId == infinityRefID {
		// Pos for unmapped reads are meaningless.  The convention is to
		// store -1 as Pos, but we don't use negative positions
		// elsewhere, so we just use 0 as a placeholder.
		addr.Pos = 0
	}
	return addr, nil
}
