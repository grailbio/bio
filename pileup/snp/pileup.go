// Copyright 2020 Grail Inc.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
package snp

import (
	"context"
	"fmt"
	"io/ioutil"
	"os"
	"runtime"
	"strconv"

	"github.com/grailbio/base/log"
	"github.com/grailbio/base/recordio"
	"github.com/grailbio/base/traverse"
	"github.com/grailbio/bio/circular"
	gbam "github.com/grailbio/bio/encoding/bam"
	"github.com/grailbio/bio/encoding/bamprovider"
	"github.com/grailbio/bio/interval"
	"github.com/grailbio/bio/pileup"
	"github.com/grailbio/hts/sam"
)

type Opts struct {
	// Commandline options.
	BedPath      string
	Region       string
	BamIndexPath string
	Clip         int
	Cols         string
	FlagExclude  int
	Mapq         int
	MaxReadLen   int
	MaxReadSpan  int
	MinBagDepth  int
	MinBaseQual  int
	Parallelism  int
	PerStrand    bool
	RemoveSq     bool
	Stitch       bool
	TempDir      string
}

var DefaultOpts = Opts{
	Clip:        0,
	FlagExclude: 0xf00,
	Mapq:        60,
	MaxReadLen:  500,
	MaxReadSpan: 511,
	MinBagDepth: 0,
	MinBaseQual: 0,
	Parallelism: 0,
	PerStrand:   false,
	RemoveSq:    false,
	Stitch:      false,
}

// Problem:
// Given a sorted BAM/PAM, we are interested in counting the number of
// read(-pair)s supporting each base at a given set of positions (usually
// indicated by a BED file).  We also want to track and report secondary
// statistics such as base-quality, fragment length, etc.  Finally, we may need
// to perform some form of read-pair stitching.
//
// Implementation strategy:
// We start with the following assumptions about the input to avoid unnecessary
// complexity:
// 1. No read is longer than -max-read-len (default 500) bases.
// 2. No read is mapped to a reference-genome interval longer than
//    -max-read-span (default 511) bases.  Note that this is incompatible with
//    RNAseq data.
// These allow us to use fixed-size buffers for most of our computation.
//
// More precisely, suppose we are processing a read starting at (0-based)
// position 1000, and we aren't stitching read-pairs together.  Then, we might
// have some intermediate stats for position 1510, but we can't possibly have
// encountered a read overlapping position 1511 yet.  So our "active interval"
// is [1000, 1511); we only need to keep track of stats for 511 consecutive
// positions at a time.  A circular buffer sized to the next-higher power of 2
// (512 in this case) is perfect for this purpose; we store the stats for
// position x in circular buffer row (x % 512).
//
// Now, let's move on to the next read; suppose it starts at position 1050.
// Before processing it, we flush pileup results for positions 1000..1049 to
// disk, since no more reads can overlap those positions.  We also clear those
// rows of the circular buffer.
//
// How does this change when we need to stitch?
//
// - We are no longer able to process the first read of each read-pair until
//   we've seen the other end (unless the two ends' alignments don't overlap,
//   anyway); instead, we save it to a "firstread-table", to be retrieved
//   later.
//
// - The "active interval" must look backwards as well as forwards.  When we're
//   processing a read starting at (0-based) position 1000, we might not be
//   able to safely report the final pileup results for position 490 yet,
//   because there might be a later read starting at position 1000 with a mate
//   at position 490, and that mate is mapped to the interval [490, 1001) so
//   we've deferred its processing until we've seen both ends of that
//   read-pair.  But we can still report the final pileup results for (0-based)
//   position 489, and anything earlier.
//
//   Then, suppose we move on to a read starting at position 1050.  Before
//   flushing results for positions 490..539, it is necessary to process all
//   unpaired entries in the firstread-table starting at position < 540.
//   (These can exist due to BAM filtering by another tool, or this tool's
//   skipping of low-MAPQ reads.)

// These constants refer to the optional output column-sets.
//   DpRef    = DP (depth) column in .ref.tsv.  Includes low-quality bases.
//   DpAlt    = DP column in .alt.tsv.
//   EndDists = Comma-separated distances from nearest read end.
//   Quals    = Comma-separated base-qualities.
//   Fraglens = Comma-separated read lengths if unstitched, fragment length
//              estimates if stitched.
//   Strands  = Comma-separated strands ('+', '-', '.').
//   HighQ    = Number of supporting reads passing the base-quality threshold.
//              Slated for renaming.
//   LowQ     = Currently an all-zero column existing for backward
//              compatibility.  Will be removed.
const (
	colBitDpRef = 1 << iota
	colBitDpAlt
	colBitEndDists
	colBitQuals
	colBitFraglens
	colBitStrands

	colBitHighQ
	colBitLowQ
)

const colPerReadMask = (colBitEndDists | colBitQuals | colBitFraglens | colBitStrands)

var colNameMap = map[string]int{
	"dpref":    colBitDpRef,
	"dpalt":    colBitDpAlt,
	"enddists": colBitEndDists,
	"quals":    colBitQuals,
	"fraglens": colBitFraglens,
	"strands":  colBitStrands,
	"highq":    colBitHighQ,
	"lowq":     colBitLowQ,
}

// Immutable (within each ref) background info needed for both the inner
// compute loop and result reporting.
type refContext struct {
	// refSeq8 is the current reference sequence, translated to
	// A=1/C=2/G=4/T=8/other=15 encoding to enable more efficient comparisons
	// against BAM/PAM Seq fields.
	// (Might change internal representation to A=0/C=1/G=2/T=3/N=4; if this is
	// done, the field will be renamed to refSeq5.)
	refSeq8 []byte
	// refID is the ID of the current reference.
	refID int
	// refName is the name of the current reference.
	refName string
}

type alignedPos struct {
	posInRef  PosType
	posInRead PosType
}

type pileupMutable struct {
	resultRingBuffer []PileupPayload // main ring buffer
	seq8Buf          []byte          // preallocated buffer to simplify stitchedPileStragglerFirstreads
	alignedBaseBufs  [2][]alignedPos // preallocated buffers for alignRelevantBases
	firstReads       firstreadSNPTable
	endMax           PosType // 1 + <last position that has a pileup entry>
	w                recordio.Writer
	writePosScanner  interval.UnionScanner
}

func newPileupMutable(nCirc PosType, maxReadLen int, stitch bool, w *os.File) (pm pileupMutable) {
	pm = pileupMutable{
		resultRingBuffer: make([]PileupPayload, nCirc),
		seq8Buf:          make([]byte, 0, maxReadLen),
		alignedBaseBufs: [2][]alignedPos{
			make([]alignedPos, 0, maxReadLen),
			nil,
		},
	}
	if stitch {
		pm.firstReads = newFirstreadSNPTable(nCirc)
		pm.alignedBaseBufs[1] = make([]alignedPos, 0, maxReadLen)
	}
	if w != nil {
		pm.w = recordio.NewWriter(w, recordio.WriterOpts{
			Marshal:      MarshalPileupRow,
			Transformers: []string{"zstd 1"},
		})
	}
	return
}

// pileupContext contains immutable data needed by the main pileup routine.
type pileupContext struct {
	bedPart       interval.BEDUnion // per-thread BED subset
	clip          int               // number of bases on ends of each read to treat as min-qual
	ignoreStrand  bool              // are we reporting strand in the output?
	minBaseQual   byte
	perReadNeeded bool           // are we reporting comma-separated per-read stats in the output, or are counts enough?
	qpt           *qualPassTable // (R1 base-qual, R2 base-qual) good enough? lookup table
	stitch        bool
}

// addBase performs a pileup update that only requires count-increments.
func (pm *pileupMutable) addBase(circPos, posInRead, isMinus PosType, seq, qual []byte, minBaseQual byte) {
	row := &pm.resultRingBuffer[circPos]
	row.Depth++
	base := pileup.Seq8ToEnumTable[seq[posInRead]]
	// Always count Ns, to preserve tsv-snp2 compatibility.
	if (qual[posInRead] >= minBaseQual) || (base == pileup.BaseX) {
		row.Counts[base][isMinus]++
	}
}

// appendBase performs a more-expensive pileup update that appends a bunch of
// per-read stats.
func (pm *pileupMutable) appendBase(circPos, posInRead, isMinus PosType, seq, qual []byte, minBaseQual, strandByte byte) {
	row := &pm.resultRingBuffer[circPos]
	row.Depth++
	base := pileup.Seq8ToEnumTable[seq[posInRead]]
	if base == pileup.BaseX {
		row.Counts[base][isMinus]++
	} else if qual[posInRead] >= minBaseQual {
		row.Counts[base][isMinus]++
		row.PerRead[base] = append(row.PerRead[base], PerReadFeatures{
			Dist5p:  uint16(posInRead),
			Fraglen: uint16(len(qual)),
			Qual:    qual[posInRead],
			Strand:  strandByte,
		})
	}
}

// addUnstitchedSegment adds an unpaired portion of a single read to the
// pileup.
func (pm *pileupMutable) addUnstitchedSegment(read *readSNP, isMinus PosType, alignedBases []alignedPos, minBaseQual byte, perReadNeeded bool) {
	mask := pm.nCirc() - 1
	qual := read.samr.Qual
	if !perReadNeeded {
		for _, ab := range alignedBases {
			pm.addBase(ab.posInRef&mask, ab.posInRead, isMinus, read.seq8, qual, minBaseQual)
		}
	} else {
		strandByte := byte(pileup.GetStrand(read.samr))
		for _, ab := range alignedBases {
			pm.appendBase(ab.posInRef&mask, ab.posInRead, isMinus, read.seq8, qual, minBaseQual, strandByte)
		}
	}
}

// alignRelevantBases checks where the read intersects loaded BED intervals,
// and 'returns' a slice of (posInRef, posInRead) tuples with those locations.
//
// The function signature requires a preallocated []alignedPos since it's an
// inner-loop operation and its runtime goes up by ~45% when we return a
// freshly allocated []alignedPos instead.
func alignRelevantBases(result *[]alignedPos, read readSNP, bedPart *interval.BEDUnion) (err error) {
	*result = (*result)[:0]
	refID := read.samr.Ref.ID()
	posInRef := PosType(read.samr.Pos)
	endPos := read.mapEnd
	relevantIntervals := bedPart.OverlapByID(refID, posInRef, endPos)
	riIdx := 0
	// Note that we already filtered out reads which don't overlap at least one
	// BED interval, so bedIntervals is guaranteed to be nonempty.
	nextIntervalStart := relevantIntervals[riIdx]
	nextIntervalEnd := relevantIntervals[riIdx+1]
	cigar := read.samr.Cigar
	posInRead := PosType(0)
	for _, co := range cigar {
		// Iterate over one CIGAR operation at a time.
		cLen := PosType(co.Len())
		switch co.Type() {
		case sam.CigarMatch:
			nextPosInRef := posInRef + cLen
			if nextPosInRef > nextIntervalStart {
				// At least one interval overlaps the current CIGAR-match region.
				cigarMatchOffset := PosType(0)
				for {
					if posInRef+cigarMatchOffset < nextIntervalStart {
						cigarMatchOffset = nextIntervalStart - posInRef
					}
					if nextIntervalEnd > nextPosInRef {
						// Interval goes past the end of the current CIGAR-match region.
						for ; cigarMatchOffset < cLen; cigarMatchOffset++ {
							(*result) = append(*result, alignedPos{
								posInRef:  posInRef + cigarMatchOffset,
								posInRead: posInRead + cigarMatchOffset,
							})
						}
						break
					}
					offsetStop := nextIntervalEnd - posInRef
					for ; cigarMatchOffset < offsetStop; cigarMatchOffset++ {
						(*result) = append(*result, alignedPos{
							posInRef:  posInRef + cigarMatchOffset,
							posInRead: posInRead + cigarMatchOffset,
						})
					}
					riIdx += 2
					if riIdx == len(relevantIntervals) {
						// We've iterated over all overlapping intervals; safe to exit now.
						return
					}
					nextIntervalStart = relevantIntervals[riIdx]
					nextIntervalEnd = relevantIntervals[riIdx+1]
				}
			}
			posInRef = nextPosInRef
			posInRead += cLen
		case sam.CigarInsertion:
			// Insertions are currently ignored.
			posInRead += cLen
		case sam.CigarSkipped:
			// Same handling as deletion for now.
			fallthrough
		case sam.CigarDeletion:
			// We don't currently report overlapping deletions; we probably want to
			// in the future.
			posInRef += cLen
			// Whenever posInRef increases, we may move past some interval(s).
			for posInRef >= nextIntervalEnd {
				riIdx += 2
				if riIdx == len(relevantIntervals) {
					return
				}
				nextIntervalEnd = relevantIntervals[riIdx+1]
			}
			nextIntervalStart = relevantIntervals[riIdx]
		case sam.CigarSoftClipped:
			// Note that we can also handle soft- and hard-clips before the main
			// loop, taking advantage of the requirement that they have to be at the
			// beginning or end.
			posInRead += cLen
		case sam.CigarHardClipped:
			// do nothing
		default:
			// There are logical ways to handle the other codes, but we currently
			// don't expect them to show up at all.  If/when we start using them, we
			// can add proper handling of those cases.
			return fmt.Errorf("alignRelevantBases: unexpected CIGAR code %v", co)
		}
	}
	return
}

// addReadPair adds a single read or read-pair to the pileup.
func (pm *pileupMutable) addReadPair(reads []readSNP, isMinus PosType, pCtx *pileupContext) (err error) {
	for i, r := range reads {
		if err = alignRelevantBases(&(pm.alignedBaseBufs[i]), r, &pCtx.bedPart); err != nil {
			return
		}
		clipQuals(r.samr, pCtx.clip)
	}
	abb0 := pm.alignedBaseBufs[0]
	abb1 := pm.alignedBaseBufs[1]
	minBaseQual := pCtx.minBaseQual
	perReadNeeded := pCtx.perReadNeeded
	if (len(reads) == 1) || (len(abb1) == 0) {
		// Empty alignedBases is possible when the read has deletions overlapping
		// all SNP positions.
		if len(abb0) != 0 {
			pm.addUnstitchedSegment(&(reads[0]), isMinus, abb0, minBaseQual, perReadNeeded)
			curEndMax := abb0[len(abb0)-1].posInRef + 1
			if pm.endMax < curEndMax {
				pm.endMax = curEndMax
			}
		}
		return
	}
	if len(abb0) == 0 {
		pm.addUnstitchedSegment(&(reads[1]), isMinus, abb1, minBaseQual, perReadNeeded)
		curEndMax := abb1[len(abb1)-1].posInRef + 1
		if pm.endMax < curEndMax {
			pm.endMax = curEndMax
		}
		return
	}
	mask := pm.nCirc() - 1
	seq0 := reads[0].seq8
	seq1 := reads[1].seq8
	qual0 := reads[0].samr.Qual
	qual1 := reads[1].samr.Qual
	idx0 := 0
	idx1 := 0
	// Loop over all relevant positions, stitching shared bases when possible,
	// until we reach the end of at least one read.
	for (len(abb0) != idx0) && (len(abb1) != idx1) {
		posInRef0 := abb0[idx0].posInRef
		posInRef1 := abb1[idx1].posInRef
		if posInRef0 == posInRef1 {
			row := &pm.resultRingBuffer[posInRef0&mask]
			row.Depth++
			posInRead0 := abb0[idx0].posInRead
			curSeq0 := seq0[posInRead0] // 'Seq0' instead of 'Base0' since this still uses BAM encoding
			posInRead1 := abb1[idx1].posInRead
			curSeq1 := seq1[posInRead1]
			if curSeq0 == curSeq1 {
				base := pileup.Seq8ToEnumTable[curSeq0]
				if !perReadNeeded {
					if pCtx.qpt.lookup2(qual0[posInRead0], qual1[posInRead1]) || (base == pileup.BaseX) {
						row.Counts[base][isMinus]++
					}
				} else {
					// dist5p/fraglen are a bit complicated in this case.  Punt for now.
					panic("stitched per-read features not yet supported")
				}
			} else {
				row.Counts[pileup.BaseX][isMinus]++
			}
			idx0++
			idx1++
		} else if posInRef0 < posInRef1 {
			if !perReadNeeded {
				pm.addBase(posInRef0&mask, abb0[idx0].posInRead, isMinus, seq0, qual0, minBaseQual)
			} else {
				panic("stitched per-read features not yet supported")
			}
			idx0++
		} else {
			if !perReadNeeded {
				pm.addBase(posInRef1&mask, abb1[idx1].posInRead, isMinus, seq1, qual1, minBaseQual)
			} else {
				panic("stitched per-read features not yet supported")
			}
			idx1++
		}
	}
	var curEndMax PosType
	if len(abb0) != idx0 {
		pm.addUnstitchedSegment(&(reads[0]), isMinus, abb0[idx0:], minBaseQual, perReadNeeded)
		curEndMax = abb0[len(abb0)-1].posInRef + 1
	} else if len(abb1) != idx1 {
		pm.addUnstitchedSegment(&(reads[1]), isMinus, abb1[idx1:], minBaseQual, perReadNeeded)
		curEndMax = abb1[len(abb1)-1].posInRef + 1
	} else {
		// Note that this is guaranteed to be equal to
		// abb1[len(abb1)-1].posInRef + 1, since we can only get here if the last
		// base was stitched, meaning abb0[idx0].posInRef was equal to
		// abb1[idx1].posInRef.
		curEndMax = abb0[len(abb0)-1].posInRef + 1
	}
	if pm.endMax < curEndMax {
		pm.endMax = curEndMax
	}
	return
}

// addOrphanReads adds every read in the firstread-table mapped to position <
// stopPos to the pileup, since their mates were filtered out or not found.
func (pm *pileupMutable) addOrphanReads(pCtx *pileupContext, stopPos PosType) (err error) {
	if !pCtx.stitch {
		return
	}
	nonempty := &pm.firstReads.nonempty
	mask := nonempty.NCirc() - 1
	var firstread [1]readSNP
	firstread[0].seq8 = pm.seq8Buf
	ignoreStrand := pCtx.ignoreStrand
	var isMinus PosType
	for pos := nonempty.FirstPos(); pos < stopPos; pos = nonempty.FirstPos() {
		circPos := pos & mask
		// & needed since this would otherwise be an array, not a slice
		bucketRow := &(pm.firstReads.buckets[circPos])
		for rs, bucketIdx := nonempty.NewRowScanner(); bucketIdx != -1; bucketIdx = rs.Next() {
			bucket := &(bucketRow[bucketIdx])
			for _, rec := range *bucket {
				convertSamr(&firstread[0], rec.samr)
				firstread[0].mapEnd = rec.mapEnd
				if !ignoreStrand {
					isMinus = PosType(pileup.GetStrand(rec.samr) - 1)
				}
				if err = pm.addReadPair(firstread[:1], isMinus, pCtx); err != nil {
					return
				}
				sam.PutInFreePool(rec.samr)
			}
			(*bucket) = (*bucket)[:0]
		}
	}
	return
}

// Key writing invariants:
// * writePos is 1 + <last 0-based position flushTo() was called with>
// * endMax is 1 + <last 0-based position with unwritten circular-buffer data>,
//   whenever there is unwritten data.
// * writePos <= firstReads.nonempty.FirstPos()

// flushToInternal flushes the contents of resultRingBuffer up to (and not
// including) 0-based position writeEnd, and then reinitializes the flushed
// rows.
// Note that this function assumes writePos is in [writeEnd - nCirc, writeEnd);
// the caller is expected to handle the other cases.
func (pm *pileupMutable) flushToInternal(rCtx *refContext, perReadNeeded bool, writeEnd PosType) (err error) {
	mask := pm.nCirc() - 1
	refID := rCtx.refID
	var start PosType
	var end PosType
	for pm.writePosScanner.Scan(&start, &end, writeEnd) {
		for pos := start; pos != end; pos++ {
			row := &pm.resultRingBuffer[pos&mask]
			if row.Depth == 0 {
				// It isn't strictly necessary to separate out this case, but it's a
				// significant performance win when zero-depth is common.
				pm.w.Append(&PileupRow{
					RefID: uint32(refID),
					Pos:   uint32(pos),
				})
			} else {
				fieldsPresent := uint32(FieldCounts)
				if !perReadNeeded {
					pm.w.Append(&PileupRow{
						FieldsPresent: fieldsPresent,
						RefID:         uint32(refID),
						Pos:           uint32(pos),
						Payload:       *row,
					})
				} else {
					fieldsPresent := FieldCounts
					// perRead contains regular slices instead of just arrays, so we need
					// to deep-copy it before clearing the ring-buffer copy.
					var perReadCopy [pileup.NBase][]PerReadFeatures
					for i := 0; i < pileup.NBase; i++ {
						if len(row.PerRead[i]) != 0 {
							fieldsPresent |= FieldPerReadA << uint(i)
							perReadCopy[i] = append([]PerReadFeatures(nil), row.PerRead[i]...)
						}
					}
					pm.w.Append(&PileupRow{
						FieldsPresent: uint32(fieldsPresent),
						RefID:         uint32(refID),
						Pos:           uint32(pos),
						Payload: PileupPayload{
							Depth:   row.Depth,
							Counts:  row.Counts,
							PerRead: perReadCopy,
						},
					})
					for i := range row.PerRead {
						row.PerRead[i] = row.PerRead[i][:0]
					}
				}
				for i := range row.Counts {
					for j := range row.Counts[i] {
						row.Counts[i][j] = 0
					}
				}
				row.Depth = 0
			}
		}
	}
	return
}

// flushTo is the main writer function.  flushEnd is
//   1 + <last position we now want to write results for>.
func (pm *pileupMutable) flushTo(rCtx *refContext, perReadNeeded bool, flushEnd PosType) (err error) {
	if pm.endMax > pm.writePosScanner.Pos() {
		writeEnd := minPosType(flushEnd, pm.endMax)
		err = pm.flushToInternal(rCtx, perReadNeeded, writeEnd)
		if (err != nil) || (pm.endMax >= flushEnd) {
			return
		}
	}
	return writeEmptyEntries(&pm.w, rCtx, flushEnd, &pm.writePosScanner)
}

type outputFormat int

const (
	formatBasestrandRio outputFormat = iota
	formatBasestrandTSV
	formatBasestrandTSVBgz
	formatTSV
	formatTSVBgz
)

type pileupSNPOpts struct {
	bedUnion         interval.BEDUnion
	clip             int
	colBitset        int
	fapath           string
	flagExclude      int
	format           outputFormat
	linearConsensus  int
	linearNosplit    bool
	mapq             int
	maxLinearBagSpan int
	maxReadLen       int
	maxReadSpan      int
	minBagDepth      int
	minBaseQual      int
	minBaseQualSum   int
	outPrefix        string
	padding          int
	parallelism      int
	provider         bamprovider.Provider
	refSeqs          [][]byte
	removeSq         bool
	shards           []gbam.Shard
	stitch           bool
	tempDir          string
}

func (pm *pileupMutable) finishRef(refIdxEnd int, opts *pileupSNPOpts, rCtx *refContext, pCtx *pileupContext) (err error) {
	if rCtx.refID != -1 {
		if err = pm.addOrphanReads(pCtx, PosTypeMax); err != nil {
			return
		}
		if err = pm.flushTo(rCtx, pCtx.perReadNeeded, PosTypeMax); err != nil {
			return
		}
	}
	// Suppose that the previous read was on chr5 and the current read is on
	// chr9.  This call writes the empty pileup entries for chr6..chr8.
	if err = pm.flushEmptyContigs(rCtx, opts.refSeqs, &pCtx.bedPart, pCtx.perReadNeeded, rCtx.refID+1, refIdxEnd); err != nil {
		return
	}
	return
}

func (pm *pileupMutable) nextRef(rCtx *refContext, newRefID int, opts *pileupSNPOpts, pCtx *pileupContext) (err error) {
	if err = pm.finishRef(newRefID, opts, rCtx, pCtx); err != nil {
		return
	}
	pm.endMax = 0
	endpoints := pCtx.bedPart.EndpointsByID(newRefID)
	pm.writePosScanner = interval.NewUnionScanner(endpoints)
	rCtx.refID = newRefID
	rCtx.refName = pCtx.bedPart.RefNames[newRefID] // only needed for error messages
	return
}

type pileupShardContext struct {
	strandReq    pileup.StrandType
	prevLimitID  int
	prevLimitPos int
	shardOverlap bool
	readPair     [2]readSNP
}

func (pm *pileupMutable) processShard(shard gbam.Shard, opts *pileupSNPOpts, rCtx *refContext, pCtx *pileupContext, psCtx *pileupShardContext) (err error) {
	iter := opts.provider.NewIterator(shard)
	defer func() {
		if e := iter.Close(); e != nil && err == nil {
			err = e
		}
	}()

	ignoreStrand := pCtx.ignoreStrand
	var isMinus PosType
	for iter.Scan() {
		curRead := iter.Record()
		if psCtx.shardOverlap {
			// The first few reads may have already been processed while iterating
			// over the previous shard, since padding > 0.  Don't reprocess them.
			notDone := true
			for (curRead.Ref.ID() == psCtx.prevLimitID) && (curRead.Pos < psCtx.prevLimitPos) {
				sam.PutInFreePool(curRead)
				notDone = iter.Scan()
				if !notDone {
					break
				}
				curRead = iter.Record()
			}
			if !notDone {
				break
			}
			psCtx.shardOverlap = false
		}
		// -flag-exclude, -mapq, and blank-read filters
		if (opts.flagExclude&int(curRead.Flags) != 0) || (opts.mapq > int(curRead.MapQ)) || (len(curRead.Cigar) == 0) {
			sam.PutInFreePool(curRead)
			continue
		}
		// -remove-sq filter
		if opts.removeSq {
			var libraryBagSize int
			if libraryBagSize, err = curRead.LibraryBagSize(); err != nil {
				return
			}
			if libraryBagSize < 2 {
				sam.PutInFreePool(curRead)
				continue
			}
		}
		// -min-bag-depth filter
		if opts.minBagDepth != 0 {
			var bagDepthFilterFail bool
			if bagDepthFilterFail, err = bagDepthFilter(curRead, opts.minBagDepth); err != nil {
				return
			}
			if bagDepthFilterFail {
				sam.PutInFreePool(curRead)
				continue
			}
		}
		strand := pileup.GetStrand(curRead)
		if psCtx.strandReq != pileup.StrandNone {
			// -per-strand filter
			if strand != psCtx.strandReq {
				sam.PutInFreePool(curRead)
				continue
			}
		} else if (strand == pileup.StrandNone) && (!ignoreStrand) {
			// We also don't need to include nonstandard-strand reads in the pileup
			// when we're only reporting (base x strand) counts.
			sam.PutInFreePool(curRead)
			continue
		}
		// Okay, this read might actually matter.
		//
		// 1. Note this read's start position, and flush as many previous positions
		//    as possible.
		if rCtx.refID != curRead.Ref.ID() {
			if err = pm.nextRef(rCtx, curRead.Ref.ID(), opts, pCtx); err != nil {
				return
			}
		}
		flushEnd := PosType(curRead.Pos)
		if opts.stitch {
			firstOrphanPos := pm.firstReads.nonempty.FirstPos()
			if firstOrphanPos != circular.FirstPosEmpty {
				// If we are stitching, and the firstread-table isn't empty, we can
				// only safely flush to
				//   max(firstOrphanPos, curRead.Pos - opts.maxReadSpan),
				// instead of curRead.Pos.
				if firstOrphanPos+PosType(opts.maxReadSpan) < flushEnd {
					flushEnd -= PosType(opts.maxReadSpan)
				} else {
					flushEnd = firstOrphanPos
				}
			}
		}
		if pm.writePosScanner.Pos() < flushEnd {
			if err = pm.addOrphanReads(pCtx, flushEnd); err != nil {
				return
			}
			if err = pm.flushTo(rCtx, pCtx.perReadNeeded, flushEnd); err != nil {
				return
			}
		}

		// 2. Note read-length and read-span, verify that they don't violate our
		//    ring-buffer sizing assumptions, and check for overlap with a BED
		//    interval.  If there's no overlap, skip this read.
		if len(curRead.Qual) > opts.maxReadLen {
			return fmt.Errorf("pileupMutable.processShard: maxReadLen is %d, but read %s at %s:%d has length %d", opts.maxReadLen, curRead.Name, rCtx.refName, curRead.Pos, len(curRead.Qual))
		}
		span, _ := curRead.Cigar.Lengths()
		if span > opts.maxReadSpan {
			return fmt.Errorf("pileupMutable.processShard: maxReadSpan is %d, but read %s at %s:%d has span %d", opts.maxReadSpan, curRead.Name, rCtx.refName, curRead.Pos, span)
		}
		mapEnd := PosType(curRead.Pos + span)
		if !pCtx.bedPart.IntersectsByID(rCtx.refID, PosType(curRead.Pos), mapEnd) {
			sam.PutInFreePool(curRead)
			continue
		}
		psCtx.readPair[0].mapEnd = mapEnd

		// 3. If stitching, look for a mate in the firstread-table.
		//    - If it's there, set nRead to 2 and save both ends to psCtx.readPair.
		//    - If it's known to be missing, set nRead to 1 and save just the
		//      current read to psCtx.readPair.
		//    - If the mate might appear later, save this read in the table for
		//      later processing, set nRead to 0, and move on to the next read in
		//      the BAM/PAM.
		nRead := 1
		if opts.stitch {
			nRead = pm.firstReads.addOrRemove(&psCtx.readPair, curRead, strand, opts.maxReadSpan)
			if nRead == 0 {
				continue
			}
		} else {
			convertSamr(&(psCtx.readPair[0]), curRead)
		}

		// 4. Add this read or read-pair to the pileup, and then recycle the memory
		//    it was using.
		if !ignoreStrand {
			isMinus = PosType(strand - 1)
		}
		if err = pm.addReadPair(psCtx.readPair[:nRead], isMinus, pCtx); err != nil {
			return
		}
		for i := 0; i != nRead; i++ {
			sam.PutInFreePool(psCtx.readPair[i].samr)
		}
	}
	return
}

func pileupSNPMain(ctx context.Context, opts *pileupSNPOpts, strandReq pileup.StrandType) (err error) {
	if (opts.format != formatTSV) && (opts.format != formatTSVBgz) && (strandReq != pileup.StrandNone) {
		err = fmt.Errorf("pileupSNPMain: single-strand mode not supported with basestrand output (strands are already tracked separately)")
		return
	}
	nShard := len(opts.shards)
	parallelism := minInt(opts.parallelism, nShard)

	// When we aren't stitching, it is always safe to flush final pileup results
	// for all positions before the current read-start; we only need to keep
	// track of the maxReadSpan positions past that point.
	// However, when stitching, we also need to track the maxReadSpan positions
	// before the current read-start, since there may be a read-pair where the
	// first read starts almost that far behind us, we haven't encountered the
	// second read yet, and there's some overlap between the two read-sides.
	nCirc := PosType(circular.NextExp2(opts.maxReadSpan))
	if opts.stitch {
		nCirc = nCirc * 2
	}

	var qpt qualPassTable
	if qpt, err = newQualPassTable(byte(opts.minBaseQual)); err != nil {
		return
	}

	if opts.tempDir != "" {
		// Note that we don't actually use the temp directory when parallelism ==
		// 1.  But may as well still force it to exist for consistency.
		if err = os.MkdirAll(opts.tempDir, 0644); err != nil {
			return
		}
	}

	tmpFiles := make([]*os.File, parallelism)
	defer func() {
		for _, f := range tmpFiles {
			if f != nil {
				if e := f.Close(); e != nil && err == nil {
					err = e
				}
			}
		}
	}()
	for jobIdx := range tmpFiles {
		if tmpFiles[jobIdx], err = ioutil.TempFile(opts.tempDir, "pileup_tmp"+strconv.Itoa(jobIdx)+"_*.rio"); err != nil {
			return
		}
	}

	log.Printf("pileupSNPMain: starting main loop (%d jobs)\n", parallelism)
	err = traverse.Each(parallelism, func(jobIdx int) error {
		startIdx := (jobIdx * nShard) / parallelism
		endIdx := ((jobIdx + 1) * nShard) / parallelism
		shardSlice := opts.shards[startIdx:endIdx]

		rCtx := refContext{
			refID: -1,
		}
		maxReadLen := opts.maxReadLen
		results := newPileupMutable(nCirc, maxReadLen, opts.stitch, tmpFiles[jobIdx])

		// We already got the header before, so it shouldn't be possible for this
		// call to generate a new error.
		header, _ := opts.provider.GetHeader()
		headerRefs := header.Refs()
		padding := PosType(opts.padding)

		// This contains information needed by some functions called by
		// pileupMutable.processShard.
		// Probable todo: set ignoreStrand to false when per-strand TSV output is
		// requested, and move the per-strand-reporting logic to the final
		// internal-format -> final format conversion; this lets up stop executing
		// the main loop twice.
		pCtx := pileupContext{
			clip:          opts.clip,
			ignoreStrand:  (opts.format == formatTSV) || (opts.format == formatTSVBgz),
			perReadNeeded: ((opts.colBitset & (colBitEndDists | colBitQuals | colBitFraglens | colBitStrands)) != 0),
			minBaseQual:   byte(opts.minBaseQual),
			stitch:        opts.stitch,
			qpt:           &qpt,
		}
		// The final concatenation step does not currently deduplicate records in
		// the overlapping region, so it's necessary to precisely split the
		// BEDUnion here.
		{
			firstCoordRange := gbam.ShardToCoordRange(shardSlice[0])
			startRefID := int(firstCoordRange.Start.RefId)
			startPos := PosType(firstCoordRange.Start.Pos)

			lastCoordRange := gbam.ShardToCoordRange(shardSlice[len(shardSlice)-1])
			limitRefID := int(lastCoordRange.Limit.RefId)
			limitPos := PosType(lastCoordRange.Limit.Pos)
			if limitRefID < 0 {
				limitRefID = len(headerRefs) - 1
				limitPos = PosType(headerRefs[limitRefID].Len())
			}
			pCtx.bedPart = opts.bedUnion.Subset(startRefID, startPos, limitRefID, limitPos)
		}

		// This contains context only needed by the top-level processShard
		// function.
		psCtx := pileupShardContext{
			strandReq:   strandReq,
			prevLimitID: -1,
		}
		psCtx.readPair[0].seq8 = make([]byte, 0, maxReadLen)
		psCtx.readPair[1].seq8 = make([]byte, 0, maxReadLen)

		for _, shard := range shardSlice {
			// May as well skip completely-nonoverlapping shards.
			if intersectionIsEmpty(&shard, headerRefs, &pCtx.bedPart) {
				continue
			}
			if e := results.processShard(shard, opts, &rCtx, &pCtx, &psCtx); e != nil {
				return e
			}
			psCtx.shardOverlap = true
			coordRange := gbam.ShardToCoordRange(shard)
			psCtx.prevLimitID = int(coordRange.Limit.RefId)
			psCtx.prevLimitPos = int(coordRange.Limit.Pos) + int(padding)
		}
		// Flush last entries, unless there were no entries at all.
		if e := results.finishRef(len(headerRefs), opts, &rCtx, &pCtx); e != nil {
			return e
		}

		return results.w.Finish()
	})
	if err != nil {
		return
	}
	log.Printf("pileupSNPMain: main loop complete")
	mainPath := opts.outPrefix
	if strandReq == pileup.StrandFwd {
		mainPath = mainPath + ".strand.fwd"
	} else if strandReq == pileup.StrandRev {
		mainPath = mainPath + ".strand.rev"
	}
	header, _ := opts.provider.GetHeader()
	var refNames []string
	for _, ref := range header.Refs() {
		refNames = append(refNames, ref.Name())
	}
	switch opts.format {
	case formatTSV:
		err = ConvertPileupRowsToTSV(ctx, tmpFiles, mainPath, opts.colBitset, false, opts.parallelism, refNames, opts.refSeqs)
	case formatTSVBgz:
		err = ConvertPileupRowsToTSV(ctx, tmpFiles, mainPath, opts.colBitset, true, opts.parallelism, refNames, opts.refSeqs)
	case formatBasestrandRio:
		err = ConvertPileupRowsToBasestrandRio(ctx, tmpFiles, mainPath, refNames)
	case formatBasestrandTSV:
		err = ConvertPileupRowsToBasestrandTSV(ctx, tmpFiles, mainPath, opts.colBitset, false, opts.parallelism, refNames, opts.refSeqs)
	case formatBasestrandTSVBgz:
		err = ConvertPileupRowsToBasestrandTSV(ctx, tmpFiles, mainPath, opts.colBitset, true, opts.parallelism, refNames, opts.refSeqs)
	}
	return
}

func Pileup(ctx context.Context, xampath, fapath, format, outPrefix string, rawOpts *Opts, refSeqs [][]byte) (err error) {
	// 1. Parse and validate command-line parameters
	// 2. Read .bam header, BED, .fa
	// 3. Construct disjoint shards with necessary padding
	var opts pileupSNPOpts
	opts.clip = rawOpts.Clip
	opts.maxReadLen = rawOpts.MaxReadLen
	if (opts.clip < 0) || (opts.clip*2 >= opts.maxReadLen) {
		return fmt.Errorf("Pileup: invalid clip= argument")
	}

	opts.fapath = fapath
	opts.flagExclude = rawOpts.FlagExclude
	opts.mapq = rawOpts.Mapq

	opts.maxReadSpan = rawOpts.MaxReadSpan
	if opts.maxReadLen > opts.maxReadSpan {
		return fmt.Errorf("Pileup: max-read-len= argument cannot be larger than max-read-span= argument")
	}

	opts.minBagDepth = rawOpts.MinBagDepth
	opts.minBaseQual = rawOpts.MinBaseQual
	opts.outPrefix = outPrefix

	opts.parallelism = rawOpts.Parallelism
	if opts.parallelism <= 0 {
		opts.parallelism = runtime.NumCPU()
	}

	opts.removeSq = rawOpts.RemoveSq
	opts.tempDir = rawOpts.TempDir
	// can look up a map instead
	if format == "basestrand-rio" {
		opts.format = formatBasestrandRio
	} else if format == "basestrand-tsv" {
		opts.format = formatBasestrandTSV
	} else if format == "basestrand-tsv-bgz" {
		opts.format = formatBasestrandTSVBgz
	} else if format == "tsv" {
		opts.format = formatTSV
	} else if format == "tsv-bgz" {
		opts.format = formatTSVBgz
	} else {
		return fmt.Errorf("Pileup: unrecognized format= argument")
	}
	colBitsetDefault := colBitDpRef | colBitHighQ | colBitLowQ
	if rawOpts.Cols != "" {
		if opts.format == formatBasestrandRio {
			return fmt.Errorf("Pileup: -cols cannot be used with basestrand-rio output")
		}
		if opts.colBitset, err = pileup.ParseCols(rawOpts.Cols, colNameMap, colBitsetDefault); err != nil {
			return err
		}
	} else {
		opts.colBitset = colBitsetDefault
	}

	dropFields := []gbam.FieldType{
		gbam.FieldTempLen,
	}
	if opts.minBagDepth == 0 {
		dropFields = append(dropFields, gbam.FieldAux)
	}
	opts.provider = bamprovider.NewProvider(xampath, bamprovider.ProviderOpts{
		Index:      rawOpts.BamIndexPath,
		DropFields: dropFields})
	defer func() {
		if e := opts.provider.Close(); e != nil && err == nil {
			err = e
		}
	}()

	var header *sam.Header
	var regionEntry interval.Entry
	if (rawOpts.Region != "") || (rawOpts.BedPath != "") {
		if header, err = opts.provider.GetHeader(); err != nil {
			return
		}
		if rawOpts.Region != "" {
			if regionEntry, err = interval.ParseRegionString(rawOpts.Region); err != nil {
				return
			}
		}
		if rawOpts.BedPath != "" {
			if rawOpts.Region != "" {
				// TODO: extend one of these functions, or the BEDUnion type, to
				// perform the necessary intersection operation.
				return fmt.Errorf("Pileup: -region and -bed flags can't be used together yet")
			}
			if opts.bedUnion, err = interval.NewBEDUnionFromPath(rawOpts.BedPath, interval.NewBEDOpts{SAMHeader: header}); err != nil {
				return
			}
		} else if rawOpts.Region != "" {
			if opts.bedUnion, err = interval.NewBEDUnionFromEntries([]interval.Entry{regionEntry}, interval.NewBEDOpts{SAMHeader: header}); err != nil {
				return
			}
		}
	} else {
		return fmt.Errorf("Pileup: either -bed or -region is currently required")
	}
	headerRefs := header.Refs()

	// padding requirement increases if we need to keep track of fragment lengths
	opts.padding = rawOpts.MaxReadSpan

	if regionEntry.RefName == "" {
		// Easiest to join the results at the end if we have each job process huge
		// disjoint (up to padding) chunks of the genome, so let's start with that
		// strategy.  Can experiment with finer-grained parallelism later.
		if opts.shards, err = opts.provider.GenerateShards(bamprovider.GenerateShardsOpts{
			Padding:   opts.padding,
			NumShards: opts.parallelism,
		}); err != nil {
			return
		}
	} else {
		// todo: add automatic sub-sharder to bamprovider
		found := false
		for _, ref := range headerRefs {
			refName := ref.Name()
			if refName == regionEntry.RefName {
				shard := gbam.Shard{
					StartRef: ref,
					EndRef:   ref,
					Start:    int(regionEntry.Start0),
					End:      int(regionEntry.End),
				}
				opts.shards = []gbam.Shard{shard}
				found = true
				break
			}
		}
		if !found {
			return fmt.Errorf("Pileup: region= contig not in BAM/PAM")
		}
	}

	if refSeqs == nil {
		if opts.refSeqs, err = pileup.LoadFa(ctx, fapath, 250000000, headerRefs); err != nil {
			return
		}
	} else {
		opts.refSeqs = refSeqs
	}

	opts.stitch = rawOpts.Stitch

	if rawOpts.PerStrand {
		// special case: run twice, filtering on different strand each time
		if err = pileupSNPMain(ctx, &opts, pileup.StrandFwd); err != nil {
			return
		}
		if err = pileupSNPMain(ctx, &opts, pileup.StrandRev); err != nil {
			return
		}
	} else {
		err = pileupSNPMain(ctx, &opts, pileup.StrandNone)
	}
	return
}
