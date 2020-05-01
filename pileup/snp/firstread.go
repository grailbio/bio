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
	gunsafe "github.com/grailbio/base/unsafe"
	"github.com/grailbio/bio/biosimd"
	"github.com/grailbio/bio/circular"
	gbam "github.com/grailbio/bio/encoding/bam"
	"github.com/grailbio/bio/pileup"
	"github.com/grailbio/hts/sam"
	"github.com/spaolacci/murmur3"
)

// We keep track of read-pairs and, when the two ends overlap, we treat the
// overlapped region as having higher quality (phred-product of individual base
// qualities).
// bamprovider.PairIterator does not provide the ordering guarantee needed by
// our circular-buffer computation strategy, so we use the ordinary iterator
// and manage read-pair lookup manually.  The -max-read-span= condition keeps
// things simpler.
//
// We don't explicitly store refID or pos since those are assumed to be known
// from other context.

// readSNP is the form of the read used for inner-loop computation.
type readSNP struct {
	// seq8 is the unpacked (1 byte per base) form of samr.Seq.
	seq8 []byte
	// samr is the original *sam.Record.  We retain a copy, and reference most
	// fields through it, since we want to call sam.PutInFreePool() on it later
	// anyway.
	samr *sam.Record
	// mapEnd is the read's end position.  We need to compute it during initial
	// read processing to decide whether/when to handle the read, and it is a
	// little bit expensive to compute so we keep it around.
	mapEnd PosType
}

// augSamr is a slightly-augmented *sam.Record.
type augSamr struct {
	// samr stores the current read.
	samr *sam.Record
	// mapEnd caches the read end position.
	mapEnd PosType
}

// convertSamr extracts the pileup-relevant fields of a *sam.Record.
func convertSamr(read *readSNP, samr *sam.Record) {
	lSeq := len(samr.Qual)
	gunsafe.ExtendBytes(&read.seq8, lSeq)
	if lSeq != 0 {
		biosimd.UnpackSeq(read.seq8, gbam.UnsafeDoubletsToBytes(samr.Seq.Seq))
	}
	read.samr = samr
}

// Size of readName hash tables.  Currently must be a power of 2, and less than
// 256 * BitsPerWord.
const readNameHtableSize = 1024

// This may be implemented as a linked list in more complex pileups, depending
// on allocation details.
type firstreadSNPBucket []augSamr

// firstreadSNPTable is a chained hash table which tracks incomplete read pairs
// relevant to the pileup computation.
type firstreadSNPTable struct {
	// buckets contains the hash table buckets.  It's two-dimensional: major
	// dimension is (pos % nCirc), minor dimension is name-hash.
	buckets [][readNameHtableSize]firstreadSNPBucket

	// nonempty tracks which buckets are nonempty.
	nonempty circular.Bitmap
}

// nCirc() returns the common size of the position-based circular buffers.  It
// is guaranteed to be a power of 2 (so that (pos % nCirc) can be computed via
// bitwise mask instead of a far slower generic integer modulus), and large
// enough to fit the pileup's "active interval".
func (frt *firstreadSNPTable) nCirc() PosType {
	return frt.nonempty.NCirc()
}

// firstreadSNPTableHashName maps a readname to a firstreadSNPTable bucket
// index.  This isn't a []firstreadSNPBucket pointer receiver since the raw
// number is also needed to update rows[].nz.
func firstreadSNPTableHashName(name string) uint32 {
	// Don't need this to be cryptographically secure.  Could try other hash
	// functions later.
	return murmur3.Sum32(gunsafe.StringToBytes(name)) % readNameHtableSize
}

func (frt *firstreadSNPTable) add(rec augSamr) {
	pos := PosType(rec.samr.Pos)
	nCirc := frt.nCirc()
	circPos := pos & (nCirc - 1)
	hashrem := firstreadSNPTableHashName(rec.samr.Name)
	bucket := &(frt.buckets[circPos][hashrem])
	if len(*bucket) == 0 {
		frt.nonempty.Set(pos, circPos, hashrem)
	}
	*bucket = append(*bucket, rec)
}

// Third parameter must change in linked list case.
func (frt *firstreadSNPTable) remove(pos PosType, hashrem uint32, i int) {
	nCirc := frt.nCirc()
	circPos := pos & (nCirc - 1)
	bucket := &(frt.buckets[circPos][hashrem])
	lenMinus1 := len(*bucket) - 1
	if i != lenMinus1 {
		(*bucket)[i] = (*bucket)[lenMinus1]
	}
	(*bucket) = (*bucket)[:lenMinus1]
	if lenMinus1 == 0 {
		frt.nonempty.Clear(pos, circPos, hashrem)
	}
}

// tryRemove tries to remove the read with the given readName, position, and
// MatePos from the table.  Returns (*sam.Record, mapEnd) on success, nil on
// failure.
func (frt *firstreadSNPTable) tryRemove(samr *sam.Record) (*sam.Record, PosType) {
	expectedMatePos := samr.Pos
	pos := PosType(samr.MatePos)
	nCirc := frt.nCirc()
	circPos := pos & (nCirc - 1)

	readName := samr.Name
	hashrem := firstreadSNPTableHashName(readName)
	bucket := &(frt.buckets[circPos][hashrem])
	for i, rec := range *bucket {
		if (rec.samr.MatePos == expectedMatePos) && (rec.samr.Name == readName) {
			// Found the read; now remove it from the table.
			frt.remove(pos, hashrem, i)
			return rec.samr, rec.mapEnd
		}
	}
	return nil, 0
}

// addOrRemove is the central function for iterating over a mix of reads and
// read-pairs in order.  More precisely:
// - If the read-pair does not have standard orientation (i.e. strand ==
//   pileup.StrandNone), or the two ends trivially don't overlap, process each
//   end separately.  readPair[0] is filled, and the return value is 1,
//   indicating that readPair[:1] should be processed now.
//   Note that, in other pileup commands, a different positional condition will
//   be checked: some non-overlapping read-pairs will still require joint
//   processing because the ends might still be close enough for fragment
//   length to be estimated.
// - Otherwise, if Pos < MatePos, store the read in the firstreadSNPTable for
//   delayed processing.  The return value is 0.
// - Otherwise, look for the mate in the firstreadSNPTable.  If the mate is
//   there, fill readPair[0] with the starting read and readPair[1] with the
//   mate, and return 2.  If it isn't in the table (due to BAM/PAM filtering,
//   or e.g. the other end having terrible MAPQ), usually fill readPair[0] and
//   return 1.  The exception is if Pos and MatePos are identical, in which
//   case the mate might still show up later; then we store the read in the
//   firstreadSNPTable and return 0.
// readPair[0].mapEnd is currently expected to be initialized to Pos + [first
// return value of samr.Cigar.Lengths()].  This may be changed to an additional
// function parameter later.
func (frt *firstreadSNPTable) addOrRemove(readPair *[2]readSNP, samr *sam.Record, strand pileup.StrandType, maxReadSpan int) int {
	matePos := PosType(samr.MatePos)
	pos := PosType(samr.Pos)
	mapEnd := readPair[0].mapEnd
	if (strand == pileup.StrandNone) || (matePos >= mapEnd) || (matePos+PosType(maxReadSpan) <= pos) {
		convertSamr(&(readPair[0]), samr)
		return 1
	}
	if pos >= matePos {
		mateSamr, mateMapEnd := frt.tryRemove(samr)
		if mateSamr != nil {
			convertSamr(&(readPair[0]), samr)
			readPair[1].mapEnd = mateMapEnd
			convertSamr(&(readPair[1]), mateSamr)
			return 2
		}
		if pos != matePos {
			// Missing mate.  Process this read on its own.
			convertSamr(&(readPair[0]), samr)
			return 1
		}
		// If we get here, pos == matePos and we haven't seen the other end.
	}
	frt.add(augSamr{
		samr:   samr,
		mapEnd: mapEnd,
	})
	return 0
}

func newFirstreadSNPTable(nCirc PosType) firstreadSNPTable {
	return firstreadSNPTable{
		buckets:  make([][readNameHtableSize]firstreadSNPBucket, nCirc),
		nonempty: circular.NewBitmap(nCirc, 1+((readNameHtableSize-1)/circular.BitsPerWord)),
	}
}
