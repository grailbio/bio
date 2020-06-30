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
	"github.com/grailbio/base/recordio"
	gbam "github.com/grailbio/bio/encoding/bam"
	"github.com/grailbio/bio/interval"
	"github.com/grailbio/hts/sam"
)

// This contains SNP-pileup support functions and definitions which are
// peripheral enough that they'd probably decrease pileup.go's overall
// readability if placed there.

func minInt(a, b int) int {
	if a < b {
		return a
	}
	return b
}

func minPosType(a, b PosType) PosType {
	if a < b {
		return a
	}
	return b
}

// Returns true if the read does not pass the filter.
func bagDepthFilter(samr *sam.Record, minBagDepth int) (bool, error) {
	ds, e := samr.BagSize()
	return (ds < minBagDepth), e
}

const minQual = 2

// clipQuals sets the given number of Qual values on both ends of the read to
// minQual.
func clipQuals(samr *sam.Record, n int) {
	if n == 0 {
		return
	}
	qual := samr.Qual
	readLen := len(qual)
	if n*2 >= readLen {
		for i := 0; i < readLen; i++ {
			qual[i] = minQual
		}
	} else {
		for i := 0; i < n; i++ {
			qual[i] = minQual
		}
		for i := readLen - n; i < readLen; i++ {
			qual[i] = minQual
		}
	}
}

func intersectionIsEmpty(shardp *gbam.Shard, headerRefs []*sam.Reference, bedPart *interval.BEDUnion) bool {
	coordRange := gbam.ShardToCoordRange(*shardp)
	startRefID := int(coordRange.Start.RefId)
	// startRefID should always >= 0 since GenerateShards was not told to
	// include unmapped reads.
	startPos := PosType(coordRange.Start.Pos)
	limitRefID := int(coordRange.Limit.RefId)
	limitPos := PosType(coordRange.Limit.Pos)
	if limitRefID < 0 {
		// We need limitRefID >= startRefID below.
		limitRefID = len(headerRefs) - 1
		limitPos = PosType(headerRefs[limitRefID].Len())
	}
	return !bedPart.Intersects(startRefID, startPos, limitRefID, limitPos)
}

// writeEmptyEntries appends empty entries to the intermediate recordio file,
// up to flushEnd.
func writeEmptyEntries(w *recordio.Writer, rCtx *refContext, flushEnd PosType, writePosScanner *interval.UnionScanner) (err error) {
	refID := rCtx.refID
	var start PosType
	var end PosType
	for writePosScanner.Scan(&start, &end, flushEnd) {
		for pos := start; pos != end; pos++ {
			(*w).Append(&PileupRow{
				RefID: uint32(refID),
				Pos:   uint32(pos),
			})
		}
	}
	return
}

// nCirc() returns the common size of the position-based circular buffers.  It
// is guaranteed to be a power of 2 (so that (pos % nCirc) can be computed via
// bitwise mask instead of a far slower generic integer modulus), and large
// enough to fit the pileup's "active interval".
func (pm *pileupMutable) nCirc() PosType {
	return PosType(len(pm.resultRingBuffer))
}

func (pm *pileupMutable) flushEmptyContigs(rCtx *refContext, bedPart *interval.BEDUnion, perReadNeeded bool, refIdx, refIdxEnd int) (err error) {
	for ; refIdx < refIdxEnd; refIdx++ {
		pm.endMax = 0
		endpoints := bedPart.EndpointsByID(refIdx)
		pm.writePosScanner = interval.NewUnionScanner(endpoints)
		rCtx.refID = refIdx
		if err = pm.flushTo(rCtx, perReadNeeded, PosTypeMax); err != nil {
			return
		}
	}
	return
}
