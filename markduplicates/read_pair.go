package markduplicates

import (
	"fmt"

	"github.com/grailbio/base/log"
	"github.com/grailbio/bio/encoding/bam"
	"github.com/grailbio/hts/sam"
)

type readPair struct {
	// When a readPair contains only one read, that read always goes
	// to left, right will be nil, and rightFileIdx will be 0.  After
	// the mate arrives, we may move left to right and put the mate in
	// left, depending on the unclipped 5' position of each read.
	//
	// We will also use readPair to store mate-unmapped reads.  In
	// those cases, the read will be in left, and index in leftFileIdx.

	left  *sam.Record // The first read of a pair according to unclipped 5' position.
	right *sam.Record // The second read of a pair according to unclipped 5' position.

	// The index into the input file of each read.  If the readPair is
	// an entry in distantMates, then the fileIdx will temporarily be
	// set to the index into the shard, but getDistantMates replaces
	// those values with a global fileIdx before finishing.
	leftFileIdx  uint64
	rightFileIdx uint64
}

func (p *readPair) String() string {
	return fmt.Sprintf("(%s,%d,%d)(%s,%d,%d)", p.left.Ref.Name(), p.left.Pos, p.leftFileIdx,
		p.right.Ref.Name(), p.right.Pos, p.rightFileIdx)
}

func (p *readPair) addRead(newRead *sam.Record, fileIdx uint64) {
	// Complete the pair, and adjust left and right order if necessary.
	if p.right != nil {
		log.Fatalf("Tried to add third read %s %d to readPair", newRead.Name, newRead.Flags)
	}

	// Order left and right by:
	//  1) refId
	//  2) unclipped position
	//  3) fileIdx
	if newRead.Ref.ID() < p.left.Ref.ID() ||
		(newRead.Ref.ID() == p.left.Ref.ID() && bam.UnclippedFivePrimePosition(newRead) < bam.UnclippedFivePrimePosition(p.left)) ||
		(newRead.Ref.ID() == p.left.Ref.ID() && bam.UnclippedFivePrimePosition(newRead) == bam.UnclippedFivePrimePosition(p.left) &&
			fileIdx < p.leftFileIdx) {
		p.right = p.left
		p.rightFileIdx = p.leftFileIdx
		p.left = newRead
		p.leftFileIdx = fileIdx
	} else {
		p.right = newRead
		p.rightFileIdx = fileIdx
	}
}
