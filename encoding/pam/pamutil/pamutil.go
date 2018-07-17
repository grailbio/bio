package pamutil

import (
	"fmt"

	"github.com/grailbio/bio/biopb"
	gbam "github.com/grailbio/bio/encoding/bam"
)

// DefaultVersion is the string embedded in ShardIndex.version.
const DefaultVersion = "PAM2"

// CoordPathString generates a string that can be used to embed in a pathname.  Use
// ParsePath() to parse such a string.
func CoordPathString(r biopb.Coord) string {
	var refStr, posStr string
	if r.RefId == biopb.InfinityRefID {
		refStr = "-"
	} else {
		refStr = fmt.Sprintf("%d", r.RefId)
	}
	if r.Pos == biopb.InfinityPos {
		posStr = "-"
	} else {
		posStr = fmt.Sprintf("%d", r.Pos)
	}
	if r.Seq == 0 {
		return fmt.Sprintf("%s:%s", refStr, posStr)
	}
	return fmt.Sprintf("%s:%s:%d", refStr, posStr, r.Seq)

}

// CoordRangePathString returns a string that can be used as part of a pathname.
func CoordRangePathString(r biopb.CoordRange) string {
	return fmt.Sprintf("%s,%s", CoordPathString(r.Start), CoordPathString(r.Limit))
}

// FieldDataPath returns the path of the file storing data for the given record
// range.
func FieldDataPath(dir string, recRange biopb.CoordRange, f gbam.FieldType) string {
	return fmt.Sprintf("%s/%s.%v", dir, CoordRangePathString(recRange), f)
}

// ShardIndexPath returns the path of shard index file.
func ShardIndexPath(dir string, recRange biopb.CoordRange) string {
	return fmt.Sprintf("%s/%s.index", dir, CoordRangePathString(recRange))
}

// BlockIntersectsRange checks if userRange and [startAddr, endAddr] intersect.
func BlockIntersectsRange(startAddr, endAddr biopb.Coord, userRange biopb.CoordRange) bool {
	// Note: We can't use biopb.CoordRange.Intersects here because
	// [b.StartAddr, b.EndAddr] is a closed section.
	if startAddr.LT(userRange.Limit) && userRange.Start.LE(endAddr) {
		return true
	}
	return false
}
