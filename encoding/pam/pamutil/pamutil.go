package pamutil

import (
	"fmt"

	"github.com/grailbio/base/backgroundcontext"
	"github.com/grailbio/base/file"
	"github.com/grailbio/bio/biopb"
	gbam "github.com/grailbio/bio/encoding/bam"
)

// DefaultVersion is the string embedded in ShardIndex.version.
const DefaultVersion = "PAM2"

// ShardIndexMagic is the value of ShardIndex.Magic.
const ShardIndexMagic = uint64(0x725c7226be794c60)

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
// range and the field.
func FieldDataPath(dir string, recRange biopb.CoordRange, field string) string {
	return fmt.Sprintf("%s/%s.%s", dir, CoordRangePathString(recRange), field)
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

func validateRecAddr(r biopb.Coord) error {
	if r.RefId < -1 || r.Pos < 0 || r.Seq < 0 {
		return fmt.Errorf("invalid record addr: %+v", r)
	}
	return nil
}

// ValidateCoordRange validates "r" and normalize its fields, if necessary. In
// particular, if the range fields are all zeros, the range is replaced by
// UniversalRange.
func ValidateCoordRange(r *biopb.CoordRange) error {
	// A Range where all values are zero is special-cased to mean "all rows"
	if r.Start.RefId == 0 && r.Start.Pos == 0 && r.Start.Seq == 0 &&
		r.Limit.RefId == 0 && r.Limit.Pos == 0 && r.Limit.Seq == 0 {
		*r = gbam.UniversalRange
		return nil
	}
	if err := validateRecAddr(r.Start); err != nil {
		return err
	}
	if err := validateRecAddr(r.Limit); err != nil {
		return err
	}
	if r.Limit.LE(r.Start) {
		return fmt.Errorf("limitref (%+v) <= startref (%+v)", r.Limit, r.Start)
	}
	return nil
}

// Remove deletes the files in the given PAM directory.  It returns an error if
// some of the existing files fails to delete.
func Remove(dir string) error {
	ctx := backgroundcontext.Get()
	err := file.RemoveAll(ctx, dir)
	file.Remove(ctx, dir) // nolint: errcheck
	return err
}
