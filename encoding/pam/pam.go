// Copyright 2018 GRAIL, Inc. All rights reserved.
// Use of this source code is governed by the Apache 2.0
// license that can be found in the LICENSE file.

// Package pam implements PAM reader and writer. PAM is a more compact and
// faster alternative to BAM.
//
// Most people, however, will want to use the bamprovider
// (https://godoc.org/github.com/grailbio/bio/encoding/bamprovider) read PAM
// instead.  The bamprovider works for both BAM and PAM files transparently.
//
// REAMDE.md (https://github.com/grailbio/bio/encoding/pam/README.md) contains
// More detailed information the PAM file format.
package pam

import (
	"fmt"
	"io"
	"regexp"
	"sort"
	"strconv"

	"github.com/grailbio/base/file"
	"github.com/grailbio/base/vcontext"
	"github.com/grailbio/bio/biopb"
	gbam "github.com/grailbio/bio/encoding/bam"
	"v.io/x/lib/vlog"
)

func validateRecAddr(r biopb.Coord) error {
	if r.RefId < -1 || r.Pos < 0 || r.Seq < 0 {
		return fmt.Errorf("Invalid record addr: %+v", r)
	}
	return nil
}

// CoordRangePathString returns a string that can be used as part of a pathname.
func CoordRangePathString(r biopb.CoordRange) string {
	return fmt.Sprintf("%s,%s", CoordPathString(r.Start), CoordPathString(r.Limit))
}

// ValidateRecRange validates "r" and normalize its fields, if necessary. In
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
		return fmt.Errorf("LimitRef (%+v) <= StartRef (%+v)", r.Limit, r.Start)
	}
	return nil
}

// ReadSeekCloser is a combination of io.ReadSeeker and io.Closer.
//
// TODO(saito) This should be moved to somewhere more generic.
type ReadSeekCloser interface {
	io.ReadSeeker
	io.Closer
}

// PathString generates a string that can be used to embed in a pathname.  Use
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

// FieldDataPath returns the path of the file storing data for the given record
// range.
func FieldDataPath(dir string, recRange biopb.CoordRange, f gbam.FieldType) string {
	return fmt.Sprintf("%s/%s.%v", dir, CoordRangePathString(recRange), f)
}

// ShardIndexPath returns the path of shard index file.
func ShardIndexPath(dir string, recRange biopb.CoordRange) string {
	return fmt.Sprintf("%s/%s.index", dir, CoordRangePathString(recRange))
}

// DefaultVersion is the string embedded in ShardIndex.version.
const DefaultVersion = "PAM2"

// ShardIndexMagic is the value of ShardIndex.Magic.
const ShardIndexMagic = uint64(0x725c7226be794c60)

// FieldIndexMagic is the value of FieldIndex.Magic.
const FieldIndexMagic = uint64(0xe360ac9026052aca)

// FileType defines the type of the file, either data or index.
type FileType int

const (
	// FileTypeUnknown is a sentinel
	FileTypeUnknown FileType = iota
	// FileTypeShardIndex represents a *.index file
	FileTypeShardIndex
	// FileTypeFieldData represents a *.<fieldname> file
	FileTypeFieldData
)

// FileInfo is the result of parsing a pathname.
type FileInfo struct {
	// Path is the value passed to ParsePath.
	Path string
	Type FileType

	// Field is the field stored in the file. Meaningful iff Type ==
	// FileTypeFieldData.
	Field gbam.FieldType
	// Dir is the directory under which the file is stored.
	Dir string
	// Range is the record range that the file stores.
	Range biopb.CoordRange
}

var basenameRe = regexp.MustCompile("^(-|\\d+):(-|\\d+)(:\\d+)?,(-|\\d+):(-|\\d+)(:\\d+)?\\.(.+)$")

func parseExtension(str string) (FileType, gbam.FieldType, bool) {
	if str == "index" {
		return FileTypeShardIndex, gbam.FieldInvalid, true
	}
	fieldType, err := gbam.ParseFieldType(str)
	if err != nil {
		return FileTypeUnknown, gbam.FieldInvalid, false
	}
	return FileTypeFieldData, fieldType, true
}

func parseRecAddr(refidstr, posstr, seqstr string) (biopb.Coord, bool) {
	mustParseText := func(s string) int {
		v, err := strconv.Atoi(s)
		if err != nil {
			panic(err)
		}
		return v
	}
	addr := biopb.Coord{biopb.InfinityRefID, biopb.InfinityPos, 0}
	if refidstr != "-" {
		addr.RefId = int32(mustParseText(refidstr))
	}
	if posstr != "-" {
		addr.Pos = int32(mustParseText(posstr))
	}
	if seqstr != "" {
		addr.Seq = int32(mustParseText(seqstr[1:]))
	}
	return addr, true
}

// ParsePath parses a PAM path into constituent parts. For example,
// ParsePath("foo:0:1,3:4.index") will result in FileInfo{Path: "foo", Type:
// FileTypeIndex, Prefix: "foo", Range: {biopb.Coord{0,1,0}, biopb.Coord{3,4,0}}}.
func ParsePath(path string) (FileInfo, error) {
	fi := FileInfo{Path: path}
	basename := file.Base(path)
	m := basenameRe.FindStringSubmatch(basename)
	if m == nil {
		return fi, fmt.Errorf("%s: Unknown file type", path)
	}
	// Dir is the all but the last component of the path, plus the first part
	// of the basename.
	fi.Dir = file.Dir(path)
	var ok bool
	if fi.Type, fi.Field, ok = parseExtension(m[7]); !ok {
		return fi, fmt.Errorf("%s: Failed to parse suffix %v", path, m[7])
	}
	if fi.Range.Start, ok = parseRecAddr(m[1], m[2], m[3]); !ok {
		return fi, fmt.Errorf("%s: Invalid range start", path)
	}
	if fi.Range.Limit, ok = parseRecAddr(m[4], m[5], m[6]); !ok {
		return fi, fmt.Errorf("%s: Invalid range limit", path)
	}
	return fi, nil
}

// Remove deletes the files in the given PAM directory.  It returns an error if
// some of the existing files fails to delete.
func Remove(dir string) error {
	ctx := vcontext.Background()
	// TODO(saito) Provide equivalent of filepath.Join that works for URLs.
	lister := file.List(ctx, dir)
	n := 0
	for lister.Scan() {
		// TODO(saito) Use grailfile once it's ready.
		if err := file.Remove(ctx, lister.Path()); err != nil {
			return err
		}
		n++
	}
	if err := lister.Err(); err != nil {
		return err
	}
	if n == 0 {
		return nil
	}
	return file.Remove(ctx, dir)
}

// ListIndexes lists shard index files found for the given PAM files.  The
// returned list will be sorted by positions.
func ListIndexes(dir string) ([]FileInfo, error) {
	ctx := vcontext.Background()
	var infos []FileInfo

	lister := file.List(ctx, dir)
	for lister.Scan() {
		fi, err := ParsePath(lister.Path())
		if err != nil {
			vlog.Infof("Ignore file %v", err)
		}
		if fi.Type == FileTypeShardIndex {
			infos = append(infos, fi)
		}
	}
	if err := lister.Err(); err != nil {
		return nil, err
	}
	if len(infos) == 0 {
		return nil, fmt.Errorf("ListIndexes %v: no index files found", dir)
	}
	// TODO(saito) Check that ranges covers the universal range.
	sort.SliceStable(infos,
		func(i, j int) bool {
			return infos[i].Range.Start.LT(infos[j].Range.Start)
		})
	return infos, nil
}

func doassert(b bool) {
	if !b {
		panic(b)
	}
}
