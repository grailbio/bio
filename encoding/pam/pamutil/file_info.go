package pamutil

import (
	"context"
	"fmt"
	"regexp"
	"sort"
	"strconv"

	"github.com/grailbio/base/file"
	"github.com/grailbio/base/log"
	"github.com/grailbio/bio/biopb"
)

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
//
// A PAM pathname looks like "dir/0:0,46:1653469.mapq" or "dir/0:0,46:1653469.index".
type FileInfo struct {
	// Path is the value passed to ParsePath.
	Path string

	// FileType is the type of the file. For "dir/0:0,46:1653469.mapq", the type
	// is FileTypeFieldData. For "dir/0:0,46:1653469.mapq", the type is
	// FileTypeFieldIndex.
	Type FileType

	// Field stores the field part of the filename. Field=="mapq" if the pathname
	// is "dir/0:0,46:1653469.mapq". It is meaningful iff Type ==
	// FileTypeFieldData.
	Field string

	// Dir is the directory under which the file is stored. Dir="dir" if the
	// pathname is "dir/0:0,46:1653469.mapq".
	Dir string
	// Range is the record range that the file stores. Range={Start:{0,0},
	// Limit:{46,1653469}} if the pathname is "dir/0:0,46:1653469.mapq".
	Range biopb.CoordRange
}

var basenameRe = regexp.MustCompile(`^(-|\d+):(-|\d+)(:\d+)?,(-|\d+):(-|\d+)(:\d+)?\.(.+)$`)

func parseExtension(str string) (FileType, string, bool) {
	if str == "index" {
		return FileTypeShardIndex, "", true
	}
	return FileTypeFieldData, str, true
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
		return fi, fmt.Errorf("parsepath %s: unknown file type", path)
	}
	// Dir is the all but the last component of the path, plus the first part
	// of the basename.
	fi.Dir = file.Dir(path)
	var ok bool
	if fi.Type, fi.Field, ok = parseExtension(m[7]); !ok {
		return fi, fmt.Errorf("parsepath %s: failed to parse suffix %v", path, m[7])
	}
	if fi.Range.Start, ok = parseRecAddr(m[1], m[2], m[3]); !ok {
		return fi, fmt.Errorf("parsepath %s: invalid range start", path)
	}
	if fi.Range.Limit, ok = parseRecAddr(m[4], m[5], m[6]); !ok {
		return fi, fmt.Errorf("parsepath %s: invalid range limit", path)
	}
	return fi, nil
}

// ListIndexes lists shard index files found for the given PAM files.  The
// returned list will be sorted by positions.
func ListIndexes(ctx context.Context, dir string) ([]FileInfo, error) {
	var infos []FileInfo

	lister := file.List(ctx, dir, true)
	for lister.Scan() {
		fi, err := ParsePath(lister.Path())
		if err != nil {
			log.Debug.Printf("Ignore file %v", err)
		}
		if fi.Type == FileTypeShardIndex {
			infos = append(infos, fi)
		}
	}
	if err := lister.Err(); err != nil {
		return nil, err
	}
	if len(infos) == 0 {
		return nil, fmt.Errorf("listindexes %s: no index files found", dir)
	}
	// TODO(saito) Check that ranges covers the universal range.
	sort.SliceStable(infos,
		func(i, j int) bool {
			return infos[i].Range.Start.LT(infos[j].Range.Start)
		})
	return infos, nil
}
