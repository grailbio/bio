// +build !cgo

package bgzf

import (
	"io"
)

// NewWriterParams fails when compiled without cgo.
func NewWriterParams(w io.Writer, level, uncompressedBlockSize, gzipStrategy, gzipXFL int) (*Writer, error) {
	panic("NewWriterParams requires cgo")
	return nil, nil
}
