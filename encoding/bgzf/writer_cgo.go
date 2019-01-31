// +build cgo

package bgzf

import (
	"fmt"
	"io"

	"github.com/yasushi-saito/zlibng"
)

// gzipCreate creates gzip WriteClosers.
type gzipFactory struct {
	level     int
	strategy  int
	gzWriter *zlibng.Writer
}

func (c *gzipFactory) create(w io.Writer) (io.WriteCloser, error) {
	var err error
	c.gzWriter, err = zlibng.NewWriter(w, zlibng.Opts{Level: c.level, Strategy: c.strategy})
	if err != nil {
		return nil, err
	}
	header := zlibng.GzipHeader{Extra: make([]byte, len(bgzfExtra))}
	copy(header.Extra[:], bgzfExtra[:])
	header.OS = 0xff // Unknown OS value
	if err = c.gzWriter.SetHeader(header); err != nil {
		c.gzWriter.Close() // nolint: errcheck
		return nil, err
	}
	return c.gzWriter, nil
}

// NewWriterParams returns a new .bgzf writer, with the given
// configuration parameters.  uncompressedBlockSize is the largest
// number of bytes to put into each .bgzf block.  gzipStrategy is a
// strategy value from gzip; possible values are DefaultStrategy,
// FilteredStrategy, HuffmanOnlyStrategy, RLEStrategy, and
// FixedStrategy.  gzipXFL will be written to the XFL gzip header
// field for each of the gzip blocks in the output; if gzipXFL is -1,
// then gzip with set XFL according to the other gzip configuration
// parameters.  Returns nil, error if there is a problem.
func NewWriterParams(w io.Writer, level, uncompressedBlockSize, gzipStrategy, gzipXFL int) (*Writer, error) {
	if uncompressedBlockSize > MaxUncompressedBlockSize {
		return nil, fmt.Errorf("uncompressedBlockSize %d is too large, max value is %d",
			uncompressedBlockSize, MaxUncompressedBlockSize)
	}
	if gzipXFL != -1 && (gzipXFL < 0 || gzipXFL > 255) {
		return nil, fmt.Errorf("gzipXFL must be -1 or in [0:255] not %d", gzipXFL)
	}

	return &Writer{
		factory:          &gzipFactory{level, gzipStrategy, nil},
		uncompressedSize: uncompressedBlockSize,
		xfl:              gzipXFL,
		w:                w,
	}, nil
}
