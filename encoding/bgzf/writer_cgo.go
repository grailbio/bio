// +build cgo

package bgzf

import (
	"fmt"
	"io"

	"github.com/youtube/vitess/go/cgzip"
)

// cgzipCreate creates cgzip WriteClosers.
type cgzipFactory struct {
	level     int
	strategy  int
	cgzWriter *cgzip.Writer
}

func (c *cgzipFactory) create(w io.Writer) (io.WriteCloser, error) {
	var err error
	c.cgzWriter, err = cgzip.NewWriterLevelFullConfig(w, c.level,
		cgzip.DEFAULT_COMPRESSED_BUFFER_SIZE,
		cgzip.DefaultWindow, cgzip.DefaultMemLevel,
		c.strategy,
	)
	if err != nil {
		return nil, err
	}
	c.cgzWriter.Header.Extra = make([]byte, len(bgzfExtra))
	copy(c.cgzWriter.Header.Extra[:], bgzfExtra[:])
	c.cgzWriter.Header.OS = 0xff // Unknown OS value

	return c.cgzWriter, nil
}

// NewWriterParams returns a new .bgzf writer, with the given
// configuration parameters.  uncompressedBlockSize is the largest
// number of bytes to put into each .bgzf block.  gzipStrategy is a
// strategy value from cgzip; possible values are DefaultStrategy,
// FilteredStrategy, HuffmanOnlyStrategy, RLEStrategy, and
// FixedStrategy.  gzipXFL will be written to the XFL gzip header
// field for each of the gzip blocks in the output; if gzipXFL is -1,
// then cgzip with set XFL according to the other cgzip configuration
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
		factory:          &cgzipFactory{level, gzipStrategy, nil},
		uncompressedSize: uncompressedBlockSize,
		xfl:              gzipXFL,
		w:                w,
	}, nil
}
