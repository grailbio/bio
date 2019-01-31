// Package bgzf includes a Writer for the .bgzf (block gzipped) file
// format.  A .bgzf file consists of one or more complete gzip blocks
// concatenated together.  Each of the gzip blocks must represent at
// most 64KB of uncompressed data, and the compressed size of the
// block must be at most 64KB.  The payload of the .bgzf file is equal
// to the uncompressed content of each block, concatenated together in
// order.  A valid .bgzf file ends with the 28 byte .bgzf terminator
// shown below; the terminator is a valid gzip block containing an
// empty payload.
//
// The .bgzf format is used by .bam files and Illumina .bcl.bgzf files
// from Nextseq instruments.
//
// For more information about the .bgzf file format, see the SAM/BAM
// spec here: https://samtools.github.io/hts-specs/SAMv1.pdf
//
// Example use with basic level parameter:
//   var bgzfFile bytes.Buffer
//   w, err := NewWriter(&bgzfFile, flate.DefaultCompression)
//   n, err := w.Write([]byte("Foo bar"))
//   err = w.Close()
//
// Example use with more configuration parameters:
//   var bgzfFile bytes.Buffer
//   w, err := NewWriterParams(
//     &bgzfFile,
//     flate.DefaultCompression,
//     DefaultUncompressedBlockSize,
//     zlibng.RLEStrategy,
//     0,
//   )
//   n, err := w.Write([]byte("Foo bar"))
//   err = w.Close()
//
// Example use with multiple compression shards:
//   // In goroutine 1
//   var shard1 bytes.Buffer
//   w, err := NewWriter(&shard1, flate.DefaultCompression)
//   n, err := w.Write([]byte("Foo bar"))
//   err = w.CloseWithoutTerminator()
//
//   // In goroutine 2
//   var shard2 bytes.Buffer
//   w, err := NewWriter(&shard2, flate.DefaultCompression)
//   n, err := w.Write([]byte(" baz!"))
//   err = w.Close()  // Terminator goes at the end of the last shard.
//
//   // Merge shards into final .bgzfFile.
//   var bgzfFile bytes.Buffer
//   _, err := io.Copy(&bgzfFile, &shard1)
//   _, err = io.Copy(&bgzfFile, &shard2)
package bgzf

import (
	"bytes"
	"fmt"
	"io"

	"github.com/grailbio/base/compress/libdeflate"
	"v.io/x/lib/vlog"
)

const (
	// DefaultUncompressedBlockSize is the default bgzf
	// uncompressedBlockSize chosen by both sambamba and biogo.  See
	// the SAM/BAM specification for details.
	DefaultUncompressedBlockSize = 0x0ff00

	// MaxUncompressedBlockSize is the largest legal value for
	// uncompressedBlockSize.  Illumina's Nextseq machines use this
	// value when creating .bcl.bgzf files.
	MaxUncompressedBlockSize = 0x10000

	// compressedBlockSize is the maximum size of the compressed data
	// for a Bgzf block.  See the SAM/BAM specification for details.
	compressedBlockSize = 0x10000
)

var (
	// bgzfExtra goes into the gzip's Extra subfield, with subfield
	// ids: 66, 67, and length 2.  See the SAM/BAM spec.
	bgzfExtra       = [...]byte{66, 67, 2, 0, 0, 0}
	bgzfExtraPrefix = [...]byte{66, 67, 2, 0}

	// terminiator is the Bgzf EOF terminator.  It belongs at the end
	// of a valid Bgzf file.  See the SAM/BAM spec.
	terminator = []byte{
		0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43,
		0x02, 0x00, 0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
	}
)

// compressFactory is an interface for creating a compressed gzip
// writer.  We use this so that we can have a cgo and non-cgo
// implementation of Writer.  The cgo version can use one of two
// factories; it uses the klauspost factory if the user uses
// NewWriter.  If the user uses NewWriterParams, then the factory
// creates a zlibng.Writer because zlibng supports more configuration
// params.  We use a factory here so that the factory can keep its own
// pointer to the libdeflate.Writer or zlibng.Writer so that the factory can
// use Reset() when possible instead of creating a new writer for each
// call to create().
type compressFactory interface {
	create(io.Writer) (io.WriteCloser, error)
}

type deflateFactory struct {
	level    int
	dfWriter *libdeflate.Writer
}

func (c *deflateFactory) create(w io.Writer) (io.WriteCloser, error) {
	if c.dfWriter == nil {
		var err error
		c.dfWriter, err = libdeflate.NewWriterLevel(w, c.level)
		if err != nil {
			return nil, err
		}
	} else {
		c.dfWriter.Reset(w)
	}
	c.dfWriter.Header.Extra = make([]byte, len(bgzfExtra))
	copy(c.dfWriter.Header.Extra[:], bgzfExtra[:])
	c.dfWriter.Header.OS = 0xff // Unknown OS value

	return c.dfWriter, nil
}

// Writer compresses data into .bgzf format.  The .bgzf format
// consists of gzip blocks concatenated together.  Each gzip block has
// an uncompressed size of at most 64KB.  The .bgzf format adds an
// Extra header field to each of the gzip headers; the Extra field
// contains the size of the uncompressed block in bytes - 1.  The
// payload data of the .bgzf file is equal to the in-order
// concatenation of all the uncompressed payloads of the gzip blocks.
// A .bgzf file also contains an EOF terminator at the end of the
// file.
type Writer struct {
	factory          compressFactory
	uncompressedSize int
	xfl              int
	w                io.Writer
	original         bytes.Buffer
	compressed       bytes.Buffer
	writer           io.WriteCloser
	coffset          uint64 // starting file position of the current gzip block
}

// NewWriter returns a new .bgzf writer with the given compression
// level.  Returns an nil, error if there is a problem.
func NewWriter(w io.Writer, level int) (*Writer, error) {
	return &Writer{
		factory:          &deflateFactory{level, nil},
		uncompressedSize: DefaultUncompressedBlockSize,
		xfl:              -1,
		w:                w,
	}, nil
}

// Writes buf to the .bgzf payload.  Returns the number of bytes
// consumed from buf and any error encountered.
func (w *Writer) Write(buf []byte) (int, error) {
	for i := 0; i < len(buf); {
		// Write one block at a time to avoid creating an entire copy of the input
		// buf.
		end := len(buf)

		// Now that libdeflate.Writer is in the picture, it's cleaner to account
		// for straggler bytes from the previous bgzf.Write() operation here.
		limit := i + w.uncompressedSize - w.original.Len()
		if limit < end {
			end = limit
		}
		n, _ := w.original.Write(buf[i:end])
		i += n
		if err := w.tryCompress(false); err != nil {
			return i, err
		}
	}
	return len(buf), nil
}

// CloseWithoutTerminator closes the current .bgzf block, but does not
// append the .bgzf terminator.  This output file is not a complete
// .bgzf file until the user calls Close().
func (w *Writer) CloseWithoutTerminator() error {
	return w.tryCompress(true)
}

// Close the current .bgzf block and also append the .bgzf terminator.
func (w *Writer) Close() error {
	if err := w.CloseWithoutTerminator(); err != nil {
		return err
	}
	_, err := w.w.Write([]byte(terminator))
	return err
}

// Removes a block from w.original, compresses the block, and
// appends the compressed block to c.output.buf.
func (w *Writer) tryCompress(compressRemainder bool) error {
	for w.original.Len() >= w.uncompressedSize || (compressRemainder && w.original.Len() > 0) {
		// Recreate gzip to start a new block
		var err error
		w.writer, err = w.factory.create(&w.compressed)
		if err != nil {
			return err
		}

		// Compress one block
		if w.original.Len() > 0 {
			_, err := w.writer.Write(w.original.Next(w.uncompressedSize))
			if err != nil {
				return err
			}
		}
		if err := w.writer.Close(); err != nil {
			return err
		}

		// Edit gzip header where necessary.
		b := w.compressed.Bytes()

		// Replace XFL value if configured.
		if w.xfl >= 0 {
			offset := 8 // This is the offset of the XFL field in the gzip header.
			b[offset] = byte(w.xfl)
		}

		// Replace bgzf BSIZE header with compressed length - 1.
		offset := 12 // This is the offset of the Extra field in the gzip header.
		bsize := w.compressed.Len() - 1
		if bsize >= compressedBlockSize {
			return fmt.Errorf("bgzf compressed block is too big: %d > %d", bsize,
				compressedBlockSize)
		}
		if w.compressed.Len() < (offset + len(bgzfExtra)) {
			vlog.Fatalf("compressed length is too short: %d < %d", w.compressed.Len(),
				offset+len(bgzfExtra))
		}
		if !bytes.Equal(b[offset:offset+len(bgzfExtraPrefix)], bgzfExtraPrefix[:]) {
			vlog.Fatalf("could not find bgzf extra prefix")
		}
		b[offset+4] = byte(bsize)
		b[offset+5] = byte(bsize >> 8)

		// Write out the compressed block.
		sz := w.compressed.Len()
		if _, err := w.compressed.WriteTo(w.w); err != nil {
			return err
		}
		w.coffset += uint64(sz)
	}
	return nil
}

// VOffset returns the virtual-offset of the next byte to be written.
func (w *Writer) VOffset() uint64 {
	return w.coffset<<16 | uint64(w.original.Len())
}
