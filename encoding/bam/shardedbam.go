package bam

import (
	"bytes"
	"fmt"
	"io"
	"sync"

	"github.com/biogo/hts/sam"
	"github.com/grailbio/base/syncqueue"
	"github.com/klauspost/compress/gzip"
	"v.io/x/lib/vlog"
)

// ShardedBAMWriter provides a way to write a BAM file as a sequence
// of shards.  Each of the shards has a sequentially increasing shard
// number starting at 0.  The ShardedBAMWriter writes these shards to
// a bam file in the order of their shard numbers.
//
// The user is responsible for creating each of the shards with the
// correct shard number (starting at 0), adding the records to the
// shards, and then closing the shards the the writers.
//
// To create a shard, the user must first create a
// ShardedBAMCompressor with GetCompressor.  The user can then ask the
// compressor to start a new shard, and then add the records to the
// shard.  When the user is finished adding the records to a shard,
// the user calls CloseShard(), which will send the shard to the
// ShardedBAMWriter.  The user can then start over with a new shard.
//
// To take advantage of multiple cores, the user can create one
// ShardedBAMCompressor for each thread, and then each thread can
// create a shard and add records to its shard independent of the
// other threads.  Compression occurs during AddRecord and CloseShard,
// so each of the cores can compress independently.
//
// Example use of ShardedBAMWriter:
//
//   f, _ := os.Create("output.bam")
//   w,err := NewShardedBAMWriter(f, gzip.DefaultCompression, 10, header)
//   c1 := w.GetCompressor()
//   c2 := w.GetCompressor()
//
//   // Each of the shards may be populated and closed concurrently.
//   c2.StartShard(0)
//   c2.AddRecord(record0)
//   c2.CloseShard()
//
//   c2.StartShard(2)
//   c2.AddRecord(record3)
//   c2.AddRecord(record4)
//   c2.CloseShard()
//
//   c1.StartShard(1)
//   c1.AddRecord(record1)
//   c1.AddRecord(record2)
//   c1.CloseShard()
//
//   w.Close()

const (
	// uncompressedBlockSize is the maximum size of the uncompressed
	// data for a Bgzf block.  Technically, the maximum size is 2^16,
	// but this threshold is slightly lower to avoid exceeding the
	// maximum size of the compressed block.  This value is the same
	// in both sambamba and biogo.  See the SAM/BAM specification for
	// details.
	uncompressedBlockSize = 0x0ff00

	// compressedBlockSize is the maximum size of the compressed data
	// for a Bgzf block.  See the SAM/BAM specification for details.
	compressedBlockSize = 0x10000
)

var (
	// bgzfExtra goes into the gzip's Extra subfield, with subfield
	// ids: 66, 67, and length 2.  See the SAM/BAM spec.
	bgzfExtra       = []byte{66, 67, 2, 0, 0, 0}
	bgzfExtraPrefix = []byte(bgzfExtra[:4])

	// magicBlock is the Bgzf EOF terminator.  It belongs at the end
	// of a valid Bgzf file.  See the SAM/BAM spec.
	magicBlock = []byte{
		0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43,
		0x02, 0x00, 0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
	}
)

// ShardedBAMCompressor contains the state of an in-progress
// compressed shard.  A caller should create a ShardedBAMCompressor
// using ShardedBAMWriter.GetCompressor().  The ShardedBAMCompressor
// will compress the records and store the compressed bytes until the
// caller is finished with the shard.  When the caller is finished
// adding records, the caller should call CloseShard().  More than one
// ShardedBAMCompressor can exist at once, and they can all compress
// records in parallel with each other.
type ShardedBAMCompressor struct {
	writer       *ShardedBAMWriter
	uncompressed bytes.Buffer
	gzWriter     *gzip.Writer
	output       *shardedBAMBuffer
}

// StartShard begins a new shard with the specified shard number.  If
// the compressor still has data from the previous shard, it will
// crash.
func (c *ShardedBAMCompressor) StartShard(shardNum int) {
	if c.output != nil {
		vlog.Fatalf("existing shard still in progress")
	}
	c.output = &shardedBAMBuffer{
		// The client sees shardNum start at 0, but internally,
		// ShardedBAMWriter needs a shard for the header, so add 1 to
		// the shardNum here, which allows the header to call
		// StartShard(-1).
		shardNum: shardNum + 1,
	}

}

// addHeader adds a sam header to the current shard.  This must be
// called on the first bam shard of the file so that the header will
// be at the correct position in the bam file; NewShardedBAMWriter
// takes care of that for users.
func (c *ShardedBAMCompressor) addHeader(h *sam.Header) error {
	if err := h.EncodeBinary(&c.uncompressed); err != nil {
		return err
	}
	if c.uncompressed.Len() > uncompressedBlockSize {
		if err := c.compress(false); err != nil {
			return err
		}
	}
	return nil
}

// AddRecord adds a sam record to the current in-progress shard.
func (c *ShardedBAMCompressor) AddRecord(r *sam.Record) error {
	if err := Marshal(r, &c.uncompressed); err != nil {
		return err
	}
	if c.uncompressed.Len() > uncompressedBlockSize {
		if err := c.compress(false); err != nil {
			return err
		}
	}
	return nil
}

// CloseShard finalizes the in-progress shard, and passes the
// compressed data to its parent ShardedBAMWriter.  It removes the
// current shard from the compressor and prepares the compressor for
// the next call to StartShard().
//
// The ShardedBAMWriter will buffer shards up to its queue size, so
// the caller must be careful about how out of order it is when
// calling CloseShard(), otherwise, calls to CloseShard() will block.
func (c *ShardedBAMCompressor) CloseShard() error {
	if c.uncompressed.Len() > 0 {
		if err := c.compress(true); err != nil {
			return err
		}
	}
	f := c.output
	c.output = nil
	return c.writer.addShard(f)
}

// Removes a block from c.uncompressed, compresses the block, and
// appends the compressed block to c.output.buf.
func (c *ShardedBAMCompressor) compress(compressRemainder bool) error {
	for c.uncompressed.Len() >= uncompressedBlockSize || (compressRemainder && c.uncompressed.Len() > 0) {
		// Setup gzip
		if c.gzWriter == nil {
			w, err := gzip.NewWriterLevel(&c.output.buf, c.writer.gzLevel)
			if err != nil {
				return err
			}
			c.gzWriter = w
		} else {
			c.gzWriter.Reset(&c.output.buf)
		}
		c.gzWriter.Header = gzip.Header{
			Extra: bgzfExtra,
			OS:    0xff, // Unknown OS value
		}

		// Compress one block
		count := 0
		originalSize := c.output.buf.Len()
		for c.uncompressed.Len() > 0 && count < uncompressedBlockSize {
			n, err := c.gzWriter.Write(c.uncompressed.Next(uncompressedBlockSize - count))
			if err != nil {
				return err
			}
			count += n
		}

		if err := c.gzWriter.Close(); err != nil {
			return err
		}
		newSize := c.output.buf.Len()

		// Fix up bgzf prefix
		sz := newSize - originalSize - 1
		if sz >= compressedBlockSize {
			return fmt.Errorf("bgzf compressed block is too big: %d > %d", sz, compressedBlockSize)
		}
		found := false
		b := c.output.buf.Bytes()
		for i := originalSize; i < newSize-len(bgzfExtra); i++ {
			if bytes.Compare(b[i:i+len(bgzfExtraPrefix)], bgzfExtraPrefix) == 0 {
				b[i+4] = byte(sz)
				b[i+5] = byte(sz >> 8)
				found = true
				break
			}
		}
		if !found {
			vlog.Fatalf("could not find bgzf extra prefix")
		}
	}
	return nil
}

// ShardedBAMBuffer represents a shard of the final output bam file.
// The shardNum should be numbered sequentially from 0 to N.
type shardedBAMBuffer struct {
	buf      bytes.Buffer
	shardNum int
}

// ShardedBAMWriter writes out ShardedBAMBuffers in the order of their
// shard numbers.
type ShardedBAMWriter struct {
	w         io.Writer
	gzLevel   int
	queue     *syncqueue.OrderedQueue
	waitGroup sync.WaitGroup
	err       error
}

// NewShardedBAMWriter creates a new ShardedBAMWriter that writes the
// output bam to w.
func NewShardedBAMWriter(w io.Writer, gzLevel, queueSize int, header *sam.Header) (*ShardedBAMWriter, error) {
	bw := ShardedBAMWriter{
		w:       w,
		gzLevel: gzLevel,
		queue:   syncqueue.NewOrderedQueue(queueSize),
	}

	c := bw.GetCompressor()
	c.StartShard(-1)
	if err := c.addHeader(header); err != nil {
		return nil, err
	}
	if err := c.CloseShard(); err != nil {
		return nil, err
	}

	bw.waitGroup.Add(1)
	go func() {
		defer bw.waitGroup.Done()
		bw.writeShards()
	}()

	return &bw, nil
}

// GetCompressor returns a child ShardedBAMCompressor.
func (bw *ShardedBAMWriter) GetCompressor() *ShardedBAMCompressor {
	return &ShardedBAMCompressor{
		writer: bw,
	}
}

// addShard inserts a shard into the ShardedBAMWriter.  The
// ShardedBAMWriter writes the shards out in sequential order by
// shardNum, so if shard is not the next shard to be written, then
// ShardedBAMWriter will buffer the shard until it has all the
// preceding shards and can write the shards in sequential order.
func (bw *ShardedBAMWriter) addShard(shard *shardedBAMBuffer) error {
	return bw.queue.Insert(shard.shardNum, shard)
}

func (bw *ShardedBAMWriter) writeShards() {
	for {
		entry, ok, err := bw.queue.Next()
		if err != nil {
			bw.err = err
			break
		}
		if !ok {
			break
		}
		shard := entry.(*shardedBAMBuffer)
		_, err = shard.buf.WriteTo(bw.w)
		if err != nil {
			bw.err = err
			bw.queue.Close(err)
			return
		}
	}
}

// Close the bam file.  This should be called only after all shards
// have been added with WriteShard.  Returns an error of failure.
func (bw *ShardedBAMWriter) Close() error {
	err := bw.queue.Close(nil)
	bw.waitGroup.Wait()
	if bw.err != nil {
		return bw.err
	}
	if err != nil {
		return err
	}

	_, err = bw.w.Write([]byte(magicBlock))
	return err
}
