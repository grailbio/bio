package sorter

import (
	"encoding/binary"
	"fmt"
	"io"
	"math"
	"sync"
	"sync/atomic"

	"github.com/golang/snappy"
	"github.com/grailbio/base/errors"
	"github.com/grailbio/base/file"
	"github.com/grailbio/base/recordio"
	"github.com/grailbio/base/vcontext"
	"github.com/grailbio/bio/biopb"
	"github.com/grailbio/hts/bam"
	"github.com/grailbio/hts/sam"
	"v.io/x/lib/vlog"
)

// Sortshard format is used by temp files used during BAM/PAM sorting and
// merging. It logically stores the same information as BAM. Compared to BAM,
// sequential reads and writes are faster. However, sortshard does not support
// random accesses, and the file size is usually larger because it uses snappy
// instead of flate. The file is a recordio, where one recordio block stores
// list of serialized BAM records, in the following format: There is no padding
// between records.
//
//   key sortEntry   // for sorting records.
//   bytes uint32  // Size of the record, in bytes.
//   data [bytes]byte  // sam.Record serialized in BAM format
//
// Each recordio block is approx. sortShardBlockSize bytes long,
// pre-compression.
//
// The record trailer stores a serialized SortShardIndex.  This proto tells
// whether the preceding data blocks are compressed or not, and it also stores
// the file offset of each block, in case one wants to seek to it directly.
//
// TODO(saito) The current codebase performs snappy compression on its own. Do
// snappy as part of recordio transformer.
type sortShardBlock []byte

// sortShardBuf stores contents of a recordio block during writes.
type sortShardBuf struct {
	buf       sortShardBlock
	remaining []byte    // part of buf[].
	lastKey   sortEntry // last key in the block, set iff nRecords>0.
	startKey  sortEntry // first key in the block, set iff nRecords>0.
	nRecords  int       // # of records stored in buf.
}

const sortShardBlockSize = 1 << 20   // size of one sortShardBlock.buf
const sortShardRecordHeaderSize = 12 // 8 byte sortEntry + 4 byte record size.

// Class for producing a SortShard file
//
// Example:
//   err := errorReporter{}
//   pool := newSortShardBlockPool()
//   w := newSortShardWriter(..., pool, &err)
//   for {
//     ...
//     // Call w.add for each record
//     w.add(key, body)
//   }
//   w.finish()
//   if err.err() != nil { panic(err) }
type sortShardWriter struct {
	rawOut          io.Writer       // File writer.
	rio             recordio.Writer // recordio wrapper for out.
	err             *errors.Once
	writeBlockIndex bool

	curBlock sortShardBuf // The block currently written to in add().
	pool     *sortShardBlockPool

	indexMu sync.Mutex
	index   biopb.SortShardIndex
}

func (w *sortShardWriter) newBuf() sortShardBuf {
	buf := w.pool.getBuf()
	return sortShardBuf{
		buf:       buf,
		remaining: buf,
	}
}

// Create a new sortShardWriter. Any error is reported through errReporter.  If
// parameters writeBlockIndex and header are set, the index block will contain
// the full block offset and bam header information.
func newSortShardWriter(out io.Writer,
	snappy, writeBlockIndex bool,
	header *sam.Header, //may be null.
	pool *sortShardBlockPool, errReporter *errors.Once) *sortShardWriter {
	w := &sortShardWriter{
		rawOut:          out,
		writeBlockIndex: writeBlockIndex,
		err:             errReporter,
		pool:            pool,
		index:           biopb.SortShardIndex{Snappy: snappy},
	}
	w.curBlock = w.newBuf()
	if header != nil {
		var err error
		if w.index.EncodedBamHeader, err = bam.MarshalHeader(header); err != nil {
			// TODO(saito) propagate error up.
			vlog.Fatalf("Failed to encode header: %v", err)
		}
	}
	w.rio = recordio.NewWriter(out, recordio.WriterOpts{
		Marshal: func(scratch []byte, v interface{}) ([]byte, error) {
			b := v.(sortShardBuf)
			return b.buf, nil
		},
		Index: func(loc recordio.ItemLocation, v interface{}) error {
			b := v.(sortShardBuf)
			if loc.Item != 0 { // This is a single-item-per-block recordio
				vlog.Fatal(loc)
			}
			if w.writeBlockIndex {
				bi := biopb.SortShardBlockIndex{
					StartKey:   uint64(b.startKey.coord),
					FileOffset: loc.Block,
					NumRecords: uint32(b.nRecords),
				}
				w.indexMu.Lock()
				w.index.Blocks = append(w.index.Blocks, bi)
				w.indexMu.Unlock()
			}
			w.pool.putBuf(b.buf)
			return nil
		},
	})
	w.rio.AddHeader(recordio.KeyTrailer, true)
	return w
}

// Add a BAM-serialized record to the buffer.
func (w *sortShardWriter) add(key sortEntry) {
	if key.coord == invalidCoord {
		vlog.Fatalf("Key: %v", key)
	}
	if key.compare(w.curBlock.lastKey) < 0 {
		vlog.Errorf("Key %v decreased, last %v", key, w.curBlock.lastKey)
		panic("key")
	}
	w.curBlock.lastKey = key
	if w.tryAdd(key) {
		return // Common case.
	}
	// For all but the last block in the shard, the buffer size is fixed at w.shardSize.
	// w.remaining = nil
	w.flush()
	vlog.VI(1).Infof("Starting new buffer at key: %+v", key)
	if !w.tryAdd(key) {
		vlog.Fatalf("Key: %v", key)
	}
}

// Flush any pending data to the file. An error is reported trhough w.err. "w"
// becomes invalid after the call.
func (w *sortShardWriter) finish() {
	w.flush()
	w.pool.putBuf(w.curBlock.buf)
	w.curBlock.buf = nil

	// Write the index block.  Index block is never compressed (the snappy flag
	// itself is embedded in the index).
	w.rio.Wait()
	indexBytes, err := w.index.Marshal()
	if err != nil {
		vlog.Fatalf("Failed to marshal index: %v", err)
	}
	w.rio.SetTrailer(indexBytes)
	w.err.Set(w.rio.Finish())
}

func (w *sortShardWriter) flush() {
	if w.curBlock.nRecords == 0 {
		return
	}
	b := w.curBlock
	w.curBlock = w.newBuf()

	bytes := b.bytes()
	if w.index.Snappy {
		compressBuf := w.pool.getBuf()
		out := snappy.Encode(compressBuf, bytes)
		w.pool.putBuf(b.buf)
		b.buf = out
	}
	w.rio.Append(b)
	w.rio.Flush()
}

// Returns a buffer that contains records added so far.
func (b *sortShardBuf) bytes() []byte {
	n := len(b.buf) - len(b.remaining)
	return b.buf[:n]
}

func (w *sortShardWriter) tryAdd(key sortEntry) bool {
	b := &w.curBlock
	if len(b.remaining) < sortShardRecordHeaderSize+len(key.body) {
		if len(b.remaining) >= sortShardRecordHeaderSize {
			binary.LittleEndian.PutUint64(b.remaining[:8], uint64(invalidCoord))
			binary.LittleEndian.PutUint32(b.remaining[8:12], 0xffffffff)
		}
		return false
	}
	binary.LittleEndian.PutUint64(b.remaining[:8], uint64(key.coord))
	binary.LittleEndian.PutUint32(b.remaining[8:12], uint32(len(key.body)))
	copy(b.remaining[sortShardRecordHeaderSize:], key.body)
	b.remaining = b.remaining[sortShardRecordHeaderSize+len(key.body):]

	if b.nRecords == 0 {
		b.startKey = key
	}
	b.nRecords++
	w.index.NumRecords++
	return true
}

// Class for extracting records in a sortShardBlock.
//
// Example:
//   for r := newSortShardBlockParser(buf); !r.done(); r.next() {
//      vlog.Infof("Key %v", r.key())
//   }
type sortShardBlockParser struct {
	curKey sortEntry // The currennt key. EOD iff curKey.key==invalidKey.
	buf    []byte    // Records that remain to be read.
}

func (r *sortShardBlockParser) reset(buf sortShardBlock) {
	r.buf = []byte(buf)
	r.next()
	if r.done() {
		vlog.Fatalf("empty buf: %v", len(buf))
	}
}

func (r *sortShardBlockParser) next() {
	if len(r.buf) <= sortShardRecordHeaderSize {
		// The header is chopped at the end.
		r.curKey = sortEntry{invalidCoord, nil}
		return
	}
	r.curKey.coord = recCoord(binary.LittleEndian.Uint64(r.buf[:8]))
	if r.curKey.coord == invalidCoord {
		return
	}
	recLen := binary.LittleEndian.Uint32(r.buf[8:12])
	if uint64(len(r.buf)) < sortShardRecordHeaderSize+uint64(recLen) {
		r.curKey.coord = invalidCoord
		return
	}
	r.curKey.body = make([]byte, recLen)
	copy(r.curKey.body, r.buf[sortShardRecordHeaderSize:sortShardRecordHeaderSize+recLen])
	r.buf = r.buf[sortShardRecordHeaderSize+recLen:]
}

func (r *sortShardBlockParser) done() bool {
	return r.curKey.coord == invalidCoord
}

func (r *sortShardBlockParser) key() sortEntry {
	if r.done() {
		vlog.Fatal(r)
	}
	return r.curKey
}

// Class for reading a SortShard file.
//
// Example:
//   err := errorReporter{}
//   pool := newSortShardBlockPool()
//   r := newSortShardReader(..., pool, &err)
//   for r.scan() {
//     use r.key()
//   }
//   if err.err() != nil { panic(err) }
type sortShardReader struct {
	path    string
	rawIn   file.File
	rio     recordio.Scanner
	index   biopb.SortShardIndex
	pool    *sortShardBlockPool
	err     *errors.Once
	lastKey sortEntry // last key read.

	parser sortShardBlockParser
	buf    []byte
	ch     chan sortShardBlock
	// draining becomes 1 on drain(). It tells asyncRead goroutine to finish
	// asap. It must be accessed via acquire-loads+release-stores.
	draining int32
}

// Read the index block in the sortshard file.
func readSortShardIndex(rio recordio.Scanner) (biopb.SortShardIndex, error) {
	index := biopb.SortShardIndex{}
	header := rio.Header()
	if !header.HasTrailer() {
		return index, fmt.Errorf("no index found in sortshard file (header: %+v, version %+v)", header, rio.Version())
	}
	if err := index.Unmarshal(rio.Trailer()); err != nil {
		return index, err
	}
	return index, nil
}

type shardReaderOpts struct {
	// If index == nil, newSortShardReader() will load the index from the end of
	// the shard file.
	index *biopb.SortShardIndex

	// Range of blocks [startOffset, limitOffset) to read from.
	startOffset int64
	limitOffset int64
}

// Create a reader for reading SortShard file "path". Any error is reported
// through errReporter.
func newSortShardReader(path string,
	pool *sortShardBlockPool,
	errReporter *errors.Once,
	optList ...shardReaderOpts) *sortShardReader {
	opts := shardReaderOpts{limitOffset: math.MaxInt64}
	if len(optList) > 1 {
		vlog.Fatalf(">1 optlist: %v", optList)
	}
	if len(optList) > 0 {
		opts = optList[0]
	}
	r := &sortShardReader{
		path: path,
		pool: pool,
		err:  errReporter,
		// The parser is initially at done() state.
		parser: sortShardBlockParser{curKey: sortEntry{invalidCoord, nil}},
		ch:     make(chan sortShardBlock),
	}

	ctx := vcontext.Background()
	cleanupOnError := func(err error) *sortShardReader {
		r.err.Set(err)
		close(r.ch)
		if r.rawIn != nil {
			r.err.Set(r.rawIn.Close(ctx))
		}
		return r
	}
	var err error
	r.rawIn, err = file.Open(ctx, path)
	if err != nil {
		return cleanupOnError(err)
	}
	r.rio = recordio.NewScanner(r.rawIn.Reader(ctx), recordio.ScannerOpts{})
	if opts.index == nil {
		r.index, err = readSortShardIndex(r.rio)
		if err != nil {
			return cleanupOnError(err)
		}
	} else {
		r.index = *opts.index
	}
	vlog.VI(0).Infof("%v: created shard reader, range=[%v,%v)", path, opts.startOffset, opts.limitOffset)
	if opts.startOffset > 0 {
		r.rio.Seek(recordio.ItemLocation{uint64(opts.startOffset), 0})
	}
	go func() {
		r.asyncRead()
		if r.rawIn != nil {
			r.err.Set(r.rawIn.Close(ctx))
		}
		close(r.ch)
	}()
	return r
}

func (r *sortShardReader) scan() bool {
	if !r.parser.done() {
		r.parser.next()
	}
	if r.parser.done() {
		if r.buf != nil {
			r.pool.putBuf(r.buf)
		}
		buf, ok := <-r.ch
		if !ok { // shard contains no data??
			r.buf = nil
			return false
		}
		r.buf = buf
		r.parser.reset(buf)
	}
	if r.parser.key().compare(r.lastKey) < 0 {
		vlog.Fatalf("Key %v decreased, last %v", r.parser.key(), r.lastKey)
	}
	r.lastKey = r.parser.key()
	return true
}

// Drain should be called when quitting reads before reaching the end of
// shard. It cleans up the reader state.  It's ok to call drain() after
// successful end of reads.
func (r *sortShardReader) drain() {
	go func() {
		n := 0
		atomic.StoreInt32(&r.draining, 1)
		for range r.ch {
			n++
		}
		vlog.VI(1).Infof("drain %v: dropped %d blocks", r.path, n)
	}()
}

// Return the key of the current record.
//
// REQUIRES: scan() returned true.
func (r *sortShardReader) key() sortEntry {
	return r.parser.key()
}

// Read a sequence of raw sortShardBlocks and send them to "r.ch".
func (r *sortShardReader) asyncRead() {
	if r.rio == nil {
		// Failed already
		return
	}
	for r.rio.Scan() && atomic.LoadInt32(&r.draining) == 0 {
		sorted := r.pool.getBuf()
		rioData := r.rio.Get().([]byte)
		if r.index.Snappy {
			var err error
			sorted, err = snappy.Decode(sorted, rioData)
			if err != nil {
				r.err.Set(err)
				break
			}
		} else {
			if len(sorted) < len(rioData) {
				// This shouldn't happen in practice, since the writer will limit the
				// bufsize to be sortShardBlockSize.
				sorted = make([]byte, len(rioData))
			}
			sorted = sorted[:len(rioData)]
			copy(sorted, rioData)
		}
		r.ch <- sorted // This may block
	}
	r.err.Set(r.rio.Err())
}

// Freepool of sortShardBlocks.
type sortShardBlockPool struct {
	sync.Pool
}

// Get a sortShardBlock from the pool. The caller should call putBuf(buf) after use.
func (p *sortShardBlockPool) getBuf() sortShardBlock {
	b := p.Get().(sortShardBlock)
	if cap(b) < sortShardBlockSize {
		b = make(sortShardBlock, sortShardBlockSize)
	} else {
		b = b[:sortShardBlockSize]
	}
	return b
}

func (p *sortShardBlockPool) putBuf(b sortShardBlock) {
	p.Put(b)
}

func newSortShardBlockPool() *sortShardBlockPool {
	return &sortShardBlockPool{sync.Pool{New: func() interface{} { return sortShardBlock{} }}}
}
