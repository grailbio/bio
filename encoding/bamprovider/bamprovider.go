package bamprovider

import (
	"fmt"
	"io"
	"sync"

	"github.com/biogo/hts/bam"
	"github.com/biogo/hts/bgzf"
	"github.com/biogo/hts/bgzf/index"
	"github.com/biogo/hts/sam"
	"github.com/grailbio/base/errorreporter"
	"github.com/grailbio/base/file"
	"github.com/grailbio/base/vcontext"
	"github.com/grailbio/bio/biopb"
	gbam "github.com/grailbio/bio/encoding/bam"
	"v.io/x/lib/vlog"
)

const (
	// DefaultBytesPerShard is the default value for GenerateShardsOpts.BytesPerShard
	DefaultBytesPerShard = int64(128 << 20)
	// DefaultMinBasesPerShard is the default value for GenerateShardsOpts.MinBasesPerShard
	DefaultMinBasesPerShard = 10000
)

// BAMProvider implements Provider for BAM files.  Both BAM and the index
// filenames are allowed to be S3 URLs, in which case the data will be read from
// S3. Otherwise the data will be read from the local filesystem.
type BAMProvider struct {
	// Path of the *.bam file. Must be nonempty.
	Path string
	// Index is the pathname of *.bam.bai file. If "", Path + ".bai"
	Index string
	err   errorreporter.T

	mu        sync.Mutex
	nActive   int
	freeIters []*bamIterator
	header    *sam.Header
}

type bamIterator struct {
	provider *BAMProvider
	in       file.File
	reader   *bam.Reader
	index    *bam.Index
	// Offset of the first record in the file.
	firstRecord bgzf.Offset
	// Half-open coordinate range to read.
	startAddr, limitAddr biopb.Coord

	active bool
	err    error
	next   *sam.Record
	done   bool
}

func (b *BAMProvider) indexPath() string {
	index := b.Index
	if index == "" {
		index = b.Path + ".bai"
	}
	return index
}

// GetHeader implements the Provider interface.
func (b *BAMProvider) GetHeader() (*sam.Header, error) {
	b.mu.Lock()
	defer b.mu.Unlock()
	if b.header != nil {
		return b.header, nil
	}

	ctx := vcontext.Background()
	reader, err := file.Open(ctx, b.Path)
	if err != nil {
		b.err.Set(err)
		return nil, err
	}
	defer reader.Close(ctx)
	bamReader, err := bam.NewReader(reader.Reader(ctx), 1)
	if err != nil {
		b.err.Set(err)
		return nil, err
	}
	defer bamReader.Close()
	b.header = bamReader.Header()
	return b.header, nil
}

// Close implements the Provider interface.
func (b *BAMProvider) Close() error {
	if b.nActive > 0 {
		vlog.Fatalf("%d iterators still active for %+v", b.nActive, b)
	}
	for _, iter := range b.freeIters {
		iter.internalClose()
	}
	b.freeIters = nil
	return b.err.Err()
}

func (b *BAMProvider) freeIterator(i *bamIterator) {
	if !i.active {
		vlog.Fatal(i)
	}
	i.active = false
	if i.Err() != nil {
		// The iter may be invalid. Don't reuse it.
		i.internalClose() // Will set b.err
		i = nil
	}
	b.mu.Lock()
	if i != nil {
		b.freeIters = append(b.freeIters, i)
	}
	b.nActive--
	if b.nActive < 0 {
		vlog.Fatalf("Negative active count for %+v", b)
	}
	b.mu.Unlock()
}

// Return an unused iterator. If b.freeIters is nonempty, this function returns
// one from freeIters. Else, it opens the BAM file, creates a BAM reader and
// returns an iterator containing them. On error, returns an iterator with
// non-nil err field.
func (b *BAMProvider) allocateIterator() *bamIterator {
	b.mu.Lock()
	b.nActive++
	if len(b.freeIters) > 0 {
		iter := b.freeIters[len(b.freeIters)-1]
		iter.active = true
		iter.err = nil
		iter.done = false
		iter.next = nil
		b.freeIters = b.freeIters[:len(b.freeIters)-1]
		b.mu.Unlock()
		return iter
	}
	b.mu.Unlock()

	iter := bamIterator{
		provider: b,
		active:   true,
	}
	ctx := vcontext.Background()
	if iter.in, iter.err = file.Open(ctx, b.Path); iter.err != nil {
		return &iter
	}

	var indexIn file.File
	if indexIn, iter.err = file.Open(ctx, b.indexPath()); iter.err != nil {
		return &iter
	}
	defer indexIn.Close(ctx)
	if iter.index, iter.err = bam.ReadIndex(indexIn.Reader(ctx)); iter.err != nil {
		return &iter
	}
	if iter.reader, iter.err = bam.NewReader(iter.in.Reader(ctx), 1); iter.err != nil {
		return &iter
	}
	iter.firstRecord = iter.reader.LastChunk().End
	return &iter
}

// GenerateShards implements the Provider interface.
func (b *BAMProvider) GenerateShards(opts GenerateShardsOpts) ([]gbam.Shard, error) {
	header, err := b.GetHeader()
	if err != nil {
		return nil, err
	}
	if opts.BytesPerShard <= 0 {
		opts.BytesPerShard = DefaultBytesPerShard
	}
	if opts.MinBasesPerShard <= 0 {
		opts.MinBasesPerShard = DefaultMinBasesPerShard
	}
	if opts.Strategy == ByteBased {
		return gbam.GetByteBasedShards(
			b.Path, b.indexPath(), opts.BytesPerShard, opts.MinBasesPerShard, opts.Padding, opts.IncludeUnmapped)
	}
	return gbam.GetPositionBasedShards(
		header, 100000, opts.Padding, opts.IncludeUnmapped)
}

// GetFileShards implements the Provider interface.
func (b *BAMProvider) GetFileShards() ([]gbam.Shard, error) {
	header, err := b.GetHeader()
	if err != nil {
		return nil, err
	}
	return []gbam.Shard{gbam.UniversalShard(header)}, nil
}

// NewIterator implements the Provider interface.
func (b *BAMProvider) NewIterator(shard gbam.Shard) Iterator {
	iter := b.allocateIterator()
	if iter.err != nil {
		return iter
	}
	if shard.StartRef.ID() != shard.EndRef.ID() {
		iter.err = fmt.Errorf("For BAMProvider, start and limit ref ID must be the same, but got %v, %v",
			shard.StartRef, shard.EndRef)
	}
	iter.reset(shard.StartRef, shard.PaddedStart(), shard.EndRef, shard.PaddedEnd())
	return iter
}

// Reset the iterator to read the range [<startRef,startPos>, <endRef, endPos>).
func (i *bamIterator) reset(startRef *sam.Reference, startPos int, endRef *sam.Reference, endPos int) {
	header := i.reader.Header()
	i.startAddr = biopb.Coord{int32(startRef.ID()), int32(startPos), 0}
	i.limitAddr = biopb.Coord{int32(endRef.ID()), int32(endPos), 0}
	if i.startAddr.GE(i.limitAddr) {
		i.err = fmt.Errorf("start coord (%v) not before limit coord (%v)", i.startAddr, i.limitAddr)
		return
	}

	// Read the index and find the file offset at which <startRef,startPos> is
	// located.
	var offset bgzf.Offset
	var err error
	ref := startRef
	for {
		if ref == nil {
			offset, err = i.findUnmappedOffset()
			break
		}
		start := 0
		if ref.ID() == startRef.ID() {
			start = startPos
		}
		end := ref.Len()
		if ref.ID() == endRef.ID() {
			end = endPos
		}
		var found bool
		found, offset, err = i.findRecordOffset(ref, start, end)
		if err != nil || found {
			break
		}
		if ref.ID() == endRef.ID() {
			// No refs in range [startRef,endRef] has any index.  There's no record to
			// read.
			i.err = io.EOF
			return
		}
		// No index is found for this ref. Try the next ref.
		ref = header.Refs()[ref.ID()+1]
	}
	if err != nil {
		i.err = err
		return
	}
	i.err = i.reader.Seek(offset)
}

// Err implements the Iterator interface.
func (i *bamIterator) Err() error {
	if i.err == io.EOF {
		return nil
	}
	return i.err
}

// Close implements the Iterator interface.
func (i *bamIterator) Close() error {
	err := i.Err()
	i.provider.freeIterator(i)
	return err
}

// Find the the file offset at which the first unmapped sequence is
// stored. This function is conservative; it may return an offset that's smaller
// than absolutely necessary.
func (i *bamIterator) findUnmappedOffset() (bgzf.Offset, error) {
	// Iterate through the endpoint of each reference to find the
	// largest offset.
	header := i.reader.Header()
	var lastOffset bgzf.Offset
	foundRefs := false
	for _, r := range header.Refs() {
		chunks, err := i.index.Chunks(r, 0, r.Len())
		if err == index.ErrInvalid {
			// There are no reads on this reference, but don't worry about it.
			continue
		}
		if err != nil {
			return lastOffset, err
		}
		foundRefs = true
		c := chunks[len(chunks)-1]
		if c.End.File > lastOffset.File ||
			(c.End.File == lastOffset.File && c.End.Block > lastOffset.Block) {
			lastOffset = c.End
		}
	}
	if !foundRefs {
		return i.firstRecord, nil
	}
	return lastOffset, nil
}

// Find the the file offset at which the first record at coordinate <ref,pos> is
// stored. This function is conservative; it may return an offset that's smaller
// than absolutely necessary.
func (i *bamIterator) findRecordOffset(ref *sam.Reference, startPos, endPos int) (bool, bgzf.Offset, error) {
	chunks, err := i.index.Chunks(ref, startPos, endPos)
	if err == index.ErrInvalid || len(chunks) == 0 {
		// No reads for this interval: return an empty iterator.
		// This needs to be handled as a special case due to the current behavior of biogo.
		// Return the same 'eofIterator' to avoid unnecessary memory allocations, this
		// is likely an artifact of micro-benchmarking which uses smaller files which
		// are likely to hit this codepath.
		return false, bgzf.Offset{}, nil
	}
	if err != nil {
		return false, bgzf.Offset{}, err
	}
	return true, chunks[0].Begin, nil
}

func (i *bamIterator) Scan() bool {
	if !i.active {
		vlog.Fatal("Reusing iterator")
	}
	if i.err != nil {
		return false
	}
	for {
		i.next, i.err = i.reader.Read()
		if i.err != nil {
			return false
		}
		recAddr := gbam.CoordFromSAMRecord(i.next, 0)
		if recAddr.LT(i.startAddr) {
			continue
		}
		return recAddr.LT(i.limitAddr)
	}
}

func (i *bamIterator) Record() *sam.Record {
	return i.next
}

func (i *bamIterator) internalClose() {
	if i.reader != nil {
		if err := i.reader.Close(); err != nil && i.err == nil {
			i.err = err
		}
		i.reader = nil
	}
	if i.in != nil {
		if err := i.in.Close(vcontext.Background()); err != nil && i.err == nil {
			i.err = err
		}
		i.in = nil
	}
	i.provider.err.Set(i.Err())
}
