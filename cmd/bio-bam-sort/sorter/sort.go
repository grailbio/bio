package sorter

import (
	"bytes"
	"crypto/sha256"
	"encoding/binary"
	"fmt"
	"io/ioutil"
	"os"
	"runtime"
	"sort"
	"sync"

	"github.com/biogo/store/llrb"
	"github.com/grailbio/base/errors"
	"github.com/grailbio/base/file"
	"github.com/grailbio/base/vcontext"
	gbam "github.com/grailbio/bio/encoding/bam"
	"github.com/grailbio/hts/bam"
	"github.com/grailbio/hts/bgzf"
	"github.com/grailbio/hts/sam"
	"v.io/x/lib/vlog"
)

// DefaultSortBatchSize is the default number of records to keep in
// memory before resorting to external sorting.
const DefaultSortBatchSize = 1 << 20

// DefaultParallelism is the default value for SortOptions.Parallelism.
const DefaultParallelism = 2

// SortOptions controls options passed to the toplevel Sort.
type SortOptions struct {
	// ShardIndex must be a number unique to this sorter, across all sorters for
	// shards that are eventually merged into one BAM or PAM file.
	//
	// ShardIndex defines the sort order of reads at the same (ref,pos), but on
	// different Sorters. If ShardIndex==0, it is set to sha(sortshardpath).
	ShardIndex uint32

	// SortBatchSize is the number of sam.Records to keep in memory before
	// resorting to external sorting.  Not for general use; the default value
	// should suffice for most applications.
	SortBatchSize int

	// MaxParallelism limits the number of background sorts. Max memory
	// consumption of the sorter grows linearly with this value. If <= 0,
	// DefaultMaxParallelism is used.
	Parallelism int

	// NoCompressTmpFiles, if false (default), compress sortshards using snappy.
	// Compression is a big win on an EC2 EBS. It will slow sort down by a minor
	// degree on fast NVMe disks.
	NoCompressTmpFiles bool

	// TmpDir defines the directory to store temp files created during merge.  ""
	// means the system default, usually /tmp.
	TmpDir string
}

// recCoord encodes reference id, alignment position, and the reverse flag.  Sort
// order of readCoord is the same as the SAM "coord" sort order.
type recCoord uint64

type sortBatch struct {
	sortTieBreaker uint64
	recs           []sortEntry
}

// sortEntry is an efficient encoding of sam.Record sort order.
type sortEntry struct {
	coord recCoord
	body  []byte // Contains the full record in bam serialized form.
}

func (k sortEntry) String() string {
	refid, pos, reverse := parseCoord(k.coord)
	return fmt.Sprintf("(%d,%d,%v,%d)", refid, pos, reverse, len(k.body))
}

// Return -1, 0, 1 if k0 < k1, k0==k1, k0 > k1, respectively.
func (k sortEntry) compare(other sortEntry) int {
	if k.coord < other.coord {
		return -1
	}
	if k.coord > other.coord {
		return 1
	}
	return bytes.Compare(k.body, other.body)
}

func coordFromRecord(rec *sam.Record) (key recCoord) {
	// This is the same as compareCoordinatesAndStrand function in sambamba.
	if rec.Ref == nil {
		key = unmappedCoord
	} else {
		key = (recCoord(rec.Ref.ID()) << 33) | recCoord(rec.Pos)<<1
	}
	if (rec.Flags & sam.Reverse) != 0 {
		key |= recCoord(1)
	}
	return
}

func (r recCoord) String() string {
	refid, pos, reverse := parseCoord(r)
	return fmt.Sprintf("(%d,%d,%v)", refid, pos, reverse)
}

func makeCoord(refid, pos int, reverse bool) (key recCoord) {
	// This is the same as compareCoordinatesAndStrand function in sambamba.
	if refid < 0 {
		key = unmappedCoord
	} else {
		key = (recCoord(refid) << 33) | recCoord(pos)<<1
	}
	if reverse {
		key |= recCoord(1)
	}
	return
}

func parseCoord(coord recCoord) (refid, pos int, reverse bool) {
	if (coord & unmappedCoord) == unmappedCoord {
		refid = -1
		pos = -1
	} else {
		refid = int(int32(coord >> 33))
		pos = int(int32((coord & 0x1ffffffff) >> 1))
	}
	if (coord & 1) != 0 {
		reverse = true
	}
	return
}

// A key that's larger than any valid key.
const invalidCoord recCoord = 0xfffffffffffffffe
const infinityCoord = invalidCoord

// A key used for all unmapped reads. Corresponds to (refid,pos)=(-1,-1)
const unmappedCoord recCoord = 0x7ffffffffffffffe

// Sorter sorts list of sam.Records and produces a sortshard file in
// "outPath". SortedShardsToBAM can be later used to merge multiple sorted shard
// files into a BAM file. "header" must contain all the references used by
// records to be added later.
//
// Sorter orders records in the following way:
//
// - Increasing reference sequence IDs, then
// - increasing alignment positions, then
// - sorts a forward read before a reverse read.
// - All else equal, sorts records the order of appearance in the input (i.e., stable sort)
//
// These criteria are the same as "samtool sort" and "sambamba sort".
//
// Example:
//   sorter := NewSorter("tmp0.sort", header)
//   for _, rec := range recordlist {
//     sorter.AddRecord(rec)
//   }
//   err := sorter.Close()
//
//   .. Similarly, produce tmp1.sort, .., tmpN.sort, possibly on
//   .. different processes or machines ..
//
//   // Merge all the sorted shards into one BAM file.
//   err := SortedShardsToBAM([]string{"tmp0.sort",..."tmpN.sort"}, "foo.bam")
type Sorter struct {
	options       SortOptions
	outPath       string
	header        *sam.Header
	smallPool     *sync.Pool          // used to serialize a single sam.Record.
	sortBlockPool *sortShardBlockPool // used to reuse buffers.
	totalRecords  uint32
	recs          []sortEntry
	err           errors.Once
	bgSorterCh    chan sortBatch

	wg     sync.WaitGroup
	mu     sync.Mutex
	shards []string // pathnames of temp sortshard files.
}

// A thin wrapper around sortShardReader to read a shard file and do
// reference-ID translations.
type mergeLeaf struct {
	// Index is a number (0,1,2..) arbitrarily assigned to distinguish mergeLeafs
	// that are merged into one destination.
	seq    int
	name   string // the path of shard file; for logging only.
	reader *sortShardReader
	done   bool // reader.scan() returned false?
	err    *errors.Once
}

func newMergeLeaf(seq int, reader *sortShardReader,
	pool *sortShardBlockPool, errorReporter *errors.Once) *mergeLeaf {
	leaf := mergeLeaf{
		seq:    seq,
		name:   reader.path,
		reader: reader,
		err:    errorReporter,
	}
	if !leaf.reader.scan() {
		return nil
	}
	return &leaf
}

func (l *mergeLeaf) Compare(c1 llrb.Comparable) int {
	l1 := c1.(*mergeLeaf)
	k0 := l.reader.key()
	k1 := l1.reader.key()
	if c := k0.compare(k1); c != 0 {
		return c
	}
	return l.seq - l1.seq
}

// Merge multiple sortShards into a new sortShard file.
func (s *Sorter) mergeShards(paths []string, header *sam.Header, outPath string) {
	ctx := vcontext.Background()
	out, err := file.Create(ctx, outPath)
	if err != nil {
		s.err.Set(err)
		return
	}
	shardReaders := make([]*sortShardReader, len(paths))
	for i, path := range paths {
		shardReaders[i] = newSortShardReader(path, s.sortBlockPool, &s.err)
	}
	writer := newSortShardWriter(out.Writer(ctx), !s.options.NoCompressTmpFiles, true, header,
		s.sortBlockPool, &s.err)
	callback := func(key sortEntry) bool {
		writer.add(key)
		return true
	}
	internalMergeShards(shardReaders, callback, s.sortBlockPool, &s.err)
	writer.finish()
	s.err.Set(out.Close(ctx))
}

// Merge sortShards. headerCallback is called once for the merged sam.Header,
// then readCallback is called sequentially for each record in sort order.  If
// readCallback returns false, this function exits immediately.
func internalMergeShards(
	shards []*sortShardReader,
	readCallback func(key sortEntry) bool,
	pool *sortShardBlockPool,
	errReporter *errors.Once) {
	// Sort all the inputs using a binary tree. This should be faster than
	// binary heap or tournament tree. The hope is that the child at the top
	// of the tree will stay at the top for many records. If that hope
	// holds, then tree will can maintain the sorted order in amortized O(1)
	// time, whereas heap always costs O(log(len(outCh)).
	leafs := llrb.Tree{}

	// Create a one-level tree.
	for i, shard := range shards {
		if c := newMergeLeaf(i, shard, pool, errReporter); c != nil {
			vlog.VI(1).Infof("Leaf %v created", c.name)
			leafs.Insert(c)
		}
	}
	vlog.VI(1).Infof("Merging %d shards, %d leafs active", len(shards), leafs.Len())

	// Do N-way merge.  readCallback will be called with increasing list of
	// records.
	done := false
	for !done && leafs.Len() > 0 {
		nthiter := 0
		// top is the smallest child. We read from top.
		// next is the 2nd smallest child, or nil if tree is the only
		// child in the tree.
		var top, next *mergeLeaf
		leafs.Do(func(item llrb.Comparable) bool {
			nthiter++
			switch nthiter {
			case 1:
				top = item.(*mergeLeaf)
				return false
			case 2:
				next = item.(*mergeLeaf)
				return true
			default:
				vlog.Fatal(nthiter)
				return false
			}
		})
		// Read records from top, until it becomes larger than next.
		for {
			if !readCallback(top.reader.key()) {
				done = true
				break
			}
			top.done = !top.reader.scan()
			if top.done || (next != nil && next.reader.key().compare(top.reader.key()) < 0) {
				break
			}
		}
		// Move top into the proper place in the tree.
		lenBefore := leafs.Len()
		leafs.DeleteMin()
		if !top.done {
			leafs.Insert(top)
			if lenAfter := leafs.Len(); lenBefore != lenAfter {
				vlog.Fatalf("Leaf size decreased from %d -> %d", lenBefore, lenAfter)
			}
			continue
		}
	}
	for _, shard := range shards {
		shard.drain()
	}
}

// NewSorter creates a Sorter object.
func NewSorter(outPath string, header *sam.Header, optList ...SortOptions) *Sorter {
	options := SortOptions{}
	if len(optList) > 0 {
		if len(optList) > 1 {
			vlog.Fatalf("More than options specified: %v", optList)
		}
		options = optList[0]
	}
	if options.ShardIndex == 0 {
		hash := sha256.Sum224([]byte(outPath))
		options.ShardIndex = binary.LittleEndian.Uint32(hash[:])
	}
	if options.SortBatchSize <= 0 {
		options.SortBatchSize = DefaultSortBatchSize
	}
	if options.Parallelism <= 0 {
		options.Parallelism = DefaultParallelism
	}
	vlog.VI(1).Infof("New Sorter: %v, %+v", outPath, options)
	sorter := &Sorter{
		options:       options,
		outPath:       outPath,
		header:        header,
		smallPool:     &sync.Pool{New: func() interface{} { return bytes.Buffer{} }},
		sortBlockPool: newSortShardBlockPool(),
		bgSorterCh:    make(chan sortBatch, options.Parallelism),
	}
	for i := 0; i < options.Parallelism; i++ {
		sorter.wg.Add(1)
		go func() {
			for batch := range sorter.bgSorterCh {
				path := sorter.sortRecords(batch.recs, batch.sortTieBreaker)
				sorter.mu.Lock()
				sorter.shards = append(sorter.shards, path)
				sorter.mu.Unlock()
			}
			sorter.wg.Done()
		}()
	}
	return sorter
}

// AddRecord adds a record to the sorter. The sorter takes ownership of
// "rec". The caller shall not read or write "rec" after the call.
func (s *Sorter) AddRecord(rec *sam.Record) {
	s.totalRecords++
	var buf bytes.Buffer
	err := bam.Marshal(rec, &buf)
	if err != nil {
		s.err.Set(err)
		return
	}
	s.recs = append(s.recs, sortEntry{coordFromRecord(rec), buf.Bytes()})
	if len(s.recs) >= s.options.SortBatchSize {
		s.startGenerateSortShard()
	}
}

func (s *Sorter) startGenerateSortShard() {
	s.bgSorterCh <- sortBatch{
		recs:           s.recs,
		sortTieBreaker: (uint64(s.options.ShardIndex) << 32) + uint64(s.totalRecords),
	}
	s.recs = nil
}

func (s *Sorter) sortRecords(records []sortEntry, sortTieBreaker uint64) string {
	vlog.VI(1).Infof("Sorting %d records, tiebreaker %x", len(records), sortTieBreaker)
	temp, err := ioutil.TempFile(s.options.TmpDir, "bamsort")
	if err != nil {
		s.err.Set(err)
		return ""
	}
	sort.SliceStable(records, func(i, j int) bool {
		return records[i].compare(records[j]) < 0
	})
	writer := newSortShardWriter(temp, !s.options.NoCompressTmpFiles, false, nil, s.sortBlockPool, &s.err)
	for _, key := range records {
		writer.add(key)
		// Note: don't call sam.PutInFreePool since rec may be produced by the
		// native biogo code.
	}
	writer.finish()
	s.err.Set(temp.Close())
	return temp.Name()
}

// Close must be called after adding all the records. It blocks the caller until
// the shard file is generated. After Close, Sorter becomes invalid.
func (s *Sorter) Close() error {
	if len(s.recs) > 0 || s.totalRecords == 0 {
		// Note: when totalRecords==0, there's no record to write but we still want
		// to create an empty shard file.
		s.startGenerateSortShard()
	}
	close(s.bgSorterCh)
	s.wg.Wait()
	if s.err.Err() == nil {
		s.mergeShards(s.shards, s.header, s.outPath)
	}
	for _, path := range s.shards {
		if err := os.Remove(path); err != nil {
			vlog.Errorf("sort %v: failed to remove sorter tmp file: %v (%v)", path, err, s.err.Err())
		}
	}
	return s.err.Err()
}

func mergeHeader(shards []*sortShardReader) (*sam.Header, error) {
	shardHeaders := make([]*sam.Header, len(shards))
	for i, shard := range shards {
		if len(shard.index.EncodedBamHeader) == 0 {
			return nil, fmt.Errorf("%s: sortshard index is missing sam.Header", shard.path)
		}
		var err error
		if shardHeaders[i], err = gbam.UnmarshalHeader(shard.index.EncodedBamHeader); err != nil {
			return nil, err
		}
	}

	header, refTranslations, err := sam.MergeHeaders(shardHeaders)
	if err != nil {
		return nil, err
	}

	// Make sure that all the refTranslations are identity translations.
	// Otherwise, the sortkey won't order records properly.
	for _, t := range refTranslations {
		for i, ref := range t {
			if ref.ID() != i {
				return nil, fmt.Errorf("cannot merge sortshards with different SAM headers")
			}
		}
	}
	header.SortOrder = sam.Coordinate
	header.GroupOrder = shardHeaders[0].GroupOrder
	return header, nil
}

// BAMFromSortShards merges a set of sortshard files into a single BAM file.
func BAMFromSortShards(paths []string, bamPath string) error {
	if len(paths) == 0 {
		return fmt.Errorf("no shards to merge")
	}
	errReporter := errors.Once{}
	pool := newSortShardBlockPool()
	shardReaders := make([]*sortShardReader, len(paths))

	for i, path := range paths {
		shardReaders[i] = newSortShardReader(path, pool, &errReporter)
	}

	if err := errReporter.Err(); err != nil {
		// TODO(saito) Close all shard readers.
		return err
	}

	mergedHeader, err := mergeHeader(shardReaders)
	if err != nil {
		// TODO(saito) Close all shard readers.
		return err
	}

	ctx := vcontext.Background()
	out, err := file.Create(ctx, bamPath)
	if err != nil {
		// TODO(saito) Close all shard readers.
		return err
	}
	gzip := bgzf.NewWriter(out.Writer(ctx), runtime.NumCPU()*4)
	writeBytes := func(bytes []byte) {
		_, e := gzip.Write(bytes)
		errReporter.Set(e)
	}

	// Write the BAM header.
	var buf bytes.Buffer
	if err := mergedHeader.EncodeBinary(&buf); err != nil {
		return err
	}
	writeBytes(buf.Bytes())

	// Write the BAM records.
	readCallback := func(key sortEntry) bool {
		writeBytes(key.body)
		return true
	}
	internalMergeShards(shardReaders, readCallback, pool, &errReporter)
	errReporter.Set(gzip.Close())
	errReporter.Set(out.Close(ctx))
	return errReporter.Err()
}
