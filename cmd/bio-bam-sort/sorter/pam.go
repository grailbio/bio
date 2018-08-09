package sorter

import (
	"fmt"
	"math"
	"sort"
	"sync"

	"github.com/biogo/hts/sam"
	"github.com/grailbio/base/errorreporter"
	"github.com/grailbio/bio/biopb"
	grailbam "github.com/grailbio/bio/encoding/bam"
	"github.com/grailbio/bio/encoding/pam"
	"github.com/grailbio/bio/encoding/pam/pamutil"
	"v.io/x/lib/vlog"
)

// Compute good PAM rowshard boundaries.  "shards" is the union of index blocks
// in all the shard files being merged. recordsPerShard is the goal # of reads
// to store in each rowshard.
//
// The resulting "bounds" are a sequence of sortKeys. The element are in an
// increasing key order. The caller should create len(bounds)+1 PAM rowshards,
// with the boundary between K'th and K+1'th shard at bounds[k].
func computePAMShardBounds(
	shards [][]biopb.SortShardBlockIndex,
	recordsPerShard int64) (bounds []recCoord) {
	clearReverseFlag := func(coord recCoord) recCoord {
		ref, pos, _ := parseCoord(coord)
		return makeCoord(ref, pos, false)
	}

	allKeys := []sortKey{}
	totalRecords := int64(0)
	totalBlocks := int64(0)
	// Sort the starting keys for all the blocks in increasing order.
	for _, shard := range shards {
		for _, block := range shard {
			allKeys = append(allKeys, sortKey{recCoord(block.StartKey), block.StartSeq})
			totalRecords += int64(block.NumRecords)
		}
		totalBlocks += int64(len(shard))
	}
	if totalRecords == 0 {
		return nil
	}
	sort.SliceStable(allKeys, func(i, j int) bool {
		return allKeys[i].compare(allKeys[j]) < 0
	})

	// Assume that the record size is uniform.
	nextGoalRecords := float64(recordsPerShard)
	recordsPerBlock := float64(totalRecords) / float64(totalBlocks)
	curNumRecords := float64(0)
	for _, key := range allKeys[1:] {
		curNumRecords += recordsPerBlock
		if curNumRecords < nextGoalRecords {
			continue
		}
		coord := clearReverseFlag(key.coord)
		// Make sure that new boundary (ref,pos) is different from the previous one.
		if len(bounds) > 0 && bounds[len(bounds)-1] == coord {
			continue
		}
		bounds = append(bounds, coord)
		nextGoalRecords = curNumRecords + float64(recordsPerShard)
	}
	return
}

// Find the location of the first recordio block to read, assuming one wants
// read records >= sortKey{coord,0}
func startFileOffset(blocks []biopb.SortShardBlockIndex, coord recCoord) int64 {
	key := sortKey{coord: coord, seq: 0}
	n := sort.Search(len(blocks), func(i int) bool {
		b := sortKey{recCoord(blocks[i].StartKey), blocks[i].StartSeq}
		return b.compare(key) >= 0
	})
	// blocks[n] is the first block >= "coord"
	if n >= len(blocks) {
		return int64(blocks[len(blocks)-1].FileOffset)
	}
	// If block[n] starts exactly at "key", then block[n] is the starting point.
	b := sortKey{recCoord(blocks[n].StartKey), blocks[n].StartSeq}
	if comp := b.compare(key); comp == 0 {
		return int64(blocks[n].FileOffset)
	} else if comp < 0 {
		vlog.Fatalf("blocks: %v coord %v", blocks[n], coord)
	}
	// Else, we need to seek one block back.
	if n >= 1 {
		return int64(blocks[n-1].FileOffset)
	}
	return 0
}

// Find location of the first recordio block that doesn't need to be read,
// assuming one wants read records < sortKey{coord,0}.
func limitFileOffset(blocks []biopb.SortShardBlockIndex, coord recCoord) int64 {
	key := sortKey{coord: coord, seq: 0}
	n := sort.Search(len(blocks), func(i int) bool {
		b := sortKey{recCoord(blocks[i].StartKey), blocks[i].StartSeq}
		return b.compare(key) >= 0
	})
	// blocks[n:] are the blocks that we can exclude.
	if n >= len(blocks) {
		return math.MaxInt64
	}
	return int64(blocks[n].FileOffset)
}

// TODO(saito) Consider unifying sortKey and RecAddr, or at least representing
// sortKey in using RecAddr.
func recCoordToRecAddr(coord recCoord) (r biopb.Coord) {
	switch coord {
	case 0:
		// zero coord
		return
	case unmappedCoord:
		// start of unmapped reads
		r.RefId = biopb.InfinityRefID
		r.Pos = 0
		return
	case invalidCoord:
		// infinity coord
		r.RefId = biopb.InfinityRefID
		r.Pos = math.MaxInt32
		return
	default:
		refid, pos, _ := parseCoord(coord)
		if refid < 0 {
			r = biopb.Coord{biopb.UnmappedRefID, 0, 0}
		} else {
			r = biopb.Coord{int32(refid), int32(pos), 0}
		}
		return
	}
}

// Generate one PAM rowshard that stores reads in range [start, limit).
func generatePAMShard(readers []*sortShardReader,
	path string,
	header *sam.Header,
	start, limit recCoord,
	pool *sortShardBlockPool,
	errReporter *errorreporter.T) {
	opts := pam.WriteOpts{
		Range: biopb.CoordRange{recCoordToRecAddr(start), recCoordToRecAddr(limit)},
	}
	vlog.VI(1).Infof("%v: Generating PAM shard %+v", path, opts)
	pamWriter := pam.NewWriter(opts, header, path)
	readCallback := func(key sortKey, data []byte) bool {
		if key.coord >= limit {
			// The rest of the records are not needed.
			return false
		}
		if key.coord >= start {
			// The first 4 bytes of data is the length field. Remove it.
			rec, err := grailbam.Unmarshal(data[4:], header)
			if err != nil {
				errReporter.Set(err)
				return false
			}
			pamWriter.Write(grailbam.CastUp(rec))
			if pamWriter.Err() != nil {
				vlog.Fatalf("ERR: %+v, opts %+v key %+v, limit %+v", pamWriter.Err(), opts, key, limit)
			}
			grailbam.PutInFreePool(rec)
		}
		return true
	}
	internalMergeShards(readers, readCallback, pool, errReporter)
	errReporter.Set(pamWriter.Close())
}

type generatePAMShardRequest struct {
	start, limit recCoord
}

// PAMFromSortShards merges a set of sortshard files into a single PAM file.
// recordsPerShard is the goal # of reads to store in each rowshard.
func PAMFromSortShards(paths []string, pamPath string, recordsPerShard int64, parallelism int) error {
	if len(paths) == 0 {
		return fmt.Errorf("No shards to merge")
	}
	vlog.VI(1).Infof("%v: Generate PAM, #recordspershard=%d", pamPath, recordsPerShard)
	// Delete existing files to avoid mixing up files from multiple generations.
	if err := pamutil.Remove(pamPath); err != nil {
		return err
	}
	errReporter := errorreporter.T{}
	pool := newSortShardBlockPool()
	baseReaders := make([]*sortShardReader, len(paths))

	allBlocks := [][]biopb.SortShardBlockIndex{}
	for i, path := range paths {
		baseReaders[i] = newSortShardReader(path, pool, &errReporter)
		// TODO(saito) baseReaders is leaked.
		allBlocks = append(allBlocks, baseReaders[i].index.Blocks)
	}
	if err := errReporter.Err(); err != nil {
		return err
	}
	mergedHeader, err := mergeHeader(baseReaders)
	if err != nil {
		return err
	}
	pamBounds := computePAMShardBounds(allBlocks, recordsPerShard)

	reqCh := make(chan generatePAMShardRequest, parallelism)
	wg := sync.WaitGroup{}
	for i := 0; i < parallelism; i++ {
		wg.Add(1)
		go func() {
			defer wg.Done()
			for req := range reqCh {
				vlog.Infof("%s: starting generating PAM shard %+v", pamPath, req)
				subReaders := make([]*sortShardReader, len(paths))
				for i, path := range paths {
					opts := shardReaderOpts{
						index:       &baseReaders[i].index,
						startOffset: startFileOffset(baseReaders[i].index.Blocks, req.start),
						limitOffset: limitFileOffset(baseReaders[i].index.Blocks, req.limit),
					}
					subReaders[i] = newSortShardReader(path, pool, &errReporter, opts)
				}
				generatePAMShard(subReaders, pamPath, mergedHeader, req.start, req.limit, pool, &errReporter)
				vlog.Infof("%s: finished generating PAM shard %+v", pamPath, req)
			}
		}()
	}

	for i := 0; i <= len(pamBounds); i++ {
		start := recCoord(0)
		if i > 0 {
			start = pamBounds[i-1]
		}
		limit := infinityCoord
		if i < len(pamBounds) {
			limit = pamBounds[i]
		}
		reqCh <- generatePAMShardRequest{start, limit}
	}
	close(reqCh)
	wg.Wait()
	return errReporter.Err()
}
