package converter

// Utility for converting BAM to PAM.

import (
	"fmt"
	"math"
	"runtime"
	"sort"
	"sync"
	"sync/atomic"

	"github.com/biogo/hts/bam"
	"github.com/biogo/hts/bgzf"
	"github.com/biogo/hts/sam"
	"github.com/grailbio/base/errors"
	"github.com/grailbio/base/file"
	"github.com/grailbio/base/traverse"
	"github.com/grailbio/base/vcontext"
	"github.com/grailbio/bio/biopb"
	gbam "github.com/grailbio/bio/encoding/bam"
	"github.com/grailbio/bio/encoding/bamprovider"
	"github.com/grailbio/bio/encoding/pam"
	"github.com/grailbio/bio/encoding/pam/pamutil"
	"github.com/klauspost/compress/gzip"
	"v.io/x/lib/vlog"
)

type bamShardBound struct {
	rec biopb.Coord
	off bgzf.Offset
}

// Given the BAM file, find shard boundaries such that each shard is roughly
// bytesPerShard bytes long.
func generateShardBoundaries(bamPath, baiPath string, bytesPerShard int64) ([]bamShardBound, error) {
	ctx := vcontext.Background()
	bamIn, err := file.Open(ctx, bamPath)
	if err != nil {
		return nil, err
	}
	defer bamIn.Close(ctx)
	bamr, err := bam.NewReader(bamIn.Reader(ctx), 1)
	if err != nil {
		return nil, err
	}

	indexIn, err := file.Open(ctx, baiPath)
	if err != nil {
		return nil, err
	}
	defer indexIn.Close(ctx)
	index, err := gbam.ReadIndex(indexIn.Reader(ctx))
	if err != nil {
		return nil, err
	}

	// Enumerate bgzf block boundaries from the index. The size of blocks
	// may be uneven.
	allChunks := []bgzf.Offset{} // guarded by mu
	for _, chunks := range index.AllOffsets() {
		allChunks = append(allChunks, chunks...)
	}
	if len(allChunks) <= 0 {
		return nil, fmt.Errorf("%v: no chunk found in the index", bamPath)
	}
	sort.SliceStable(allChunks, func(i, j int) bool {
		c0 := allChunks[i]
		c1 := allChunks[j]
		if c0.File != c1.File {
			return c0.File < c1.File
		}
		return c0.Block < c1.Block
	})

	// Group the bgzf blocks so that each group is roughly bytesPerShard
	// long.
	bamStat, err := bamIn.Stat(ctx)
	if err != nil {
		return nil, err
	}
	bounds := []bamShardBound{
		bamShardBound{
			rec: biopb.Coord{0, 0, 0},
			off: bgzf.Offset{0, 0},
		},
	}
	if bytesPerShard <= 0 {
		return bounds, nil
	}
	nextGoalOff := bytesPerShard
	for _, chunk := range allChunks[1:] {
		if chunk.File >= nextGoalOff {
			coord, err := gbam.GetCoordAtOffset(bamr, chunk)
			if err != nil {
				vlog.Panic(err)
				continue
			}
			if len(bounds) > 0 {
				cmp := coord.Compare(bounds[len(bounds)-1].rec)
				if cmp < 0 {
					vlog.Panicf("BAM record not sorted properly: %+v %+v", coord, bounds)
				}
				if cmp == 0 {
					continue
				}
			}
			bounds = append(bounds, bamShardBound{rec: coord, off: chunk})
			nextGoalOff = chunk.File + bytesPerShard
		}
	}
	vlog.Infof("Found %d chunk bounds (file size: %d, goal shard size: %d)", len(bounds), bamStat.Size(), bytesPerShard)
	return bounds, nil
}

// convertShard converts range [startShard, limitShard) of the BAM file to
// PAM. Returns the number of sam.Records converted.
func convertShard(opts pam.WriteOpts, pamPath string, in bamprovider.Provider, startShard, limitShard bamShardBound) (int64, error) {
	header, e := in.GetHeader()
	if e != nil {
		return 0, e
	}

	getRef := func(id int32) *sam.Reference {
		if id == biopb.UnmappedRefID {
			return nil
		}
		return header.Refs()[id]
	}

	iter := in.NewIterator(gbam.Shard{
		StartRef: getRef(startShard.rec.RefId),
		Start:    int(startShard.rec.Pos),
		EndRef:   getRef(limitShard.rec.RefId),
		End:      int(limitShard.rec.Pos),
	})
	opts.Range = biopb.CoordRange{Start: startShard.rec, Limit: limitShard.rec}
	w := pam.NewWriter(opts, header, pamPath)
	err := errors.Once{}
	var nRecs int64
	for iter.Scan() {
		rec := iter.Record()
		w.Write(rec)
		if w.Err() != nil {
			break
		}
		nRecs++
		gbam.PutInFreePool(gbam.CastDown(rec))
	}
	err.Set(iter.Close())
	err.Set(w.Close())
	vlog.Infof("%v: Finished converting [%+v, %+v) with %d recs read: %v",
		pamPath, startShard, limitShard, nRecs, err.Err())
	return nRecs, err.Err()
}

// ConvertToPAM copies BAM to PAM.
//
// bytesPerShard is specified in term of the input BAM file. For example, if you
// have a 100GB BAM file and set bytesPerShard to 20GB, this function will
// create roughly five (=100/20) PAM shards.
func ConvertToPAM(opts pam.WriteOpts, pamPath, bamPath, baiPath string, bytesPerShard int64) error {
	if baiPath == "" {
		baiPath = bamPath + ".bai"
	}
	if pamPath == "" {
		return fmt.Errorf("Empty pam path")
	}
	if bytesPerShard <= 0 {
		return fmt.Errorf("Negative bytesPerShard: %v", bytesPerShard)
	}
	shards, e := generateShardBoundaries(bamPath, baiPath, bytesPerShard)
	if e != nil {
		return e
	}
	if e := pamutil.ValidateCoordRange(&opts.Range); e != nil {
		return e
	}
	if !opts.Range.EQ(gbam.UniversalRange) {
		return fmt.Errorf("WriteOpts.Range to ConvertFromBAM must be a universal range, but found %+v", opts)
	}
	vlog.Infof("%v: Creating %d shards: %+v", pamPath, len(shards), shards)
	// Delete existing files to avoid mixing up files from multiple generations.
	if e := pamutil.Remove(pamPath); e != nil {
		return e
	}

	var totalRecs int64
	bam := bamprovider.BAMProvider{Path: bamPath, Index: baiPath}
	err := traverse.Each(len(shards), func(i int) error {
		shard := shards[i]
		nextShard := bamShardBound{
			rec: biopb.Coord{biopb.InfinityRefID, biopb.InfinityPos, 0},
			off: bgzf.Offset{math.MaxInt64, math.MaxUint16},
		}
		if i < len(shards)-1 {
			nextShard = shards[i+1]
		}
		nRecs, err := convertShard(opts, pamPath, &bam, shard, nextShard)
		atomic.AddInt64(&totalRecs, nRecs)
		return err
	})
	if e := bam.Close(); e != nil && err == nil {
		err = e
	}
	vlog.Infof("%v: Finished converting, written %d records, error %v", pamPath, totalRecs, err)
	return err
}

type convertRequest struct {
	shardIdx int
	records  []*sam.Record
}

// ConvertToBAM copies "provider" to a BAM file. Existing contents of "bamPath",
// if any, are destroyed.
func ConvertToBAM(bamPath string, provider bamprovider.Provider) error {
	const recordsPerShard = 128 << 10
	parallelism := runtime.NumCPU()

	ctx := vcontext.Background()
	header, e := provider.GetHeader()
	if e != nil {
		return e
	}
	out, e := file.Create(ctx, bamPath)
	if e != nil {
		return e
	}
	w, e := gbam.NewShardedBAMWriter(out.Writer(ctx), gzip.DefaultCompression, parallelism*4, header)
	if e != nil {
		return e
	}

	wg := sync.WaitGroup{}
	reqCh := make(chan convertRequest, parallelism)
	var err errors.Once
	vlog.Infof("Creating %d threads", parallelism)
	for wi := 0; wi < parallelism; wi++ {
		wg.Add(1)
		go func() {
			c := w.GetCompressor()
			for req := range reqCh {
				if e := c.StartShard(req.shardIdx); e != nil {
					err.Set(e)
					break
				}
				for _, r := range req.records {
					c.AddRecord(r)
					gbam.PutInFreePool(gbam.CastDown(r))
				}
				err.Set(c.CloseShard())
			}
			wg.Done()
		}()
	}

	iter := provider.NewIterator(gbam.UniversalShard(header))
	req := convertRequest{
		records:  make([]*sam.Record, 0, recordsPerShard),
		shardIdx: 0,
	}
	for iter.Scan() {
		req.records = append(req.records, iter.Record())
		if len(req.records) >= recordsPerShard {
			reqCh <- req
			req.records = make([]*sam.Record, 0, recordsPerShard)
			req.shardIdx++
		}
	}
	if len(req.records) > 0 {
		reqCh <- req
	}
	close(reqCh)

	err.Set(iter.Close())
	wg.Wait()
	err.Set(w.Close())
	err.Set(out.Close(ctx))
	return err.Err()
}
