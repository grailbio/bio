package bamprovider

import (
	"fmt"
	"sync"

	"github.com/grailbio/base/errors"
	"github.com/grailbio/base/file"
	"github.com/grailbio/base/vcontext"
	"github.com/grailbio/bio/biopb"
	gbam "github.com/grailbio/bio/encoding/bam"
	"github.com/grailbio/bio/encoding/pam"
	"github.com/grailbio/bio/encoding/pam/pamutil"
	"github.com/grailbio/hts/sam"
)

// PAMProvider reads PAM files.  The path can be S3 URLs, in which case the data
// will be read from S3. Otherwise the data will be read from the local
// filesystem.
type PAMProvider struct {
	// Path prefix. Must be nonempty.
	Path string
	// Opts is passed to pam.NewReader.
	Opts pam.ReadOpts
	err  errors.Once

	mu      sync.Mutex
	header  *sam.Header        // extracted from <dir>/<range>.index.
	info    FileInfo           // extracted from <dir>/<range>.index.
	indexes []pamutil.FileInfo // files found in the pam directory.
}

// pamIterator implements the Iterator interface.
type pamIterator struct {
	provider *PAMProvider
	reader   *pam.Reader
}

func (p *PAMProvider) initInfo() {
	p.mu.Lock()
	defer p.mu.Unlock()
	if p.header != nil || p.err.Err() != nil {
		return
	}

	ctx := vcontext.Background()
	if len(p.indexes) == 0 {
		indexes, err := pamutil.ListIndexes(ctx, p.Path)
		if err != nil {
			p.err.Set(err)
			return
		} else if len(indexes) == 0 {
			p.err.Set(fmt.Errorf("pamprovider %v: no pam file found for range %+v", p.Path, p.Opts))
			return
		}
		p.indexes = indexes
	}

	indexPath := pamutil.ShardIndexPath(p.Path, p.indexes[0].Range)
	in, err := file.Open(ctx, indexPath)
	if err != nil {
		p.err.Set(err)
		return
	}
	defer in.Close(ctx) // nolint: errcheck
	info, err := in.Stat(ctx)
	if err != nil {
		p.err.Set(err)
		return
	}
	p.info = FileInfo{ModTime: info.ModTime(), Size: info.Size()}
	index, err := pamutil.ReadShardIndex(ctx, p.Path, p.indexes[0].Range)
	if err != nil {
		p.err.Set(err)
		return
	}
	p.header, err = gbam.UnmarshalHeader(index.EncodedBamHeader)
	if err != nil {
		p.err.Set(err)
		return
	}
}

// FileInfo implements the Provider interface.
func (p *PAMProvider) FileInfo() (FileInfo, error) {
	p.initInfo()
	// p.info is constant after initInfo, so it's ok to read it unlocked.
	return p.info, p.err.Err()
}

// GetHeader implements the Provider interface.
func (p *PAMProvider) GetHeader() (*sam.Header, error) {
	p.initInfo()
	// p.header is constant after initInfo, so it's ok to read it unlocked.
	return p.header, p.err.Err()
}

// Close implements the Provider interface.
func (p *PAMProvider) Close() error {
	return p.err.Err()
}

// GenerateShards implements the Provider interface.
func (p *PAMProvider) GenerateShards(opts GenerateShardsOpts) ([]gbam.Shard, error) {
	if opts.Strategy != Automatic && opts.Strategy != ByteBased {
		return nil, fmt.Errorf("GenerateShards: strategy %v not supported", opts.Strategy)
	}
	if (opts.SplitMappedCoords || opts.SplitUnmappedCoords) && (opts.Padding != 0) {
		// We might want to support this: an operation which doesn't care about
		// distant mates could be parallelized by making each goroutine responsible
		// for the read-pairs where the first read lands inside a given shard
		// slice, and the shards might be more even if SplitMappedCoords is
		// specified.
		// However, this requires the caller to have access to each read's Seq
		// value within the mapped position; as of this writing, we don't expose
		// it.  And it's plausible that practically all padding use cases are best
		// handled without coordinate-splitting.  So just prohibit it unless/until
		// we run into a performance problem that this solves.
		return nil, fmt.Errorf("GenerateShards: nonzero Padding cannot be specified with Split*Coords")
	}
	header, err := p.GetHeader()
	if err != nil {
		return nil, err
	}
	popts := pamutil.GenerateReadShardsOpts{
		Range:                              gbam.UniversalRange,
		SplitMappedCoords:                  opts.SplitMappedCoords,
		SplitUnmappedCoords:                opts.SplitUnmappedCoords,
		AlwaysSplitMappedAndUnmappedCoords: opts.AlwaysSplitMappedAndUnmappedCoords,
		BytesPerShard:                      opts.BytesPerShard,
		NumShards:                          opts.NumShards,
	}
	if !opts.IncludeUnmapped {
		popts.Range = gbam.MappedRange
	}
	ctx := vcontext.Background()
	pamShardIndexes, err := pamutil.ReadIndexes(ctx, p.Path, popts.Range, gbam.FieldNames)
	if err != nil {
		return nil, err
	}
	pamShards, err := pamutil.GenerateReadShards(popts, pamShardIndexes)
	if err != nil {
		return nil, err
	}
	bamShards := make([]gbam.Shard, len(pamShards))
	for index, r := range pamShards {
		bamShards[index] = gbam.CoordRangeToShard(header, r, opts.Padding, index)
	}
	return bamShards, nil
}

// GetFileShards implements the Provider interface.
func (p *PAMProvider) GetFileShards() ([]gbam.Shard, error) {
	header, err := p.GetHeader()
	if err != nil {
		return nil, err
	}
	p.mu.Lock()
	defer p.mu.Unlock()
	if len(p.indexes) == 0 {
		panic(p)
	}
	bamShards := make([]gbam.Shard, len(p.indexes))
	for index, f := range p.indexes {
		bamShards[index] = gbam.CoordRangeToShard(header, f.Range, 0, index)
	}
	return bamShards, nil
}

// NewIterator implements Provider.GetIndexedReader.
func (p *PAMProvider) NewIterator(shard gbam.Shard) Iterator {
	opts := p.Opts
	// This assumes that either padding is zero and/or Split*Coords isn't
	// specified.
	opts.Range.Start = biopb.Coord{int32(shard.StartRef.ID()), int32(shard.PaddedStart()), int32(shard.StartSeq)}
	opts.Range.Limit = biopb.Coord{int32(shard.EndRef.ID()), int32(shard.PaddedEnd()), int32(shard.EndSeq)}
	return &pamIterator{
		provider: p,
		reader:   pam.NewReader(opts, p.Path),
	}
}

func (i *pamIterator) Scan() bool          { return i.reader.Scan() }
func (i *pamIterator) Record() *sam.Record { return i.reader.Record() }
func (i *pamIterator) Err() error          { return i.reader.Err() }

func (i *pamIterator) Close() error {
	err := i.reader.Close()
	if err != nil {
		i.provider.err.Set(err)
	}
	return err
}
