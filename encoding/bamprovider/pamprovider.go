package bamprovider

import (
	"sync"

	"github.com/biogo/hts/sam"
	"github.com/grailbio/base/errorreporter"
	"github.com/grailbio/bio/biopb"
	gbam "github.com/grailbio/bio/encoding/bam"
	"github.com/grailbio/bio/encoding/pam"
	"github.com/pkg/errors"
)

// PAMProvider reads PAM files.  The path can be S3 URLs, in which case the data
// will be read from S3. Otherwise the data will be read from the local
// filesystem.
type PAMProvider struct {
	// Path prefix. Must be nonempty.
	Path string
	// Opts is passed to pam.NewReader.
	Opts pam.ReadOpts
	err  errorreporter.T

	mu      sync.Mutex
	header  *sam.Header    // extracted from <dir>/<range>.index.
	indexes []pam.FileInfo // files found in the pam directory.
}

// pamIterator implements the Iterator interface.
type pamIterator struct {
	provider *PAMProvider
	reader   *pam.Reader
}

// GetHeader implements the Provider interface.
func (p *PAMProvider) GetHeader() (*sam.Header, error) {
	p.mu.Lock()
	defer p.mu.Unlock()
	if p.header != nil {
		return p.header, nil
	}
	if len(p.indexes) == 0 {
		indexes, err := pam.ListIndexes(p.Path)
		if err != nil {
			p.err.Set(err)
			return nil, err
		}
		if len(indexes) == 0 {
			panic(p)
		}
		p.indexes = indexes
	}
	index, err := pam.ReadShardIndex(p.Path, p.indexes[0].Range)
	if err != nil {
		p.err.Set(err)
		return nil, err
	}
	p.header, err = gbam.UnmarshalHeader(index.EncodedBamHeader)
	if err != nil {
		p.err.Set(err)
		return nil, err
	}
	return p.header, nil
}

// Close implements the Provider interface.
func (p *PAMProvider) Close() error {
	return p.err.Err()
}

// GenerateShards implements the Provider interface.
func (p *PAMProvider) GenerateShards(opts GenerateShardsOpts) ([]gbam.Shard, error) {
	if opts.Strategy != Automatic && opts.Strategy != ByteBased {
		return nil, errors.Errorf("GenerateShards: strategy %v not supported", opts.Strategy)
	}
	header, err := p.GetHeader()
	if err != nil {
		return nil, err
	}
	popts := pam.GenerateReadShardsOpts{
		Range:                              gbam.UniversalRange,
		SplitMappedCoords:                  opts.SplitMappedCoords,
		SplitUnmappedCoords:                opts.SplitUnmappedCoords,
		AlwaysSplitMappedAndUnmappedCoords: opts.AlwaysSplitMappedAndUnmappedCoords,
		BytesPerShard:                      opts.BytesPerShard,
	}
	if !opts.IncludeUnmapped {
		popts.Range = gbam.MappedRange
	}
	pamShards, err := pam.GenerateReadShards(popts, p.Path)
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
