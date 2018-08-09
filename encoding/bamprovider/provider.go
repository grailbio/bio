package bamprovider

import (
	"strings"

	"github.com/biogo/hts/sam"
	"github.com/grailbio/base/vcontext"
	gbam "github.com/grailbio/bio/encoding/bam"
	"github.com/grailbio/bio/encoding/pam"
	"github.com/grailbio/bio/encoding/pam/pamutil"
	"v.io/x/lib/vlog"
)

// ProviderOpts defines options for NewProvider.
type ProviderOpts struct {
	// Index specifies the name of the BAM inde file. This field is meaningful
	// only for BAM files. If Index=="", it defaults to path + ".bai".
	Index string

	// DropFields causes the listed fields not to be filled in sam.Record. This
	// option is recognized only by the PAM reader.
	DropFields []gbam.FieldType
}

// ShardingStrategy defines algorithms used by Provider.GenerateShards.
type ShardingStrategy int

const (
	// Automatic picks some good strategy. In practice, it means ByteBased for
	// PAM, PositionBased for BAM.
	Automatic ShardingStrategy = iota
	// ByteBased strategy partitions the file so that each shard has roughly equal
	// number of bytes.
	ByteBased
	// PositionBased strategy partitions the file so that each shard covers a
	// genomic coordinate range ([<startref,startpos>, <limitref,limitpos>) with
	// uniform width - i.e., value of (limitpos - startpos) is uniform across
	// shards.
	PositionBased
)

// GenerateShardsOpts defines behavior of Provider.GenerateShards.
type GenerateShardsOpts struct {
	// Strategy defines sharding strategy.
	Strategy ShardingStrategy

	Padding int
	// IncludeUnmapped causes GenerateShards() to produce shards for the
	// unmapped && mate-unmapped reads.
	IncludeUnmapped bool

	// SplitUnmappedCoords allows GenerateShards() to split unmapped
	// reads into multiple shards. Setting this flag true will cause shard
	// size to be more even, but the caller must be able to handle split
	// unmapped reads.
	SplitUnmappedCoords bool

	// SplitMappedCoords allows GenerateShards() to split mapped reads of
	// the same <refid, alignment position> into multiple shards. Setting
	// this flag true will cause shard size to be more even, but the caller
	// must be able to handle split reads.
	SplitMappedCoords bool

	// AlwaysSplitMappedAndUnmappedCoords causes GenerateShard always to split
	// shards at the boundary of mapped and unmapped reads.
	AlwaysSplitMappedAndUnmappedCoords bool

	// BytesPerShard is the target shard size, in bytes. This is consulted only in
	// ByteBased sharding strategy.
	BytesPerShard int64

	// NumShards is the target shard count. It is consulted by the ByteBased sharding
	// strategy, and is ignored if BytesPerShard is set.
	NumShards int

	// MinBasesPerShard defines the nimimum number of bases in each shard. This is
	// consulted only in ByteBased sharding strategy.
	MinBasesPerShard int
}

// Provider allows reading BAM or PAM file in parallel. Thread safe.
type Provider interface {
	// GetHeader returns the header for the provided BAM data.  The callee
	// must not modify the returned header object.
	//
	// REQUIRES: Close has not been called.
	GetHeader() (*sam.Header, error)

	// GenerateShards prepares for parallel reading of genomic data.
	//
	// The Shards split the BAM data from the given provider into
	// contiguous, non-overlapping genomic intervals. A SAM record is
	// associated with a shard if its alignment start position is within the
	// given padding distance of the shard. This means reads near shard
	// boundaries may be associated with more than one shard.
	//
	// Use NewIterator to read records in a shard.
	//
	// REQUIRES: Close has not been called.
	GenerateShards(opts GenerateShardsOpts) ([]gbam.Shard, error)

	// GetFileShards describes how records are split into files. For BAM, this
	// function just returns a UniversalShard since all records are in one
	// file. For PAM, this function returns one entry per fileshard.
	//
	// REQUIRES: Close has not been called.
	GetFileShards() ([]gbam.Shard, error)

	// NewIterator returns an iterator over record contained in the shard.  The
	// "shard" parameter is usually produced by GenerateShards, but the caller may
	// also manually construct it.
	//
	// REQUIRES: Close has not been called.
	NewIterator(shard gbam.Shard) Iterator

	// Close must be called exactly once. It returns any error encountered
	// by the provider, or any iterator created by the provider.
	//
	// REQUIRES: All the iterators created by NewIterator have been closed.
	Close() error
}

// Iterator iterates over sam.Records in a particular genomic range, in
// coordinate order. Thread compatible.
type Iterator interface {
	// Scan returns where there are any records remaining in the iterator,
	// and if so, advances the iterator to the next record. If the iterator
	// reaches the end of its range, Scan() returns false.  If an error
	// occurs, Scan() returns false and the error can be retrieved by
	// calling Error().
	//
	// Scan and Record always yield records in the ascending coordinate
	// (refid,position) order.
	//
	// REQUIRES: Close has not been called.
	Scan() bool

	// Record returns the current record in the iterator. This must be
	// called only after a call to Scan() returns true.
	//
	// REQUIRES: Close has not been called.
	Record() *sam.Record

	// Err returns the error encoutered during iteration, or nil if no error
	// occurred.  An io.EOF error will be translated to nil.
	Err() error

	// Close must be called exactly once. It returns the value of Err().
	Close() error
}

// FileType represents the type of a BAM-like file.
type FileType int

const (
	// Unknown is a sentinel.
	Unknown FileType = iota
	// BAM file
	BAM
	// PAM file
	PAM
)

// ParseFileType parses the file type string. "bam" returns bamprovider.BAM, for
// example. On error, it returns Unknown.
func ParseFileType(name string) FileType {
	switch name {
	case "bam":
		return BAM
	case "pam":
		return PAM
	default:
		return Unknown
	}
}

// GuessFileType returns the file type from the pathname and/or
// contents. Returns Unknown on error.
func GuessFileType(path string) FileType {
	if strings.HasSuffix(path, ".bam") {
		return BAM
	}
	if strings.Contains(path, ".pam") {
		return PAM
	}
	ctx := vcontext.Background()
	if _, err := pamutil.ListIndexes(ctx, path); err == nil {
		return PAM
	}
	vlog.VI(1).Infof("%v: could not detect file type.", path)
	return Unknown
}

func mergeOpts(optList []ProviderOpts) ProviderOpts {
	opts := ProviderOpts{}
	for _, o := range optList {
		if o.Index != "" {
			opts.Index = o.Index
		}
		opts.DropFields = append(opts.DropFields, o.DropFields...)
	}
	return opts
}

// NewProvider creates a Provider object that can handle BAM or PAM file of
// "path". The file type is autodetected from the path.
func NewProvider(path string, optList ...ProviderOpts) Provider {
	opts := mergeOpts(optList)
	switch GuessFileType(path) {
	case BAM, Unknown:
		return &BAMProvider{Path: path, Index: opts.Index}
	case PAM:
		return &PAMProvider{Path: path, Opts: pam.ReadOpts{DropFields: opts.DropFields}}
	}
	panic("shouldn't reach here")
}
