package markduplicates

import (
	"fmt"

	"github.com/grailbio/bio/encoding/bamprovider"
)

func validate(opts *Opts) error {
	if opts.BamFile == "" {
		return fmt.Errorf("you must specify a bam file with --bam")
	}
	if opts.ShardSize <= 0 {
		return fmt.Errorf("shard-size must be non-zero")
	}
	if opts.Padding < 0 {
		return fmt.Errorf("padding must be non-negative")
	}
	if opts.Padding >= opts.ShardSize {
		return fmt.Errorf("padding must be less than shard-size")
	}
	if opts.MinBases <= 0 {
		return fmt.Errorf("min-bases should be positive")
	}
	if opts.IndexFile == "" {
		opts.IndexFile = opts.BamFile + ".bai"
	}
	if len(opts.UmiFile) > 0 && !opts.UseUmis {
		return fmt.Errorf("umi-file is set, but use-umis is false")
	}
	if opts.ScavengeUmis > -1 && !opts.UseUmis {
		return fmt.Errorf("scavenge-umis is set, but use-umis is false")
	}
	if opts.ScavengeUmis > -1 && opts.UmiFile == "" {
		return fmt.Errorf("scavenge-umis is set, but umi-file is empty")
	}
	if bamprovider.ParseFileType(opts.Format) == bamprovider.Unknown {
		return fmt.Errorf("unknown outputformat %s", opts.Format)
	}
	return nil
}
