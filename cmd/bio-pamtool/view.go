package main

import (
	"fmt"
	"math"
	"regexp"
	"strconv"
	"strings"

	"github.com/biogo/hts/sam"
	gbam "github.com/grailbio/bio/encoding/bam"
	"github.com/grailbio/bio/encoding/bamprovider"
)

type viewRegion struct {
	// 0-based, half-open interval within refName.
	startRefName       string
	startPos, startSeq int

	limitRefName       string
	limitPos, limitSeq int
}

// Parse flag of form "chr0:beg0-end0,...,chrk:begk-endk".
func parseRegionsFlag(flag string) ([]viewRegion, error) {
	regions := []viewRegion{}
	// samtools-compatible format
	re0 := regexp.MustCompile("([^:]*):(\\d+)-(\\d+)")
	// ref:pos:seq-ref:pos:seq
	re1 := regexp.MustCompile("([^:]*):(\\d+)(:\\d+)?-([^:]*):(\\d+)(:\\d+)?")
	for _, val := range strings.Split(flag, ",") {
		var r viewRegion
		if matches := re0.FindStringSubmatch(val); matches != nil {
			begin, err := strconv.ParseInt(matches[2], 0, 64)
			if err != nil {
				return nil, err
			}
			end, err := strconv.ParseInt(matches[3], 0, 64)
			if err != nil {
				return nil, err
			}
			r = viewRegion{
				startRefName: matches[1],
				startPos:     int(begin - 1), // convert to 0-based closed bound
				startSeq:     0,
				limitRefName: matches[1],
				limitPos:     int(end), // convert to 0-based open bound.
				limitSeq:     0,
			}
		} else if matches := re1.FindStringSubmatch(val); matches != nil {
			parseCoord := func(matches []string) (int, int, error) {
				var pos, seq int64
				var err error
				if pos, err = strconv.ParseInt(matches[0], 0, 64); err != nil {
					return -1, -1, err
				}
				if matches[1] != "" {
					if seq, err = strconv.ParseInt(matches[1][1:], 0, 64); err != nil {
						return -1, -1, err
					}
				}
				return int(pos), int(seq), nil
			}
			r.startRefName = matches[1]
			r.limitRefName = matches[4]
			var err error
			if r.startPos, r.startSeq, err = parseCoord(matches[2:4]); err != nil {
				return nil, err
			}
			if r.limitPos, r.limitSeq, err = parseCoord(matches[5:7]); err != nil {
				return nil, err
			}
		} else {
			return nil, fmt.Errorf("%s: must be of form 'chr:beg-end' or 'chr:beginpos:beginseq-chr:limitpos:limitseq", val)
		}
		regions = append(regions, r)
	}
	return regions, nil
}

func readShard(provider bamprovider.Provider, shard gbam.Shard, filter *filterExpr) error {
	iter := provider.NewIterator(shard)
	for iter.Scan() {
		if filter == nil || evaluateFilterExpr(filter, iter.Record()) {
			s, err := iter.Record().MarshalText()
			if err != nil {
				return err
			}
			fmt.Print(string(s), "\n")
		}
	}
	return iter.Close()
}

func viewAll(provider bamprovider.Provider, filter *filterExpr) error {
	header, err := provider.GetHeader()
	if err != nil {
		return err
	}
	shard := gbam.Shard{
		StartRef: header.Refs()[0],
		EndRef:   nil, // unmapped
		Start:    0,
		End:      math.MaxInt32,
	}
	return readShard(provider, shard, filter)
}

func viewSubregion(provider bamprovider.Provider, region viewRegion, filter *filterExpr) error {
	header, err := provider.GetHeader()
	if err != nil {
		return err
	}
	findRef := func(name string) (*sam.Reference, error) {
		if name == "" || name == "*" {
			return nil, nil
		}
		for _, r := range header.Refs() {
			if r.Name() == name {
				return r, nil
			}
		}
		return nil, fmt.Errorf("Reference %v not found in header", name)
	}
	shard := gbam.Shard{
		Start:    region.startPos,
		StartSeq: region.startSeq,
		End:      region.limitPos,
		EndSeq:   region.limitSeq,
	}
	if shard.StartRef, err = findRef(region.startRefName); err != nil {
		return err
	}
	if shard.EndRef, err = findRef(region.limitRefName); err != nil {
		return err
	}
	return readShard(provider, shard, filter)
}

type viewFlags struct {
	bamIndex   *string
	withHeader *bool
	headerOnly *bool
	regions    *string
	filter     *string
}

// TODO(saito) Currently this function only dumps the index info.  Add feature
// to read data sections too.
func view(flags viewFlags, path string) error {
	regions := []viewRegion{}
	if *flags.regions != "" {
		var err error
		regions, err = parseRegionsFlag(*flags.regions)
		if err != nil {
			return err
		}
	}
	var filter *filterExpr
	if *flags.filter != "" {
		var err error
		filter, err = parseFilterExpr(*flags.filter)
		if err != nil {
			return err
		}
	}
	provider := bamprovider.NewProvider(path, bamprovider.ProviderOpts{Index: *flags.bamIndex})
	if *flags.headerOnly || *flags.withHeader {
		header, err := provider.GetHeader()
		if err != nil {
			return err
		}
		h, err := header.MarshalText()
		if err != nil {
			return err
		}
		fmt.Print(string(h))
		if *flags.headerOnly {
			return nil
		}
	}
	if len(regions) > 0 {
		for _, region := range regions {
			if err := viewSubregion(provider, region, filter); err != nil {
				return err
			}
		}
	} else {
		if err := viewAll(provider, filter); err != nil {
			return err
		}
	}
	return nil
}
