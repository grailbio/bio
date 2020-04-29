// Copyright 2020 Grail Inc.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
package main

/*
bio-pileup is a variant calling tool which reports the number of reads in a
BAM/PAM supporting each allele at each genomic position.
*/

import (
	"flag"
	"fmt"
	"os"
	"strings"

	"github.com/grailbio/base/grail"
	"github.com/grailbio/base/log"
	"github.com/grailbio/base/vcontext"
	"github.com/grailbio/bio/cmd/bio-pileup/snp"
)

var (
	bedPath      = flag.String("bed", "", "Input BED path; this xor -region required")
	region       = flag.String("region", "", "Restrict pileup computation to the specified region. Format as <contig ID>:<1-based first pos>-<last pos>, <contig ID>:<1-based pos>, or just <contig ID>; this xor -bed required")
	bamIndexPath = flag.String("index", "", "Input BAM index path. Defaults to bampath + .bai")
	clip         = flag.Int("clip", 0, "Number of bases on end of each read to treat as minimum-quality")
	cols         = flag.String("cols", "", "Output TSV column sets. #CHROM/POS/REF(/ALT) are always present. Currently supported optional sets are 'dpref', 'dpalt', 'enddists', 'quals', 'fraglens', 'strands', 'highq', and 'lowq'; default is \"dpref,highq,lowq\"")
	flagExclude  = flag.Int("flag-exclude", 0xf00, "Reads with a FLAG bit intersecting this value are skipped")
	format       = flag.String("format", "tsv", "Output format; 'basestrand-rio', 'basestrand-tsv', 'basestrand-tsv-bgz', 'tsv', and 'tsv-bgz' supported")
	mapq         = flag.Int("mapq", 60, "Reads with MAPQ below this level are skipped")
	maxReadLen   = flag.Int("max-read-len", 500, "Upper bound on individual read length")
	maxReadSpan  = flag.Int("max-read-span", 511, "Upper bound on size of reference-genome region a read maps to")
	minBagDepth  = flag.Int("min-bag-depth", 0, "Lower bound on bag depth (DS aux tag value")
	minBaseQual  = flag.Int("min-base-qual", 0, "Lower bound on base quality in a single read")
	outPrefix    = flag.String("out", "bio-pileup", "Output path prefix")
	parallelism  = flag.Int("parallelism", 0, "Maximum number of simultaneous (local) pileup jobs to launch; 0 = runtime.NumCPU()")
	perStrand    = flag.Bool("per-strand", false, "Generate two pairs of output files, one for each strand")
	removeSq     = flag.Bool("remove-sq", false, "Remove sequencing duplicates (no DL aux tag with value > 1)")
	stitch       = flag.Bool("stitch", false, "Stitch read-pairs")
	tempDir      = flag.String("temp-dir", "", "Directory to write temporary files to (default os.TempDir())")
)

func bioPileupUsage() {
	fmt.Printf("Usage: %s [OPTIONS] {b,p}ampath fapath\n", os.Args[0])
	fmt.Printf("Other options:\n")
	flag.PrintDefaults()
}

func main() {
	flag.Usage = bioPileupUsage
	shutdown := grail.Init()
	defer shutdown()

	allArgs := flag.Args()
	nPositionalArgs := flag.NArg()
	positionalArgs := allArgs[len(allArgs)-nPositionalArgs:]
	if nPositionalArgs != 2 {
		if nPositionalArgs < 2 {
			log.Fatalf("Missing positional arguments ({b,p}ampath and fapath required); please check flag syntax: '%s'", strings.Join(positionalArgs, " "))
		} else {
			log.Fatalf("Too many positional arguments (only {b,p}ampath and fapath expected); please check flag syntax: '%s'", strings.Join(positionalArgs, " "))
		}
	}
	ctx := vcontext.Background()
	opts := snp.Opts{
		BedPath:      *bedPath,
		Region:       *region,
		BamIndexPath: *bamIndexPath,
		Clip:         *clip,
		Cols:         *cols,
		FlagExclude:  *flagExclude,
		Mapq:         *mapq,
		MaxReadLen:   *maxReadLen,
		MaxReadSpan:  *maxReadSpan,
		MinBagDepth:  *minBagDepth,
		MinBaseQual:  *minBaseQual,
		Parallelism:  *parallelism,
		PerStrand:    *perStrand,
		RemoveSq:     *removeSq,
		Stitch:       *stitch,
		TempDir:      *tempDir,
	}
	if err := snp.Pileup(ctx, positionalArgs[0], positionalArgs[1], *format, *outPrefix, &opts); err != nil {
		log.Panicf("%v", err)
	}
	log.Debug.Printf("exiting")
}
