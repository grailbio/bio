/*
This is the main package for parsegencode. It accepts a gencode GTF, a genome fasta and prints
transcript fasta records to a user-specified (or default) outfile, optionally padding the exons by
any number of bases.

To run parsegencode minimally:

parse_gencode /path/to/gencode_gtf /path/to/genome_fasta

This will create /path/to/genome_fasta_parsed.fa

To pad exons by 20bp
parse_gencode /path/to/gencode_gtf /path/to/genome_fasta -exon_padding 20

See `func Usage` for more details.
*/
package cmd

import (
	"context"
	"io"
	"os"

	"github.com/grailbio/base/file"
	"github.com/grailbio/base/log"
	"github.com/grailbio/bio/encoding/fasta"
	"github.com/grailbio/bio/fusion/parsegencode"
)

type gencodeFlags struct {
	exonPadding                int
	output                     string // the destination FASTA file. If empty, stdout.
	codingOnly                 bool
	separateJns                bool
	retainedExonBases          int
	wholeGenes                 bool
	collapseTranscripts        bool
	keepMitochondrialGenes     bool
	keepReadthroughTranscripts bool
	keepPARYLocusTranscripts   bool
	keepVersionedGenes         bool
}

// GenerateTranscriptome generates the AF4 transcriptome FASTA file from the
// given inputs. gtfPath is gencode comprehenve antotation file (e.g.,
// gencode.v26.annotation.gtf), and fastaPath is the reference genome (e.g.,
// hg38.fa). gtfPath may be compressed, but fastaPath must be uncompressed.
func GenerateTranscriptome(ctx context.Context, gtfPath, fastaPath string, flags gencodeFlags) {
	if flags.exonPadding < 0 {
		log.Fatal("Pad cannot be negative.")
	}
	var out io.Writer = os.Stdout
	if flags.output != "" {
		var closer func()
		out, closer = createFile(ctx, flags.output)
		defer closer()
	}
	n := 0
	if flags.separateJns {
		n++
		if flags.retainedExonBases < 0 {
			log.Fatal("-retained_exon_bases cannot be negative")
		}
		if flags.exonPadding < 1 {
			log.Fatal("-exon_padding cannot be non-zero if separate junctions are requested")
		}
	}
	if flags.wholeGenes {
		if flags.separateJns || flags.collapseTranscripts {
			log.Fatalf("illegal flag combinations (-whole-genes): %+v", flags)
		}
		n++
	}
	if flags.collapseTranscripts {
		if flags.separateJns || flags.wholeGenes {
			log.Fatalf("illegal flag combinations (-collapse-transcripts): %+v", flags)
		}
		n++
	}
	if n > 1 {
		log.Fatal("-separate_junctions, -whole_genes, and -collapse_transcripts are mutually " +
			"exclusive. Please specify just one")
	}
	records := parsegencode.ReadGTF(ctx, gtfPath,
		// flags.codingOnly is used for exons and transcripts but not genes. When
		// flag.wholeGenes is set, transcript entries need not be inspected in the
		// gtf, only the start and end of gene record will be used; when
		// flags.collapseTranscripts is set, all transcripts that overlap are
		// collapsed into one single final record"
		flags.codingOnly && !flags.wholeGenes && !flags.collapseTranscripts,
		flags.exonPadding,
		flags.separateJns,
		flags.retainedExonBases)

	fastaIn, err := file.Open(ctx, fastaPath)
	if err != nil {
		log.Panic(err)
	}
	defer func() {
		if err := fastaIn.Close(ctx); err != nil {
			log.Panic(err)
		}
	}()
	fasta, err := fasta.New(fastaIn.Reader(ctx))
	if err != nil {
		log.Panic(err)
	}
	switch {
	case flags.wholeGenes:
		parsegencode.PrintWholeGenes(out, fasta, records, flags.codingOnly, flags.exonPadding)
	case flags.collapseTranscripts:
		parsegencode.PrintCollapsedTranscripts(out, fasta, records, flags.codingOnly)
	default:
		parsegencode.PrintParsedGTFRecords(
			out, fasta, records, flags.codingOnly, flags.separateJns,
			false, /* new format */
			flags.keepMitochondrialGenes,
			flags.keepReadthroughTranscripts,
			flags.keepPARYLocusTranscripts,
			flags.keepVersionedGenes)
	}
}
