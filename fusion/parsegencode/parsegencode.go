// Package parsegencode contains the required methods for parsing a gencode GTF
// annotation and printing out the transcripts to an output file. It has the
// added functionality of allowing the user to pad introns by a specified number
// of base pairs.
package parsegencode

import (
	"bufio"
	"context"
	"io"
	"sort"
	"strings"

	"github.com/grailbio/base/file"
	"github.com/grailbio/base/log"
	"github.com/grailbio/base/tsv"
	"github.com/grailbio/bio/encoding/fasta"
)

type genomicRange struct {
	start, stop int // both ends are closed.
}

type genomicRanges []genomicRange

// Append a genomicRange to g. If last item in g overlaps with r, r is merged.
func (g *genomicRanges) merge(r genomicRange) {
	if overlap, newRange := (*g).overlaps(r); overlap {
		//overlap
		(*g)[len(*g)-1] = newRange
		return
	}
	// No overlap
	*g = append(*g, r)
}

// If the new genomicRange overlaps with the last range of g, return true and update
// the last range to be the union.
func (g genomicRanges) overlaps(gRange2 genomicRange) (bool, genomicRange) {
	if len(g) == 0 {
		return false, genomicRange{}
	}
	gRange1 := g[len(g)-1]
	if gRange1.start <= gRange2.start {
		if gRange2.start <= gRange1.stop+1 {
			// Includes full overlap of gRange2 by gRange1. +1 accounts for touching ranges
			// gRange1:                |----------|
			// gRange2:                           |---...
			// gRange2:                     |---...
			// gRange2:                |---...
			return true, genomicRange{gRange1.start, max(gRange1.stop, gRange2.stop)}
		}
		//else
		// gRange1:                |----------|
		// gRange2:                             |----------|
		return false, genomicRange{}
	} else if gRange1.start-1 <= gRange2.stop {
		// Includes full overlap of gRange1 by gRange2. -1 accounts for touching ranges
		// gRange1:                |----------|
		// gRange2:         ...----|
		// gRange2:              ...---|
		// gRange2:              ...----------|
		return true, genomicRange{gRange2.start, max(gRange1.stop, gRange2.stop)}
	}
	// else
	// gRange1:                |----------|
	// gRange2:  |----------|
	return false, genomicRange{}
}

// Collapse all overlapping ranges in a vector of ranges
func (g *genomicRanges) collapse() {
	newGR := genomicRanges{}
	for _, gr := range *g {
		newGR.merge(gr)
	}
	*g = newGR
}

// Reverse genomicRanges
func (g *genomicRanges) reverse() {
	for i, j := 0, len(*g)-1; i < j; i, j = i+1, j-1 {
		(*g)[i], (*g)[j] = (*g)[j], (*g)[i]
	}
}

// gencodeTranscript will store a gencode transcript record
type gencodeTranscript struct {
	geneID           string
	transcriptID     string
	start            int
	stop             int
	transcriptType   string
	transcriptName   string
	havanaTranscript string
	exons            genomicRanges
	junctions        genomicRanges
}

// GencodeGene will store a gencode gene record
type GencodeGene struct {
	geneID      string // gene name
	transcripts map[string]*gencodeTranscript
	chrom       string
	start       int
	stop        int
	strand      string
	havanaGene  string
	geneName    string
	geneType    string
}

// Return gene.transcripts in ascending order of transcriptIDs.
func (gene GencodeGene) sortedTranscripts() []*gencodeTranscript {
	t := make([]*gencodeTranscript, 0, len(gene.transcripts))
	for k, v := range gene.transcripts {
		if k != v.transcriptID {
			panic(k)
		}
		if v.geneID != gene.geneID {
			panic(gene)
		}
		t = append(t, v)
	}
	sort.Slice(t, func(i, j int) bool {
		return t[i].transcriptID < t[j].transcriptID
	})
	return t
}

// ParseInfoFields parses the "INFO" field of the record to yield a map of key,value pairs.
func parseInfoFields(parsedInfo *map[string]string, info string) {
	for k := range *parsedInfo {
		delete(*parsedInfo, k)
	}
	fields := strings.Split(strings.TrimSpace(info), ";")
	for _, field := range fields {
		field = strings.TrimSpace(field)
		if field == "" {
			continue
		}
		pair := strings.Split(field, " ")
		(*parsedInfo)[pair[0]] = strings.Trim(pair[1], "\"")
	}
}

// ReverseComplement a nucleotide sequence
func reverseComplement(seq string) string {
	var revcomp strings.Builder
	seqSize := len(seq)
	for i := range seq {
		switch x := seq[seqSize-i-1]; x {
		case 'A':
			revcomp.WriteRune('T')
		case 'C':
			revcomp.WriteRune('G')
		case 'T':
			revcomp.WriteRune('A')
		case 'G':
			revcomp.WriteRune('C')
		case 'N':
			revcomp.WriteRune('N')
		default:
			log.Fatal("Unrecognized nucleotide : " + string(x))
		}
	}
	return revcomp.String()
}

func max(x, y int) int {
	if x > y {
		return x
	}
	return y
}

// ReadGTF will read a GTF file line by line into a list of GencodeGenes.  The
// genes are sorted in the increasing order of geneIDs.
func ReadGTF(ctx context.Context,
	pathname string,
	codingOnly bool,
	exonPadding int,
	separateJns bool,
	retainedExonBases int) []*GencodeGene {
	log.Print("Parsing gtf at " + pathname)

	// to retain 5 exonic basepairs, you need to add or subtract 4 from the start and end coordinates
	retainedExonBases--

	var totalGenes, totalTranscripts, retainedTranscripts, totalExons, retainedExons int
	// Define a map to store the data
	records := make(map[string]*GencodeGene)
	in, err := file.Open(ctx, pathname)
	if err != nil {
		log.Fatal(err)
	}
	// create a new scanner and apply the function to the scanner
	scanner := tsv.NewReader(in.Reader(ctx))
	scanner.Comment = '#'
	scanner.LazyQuotes = true
	fields := map[string]string{}
	for {
		// gtfRecord will store data read from one line of the gencode file
		type gtfRecord struct {
			Chrom    string
			Source   string
			Molecule string
			Start    int
			Stop     int
			Score    string // unused floating point value, but may be "."
			Strand   string
			Frame    string
			Fields   string
		}
		var gtfLine gtfRecord
		err := scanner.Read(&gtfLine)
		if err != nil {
			if err != io.EOF {
				log.Panicf("%s: %v", pathname, err)
			}
			break
		}
		parseInfoFields(&fields, gtfLine.Fields)
		if gtfLine.Molecule == "gene" {
			gene := GencodeGene{
				chrom:       gtfLine.Chrom,
				start:       gtfLine.Start,
				stop:        gtfLine.Stop,
				strand:      gtfLine.Strand,
				geneName:    fields["gene_name"],
				geneType:    fields["gene_type"],
				geneID:      fields["gene_id"],
				havanaGene:  fields["havana_gene"],
				transcripts: make(map[string]*gencodeTranscript)}
			records[gene.geneID] = &gene
			totalGenes++
			continue
		}
		if gtfLine.Molecule == "transcript" {
			geneID := fields["gene_id"]
			gene, ok := records[geneID]
			if !ok {
				log.Fatal("GTF not sorted properly.")
			}
			totalTranscripts++
			if codingOnly && !isCodingBiotype(fields["transcript_type"]) {
				continue
			}
			transcript := gencodeTranscript{
				start:            gtfLine.Start,
				stop:             gtfLine.Stop,
				geneID:           fields["gene_id"],
				transcriptName:   fields["transcript_name"],
				transcriptType:   fields["transcript_type"],
				transcriptID:     fields["transcript_id"],
				havanaTranscript: fields["havana_transcript"],
				exons:            make([]genomicRange, 0),
				junctions:        make([]genomicRange, 0)}
			if transcript.havanaTranscript == "" {
				transcript.havanaTranscript = "-"
			}
			gene.transcripts[transcript.transcriptID] = &transcript
			retainedTranscripts++
			continue
		}
		if gtfLine.Molecule == "exon" {
			totalExons++
			if codingOnly && !isCodingBiotype(fields["transcript_type"]) {
				continue
			}
			geneID := fields["gene_id"]
			transcriptID := fields["transcript_id"]
			transcript, ok := records[geneID].transcripts[transcriptID]
			if !ok {
				log.Fatal("GTF not sorted properly.")
			}

			var gRange genomicRange
			gRange.start = gtfLine.Start
			gRange.stop = gtfLine.Stop

			if !separateJns {
				if exonPadding > 0 {
					gRange.start -= exonPadding
					if gRange.start < 0 {
						gRange.start = 0
					}
					gRange.stop += exonPadding
				}
				transcript.exons.merge(gRange)
			} else {
				// Add the exon to the list of exons
				transcript.exons = append(transcript.exons, gRange)
				// We want separate junctions. Store the junction ranges in gencodeTranscript.junctions
				// as genomicRanges
				exonStartJn := genomicRange{gRange.start - exonPadding, gRange.start + retainedExonBases}
				exonStopJn := genomicRange{gRange.stop - retainedExonBases, gRange.stop + exonPadding}
				if records[geneID].strand == "-" {
					// do this is reverse order so the junctions are printed out in the correct order
					tmpJn := exonStartJn
					exonStartJn = exonStopJn
					exonStopJn = tmpJn
				}

				// If there is a previous juntion, see if it merges with the start of this one.
				// These merges are different from the non-separate-junctions ones . E.g.
				//                     ...-----=|------||
				//                          ||------|=-------
				//              =>             =|---|=             (padding + intron overlap)
				//             !=>          ||-=-----=-||          (largest span)
				if overlap, _ := transcript.junctions.overlaps(
					exonStartJn); overlap {
					// overlaps
					idx := len(transcript.junctions) - 1
					if records[geneID].strand == "+" {
						// Need to replace only the stop
						transcript.junctions[idx].stop = exonStartJn.stop
					} else {
						// Need to replace only the start
						transcript.junctions[idx].start = exonStartJn.start
					}
				} else {
					// no overlap. append start jn
					transcript.junctions = append(transcript.junctions, exonStartJn)
				}
				// Append the stop jn
				transcript.junctions = append(transcript.junctions, exonStopJn)
			}
			retainedExons++
		}
	}
	log.Printf("Successfully parsed %d genes, %d transcripts and %d exons\n", totalGenes,
		totalTranscripts, totalExons)
	log.Printf("Retained %d transcripts and %d exons\n", retainedTranscripts, retainedExons)
	if err := in.Close(ctx); err != nil {
		log.Panic(err)
	}

	genes := make([]*GencodeGene, 0, len(records))
	for _, gene := range records {
		genes = append(genes, gene)
	}
	sort.Slice(genes, func(i, j int) bool {
		return genes[i].geneID < genes[j].geneID
	})
	return genes
}

// Obtained from https://www.gencodegenes.org/gencode_biotypes.html
var codingBiotypes = []string{
	"protein_coding",
	"nonsense_mediated_decay",
	"non_stop_decay",
	"IG_C_gene",
	"IG_D_gene",
	"IG_J_gene",
	"IG_LV_gene",
	"IG_V_gene",
	"TR_C_gene",
	"TR_J_gene",
	"TR_V_gene",
	"TR_D_gene",
	"polymorphic_pseudogene"}

func isCodingBiotype(query string) bool {
	for _, q := range codingBiotypes {
		if query == q {
			return true
		}
	}
	return false
}

func write(out *bufio.Writer, s string) {
	if _, err := out.WriteString(s); err != nil {
		log.Panic(err)
	}
}

func flush(out *bufio.Writer) {
	if err := out.Flush(); err != nil {
		log.Panic(err)
	}
}

// PrintParsedGTFRecords will print transcript fasta records to an output file given a map of fasta
// records and a map of gencode genes.
func PrintParsedGTFRecords(
	out io.Writer,
	fasta fasta.Fasta,
	newGTFRecords []*GencodeGene,
	codingOnly bool,
	separateJns bool) {
	w := bufio.NewWriter(out)
	defer flush(w)
	for geneID := range newGTFRecords {
		gene := *newGTFRecords[geneID]
		for _, transcript := range gene.sortedTranscripts() {
			write(w, ">"+
				strings.Join([]string{transcript.transcriptID,
					gene.geneID,
					gene.havanaGene,
					transcript.havanaTranscript,
					transcript.transcriptName,
					gene.geneName,
					"NA",
					"NA",
					"\n"}, "|"))
			for _, exon := range transcript.exons {
				seq, err := fasta.Get(gene.chrom, uint64(exon.start-1), uint64(exon.stop))
				if err != nil {
					log.Panic(err)
				}
				if gene.strand == "+" {
					write(w, seq)
				} else {
					write(w, reverseComplement(seq))
				}
			}
			for _, junction := range transcript.junctions {
				write(w, "|")
				seq, err := fasta.Get(gene.chrom, uint64(junction.start-1), uint64(junction.stop))
				if err != nil {
					log.Panic(err)
				}
				if gene.strand == "+" {
					write(w, seq)
				} else {
					write(w, reverseComplement(seq))
				}
			}
			write(w, "\n")
		}
	}
}

// PrintWholeGenes will print whole gene records to an output file given a map of fasta
// records and a map of gencode genes.
func PrintWholeGenes(
	out io.Writer,
	fasta fasta.Fasta,
	newGTFRecords []*GencodeGene,
	codingOnly bool,
	genePadding int) {
	w := bufio.NewWriter(out)
	defer flush(w)
	for geneID := range newGTFRecords {
		gene := *newGTFRecords[geneID]
		if codingOnly && !isCodingBiotype(gene.geneType) {
			continue
		}
		write(w, ">"+
			strings.Join([]string{gene.geneID,
				gene.havanaGene,
				gene.geneName,
				gene.geneType,
				"\n"}, "|"))
		seq, err := fasta.Get(gene.chrom, uint64(gene.start-genePadding-1), uint64(gene.stop+genePadding))
		if err != nil {
			log.Panic(err)
		}
		if gene.strand == "+" {
			write(w, seq)
		} else {
			write(w, reverseComplement(seq))
		}
		write(w, "\n")
	}
}

// PrintCollapsedTranscripts will print all exons (with or without padding) for gene records to an
// output file given a map of fasta records and a map of gencode genes. Overlapping exons will be
// collapsed into single entries.
func PrintCollapsedTranscripts(
	out io.Writer,
	fasta fasta.Fasta,
	newGTFRecords []*GencodeGene,
	codingOnly bool) {
	w := bufio.NewWriter(out)
	defer flush(w)
	for geneID := range newGTFRecords {
		gene := *newGTFRecords[geneID]
		if codingOnly && !isCodingBiotype(gene.geneType) {
			continue
		}
		write(w, ">"+
			strings.Join([]string{gene.geneID,
				gene.havanaGene,
				gene.geneName,
				gene.geneType,
				"\n"}, "|"))
		var geneExons genomicRanges
		for _, transcript := range gene.sortedTranscripts() {
			geneExons = append(geneExons, transcript.exons...)
		}
		sort.SliceStable(geneExons, func(i, j int) bool {
			return (geneExons[i].start < geneExons[j].start) ||
				(geneExons[i].start == geneExons[j].start && geneExons[i].stop < geneExons[j].stop)
		})
		geneExons.collapse()
		if gene.strand == "-" {
			geneExons.reverse()
		}
		for i, exon := range geneExons {
			if i > 0 {
				write(w, "|")
			}
			seq, err := fasta.Get(gene.chrom, uint64(exon.start-1), uint64(exon.stop))
			if err != nil {
				log.Panic(err)
			}
			if gene.strand == "+" {
				write(w, seq)
			} else {
				write(w, reverseComplement(seq))
			}
		}
		write(w, "\n")
	}
}
