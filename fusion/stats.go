package fusion

// Stats represents high-level statistics during the run of the stage1 of AF4.
type Stats struct {
	// LowComplexityReads2 is the # of readpairs where both reads are
	// found to have low complexity.
	LowComplexityReads2 int
	// LowComplexityReads2 is the # of readpairs where one of the reads are found
	// to have low complexity.
	LowComplexityReads1 int
	// LowComplexityStitched is the # of readpairs that were successfully
	// stitched, but then found to have low complexity.
	LowComplexityReadsStitched int
	// Stitched is the # of reads successuflly stitched
	Stitched int
	// RawGenes is the total genes found during kmer lookup.
	RawGenes int
	// Genes is the total genes found during kmer lookup, after
	// Opts.MaxGeneCandidatesPerFragment cutoff.
	Genes int
	// Fragments counts the total number of fragments processed.
	Fragments int
	// FragmentsWithMatchingGenes[k] (0<=k<4) counts the total # of fragments that
	// are found to have k genes matching one of its kmers. The last element in
	// this array counts all the fragmenst with >=4 matching genes.
	FragmentsWithMatchingGenes [5]int
	// Ranges is the total # of ranges covered by any gene.
	RawRanges int
	// Ranges is the total # of ranges covered by any gene, after
	// Opts.MaxGeneCandidatesPerFragment cutoff.
	Ranges int
}

// Merge adds the field values of the two Stats objects and creates new Stats.
func (s Stats) Merge(o Stats) Stats {
	s.LowComplexityReads2 += o.LowComplexityReads2
	s.LowComplexityReads1 += o.LowComplexityReads1
	s.LowComplexityReadsStitched += o.LowComplexityReadsStitched
	s.Stitched += o.Stitched
	s.Genes += o.Genes
	s.RawGenes += o.RawGenes
	for i, n := range o.FragmentsWithMatchingGenes {
		s.FragmentsWithMatchingGenes[i] += n
	}
	s.Fragments += o.Fragments
	s.Ranges += o.Ranges
	s.RawRanges += o.RawRanges
	return s
}
