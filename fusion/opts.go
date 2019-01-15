package fusion

type Opts struct {
	UMIInRead bool
	UMIInName bool
	// LowComplexityFraction determines whether a fragment (or read) should be
	// dropped because it contains too many repetition of the same base types.  If
	// LowComplexityFraction of bases in a fragment are such repetitions, it is
	// dropped without further analyses.
	LowComplexityFraction float64
	// KmerLength is the length of kmer used to match DNA sequences.
	KmerLength int
	Denovo     bool
	// Is this an unstranded RNA-seq library? Default is stranded.
	UnstrandedPrep bool

	// maxGap specifies the max gap allowed between
	// two consecutive ranges (this is to tolerate sequence errors).
	MaxGap int

	// MaxHomology is the max overlap allowed b/w genes in a fusion
	MaxHomology int
	// Min base evidence for a gene in the fusion.
	MinSpan int

	// MaxKmerFrequency the number of genes a kmer belongs to, default 5. Used to
	// be --cap flag in the C++ code.
	MaxGenesPerKmer int

	// MaxGeneCandidatesPerFragment caps the number of genes considered for each
	// fragment. Doing so will reduce the number of pairs that will be processed
	// for every read.
	//
	// TODO(saito,xyang) report read through events but differentiating this from
	// the overlapping genes
	MaxGeneCandidatesPerFragment int

	// MaxProximityDistance is the distance cutoff below which a candidate will be
	// rejected as a readthrough event
	MaxProximityDistance int

	// MaxProximityGenes is number of genes separating a gene pair (If on the
	// same chromsosome) below which they will be flagged as read-through events.
	MaxProximityGenes int

	// MaxGenePartners x capping number of partners a gene can have
	// this is used in the filtering stage.
	MaxGenePartners int

	// Minimum number of supporting reads required
	// to consider a fusion
	MinReadSupport int
}

// DefaultOpts sets the default values to Opts.
var DefaultOpts = Opts{
	UMIInRead:                    false,
	UMIInName:                    false,
	KmerLength:                   19,
	UnstrandedPrep:               true,
	MaxGap:                       9, //kmerLength/2
	MaxHomology:                  15,
	MinSpan:                      25,
	LowComplexityFraction:        0.9,
	MaxGenesPerKmer:              2,
	MaxGeneCandidatesPerFragment: 5,

	MaxProximityDistance: 1000,
	MaxProximityGenes:    0,
	MaxGenePartners:      5,
	MinReadSupport:       2,
}
