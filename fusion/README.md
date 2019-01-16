## Description

AF4 is an alignment free fusion detector that can work with RNA-seq data, cfRNA data and cfDNA data. It consists of two stages. In stage one, it identifies candidate fusion fragments and the corresponding fusion pairs. In the second stage, a list of filtering criteria are applied to narrow down the fusion candidates. 

Method: AF4 decomposes each fragment into kmers and identifies for each kmer the corresponding genes. Then it infers the max-cover of the fragment that is best explained by a single gene or gene-pair. This effectively removes most of the fragments that are not supporting any fusions while kept real fusion fragments. These fragments are then filtered in the second stage using criteria like the minimum supporting fragments per fusion event, genomic distance of fusion genes etc. 

## Compile/Install 

**to revise** 

## How to run 

**to revise** 

- Denovo mode 
- Target mode 
- First stage alone
- Second stage alone 

## Reference files

- Optional target gene pair files consisting of all COSMIC fusion genes and curated fusion gene pairs. This file can be modified to fit specific needs as long as it is in COSMIC fusion file format. **to provide file used**

- The reference transcriptome file for AF4 required format **to provide file used; cfRNA and cfDNA**

		# This will parse the file to add the junction sequences but will maintain the same fasta header
		
		parse_gencode
			-exon_padding 250 \  # Optional: Pad exons with 250 BP of intron sequence
			-retained_exon_bases 18 \  # Optional: Use 18 exonic bases at junctions (This is the default value)
			-separate_junctions \  # Optional: Add junction info at the end of the fasta record
			-o gencode.v26.250padded_separate_jns.fa \  # Output file name (default is STDOUT)
          /path/to/gencode.v26.annotation.gtf \  # Required positional: Gencode GTF
          /path/to/gencode.v26.transcripts.fa  # Required positional: Gencode Transcript Fasta

		# The following command will process the fasta into the af4-ready format, It will also name the output gencode.v26.250padded_separate_jns_transcripts_parsed_no_mt_no_overlap_no_pary_no_versioned.fa by default
	
		python create_af4ready_fasta.py \
			-F  gencode.v26.250padded_separate_jns.fa \  # Required: Gencode-like fasta. This is th e output from parse_gencode
			-G /path/to/gencode.v26.annotation.gtf  # Required: Gencode GTF
                               
	
	Users can provide their own reference input as long as it is compatible with the AF4-ready format, which is FASTA file with each entry formatted as: 
	
		>`transcript_name`|`gene_name`|`genomic position with index`|`exon lengths`
		Gene/transcript sequence on one line

		E.g.
		>ENST00000420190.6|SAMD11|chr1:923928-944581:49|1021,92,182,51,125,90,17
		CGGAGTCTCCCAAGTCCCCGCCGGGCGGGCGCGCGCCAGTGGACGCGGGTGCACGACTGACGCGGCCCGGGCGGCGGGGCGGGGGC ...

**NOTE**
**`$GRAIL/go/src/grail.com/bio/rna/parse_gencode`**
**`$GRAIL/bio/rna/fusion/analysis/generate_fasta/create_af4ready_fasta.py`**

## Output

The output consists of `all.fa` and `filtered.fa`. `all.fa` is the output after the first stage of the program. The corresponding read pairs can then be pulled from the input and serve as a small input to other fusion calling programs. `filtered.fa` is the final output that went through built-in filters in AF4.

	Each sequence in the output is represented as a FASTA format. 
	The header consists of `|` delimited fields:
	
	- Read Name
	- Fusion Gene Pair: G1/G2
	- Chr:coordinates:gene number on Chr/Chr:coordinates:gene number of Chr
	- number of bases supporting G1/G2
	- range of bases on stitched fragment supporting G1/G2
	
	For example, 
	
	`>E00481:58:H53VWALXX:1:1101:28574:38754:AAATCC+CTATAC 1:N:0:CTGAAGCT+ACGTCCTG|BORCS8/MEF2B|chr19:19176903-19192591:840/chr19:19145568-19170289:839|200/27|1:257/258:284`
		
	`E00481:58:H53VWALXX:1:1101:28574:38754:AAATCC+CTATAC 1:N:0:CTGAAGCT+ACGTCCTG`: read name
	`BORCS8/MEF2B`: fusion gene pairs
	`chr19:19176903-19192591:840`: BORCS8 is on chr19, position 19176903-19192591, and 840th gene from 5' region
	`200/27`: 200 and 27 bases supporting BORCS8 and MEF2B, respectively, on the stitched fragment
	`1:257/258:284` : 1:257bp supports BORCS8 and 258:284bp supports MEF2B.

## Miscellaneous

- PCR duplicates, UMIs

	By default, AF4 deals with 6bp dual UMIs and collapse sequences when their UMIs are within Hamming distance of two. However, when UMIs are absent, PCR duplicates are identified based on sequence similarity, where sequneces are collapsed if they are highly similar. 