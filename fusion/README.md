# AF4

AF4 is an alignment free fusion detector that can work with RNA-seq data, cfRNA
data and cfDNA data. It consists of two stages. In stage one, it identifies
candidate fusion fragments and the corresponding fusion pairs. In the second
stage, a list of filtering criteria are applied to narrow down the fusion
candidates.

Method: AF4 decomposes each fragment into kmers and identifies for each kmer the
corresponding genes. Then it infers the max-cover of the fragment that is best
explained by a single gene or gene-pair. This effectively removes most of the
fragments that are not supporting any fusions while kept real fusion
fragments. These fragments are then filtered in the second stage using criteria
like the minimum supporting fragments per fusion event, genomic distance of
fusion genes etc.

Contact: saito@grail.com, xyang@grail.com

## Installation

Currently, AF4 only supports UNIX platforms with ivy bridge CPU or newer.  Our reference machines are Ubuntu16.04 with 256GiB memory and 56 to 64 CPUs.

- If you have go environment set up, please use the following command directly

    	go install github.com/grailbio/bio/cmd/bio-fusion

	The binary `bio-fusion` should be ready to use, for usage, type

		bio-fusion -h

- Setup go then install

	- Assuming your current working folder is `$HOME/go_install`
	- Obtain pre-compiled go package (this works for linux only, please check sha256 to make sure the `.tar.gz` file is properly downloaded)

			wget https://dl.google.com/go/go1.11.5.linux-amd64.tar.gz
			tar -C $HOME/go_install -xvzf go1.11.5.linux-amd64.tar.gz

		The go binary path should be at: `$HOME/go_install/go/bin/go`

	- Set up go env

			export GOROOT=${HOME}/go_install/go
			export PATH="${HOME}/go_install/go/bin:${PATH}"

	- Set up go work space & install fusion package

			# create workspace assuming you want to use $HOME/workdir folder
			mkdir -p $HOME/workdir
			cd $HOME/workdir
			go mod init playground/
			go install github.com/grailbio/bio/cmd/bio-fusion

	- The binary `bio-fusion` should be ready to use, for usage, type

			bio-fusion -h



## Running

    bio-fusion -r1=a0.fastq.gz,b0.fastq.gz,... -r2=a1.fastq.gz,b1.fastq.gz,... -transcript=transcipt.fa [-cosmic-fusion=fusion.txt] [-fasta-output=all.fa] [-filtered-output=filtered.fa]

- Flags `r1=...` and `r2=...` specify comma-separated lists of FASTQ files. In
  the above example, files `a0.fastq.gz` and `a1.fastq.gz` should contain paired
  reads from R1 and R2, respectively.  The two lists must contain the same
  number of path names, and file-pair in the lists must contain exactly the same
  number of reads.  If a path ends with `.gz` or `.bz2`, they will be
  decompressed by gzip and bz2, respectively.

- Flag `-transcript` specifies the transcriptome. The next section describes the
  format of this file in more detail.

- Flag `-cosmic-fusion` is used only in target mode. It provides curated fusion
  gene pairs. This file should be in COSMIC TSV format.  The first line is a
  header, which is ignored by `bio-fusion`.  The rest of file should contain at
  least one TSV column, in form "gene1/gene2", per line. The rest of the columns
  are ignored by AF4. For example:

```
Genes	Samples	Mutations	Papers
AAGAB/FOSL2	10	1	1
ABTB1/SDC2	10	1	1
ACAD10/SELP	10	1	1
ACBD6/RRP15	2	2	1
```


- Flag `-fasta-output` specifies the path name of the 1st-stage output. The
  default value is `all.fa`. The output file format is described in a later
  section.

- Flag `-filtered-output` specifies the path name of the 2nd-stage output. It is
  a subset of the `fasta-output`. The default value is `filtered.fa`.

- Passing `-h` will show more minor flags supported by `bio-fusion`.

### PCR duplicates, UMIs

By default, AF4 deals with 6bp dual UMIs and collapse sequences when their UMIs
are within Hamming distance of two. However, when UMIs are absent, PCR
duplicates are identified based on sequence similarity, where sequneces are
collapsed if they are highly similar.

## Reference transcriptome

The transcriptome is a FASTA file. Each key should be of form
    >transcriptname|genename|chr:start-end:index|len0,len1,...,lenN

The sequence should list all the exons without any separator. E.g.:

    >ENST00000420190.6|SAMD11|chr1:923928-944581:49|1021,92,182,51,125,90,17
    CGGAGTCTCCCAAGTCCCCGCCGGGCGGGCGCGCGCCAGTGGACGCGGGTGCACGACTGACGCGGCCCGGGCGGCGGGGCGGGGGC ...

`start`, `end` are 1-based positions in the chromosome, both closed. `index` is
the rank of the transcript's `start` position within the chromosome. That is,
the transcript with the smallest start position for a particular chromosome has
index of 0, the 2nd smallest start position will have index of 1, so on.
`len0`, `len1`, ... are length of each exon, but the exon lengths are currently
unused by AF4.

## Generating reference transcriptome from Gencode annotations

The AF4 transcriptome file can be generated from [Gencode annotation GTF
files](https://www.gencodegenes.org/human/release_26.html) and a reference human
genome, usually hg38.  The reference files can be found in the following
locations:

```
ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_26/gencode.v26.annotation.gtf.gz
http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz # you must uncompress the file before feeding it to bio-fusion
```

To generate
`gencode.v26.250padded_separate_jns_transcripts_parsed_no_mt_no_overlap_no_pary_no_versioned.fa`,
used for cfRNA datasets, run the following:

```
bio-fusion -generate-transcriptome -exon-padding 250 -retained-exon-bases 18  -separate-junctions -output gencode.v26.250padded_separate_jns_transcripts_parsed_no_mt_no_overlap_no_pary_no_versioned.fa /path/to/gencode.v26.annotation.gtf /path/to/hg38.fa
```

For non-padded transcrptome used for other datasets, run the folowing:

```
bio-fusion -generate-transcriptome -output gencode.v26.whole_genes.fa -keep-mitochondrial-genes -keep-readthrough-transcripts -keep-pary-locus-transcripts -keep-versioned-genes -/path/to/gencode.v26.annotation.gtf /path/to/hg38.fa
```

Pre-generated transcriptome files used in our benchmark are found in
[s3://grail-publications/2019-ISMB/references](https://grail-publications.s3-us-west-2.amazonaws.com/2019-ISMB/list.html).

## Output format

AF4 produces two outputs, `all.fa` and `filtered.fa`. File `all.fa` is the output
after the first stage of the program.  `filtered.fa` is the final output that
went through the built-in filters in AF4.

Each sequence in the output is represented as a FASTA format.
The header consists of `|`-delimited fields:

- Read name: the R1 name is reported here.

- Comma-separated fusion gene pairs. Each genepair is in form `Ga0/Ga1`.  Most often,
  this column contains just one gene pair, `Ga0/Ga1`, which means that the left
  half of the fragment matches `Ga0`, the right half of the fragment matches
  `Ga1`.  Only when one fragment matches multiple genepairs equially well, more
  than one pairs are shown in this column.

- Comma-separated genepair positions. Each position-pair is in form
  `chr0:start0-end0:index0/chr1:start1-end1:index1`. The positions is copied
  from the transcriptome.  See the previous section for the meaning of
  `start[01]`, `end[01]`, and `index[01]`. The number of positionpairs in this
  column always matches the number of genepairs in the previous column.

- Comma-separated numbers of bases supporting G1/G2. The number of entries in
  this column always matches the number of genepairs.

- Range of bases on stitched fragment supporting G1/G2

The sequence is made up of the stitched readpair. When the readpair cannot be
stitched (e.g., because they don't overlap), then the line contains the R1
sequence, followed by `|`, followed dy the reverse complement of the R2
sequence.

For example,

    >E00481:58:H53VWALXX:1:1101:28574:38754:AAATCC+CTATAC 1:N:0:CTGAAGCT+ACGTCCTG|BORCS8/MEF2B|chr19:19176903-19192591:840/chr19:19145568-19170289:839|200/27|1:257/258:284

- `E00481:58:H53VWALXX:1:1101:28574:38754:AAATCC+CTATAC 1:N:0:CTGAAGCT+ACGTCCTG`: read name
- `BORCS8/MEF2B`: fusion gene pairs
- `chr19:19176903-19192591:840`: BORCS8 is on chr19, position 19176903-19192591, and 840th gene from 5' region
- `200/27`: 200 and 27 bases supporting BORCS8 and MEF2B, respectively, on the stitched fragment
- `1:257/258:284` : 1:257bp supports BORCS8 and 258:284bp supports MEF2B.

## Running the benchmarks

References and FASTQ files used in the ISMB paper are in
[s3://grail-publications/2019-ISMB](https://grail-publications.s3-us-west-2.amazonaws.com/2019-ISMB/list.html).
Directory [benchmark] also contains helper scripts to run the benchmarks.

### Simulated dataset

[Simulated fusion FASTQ files](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4797269/)
are in [s3://grail-publications/2019-ISMB/simulated_benchmark](https://grail-publications.s3-us-west-2.amazonaws.com/2019-ISMB/list.html). First download the [reference files](https://grail-publications.s3-us-west-2.amazonaws.com/2019-ISMB/list.html) and the FASTQ files in a local directory, then run:

```
bio-fusion -r1=path_to/r1.fastq.gz -r2=path_to/r2.fastq.gz -transcript=path_to/gencode.v26.250padded_separate_jns_transcripts_parsed_no_mt_no_overlap_no_pary_no_versioned.fa -umi-in-name -max-genes-per-kmer=2 -max-proximity-distance=1000 -max-proximity-genes=0 [-cosmic_fusion=path_to/all_pair_art_lod_gpair_merged.txt]
```

The results will be created in `./all.fa` and `./filtered.fa`.

### CfRNA dataset

CfRNA files are in [s3://grail-publications/2019-ISMB/rna_benchmark](https://grail-publications.s3-us-west-2.amazonaws.com/2019-ISMB/list.html). First download the [reference files](https://grail-publications.s3-us-west-2.amazonaws.com/2019-ISMB/list.html) and the FASTQ files in a local directory, then run:

```
bio-fusion -r1=path_to/r1.fastq.gz,... -r2=path_to/r2.fastq.gz,... -transcript=path_to/gencode.v26.250padded_separate_jns_transcripts_parsed_no_mt_no_overlap_no_pary_no_versioned.fa -umi-in-name -max-genes-per-kmer=2 -max-proximity-distance=1000 -max-proximity-genes=0 [-cosmic_fusion=path_to/all_pair_art_lod_gpair_merged.txt]
```

The results will be created in `./all.fa` and `./filtered.fa`.

### CfDNA dataset

CfDNA files are in [s3://grail-publications/2019-ISMB/titration_benchmark](https://grail-publications.s3-us-west-2.amazonaws.com/2019-ISMB/list.html). First download the [reference files](https://grail-publications.s3-us-west-2.amazonaws.com/2019-ISMB/list.html) and the FASTQ files in a local directory, then run:

```
bio-fusion -r1=path_to/r1.fastq.gz,... -r2=path_to/r2.fastq.gz,... -transcript=path_to/gencode.v26.whole_genes.fa -umi-in-name -max-genes-per-kmer=2 -max-proximity-distance=1000 -max-proximity-genes=0 [-cosmic_fusion=path_to/all_pair_art_lod_gpair_merged.txt]
```

The results will be created in `./all.fa` and `./filtered.fa`.
