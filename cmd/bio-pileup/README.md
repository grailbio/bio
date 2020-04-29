# bio-pileup

Given a BAM or PAM, and a BED file describing genomic positions of interest,
bio-pileup reports the number of reads supporting each allele at each position.
This command is similar to "bcftools mpileup".

There are options for "collapsing" the two ends of a read-pair together, or
entire duplicate-sets ("bags") identified by doppelmark.

Only SNPs are currently reported, but indel support is very likely to be added
in the future.

Sample usage:
bio-pileup \
    --bed my-regions.bed \
    --out output-prefix \
    my.bam \
    ref.fa

Run "bio-pileup --help" for more details.
