# bio-bam-sort

bio-bam-sort is a tool for producing BAM or PAM files from aligner outputs. It
can be used as a backend for [bwa](http://bio-bwa.sourceforge.net/).

Example usage:

    bwa ..... | bio-bam-sort -sam out1.shard
    bwa ..... | bio-bam-sort -sam out2.shard
    bio-bam-sort -pam foo.pam out1.shard out2.shard
