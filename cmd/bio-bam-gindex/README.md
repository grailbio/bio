bio-bam-gindex
==============

## Background

A [gindex](https://github.com/grailbio/bio/tree/master/encoding/bam/gindex.go) (.gbai)
file is an alternate index format for BAM files.  The gindex format
allows readers to seek into a BAM file more efficiently because the
gindex can contain pointers into the BAM file at a finer granularity.

The [encoding/bamprovider](https://github.com/grailbio/bio/tree/master/encoding/bamprovider)
interface can use a .gbai file as a drop-in replacement for a .bai file.

## Usage

Command bio-bam-gindex reads a .bam file and writes a .gbai index file.
bio-bam-gindex expects the bam file to arrive on stdin, and
writes to stdout.  It has a single parameter --shard-size which is
the approximate shard size in bytes.

Example usage:

    cat foo.bam | bio-bam-gindex --shard-size=65536 > foo.bam.gbai

