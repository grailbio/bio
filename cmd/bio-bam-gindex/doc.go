/*Command bio-bam-gindex reads a .bam file and writes a .gbai index
  file.  bio-bam-gindex expects the bam file to arrive on stdin, and
  writes to stdout.  It has a single parameter --shard-size which is
  the approximate shard size in bytes.

  Usage: cat foo.bam | bio-bam-gindex --shard-size=65536 > foo.bam.gbai
*/
package main
