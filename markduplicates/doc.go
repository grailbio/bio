/*Command bio-mark-duplicates marks or removes duplicates from .bam
  files.

  This tool is meant to replicate the behavior of picard MarkDuplicates
  in a more scalable, and efficient way.  The intention is to make
  bio-mark-duplicates scale well with machines exceeding 16 cores.

  Duplicate Marking Concepts:

  At the conceptual level, this tool considers two reads A and B as
  duplicates (isDuplicate(A, B)) if their:
    1) reference
    2) unclipped 5' position
    3) read direction (orientation)
  are ALL identical.

  Two pairs P1 and P2 are considered duplicates of each other, if
  isDuplicate(P1.leftRead, P2.leftRead) and isDuplicate(P1.rightRead,
  P2.rightRead).  Left vs right is determined by the unclipped 5'
  position of each read in the pair.

  Mapped pairs vs. Mapped-Unmapped pairs: For some read pairs, both
  reads will be mapped (mapped pairs).  For other read pairs, only one
  of the reads will be mapped (mapped-unmapped pairs).  A mapped pair
  can be a duplicate of another mapped pair, but a mapped pair P1 may
  NOT be a duplicate of a mapped-unmapped pair P2 because one read of
  P2 will have no alignment position, and thus cannot be equal to one
  of the mapped reads of P1.

  However, the mapped read of a mapped-unmapped pair can be considered
  a duplicate of one read on a mapped pair.  So in this example, P2.left
  could be a duplicate of P1.left.  We call P2.left a "mate-unmapped read".

    P1: left(chr1, 1020, F) right(chr1, 1040, R)
    P2: left(chr1, 1020, F) right(chr1, 0, ?)

    P1 is not a duplicate of P2, but P2.left is a duplicate of P1.left.

  After identifying the duplicates, this tool will select a primary
  pair or read for each set of duplicates.  The primary will be the
  duplicate with the highest score based on the sum of its base
  qualities.  To break ties, a higher priority is given to reads that
  appear earlier in the bam input.

  In choosing a primary, pairs are given priority over mate-unmapped
  reads.  So if a mate-unmapped read is found to be a duplicate of one
  read in a mapped pair, the pair would always be the primary, and the
  mate-unmapped read would always be the duplicate.  If two
  mate-unmapped reads are duplicates of each other, but not duplicates
  of a mapped pair, then the primary is chosen from the two
  mate-unmapped reads.

  After identifying the primary and the duplicates, this tool can be
  configured to mark each read with the duplicate flag 1024, or to
  remove each of the duplicate reads.

  Tagging:

  If the caller specifies the "tag-duplicates" parameter, the tool
  will attach auxiliary tags DI, DL, DS, and DT to the output.

  DI is the duplicate index of a duplicate set.  All pairs in a
  duplicate set, including the primary, will have the same value for
  DI; the value for DI is the file index of the left-most read of the
  primary duplicate pair.  DI is not set for mate-unmapped reads.

  DL is the number of library (LB aka PCR) duplicate pairs in the
  duplicate set. This is the DS value minus the number of "SQ"
  duplicates pairs in the duplicate set.

  DS is the number of pairs in the duplicate set.  DS is not set on
  mate-unmapped reads, and it also does not count mate-unmapped
  duplicates.

  DT is set on duplicate pairs (not the primary) and mate-unmapped
  reads.  It is set to "SQ" for optical duplicates, and "LB" for all
  other duplicates.

  Implementation:

  The implementation splits the input bam file into non-overlapping
  shards, and processes each of those shards in separate workers, and
  possibly in parallel with the other shards.  Each worker is
  responsible for taking the reads in one shard, and outputting the
  same reads in the same order (except for marking or removing
  duplicates).

  Matching up pairs:

  To determine which pairs are duplicates, a worker must read both
  reads in a pair to determine each read's 5' position and
  orientation.  Each read record contains a pointer to its mate, but
  it does not contain the orientation or the unclipped position which
  is in the mate's flags and cigar respectively.

  Matching up both reads in a pair is somewhat tricky.  For most of
  the reads that live in a particular shard, the read's mate will
  appear within that shard, and a worker can easily match up those
  pairs.  However, some pairs will span from inside the shard to
  outside the shard.  Most of these read pairs will have their reads
  relatively close together, so the a worker will read an additional
  pair-padding distance on either end of the shard for matching up
  pairs, but even that is not enough.  Sometimes, a pair spans from
  inside the shard to beyond the pair-padding; we'll call these
  "distant-pairs".  To deal with distant-pairs, this tool pre-scans
  the entire bam input, to identify and save the distant-pairs to
  memory so that the workers can lookup the distant pair reads without
  needing to seek, decompress, and read them from the bam file.  Note
  that the distant-pairs allow all pairs to be resolved, so the
  pair-padding is a memory optimization to avoid storing mates that
  live in the pair-padding in memory for the duration of the entire
  run.  Pair-padding should be chosen so that most mates will fall
  into the pair-padding.

  Matching up duplicates:

  Duplicates are matched up using their 5' positions, and since a 5'
  position can differ from the alignment's start position, each shard
  needs some fuzziness at its boundaries to ensure all potential
  duplicates will be compared against each other.  For example:

           shard1                  shard2
   |---------------------|-------------------------|
                      |cccc|--------------|         read1 (with clipping)
                      5    S              E         5', Begin, End

                      |c|-----------------|         read2 (with clipping)
                      5 S                 E         5', Begin, End

                clip-pad           shard2            clip-pad
             |-----------|-------------------------|-----------|

  In this example, read1 and read2 are duplicates according to their 5'
  position, but the two reads will reside in different shards based on
  their alignment Begin positions.

  To handle these overlap cases, each shard has "clip-padding" on each
  side of the shard.  When searching for duplicates in each shard, the
  worker considers each read that is either in the clip-padding or in
  the shard.  So in the example, read1 is in shard2, and read2 is in
  shard2's clip-padding, but the shard2 worker still compares the two
  reads with each other, decides they are duplicates, and marks one as
  a duplicate based on their scores.  Note that the worker for shard1
  also compares read1 and read2 and also decides they are duplicates,
  and marks one based on their scores.  Since the scoring is
  deterministic, shard1 and shard2 agree on which read to mark as
  duplicate.

  Clip-padding and pair-padding serve different purposes.
  Clip-padding is for correctness and must exceed the largest clip
  distance in the input file.  Pair-padding is a memory optization.
  The complete shard diagram looks like this:

   shard-pad  clip-pad            shard1            clip-pad   shard-pad
  |---------|-----------|-------------------------|-----------|---------|


  Output ordering:

  As the workers complete each shard, they output the shard's marked
  reads to an output queue that preserves the original order of the
  shards.  A writer reads shards from the queue in order, and writes the
  records to the bam file output.  The output queue has a maximum size,
  and when full only allows a worker to insert the shard that the writer
  currently needs.  This ensures that the queue does not grow too long
  when a worker takes a long time to process a particular shard.
*/
package markduplicates
