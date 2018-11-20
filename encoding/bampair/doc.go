/*Package bampair provides a way to get the mate of each read when
  reading a BAM/PAM file that is sorted by position.  Package bampair
  assumes that the user is reading the BAM/PAM file in a sharded way,
  and that while reading a shard, it is possible to store all the
  shard's records in memory while processing the shard.  Package
  bampair makes it possible for a user to process each record and its
  mate in file order, making each shard's processing deterministic.

  To use bampair, the user first calls GetDistantMates() with a
  bamprovider and a list of shards, which returns a DistantMateTable.
  The DistantMateTable contains the mate for each read who's mate is
  *not* in the same shard.  For example, if R1 and R2 are mates, and
  R1 is in shard2 and R2 is in shard4, then the DistantMateTable will
  contain both R1 and R2.  On the other hand, if R1 and R2 are both in
  shard3, then the DistantMateTable will contain neither R1 nor R2.

  After calling GetDistantMates(), the user can then open each shard
  and find the mate for reach record in the shard using the following
  procedure: For a record who's mate is in the same shard, the user
  must store the record in memory and continue reading the shard until
  the user encounters the mate.  For a record who's mate is not in the
  same shard, the user can call DistantMateTable.GetMate() right away
  to retrieve the mate.  For an usage example, see
  ExampleResolvePairs() in distant_mates_test.go.

  Some applications may need to add padding to beginning and end of
  each shard.  In this case, if R1 and R2 are in the same padded
  shard, then neither will be in distant mates.  If R1 is in the
  padded shard, and R2 is not, then R2 will be in distant mates.
*/
package bampair
