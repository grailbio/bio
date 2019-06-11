/*Package interval implements interval-union operations in a manner optimized
  for sets of genomic coordinates represented by BED files.
  (Note the 'union'.  Overlapping intervals are merged, not tracked
  separately; it is currently necessary to use another package when that is not
  the desired behavior.)
  It assumes every position fits in a PosType, which is currently defined as
  int32 since that's what BAM files are limited to.
*/
package interval
