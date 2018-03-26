# PAM file format

PAM replaces BAM. PAM can be translated from BAM and vice versa without
information loss. Compared to BAM, PAM is smaller and faster. Major differences
are the following.

- It uses larger compression blocks, by default 8MB. BAM uses 16KB blocks.

- It stores a field in a separate file. This layout has three benefits.

  - It results in better compression, since values in a field are more similar to each other.
  - An application does not need these fields during reads, it can skip the blocks altogether.
  - An application that wants to modify some, but not all the fields can do that
    efficiently.  markduplicates is one such application.

- PAM supports file sharding. It allows creating a separate file for a subrange
  of alignment positions, [<start-refid,start-alignpos>,
  <limit-refid,limit-alignpos>).  Row sharding improves parallelism during
  writes, although it doesn't help during reads (it doesn't hurt either).


# File-system layout

A single BAM file is translated into a set of files in a directory. These Below
is an example:

        foo/0:0,46:1653469.index
        foo/0:0,46:1653469.aux
        foo/0:0,46:1653469.coord
        foo/0:0,46:1653469.cigar
        foo/0:0,46:1653469.flags
        foo/0:0,46:1653469.mapq
        foo/0:0,46:1653469.matepos
        foo/0:0,46:1653469.materefid
        foo/0:0,46:1653469.name
        foo/0:0,46:1653469.qual
        foo/0:0,46:1653469.seq
        foo/0:0,46:1653469.templen
        foo/46:1653469,-:-.index
        foo/46:1653469,-:-.aux
        foo/46:1653469,-:-.coord
        foo/46:1653469,-:-.cigar
        foo/46:1653469,-:-.flags
        foo/46:1653469,-:-.mapq
        foo/46:1653469,-:-.matepos
        foo/46:1653469,-:-.materefid
        foo/46:1653469,-:-.name
        foo/46:1653469,-:-.qual
        foo/46:1653469,-:-.seq
        foo/46:1653469,-:-.templen

`foo/` is the prefix. The `0:0,46:1653469` part defines the coordinate range,
(reference sequence IDs : alignment positions), that is stored in the files.  In
this example, files `foo/0:0,46:1653469.*` store records whose (refid:pos) are
in range [(0:0), (46:1653469)).  Files `foo/46:1653469,-:-.*` store records in
coord range [(46:1653469), (∞:∞)). The refid of "-" means "unmapped", and position
of "-" means ∞.

    Note: `[a,b)' means a half-open range, inclusive of `a`, exclusive of `b`.

PAM always sorts reads according to BAM's "coordinate" sort order. That is,
records are reference ID, then position, then forward reads before reverse
reads.  Ordering of records with the same (refid:position) is unspecified.
There is no option to sort by sequence names.  PAM file-shard boundary is always
at a position boundary. That is, a set of records with the same (refid:position)
are always in one shard. In particular, unmapped sequences are always in the
last shard.

With a well-formed set of PAM files, the coordranges of the shardfiles cover the
universal range `0:0,-:-` without a gap or an overlap.

# Index file

File `foo/0:0,46:1653469.index` is an index file for the coordinate
range [(0:0), (46:1653469)).
practice, index stores an equivalent of a `sam.Header`:

```
// Contents of the fileshard index (*.index) file. It is stored in a flate-compressed,
// recordio block.
message ShardIndex {
  // A magic number. Always 0x725c7226be794c60
  fixed64 magic = 1;

  // As of 2017-10-18, always "PAM1"
  string version = 3;

  // Range of records. This is the same as the record range encoded in the
  // filename. The range of records actually stored the data files may be a
  // subset of this range.
  RecRange range = 4 [(gogoproto.nullable) = false];

  // sam.Header encoded in BAM format.
  bytes encoded_bam_header = 15
      [(gogoproto.nullable) = false, (gogoproto.customtype) = "SAMHeader"];
}
```

# Field data file

Files of form `dir/coordrange.field` store data for the given field. There are
eleven fields: coord, flags, mapq, cigar, materefid, matepos, templen, name,
seq, qual, aux.  The "coord" file stores the (refid,pos) of each record, and
other field files store the values of the field with the same name in
sam.Record.  Each data file is a recordio file
(https://github.com/grailbio/base/tree/master/recordio), with 8MB
pre-compression block size, and using zstd for compression.

The field data files for a given coordinate range always store exactly the same
number of records. However, the recordio block boundaries aren't necessarily
aligned across fields, because values for some fields (e.g., `seq` or `qual`)
are larger than others' (e.g., `mapq` or `templen`).  For example, (say) 10th
recordio block for `foo/0:0,46:1653469.mapq` may store reads in range
[(0:10000),(0:12345)], but the 10th recordio block for `foo/0:0,46:1653469.seq`
may store reads in range [(0:1234),(0:4567)].

## Field data file format

A field data file (e.g., `foo.pam.1:0,46:1653469.mapq.data`) is stored as a recordio file
(https://github.com/grailbio/base/tree/master/recordio).

As the PAM writer receives records to write, it extracts their field values and
appends them to the eleven buffers, one for each field.  Once the size of a
field buffer a limit (8MB, pre-compression), the buffer is compressed and
written to the corresponding recordio file. Each (8MB-pre-compression) buffer
becomes one recordio block.  Each recordio block has the following layout:

             Block header length (4 bytes varint)
             Block header        (serialized BlockHeader proto)
             Subfield 0
             Subfield 1     // optional

The block header is a serialized `BlockHeader` as defined below.  It encodes the
location of each field within the recordio block.

```
message BlockHeader {
  // Field type stored in the block.
  int32 field = 1 [(gogoproto.casttype) = "FieldType"];

  // Locations of the default and blob subfields. The values are the number of
  // bytes from the end of the BlockHeader, pre compression.
  uint32 offset = 2;
  uint32 blob_offset = 3;
}
```

### Subfields

The rest of the block contains values for the field.  Some field value may be
split into up to two components called _subfields_, and are stored in
subfield-major order. More specifically:

- A simple numeric field, such as `sam.Record.Flags` or `sam.Record.Mapq` are
  stored in field-major order naively.

- A variable-length field, such as `sam.Record.Qual` or `sam.Record.Aux`, is
  split into two subfields.  One subfield, the *default subfield*, stores the
  length of the field. Another subfield, *blob0* stores the actual quality data
  as a byte sequence.  The default subfield is stored before the blob0 in the
  recordio block.

  > Note: the default subfield is used to store simple numeric field values.

- The coord field stores sam.Record.Ref.ID() is the default field, and
  sam.Record.Pos in the blob0 subfield.

To give a simple example, assume that a block contains three fields, RefID, Pos
and Qual from two records. We ignore other fields, such as Seq or Aux in this
example.

    R0: {RefID=1, Pos=100, Qual=[41,39,41,41])
    R1: {RefID=1, Pos=101, Qual=[41,41,39,39])

Data file for coord field:

    1 1 // refid
    100 101  // pos

Data file for qual field:

    4 4  // length of qual fields
    41 39 41 41 41 41 39 39  // concat of qual field values

The subfield-major layout results in a highly repetitive data pattern
pre-compression, allowing the compressor to do a more effective job.

### Encoding of field types

- RefID, MateRefID, MatePos: A varint diff from the previous record.  The diff
  is reset at a recordio block boundary. The value is stored in the
  default subfield.

  Imagine we store records r0, r1, r2 in one recordio block. Then, for
  MateRefID, we store values r0.MateRefID, (r1.MateRefID-r0.MateRefID),
  (r2.MateRefID-r1.MateRefID) in the materef data field file. The reason we
  store diffs instead of absolute values is that the diffs are often zero
  or otherwise small value, and they will encode smaller and compress better.

- Pos: A varint diff from the previous record.  The diff is reset at a recordio
  block boundary. The value is stored in the blob subfield.  This is because
  (RefID, Pos) are stored in one file, with RefID going to the default subfield
  and pos going to the blob subfield.

- TemplateLen: varint (default subfield)
- Flags: fixed uint16 (default subfield)
- MapQ: fixed uint8 (default subfield)

- Name: prefix encoded in the following fashion:
  * default subfield: length of prefix shared with the previous record
  * blob0: the suffix that differs from the previous record. The diff is reset
    at a recordio(column) block boundary.

- Cigar:
  * default subfield: cigar length
  * blob0:  uvarint encoding of (length<<4 |  optype)

- Seq:
  * default subfield: sequence length
  * blob0:  sequence of bytes (4 bits per base)

- Qual:
  * default subfield: qual length
  * blob0:  sequence of bytes (8 bits per base)

- Aux:
  * Encode # of aux tags as a varint in the default subfield
  * For each aux tag:
      - Encode the first three bytes (the two-byte tag type + one-byte type) of
        the tag in the blob subfield.  We use the type fields defined in
        github.com/biogo/hts/sam/auxtags.go. It is an extension of SAM's.

  * For each aux tag:
      - For a fixed-length field, such as 'i', 'f', etc, encode value in binary
        in the blob field.  For a variable-length field, i.e., 'Z' and 'H',
        encode the value length (excluding the first three bytes) in the default
        subfield, then encode the value in the blob field.

  For example, tags `AS:i:96 RG:Z:12345` will be encoded in the following fashion:

```
  default subfield:
    2  // number of tags, as varint
    5  // length of the RG:Z payload, as varint

  blob0 subfield:
    'A', 'S', 'i', 'R', 'G', 'Z'  // sequence of tag names+types
    96                            // "AS:i" payload as binary byte
    '1', '2', '3', '4', '5'       // payload for RG:Z tag
```

### Field data index

Each field-data file stores an index in the recordio trailer
(https://github.com/grailbio/base/tree/master/recordio) part.

    Note: This index is
    unrelated to the fileshard index file (foo/coord.index) described in the
    previous section.

The index is a serialized `FieldIndex` as defined in `pam.proto`:

```
// BlockIndexEntry summarizes contents of one block. It is stored in the
// ShardIndex.
message BlockIndexEntry {
  // Offset of the start of the block within the recordio data file.  It is
  // passed directly to File.Seek().
  uint64 file_offset = 1;
  // Pre-compression size of the recordio block.
  uint32 uncompressed_bytes = 2;

  uint32 num_records = 3;  // Number of sam.Records stored

  // The refid and alignment position of the first record stored in the block.
  RecAddr start_addr = 4 [(gogoproto.nullable) = false];

  // The refid and alignment position of the last record stored in the block.
  // Note: unlike RecRange, [startAddr, endAddr] is closed, because we don't
  // know the open limit addr when flushing a recordioblock. Use
  // blockIntersectsRange to check if [startAddr,endAddr] intersects a RecAddr.
  RecAddr end_addr = 5 [(gogoproto.nullable) = false];
}

// RecAddr uniquely identifies a sam.Record in a PAM file.
//
// Field seq is used to distinguish reads at the same coordinate (ref_id,
// pos). The first PAM record at a particular (ref_id, pos) has seq=0, the
// second record at the same (ref_id, pos) has seq=1, and so on.  sam.Record
// does not store the seq value. It is computed by PAM reader and writer during
// I/O.
//
// The total ordering of RecAddr is defined (1) by increasing order of ref_id,
// then (2) by increasing order of pos, then (3) by increasing order of seq.
// However, there is one exception: records with ref_id=-1 are sorted after any
// other records.  Records with ref_id=-1 are by SAM/BAM convention unmapped,
// and they are stored after mapped reads.
//
// For unmapped reads, pos is meaningless. We use address ref_id=-1,pos=0 for
// any unmapped records (see RecAddrFromSAMRecord).  For a RecRange for unmapped
// reads, we use {{-1,0,0},{-1,int32max,int32max}}.
//
// TODO(sits) This convention is a bit hairy. We could instead define the
// universal range to be {{0,int32min},{-1,int32max}}. This has a downside that
// using 0 as the min position causes the code to misbehave.
message RecAddr {
  int32 ref_id = 1;
  int32 pos = 2;
  int32 seq = 3;
}

// RecRange is a half-open range of sam.Records. The start bound is closed, the
// limit is open.
//
// Examples:
//   RecRange{{0,0,0},{InfinityRefID, InfinityPos, 0}} : covers all possible
//   sequences.
//   RecRange{{0,0,0},{MaxValidRefID, InfinityPos, 0}} : covers all mapped sequences.
//   RecRange{{UnmappedRefID,0,0},{UnmappedRefID, InfinityPos,0}} : covers all
//   unmapped sequences.
//
// INVARIANT: [(start.ref_id, start.pos), (limit.ref_id, limit.pos)) must
// represent an nonempty interval.
message RecRange {
  RecAddr start = 1 [(gogoproto.nullable) = false];
  RecAddr limit = 2 [(gogoproto.nullable) = false];
}

// Contents of the field index (*.<fieldname>.index) file. It is stored in a
// flate-compressed, recordio block.
message FieldIndex {
  // A magic number. Always 0x725c7226be794c60
  fixed64 magic = 1;

  // As of 2017-10-18, always "PAM1"
  string version = 3;

  int32 field = 4 [(gogoproto.casttype) = "FieldType"];

  // Stores one entry per recordio block. Sorted by RecAddrs.
  repeated BlockIndexEntry blocks = 16 [(gogoproto.nullable) = false];
}
```

# Comparison to BAM

## File size

As an example a 113GB BAM file,
s3://grail-avalon/pipeline/clinical/v2/bams/CNVS-NORM-110033752-cfDNA-WGBS-Rep1.bam,
becomes 69GB, or 57%, in PAM. A later section discusses effects of different
compression algorithms.

+ BAM: 122012650919 bytes

- PAM with default compression (zstd 3): 69920522240 bytes (57% of BAM)

    Cost of conversion from BAM⟶PAM:
    real	16m18.042s
    user	301m7.572s
    sys	16m29.188s

- PAM with "zstd 6" compression: 67525312512 bytes (55% of BAM)

    Cost of conversion from BAM⟶PAM:
    real	17m29.566s
    user	374m55.164s
    sys	16m54.620s

- PAM with "zstd 10" compression: 67525312512 bytes (53% of BAM)

    Cost of conversion from BAM⟶PAM:
    real	23m20.641s
    user	673m42.464s
    sys	17m51.464s

- PAM with "zstd 20" compression: 57791840256 bytes (47% of BAM)

    Cost of conversion from BAM⟶PAM:
    real	150m36.767s
    user	6637m9.440s
    sys	624m5.760s

## Performance

Creating PAM files for CNVS-NORM-110033752-cfDNA-WGBS-Rep1.bam takes about 980s,
whereas creating the original BAM file using the biogo BAM writer interface
would have taken about 1800-2400s (this is guesstimate).

Reading the same BAM/PAM files (on ubuntu02.mpk, with 56 cpus):

    - BAM sharded reader in Go: 900s
    - PAM sharded reader: 470s
    - PAM sharded reader, but dropping Name and Qual fields: 330s

> Note: This result looks a bit strange: it takes longer to read a BAM than to
> convert PAM to BAM, which involves reading the same BAM file.  This happens
> because the BAM sharded reader creates very uneven shards.  We could use the
> same internal sharder used by the BAM-to-PAM converter and the BAM’s
> performance would improve to something much closer to PAM.

Reading the same BAM/PAM files excluding unmapped reads (on ubuntu02.mpk):

    - BAM sharded reader in Go: 310s
    - PAM sharded reader: 130s
    - PAM sharded reader, but dropping Name and Qual fields: ???s

> Note: The huge difference between unmapped reads and no-unmapped reads is that
> reading of unmapped segments is non-parallelizable in the current
> implementation. At least in PAM, we have fine-grained index for the unmapped
> segment, so we can potentially parallelize reading it.

## Downsides

Random seeking to a specific PAM record is slow. To read one record, the
application needs to open eleven data files, one per field. For each data file,
it needs to seek to the recordio block containing the record, read and
uncompress the block and extract one value.

In comparison, BAM requires opening one file, reading and gunzipping one 16KB
block.

# Go API

We currently only offer Go API. We might offer a C++ API if a need arises.

## pam.Writer

```
package pam

type WriteOpts struct {
	// Compression block size, by default 8MB
    MaxBufSize          int
    // Number of write ahead allowed, by default 16.
	MaxParallelism int

	// WriteParallelism limits the max number of pending recordio flushes
	// allowed. If <= 0, DefaultWriteParallelism is used.
	WriteParallelism int

	// DropFields causes the writer not to write the specified fields to file.
	DropFields []FieldType

	// Transformers defines the recordio block transformers. It can be used to
	// change the compression algorithm, for example. The value is passed to
	// recordio.WriteOpts.Transformers. If empty, {"zstd"} is used.
	Transformers []string

	// Range defines the range of records that can be stored in the PAM
	// file.  The range will be encoded in the path name. Also, Write() will
	// cause an error if it sees a record outside the range. An empty range
	// (default) means UniversalRange.
	//
	// The range bound is closed at the start, open at the limit.
	Range RecRange
}


// NewWriter creates a new PAM writer. Files are created in "dir". If "dir" does
// not exist already it is created. Existing contents of "dir", if any, are
// deleted.
func NewWriter(wo WriteOpts,
    samHeader *sam.Header,
    dir string) *Writer

// Write a new record.
//
// REQUIRES: the record must be in WriteOpts.Range.
func (w *Writer) Write(r *sam.Record)

// Report any error encountered so far
func (w *Writer) Err() error

// Close all the files. Must be called exactly once.
func (w *Writer) Close() error
```

## pam.Reader

Example:

```
   r := pam.NewReader(...)
   for r.Scan() {
      rec := r.Record()
      ... use the rec ...
   }
   err := r.Close()
```

```
package pam

type ReadOpts struct {
	// DropFields causes the listed fields not to be filled in Read().
	DropFields []FieldType

    // Coordinate range to read.
	Range RecRange
}

// Create a new reader. The reader is thread compatible.
func NewReader(opts ReadOpts, pathPrefix string) *Reader

// Scan reads the next record. It returns true if a record is found. It returns
// false on EOF or an error. Call Record to get the record, and Err or Close to
// obtain the error code.
func (r *Reader) Scan() bool

// Record returns the most recent record read by Scan.
//
// REQUIRES: Scan() returned true.
func (r *Reader) Record() *sam.Record

// Err returns any error encountered so far.
//
// Note: Err never returns io.EOF. On EOF, Scan() returns false, and Err()
// returns nil.
func (r *Reader) Err() error

// Close all the files. Must be called exactly once.
func (r *Reader) Close() error

// GenerateReadShards returns a list of RecRanges. The RecRanges can be passed
// to NewReader for parallel, sharded record reads. The returned list satisfies
// the following conditions.
//
// 1. The ranges in the list fill opts.Range (or the UniversalRange if not set)
//    exactly, without an overlap or a gap.
//
// 2. Length of the list is at least nShards. The length may exceed nShards
//    because this function tries to split a range at a fileshard boundary.
//
// 3. The bytesize of the file region(s) that covers each RecRange is roughly
// the same.
//
// 4. The ranges are sorted in an increasing order of RecAddr.
func GenerateReadShards(opts ReadOpts, path string, nShards int) ([]RecRange, error) {
```

## Future extensions

### Adding annotations

We want to allow attaching extra aux tags, or possibly fields to records in an
existing PAM file.  An annotation is stored as a small file that shares path
prefix with the original PAM file. The original PAM file remains intact. Tools
such as [bio-mark-duplicates](../../../cmd/bio-mark-duplicates) will create
small annotation files instead of duplicating the whole PAM file.

The below snippet shows an example of adding a new Aux field, `D`, to every
record in "foo.pam".

```
r := pam.NewReader(opts, "foo.pam")
w := pam.NewAnnotationWriter("foo.pam", "duplicates")
for r.Scan() {
  // Add an annotation to "rec".
  ann := pam.CloneRequiredFields(r.Record())
  ann.Aux = []sam.Aux{"D:I:123"}
  w.Write(ann)
}
if err := r.Close(); err != nil { panic(err) }
if err := w.Close(); err != nil { panic(err) }
```

Below is a snippet for reading the new annotation.

```
r := pam.NewReader(pam.ReadOpts{Annotations: []string{"duplicates"}}, "foo.pam")
for r.Scan() {
  rec := r.Record()
  // rec.Aux will have "D:I:123", along with any other annotation found in the
  // base file.
}
if err := r.Close(); err != nil { panic(err) }
if err := w.Close(); err != nil { panic(err) }
```

An annotation will be encoded mostly like a regular PAM file. When reading the
annotation, its contents are merged into records read from the main PAM
file. More specifically,

- The records stored in an annotation may not have all the sam.Record fields.
  Typically, the file only stores fields (RefID, Pos, Name, Flag, Aux). The
  first four are necessary to uniquely identify and sort the records. The Aux
  field will be appended to the Aux field in the original record.

  While Aux is the field most often stored in the annotation, it can store other
  fields to augment or replace fields in the base PAM file.

- The annotation file need not list all the records found in the base PAM file.
  The recordio block boundaries need not match those in the base PAM file.

- We may provide a way to delete a field or record from the base PAM file, but
  this is a lower-priority feature.

### Adding new fields

The same idea can be extended to add a field to records in the base file.  The
internal machinery is virtually the same as Aux fields, but the annotation will
be stored in a new field in sam.Record. This may lead to more readable or
compact data and code. The downside of adding a new field is that we need to add
a field in sam.Record apriori, and have people agree on the name and type of it.

### Paired reads

Many applications that consume BAM/PAM want to read sequences in pairs
(R1+R2). Currently, we have many different abstractions, such as fragment
files, `PairIterators`, to create pairs from a BAM/PAM file. We
want to embed pairs more explicitly in the PAM file.

The code for writing pairs will look like below.

```
w := pam.NewWriter(pam.WriteOpts{}, "foo.pam")
for somecondition {
  r1 := record to write
  r2 := mate record
  w.WritePair(r1, r2)
}
```

The code for reading pairs will look like below.

```
r := pam.NewReader(pam.ReadOpts{}, "foo.pam")
for r.Scan() {
  r1, r2 := r.Pair()
  .. use r1 and r2 ..
}
```

`r.Pair` can be used only when the PAM file was written using `WritePair`.  On
the other hand, it's ok to use non-paired reads for a file written using
`WritePair`.

Internally, the writer will perform the following things:

- Sort R1 and R2 in increasing order of positions.
- If both R1 and R2 fit in one record block, store them normally.
- If they do not fit in one recordio block, the records are replicated, with a
  special flag value to indicate that they are replicas.
- If the mates span across rowshard boundaries, the caller must call WritePair
  in both shards.
