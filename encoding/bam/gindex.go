package bam

import (
	"bytes"
	"encoding/binary"
	"fmt"
	"io"
	"sort"

	"github.com/biogo/hts/bgzf"
	"github.com/biogo/hts/sam"
	"github.com/klauspost/compress/gzip"
)

const (
	maxRecordSize = 0xffffff
)

// GIndex is an alternate .bam file index format that uses the .gbai
// file extension.  The .gbai file format contains mappings from
// genomic position to .bam file voffset.  This index format is
// simpler than the legacy style .bai file, but allows a user to seek
// into a .bam file much more efficiently for some genomic positions.
//
// The .gbai format exists because the .bai format can point into a
// .bam file with a minimum genomic spacing of 16 kbp.  The problem
// with this minimum spacing is that if there are many alignments in
// the .bam file within a 16 kbp region, then seeking to a target
// genomic position within the 16 kbp region requires the reader to
// seek to the beginning for the 16 kbp region and then scan through
// bam records until reaching the target genomic position.  This
// scanning requires unnecessary IO and CPU time for reading and
// decompressing records that come before the target position.
//
// The .gbai file format contains a set of mappings from (genomic
// position, and record number at that position) to the voffset in the
// bam file where the record begins.  In typical use, the spacing
// between the genomic positions in the .gbai file are chosen so that
// the spacing between voffsets in the .bam file are uniform and
// relatively small.  This allows a user to divide the .bam file into
// uniform sized shards.  For example, 64 KBytes is a reasonable
// default spacing between voffsets.  This spacing allows a reader to
// seek directly to within 64 KBytes of any target genomic position.
//
// The on disk .gbai format is a header followed by a sequence of
// entries.  The header consists of the magic byte sequence
// {0x47, 0x42, 0x41, 0x49, 0x01, 0xf1, 0x78, 0x5c,
//  0x7b, 0xcb, 0xc1, 0xba, 0x08, 0x23, 0xb1, 0x19}
// which is "GBAI1" followed by 11 random bytes.
//
// Each entry consists of 4 values, each in little-endian byte order:
//   1) int32 RefID to match the .bam file RefIDs.  The unmapped
//      records at the end of the .bam have RefID equal to -1.
//   2) int32 Position to match the .bam file Positions
//   3) uint32 Sequence number of the record at the particular (RefID,
//      Position) pair.  If the record is the first record with this
//      (RefID, Position) pair, then Sequence will be 0.  If the
//      record is the second, then Sequence will be 1, and so on.
//   4) uint64 VOffset of the record in the .bam file as described in
//      the .bam specification.
//
// The .gbai index entries are sorted in ascending order using the key
// (RefID, Position, Sequence) and the .gbai index requires that the
// corresponding .bam file is also sorted by position.
//
// If the bam file contains a bam record for a given RefID, then the
// gindex contains an entry for the first bam record with the given
// RefID.  This implies that the first entry in the gindex points to
// the first record in the bam file.  If there are no bam records with
// RefID R, then there will be no entries in the gindex with RefID R.
//
// The series of index entries is then compressed with gzip before
// writing to the .gbai file.
type GIndex []GIndexEntry

var gbaiMagic = []byte{
	'G', 'B', 'A', 'I', 0x01, 0xf1, 0x78, 0x5c,
	0x7b, 0xcb, 0xc1, 0xba, 0x08, 0x23, 0xb1, 0x19,
}

// GIndexEntry is one entry of the .gbai index.
type GIndexEntry struct {
	RefID   int32
	Pos     int32
	Seq     uint32
	VOffset uint64
}

// RecordOffset returns a voffset into the bam from which, reading
// forward will eventually read records at the target position.  When
// reading from the returned voffset, if the bam record's (refid,
// position) is greater than the target (refid, position), then the
// target position is not present in the bam file.
func (idx *GIndex) RecordOffset(refID, pos int32, seq uint32) bgzf.Offset {
	if len(*idx) < 1 {
		panic("GIndex must have at least one entry")
	}
	target := GIndexEntry{refID, pos, seq, 0}
	x := sort.Search(len(*idx), func(i int) bool {
		return comparePos(&(*idx)[i], &target) >= 0
	})

	if x == len(*idx) {
		return ToBGZFOffset((*idx)[x-1].VOffset)
	}

	// If search returned an entry that is larger than target, then
	// try to back up one entry.
	if comparePos(&(*idx)[x], &target) > 0 {
		if x > 0 {
			x--
		}
	}
	return ToBGZFOffset((*idx)[x].VOffset)
}

// UnmappedOffset returns a voffset at or before the first read in the
// .bam's unmapped section.
func (idx *GIndex) UnmappedOffset() bgzf.Offset {
	return idx.RecordOffset(-1, 0, 0)
}

// ToBGZFOffset takes a uint64 voffset and returns a bgzf.Offset.
func ToBGZFOffset(voffset uint64) bgzf.Offset {
	return bgzf.Offset{int64(voffset >> 16), uint16(voffset & 0xffff)}
}

// toVOffset takes a bgzf.Offset and returns a uint64 voffset.
func toVOffset(offset bgzf.Offset) uint64 {
	return uint64(offset.File)<<16 | uint64(offset.Block)
}

// gIndexWriter writes a .gbai index file.
type gIndexWriter struct {
	gz *gzip.Writer
}

// newGIndexWriter creates a new gIndexWriter.
func newGIndexWriter(w io.Writer) *gIndexWriter {
	return &gIndexWriter{
		gz: gzip.NewWriter(w),
	}
}

func (w *gIndexWriter) writeHeader() error {
	n, err := w.gz.Write(gbaiMagic)
	if err != nil {
		return err
	}
	if n != len(gbaiMagic) {
		return fmt.Errorf("short write to gbai header: %d should be %d", n, len(gbaiMagic))
	}
	return nil
}

func (w *gIndexWriter) append(entry *GIndexEntry) error {
	return binary.Write(w.gz, binary.LittleEndian, entry)
}

func (w *gIndexWriter) close() error {
	return w.gz.Close()
}

func comparePos(x, y *GIndexEntry) int {
	if x.RefID != y.RefID {
		if x.RefID < 0 && y.RefID >= 0 {
			return 1
		} else if x.RefID >= 0 && y.RefID < 0 {
			return -1
		}
		return int(x.RefID) - int(y.RefID)
	}
	if x.Pos > y.Pos {
		return 1
	} else if x.Pos < y.Pos {
		return -1
	}

	if x.Seq > y.Seq {
		return 1
	} else if x.Seq < y.Seq {
		return -1
	}
	return 0
}

func compareFilePos(x, y *GIndexEntry) int {
	return int(int64(x.VOffset) - int64(y.VOffset))
}

// WriteGIndex reads a .bam file from r, and writes a .gbai file to w.
// The spacing between voffset file locations will be approximately
// byteInterval, and parallelism controls the .bam file read
// parallelism.  Currently, WriteGIndex will not create two index
// entries for a given (RefID, Pos) pair, i.e. Seq will always be
// zero.  That means there will be only one entry for the entire
// unmapped region.
func WriteGIndex(w io.Writer, r io.Reader, byteInterval, parallelism int) error {
	bgzfReader, err := bgzf.NewReader(r, parallelism)
	if err != nil {
		return err
	}
	header, err := sam.NewHeader(nil, nil)
	if err != nil {
		return err
	}
	if err := header.DecodeBinary(bgzfReader); err != nil {
		return err
	}
	gindex := newGIndexWriter(w)
	if err := gindex.writeHeader(); err != nil {
		return err
	}

	// Read through all alignment records, and output voffsets at shard boundaries.
	prevRefID := int32(0)
	prevPos := int32(0)
	prevFileOffset := uint64(0)

	sizeBuf := make([]byte, 4)
	buf := make([]byte, maxRecordSize)
	firstRecord := true

	for {
		// Read the record size.
		_, err := io.ReadFull(bgzfReader, sizeBuf)
		if err == io.EOF {
			break
		} else if err != nil {
			return err
		}
		recordVOffset := bgzfReader.LastChunk().Begin
		sz := int(binary.LittleEndian.Uint32(sizeBuf))
		if sz > maxRecordSize {
			return fmt.Errorf("bam record exceeds max: %d", sz)
		}
		_, err = io.ReadFull(bgzfReader, buf[0:sz])
		if err == io.EOF {
			return fmt.Errorf("could not read full bam record")
		} else if err != nil {
			return err
		}

		// Parse the relevant fields
		refID := int32(binary.LittleEndian.Uint32(buf[0:4]))
		pos := int32(binary.LittleEndian.Uint32(buf[4:8]))

		// Always add an entry for the first record of a new RefID.
		if firstRecord || refID != prevRefID {
			entry := GIndexEntry{
				RefID:   refID,
				Pos:     pos,
				Seq:     0,
				VOffset: toVOffset(recordVOffset),
			}
			if err := gindex.append(&entry); err != nil {
				return err
			}
			prevRefID = refID
			prevPos = pos
			prevFileOffset = uint64(recordVOffset.File)
			firstRecord = false
			continue
		}

		firstOccurrence := false
		if pos != prevPos {
			prevPos = pos
			firstOccurrence = true
		}
		if firstOccurrence && (uint64(recordVOffset.File)-prevFileOffset) >= uint64(byteInterval) {
			entry := GIndexEntry{
				RefID:   refID,
				Pos:     pos,
				Seq:     0,
				VOffset: toVOffset(recordVOffset),
			}
			if err := gindex.append(&entry); err != nil {
				return err
			}
			prevFileOffset = uint64(recordVOffset.File)
		}
	}
	return gindex.close()
}

// ReadGIndex expects a .gbai file as r, and returns the parsed GIndex
// and any errors encountered while reading and unmarshalling r.
func ReadGIndex(r io.Reader) (gindex *GIndex, err error) {
	var gz *gzip.Reader
	gz, err = gzip.NewReader(r)
	if err != nil {
		return nil, err
	}
	defer func() {
		if cerr := gz.Close(); cerr != nil && err != nil {
			err = cerr
		}
	}()

	buf := make([]byte, len(gbaiMagic))
	if _, err = io.ReadFull(gz, buf); err != nil {
		return nil, err
	}

	if !bytes.Equal(gbaiMagic, buf) {
		return nil, fmt.Errorf("Unexpected gbai magic: %v should be %v", buf, gbaiMagic)
	}

	index := make(GIndex, 0)
	for i := 0; ; i++ {
		entry := GIndexEntry{}
		if err = binary.Read(gz, binary.LittleEndian, &entry); err == io.EOF {
			break
		} else if err != nil {
			return nil, err
		}

		if i > 0 {
			prev := index[i-1]
			if comparePos(&prev, &entry) >= 0 {
				return nil, fmt.Errorf("Index positions are out of order %v must be less than %v", prev, entry)
			}
			if compareFilePos(&prev, &entry) >= 0 {
				return nil, fmt.Errorf("Voffsets are out of order %v must be less than %v", prev, entry)
			}
		}

		// It is possible for some reference IDs to be skipped if there are no records in that Ref.
		index = append(index, entry)
	}
	return &index, nil
}
