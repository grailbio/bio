package bam

import (
	"fmt"

	"github.com/grailbio/hts/sam"
)

// IsPaired returns true if record is paired.
func IsPaired(record *sam.Record) bool {
	return record.Flags&sam.Paired != 0
}

// IsProperPair returns true if record is properly aligned.
func IsProperPair(record *sam.Record) bool {
	return record.Flags&sam.ProperPair != 0
}

// IsUnmapped returns true if record is unmapped.
func IsUnmapped(record *sam.Record) bool {
	return record.Flags&sam.Unmapped != 0
}

// IsMateUnmapped returns true if mate of record is unmapped.
func IsMateUnmapped(record *sam.Record) bool {
	return record.Flags&sam.MateUnmapped != 0
}

// IsReverse returns true if record maps to reverse strand.
func IsReverse(record *sam.Record) bool {
	return record.Flags&sam.Reverse != 0
}

// IsMateReverse returns true if mate of record maps to reverse strand.
func IsMateReverse(record *sam.Record) bool {
	return record.Flags&sam.MateReverse != 0
}

// IsRead1 returns true if record is first in pair.
func IsRead1(record *sam.Record) bool {
	return record.Flags&sam.Read1 != 0
}

// IsRead2 returns true if record is second in pair.
func IsRead2(record *sam.Record) bool {
	return record.Flags&sam.Read2 != 0
}

// IsSecondary returns true if record is a secondary alignment.
func IsSecondary(record *sam.Record) bool {
	return record.Flags&sam.Secondary != 0
}

// IsQCFail returns true if record does not pass quality control filters.
func IsQCFail(record *sam.Record) bool {
	return record.Flags&sam.QCFail != 0
}

// IsDuplicate returns true if record is a duplicate.
func IsDuplicate(record *sam.Record) bool {
	return record.Flags&sam.Duplicate != 0
}

// IsLinearDuplicate returns true if record is a linear duplicate.
func IsLinearDuplicate(record *sam.Record) bool {
	dupType, err := record.LinearDup()
	if err != nil {
		panic(err)
	}
	return dupType == sam.LinearDuplicate
}

// IsSupplementary returns true if record is a supplementary alignment.
func IsSupplementary(record *sam.Record) bool {
	return record.Flags&sam.Supplementary != 0
}

// IsPrimary returns true if record is a primary alignment.
func IsPrimary(record *sam.Record) bool {
	return record.Flags&sam.Secondary == 0 && record.Flags&sam.Supplementary == 0
}

// HasNoMappedMate returns true if record is unpaired or has an unmapped mate.
func HasNoMappedMate(record *sam.Record) bool {
	return (record.Flags&sam.Paired) == 0 || (record.Flags&sam.MateUnmapped) != 0
}

// ClearAuxTags removes all instances of the tags in tagsToRemove[] from r.
// (Current implementation is not designed for very large len(tagsToRemove); at
// some point map lookups become better.)
func ClearAuxTags(r *sam.Record, tagsToRemove []sam.Tag) {
	found := 0
	for i := range r.AuxFields {
		tag := r.AuxFields[i].Tag()
		for _, toRemove := range tagsToRemove {
			if tag == toRemove {
				found++
			}
		}
	}

	if found > 0 {
		newAux := make([]sam.Aux, len(r.AuxFields)-found)
		size := 0
		for i := range r.AuxFields {
			tag := r.AuxFields[i].Tag()
			keepAux := true
			for _, toRemove := range tagsToRemove {
				if tag == toRemove {
					keepAux = false
					break
				}
			}
			if keepAux {
				newAux[size] = r.AuxFields[i]
				size++
			}
		}
		r.AuxFields = newAux
	}
}

type StrandType int

const (
	// StrandNone denotes an unmapped read.
	StrandNone StrandType = iota
	// StrandFwd denotes a read mapped to the R1+ strand.
	StrandFwd
	// StrandRev denotes a read mapped to the R1- strand.
	StrandRev
)

// GetStrand returns whether the current read is mapped to the R1+ strand, the
// R1- strand, or neither, in a manner that ignores all flags associated with
// the other read end.
func GetStrand(flags sam.Flags) StrandType {
	maskedFlag := flags & (sam.Unmapped | sam.Reverse | sam.Read1 | sam.Read2)
	if (maskedFlag == sam.Read1) || (maskedFlag == (sam.Reverse | sam.Read2)) {
		return StrandFwd
	}
	if (maskedFlag == sam.Read2) || (maskedFlag == (sam.Reverse | sam.Read1)) {
		return StrandRev
	}
	return StrandNone
}

// IsReversedRead returns true if the reverse flag is set on record.
func IsReversedRead(record *sam.Record) bool {
	return (record.Flags & sam.Reverse) != 0
}

// IsQCFailed returns true if the QC failed flag is set on record.
func IsQCFailed(record *sam.Record) bool {
	return (record.Flags & sam.QCFail) != 0
}

// LeftClipDistance returns the total amount of clipping (both hard
// and soft) on the left-most side of record.
func LeftClipDistance(record *sam.Record) int {
	total := 0
	for _, op := range record.Cigar {
		if op.Type() == sam.CigarSoftClipped || op.Type() == sam.CigarHardClipped {
			total += op.Len()
		} else {
			break
		}
	}
	return total
}

// RightClipDistance returns the total amount of clipping (both hard
// and soft) on the right-most side of record.
func RightClipDistance(record *sam.Record) int {
	total := 0
	for i := len(record.Cigar) - 1; i >= 0; i-- {
		op := record.Cigar[i]
		if op.Type() == sam.CigarSoftClipped || op.Type() == sam.CigarHardClipped {
			total += op.Len()
		} else {
			break
		}
	}
	return total
}

// FivePrimeClipDistance returns the total amount of clipping (both
// hard and soft) on the 5' side of record.
func FivePrimeClipDistance(record *sam.Record) int {
	if IsReversedRead(record) {
		return RightClipDistance(record)
	}
	return LeftClipDistance(record)
}

// UnclippedStart returns the unclipped left-most position of record,
// regardless of record's read direction.
func UnclippedStart(record *sam.Record) int {
	return record.Start() - LeftClipDistance(record)
}

// UnclippedStart returns the unclipped right-most position of record,
// regardless of record's read direction.
func UnclippedEnd(record *sam.Record) int {
	unclippedEnd := record.End() + RightClipDistance(record) - 1
	if unclippedEnd < record.Start() {
		panic(fmt.Sprintf("unclippedEnd is less than start: %d < %d for %s",
			unclippedEnd, record.Start(), record.Name))
	}
	return unclippedEnd
}

// UnclippedFivePrimePosition returns the unclipped 5' position of
// record.
func UnclippedFivePrimePosition(record *sam.Record) int {
	if IsReversedRead(record) {
		return UnclippedEnd(record)
	}
	return UnclippedStart(record)
}

// indexAtPos returns an index into the sequence and quality string in
// the record that matches the pos (0 based) argument.  If the read
// overlaps pos, and the cigar maps a base to the reference position,
// then return (sequenceIndex, true).  If the read overlaps the
// position, but the cigar does not map a base to the reference
// position, then return (-1, true).  If the read does not overlap the
// reference pos, then return (-1, false).
func indexAtPos(record *sam.Record, pos int) (int, bool) {
	if pos < record.Start() || pos >= record.End() {
		return -1, false
	}

	refPos := record.Pos
	seqPos := 0

	for _, co := range record.Cigar {
		t := co.Type()
		if t == sam.CigarBack {
			panic("No support for CigarBack")
		}

		consume := t.Consumes()
		lenref := co.Len()

		if refPos <= pos && refPos+lenref > pos {
			if consume.Query == 1 && consume.Reference == 1 {
				return seqPos + (pos - refPos), true
			} else if consume.Query == 1 {
				seqPos += lenref
			} else if consume.Reference == 1 {
				return -1, true // pos on the reference was skipped.
			}
		} else {
			if consume.Query == 1 {
				seqPos += lenref
			}
			if consume.Reference == 1 {
				refPos += lenref
			}
		}
	}
	panic("Unexpected exit of loop")
	return -1, false
}

// BaseAtPos returns the base at reference pos (0 based) from record,
// if the mapped part of the read overlaps pos.  If not, return (0,
// false).  If pos is in the mapped portion of the read, but the
// reference base was skipped, then the returned value will be (0,
// true).
func BaseAtPos(record *sam.Record, pos int) (byte, bool) {
	index, overlap := indexAtPos(record, pos)
	if !overlap {
		return 0, false
	} else if index >= 0 {
		return record.Seq.BaseChar(index), true
	}
	return 0, true
}

// QualAtPos returns the base quality byte at reference pos (0 based)
// from record, if the mapped part of the read overlaps pos.  If not,
// return (0, false).  If pos is in the mapped portion of the read,
// but the reference base was skipped, then the returned value will be
// (0, true).
func QualAtPos(record *sam.Record, pos int) (byte, bool) {
	index, overlap := indexAtPos(record, pos)
	if !overlap {
		return 0, false
	} else if index >= 0 {
		return record.Qual[index], true
	}
	return 0, true
}
