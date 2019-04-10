package bam

import "github.com/grailbio/hts/sam"

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
