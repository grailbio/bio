package bam

import "github.com/grailbio/hts/sam"

// HasNoMappedMate returns true if record is unpaired or has an unmapped mate.
func HasNoMappedMate(record *sam.Record) bool {
	return (record.Flags&sam.Paired) == 0 || (record.Flags&sam.MateUnmapped) != 0
}
