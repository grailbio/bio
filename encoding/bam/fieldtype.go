package bam

import (
	"fmt"
)

// FieldType defines a sam.Record field. Each field is stored in a separate
// file.
type FieldType uint8

const (
	// FieldCoord combines <sam.Reference.ID(), sam.Record.Pos>. They need to be
	// read together when seeking to a specific coordinate.
	FieldCoord FieldType = iota
	// Rest of Field* stands for the sam.Record field with the same name.
	FieldFlags
	FieldMapq
	FieldCigar
	FieldMateRefID
	FieldMatePos
	FieldTempLen
	FieldName
	FieldSeq
	FieldQual
	FieldAux

	// FieldInvalid is a sentinel
	FieldInvalid
	MinField   = FieldCoord
	NumFields  = int(FieldInvalid)
	MaxNumCols = NumFields
)

var fieldNames = []string{
	"coord",
	"flags",
	"mapq",
	"cigar",
	"materefid",
	"matepos",
	"templen",
	"name",
	"seq",
	"qual",
	"aux",
}

func (f FieldType) String() string {
	if f >= 0 && int(f) < len(fieldNames) {
		return fieldNames[f]
	}
	return fmt.Sprintf("Field%d", f)
}

// ParseFieldType converts a string to FieldType. For example, "Cigar" will
// return FieldCigar.
func ParseFieldType(v string) (FieldType, error) {
	for f, name := range fieldNames {
		if name == v {
			return FieldType(f), nil
		}
	}
	return FieldAux, fmt.Errorf("%v: invalid PAM field type", v)
}
