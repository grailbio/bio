// Copyright 2018 GRAIL, Inc. All rights reserved.
// Use of this source code is governed by the Apache 2.0
// license that can be found in the LICENSE file.

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

// FieldNames lists all the bam Field names.
var FieldNames = []string{
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

// String returns the name of the type.  The name is used as part of the PAM
// filenames, so it shall not be changed.
func (f FieldType) String() string {
	if int(f) < len(FieldNames) {
		return FieldNames[f]
	}
	return fmt.Sprintf("Field%d", f)
}

// ParseFieldType converts a string to FieldType. For example, "cigar" will
// return FieldCigar.
func ParseFieldType(v string) (FieldType, error) {
	for f, name := range FieldNames {
		if name == v {
			return FieldType(f), nil
		}
	}
	return FieldAux, fmt.Errorf("%v: invalid PAM field type", v)
}
