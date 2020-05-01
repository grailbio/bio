// Copyright 2020 Grail Inc.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
package pileup

import (
	"bufio"
	"context"
	"fmt"
	"io"
	"strings"

	"github.com/grailbio/base/file"
	"github.com/grailbio/base/fileio"
	gunsafe "github.com/grailbio/base/unsafe"
	"github.com/grailbio/bio/biosimd"
	"github.com/grailbio/bio/interval"
	"github.com/grailbio/hts/sam"
	"github.com/klauspost/compress/gzip"
)

// Common pileup components.

// PosType is the integer type used to represent genomic positions.
type PosType = interval.PosType

// PosTypeMax is the maximum value that can be represented by a PosType.
const PosTypeMax = interval.PosTypeMax

// These constants have two relevant meanings:
// 1. In the .bam seq[] encoding (sam.BaseA, sam.BaseC, etc.), it's the
//    position of A's set bit.  This is relevant when using
//    __builtin_ctzl([read seq word] & (~[ref seq word])) to quickly iterate
//    over differences-from-reference.
// 2. It's the natural value for A/C/G/T in a packed 2-bit representation
//    (useful anywhere we don't have to worry about Ns).

const (
	// BaseA represents an A base.
	BaseA byte = iota
	// BaseC represents an C base.
	BaseC
	// BaseG represents an G base.
	BaseG
	// BaseT represents an T base.
	BaseT
	// BaseX is a catch-all.
	BaseX
)

const (
	// NBase is the number of regular base types.
	NBase = 4
	// NBaseEnum counts BaseX as well as the regular base types.
	NBaseEnum = 5
)

// Seq8ToEnumTable is the .bam seq nibble -> A/C/G/T/X enum mapping.
var Seq8ToEnumTable = [...]byte{BaseX, BaseA, BaseC, BaseX, BaseG, BaseX, BaseX, BaseX, BaseT, BaseX, BaseX, BaseX, BaseX, BaseX, BaseX, BaseX}

// EnumToASCIITable is the A/C/G/T/X -> ASCII mapping, with X rendered as 'N'.
var EnumToASCIITable = [...]byte{'A', 'C', 'G', 'T', 'N'}

// Seq8ToASCIITable is the .bam seq nibble -> ASCII mapping.
var Seq8ToASCIITable = [...]byte{'=', 'A', 'C', 'M', 'G', 'R', 'S', 'V', 'T', 'W', 'Y', 'H', 'K', 'D', 'B', 'N'}

// StrandType describes which strand a read-pair is aligned to.
type StrandType int

const (
	// StrandNone means either no strand restriction, or (when returned by
	// GetStrand) undefined-strand (read ends on different chromosomes, or appear
	// to be part of an inversion).
	StrandNone StrandType = iota
	// StrandFwd means that the read-pair's start is on the 5' side and the end
	// is on the 3' side of the same chromosome.
	StrandFwd
	// StrandRev means that the read-pair's start is on the 3' side and the end
	// is on the 5' side of the same chromosome.
	StrandRev
)

// StrandTypeToASCIITable is the StrandType -> ASCII mapping.
var StrandTypeToASCIITable = [...]byte{'.', '+', '-'}

// GetStrand returns the strand the read-pair is aligned to.
func GetStrand(samr *sam.Record) StrandType {
	if samr.Ref != samr.MateRef {
		return StrandNone
	}
	flagStrand := samr.Flags & (sam.Reverse | sam.MateReverse | sam.Read1 | sam.Read2)
	if (flagStrand == (sam.MateReverse | sam.Read1)) || (flagStrand == (sam.Reverse | sam.Read2)) {
		return StrandFwd
	} else if (flagStrand == (sam.Reverse | sam.Read1)) || (flagStrand == (sam.MateReverse | sam.Read2)) {
		return StrandRev
	}
	if samr.Flags&sam.MateUnmapped == sam.MateUnmapped {
		// Support an alternate encoding emitted by some 'collapser' programs.
		flagStrand &= sam.Reverse | sam.MateReverse
		if flagStrand == 0 {
			return StrandFwd
		} else if flagStrand == (sam.Reverse | sam.MateReverse) {
			return StrandRev
		}
	}
	return StrandNone
}

// ParseCols parses a column-set-descriptor string given on the command line
// (colsParam) into a 64-bit integer bitset for internal use.
func ParseCols(colsParam string, colNameMap map[string]int, defaultColBitset int) (colBitset int, err error) {
	if colsParam == "" {
		return defaultColBitset, nil
	}

	colsParamParts := strings.Split(colsParam, ",")
	// Two cases:
	// 1. Each part has a '+' or a '-' in front.  Treat these as patches to the
	//    default column set.
	// 2. No part has a '+' or a '-' in front.  Ignore the default and treat this
	//    as the full set.
	firstChar := colsParamParts[0][0]
	if (firstChar == '+') || (firstChar == '-') {
		colBitset = defaultColBitset
		for _, part := range colsParamParts {
			firstChar = part[0]
			if (firstChar != '+') && (firstChar != '-') {
				err = fmt.Errorf("parseCols: either all terms in column set descriptor must be preceded by +/-, or none can be")
				return
			}
			v := colNameMap[part[1:]]
			if v == 0 {
				err = fmt.Errorf("parseCols: %v not found", part[1:])
				return
			}
			// can also check for duplicates
			if firstChar == '+' {
				colBitset |= v
			} else {
				colBitset &= ^v
			}
		}
	} else {
		for _, part := range colsParamParts {
			firstChar = part[0]
			if (firstChar == '+') || (firstChar == '-') {
				err = fmt.Errorf("parseCols: either all terms in column set descriptor must be preceded by +/-, or none can be")
				return
			}
			v := colNameMap[part]
			if v == 0 {
				err = fmt.Errorf("parseCols: %v not found", part)
				return
			}
			colBitset |= v
		}
	}
	return colBitset, nil
}

// LoadFa loads a .fa file, and currently translates the contents to
// A=1/C=2/G=4/T=8/other=15 encoding ("seq8") to facilitate efficient
// comparisons against BAM/PAM Seq data.
// An option may be added in the future for conversion to
// A=0/C=1/G=2/T=3/other=4 enum encoding; that may be better for some purposes.
// TODO(cchang): Use grailbio/bio/encoding/fasta instead, updating that package
// as needed.
func LoadFa(ctx context.Context, fapath string, maxline int, headerRefs []*sam.Reference) (refSeqs [][]byte, err error) {
	var infile file.File
	if infile, err = file.Open(ctx, fapath); err != nil {
		return
	}
	defer func() {
		if e := infile.Close(ctx); e != nil && err == nil {
			err = e
		}
	}()
	reader := io.Reader(infile.Reader(ctx))
	switch fileio.DetermineType(fapath) {
	case fileio.Gzip:
		if reader, err = gzip.NewReader(reader); err != nil {
			return
		}
	}
	scanner := bufio.NewScanner(reader)
	// Provide a way to accept .fa files which put an entire reference on a
	// single line.
	startSize := bufio.MaxScanTokenSize
	if startSize > maxline {
		startSize = maxline
	}
	buf := make([]byte, startSize, maxline)
	scanner.Buffer(buf, maxline)

	bamRefMap := make(map[string]int)
	for i, curRef := range headerRefs {
		name := curRef.Name()
		bamRefMap[name] = i
		// possible todo: tolerate '1' vs. 'chr1', etc.
	}
	refSeqs = make([][]byte, len(headerRefs))

	lineIdx := 0
	refIdx := 0
	keepRef := false
	refSeq := []byte{}
	refPos := 0
	for scanner.Scan() {
		lineIdx++
		curLine := scanner.Bytes()
		nByte := len(curLine)
		if nByte == 0 {
			continue
		}
		if curLine[0] == '>' {
			if keepRef {
				if refPos != len(refSeq) {
					err = fmt.Errorf("loadFa: inconsistent lengths for contig %s (%d in .bam header, %d in .fa)", headerRefs[refIdx].Name(), len(refSeq), refPos)
					return
				}
				refSeqs[refIdx] = refSeq
			}
			refIdx, keepRef = bamRefMap[gunsafe.BytesToString(curLine[1:])]
			if keepRef {
				newRef := headerRefs[refIdx]
				refSeq = make([]byte, newRef.Len())
				refPos = 0
			}
			continue
		}
		if !keepRef {
			continue
		}
		biosimd.ASCIIToSeq8(refSeq[refPos:refPos+nByte], curLine)
		refPos += nByte
	}
	if err = scanner.Err(); err != nil {
		return
	}
	if keepRef {
		if refPos != len(refSeq) {
			err = fmt.Errorf("loadFa: inconsistent lengths for ref %s (%d in .bam header, %d in .fa)", headerRefs[refIdx].Name(), len(refSeq), refPos)
			return
		}
		refSeqs[refIdx] = refSeq
	}
	return
}
