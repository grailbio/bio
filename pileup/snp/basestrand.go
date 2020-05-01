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
package snp

import (
	"bytes"
	"context"
	"encoding/binary"
	"fmt"
	"io"
	"strconv"
	"strings"
	"sync"

	"github.com/grailbio/base/file"
	"github.com/grailbio/base/log"
	"github.com/grailbio/base/recordio"
	"github.com/grailbio/base/recordio/recordiozstd"
	"github.com/grailbio/base/tsv"
	"github.com/grailbio/bio/pileup"
)

const (
	refNamesHeader = "RefNames"
	trailerVersion = 1
)

// BaseStrandPile represents a single pileup entry with a count for every
// (base, strand) tuple.
//
// - Pos is zero-based; it is necessary to add 1 when converting to most text
//   formats (but not BED).
// - In Counts[][], base is the major dimension, with pileup.BaseA=0, C=1, G=2,
//   T=3.  Strand is the minor dimension, with strandFwd=0 and strandRev=1.
//   TODO(cchang): strandFwd=0, strandRev=1 is inconsistent with bio/pileup's
//   internal representation (which has None=0).  We have enough other code at
//   this point with Fwd=0, Rev=1 that it's probably time to change
//   bio/pileup's representation to match that.
type BaseStrandPile struct {
	RefID  uint32
	Pos    uint32
	Counts [pileup.NBase][2]uint32
}

// Minimal implementation; can optimize later if we're using this on large
// pileups.

// WriteBaseStrandsRio writes the given BaseStrand-pileup entries to the given
// writer, using recordio.
func WriteBaseStrandsRio(piles []BaseStrandPile, refNames []string, out io.Writer) error {
	// recordiozstd.Init() is called in singleton.go's init().
	recordWriter := recordio.NewWriter(out, recordio.WriterOpts{
		Marshal:      marshalBaseStrand,
		Transformers: []string{recordiozstd.Name},
	})
	// could error out if refNames is empty
	recordWriter.AddHeader(refNamesHeader, strings.Join(refNames, "\000"))
	recordWriter.AddHeader(recordio.KeyTrailer, true)
	for i := range piles {
		recordWriter.Append(&piles[i])
	}
	recordWriter.SetTrailer(baseStrandsRioTrailer(len(piles)))
	return recordWriter.Finish()
}

func baseStrandsRioTrailer(numPiles int) []byte {
	var buffer bytes.Buffer
	if err := binary.Write(&buffer, binary.LittleEndian, int64(trailerVersion)); err != nil {
		panic("couldn't write trailer version")
	}
	if err := binary.Write(&buffer, binary.LittleEndian, int64(numPiles)); err != nil {
		panic("couldn't write numPiles to trailer")
	}
	return buffer.Bytes()
}

func parseBaseStrandsTrailer(trailer []byte) (int64, error) {
	r := bytes.NewReader(trailer)
	var version, numPiles int64
	if err := binary.Read(r, binary.LittleEndian, &version); err != nil {
		return 0, err
	}
	if version != trailerVersion {
		return 0, fmt.Errorf("unrecognized trailer version: got %d, want %d", version, trailerVersion)
	}
	if err := binary.Read(r, binary.LittleEndian, &numPiles); err != nil {
		return 0, err
	}
	return numPiles, nil
}

func marshalBaseStrand(scratch []byte, p interface{}) ([]byte, error) {
	t := scratch
	if len(t) < 40 {
		t = make([]byte, 40)
	}

	// As of Go 1.12, the bounds-check-eliminator is not smart enough to
	// determine that len(t) is always 40 if we move t = scratch[:40] in an
	// else-branch above.
	t = t[:40]

	pile := p.(*BaseStrandPile)
	binary.LittleEndian.PutUint32(t[:4], pile.RefID)
	binary.LittleEndian.PutUint32(t[4:8], pile.Pos)
	// I tried a few different ways to write a double-loop here, but the
	// bounds-check-eliminator did not handle any of them well (on every
	// inner-loop iteration, there were three separate t[] bounds-checks:
	// offset < 0?  offset > offset+4?  offset+4 > len(t)?)
	//
	// The code below is handled decently.
	//
	// The trailing 0/1 indexes will be replaced with enums when bio/pileup as a
	// whole has been switched over to strandFwd=0, strandRev=1, strandNone=2.
	binary.LittleEndian.PutUint32(t[8:12], pile.Counts[pileup.BaseA][0])
	binary.LittleEndian.PutUint32(t[12:16], pile.Counts[pileup.BaseA][1])
	binary.LittleEndian.PutUint32(t[16:20], pile.Counts[pileup.BaseC][0])
	binary.LittleEndian.PutUint32(t[20:24], pile.Counts[pileup.BaseC][1])
	binary.LittleEndian.PutUint32(t[24:28], pile.Counts[pileup.BaseG][0])
	binary.LittleEndian.PutUint32(t[28:32], pile.Counts[pileup.BaseG][1])
	binary.LittleEndian.PutUint32(t[32:36], pile.Counts[pileup.BaseT][0])
	binary.LittleEndian.PutUint32(t[36:40], pile.Counts[pileup.BaseT][1])
	return t, nil
}

// ReadBaseStrandsRio reads BaseStrand piles from a recordio file written by
// WriteBaseStrandsRio.
func ReadBaseStrandsRio(rs io.ReadSeeker) (piles []BaseStrandPile, refNames []string, err error) {
	var unmarshaller BaseStrandUnmarshaller
	scanner := recordio.NewScanner(rs, recordio.ScannerOpts{
		Unmarshal: unmarshaller.UnmarshalBaseStrand,
	})
	if len(scanner.Trailer()) != 0 {
		var numPiles int64
		if numPiles, err = parseBaseStrandsTrailer(scanner.Trailer()); err != nil {
			return
		}
		unmarshaller.init(numPiles)
	}

	hdr := scanner.Header()
	for _, kv := range hdr {
		switch kv.Key {
		case refNamesHeader:
			packedRefNames := kv.Value.(string)
			refNames = strings.Split(packedRefNames, "\000")
			// Cannot return an error on unrecognized key since recordio can write its own.
		}
	}

	for scanner.Scan() {
		piles = append(piles, *scanner.Get().(*BaseStrandPile))
	}
	err = scanner.Err()
	return
}

// BaseStrandUnmarshaller is used to allocate memory in large blocks during
// unmarshalling, to prevent contention with other goroutines.
type BaseStrandUnmarshaller struct {
	piles  []BaseStrandPile
	offset int
}

func (b *BaseStrandUnmarshaller) init(size int64) {
	if b.piles != nil {
		panic("tried to initialize when already initialized")
	}
	b.piles = make([]BaseStrandPile, size)
}

func (b *BaseStrandUnmarshaller) UnmarshalBaseStrand(in []byte) (out interface{}, err error) {
	in = in[:40] // help the bounds-checker
	if b.offset == len(b.piles) {
		b.piles = append(b.piles, BaseStrandPile{})
	}
	pile := &b.piles[b.offset]
	b.offset++
	pile.RefID = binary.LittleEndian.Uint32(in[:4])
	pile.Pos = binary.LittleEndian.Uint32(in[4:8])
	pile.Counts[pileup.BaseA][0] = binary.LittleEndian.Uint32(in[8:12])
	pile.Counts[pileup.BaseA][1] = binary.LittleEndian.Uint32(in[12:16])
	pile.Counts[pileup.BaseC][0] = binary.LittleEndian.Uint32(in[16:20])
	pile.Counts[pileup.BaseC][1] = binary.LittleEndian.Uint32(in[20:24])
	pile.Counts[pileup.BaseG][0] = binary.LittleEndian.Uint32(in[24:28])
	pile.Counts[pileup.BaseG][1] = binary.LittleEndian.Uint32(in[28:32])
	pile.Counts[pileup.BaseT][0] = binary.LittleEndian.Uint32(in[32:36])
	pile.Counts[pileup.BaseT][1] = binary.LittleEndian.Uint32(in[36:40])
	return pile, nil
}

// WriteBaseStrandToTSV writes a []BaseStrandPile as a TSV.
func WriteBaseStrandToTSV(piles []BaseStrandPile, refNames []string, w io.Writer) (err error) {
	// Note that this is slightly different from the .basestrand.tsv format:
	// there is no REF base column, or '#' marking the header line.
	outTSV := tsv.NewWriter(w)
	outTSV.WriteString("CHROM\tPOS\tA+\tA-\tC+\tC-\tG+\tG-\tT+\tT-")
	if err = outTSV.EndLine(); err != nil {
		return
	}
	for _, pile := range piles {
		outTSV.WriteString(refNames[pile.RefID])
		outTSV.WriteUint32(pile.Pos)
		outTSV.WriteUint32(pile.Counts[pileup.BaseA][0])
		outTSV.WriteUint32(pile.Counts[pileup.BaseA][1])
		outTSV.WriteUint32(pile.Counts[pileup.BaseC][0])
		outTSV.WriteUint32(pile.Counts[pileup.BaseC][1])
		outTSV.WriteUint32(pile.Counts[pileup.BaseG][0])
		outTSV.WriteUint32(pile.Counts[pileup.BaseG][1])
		outTSV.WriteUint32(pile.Counts[pileup.BaseT][0])
		outTSV.WriteUint32(pile.Counts[pileup.BaseT][1])
		if err = outTSV.EndLine(); err != nil {
			return
		}
	}
	if err = outTSV.Flush(); err != nil {
		return
	}
	return
}

// WriteBaseStrandsRioAsTSV converts the given recordio pileup to TSV.
func WriteBaseStrandsRioAsTSV(ctx context.Context, path string, w io.Writer) error {
	file, err := file.Open(ctx, path)
	if err != nil {
		return err
	}
	piles, refNames, err2 := ReadBaseStrandsRio(file.Reader(ctx))
	if err2 != nil {
		return err2
	}
	return WriteBaseStrandToTSV(piles, refNames, w)
}

// BaseStrandTsvRow represents a single row of a basestrand.tsv file.
type BaseStrandTsvRow struct {
	Chr  string `tsv:"#CHROM"` // Chromosome
	Pos  int64  `tsv:"POS"`    // Position in chromosome
	Ref  string `tsv:"REF"`    // Reference base
	FwdA int64  `tsv:"A+"`     // A count on the forward strand
	RevA int64  `tsv:"A-"`     // A count on the reverse strand
	FwdC int64  `tsv:"C+"`     // C count on the forward strand
	RevC int64  `tsv:"C-"`     // C count on the reverse strand
	FwdG int64  `tsv:"G+"`     // G count on the forward strand
	RevG int64  `tsv:"G-"`     // G count on the reverse strand
	FwdT int64  `tsv:"T+"`     // T count on the forward strand
	RevT int64  `tsv:"T-"`     // T count on the reverse strand
}

// ReadBaseStrandTsv reads a basestrand.tsv file from the given io.Reader.
func ReadBaseStrandTsv(r io.Reader) ([]BaseStrandTsvRow, error) {
	tsvReader := tsv.NewReader(r)
	tsvReader.Comment = '#'

	rows := make([]BaseStrandTsvRow, 0)
	for {
		var row BaseStrandTsvRow
		if err := tsvReader.Read(&row); err != nil {
			if err == io.EOF {
				break
			} else {
				return nil, err
			}
		}
		rows = append(rows, row)
	}
	return rows, nil
}

// WriteBaseStrandTsv writes a basestrand.tsv file to the given writer
func WriteBaseStrandTsv(rows []BaseStrandTsvRow, writer io.Writer) error {
	tsvWriter := tsv.NewRowWriter(writer)
	for _, row := range rows {
		if err := tsvWriter.Write(&row); err != nil {
			return err
		}
	}
	return tsvWriter.Flush()
}

// singleStrandBaseStrandTsvRow represents a single row of a strand.<fwd/rev>.tsv file.
type singleStrandBaseStrandTsvRow struct {
	Chr   string `tsv:"CHROM"` // Chromosome
	Pos   int64  `tsv:"POS"`   // Position in chromosome
	Depth int64  `tsv:"DEPTH"` // Depth
	Ref   string `tsv:"REF"`   // Reference base
	A     int64  `tsv:"A"`     // A count
	C     int64  `tsv:"C"`     // C count
	G     int64  `tsv:"G"`     // G count
	T     int64  `tsv:"T"`     // T count
	N     int64  `tsv:"N"`     // N count
	Ins   int64  `tsv:"INS"`   // Insertion
	Del   int64  `tsv:"DEL"`   // Deletion
}

// ReadSingleStrandBaseStrandTsv reads strand specific strand.<fwd/rev>.snp.tsv
// files from the given io.Reader.
func ReadSingleStrandBaseStrandTsv(forward, reverse io.Reader) ([]BaseStrandTsvRow, error) {
	fwdReader := tsv.NewReader(forward)
	fwdReader.HasHeaderRow = true
	fwdReader.UseHeaderNames = true

	revReader := tsv.NewReader(reverse)
	revReader.HasHeaderRow = true
	revReader.UseHeaderNames = true

	rows := make([]BaseStrandTsvRow, 0)
	for {
		var fwdRow, revRow singleStrandBaseStrandTsvRow

		// Read in the forward strand.
		if err := fwdReader.Read(&fwdRow); err != nil {
			if err == io.EOF {
				break
			} else {
				return nil, err
			}
		}

		// Read in the reverse strand.
		if err := revReader.Read(&revRow); err != nil {
			if err == io.EOF {
				break
			} else {
				return nil, err
			}
		}

		// Intersect the rows.
		// Check that the chromosome position, and reference are equal.
		if fwdRow.Chr != revRow.Chr || fwdRow.Pos != revRow.Pos {
			return nil, fmt.Errorf("expected equal chromosome position but got fwd: %s %d and rev: %s %d",
				fwdRow.Chr, fwdRow.Pos, revRow.Chr, revRow.Pos)
		}
		if fwdRow.Ref != revRow.Ref {
			return nil, fmt.Errorf("expected equal reference alleles but got fwd: %s and rev: %s", fwdRow.Ref, revRow.Ref)
		}

		row := BaseStrandTsvRow{
			Chr:  fwdRow.Chr,
			Pos:  fwdRow.Pos,
			Ref:  fwdRow.Ref,
			FwdA: fwdRow.A,
			RevA: revRow.A,
			FwdC: fwdRow.C,
			RevC: revRow.C,
			FwdG: fwdRow.G,
			RevG: revRow.G,
			FwdT: fwdRow.T,
			RevT: revRow.T,
		}
		rows = append(rows, row)
	}

	return rows, nil
}

// ReadBaseStrandTsvIntoChannel reads a basestrand.tsv file from the given tsv.Reader into the given channel.
func ReadBaseStrandTsvIntoChannel(reader *tsv.Reader, c chan []BaseStrandTsvRow, bufferLen int, fileName string, wg *sync.WaitGroup) {
	wg.Add(1)
	index := 0
	buffer := make([]BaseStrandTsvRow, 0, bufferLen)
	var row BaseStrandTsvRow
	for {
		if err := reader.Read(&row); err != nil {
			if err == io.EOF {
				c <- buffer
				break
			} else {
				log.Panicf("could not parse BaseStrandTsvRow file error %v", err.Error())
			}
		}
		if len(buffer) == cap(buffer) {
			c <- buffer
			buffer = make([]BaseStrandTsvRow, 0, bufferLen)
		}
		buffer = append(buffer, row)
		index++
		if index%1000000 == 0 {
			log.Info.Printf("Processed %v Lines in %v", index, fileName)
		}
	}
	close(c)
	wg.Done()
	log.Info.Printf("Finished reading %v", fileName)
}

// ChrId returns the index for the given chromosome string in bio/pileup
// format.
func ChrId(chr string) (int, error) {
	// In our current set of reference files, the chromosome ID is prefixed with
	// "chr".
	chrId := chr[3:]
	switch chrId {
	case "X":
		return 23, nil
	case "Y":
		return 24, nil
	case "M":
		return 25, nil
	default:
		i, err := strconv.Atoi(chrId)
		if err != nil {
			return 0, err
		}
		if i <= 0 || i >= 23 {
			return 0, fmt.Errorf("%d is not a valid chromosome id", i)
		}
		return i, nil
	}
}
