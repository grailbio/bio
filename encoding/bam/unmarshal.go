package bam

import (
	"encoding/binary"
	"fmt"
	"reflect"
	"unsafe"

	"bytes"
	"errors"
	"github.com/biogo/hts/sam"
)

const sizeofSliceHeader = int(unsafe.Sizeof(reflect.SliceHeader{}))
const bamFixedBytes = 32

// parseAux examines the data of a SAM record's OPT fields,
// returning a slice of sam.Aux that are backed by the original data.
func parseAux(aux []byte, aa []sam.Aux) {
	naa := 0
	for i := 0; i+2 < len(aux); {
		t := aux[i+2]
		switch j := jumps[t]; {
		case j > 0:
			j += 3
			aa[naa] = sam.Aux(aux[i : i+j : i+j])
			naa++
			i += j
		case j < 0:
			switch t {
			case 'Z', 'H':
				var (
					j int
					v byte
				)
				for j, v = range aux[i:] {
					if v == 0 { // C string termination
						break // Truncate terminal zero.
					}
				}
				aa[naa] = sam.Aux(aux[i : i+j : i+j])
				naa++
				i += j + 1
			case 'B':
				length := binary.LittleEndian.Uint32(aux[i+4 : i+8])
				j = int(length)*jumps[aux[i+3]] + int(unsafe.Sizeof(length)) + 4
				aa[naa] = sam.Aux(aux[i : i+j : i+j])
				naa++
				i += j
			}
		default:
			panic(fmt.Sprintf("bam: unrecognised optional field type: %q", t))
		}
	}
}

// Round "off" up so that it is a multiple of 8. Used when storing a pointer in
// []byte.  8-byte alignment is sufficient for all CPUs we care about.
func alignOffset(off int) int {
	const pointerSize = 8
	return ((off-1)/pointerSize + 1) * pointerSize
}

var jumps = [256]int{
	'A': 1,
	'c': 1, 'C': 1,
	's': 2, 'S': 2,
	'i': 4, 'I': 4,
	'f': 4,
	'Z': -1,
	'H': -1,
	'B': -1,
}

var (
	errCorruptAuxField = errors.New("Corrupt aux field")
	errRecordTooShort  = errors.New("Record too short")
)

// countAuxFields examines the data of a SAM record's OPT field to determine
// the number of auxFields there are.
func countAuxFields(aux []byte) (int, error) {
	naux := 0
	for i := 0; i+2 < len(aux); {
		t := aux[i+2]
		switch j := jumps[t]; {
		case j > 0:
			j += 3
			i += j
			naux++
		case j < 0:
			switch t {
			case 'Z', 'H':
				var (
					j int
					v byte
				)
				for j, v = range aux[i:] {
					if v == 0 { // C string termination
						break // Truncate terminal zero.
					}
				}
				i += j + 1
				naux++
			case 'B':
				if len(aux) < i+8 {
					return -1, errCorruptAuxField
				}
				length := binary.LittleEndian.Uint32(aux[i+4 : i+8])
				j = int(length)*jumps[aux[i+3]] + int(unsafe.Sizeof(length)) + 4
				i += j
				naux++
			}
		default:
			return -1, errCorruptAuxField
		}
	}
	return naux, nil
}

// Unmarshal a serialized BAM record.
func Unmarshal(b []byte, header *sam.Header) (*Record, error) {
	if len(b) < bamFixedBytes {
		return nil, errRecordTooShort
	}
	// Need to use int(int32(uint32)) to ensure 2's complement extension of -1.
	rec := GetFromFreePool()
	refID := int(int32(binary.LittleEndian.Uint32(b)))
	rec.Pos = int(int32(binary.LittleEndian.Uint32(b[4:])))
	nLen := int(b[8])
	rec.MapQ = b[9]
	nCigar := int(binary.LittleEndian.Uint16(b[12:]))
	rec.Flags = sam.Flags(binary.LittleEndian.Uint16(b[14:]))
	lSeq := int(binary.LittleEndian.Uint32(b[16:]))
	nextRefID := int(int32(binary.LittleEndian.Uint32(b[20:])))
	rec.MatePos = int(int32(binary.LittleEndian.Uint32(b[24:])))
	rec.TempLen = int(int32(binary.LittleEndian.Uint32(b[28:])))

	// Read variable length data.
	srcVariableBytes := len(b) - bamFixedBytes

	nDoubletBytes := (lSeq + 1) >> 1
	srcAuxOffset := bamFixedBytes + nLen + (nCigar * 4) + nDoubletBytes + lSeq
	if len(b) < srcAuxOffset {
		return nil, fmt.Errorf("Corrupt BAM aux record: len(b)=%d, auxoffset=%d", len(b), srcAuxOffset)
	}
	nAuxFields, err := countAuxFields(b[srcAuxOffset:])
	if err != nil {
		return nil, err
	}

	shadowCigarOffset := alignOffset(srcVariableBytes)               // store the cigar int32s here
	shadowAuxOffset := alignOffset(shadowCigarOffset + (nCigar * 4)) // store the AuxFields here
	shadowSize := shadowAuxOffset + (nAuxFields * sizeofSliceHeader)

	// shadowBuf is used as an 'arena' from which all objects/slices
	// required to store the result of parsing the bam alignment record.
	// This reduces the load on GC and consequently allows for better
	// scalability with the number of cores used by clients of this package.
	ResizeScratch(&rec.Scratch, shadowSize)
	shadowBuf := rec.Scratch
	copy(shadowBuf, b[bamFixedBytes:])

	bufHdr := (*reflect.SliceHeader)(unsafe.Pointer(&shadowBuf))

	// Note that rec.Name now points to the shadow buffer
	hdr := (*reflect.StringHeader)(unsafe.Pointer(&rec.Name))
	hdr.Data = uintptr(unsafe.Pointer(bufHdr.Data))
	hdr.Len = nLen - 1 // drop trailing '\0'
	shadowOffset := nLen

	var sliceHdr *reflect.SliceHeader
	if nCigar > 0 {
		for i := 0; i < nCigar; i++ {
			*(*uint32)(unsafe.Pointer(&shadowBuf[shadowCigarOffset+(i*4)])) = uint32(binary.LittleEndian.Uint32(shadowBuf[shadowOffset+(i*4):]))
		}
		rec.Cigar = UnsafeBytesToCigar(shadowBuf[shadowCigarOffset : shadowCigarOffset+nCigar*4])
		shadowOffset += nCigar * 4
	} else {
		rec.Cigar = nil
	}

	rec.Seq.Length = int(lSeq)
	rec.Seq.Seq = UnsafeBytesToDoublets(shadowBuf[shadowOffset : shadowOffset+nDoubletBytes])
	shadowOffset += nDoubletBytes

	rec.Qual = shadowBuf[shadowOffset : shadowOffset+lSeq]
	shadowOffset += lSeq
	if nAuxFields > 0 {
		// Clear the array before updating rec.AuxFields. GC will be
		// confused otherwise.
		for i := shadowAuxOffset; i < shadowAuxOffset+nAuxFields*sizeofSliceHeader; i++ {
			shadowBuf[i] = 0
		}
		sliceHdr = (*reflect.SliceHeader)(unsafe.Pointer(&rec.AuxFields))
		sliceHdr.Data = uintptr(unsafe.Pointer(bufHdr.Data + uintptr(shadowAuxOffset)))
		sliceHdr.Len = nAuxFields
		sliceHdr.Cap = sliceHdr.Len
		parseAux(shadowBuf[shadowOffset:srcVariableBytes], rec.AuxFields)
	}

	refs := len(header.Refs())
	if refID != -1 {
		if refID < -1 || refID >= refs {
			return nil, fmt.Errorf("bam: reference id %v out of range", refID)
		}
		rec.Ref = header.Refs()[refID]
	}
	if nextRefID != -1 {
		if refID == nextRefID {
			rec.MateRef = rec.Ref
			return rec, nil
		}
		if nextRefID < -1 || nextRefID >= refs {
			return nil, fmt.Errorf("bam: mate reference id %v out of range", nextRefID)
		}
		rec.MateRef = header.Refs()[nextRefID]
	}
	return rec, nil
}

// UnmarshalHeader parses a sam.Header encoded in BAM binary format.
func UnmarshalHeader(buf []byte) (*sam.Header, error) {
	header, err := sam.NewHeader(nil, nil)
	if err != nil {
		return nil, err
	}
	hr := bytes.NewReader(buf)
	if err := header.DecodeBinary(hr); err != nil {
		return nil, err
	}
	if hr.Len() > 0 {
		return nil, fmt.Errorf("%d byte junk at the end of SAM header", hr.Len())
	}
	return header, nil
}
