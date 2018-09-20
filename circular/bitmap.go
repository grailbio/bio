package circular

import (
	"github.com/grailbio/base/bitset"
	"github.com/grailbio/base/log"
	"github.com/grailbio/base/simd"
	bi "github.com/grailbio/bio/interval"
)

// BitsPerWord is the number of bits per machine word.  (Don't want to import
// base/simd or base/bitset in files where we only need this constant.)
const BitsPerWord = simd.BitsPerWord

// indexNonzeroRev returns the index of the first encountered nonzero byte in
// arr[], treating it as a circular buffer, and searching backward from start.
// If no nonzero byte is encountered before stop, stop is returned.
// May export this later, but let's wait until there's a known external
// consumer.
func indexNonzeroRev(arr []byte, start, stop int) int {
	nCirc := len(arr)
	mask := nCirc - 1
	offset := start &^ mask
	circPos := start & mask
	circStop := stop & mask
	// Probable todo: add simd.LastGreater8Unsafe(), since it's a straightforward
	// variation of simd.FirstGreater8Unsafe().
	circPause := circStop
	if circStop > circPos {
		circPause = 0
	}
	for ; ; circPos-- {
		if arr[circPos] != 0 {
			return circPos + offset
		}
		if circPos == circPause {
			if circPos == circStop {
				return stop
			}
			circPos = nCirc
			offset -= nCirc
		}
	}
}

// FirstPosEmpty is an empty-bitmap sentinel value.  It must be larger than any
// real coordinate.
const FirstPosEmpty = bi.PosTypeMax

// Bitmap is a 2-dimensional bitmap, with circular major dimension.
// Variable and type names currently make the assumption that the major
// dimension corresponds to position; this can be abstracted out.
type Bitmap struct {
	// bits stores the raw bits.  Logical row n of the bitmap is
	// bits[n*rowWidth:(n+1)*rowWidth].
	bits []uintptr
	// wordPops[i] stores the number of nonzero words in bits[i].  This lets us
	// find the next nonempty bitarray with a (circular) find-next-nonzero-byte
	// operation.
	wordPops []byte
	// firstPos stores the position of the first entry in the table (high bits
	// preserved), or FirstPosEmpty when the table is empty.
	// (We may want the PosType definition to be owned by another package.)
	firstPos bi.PosType
	// lastPos stores the position of the last entry in the table (high bits
	// preserved), or -1 when the table is empty.
	lastPos bi.PosType
	// rowWidth stores the number of words in each logical bitmap row.  Width is
	// currently limited to 255.
	rowWidth bi.PosType
}

// NewBitmap creates an empty Bitmap.
func NewBitmap(nCirc, rowWidth bi.PosType) (b Bitmap) {
	if rowWidth > 255 {
		// wordPops would need to be []uint16 instead, or we'd need to switch to a
		// different data structure.
		log.Panicf("circular.Bitmap does not support rowWidth > 255")
	}
	if (nCirc & (nCirc - 1)) != 0 {
		log.Panicf("circular.Bitmap requires nCirc to be a power of two")
	}
	b.bits = make([]uintptr, nCirc*rowWidth)
	capacity := nCirc
	if capacity < bi.PosType(simd.BytesPerVec()) {
		capacity = bi.PosType(simd.BytesPerVec())
	}
	b.wordPops = make([]byte, nCirc, capacity)
	b.firstPos = FirstPosEmpty
	b.lastPos = -1
	b.rowWidth = rowWidth
	return
}

// NCirc returns the major dimension size.
func (b *Bitmap) NCirc() bi.PosType {
	return bi.PosType(len(b.wordPops))
}

// FirstPos returns the position of the first table entry, or FirstPosEmpty
// when the table is empty.
func (b *Bitmap) FirstPos() bi.PosType {
	return b.firstPos
}

// Row returns a []uintptr corresponding to a single row of the bitmap.
func (b *Bitmap) Row(circPos bi.PosType) []uintptr {
	base := circPos * b.rowWidth
	return b.bits[base : base+b.rowWidth]
}

// Set sets a single bit of the bitmap.  (Nothing bad happens if the bit was
// already set.)
func (b *Bitmap) Set(pos, circPos bi.PosType, colIdx uint32) {
	row := b.Row(circPos)
	colWordIdx := colIdx / BitsPerWord
	curWord := row[colWordIdx]
	if curWord == 0 {
		b.wordPops[circPos]++
	}
	row[colWordIdx] = curWord | (uintptr(1) << (colIdx % BitsPerWord))
	if b.firstPos > pos {
		b.firstPos = pos
	}
	if b.lastPos < pos {
		b.lastPos = pos
	}
}

// firstNonemptyPos returns the position of the first nonempty bitmap row in
// [pos, stopPos), or stopPos if there aren't any.
func (b *Bitmap) firstNonemptyPos(pos, stopPos bi.PosType) bi.PosType {
	arr := b.wordPops
	nCirc := len(arr)
	mask := nCirc - 1
	offset := int(pos) &^ mask
	circPos := int(pos) & mask
	circStop := int(stopPos) & mask
	if circStop < circPos {
		result := simd.FirstGreater8Unsafe(arr, 0, circPos)
		if result != nCirc {
			return bi.PosType(result + offset)
		}
		circPos = 0
		offset += nCirc
	}
	return bi.PosType(simd.FirstGreater8Unsafe(arr[:circStop], 0, circPos) + offset)
}

// Clear clears a single bit of the bitmap.  (Nothing bad happens if the bit
// was already clear.)
func (b *Bitmap) Clear(pos, circPos bi.PosType, colIdx uint32) {
	row := b.Row(circPos)
	colWordIdx := colIdx / BitsPerWord
	curWord := row[colWordIdx] &^ (uintptr(1) << (colIdx % BitsPerWord))
	row[colWordIdx] = curWord
	if curWord == 0 {
		wPop := b.wordPops[circPos] - 1
		b.wordPops[circPos] = wPop
		if wPop == 0 {
			if pos == b.firstPos {
				if pos == b.lastPos {
					b.firstPos = FirstPosEmpty
					b.lastPos = -1
				} else {
					b.firstPos = b.firstNonemptyPos(pos+1, b.lastPos)
				}
			} else if pos == b.lastPos {
				b.lastPos = bi.PosType(indexNonzeroRev(b.wordPops, int(pos-1), int(b.firstPos)))
			}
		}
	}
}

// NewRowScanner returns a bitset.NonzeroWordScanner for the first nonempty
// table row, along with the position of the first set bit in the row.  It is
// assumed that the row will be fully scanned before any other changes are made
// to the table.
func (b *Bitmap) NewRowScanner() (bitset.NonzeroWordScanner, int) {
	pos := b.firstPos
	if pos == FirstPosEmpty {
		log.Panicf("Bitmap.NewRowScanner() called on an empty table.")
	}
	mask := b.NCirc() - 1
	circPos := pos & mask
	if pos < b.lastPos {
		b.firstPos = b.firstNonemptyPos(pos+1, b.lastPos)
	} else {
		b.firstPos = FirstPosEmpty
		b.lastPos = -1
	}
	nzwPop := int(b.wordPops[circPos])
	b.wordPops[circPos] = 0
	return bitset.NewNonzeroWordScanner(b.Row(circPos), nzwPop)
}

// CheckPanic verifies the following invariants for a Bitmap, panicking on
// failure:
// * If either firstPos or lastPos is set to the "empty" value
//   (FirstPosEmpty and -1, respectively), the other is as well.
// * If the bitmap isn't empty:
//   * firstPos is in [lastPos, lastPos + nCirc).
//   * wordPops[firstPos & (nCirc - 1)] and wordPops[lastPos & (nCirc - 1)] are
//     nonzero.
// * For each i in [0, nCirc), wordPops[i] = # of nonzero words in Row(i).
func (b *Bitmap) CheckPanic(tag string) {
	nCirc := b.NCirc()
	mask := nCirc - 1
	if b.lastPos == -1 {
		if b.firstPos != FirstPosEmpty {
			log.Panicf("b.firstPos = %d, b.lastPos = -1, tag: %s", b.firstPos, tag)
		}
	} else {
		if b.firstPos > b.lastPos {
			log.Panicf("b.firstPos = %d, b.lastPos = %d, tag: %s", b.firstPos, b.lastPos, tag)
		}
		if b.firstPos+mask <= b.lastPos {
			log.Panicf("b.firstPos = %d, b.lastPos = %d, tag: %s", b.firstPos, b.lastPos, tag)
		}
		circFirst := b.firstPos & mask
		circLast := b.lastPos & mask
		if b.wordPops[circFirst] == 0 {
			log.Printf("b.firstPos = %d, b.wordPops[] zero", b.firstPos)
			for i := b.firstPos + 1; i <= b.lastPos; i++ {
				log.Printf("b.wordPops[%d]: %d", i, b.wordPops[i&mask])
			}
			log.Panicf("tag: %s", tag)
		}
		if b.wordPops[circLast] == 0 {
			for i := b.firstPos; i < b.lastPos; i++ {
				log.Printf("b.wordPops[%d]: %d", i, b.wordPops[i&mask])
			}
			log.Panicf("b.lastPos = %d, b.wordPops[] zero, tag: %s", b.lastPos, tag)
		}
	}
	for i := bi.PosType(0); i != nCirc; i++ {
		nz := 0
		row := b.Row(i)
		for j := bi.PosType(0); j != b.rowWidth; j++ {
			if row[j] != 0 {
				nz++
			}
		}
		if byte(nz) != b.wordPops[i] {
			log.Panicf("b.firstPos = %d, b.lastPos = %d; b.wordPops[%d] out of sync (%d, %d expected); tag: %s", b.firstPos, b.lastPos, i, b.wordPops[i], nz, tag)
		}
	}
}
