package fasta

import (
	"bufio"
	"fmt"
	"io"
	"sort"
	"strconv"
	"sync"
)

type indexEntry struct {
	length    uint64
	offset    uint64
	lineBase  uint64
	lineWidth uint64
}

type indexedFasta struct {
	seqs      map[string]indexEntry
	seqNames  []string // returned by SeqNames()
	reader    io.ReadSeeker
	bufOff    int64
	buf       []byte // caches file contents starting at bufOff.
	resultBuf []byte // temp for concatenating multi-line sequences.
	mutex     sync.Mutex
}

// NewIndexed creates a new Fasta that can perform efficient random lookups
// using the provided index, without reading the data into memory.
func NewIndexed(fasta io.ReadSeeker, index io.Reader) (Fasta, error) {
	f := &indexedFasta{seqs: make(map[string]indexEntry), reader: fasta}
	scanner := bufio.NewScanner(index)
	scanner.Split(bufio.ScanLines)
	for scanner.Scan() {
		matches := indexRegExp.FindStringSubmatch(scanner.Text())
		if len(matches) != 6 {
			return nil, fmt.Errorf("Invalid index line: %s", scanner.Text())
		}
		ent := indexEntry{}
		ent.length, _ = strconv.ParseUint(matches[2], 10, 64)
		ent.offset, _ = strconv.ParseUint(matches[3], 10, 64)
		ent.lineBase, _ = strconv.ParseUint(matches[4], 10, 64)
		ent.lineWidth, _ = strconv.ParseUint(matches[5], 10, 64)
		f.seqs[matches[1]] = ent
		f.seqNames = append(f.seqNames, matches[1])
	}
	sort.SliceStable(f.seqNames, func(i, j int) bool {
		return f.seqs[f.seqNames[i]].offset < f.seqs[f.seqNames[j]].offset
	})
	return f, nil
}

// FaiToReferenceLengths reads in a fasta fai file and returns a map of
// reference name to reference length. This doesn't require reading in the fasta
// itself.
func FaiToReferenceLengths(index io.Reader) (map[string]uint64, error) {
	newIndexed, err := NewIndexed(nil, index)
	if err != nil {
		return nil, err
	}
	newMap := make(map[string]uint64)
	for _, ref := range newIndexed.SeqNames() {
		refLength, err := newIndexed.Len(ref)
		if err != nil {
			return nil, err
		}
		newMap[ref] = refLength
	}
	return newMap, nil
}

// Len implements Fasta.Len().
func (f *indexedFasta) Len(seqName string) (uint64, error) {
	ent, ok := f.seqs[seqName]
	if !ok {
		return 0, fmt.Errorf("sequence not found in index: %s", seqName)
	}
	return ent.length, nil
}

// Read range [off, off+n) from the underlying fasta file.
func (f *indexedFasta) read(off int64, n int) ([]byte, error) {
	limit := off + int64(n)
	if off < f.bufOff || limit > f.bufOff+int64(len(f.buf)) {
		if newOffset, err := f.reader.Seek(off, io.SeekStart); err != nil || newOffset != off {
			return nil, fmt.Errorf("failed to seek to offset %d: %d, %v", off, newOffset, err)
		}
		bufSize := 8192
		if bufSize < n {
			bufSize = n
		}
		f.resizeBuf(&f.buf, bufSize)
		bytesRead, err := f.reader.Read(f.buf)
		if bytesRead < n {
			return nil, fmt.Errorf("encountered unexpected end of file (bad index? file doesn't end in newline?)")
		}
		if err != nil && err != io.EOF {
			return nil, err
		}
		f.bufOff = off
		f.buf = f.buf[:bytesRead]
		if off < f.bufOff || limit > f.bufOff+int64(len(f.buf)) {
			panic(off)
		}
	}
	return f.buf[off-f.bufOff : limit-f.bufOff], nil
}

func (f *indexedFasta) resizeBuf(buf *[]byte, n int) {
	if cap(*buf) < n {
		*buf = make([]byte, n)
	} else {
		*buf = (*buf)[0:n]
	}
}

// Get implements Fasta.Get().
func (f *indexedFasta) Get(seqName string, start uint64, end uint64) (string, error) {
	f.mutex.Lock()
	defer f.mutex.Unlock()

	if end <= start {
		return "", fmt.Errorf("start must be less than end")
	}
	ent, ok := f.seqs[seqName]
	if !ok {
		return "", fmt.Errorf("sequence not found in index: %s", seqName)
	}
	if end > ent.length {
		return "", fmt.Errorf("end is past end of sequence %s: %d", seqName, ent.length)
	}

	// Start the read at a byte offset allowing for the presence of newline
	// characters.
	charsPerNewline := ent.lineWidth - ent.lineBase
	offset := ent.offset + start + charsPerNewline*(start/ent.lineBase)

	// Figure out how many characters (including newlines) we should read,
	// and read them.
	firstLineBases := ent.lineBase - (start % ent.lineBase)
	newlinesToRead := uint64(0)
	if end-start > firstLineBases {
		newlinesToRead = 1 + (end-start-firstLineBases)/ent.lineBase
	}
	capacity := end - start + newlinesToRead*charsPerNewline

	buffer, err := f.read(int64(offset), int(capacity))
	if err != nil && err != io.EOF {
		return "", err
	}

	// Traverse the bytes we just read and copy the non-newline characters
	// to the result.
	f.resizeBuf(&f.resultBuf, int(end-start))
	linePos := (offset - ent.offset) % ent.lineWidth
	resultPos := 0
	for i := range buffer {
		if linePos < ent.lineBase {
			f.resultBuf[resultPos] = buffer[i]
			resultPos++
		}
		linePos++
		if linePos == ent.lineWidth {
			linePos = 0
		}
	}
	return string(f.resultBuf), nil
}

// SeqNames implements Fasta.SeqNames().
func (f *indexedFasta) SeqNames() []string {
	return f.seqNames
}
