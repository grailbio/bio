package fasta

import (
	"bufio"
	"fmt"
	"io"

	"github.com/grailbio/base/unsafe"
	"github.com/grailbio/bio/biosimd"
)

func newEagerIndexed(fastaR io.Reader, index []indexEntry, parsedOpts opts) (Fasta, error) {
	var entireLen uint64
	entireSeqStarts := make([]uint64, len(index))
	for e, entry := range index {
		entireSeqStarts[e] = entireLen
		entireLen += entry.length
	}
	entire := make([]byte, entireLen)

	var fileOffset uint64
	bufR := bufio.NewReaderSize(fastaR, bufferInitSize)
	for e, entry := range index {
		n, err := bufR.Discard(int(entry.offset - fileOffset))
		fileOffset += uint64(n)
		if err != nil {
			return nil, fmt.Errorf("fasta.NewIndexedAll: seeking seq: %v", err)
		}
		var basesRead uint64
		for basesRead < entry.length {
			// Compute length of the next line (may be partial if it's the last).
			nextBasesRead := basesRead + entry.lineBase
			if nextBasesRead > entry.length {
				nextBasesRead = entry.length
			}
			lineBases := nextBasesRead - basesRead

			// Read the line.
			entireLineStart := entireSeqStarts[e] + basesRead
			entireLine := entire[entireLineStart : entireLineStart+lineBases]
			n, err := io.ReadFull(bufR, entireLine)
			fileOffset += uint64(n)
			if err != nil {
				return nil, fmt.Errorf("fasta.NewIndexedAll: reading: %v", err)
			}
			basesRead += lineBases

			// Skip line terminator(s) unless we're at the end of the sequence.
			if basesRead < entry.length {
				n, err := bufR.Discard(int(entry.lineWidth - entry.lineBase))
				fileOffset += uint64(n)
				if err != nil {
					return nil, fmt.Errorf("fasta.NewIndexedAll: seeking line: %v", err)
				}
			}
		}
	}

	if parsedOpts.Enc == CleanASCII {
		biosimd.CleanASCIISeqInplace(entire)
	} else if parsedOpts.Enc == Seq8 {
		biosimd.ASCIIToSeq8Inplace(entire)
	}

	fa := fasta{
		seqs:     make(map[string]string, len(index)),
		seqNames: make([]string, 0, len(index)),
	}
	for e, entry := range index {
		seqBytes := entire[entireSeqStarts[e] : entireSeqStarts[e]+entry.length]
		fa.seqs[entry.name] = unsafe.BytesToString(seqBytes)
		fa.seqNames = append(fa.seqNames, entry.name)
	}
	return &fa, nil
}
