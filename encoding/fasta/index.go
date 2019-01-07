package fasta

import (
	"bufio"
	"bytes"
	"io"
	"strings"

	"github.com/grailbio/base/errors"
	"github.com/grailbio/base/tsv"
)

// GenerateIndex generates an index (*.fai) from FASTA.  The index can be later
// passed to NewIndexed() to random-access the FASTA file quickly.
//
// The index format is defined by "samtool faidx"
// (http://www.htslib.org/doc/faidx.html).
func GenerateIndex(out io.Writer, in io.Reader) (err error) {
	var (
		tsvOut      = tsv.NewWriter(out)
		r           = bufio.NewReader(in)
		seqName     string
		seqStartOff int64
		totalBases  int
		lineBases   int
		lineWidth   int
		cumByte     int64
		eof         bool
	)

	setErr := func(e error) {
		if e != nil && err == nil {
			err = e
		}
	}
	flush := func() {
		tsvOut.WriteString(seqName)
		tsvOut.WriteInt64(int64(totalBases))
		tsvOut.WriteInt64(seqStartOff)
		tsvOut.WriteInt64(int64(lineBases))
		tsvOut.WriteInt64(int64(lineWidth))
		setErr(tsvOut.EndLine())
	}
	for !eof && err == nil {
		fullLine, e := r.ReadBytes('\n')
		if e == io.EOF { // Process fullLine, then exit the loop
			eof = true
		} else if e != nil {
			setErr(e)
		}
		cumByte += int64(len(fullLine))
		line := bytes.TrimRight(fullLine, "\r\n")
		if len(line) == 0 {
			continue
		}
		if line[0] == '>' { // Start a new sequence.
			if lineWidth != 0 {
				if seqName == "" {
					setErr(errors.E("malformed FASTA file"))
				}
				flush()
			}
			seqName = strings.Split(string(line[1:]), " ")[0]
			seqStartOff = cumByte
			lineWidth = 0
			lineBases = 0
			totalBases = 0
			continue
		}
		if lineWidth == 0 {
			lineWidth = len(fullLine)
			lineBases = len(line)
		}
		totalBases += len(line)
	}
	flush()
	setErr(tsvOut.Flush())
	if cumByte == 0 {
		setErr(errors.E("empty FASTA file"))
	}
	return
}
