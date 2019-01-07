package fastq

import "io"

var newline = []byte{'\n'}

// Writer is a FASTQ file writer.
type Writer struct {
	w   io.Writer
	err error
}

// NewWriter constructs a new FASTQ writer
// that writes reads to the underlying writer w.
func NewWriter(w io.Writer) *Writer {
	return &Writer{w: w}
}

// Write writes the read r in FASTQ format.
// An error is returned if the write failed.
func (w *Writer) Write(r *Read) error {
	w.writeln(r.ID)
	w.writeln(r.Seq)
	w.writeln(r.Unk)
	w.writeln(r.Qual)
	return w.err
}

func (w *Writer) writeln(line string) {
	if w.err != nil {
		return
	}
	_, w.err = io.WriteString(w.w, line)
	if w.err == nil {
		_, w.err = w.w.Write(newline)
	}
}
