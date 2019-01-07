package fastq

import (
	"bufio"
	"errors"
	"io"
)

var (
	// ErrShort is returned when a truncated FASTQ file is encountered.
	ErrShort = errors.New("short FASTQ file")
	// ErrInvalid is returned when an invalid FASTQ file is encountered.
	ErrInvalid = errors.New("invalid FASTQ file")
	// ErrDiscordant is returned when two underlying FASTQ files are discordant.
	ErrDiscordant = errors.New("discordant FASTQ pairs")
)

// A Read is a FASTQ read, comprising an ID, sequence, line 3
// ("unknown"), and a quality string.
type Read struct {
	ID, Seq, Unk, Qual string
}

// Trim cuts the read and quality lengths to at most n.
func (r *Read) Trim(n int) {
	r.Seq = r.Seq[:n]
	r.Qual = r.Qual[:n]
}

var errEOF = errors.New("eof")

// Scanner provides a convenient interface for reading FASTQ read
// data. The Scan method returns the next read, returning a boolean
// indicating whether the read succeeded. Scanners are not
// threadsafe.
//
// Scanner performs some validation: it requires ID lines to begin
// with "@" and that line 3 begins with "+", but does not perform
// further validation (e.g., seq/qual being of equal length,
// containing only data in range, etc.)
type Scanner struct {
	b      *bufio.Scanner
	err    error
	fields Field
}

// Field enumerates FASTQ fields. It is used to specify fields to read in
// NewScanner.
type Field uint

const (
	// ID causes the Read.ID field to be filled
	ID Field = 1 << iota
	// Seq causes the Read.Seq field to be filled
	Seq
	// Unk causes the Read.Unk field to be filled
	Unk
	// Qual causes the Read.Unk field to be filled
	Qual
	// All equals ID|Seq|Unk|Qual.
	All = ID | Seq | Unk | Qual
)

// NewScanner constructs a new Scanner that reads raw FASTQ data from the
// provided reader. Fields is a bitset of the fields to read. A typical value
// would be All or ID|Seq|Qual.
func NewScanner(r io.Reader, fields Field) *Scanner {
	return &Scanner{b: bufio.NewScanner(r), fields: fields}
}

// Scan the next read into the provided read. Scan returns a boolean
// indicating whether the scan succeeded. Once Scan returns false, it
// never returns true again. Upon completion, the user should check
// the Err method to determine whether scanning stopped because of an
// error or because the end of the stream was reached.
func (f *Scanner) Scan(read *Read) bool {
	if f.err != nil {
		return false
	}
	if !f.b.Scan() {
		if f.err = f.b.Err(); f.err == nil {
			f.err = errEOF
		}
		return false
	}
	id := f.b.Bytes()
	if len(id) == 0 || id[0] != '@' {
		f.err = ErrInvalid
		return false
	}
	if f.fields&ID != 0 {
		read.ID = string(id)
	}
	if !f.scan() {
		return false
	}
	if f.fields&Seq != 0 {
		read.Seq = f.b.Text()
	}
	if !f.scan() {
		return false
	}
	unk := f.b.Bytes()
	if len(unk) == 0 || unk[0] != '+' {
		f.err = ErrInvalid
		return false
	}
	if f.fields&Unk != 0 {
		read.Unk = string(unk)
	}
	if !f.scan() {
		return false
	}
	if f.fields&Qual != 0 {
		read.Qual = f.b.Text()
	}
	return true
}

func (f *Scanner) scan() bool {
	ok := f.b.Scan()
	if !ok {
		if f.err = f.b.Err(); f.err == nil {
			f.err = ErrShort
		}
	}
	return ok
}

// Err returns the scanning error, if any.
func (f *Scanner) Err() error {
	if f.err == errEOF {
		return nil
	}
	return f.err
}

// PairScanner composes a pair of scanners to scan a pair of FASTQ
// streams.
type PairScanner struct {
	r1, r2 *Scanner
	err    error
}

// NewPairScanner creates a new FASTQ pair scanner from the provided
// R1 and R2 readers.
func NewPairScanner(r1, r2 io.Reader, fields Field) *PairScanner {
	return &PairScanner{
		r1: NewScanner(r1, fields),
		r2: NewScanner(r2, fields),
	}
}

// Scan scans the next read pair into r1, r2. Scan returns a boolean
// indicating whether the scan succeeded. Once Scan returns false, it
// never returns true again. Upon completion, the user should check
// the Err method to determine whether scanning stopped because of an
// error or because the end of the stream was reached.
func (p *PairScanner) Scan(r1, r2 *Read) bool {
	ok1 := p.r1.Scan(r1)
	ok2 := p.r2.Scan(r2)
	if ok1 != ok2 {
		p.err = ErrDiscordant
	}
	return ok1 && ok2
}

// Err returns the scanning error, if any. It should be checked
// after Scan returns false.
func (p *PairScanner) Err() error {
	if err := p.r1.Err(); err != nil {
		return err
	}
	if err := p.r2.Err(); err != nil {
		return err
	}
	return p.err
}
