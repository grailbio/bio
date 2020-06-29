// Package fasta contains code for parsing (optionally indexed) FASTA files.
// See http://www.htslib.org/doc/faidx.html.  Briefly, FASTA files consist of a
// number of named sequences that may be interrupted by newlines.  For example:
//
// >chr7
// ACGTAC
// GAGGAC
// GCG
// >chr8
// ACGT
//
// Note: Sequence names are defined to be the stretch of characters excluding
// spaces immediately after '>'.  Any text appear after a space are ignored.
// For example, '>chr1 A viral sequence' becomes 'chr1'.
package fasta

import (
	"bufio"
	"bytes"
	"fmt"
	"io"

	"github.com/grailbio/bio/biosimd"
	"github.com/pkg/errors"
)

const (
	mib            = 1024 * 1024
	bufferInitSize = 300 * mib
)

// Fasta represents FASTA-formatted data, consisting of a set of named
// sequences.
type Fasta interface {
	// Get returns a substring of the given sequence name at the given
	// coordinates, which are treated as a 0-based half-open interval
	// [start, end). Get is thread-safe.
	Get(seqName string, start, end uint64) (string, error)

	// Len returns the length of the given sequence.
	Len(seqName string) (uint64, error)

	// SeqNames returns the names of all sequences, in the order of appearance in
	// the FASTA file.
	SeqNames() []string
}

type Encoding byte

const (
	// RawASCII encoding preserves the original bytes, including case.
	RawASCII Encoding = iota
	// CleanASCII encoding capitalizes all lowercase 'a'/'c'/'g'/'t', and
	// converts all non-ACGT characters to 'N'.
	CleanASCII
	// Seq8 encoding is 'A'/'a' = 1, 'C'/'c' = 2, 'G'/'g' = 4, 'T'/'t' = 8,
	// anything else = 15.  This plays well with BAM/PAM files.
	Seq8
	// TODO(cchang): Add 'Base5' encoding, where 'A'/'a' = 0, 'C'/'c' = 1,
	// 'G'/'g' = 2, 'T'/'t' = 3, anything else = 4.
	EncodingLimit
)

type opts struct {
	Enc   Encoding
	Index []byte
}

// Opt is an optional argument to New, NewIndexed.
type Opt func(*opts)

// OptClean specifies returned FASTA sequences should be cleaned as described
// in biosimd.CleanASCIISeq*.  It is equivalent to OptEncoding(CleanASCII).
func OptClean(o *opts) {
	if o.Enc != RawASCII {
		panic("fasta.OptClean: multiple encodings specified")
	}
	o.Enc = CleanASCII
}

// OptEncoding specifies the encoding of the in-memory FASTA sequences.
func OptEncoding(e Encoding) Opt {
	return func(o *opts) {
		if o.Enc != RawASCII {
			panic("fasta.OptEncoding: multiple encodings specified")
		}
		if o.Enc >= EncodingLimit {
			panic("fasta.OptEncoding: invalid encoding value")
		}
		o.Enc = e
	}
}

// OptIndex makes New read FASTA file with a provided index, like NewIndexed.
// Unlike NewIndexed, New with OptIndex is optimized for reading all sequences
// in the FASTA file rather than a small, random subset. Callers that plan to
// read many or all FASTA sequences should use this (though as always, profile
// in your application).
func OptIndex(index []byte) Opt {
	return func(o *opts) {
		o.Index = index
	}
}

func makeOpts(userOpts ...Opt) opts {
	var parsedOpts opts
	for _, userOpt := range userOpts {
		userOpt(&parsedOpts)
	}
	return parsedOpts
}

type fasta struct {
	seqs     map[string]string
	seqNames []string
}

// New creates a new Fasta that holds all the FASTA data from the given reader
// in memory. Pass OptIndex, if possible, to read much faster.
func New(r io.Reader, opts ...Opt) (Fasta, error) {
	parsedOpts := makeOpts(opts...)
	if len(parsedOpts.Index) == 0 {
		return newEagerUnindexed(r, parsedOpts)
	}
	index, err := parseIndex(bytes.NewReader(parsedOpts.Index))
	if err != nil {
		return nil, err
	}
	return newEagerIndexed(r, index, parsedOpts)
}

func newEagerUnindexed(r io.Reader, parsedOpts opts) (Fasta, error) {
	f := &fasta{seqs: make(map[string]string)}
	scanner := bufio.NewScanner(r)
	scanner.Buffer(nil, bufferInitSize)
	var seqName string
	// We don't use strings.Builder here, since that would force us to perform
	// unsafe string -> []byte -> string conversions at the end.
	seqBuf := make([]byte, 0, bufferInitSize)
	for scanner.Scan() {
		line := scanner.Bytes()
		if len(line) == 0 {
			continue
		}
		if line[0] == '>' { // Start a new sequence.
			if len(seqBuf) != 0 { // We need to store the previous sequence first.
				if seqName == "" {
					return nil, errors.Errorf("malformed FASTA file")
				}
				if parsedOpts.Enc == CleanASCII {
					biosimd.CleanASCIISeqInplace(seqBuf)
				} else if parsedOpts.Enc == Seq8 {
					biosimd.ASCIIToSeq8Inplace(seqBuf)
				}
				f.seqs[seqName] = string(seqBuf)
				f.seqNames = append(f.seqNames, seqName)
				seqBuf = seqBuf[:0]
			}
			seqName = string(bytes.SplitN(line[1:], []byte{' '}, 2)[0])
		} else {
			seqBuf = append(seqBuf, line...)
		}
	}
	if scanner.Err() != nil {
		return nil, errors.Wrap(scanner.Err(), "couldn't read FASTA data")
	}
	if parsedOpts.Enc == CleanASCII {
		biosimd.CleanASCIISeqInplace(seqBuf)
	} else if parsedOpts.Enc == Seq8 {
		biosimd.ASCIIToSeq8Inplace(seqBuf)
	}
	f.seqs[seqName] = string(seqBuf)
	f.seqNames = append(f.seqNames, seqName)
	return f, nil
}

// Get implements Fasta.Get().
func (f *fasta) Get(seqName string, start, end uint64) (string, error) {
	s, ok := f.seqs[seqName]
	if !ok {
		return "", errors.Errorf("sequence not found: %s", seqName)
	}
	if end <= start {
		return "", fmt.Errorf("start must be less than end")
	}
	if start < 0 || end > uint64(len(s)) {
		return "", errors.Errorf("invalid query range %d - %d for sequence %s with length %d",
			start, end, seqName, len(s))
	}
	return s[start:end], nil
}

// Len implements Fasta.Len().
func (f *fasta) Len(seq string) (uint64, error) {
	s, ok := f.seqs[seq]
	if !ok {
		return 0, errors.Errorf("sequence not found: %s", seq)
	}
	return uint64(len(s)), nil
}

// SeqNames implements Fasta.SeqNames().
func (f *fasta) SeqNames() []string {
	return f.seqNames
}
