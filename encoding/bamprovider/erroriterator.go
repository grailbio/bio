package bamprovider

import (
	"github.com/biogo/hts/sam"
)

type errorIterator struct {
	err error
}

func (i *errorIterator) Scan() bool          { return false }
func (i *errorIterator) Record() *sam.Record { panic("shall not be called") }
func (i *errorIterator) Err() error          { return i.err }
func (i *errorIterator) Close() error        { return i.err }

// NewErrorIterator creates an Iterator that yields no record and returns "err"
// in Err and Close.
func NewErrorIterator(err error) Iterator {
	return &errorIterator{err: err}
}
