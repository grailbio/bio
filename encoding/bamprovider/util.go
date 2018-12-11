package bamprovider

import (
	"fmt"

	gbam "github.com/grailbio/bio/encoding/bam"
	"github.com/grailbio/hts/sam"
)

// RefByName finds a sam.Reference with the given name. It returns nil if a
// reference is not found.
func RefByName(h *sam.Header, refName string) *sam.Reference {
	for _, ref := range h.Refs() {
		if ref.Name() == refName {
			return ref
		}
	}
	return nil
}

// NewRefIterator creates an iterator for half-open range [refName:start,
// refName:limit). Start and limit are both base zero.  The iterator will yield
// reads whose start positions are in the given range.
func NewRefIterator(p Provider, refName string, start, limit int) Iterator {
	h, err := p.GetHeader()
	if err != nil {
		return NewErrorIterator(err)
	}
	ref := RefByName(h, refName)
	if ref == nil {
		return NewErrorIterator(fmt.Errorf("bamprovider.NewRefIterator: reference '%s' not found", refName))
	}
	shard := gbam.Shard{
		StartRef: ref,
		EndRef:   ref,
		Start:    start,
		End:      limit,
	}
	return p.NewIterator(shard)
}
