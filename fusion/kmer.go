package fusion

import (
	"github.com/grailbio/base/simd"
	gunsafe "github.com/grailbio/base/unsafe"
	"github.com/grailbio/bio/biosimd"
)

const (
	invalidKmerBits = uint8(255)
)

var (
	asciiToKmerMap                  [256]uint8
	asciiToReverseComplementKmerMap [256]uint8
)

func init() {
	for i := range asciiToKmerMap {
		asciiToKmerMap[i] = invalidKmerBits
		asciiToReverseComplementKmerMap[i] = invalidKmerBits
	}
	asciiToKmerMap['A'] = 0
	asciiToKmerMap['a'] = 0
	asciiToKmerMap['C'] = 1
	asciiToKmerMap['c'] = 1
	asciiToKmerMap['G'] = 2
	asciiToKmerMap['g'] = 2
	asciiToKmerMap['T'] = 3
	asciiToKmerMap['t'] = 3

	asciiToReverseComplementKmerMap['A'] = 3
	asciiToReverseComplementKmerMap['a'] = 3
	asciiToReverseComplementKmerMap['C'] = 2
	asciiToReverseComplementKmerMap['c'] = 2
	asciiToReverseComplementKmerMap['G'] = 1
	asciiToReverseComplementKmerMap['g'] = 1
	asciiToReverseComplementKmerMap['T'] = 0
	asciiToReverseComplementKmerMap['t'] = 0
}

// Kmer is a compact encoding of a sequence of ACGT, up to 32bases.
type Kmer uint64

// invalidKmer is a sentinel kmer.
const invalidKmer = Kmer(0xffffffffffffffff)

type kmersAtPos struct {
	// Pos is the position in the fragment. its readType may be R1 or R2.
	pos Pos
	// Forward and R-C kmer that encodes the subsequence
	// [Pos,Pos+opts.KmerLength).
	forward, reverseComplement Kmer
}

func (km kmersAtPos) minKmer() Kmer {
	if km.forward < km.reverseComplement {
		return km.forward
	}
	return km.reverseComplement
}

type kmerizer struct {
	kmerLength int
	tmpSeq     []byte
	mask       Kmer // ~0 << (2*kmerLength)

	seq string
	si  int
	cur kmersAtPos
}

func newKmerizer(kmerLength int) *kmerizer {
	return &kmerizer{
		kmerLength: kmerLength,
		mask:       ^(Kmer(0xffffffffffffffff) << Kmer(kmerLength*2 /*2==#bits per base*/)),
	}
}

func asciiToKmer(seq string) Kmer {
	var k Kmer
	for _, ch := range []byte(seq) {
		b := asciiToKmerMap[ch]
		if b == invalidKmerBits {
			return invalidKmer
		}
		k = (k << 2) | Kmer(b)
	}
	return k
}

func nextAmbiguousPosition(seq string, si int) int {
	for i := si; si < len(seq); i++ {
		if asciiToKmerMap[seq[i]] == invalidKmerBits {
			return i
		}
	}
	return len(seq)
}

func (k *kmerizer) Reset(seq string) {
	k.seq = seq
	k.si = 0
}

func (k *kmerizer) Scan() bool {
	if k.si > 0 /*k.cur is set*/ && k.si+k.kmerLength <= len(k.seq) {
		nextCh := k.seq[k.si+k.kmerLength-1]
		if bits := asciiToKmerMap[nextCh]; bits != invalidKmerBits {
			// Fast path. Directly add the 2-bit encoding of "nextCh" to k.cur.forward
			// and k.cur.reverseComplement.
			k.cur.pos = Pos(k.si)
			k.cur.forward = ((k.cur.forward << 2) | Kmer(bits)) & k.mask
			shift := (Kmer(k.kmerLength) - 1) * 2
			k.cur.reverseComplement = (k.cur.reverseComplement >> 2) | (Kmer(asciiToReverseComplementKmerMap[nextCh]) << shift)
			k.si++
			return true
		}
		// Fall through
	}

	for k.si+k.kmerLength <= len(k.seq) {
		forwardStr := k.seq[k.si : k.si+k.kmerLength]
		var forwardKmer, reverseKmer Kmer
		if forwardKmer = asciiToKmer(forwardStr); forwardKmer == invalidKmer {
			k.si = nextAmbiguousPosition(k.seq, k.si) + 1
			continue
		}
		simd.ResizeUnsafe(&k.tmpSeq, k.kmerLength)
		biosimd.ReverseComp8NoValidate(k.tmpSeq, gunsafe.StringToBytes(forwardStr))
		if reverseKmer = asciiToKmer(gunsafe.BytesToString(k.tmpSeq)); reverseKmer == invalidKmer {
			panic("shoulnd't happen")
		}
		k.cur = kmersAtPos{pos: Pos(k.si), forward: forwardKmer, reverseComplement: reverseKmer}
		k.si++
		return true
	}
	return false
}

func (k *kmerizer) Get() kmersAtPos { return k.cur }
