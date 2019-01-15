package fusion

// Stitcher stitches two reads (R1 and R2) and produces a Fragment.  It can be
// used to stitch multiple read pairs. Thread compatible.
type Stitcher struct {
	// config params
	kmerLength            int
	lowComplexityFraction float64

	// temp states
	kmerizer  *kmerizer
	r2KmerMap map[Kmer]Pos
	freePool  []Fragment
}

// NewStitcher creates a new stitcher. kmerLength and lowComplexityFraction
// should be copied from the counterparts in Opts.
func NewStitcher(kmerLength int, lowComplexityFraction float64) *Stitcher {
	return &Stitcher{
		kmerLength:            kmerLength,
		lowComplexityFraction: lowComplexityFraction,
		kmerizer:              newKmerizer(kmerLength),
		r2KmerMap:             map[Kmer]Pos{},
	}
}

// FreeFragment puts the fragment in a freepool. The caller must not retain any
// reference to the fragment after the call.  The future calls to Stitch will
// use fragments in the freepool.
func (s *Stitcher) FreeFragment(f Fragment) {
	if len(s.freePool) > 4 {
		return
	}
	f = Fragment{kmers: f.kmers[:0]}
	s.freePool = append(s.freePool, f)
}

func (s *Stitcher) allocFragment() Fragment {
	if l := len(s.freePool); l > 0 {
		f := s.freePool[l-1]
		s.freePool = s.freePool[:l-1]
		return f
	}
	return Fragment{}
}

// Stitch tries to find a segment shared between the readpair and combine them
// into once sequence.  Name should be the R1 name from the FASTQ.
func (s *Stitcher) Stitch(name, r1Seq, r2Seq string, stats *Stats) Fragment {
	frag := s.allocFragment()
	frag.Name = name
	stitchedSeq, ok := s.tryStitch(r1Seq, r2Seq)
	switch {
	case !ok:
		frag.R1Seq = r1Seq
		frag.R2Seq = r2Seq
	case IsLowComplexity(stitchedSeq, s.lowComplexityFraction):
		stats.LowComplexityReadsStitched++
		stats.Stitched++
		frag.R1Seq = "N"
	default:
		stats.Stitched++
		frag.R1Seq = stitchedSeq
	}
	if frag.R1Seq == "N" && frag.R2Seq == "" {
		// shortcurcuit for performance
		return frag
	}
	s.kmerizer.Reset(frag.R1Seq)
	for s.kmerizer.Scan() {
		k := s.kmerizer.Get()
		if k.pos.ReadType() == R2 {
			panic(k)
		}
		// TODO(saito) fix memory consumption
		frag.kmers = append(frag.kmers, k)
	}
	if len(frag.R2Seq) > 0 {
		s.kmerizer.Reset(frag.R2Seq)
		for s.kmerizer.Scan() {
			k := s.kmerizer.Get()
			k.pos = newR2Pos(k.pos)
			frag.kmers = append(frag.kmers, k)
		}
	}
	return frag
}

func (s *Stitcher) tryStitch(r1Seq, r2Seq string) (string, bool) {
	if r1Seq == "" || r2Seq == "" {
		return "", false
	}
	for k := range s.r2KmerMap {
		delete(s.r2KmerMap, k)
	}
	s.kmerizer.Reset(r2Seq)
	for s.kmerizer.Scan() {
		k := s.kmerizer.Get()
		s.r2KmerMap[k.forward] = k.pos
	}

	s.kmerizer.Reset(r1Seq)
	for s.kmerizer.Scan() {
		k1 := s.kmerizer.Get()
		r1Pos := k1.pos
		r2Pos, ok := s.r2KmerMap[k1.forward]
		if !ok {
			continue
		}
		// Normalize the positions by setting origin at the start of the shared kmer.
		//
		//		 r1start 			 	r1pos					 	r1end
		//     |              |               |
		// R1  ===============|--kmer-|=======>
		//
		// R2          <======|--kmer-|==============
		//              |     |										 	 |
		//							|			r2pos									 r2end
		//              |
		//							r2start
		// That is, the coordinates used in this block is adjusted so that r1Pos and
		// r2Pos become zero.
		r1Start, r1End := Pos(0)-r1Pos, Pos(len(r1Seq))-r1Pos
		r2Start, r2End := Pos(0)-r2Pos, Pos(len(r2Seq))-r2Pos

		overlapStart := maxPos(r1Start, r2Start)
		overlapEnd := minPos(r1End, r2End)
		overlap := posSpan(overlapEnd, overlapStart)
		if overlap < s.kmerLength {
			panic("impossible")
		}
		// Check that the overlapped bases are similar enough.
		r1Off, r2Off := int(overlapStart+r1Pos), int(overlapStart+r2Pos)
		hd := hammingDistance(r1Seq[r1Off:r1Off+overlap], r2Seq[r2Off:r2Off+overlap])
		if float64(hd)/float64(overlap) > 0.1 {
			continue
		}

		// Create a single sequence that consists of the overlapped parts.
		stitchStart := minPos(r1Start, r2Start)
		stitchEnd := maxPos(r1End, r2End)
		switch {
		case stitchStart == r1Start && stitchEnd == r1End:
			// (1) R1 covers R2.
			// R1  ===============|--kmer-|==============>
			// R2           <=====|--kmer-|=======
			//     <---------------stitch-------->
			return r1Seq[r1Start+r1Pos : r2End+r1Pos], true
		case stitchStart == r2Start && stitchEnd == r2End: // r2 covers r1
			// (2) R2 covers R1.
			// R1           ======|--kmer-|=======>
			// R2  <==============|--kmer-|==============
			//              <------stitch--------------->
			return r2Seq[r1Start+r2Pos : r2End+r2Pos], true
		case stitchStart == r1Start:
			// (3)
			// R1  ===============|--kmer-|====>
			// R2           <=====|--kmer-|============
			//     <------stitch---------------------->
			return r1Seq + r2Seq[r1End+r2Pos:r2End+r2Pos], true
		default:
			// (4)
			// R1           ======|--kmer-|============>
			// R2  <==============|--kmer-|====
			//              <------stitch----->
			return r1Seq[:r2End+r1Pos], true
		}
	}
	return "", false
}

func hammingDistance(s1, s2 string) int {
	l := len(s1)
	if l != len(s2) {
		panic("hamming")
	}
	d := 0
	for i := 0; i < l; i++ {
		if s1[i] != s2[i] {
			d++
		}
	}
	return d
}
