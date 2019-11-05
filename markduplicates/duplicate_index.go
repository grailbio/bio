package markduplicates

import (
	"fmt"
	"regexp"
	"strings"

	"github.com/grailbio/base/log"
	"github.com/grailbio/bio/encoding/bam"
	"github.com/grailbio/bio/umi"
	"github.com/grailbio/bio/util"
	"github.com/grailbio/hts/sam"
)

var umiRe = regexp.MustCompile(`([ACGTNacgtn]+)\+([ACGTNacgtn]+)`)

// If the set has any pairs, the primary will be in pairs[0],
// otherwise, the primary will be in singles[0].  Each name in
// opticals will also be in pairs.  This is the externally visible
// data structure.
type duplicateSet struct {
	pairs     []string
	singles   []string
	opticals  []string
	corrected map[string]string
}

type DuplicateEntry interface {
	Name() string
	BaseQScore() int
	FileIdx() uint64
}

type IndexedSingle struct {
	R        *sam.Record
	FileIdx_ uint64
}

// Use this to order two reads in a read pair.  If the refid, pos, and
// orientation all match, then R1 is less than R2.  If everything
// matches except for the read number, then the order does not matter
// for comparing potential positional duplicate pairs because only
// ref, pos, and orientation are compared for determining positional
// duplicates.
func (s *IndexedSingle) lessThan(other IndexedSingle) bool {
	sPos := bam.UnclippedFivePrimePosition(s.R)
	otherPos := bam.UnclippedFivePrimePosition(other.R)
	sOrientation := orientationByteSingle(bam.IsReversedRead(s.R))
	otherOrientation := orientationByteSingle(bam.IsReversedRead(other.R))

	return s.R.Ref.ID() < other.R.Ref.ID() ||
		(s.R.Ref.ID() == other.R.Ref.ID() && sPos < otherPos) ||
		(s.R.Ref.ID() == other.R.Ref.ID() && sPos == otherPos && sOrientation < otherOrientation) ||
		(s.R.Ref.ID() == other.R.Ref.ID() && sPos == otherPos && sOrientation == otherOrientation && s.R.Flags&sam.Read1 != 0)
}

func (s IndexedSingle) Name() string {
	return s.R.Name
}

func (s IndexedSingle) BaseQScore() int {
	return baseQScore(s.R)
}

func (s IndexedSingle) FileIdx() uint64 {
	return s.FileIdx_
}

type IndexedPair struct {
	Left  IndexedSingle
	Right IndexedSingle
}

func (p IndexedPair) Name() string {
	return p.Left.R.Name
}

func (p IndexedPair) BaseQScore() int {
	score := baseQScore(p.Left.R)
	if p.Right.R != nil {
		score += baseQScore(p.Right.R)
	}
	return score
}

func (p IndexedPair) FileIdx() uint64 {
	return p.Left.FileIdx_
}

func (p IndexedPair) GetR1R2() (r1 *sam.Record, r2 *sam.Record) {
	if bam.IsRead1(p.Left.R) {
		return p.Left.R, p.Right.R
	}
	return p.Right.R, p.Left.R
}

type IntermediateDuplicateSet struct {
	Pairs     []DuplicateEntry
	Singles   []DuplicateEntry
	Corrected map[string]string // Maps read name to corrected UMI pair: "GAC+GAG"
}

type umiKey struct {
	leftRefId   int
	leftPos     int
	rightRefId  int
	rightPos    int
	Orientation Orientation
	Strand      strand
	leftUmi     string
	rightUmi    string
}

func (k *umiKey) isSingle() bool {
	return k.Orientation == f || k.Orientation == r
}

func (k *umiKey) distance(other *umiKey) int {
	if k.isSingle() != other.isSingle() {
		log.Fatalf("Compared single key with pair key %v %v", k, other)
	}
	dist := 0
	if len(k.leftUmi) > 0 {
		dist += util.Levenshtein(k.leftUmi, other.leftUmi, "", "")
	}
	if len(k.rightUmi) > 0 {
		dist += util.Levenshtein(k.rightUmi, other.rightUmi, "", "")
	}
	return dist
}

// duplicateIndex contains the logic used to resolve duplicates.
type duplicateIndex struct {
	worker           int
	entries          map[duplicateKey][]DuplicateEntry
	readGroupLibrary map[string]string
	queue            []*duplicateSet
	umiCorrector     *umi.SnapCorrector
	opts             *Opts
	bagProcessors    []BagProcessor
	startedRemoving  bool
}

// newDuplicateIndex returns a duplicateIndex with the given
// parameters.  The duplicateIndex prints the worker id parameter when
// logging, so that we can separate the output from each of the
// concurrent worker threads.
func newDuplicateIndex(
	worker int,
	header *sam.Header,
	readGroupLibrary map[string]string,
	opts *Opts,
	umiCorrector *umi.SnapCorrector) *duplicateIndex {
	di := &duplicateIndex{
		worker:           worker,
		entries:          make(map[duplicateKey][]DuplicateEntry),
		readGroupLibrary: readGroupLibrary,
		queue:            make([]*duplicateSet, 0),
		umiCorrector:     umiCorrector,
		opts:             opts,
	}

	for i := range opts.BagProcessorFactories {
		di.bagProcessors = append(di.bagProcessors, opts.BagProcessorFactories[i].Create())
	}
	return di
}

// insert a record that is mate-unmapped, sometimes called a singleton.
func (d *duplicateIndex) insertSingleton(r *sam.Record, fileIdx uint64) {
	if d.startedRemoving {
		log.Fatalf("cannot insert after started removing")
	}

	fivePosition := bam.UnclippedFivePrimePosition(r)
	orientation := orientationByteSingle(bam.IsReversedRead(r))
	var strand strand
	if d.opts.StrandSpecific {
		strand = r1Strand(r)
	}
	key := duplicateKey{r.Ref.ID(), fivePosition, -1, -1, orientation, strand}
	d.entries[key] = append(d.entries[key], IndexedSingle{r, fileIdx})
}

// insert a read pair.  a and b need not be in any particular order;
// insertPair will order them in a canonical order internally.
func (d *duplicateIndex) insertPair(a, b *sam.Record, aFileIdx, bFileIdx uint64) {
	if d.startedRemoving {
		log.Fatalf("cannot insert after started removing")
	}

	aIndexed := IndexedSingle{a, aFileIdx}
	bIndexed := IndexedSingle{b, bFileIdx}
	var left, right IndexedSingle
	if aIndexed.lessThan(bIndexed) {
		left = IndexedSingle{a, aFileIdx}
		right = IndexedSingle{b, bFileIdx}
	} else {
		left = IndexedSingle{b, bFileIdx}
		right = IndexedSingle{a, aFileIdx}
	}

	// Update duplicate set.
	var strand strand
	if d.opts.StrandSpecific {
		strand = r1Strand(a)
	}
	key := duplicateKey{
		left.R.Ref.ID(), bam.UnclippedFivePrimePosition(left.R),
		right.R.Ref.ID(), bam.UnclippedFivePrimePosition(right.R),
		orientationBytePair(bam.IsReversedRead(left.R), bam.IsReversedRead(right.R)),
		strand,
	}
	d.entries[key] = append(d.entries[key], IndexedPair{left, right})
}

func ChoosePrimary(entries []DuplicateEntry) int {
	bestIndex := -1
	bestScore := -1
	bestFileIdx := uint64(0)
	for i, entry := range entries {
		currentScore := entry.BaseQScore()
		// Choose primary using score, and break ties using the fileIdx of left.
		if bestIndex < 0 || currentScore > bestScore || (currentScore == bestScore && entry.FileIdx() < bestFileIdx) {
			bestIndex = i
			bestScore = currentScore
			bestFileIdx = entry.FileIdx()
		}
	}
	return bestIndex
}

// The user should call computeDupSets() after inserting all
// singletons and pairs with insertSingle() or insertPair(), and
// before calling nextDupSet().  Do not call insertSingle() or
// insertPair() after calling removeDupSet().
//
//  1) Create an intermediate IntermediateDuplicateSet which contains pairs and singles.
//     Currently this may contain
//       a) exact position matches
//       b) exact position matches + exact match umi.
//     In the future, this may contain matches like fuzzy umi matches.
//  2) Decides the primary, and computes opticals based on the IntermediateDuplicateSet groups.
func (d *duplicateIndex) computeDupSets(metrics *MetricsCollection) {
	d.startedRemoving = true

	// Create groups according to opts.
	var groups []*IntermediateDuplicateSet
	if d.opts.UseUmis {
		groups = d.groupByPositionAndUmi()
	} else {
		groups = d.groupByPosition()
	}

	for _, bagProcessor := range d.bagProcessors {
		groups = bagProcessor(groups)
	}

	// Choose primary & compute opticals.
	for _, g := range groups {
		set := duplicateSet{
			corrected: g.Corrected,
		}

		if len(g.Pairs) > 0 {
			bestIndex := ChoosePrimary(g.Pairs)
			set.pairs = append(set.pairs, g.Pairs[bestIndex].(IndexedPair).Left.R.Name)
			for i, pair := range g.Pairs {
				if i != bestIndex {
					set.pairs = append(set.pairs, pair.(IndexedPair).Left.R.Name)
				}
			}
			for _, single := range g.Singles {
				set.singles = append(set.singles, single.(IndexedSingle).R.Name)
			}
			if d.opts.OpticalDetector != nil {
				set.opticals = d.opts.OpticalDetector.Detect(d.readGroupLibrary, g.Pairs, bestIndex)
			}
			if len(d.opts.OpticalHistogram) > 0 {
				addOpticalDistances(d.opts, d.readGroupLibrary, g.Pairs, metrics)
			}
		} else {
			bestIndex := ChoosePrimary(g.Singles)
			set.singles = append(set.singles, g.Singles[bestIndex].(IndexedSingle).R.Name)
			for i, single := range g.Singles {
				if i != bestIndex {
					set.singles = append(set.singles, single.(IndexedSingle).R.Name)
				}
			}
		}

		d.queue = append(d.queue, &set)
	}
}

func (d *duplicateIndex) groupByPosition() []*IntermediateDuplicateSet {
	getDupSingles := func(refId, pos int, orientation Orientation, strand strand) []DuplicateEntry {
		k := duplicateKey{refId, pos, -1, -1, orientation, strand}
		singles, ok := d.entries[k]
		if ok {
			delete(d.entries, k)
			return singles
		}
		return []DuplicateEntry{}
	}

	groups := make([]*IntermediateDuplicateSet, 0)

	for k, duplicates := range d.entries {
		if !k.isSingle() {
			singles := make([]DuplicateEntry, 0)
			if !d.opts.SeparateSingletons {
				singles = append(getDupSingles(k.leftRefId, k.leftPos, leftOrientation(k.Orientation), k.Strand),
					getDupSingles(k.rightRefId, k.rightPos, rightOrientation(k.Orientation), k.Strand)...)
			}

			groups = append(groups, &IntermediateDuplicateSet{
				Pairs:   duplicates,
				Singles: singles,
			})
			delete(d.entries, k)
		}
	}

	for k, duplicates := range d.entries {
		if k.isSingle() {
			groups = append(groups, &IntermediateDuplicateSet{
				Singles: duplicates,
			})
			delete(d.entries, k)
		}
	}
	return groups
}

// Note: a singleton will match against a pair if just the singleton's
// one umi matches the relevant read in the pair, even if the
// singleton's read name contains two umis.
func (d *duplicateIndex) groupByPositionAndUmi() []*IntermediateDuplicateSet {

	scavenge := func(scavengeCandidates, knownUmis map[umiKey]bool, umiToGroup map[umiKey][]DuplicateEntry) {
		for key := range scavengeCandidates {
			numCloseEnough := 0
			var closeEnough umiKey
			for knownKey := range knownUmis {
				if key.distance(&knownKey) <= d.opts.ScavengeUmis {
					closeEnough = knownKey
					numCloseEnough++
					if numCloseEnough > 1 {
						break
					}
				}
			}

			// If there is exactly one knownUmi bag that is within
			// the scavenge distance, then combine those two bags.
			if numCloseEnough == 1 {
				log.Debug.Printf("scavenge success for %v to %v", key, closeEnough)
				umiToGroup[closeEnough] = append(umiToGroup[closeEnough], umiToGroup[key]...)
				delete(umiToGroup, key)
			} else {
				// Note that this does not attempt to error correct a scavengeCandiate against another scavengeCandiate.
				// We could add that later if we think it would be helpful.
				if log.At(log.Debug) {
					for _, s := range umiToGroup[key] {
						log.Debug.Printf("could not scavenge %s", s.(IndexedSingle).R.Name)
					}
				}
			}
		}
	}

	// For each position-based group, further split pairs and singles by umi.
	umiToGroup := map[umiKey][]DuplicateEntry{}

	for k, entries := range d.entries {
		scavengeCandidates := map[umiKey]bool{}
		knownUmis := map[umiKey]bool{}

		for _, e := range entries {
			leftUmi, rightUmi, fullyCorrected, correctedSome := d.tryCorrectUmis(e)
			// If the resulting UMIs are both known umis, then save the corrected umi values.
			if d.opts.TagDups && fullyCorrected && correctedSome {
				log.Debug.Printf("snap correcting %s", e.Name())
			}

			// Put each pair into the duplicate umi map.
			key := umiKey{k.leftRefId, k.leftPos, k.rightRefId, k.rightPos, k.Orientation,
				k.Strand, leftUmi, rightUmi}
			umiToGroup[key] = append(umiToGroup[key], e)

			// remember which keys were not fully corrected.
			if !fullyCorrected {
				scavengeCandidates[key] = true
			} else {
				knownUmis[key] = true
			}
		}

		if d.opts.ScavengeUmis > 0 {
			// Attempt to match scavengeCandidates against bags that have known umis.
			scavenge(scavengeCandidates, knownUmis, umiToGroup)
		}
		delete(d.entries, k)
	}

	getDupSingles := func(refId, pos int, orientation Orientation, strand strand, umi string) []DuplicateEntry {
		k := umiKey{refId, pos, -1, -1, orientation, strand, umi, ""}
		singles, ok := umiToGroup[k]
		if ok {
			delete(umiToGroup, k)
		}
		return singles
	}

	// attach the relevant corrections and return a IntermediateDuplicateSet.
	createDupSetInternal := func(key umiKey, pairs []DuplicateEntry, singles []DuplicateEntry) *IntermediateDuplicateSet {
		corrected := map[string]string{}
		if d.opts.TagDups {
			for _, p := range pairs {
				left, right, swapped := getCanonicalUmis(p.(IndexedPair))
				if left != key.leftUmi || right != key.rightUmi {
					if swapped {
						corrected[p.Name()] = fmt.Sprintf("%s+%s", key.rightUmi, key.leftUmi)
					} else {
						corrected[p.Name()] = fmt.Sprintf("%s+%s", key.leftUmi, key.rightUmi)
					}
				}
			}
			for _, single := range singles {
				s := single.(IndexedSingle)
				umi, mateUmi, swapped := getCanonicalUmi(s)

				if s.R.Ref.ID() == key.leftRefId && s.R.Pos == key.leftPos &&
					((key.isSingle() && orientationByteSingle(bam.IsReversedRead(s.R)) == key.Orientation) ||
						!key.isSingle() && orientationByteSingle(bam.IsReversedRead(s.R)) == leftOrientation(key.Orientation)) &&
					umi != key.leftUmi {
					// key.leftUmi is the corrected value.
					if swapped {
						corrected[s.Name()] = fmt.Sprintf("%s+%s", mateUmi, key.leftUmi)
					} else {
						corrected[s.Name()] = fmt.Sprintf("%s+%s", key.leftUmi, mateUmi)
					}
				} else if s.R.Ref.ID() == key.rightRefId && s.R.Pos == key.rightPos &&
					((key.isSingle() && orientationByteSingle(bam.IsReversedRead(s.R)) == key.Orientation) ||
						!key.isSingle() && orientationByteSingle(bam.IsReversedRead(s.R)) == rightOrientation(key.Orientation)) &&
					umi != key.rightUmi {
					// key.rightUmi is the corrected value.
					if swapped {
						corrected[s.Name()] = fmt.Sprintf("%s+%s", mateUmi, key.rightUmi)
					} else {
						corrected[s.Name()] = fmt.Sprintf("%s+%s", key.rightUmi, mateUmi)
					}
				}
			}
		}
		return &IntermediateDuplicateSet{
			Pairs:     pairs,
			Singles:   singles,
			Corrected: corrected,
		}
	}

	// The following code is the same as groupByPosition() except for
	// the addition of umis to the lookup key.  It would be nice to
	// use a common piece of code for this.
	groups := make([]*IntermediateDuplicateSet, 0)
	for k, pairs := range umiToGroup {
		if k.isSingle() {
			continue
		}

		// Find singles that match on position and umi.
		singles := make([]DuplicateEntry, 0)
		// Find singles that match on position and umi.
		if !d.opts.SeparateSingletons {
			// Collect matching singles for each read who's umi lacks N.
			if !strings.ContainsAny(k.leftUmi, "Nn") {
				singles = append(singles, getDupSingles(k.leftRefId, k.leftPos, leftOrientation(k.Orientation),
					k.Strand, k.leftUmi)...)
			}
			if !strings.ContainsAny(k.rightUmi, "Nn") {
				singles = append(singles, getDupSingles(k.rightRefId, k.rightPos, rightOrientation(k.Orientation),
					k.Strand, k.rightUmi)...)
			}
		}

		// If either umi contains N, split each pair into a separate
		// group.  The first group should contain any singletons that
		// matched this k.
		if strings.ContainsAny(k.leftUmi, "Nn") || strings.ContainsAny(k.rightUmi, "Nn") {
			for i, p := range pairs {
				groups = append(groups, createDupSetInternal(k, []DuplicateEntry{p}, singles))
				if i == 0 {
					// We attach the singles only to pairs[0].
					singles = []DuplicateEntry{}
				}
			}
		} else {
			groups = append(groups, createDupSetInternal(k, pairs, singles))
		}
		delete(umiToGroup, k)
	}

	for k, singles := range umiToGroup {
		if strings.ContainsAny(k.leftUmi, "Nn") {
			for _, s := range singles {
				groups = append(groups, createDupSetInternal(k, []DuplicateEntry{}, []DuplicateEntry{s}))
			}
		} else {
			groups = append(groups, createDupSetInternal(k, []DuplicateEntry{}, singles))
		}
		delete(umiToGroup, k)
	}
	return groups
}

func (d *duplicateIndex) tryCorrectUmis(e DuplicateEntry) (leftUmi, rightUmi string, fullyCorrected, correctedSome bool) {
	switch v := e.(type) {
	case IndexedPair:
		leftUmi, rightUmi, _ = getCanonicalUmis(v)
		if d.umiCorrector != nil {
			correctedLeftUmi, leftDist, correctedLeft := d.umiCorrector.CorrectUMI(leftUmi)
			correctedRightUmi, rightDist, correctedRight := d.umiCorrector.CorrectUMI(rightUmi)

			leftUmi = correctedLeftUmi
			rightUmi = correctedRightUmi
			fullyCorrected = (leftDist >= 0 && rightDist >= 0)
			correctedSome = (correctedLeft || correctedRight)
		} else {
			fullyCorrected = false
			correctedSome = false
		}
	case IndexedSingle:
		leftUmi, _, _ = getCanonicalUmi(v)
		if d.umiCorrector != nil {
			correctedUmi, dist, corrected := d.umiCorrector.CorrectUMI(leftUmi)

			leftUmi = correctedUmi
			rightUmi = ""
			fullyCorrected = dist >= 0
			correctedSome = corrected
		} else {
			fullyCorrected = false
			correctedSome = false
		}
	}
	return
}

func getUmiField(name string) string {
	idx := strings.LastIndexByte(name, ':')
	if idx < 0 {
		log.Fatalf("Could not parse UMI in qname: %s", name)
	}
	return name[idx:]
}

// getCanonicalUmis returns the 'left' and 'right' umis for a given
// pair.  Even though the pair has a left and right, those left and
// right are not always ordered in a canonical way because that sort
// order relies on R1 and R2 to break the tie when the ref, pos, and
// orientations are equal for both reads in a pair.  In those cases,
// getCanonicalUmis must order the umis canonically, and it does so
// based on this criteria: (refid, pos, orientation, umi) which
// ignores the R1 and R2 flags.  Also returns a boolean that is true
// if leftUmi came from R2.
func getCanonicalUmis(pair IndexedPair) (leftUmi string, rightUmi string, swapped bool) {
	umis := umiRe.FindStringSubmatch(getUmiField(pair.Left.R.Name))
	if umis == nil {
		log.Fatalf("Could not parse UMI in qname: %s", pair.Left.R.Name)
	}

	// If it's a tie based on ref, pos, and orientation, then order by umi value.
	if pair.Left.R.Ref.ID() == pair.Right.R.Ref.ID() &&
		bam.UnclippedFivePrimePosition(pair.Left.R) == bam.UnclippedFivePrimePosition(pair.Right.R) &&
		bam.IsReversedRead(pair.Left.R) == bam.IsReversedRead(pair.Right.R) {
		if strings.Compare(umis[1], umis[2]) < 0 {
			return umis[1], umis[2], false
		}
		return umis[2], umis[1], true
	}

	// Otheriwse keep the left/right order as given by the pair.
	if (pair.Left.R.Flags & sam.Read1) != 0 {
		return umis[1], umis[2], false
	}
	return umis[2], umis[1], true
}

// getCanonicalUmi returns the UMI associated with read, and also the
// UMI associated with the read's mate.  The third return value is
// true if umi is from R2.
func getCanonicalUmi(read IndexedSingle) (umi string, mateUmi string, swapped bool) {
	umis := umiRe.FindStringSubmatch(getUmiField(read.R.Name))
	if umis == nil {
		log.Fatalf("Could not parse UMI in qname: %s", read.R.Name)
	}
	if (read.R.Flags & sam.Read1) != 0 {
		return umis[1], umis[2], false
	}
	return umis[2], umis[1], true
}

// This is the method for outside users.  This will remove and return
// one set of duplicates.  The duplicateSet might be based on a pair
// or a singleton.  If there are no more duplicateSets, returns (nil,
// false).
func (d *duplicateIndex) nextDupSet() (*duplicateSet, bool) {
	if len(d.queue) > 0 {
		var dupSet *duplicateSet
		dupSet, d.queue = d.queue[0], d.queue[1:]
		return dupSet, true
	}
	return nil, false
}
