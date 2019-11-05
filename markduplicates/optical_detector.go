package markduplicates

import (
	"sort"
	"strings"

	"github.com/grailbio/base/log"
	"github.com/grailbio/bio/encoding/bam"
	"github.com/grailbio/bio/encoding/bampair"
)

type sortingEntry struct {
	library      string
	leftRefId    int
	left5Pos     int
	orientation  Orientation
	rightRefId   int
	right5Pos    int
	leftFileIdx  uint64
	rightFileIdx uint64
	pair         IndexedPair

	// These are not used for sorting.
	location  PhysicalLocation
	duplicate bool
}
type sortingTable []sortingEntry

func (t sortingTable) Swap(i, j int) {
	t[i], t[j] = t[j], t[i]
}
func (t sortingTable) Len() int {
	return len(t)
}

// Use the same sort order as picard.
func (t sortingTable) Less(i, j int) bool {
	diff := strings.Compare(t[i].library, t[j].library)
	if diff == 0 {
		diff = t[i].leftRefId - t[j].leftRefId
	}
	if diff == 0 {
		diff = t[i].left5Pos - t[j].left5Pos
	}
	if diff == 0 {
		diff = int(t[i].orientation - t[j].orientation)
	}
	if diff == 0 {
		diff = t[i].rightRefId - t[j].rightRefId
	}
	if diff == 0 {
		diff = t[i].right5Pos - t[j].right5Pos
	}
	if diff == 0 {
		if t[i].leftFileIdx > t[j].leftFileIdx {
			diff = 1
		} else if t[i].leftFileIdx < t[j].leftFileIdx {
			diff = -1
		} else {
			diff = 0
		}
	}
	if diff == 0 {
		if t[i].rightFileIdx > t[j].rightFileIdx {
			diff = 1
		} else if t[i].rightFileIdx < t[j].rightFileIdx {
			diff = -1
		} else {
			diff = 0
		}
	}
	return diff < 0
}

// TileOpticalDetector detects optical duplicates with a tile. For two
// reads to be optical duplicates, their tile, lane, surface, library,
// and read orientations must be identical
type TileOpticalDetector struct {
	OpticalDistance int
}

// GetRecordProcessor implements OpticalDetector.
func (t *TileOpticalDetector) GetRecordProcessor() bampair.RecordProcessor {
	return nil
}

// RecordProcessorsDone implements OpticalDetector.
func (t *TileOpticalDetector) RecordProcessorsDone() {
}

// Detect implements OpticalDetector.
func (t *TileOpticalDetector) Detect(readGroupLibrary map[string]string, duplicates []DuplicateEntry, bestIndex int) []string {
	// Split duplicates by tile number into batches before marking the
	// optical duplicates.  We split by tile to reduce the cost of
	// comparing each pair against the other pairs.
	type batchKey struct {
		lane            int
		tile            int
		readGroup       string
		readGroupFound  bool
		r1R2Orientation Orientation
	}

	batches := make(map[batchKey]sortingTable)
	var bestBatchKey batchKey
	bestName := ""
	duplicateNames := make([]string, 0)
	for i, pair := range duplicates {
		p := pair.(IndexedPair)
		location := ParseLocation(pair.Name())
		readGroup, readGroupFound := getReadGroup(p.Left.R)
		key := batchKey{
			lane:            location.Lane,
			tile:            location.TileName,
			readGroup:       readGroup,
			readGroupFound:  readGroupFound,
			r1R2Orientation: GetR1R2Orientation(&p),
		}

		if i == bestIndex {
			bestBatchKey = key
			bestName = pair.Name()
		}

		if _, found := batches[key]; !found {
			batches[key] = make([]sortingEntry, 0)
		}
		batches[key] = append(batches[key],
			sortingEntry{
				library:      GetLibrary(readGroupLibrary, p.Left.R),
				leftRefId:    p.Left.R.Ref.ID(),
				left5Pos:     bam.UnclippedFivePrimePosition(p.Left.R),
				orientation:  orientationBytePair(bam.IsReversedRead(p.Left.R), bam.IsReversedRead(p.Right.R)),
				rightRefId:   p.Right.R.Ref.ID(),
				right5Pos:    bam.UnclippedFivePrimePosition(p.Right.R),
				leftFileIdx:  p.Left.FileIdx_,
				rightFileIdx: p.Right.FileIdx_,
				pair:         p,
				location:     location,
				duplicate:    false,
			})
	}

	// Mark optical duplicates for each tile at a time.
	for key, batch := range batches {
		if log.At(log.Debug) && len(batch) > 1 {
			log.Debug.Printf("optical batch size: %d, %v", len(batch), key)
		}
		sort.Sort(batch)
		bestIdx := -1
		foundOptical := false
		if key == bestBatchKey {
			// If this batch contains the primary pair, then compare
			// all pairs against the primary first.
			for i := range batch {
				if batch[i].pair.Left.R.Name == bestName {
					bestIdx = i
					break
				}
			}
			for i := range batch {
				if bestIdx == i {
					continue
				}
				if isOpticalDup(t.OpticalDistance, &batch[bestIdx].location, &batch[i].location) {
					foundOptical = true
					batch[i].duplicate = true
					duplicateNames = append(duplicateNames, batch[i].pair.Left.R.Name)
					if log.At(log.Debug) {
						log.Debug.Printf("optical dups: %s %s (dup)", batch[bestIdx].pair.Left.R.Name,
							batch[i].pair.Left.R.Name)
					}
				}
			}
		}

		// Next, compare each pair with each other pair.
		for i := 0; i < len(batch); i++ {
			if i == bestIdx {
				continue
			}
			for j := i + 1; j < len(batch); j++ {
				if j == bestIdx {
					continue
				}
				if batch[i].duplicate && batch[j].duplicate {
					continue
				}
				if isOpticalDup(t.OpticalDistance, &batch[i].location, &batch[j].location) {
					if batch[j].duplicate {
						foundOptical = true
						batch[i].duplicate = true
						duplicateNames = append(duplicateNames, batch[i].pair.Left.R.Name)
						if log.At(log.Debug) {
							log.Debug.Printf("optical dups: %s %s (dup)", batch[j].pair.Left.R.Name,
								batch[i].pair.Left.R.Name)
						}
					} else {
						foundOptical = true
						batch[j].duplicate = true
						duplicateNames = append(duplicateNames, batch[j].pair.Left.R.Name)
						if log.At(log.Debug) {
							log.Debug.Printf("optical dups: %s %s (dup)", batch[i].pair.Left.R.Name,
								batch[j].pair.Left.R.Name)
						}
					}
				}
			}
		}
		if log.At(log.Debug) && foundOptical {
			log.Debug.Printf("duplicate group:")
			for i, e := range batch {
				log.Debug.Printf("  names[%d] %s optical dup: %v, best: %v, entry: %v",
					i, e.pair.Left.R.Name, e.duplicate, i == bestIdx, e)
			}
		}
	}
	return duplicateNames
}

func isOpticalDup(opticalDistance int, a, b *PhysicalLocation) bool {
	return abs(a.X-b.X) <= opticalDistance && abs(a.Y-b.Y) <= opticalDistance
}
