package markduplicates

import (
	"math"
	"math/rand"
	"sort"
	"strconv"
	"strings"

	"github.com/grailbio/base/log"
)

// PhysicalLocation describes a read's physical location on the flow
// cell. Lane, Surface, Swatch, Section, and TileNumber together
// specify which flowcell tile the read was found in. TileName is the
// 4 or 5 digit representation of the tile, e.g. 1203 means surface 1,
// swath 2 and tile 3. 12304 means surface 1, swath 2, section 3, and
// tile 4. X and Y describe the X and Y coordinates of the well within
// the tile.
type PhysicalLocation struct {
	Lane       int
	Surface    int
	Swath      int
	Section    int
	TileNumber int
	TileName   int
	X          int
	Y          int
}

const (
	// Illumina read names come in 3 varieties: 5, 7, and 8 columns.
	// For 5 and 7 field read names, the last three fields are:
	// tileName, X and Y. For 8 field read names, the last four fields
	// are tileName, X, Y, and UMI. These constants help keep track of
	// which fields are what.

	// IlluminaReadName5Fields is the number of columns in a 5 field read name.
	IlluminaReadName5Fields = 5
	// IlluminaReadName5FieldsTileField is 0-based field number that
	// contains the tileName for 5 field read names.
	IlluminaReadName5FieldsTileField = 2

	// IlluminaReadName7Fields is the number of columns in a 5 field read name.
	IlluminaReadName7Fields = 7
	// IlluminaReadName7FieldsTileField is 0-based field number that
	// contains the tileName for 7 field read names.
	IlluminaReadName7FieldsTileField = 4

	// IlluminaReadName8Fields is the number of columns in a 5 field read name.
	IlluminaReadName8Fields = 8
	// IlluminaReadName8FieldsTileField is 0-based field number that
	// contains the tileName for 8 field read names.
	IlluminaReadName8FieldsTileField = 4
)

// addOpticalDistances adds the optical distances between readpairs in
// duplicates to metrics. If opts.OpticalHistogramMax is >= 0, then
// limit to the first opts.OpticalHistogramMax readpairs after sorting
// by fileidx.
func addOpticalDistances(opts *Opts, readGroupLibrary map[string]string,
	originalDuplicates []DuplicateEntry, metrics *MetricsCollection) {

	// First sort pairs by fileidx to ensure deterministic behavior.
	duplicates := make([]DuplicateEntry, len(originalDuplicates))
	copy(duplicates, originalDuplicates)
	sort.Slice(duplicates, func(i, j int) bool {
		return duplicates[i].FileIdx() < duplicates[j].FileIdx()
	})

	// If we are capping the number of duplicate readpairs in the
	// optical histogram, then shuffle the reads so that the histogram
	// has a random sampling of the flow cell positions.
	if opts.OpticalHistogramMax >= 0 {
		r := rand.New(rand.NewSource(int64(duplicates[0].FileIdx())))
		r.Shuffle(len(duplicates), func(i, j int) {
			duplicates[i], duplicates[j] = duplicates[j], duplicates[i]
		})
	}

	if len(opts.OpticalHistogram) > 0 {
		type key struct {
			lane           int
			readGroup      string
			readGroupFound bool
			orientation    Orientation
		}
		m := map[key][]PhysicalLocation{}
		for _, dup := range duplicates {
			pair := dup.(IndexedPair)
			location := ParseLocation(dup.Name())
			readGroup, readGroupFound := getReadGroup(pair.Left.R)
			orientation := GetR1R2Orientation(&pair)

			k := key{
				lane:           location.Lane,
				readGroup:      readGroup,
				readGroupFound: readGroupFound,
				orientation:    orientation,
			}
			m[k] = append(m[k], location)
		}
		for _, locations := range m {
			for i := 0; i < len(locations) &&
				(opts.OpticalHistogramMax < 0 || i < opts.OpticalHistogramMax); i++ {
				for j := i + 1; j < len(locations) &&
					(opts.OpticalHistogramMax < 0 || j < opts.OpticalHistogramMax); j++ {
					metrics.AddDistance(len(duplicates),
						opticalDistance(&locations[i], &locations[j]))
				}
			}
		}
	}
}

func opticalDistance(a, b *PhysicalLocation) int {
	return int(math.Sqrt(math.Pow(float64(a.X-b.X), 2.0) + math.Pow(float64(a.Y-b.Y), 2.0)))
}

// ParseLocation returns a physical location given an Illumina style
// read name. The read name must have 5, 7, or 8 fields separated by
// ':'. When there are 5 or 7 fields, the last three fields are
// tileName, X and Y.  When there are 8 fields, the last four fields
// are tileName, X, Y, and UMI.
//
// The tileName be formatted as a 4 or 5 digit Illumina tileName.
// For a description of 4 digit tile numbers, see Appendix B, section Tile Numbering in
//  http://support.illumina.com.cn/content/dam/illumina-support/documents/documentation/system_documentation/hiseqx/hiseq-x-system-guide-15050091-e.pdf
//
// For a description of 5 digit tile numbers, see Appendix C, section Tile Numbering in
//   https://support.illumina.com/content/dam/illumina-support/documents/documentation/system_documentation/nextseq/nextseq-550-system-guide-15069765-05.pdf
func ParseLocation(qname string) PhysicalLocation {
	fields := strings.Split(qname, ":")
	var tileIdx int
	switch len(fields) {
	case IlluminaReadName5Fields:
		tileIdx = IlluminaReadName5FieldsTileField
	case IlluminaReadName7Fields:
		tileIdx = IlluminaReadName7FieldsTileField
	case IlluminaReadName8Fields:
		tileIdx = IlluminaReadName8FieldsTileField
	default:
		log.Fatalf("Could not parse name: %s, expected 5, 7, or 8 fields separated by ':'", qname)
	}

	var (
		location PhysicalLocation
		err      error
	)
	location.Lane, err = strconv.Atoi(fields[tileIdx-1])
	if err != nil {
		log.Fatalf("Could not parse name: %s, could not convert lane to integer: %v",
			qname, err)
	}
	location.TileName, err = strconv.Atoi(fields[tileIdx])
	if err != nil {
		log.Fatalf("Could not parse name: %s, could not convert tile to integer: %v",
			qname, err)
	}
	location.X, err = strconv.Atoi(fields[tileIdx+1])
	if err != nil {
		log.Fatalf("Could not parse name: %s, could not convert x to integer: %v",
			qname, err)
	}
	location.Y, err = strconv.Atoi(fields[tileIdx+2])
	if err != nil {
		log.Fatalf("Could not parse name: %s, could not convert y to integer: %v",
			qname, err)
	}

	if location.TileName > 99999 {
		log.Fatalf("Could not parse name: %s, unexpected tile name %d, expected 4 or 5 digits",
			qname, location.TileName)
	} else if location.TileName > 9999 {
		location.Surface = location.TileName / 10000
		location.Swath = (location.TileName % 10000) / 1000
		location.Section = (location.TileName % 1000) / 100
		location.TileNumber = location.TileName % 100
	} else {
		location.Surface = location.TileName / 1000
		location.Swath = (location.TileName % 1000) / 100
		location.TileNumber = location.TileName % 100
	}
	return location
}
