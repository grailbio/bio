package markduplicates

import (
	"context"
	"fmt"
	"os"
	"sync"

	"github.com/grailbio/base/errors"
	"github.com/grailbio/base/log"
)

// Metrics contains metrics from mark duplicates.
type Metrics struct {
	// Implement the metrics reported by picard

	// UnpairedReads is the number of mapped reads examined which did
	// not have a mapped mate pair, either because the read is
	// unpaired, or the read is paired to an unmapped mate.
	UnpairedReads int

	// ReadPairsExamined is the number of mapped read pairs
	// examined. (Primary, non-supplemental).
	ReadPairsExamined int

	// SecondarySupplementary is the number of reads that were either
	// secondary or supplementary.
	SecondarySupplementary int

	// UnmappedReads is the total number of unmapped reads
	// examined. (Primary, non-supplemental).
	UnmappedReads int

	// UnpairedDups is the number of fragments that were marked as duplicates.
	UnpairedDups int

	// ReadPairDups is the number of read pairs that were marked as duplicates.
	ReadPairDups int

	// ReadPairOpticalDups is the number of read pairs duplicates that
	// were caused by optical duplication. Value is always <
	// READ_PAIR_DUPLICATES, which counts all duplicates regardless of
	// source.
	ReadPairOpticalDups int
}

// String returns a string representation of the metrics contained in
// m. The string can be used as metrics file output.
func (m *Metrics) String() string {
	librarySizeStr := "0"
	a := uint64((m.ReadPairsExamined / 2) - (m.ReadPairOpticalDups / 2))
	b := uint64((m.ReadPairsExamined / 2) - (m.ReadPairDups / 2))
	librarySize, err := estimateLibrarySize(a, b)
	if err == nil {
		librarySizeStr = fmt.Sprintf("%v", librarySize)
	} else {
		log.Error.Printf("error in estimateLibrarySize(%v, %v): %v, ", a, b, err)
	}

	return fmt.Sprintf("%d\t%d\t%d\t%d\t%d\t%d\t%d\t%0.6f\t%v", m.UnpairedReads, m.ReadPairsExamined/2,
		m.SecondarySupplementary, m.UnmappedReads, m.UnpairedDups,
		m.ReadPairDups/2, m.ReadPairOpticalDups/2,
		100*(float64(m.UnpairedDups+m.ReadPairDups)/float64(m.UnpairedReads+m.ReadPairsExamined)),
		librarySizeStr)
}

// Add adds the metrics in other to m.
func (m *Metrics) Add(other *Metrics) {
	m.UnpairedReads += other.UnpairedReads
	m.ReadPairsExamined += other.ReadPairsExamined
	m.SecondarySupplementary += other.SecondarySupplementary
	m.UnmappedReads += other.UnmappedReads
	m.UnpairedDups += other.UnpairedDups
	m.ReadPairDups += other.ReadPairDups
	m.ReadPairOpticalDups += other.ReadPairOpticalDups
}

// MetricsCollection contains metrics computed by Mark.
type MetricsCollection struct {
	// Global metrics
	maxAlignDist int

	// OpticalDistance stores the number of duplicate read pairs that
	// have the given Euclidean distance.
	OpticalDistance [][]int64

	// LibraryMetrics contains per-library metrics.
	LibraryMetrics map[string]*Metrics

	mutex sync.Mutex
}

func newMetricsCollection() *MetricsCollection {
	mc := &MetricsCollection{
		LibraryMetrics:  make(map[string]*Metrics),
		OpticalDistance: make([][]int64, 4),
	}
	for i := range mc.OpticalDistance {
		mc.OpticalDistance[i] = make([]int64, 60000)
	}
	return mc
}

// Get returns Metrics for the given library. If there is no Metrics
// for library yet, create one and return it.
func (mc *MetricsCollection) Get(library string) *Metrics {
	m, found := mc.LibraryMetrics[library]
	if found {
		return m
	}
	m = &Metrics{}
	mc.LibraryMetrics[library] = m
	return m
}

// Merge per-library and optical distance metrics from other
// into mc.
func (mc *MetricsCollection) Merge(other *MetricsCollection) {
	mc.mutex.Lock()
	defer mc.mutex.Unlock()

	for library, otherMetrics := range other.LibraryMetrics {
		existing, found := mc.LibraryMetrics[library]
		if found {
			existing.Add(otherMetrics)
		} else {
			// Make a copy to be owned by m.
			new := *otherMetrics
			mc.LibraryMetrics[library] = &new
		}
	}
	for i := range mc.OpticalDistance {
		if len(mc.OpticalDistance[i]) < len(other.OpticalDistance[i]) {
			temp := make([]int64, len(other.OpticalDistance[i]))
			copy(temp, mc.OpticalDistance[i])
			mc.OpticalDistance[i] = temp
		}
		for j := range other.OpticalDistance[i] {

			mc.OpticalDistance[i][j] += other.OpticalDistance[i][j]
		}
	}
}

// AddDistance increments the histogram counter for the given bagsize
// and distance.
func (mc *MetricsCollection) AddDistance(bagSize, distance int) {
	if distance >= len(mc.OpticalDistance[0]) {
		for i := range mc.OpticalDistance {
			temp := make([]int64, distance+1)
			copy(temp, mc.OpticalDistance[i])
			mc.OpticalDistance[i] = temp
		}
	}

	if bagSize <= 2 {
		mc.OpticalDistance[0][distance]++
	} else if bagSize >= 3 && bagSize <= 4 {
		mc.OpticalDistance[1][distance]++
	} else if bagSize >= 5 && bagSize <= 7 {
		mc.OpticalDistance[2][distance]++
	} else if bagSize >= 8 {
		mc.OpticalDistance[3][distance]++
	}
}

func writeMetrics(ctx context.Context, opts *Opts, globalMetrics *MetricsCollection) (err error) {
	var f *os.File
	f, err = os.Create(opts.MetricsFile)
	if err != nil {
		return errors.E(err, "Couldn't create metrics file:", opts.MetricsFile)
	}
	defer func() {
		if err2 := f.Close(); err == nil && err2 != nil {
			err = err2
		}
	}()

	s := "# bio-mark-duplicates\n" +
		"# maximum 5' alignment distance: " + fmt.Sprintf("%d", globalMetrics.maxAlignDist) + "\n" +
		"LIBRARY\tUNPAIRED_READS_EXAMINED\tREAD_PAIRS_EXAMINED\t" +
		"SECONDARY_OR_SUPPLEMENTARY_RDS\tUNMAPPED_READS\tUNPAIRED_READ_DUPLICATES\t" +
		"READ_PAIR_DUPLICATES\tREAD_PAIR_OPTICAL_DUPLICATES\tPERCENT_DUPLICATION\t" +
		"ESTIMATED_LIBRARY_SIZE\n"

	for library, metrics := range globalMetrics.LibraryMetrics {
		s += library + "\t" + metrics.String() + "\n"
	}
	if _, err = f.Write([]byte(s)); err != nil {
		return errors.E(err, "error writing to metrics file:", opts.MetricsFile)
	}
	return nil
}

func writeOpticalHistogram(ctx context.Context, opts *Opts, globalMetrics *MetricsCollection) (err error) {
	var f *os.File
	f, err = os.Create(opts.OpticalHistogram)
	if err != nil {
		return errors.E(err, "Couldn't create optical histogram file:", opts.OpticalHistogram)
	}
	defer func() {
		if err2 := f.Close(); err == nil && err2 != nil {
			err = err2
		}
	}()

	if _, err = fmt.Fprintf(f, "#bag_size_range\toptical_dist\tcount\n"); err != nil {
		return errors.E(err, "error writing to optical histogram file:", opts.OpticalHistogram)
	}
	for i, prefix := range []string{"bagsize-2", "bagsize3-4", "bagsize5-7", "bagsize8-"} {
		for dist, count := range globalMetrics.OpticalDistance[i] {
			if _, err = fmt.Fprintf(f, "%s\t%d\t%d\n", prefix, dist, count); err != nil {
				return errors.E(err, "error writing to optical histogram file:", opts.OpticalHistogram)
			}
		}
	}
	return nil
}
