package main

import (
	"flag"
	"fmt"
	"os"
	"runtime"
	"runtime/debug"
	"sync"
	"time"

	"github.com/grailbio/base/grail"
	gbam "github.com/grailbio/bio/encoding/bam"
	"github.com/grailbio/bio/encoding/bamprovider"
	"github.com/grailbio/hts/sam"
	"v.io/x/lib/vlog"
)

var (
	bam                 string
	index               string
	parallelism         int
	shardSize           int
	failOnError         = false
	verbose             = true
	gcPercent           int
	printProgressAfterN int
	pairIterator        bool
)

func printProgress(i int64) bool {
	return i > 0 && int(i)%printProgressAfterN == 0
}

func init() {
	flag.StringVar(&bam, "bam", "testdata/test.bam", "BAM filename")
	flag.StringVar(&index, "index", "", "BAM index filename")
	flag.IntVar(&parallelism, "parallelism", runtime.NumCPU(), "Number of parallel computations to run")
	flag.IntVar(&shardSize, "shard-size", 1000000, "shard size in bp")
	flag.IntVar(&printProgressAfterN, "print-progress", 1000000, "print progress after this many records")
	flag.IntVar(&gcPercent, "gc-percent", 100, "set debug.SetGCPercent")
}

func main() {
	shutdown := grail.Init()
	defer shutdown()
	debug.SetGCPercent(gcPercent)
	if gcPercent != 100 {
		fmt.Printf("GC Percent set to %v\n", gcPercent)
	}
	if err := countPairs(parallelism, shardSize); err != nil {
		fmt.Fprintf(os.Stderr, "countPairs failed: %v", err)
	}
}

// Per-thread counter
type counter struct {
	sameReference      int64
	differentReference int64
	nerrs              int64
	npairs             int64
}

type stats struct {
	start  time.Time
	counts []counter // one per thread
}

func handleRecord(threadID int, r gbam.Pair, s *counter) {
	if r.Err != nil {
		s.nerrs++
		if failOnError {
			vlog.Fatalf("error getting records: %v", r.Err)
		} else {
			if verbose {
				vlog.Infof("error getting records: %v", r.Err)
			}
		}
		return
	}
	s.npairs++
	if r.R1.RefID() == r.R2.RefID() {
		s.sameReference++
	} else {
		s.differentReference++
	}
	if printProgress(s.npairs) {
		fmt.Fprintf(os.Stderr, "Counter %d saw %d same ref, %d different ref so far\n",
			threadID, s.sameReference, s.differentReference)
	}
	sam.PutInFreePool(r.R1)
	sam.PutInFreePool(r.R2)
}

func countPairs(parallelism int, shardSize int) error {
	stats := stats{
		start:  time.Now(),
		counts: make([]counter, parallelism),
	}
	vlog.Error("Using pair iterators")
	provider := bamprovider.NewProvider(bam, bamprovider.ProviderOpts{Index: index})
	iters, err := bamprovider.NewPairIterators(provider, false)
	if err != nil {
		return err
	}
	var waitGroup sync.WaitGroup
	waitGroup.Add(len(iters))
	for i, iter := range iters {
		go func(id int, iter *bamprovider.PairIterator, c *counter) {
			for iter.Scan() {
				recordOrError := iter.Record()
				handleRecord(id, recordOrError, c)
			}
			waitGroup.Done()
		}(i, iter, &stats.counts[i])
	}
	waitGroup.Wait()
	if err := provider.Close(); err != nil {
		return err
	}
	if verbose {
		var total counter
		fmt.Printf("Took: %s\n", time.Now().Sub(stats.start))
		for i, count := range stats.counts {
			fmt.Printf("Counter %d saw %d total pairs\n", i, count.differentReference+count.sameReference)
			total.differentReference += count.differentReference
			total.sameReference += count.sameReference
			total.npairs += count.npairs
			total.nerrs += count.nerrs
		}
		fmt.Printf("%d pairs with same reference, %d pairs with different reference\n",
			total.sameReference, total.differentReference)
		fmt.Printf("%v pairs, %v errors\n", total.npairs, total.nerrs)
	}
	return nil
}
