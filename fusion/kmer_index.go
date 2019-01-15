package fusion

import (
	"math"
	"sort"
	"unsafe"

	farm "github.com/dgryski/go-farm"
	"github.com/grailbio/base/log"
	"golang.org/x/sys/unix"
)

// This file implements a singleton kmer -> genelist map.  This map is
// physically sharded 256-ways, using upper 8 bits of farmhash(kmer) to pick the
// shard. Within one shard, we use the lower bits of the same farmhash(kmer) to
// implement a vanilla linear-probing hashtable.

const (
	nKmerIndexShard    = 256 // # of shards in the hash table.
	maxCollisions      = 64  // max# of collisions allowed in a hash table for a single lookup.
	kmerIndexEntrySize = unsafe.Sizeof(kmerIndexEntry{})
)

// kmerIndex is the central kmer->genelist map.  It is logically equivalent to
// map[Kmer][]geneID, but it is hand-coded to reduce memory and
// computational overhead.
//
// kmerIndex is sharded 4096 ways.  The lower 12 bits of farmhash(kmer) is used
// to select the shard. The upper bits of the same hash is used to pick a bucket
// in the shard.
type kmerIndex [nKmerIndexShard]kmerIndexShard

// One shard of kmerIndex. It is a vanilla linear-probing hash table.
type kmerIndexShard struct {
	nShift uint32 // == ceil(log2(#-of-kmers))

	// The hash table is logically [n]kmerIndexEntry, where n=#-of-kmers, but it
	// is created in an anon-mapped memory region, with madvise(MADV_HUGEPAGE) to
	// reduce TLB misses.  tableStart and tableLimit are the start and the limit
	// of the (logical) [n]kmerIndexEntry.
	tableStart unsafe.Pointer
	tableLimit unsafe.Pointer

	// outlined is an arrray of geneID. It is used to store genelist when one
	// kmer has >2 entries. See the kmerIndexEntry doc for more details.
	outlined unsafe.Pointer
}

// kmerIndexEntry is the entry for the central kmer->genelist map.  If an entry
// contains <= 2 genepos, they are stored in the "inlined" field. Else, genepos
// are stored in slice kmerIndexShard.outlined[-inlined[0]:-inlined[1]].
// One can tell whether the genepos are inlined or outlined by checking
// inlined[0] > 0.
//
// This layout reduces the memory consumption and avoids cache misses for the
// common case where there are <= 2 genes for the kmer. In addition, this layout
// contains only two pointers per index shard (kmerIndexShard.{table,outlined}),
// so it minimizes the GC scan overhead.
//
// Genes in kmerIndexEntry are sorted in ascending order of geneID, then pos.
type kmerIndexEntry struct {
	kmer    Kmer
	inlined [2]GeneID
}

func hashKmer(k Kmer) uint64 {
	// TODO(saito) consider switching to xxhash. Should give about 6% speedup.
	return farm.Hash64WithSeed(nil, uint64(k))
}

// initShard fills one shard of of kmerIndex. maxGenesPerKmer should be a copy of
// Opts.MaxGenesPerKmer.The caller must guarantee that all the kmers in the
// input are in the given shard. Thread compatible.
func (idx *kmerIndex) initShard(shard int, input map[Kmer]*[]GeneID, maxGenesPerKmer int) {
	const (
		hugePageSize = 2 << 20 // size of Linux transparent hugetlb.
		loadFactor   = 4       // hashtable load factor
	)
	minSize := int((float64(len(input) + 1)) * loadFactor)
	// Compute shift = ceil(log2(minSize)), size = 2^shift
	size := 1
	shift := 0
	for size < minSize {
		if size*2 < size {
			panic(minSize)
		}
		size *= 2
		shift++
	}

	// Use the upper "shift" bits of the hash to select a bucket.
	sizeShift := 64 - shift

	// Set up transparent hugepages.  Ubuntu, by default, activates THPs only for
	// mavdised regions, so we bypass Go's standard memory allocator.
	//
	// For more details, see:
	// https://www.kernel.org/doc/Documentation/vm/transhuge.txt.
	tableData, err := unix.Mmap(-1, 0, size*int(kmerIndexEntrySize)+hugePageSize,
		unix.PROT_READ|unix.PROT_WRITE, unix.MAP_ANON|unix.MAP_PRIVATE)
	if err != nil {
		log.Panic(err)
	}
	if err := unix.Madvise(tableData, unix.MADV_HUGEPAGE); err != nil {
		log.Panic(err)
	}
	// Round the page up to a hugePageSize boundary. It's not clear if this helps,
	// but at worst, it is a noop.
	tableStart := ((uintptr(unsafe.Pointer(&tableData[0]))-1)/hugePageSize + 1) * hugePageSize
	tableLimit := (tableStart + uintptr(size)*kmerIndexEntrySize)

	// At this point, memory range [tableStart,tableLimit) covers the hashtable.
	// First initialize the hash table.
	for i := 0; i < size; i++ {
		ent := (*kmerIndexEntry)(unsafe.Pointer(tableStart + kmerIndexEntrySize*uintptr(i)))
		ent.kmer = invalidKmer
	}

	// Add the entries into the hash table.
	var (
		outlined []GeneID // accumulates outlined index entries.
	)
	for kmer, genesPtr := range input {
		genes := *genesPtr
		// Sort then dedup the genes
		sort.SliceStable(genes, func(i, j int) bool {
			return genes[i] < genes[j]
		})
		nGene := 1
		for i := 1; i < len(genes); i++ {
			if genes[nGene-1] != genes[i] {
				nGene++
				genes[nGene-1] = genes[i]
			}
		}
		genes = genes[:nGene]
		if len(genes) > maxGenesPerKmer {
			continue
		}
		h := hashKmer(kmer)
		if h&(nKmerIndexShard-1) != uint64(shard) {
			panic(kmer)
		}
		entPtr := tableStart + kmerIndexEntrySize*uintptr(h>>uint(sizeShift))
		var ent *kmerIndexEntry
		// Linear-probe to find the place for this kmer.
		for iter := 0; ; iter++ {
			ent = (*kmerIndexEntry)(unsafe.Pointer(entPtr))
			if ent.kmer == invalidKmer {
				break
			}
			if iter > maxCollisions {
				log.Panicf("size %d, shift %d", size, shift)
			}
			entPtr += kmerIndexEntrySize
			if entPtr >= tableLimit {
				entPtr = tableStart
			}
		}

		ent.kmer = kmer
		switch len(genes) {
		case 1:
			ent.inlined[0] = genes[0]
		case 2:
			ent.inlined[0] = genes[0]
			ent.inlined[1] = genes[1]
		default:
			if len(genes) > math.MaxInt32 {
				panic("too many genes at a kmer")
			}
			ent.inlined[0] = -GeneID(len(outlined))
			ent.inlined[1] = -GeneID(len(outlined) + len(genes))
			outlined = append(outlined, genes...)
		}
	}
	var outlinedPtr unsafe.Pointer
	if len(outlined) > 0 {
		outlinedPtr = unsafe.Pointer(&outlined[0])
	}

	idx[shard] = kmerIndexShard{
		nShift:     uint32(sizeShift),
		tableStart: unsafe.Pointer(tableStart),
		tableLimit: unsafe.Pointer(tableLimit),
		outlined:   outlinedPtr,
	}
}

// kmerIndexIterator lists geneIDs found by kmerIndex.get().
type kmerIndexIterator struct {
	ent      *kmerIndexEntry // nil if kmer is not found.
	outlined unsafe.Pointer  // copy of kmerIndexShard.outlined.
}

// Get finds an entry for the given kmer.  Thread safe.
//
// Example:
//   iter := idx.get(kmer)
//   for i := 0; ; i++ {
//      pos, ok := iter.get(i)
//      if !ok { break }
//      .. use pos ..
//   }
func (idx *kmerIndex) get(kmer Kmer) kmerIndexIterator {
	h := hashKmer(kmer)
	shard := (*idx)[h&(nKmerIndexShard-1)]

	tableStart := uintptr(shard.tableStart)
	tableLimit := uintptr(shard.tableLimit)
	entPtr := tableStart + kmerIndexEntrySize*uintptr(h>>shard.nShift)
	for iter := 0; iter <= maxCollisions; iter++ {
		ent := (*kmerIndexEntry)(unsafe.Pointer(entPtr))
		if ent.kmer == kmer {
			return kmerIndexIterator{ent: ent, outlined: shard.outlined}
		}
		if ent.kmer == invalidKmer {
			return kmerIndexIterator{}
		}
		entPtr += kmerIndexEntrySize
		if entPtr >= tableLimit {
			entPtr = tableStart
		}
	}
	return kmerIndexIterator{}
}

// Get returns the i'th geneID. It returns invalidGeneID if i >= number of genes
// found for the kmer.
func (iter kmerIndexIterator) get(i int) GeneID {
	if iter.ent == nil {
		return invalidGeneID
	}
	if iter.ent.inlined[0] > 0 { // inlined repl
		if i == 0 {
			return iter.ent.inlined[0]
		}
		if i == 1 && iter.ent.inlined[1] > 0 {
			return iter.ent.inlined[1]
		}
		return invalidGeneID
	}
	// outlined.
	start := -int(iter.ent.inlined[0])
	limit := -int(iter.ent.inlined[1])
	if start >= limit-2 {
		panic(*iter.ent)
	}
	if start+i >= limit {
		return invalidGeneID
	}
	const geneIDSize = 4 // sizeof(GeneID)
	p := unsafe.Pointer((uintptr)(iter.outlined) + geneIDSize*uintptr(start+i))
	return *(*GeneID)(p)
}
