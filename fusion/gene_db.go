package fusion

import (
	"context"
	"fmt"
	"io"
	"io/ioutil"
	"os"
	"regexp"
	"runtime"
	"strconv"
	"strings"
	"sync"

	"github.com/grailbio/base/file"
	"github.com/grailbio/base/log"
	"github.com/grailbio/base/tsv"
	"github.com/grailbio/bio/encoding/fasta"
)

// GeneID is a dense sequence number (1, 2, 3, ...) assigned to a gene (e.g.,
// "MAPK10").  IDs are valid only within one process invocation.
//
// Caution: this type must be signed. kmer_index uses negative geneids to indicate
// outlined slices.
type GeneID int32

const invalidGeneID = GeneID(0)

type fusionEventPair struct{ gene1, gene2 GeneID }

// GeneInfo stores the info about a gene. It is parsed out from the
// transcriptome fasta key.
//
// Transcriptome key example: "ENST00000279783.3|OR8K1|chr11:56346039-56346998:1051|960"
type GeneInfo struct {
	// ID is a dense sequence number (1, 2, ...). It is valid only during the current run.
	ID GeneID
	// EnsemblID is parsed from the transcriptome FASTA key. E.g., "ENST00000279783.3"
	EnsemblID string
	// Gene is parsed from the transcriptome FASTA key. E.g., "OR8K1"
	Gene string
	// Chrom is parsed from the transcriptome FASTA key. E.g., "chr11"
	Chrom string
	// Start is parsed from the transcriptome FASTA key. E.g., 56346039
	Start int
	// End is parsed from the transcriptome FASTA key. E.g., 56346998
	End int
	// Index is the rank of this gene in Chrom. Gene with the smallest Start in
	// the given Chrom will have Index of zero.
	Index int
	// FusionEvent is true if this gene appears in the cosmic TSV file added via
	// ReadFusionEvents.
	FusionEvent bool
}

var transcriptomeRefRE = regexp.MustCompile(`^([^|]+)\|([^|]+)\|([^:]+):(\d+)-(\d+):(\d+)`)

// ParseTranscriptomeKey parses a transcriptome fasta key.
//
// Transcriptome key example: "ENST00000279783.3|OR8K1|chr11:56346039-56346998:1051|960"
func ParseTranscriptomeKey(seqName string) (ensemblID, gene, chrom string, start, end, index int, err error) {
	// Key encoding:
	// ENST00000279783.3|OR8K1|chr11:56346039-56346998:1051|960
	//
	// TODO(saito) make this faster.
	m := transcriptomeRefRE.FindStringSubmatch(seqName)
	if m == nil {
		err = fmt.Errorf("ParseTranscriptome: failed to parse '%s'", seqName)
		return
	}
	ensemblID = m[1]
	gene = m[2]
	chrom = m[3]

	parseInt := func(s string) int {
		v, err := strconv.ParseInt(s, 0, 32)
		if err != nil {
			log.Panicf("ParseTranscrptome %s: %v", seqName, err)
		}
		return int(v)
	}
	// TODO(saito) is this base zero or one?
	start = parseInt(m[4])
	end = parseInt(m[5])
	index = parseInt(m[6])
	return
}

// GeneDB is a singleton object that stores transcriptomes, kmers generated from
// the transcripts, and candidate fusion event pairs. Thread compatible.
type GeneDB struct {
	opts Opts

	// Names maps gene names to dense IDs.  Names contains all the genes listed in
	// fusion events file and transcriptome.
	names           map[string]GeneID
	namesFrozen     bool
	hasFusionEvents bool

	genes []*GeneInfo // indexed by geneID
	pairs map[fusionEventPair]struct{}

	// kmerIndex maps kmer -> list of genes that have the kmer
	kmerIndex kmerIndex
}

// GeneIDRange returns the range of gene IDs registered in this object.  The low
// end is closed, high end is open. For example, the return value of (1, 95)
// means this DB holds 94 genes, IDs from 1 to 94. You can use GeneInfo() to get
// the information about the gene.
func (m *GeneDB) GeneIDRange() (GeneID, GeneID) { return 1, GeneID(len(m.genes)) }

// GeneInfo gets the GeneInfo given an ID. It always returns a non-nil info.
//
// REQUIRES: ID is valid.
func (m *GeneDB) GeneInfo(id GeneID) *GeneInfo {
	if id == invalidGeneID {
		panic(id)
	}
	return m.genes[id]
}

// GeneInfoByName gets GeneINfo given a gene name. It returns nil if the gene is
// not registered.
func (m *GeneDB) GeneInfoByName(name string) *GeneInfo {
	id := m.names[name]
	if id == invalidGeneID {
		return nil
	}
	return m.genes[id]
}

// FindByKmer finds info about the given kmer.
func (m *GeneDB) findByKmer(k Kmer) kmerIndexIterator {
	return m.kmerIndex.get(k)
}

// IsFusionPair checks if the given pair of genes appear in the cosmic TSV file
// added by ReadFusionEvents.
//
// REQUIRES: IDs are valid.
func (m *GeneDB) IsFusionPair(g1, g2 GeneID) bool {
	_, ok := m.pairs[fusionEventPair{g1, g2}]
	return ok
}

// PrepopulateGenes x assigns gene IDs to genes. NOT FOR GENERAL USE. It is used only to change
// the genename <-> geneID assignments to reproduce the behavior of the C++
// code.
func (m *GeneDB) PrepopulateGenes(names []string) {
	for _, name := range names {
		m.internGene(name)
	}
	m.namesFrozen = true
}

// PrepopulateGeneInfo x fills geneinfo in batch. NOT FOR GENERAL USE. It is
// used to populate gene DB from recordio dump.
func (m *GeneDB) PrepopulateGeneInfo(genes []GeneInfo) {
	for _, info := range genes {
		gid := m.internGene(info.Gene)
		if gid != info.ID {
			panic(info)
		}
		gi := m.GeneInfo(gid)
		*gi = info
	}
	m.namesFrozen = true
}

// InternGene finds or assigns an ID to the gene with the given name.
func (m *GeneDB) internGene(name string) GeneID {
	if id, ok := m.names[name]; ok {
		return id
	}
	if m.namesFrozen {
		panic(name)
	}
	id := GeneID(len(m.genes))
	m.names[name] = id
	m.genes = append(m.genes, &GeneInfo{ID: id, Gene: name})
	return id
}

func (m *GeneDB) registerFusionPair(id0, id1 GeneID) {
	m.genes[id0].FusionEvent = true
	m.genes[id1].FusionEvent = true
	m.pairs[fusionEventPair{id0, id1}] = struct{}{}
}

// geneID retrieves the gene ID given a name. It returns invalidGeneID if not
// found.
func (m *GeneDB) geneID(name string) GeneID {
	return m.names[name]
}

// NewGeneDB creates an empty GeneDB.
func NewGeneDB(opts Opts) *GeneDB {
	return &GeneDB{
		opts:  opts,
		names: map[string]GeneID{},
		genes: []*GeneInfo{&GeneInfo{Gene: "invalid"}},
		pairs: map[fusionEventPair]struct{}{},
	}
}

// ReadFusionEvents reads from a Cosmic TSV file the names of gene pairs that
// form fusions. The first column of each line must be of form "gene1/gene2",
// for example "ACSL3/ETV1".
func (m *GeneDB) ReadFusionEvents(ctx context.Context, path string) {
	m.hasFusionEvents = true
	in, err := file.Open(ctx, path)
	if err != nil {
		log.Panicf("open %s: %v", path, err)
	}
	r := tsv.NewReader(in.Reader(ctx))
	r.HasHeaderRow = true
	r.ValidateHeader = true

	row := struct{ Genes string }{} // Rest of the fields are ignored
	nLine := 0
	for {
		if err := r.Read(&row); err != nil {
			if err == io.EOF {
				break
			}
			log.Panic(err)
		}
		genes := strings.Split(row.Genes, "/")
		if len(genes) != 2 {
			log.Panicf("read tsv %s:%d: expect 'gene1/gene2', but found %+v", path, nLine, row)
		}
		m.registerFusionPair(m.internGene(genes[0]), m.internGene(genes[1]))
		nLine++
	}
	if err := in.Close(ctx); err != nil {
		log.Panicf("close %s: %v", path, err)
	}
}

// ReadTranscriptome reads a transcriptome reference fasta file.
func (m *GeneDB) ReadTranscriptome(ctx context.Context, fastaPath string, filter bool) {
	if filter != m.hasFusionEvents {
		panic("filter")
	}
	generateIndex := func() (string, func()) {
		index, err := ioutil.TempFile("", "")
		if err != nil {
			log.Panicf("tempfile: %v", err)
		}

		in, err := file.Open(ctx, fastaPath)
		if err != nil {
			log.Panicf("generateIndex %s: %v", fastaPath, err)
		}
		if err = fasta.GenerateIndex(index, in.Reader(ctx)); err != nil {
			log.Panicf("generateIndex %s: %v", fastaPath, err)
		}
		if err = in.Close(ctx); err != nil {
			log.Panicf("generateIndex close %s: %v", fastaPath, err)
		}
		if err = index.Close(); err != nil {
			log.Panicf("generateIndex close %s: %v", index.Name(), err)
		}
		indexPath := index.Name()
		return indexPath, func() {
			if err := os.Remove(indexPath); err != nil {
				log.Panicf("remove %s: %v", indexPath, err)
			}
		}
	}

	type fa struct {
		in, idxIn file.File
		fa        fasta.Fasta
	}

	// TODO(saito) Use a preexisting index if provided.
	indexPath, cleanup := generateIndex()
	defer cleanup()

	openFASTA := func() *fa {
		fa := fa{}
		var err error
		if fa.in, err = file.Open(ctx, fastaPath); err != nil {
			log.Panicf("open %s: %v", fastaPath, err)
		}
		if fa.idxIn, err = file.Open(ctx, indexPath); err != nil {
			log.Panicf("open %s: %v", indexPath, err)
		}
		if fa.fa, err = fasta.NewIndexed(fa.in.Reader(ctx), fa.idxIn.Reader(ctx)); err != nil {
			log.Panicf("fasta.NewIndexed %s,%s: %v", fastaPath, indexPath, err)
		}
		return &fa
	}

	closeFASTA := func(fa *fa) {
		if err := fa.in.Close(ctx); err != nil {
			log.Panicf("close %s: %v", fastaPath, err)
		}
		if err := fa.idxIn.Close(ctx); err != nil {
			log.Panicf("close %s: %v", fastaPath, err)
		}
	}

	if filter {
		log.Printf("Reading transcriptome %s filtered", fastaPath)
	} else {
		log.Printf("Reading transcriptome %s denovo", fastaPath)
	}
	type kmerShard struct {
		mu    sync.Mutex
		kmers map[Kmer]*[]GeneID
	}
	kmerShards := [nKmerIndexShard]kmerShard{}

	registerKmer := func(geneID GeneID, kmPos kmersAtPos) {
		k := kmPos.minKmer()
		shard := &kmerShards[hashKmer(k)%nKmerIndexShard]
		shard.mu.Lock()
		if shard.kmers == nil {
			shard.kmers = map[Kmer]*[]GeneID{}
		}
		ent, ok := shard.kmers[k]
		if !ok {
			ent = &[]GeneID{}
			shard.kmers[k] = ent
		}
		*ent = append(*ent, geneID)
		shard.mu.Unlock()
	}

	// Produce kmers for genes in parallel.
	var nAdded, nSkipped int
	{
		type req struct {
			geneID  GeneID
			seqName string
		}
		reqCh := make(chan req, 1024)
		wg := sync.WaitGroup{}
		for i := 0; i < runtime.NumCPU(); i++ {
			wg.Add(1)
			go func() {
				km := newKmerizer(m.opts.KmerLength)
				fa := openFASTA()
				defer closeFASTA(fa)
				for req := range reqCh {
					seqLen, err := fa.fa.Len(req.seqName)
					if err != nil {
						panic(err)
					}
					seq, err := fa.fa.Get(req.seqName, 0, seqLen)
					if err != nil {
						panic(err)
					}
					km.Reset(seq)
					for km.Scan() {
						registerKmer(req.geneID, km.Get())
					}
				}
				wg.Done()
			}()
		}

		fa := openFASTA()
		for _, seqName := range fa.fa.SeqNames() {
			ensemblID, gene, chrom, start, end, index, err := ParseTranscriptomeKey(seqName)
			if err != nil {
				log.Panic(err)
			}
			var gi *GeneInfo
			if filter {
				gi = m.GeneInfoByName(gene)
				if gi == nil || !gi.FusionEvent {
					nSkipped++
					continue
				}
			} else {
				gi = m.GeneInfo(m.internGene(gene))
				gi.FusionEvent = true
			}
			// Update the detailed info. In case the fasta file has multiple entries
			// for one gene, they should have the same GeneInfo.
			//
			// TODO(saito) This may overwrite the ensembl ID.
			gi.EnsemblID, gi.Gene, gi.Chrom, gi.Start, gi.End, gi.Index = ensemblID, gene, chrom, start, end, index
			nAdded++
			reqCh <- req{gi.ID, seqName}
		}
		closeFASTA(fa)
		close(reqCh)
		wg.Wait()
	}

	log.Printf("Finished reading transcriptome %s, added %d genes, skipped %d genes", fastaPath, nAdded, nSkipped)

	// Post-process each kmer - dedup the list of genes and sort them.
	{
		// List of raw info for kmers
		type req struct {
			shardIndex int
			kmers      map[Kmer]*[]GeneID
		}
		reqCh := make(chan req, nKmerIndexShard)
		producerWg := sync.WaitGroup{}
		for i := 0; i < runtime.NumCPU(); i++ {
			producerWg.Add(1)
			go func() {
				for req := range reqCh {
					m.kmerIndex.initShard(req.shardIndex, req.kmers, m.opts.MaxGenesPerKmer)
				}
				producerWg.Done()
			}()
		}

		for si := range kmerShards {
			reqCh <- req{shardIndex: si, kmers: kmerShards[si].kmers}
			kmerShards[si].kmers = nil
		}
		close(reqCh)
		producerWg.Wait()
	}
	log.Printf("Finished postprocessing transcriptome %s", fastaPath)
}
