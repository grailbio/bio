package markduplicates

import (
	"compress/gzip"
	"context"
	"fmt"
	"io"
	"io/ioutil"
	"os"
	"strconv"
	"sync"
	"time"

	"github.com/grailbio/base/errors"
	"github.com/grailbio/base/file"
	"github.com/grailbio/base/log"
	"github.com/grailbio/base/vcontext"
	"github.com/grailbio/bio/biopb"
	"github.com/grailbio/bio/encoding/bam"
	"github.com/grailbio/bio/encoding/bampair"
	"github.com/grailbio/bio/encoding/bamprovider"
	"github.com/grailbio/bio/encoding/pam"
	"github.com/grailbio/bio/umi"
	"github.com/grailbio/hts/sam"
)

// BagProcessor takes the set of bags from a particular shard, and
// returns the same set of reads, but the reads may now be bagged
// differently. The new set of bags may contain more bags or fewer
// bags than the original set of bags.
type BagProcessor func([]*IntermediateDuplicateSet) []*IntermediateDuplicateSet

// BagProcessorFactory creates BagProcessors. Mark-duplciates creates
// one BagProcessor per goroutine.
type BagProcessorFactory interface {
	Create() BagProcessor
}

// OpticalDetector is a general interface for optical duplicate detection.
type OpticalDetector interface {
	// GetRecordProcessor returns a RecordProcessor that sees every
	// record in the bam input before any calls to Detect happen. The
	// OpticalDetector can use this to calculate statistics that
	// influence optical detection.
	GetRecordProcessor() bampair.RecordProcessor

	// RecordProcessorsDone should be called after the RecordProcessors
	// have seen all the input records.
	RecordProcessorsDone()

	// Detect identifies the optical duplicates in pairs and returns
	// their names in a slice. readGroupLibrary maps readGroup to
	// library name. pairs contains all the readpairs in the bag, and
	// bestIndex is an index into pairs that points to the bag's
	// primary readpair.
	Detect(readGroupLibrary map[string]string, pairs []DuplicateEntry, bestIndex int) []string
}

// Opts for mark-duplicates.
type Opts struct {
	// Commandline options.
	BamFile              string
	IndexFile            string
	MetricsFile          string
	Format               string
	ShardSize            int
	MinBases             int
	Padding              int
	DiskMateShards       int
	ScratchDir           string
	Parallelism          int
	QueueLength          int
	ClearExisting        bool
	RemoveDups           bool
	TagDups              bool
	IntDI                bool
	UseUmis              bool
	UmiFile              string
	ScavengeUmis         int
	EmitUnmodifiedFields bool
	SeparateSingletons   bool
	OutputPath           string
	StrandSpecific       bool
	OpticalHistogram     string
	OpticalHistogramMax  int

	// Data and operators derived from commandline options.
	BagProcessorFactories []BagProcessorFactory
	OpticalDetector       OpticalDetector
	KnownUmis             []byte
}

type duplicateMatcher interface {
	insertSingleton(r *sam.Record, fileIdx uint64)
	insertPair(a, b *sam.Record, aFileIdx, bFileIdx uint64)
	computeDupSets(*MetricsCollection)
	nextDupSet() (*duplicateSet, bool)
}

type maxAlignDistCheck struct {
	clearExisting      bool
	padding            int
	maxAlignDist       int
	maxX               int
	maxY               int
	globalMaxAlignDist *int
	mutex              *sync.Mutex
}

func (m *maxAlignDistCheck) Process(r *sam.Record) error {
	if m.clearExisting {
		clearDupFlagTags(r)
	}

	d := r.Pos - bam.UnclippedFivePrimePosition(r)
	if d < 0 {
		d = -d
	}
	if d > m.padding {
		return fmt.Errorf("5' alignment distance(%d) exceeds padding(%d) on read: %v", d, m.padding, r.Name)
	}
	if d > m.maxAlignDist {
		m.maxAlignDist = d
	}
	return nil
}

func (m *maxAlignDistCheck) Close() {
	log.Debug.Printf("maximum alignment distance: %d", m.maxAlignDist)
	m.mutex.Lock()
	defer m.mutex.Unlock()
	if m.maxAlignDist > *m.globalMaxAlignDist {
		*m.globalMaxAlignDist = m.maxAlignDist
	}
}

// MarkDuplicates implements duplicate marking.
type MarkDuplicates struct {
	Provider           bamprovider.Provider
	Opts               *Opts
	shardList          []bam.Shard
	readGroupLibrary   map[string]string
	umiCorrector       *umi.SnapCorrector
	distantMates       *bampair.DistantMateTable
	shardInfo          *bampair.ShardInfo
	globalMetrics      *MetricsCollection
	globalMaxAlignDist int
	mutex              sync.Mutex
}

// Mark marks the duplicates, and returns metrics, and an error if encountered.
func (m *MarkDuplicates) Mark(shards []bam.Shard) (*MetricsCollection, error) {
	header, err := m.Provider.GetHeader()
	if err != nil {
		return nil, err
	}

	if shards == nil {
		m.shardList, err = m.Provider.GenerateShards(bamprovider.GenerateShardsOpts{
			Strategy:                           bamprovider.ByteBased,
			Padding:                            m.Opts.Padding,
			IncludeUnmapped:                    true,
			BytesPerShard:                      int64(m.Opts.ShardSize),
			MinBasesPerShard:                   m.Opts.MinBases,
			SplitUnmappedCoords:                false,
			SplitMappedCoords:                  false,
			AlwaysSplitMappedAndUnmappedCoords: true,
		})
	} else {
		m.shardList = shards
	}
	if err != nil {
		return nil, err
	}
	// Collect some info from the bam header
	m.readGroupLibrary = make(map[string]string)
	for _, readGroup := range header.RGs() {
		m.readGroupLibrary[readGroup.Name()] = readGroup.Library()
	}

	// Create umi corrector.
	if m.Opts.KnownUmis != nil {
		m.umiCorrector = umi.NewSnapCorrector(m.Opts.KnownUmis)
	}

	m.globalMetrics = newMetricsCollection()

	// Scan the file once to find each distant mate, and save them to distantMates.
	log.Debug.Printf("Scanning %d shards", len(m.shardList))
	distantMatesOpts := &bampair.Opts{
		Parallelism: m.Opts.Parallelism,
		DiskShards:  m.Opts.DiskMateShards,
		ScratchDir:  m.Opts.ScratchDir,
	}

	recordProcessors := []func() bampair.RecordProcessor{
		func() bampair.RecordProcessor {
			return &maxAlignDistCheck{
				clearExisting:      m.Opts.ClearExisting,
				padding:            m.Opts.Padding,
				globalMaxAlignDist: &m.globalMaxAlignDist,
				mutex:              &m.mutex,
			}
		},
	}
	if m.Opts.OpticalDetector != nil {
		recordProcessors = append(recordProcessors, m.Opts.OpticalDetector.GetRecordProcessor)
	}

	distantMates, shardInfo, err := bampair.GetDistantMates(m.Provider, m.shardList,
		distantMatesOpts, recordProcessors)
	if err != nil {
		return nil, fmt.Errorf("failed while scanning for distant mates: %v", err)
	}
	m.distantMates = distantMates
	m.shardInfo = shardInfo
	m.globalMetrics.maxAlignDist = m.globalMaxAlignDist
	if m.Opts.OpticalDetector != nil {
		m.Opts.OpticalDetector.RecordProcessorsDone()
	}

	for i := 0; i < m.shardInfo.Len(); i++ {
		log.Debug.Printf("shard[%d] info: %v", i, m.shardInfo.GetInfoByIdx(i))
	}

	switch bamprovider.ParseFileType(m.Opts.Format) {
	case bamprovider.BAM:
		err = m.generateBAM()
	case bamprovider.PAM:
		err = m.generatePAM()
	}
	if err != nil {
		return nil, err
	}
	return m.globalMetrics, nil
}

type pamOutputShard struct {
	index     int // 0, 1, ...
	fileShard bam.Shard
	fileRange biopb.CoordRange
	remaining []bam.Shard
}

func newPAMShardsWriter(header *sam.Header, fileShards []bam.Shard, readShards []bam.Shard) ([]*pamOutputShard, error) {
	s := make([]*pamOutputShard, len(fileShards))
	j := 0
	for i := range fileShards {
		ps := &pamOutputShard{fileShard: fileShards[i]}
		s[i] = ps
		ps.index = i
		ps.fileRange = bam.ShardToCoordRange(fileShards[i])

		// Collect read shards in the filerange. The read range boundaries are
		// always aligned at file shard boundaries.
		r := []bam.Shard{}
		for j < len(readShards) {
			readRange := bam.ShardToCoordRange(readShards[j])
			if readRange.Start.GE(ps.fileRange.Limit) {
				break
			}
			if !ps.fileRange.ContainsRange(readRange) {
				log.Fatalf("fileRange %v, readrange %v", ps.fileRange, readRange)
			}
			r = append(r, readShards[j])
			j++
		}
		if len(r) == 0 {
			return nil, fmt.Errorf("empty fileRange %v", s[i])
		}
		ps.remaining = r
	}
	if j != len(readShards) {
		log.Fatalf("fileShards %v does not cover the entire readshards range %v", fileShards, readShards)
	}
	return s, nil
}

func (m *MarkDuplicates) generatePAM() error {
	header, err := m.Provider.GetHeader()
	if err != nil {
		return err
	}
	fileShards, err := m.Provider.GetFileShards()
	if err != nil {
		return err
	}
	outputShards, err := newPAMShardsWriter(header, fileShards, m.shardList)
	if err != nil {
		return err
	}

	e := errors.Once{}
	wg := sync.WaitGroup{}

	outShardCh := make(chan *pamOutputShard, len(outputShards))
	nShards := len(outputShards)
	outShardCh <- outputShards[nShards-1]
	for i := 0; i < nShards-1; i++ {
		outShardCh <- outputShards[i]
	}
	close(outShardCh)
	for wi := 0; wi < m.Opts.Parallelism; wi++ {
		wg.Add(1)
		go func() {
			defer wg.Done()
			for outShard := range outShardCh {
				opts := pam.WriteOpts{
					Range: outShard.fileRange,
				}
				if !m.Opts.EmitUnmodifiedFields {
					opts.DropFields = []bam.FieldType{
						bam.FieldMapq,
						bam.FieldMateRefID, bam.FieldMatePos,
						bam.FieldTempLen,
						bam.FieldName,
						bam.FieldSeq,
						bam.FieldQual}
				}
				writer := pam.NewWriter(opts, header, m.Opts.OutputPath)
				for len(outShard.remaining) > 0 {
					bs := outShard.remaining[0]
					outShard.remaining = outShard.remaining[1:]
					log.Debug.Printf("file %d: starting shard %s, %d remaining", outShard.index, bs.String(), len(outShard.remaining))
					iter := m.Provider.NewIterator(bs)
					m.processShard(iter, bs, outShard.index, func(r *sam.Record) {
						writer.Write(r)
						sam.PutInFreePool(r)
					})
					e.Set(iter.Close())
					log.Debug.Printf("file %d: finished shard %s, %d remaining", outShard.index, bs.String(), len(outShard.remaining))
				}
				e.Set(writer.Close())
				log.Debug.Printf("file %d: all done", outShard.index)
			}
		}()
	}
	wg.Wait()
	return e.Err()
}

func (m *MarkDuplicates) generateBAM() error {
	ctx := vcontext.Background()
	// Prepare outputs.
	var outputStream io.Writer
	if m.Opts.OutputPath == "" {
		outputStream = os.Stdout
	} else {
		out, err := file.Create(ctx, m.Opts.OutputPath)
		if err != nil {
			log.Fatalf("Couldn't create output file %s: %v", m.Opts.OutputPath, err)
		}
		defer func() {
			if err := out.Close(ctx); err != nil {
				log.Fatalf("close %s: %v", m.Opts.OutputPath, err)
			}
		}()
		outputStream = out.Writer(ctx)
	}
	header, err := m.Provider.GetHeader()
	if err != nil {
		log.Fatalf("Could not read header from provider %s: %s", m.Provider, err)
	}
	var writer *bam.ShardedBAMWriter
	if writer, err = bam.NewShardedBAMWriter(outputStream, gzip.DefaultCompression,
		m.Opts.QueueLength, header); err != nil {
		log.Fatalf("Couldn't create bam writer for %s: %v", m.Opts.OutputPath, err)
	}

	// Create workers to process shards off the shardChannel.
	t0 := time.Now()
	var workerGroup sync.WaitGroup
	shardChannel := make(chan bam.Shard, len(m.shardList))
	// The last shard is the unmapped (which can be very large), so
	// move it to the front to process it first.
	unmappedShard := m.shardList[len(m.shardList)-1]
	m.shardList = m.shardList[0 : len(m.shardList)-1]
	if unmappedShard.EndRef != nil {
		log.Fatalf("expected unmapped shard to be last, instead got %v", unmappedShard)
	}
	shardChannel <- unmappedShard
	for _, shard := range m.shardList {
		shardChannel <- shard
	}
	close(shardChannel)

	log.Debug.Printf("Creating %d workers", m.Opts.Parallelism)
	for i := 0; i < m.Opts.Parallelism; i++ {
		workerGroup.Add(1)
		go func(worker int) {
			defer workerGroup.Done()
			compressor := writer.GetCompressor()
			for {
				shard, ok := <-shardChannel
				if !ok {
					break
				}
				log.Debug.Printf("starting shard %s", shard.String())
				if err := compressor.StartShard(shard.ShardIdx); err != nil {
					log.Fatalf("could not create bam shard: %v", err)
				}
				iter := m.Provider.NewIterator(shard)
				m.processShard(iter, shard, worker, func(r *sam.Record) {
					if err := compressor.AddRecord(r); err != nil {
						panic(err)
					}
				})
				if err := iter.Close(); err != nil {
					log.Fatalf("close shard %d: %s", shard.ShardIdx, err)
				}
				// Close the shard (this will block if the queue is full)
				if err := compressor.CloseShard(); err != nil {
					log.Fatalf("close shard compressor %d: %v", shard.ShardIdx, err)
				}
			}
		}(i)
	}
	workerGroup.Wait()
	t1 := time.Now()
	log.Debug.Printf("workers all done in %v", t1.Sub(t0))

	// Close distantMates to clean up any files it may have created.
	if err := m.distantMates.Close(); err != nil {
		log.Fatalf("Error while closing distant mates: %v", err)
	}

	// Wait for the writer to finish writing and then close.
	if err := writer.Close(); err != nil {
		log.Fatalf("Error while closing bam: %v", err)
	}
	t2 := time.Now()
	log.Debug.Printf("closed writer in %v ms", t2.Sub(t1))

	return nil
}

func updateMetrics(readGroupLibrary map[string]string, MetricsCollection *MetricsCollection, record *sam.Record) {
	library := GetLibrary(readGroupLibrary, record)
	metrics := MetricsCollection.Get(library)

	if (record.Flags & sam.Unmapped) != 0 {
		metrics.UnmappedReads++
	} else if bam.HasNoMappedMate(record) &&
		(record.Flags&sam.Secondary) == 0 && (record.Flags&sam.Supplementary) == 0 {
		metrics.UnpairedReads++
	}

	if (record.Flags&sam.Paired) != 0 &&
		(record.Flags&sam.Unmapped) == 0 && (record.Flags&sam.MateUnmapped) == 0 &&
		(record.Flags&sam.Secondary) == 0 && (record.Flags&sam.Supplementary) == 0 {
		metrics.ReadPairsExamined++
	}
	if (record.Flags&sam.Secondary) != 0 || (record.Flags&sam.Supplementary) != 0 {
		metrics.SecondarySupplementary++
	}
}

func (m *MarkDuplicates) processShard(
	iter bamprovider.Iterator,
	shard bam.Shard,
	worker int,
	writeCallback func(*sam.Record)) {

	header, err := m.Provider.GetHeader()
	if err != nil {
		log.Fatalf("error getting header: %v", err)
	}

	if err := m.distantMates.OpenShard(shard.ShardIdx); err != nil {
		log.Fatalf("error opening distant mate shard: %v", err)
	}
	defer m.distantMates.CloseShard(shard.ShardIdx)
	t0 := time.Now()
	orderedReads := []*sam.Record{}
	pairsByName := make(map[string]*readPair)
	singlesByName := make(map[string]*readPair)

	var matcher duplicateMatcher = newDuplicateIndex(worker, header, m.readGroupLibrary, m.Opts, m.umiCorrector)
	MetricsCollection := newMetricsCollection()
	pending := make(map[string]bool)
	readCount := 0

	// readIdx is the index of each read, zeroed at the start of
	// the padded shard.  We add readIdx to
	// shardInfo.PaddingStartFileIdx to calculate the global file
	// index of each read.
	readIdx := uint64(0)

	for iter.Scan() {
		record := iter.Record()
		if m.Opts.ClearExisting {
			clearDupFlagTags(record)
		}

		// In the unmapped shard (record.Ref == nil), all records are in the shard.
		if shard.RecordInShard(record) {
			updateMetrics(m.readGroupLibrary, MetricsCollection, record)
		}

		// Compress reads in the unmapped shard right away instead
		// of storing in orderedReads to limit memory consumption.
		if record.Ref == nil && shard.RecordInShard(record) {
			writeCallback(record)
			continue
		}
		orderedReads = append(orderedReads, record)

		if (record.Flags&sam.Secondary) != 0 || (record.Flags&sam.Supplementary) != 0 {
			log.Debug.Printf("Ignoring secondary or supplementary read: %s", record.Name)
		} else if (record.Flags & sam.Unmapped) != 0 {
			// Pass through Secondary alignments and Unmapped records.
			log.Debug.Printf("Ignoring unmapped read: %s", record.Name)
		} else if !shard.RecordInPaddedShard(record) &&
			!mateInPaddedShard(&shard, record) {
			log.Debug.Printf("Ignoring read outside of padding: %s", record.Name)
		} else if bam.HasNoMappedMate(record) {
			// Handle reads with an unmapped mate differently.
			info := m.shardInfo.GetInfoByShard(&shard)
			singlesByName[record.Name] = &readPair{
				left:        record,
				leftFileIdx: readIdx + info.PaddingStartFileIdx,
			}

			matcher.insertSingleton(record, readIdx+info.PaddingStartFileIdx)
			record = nil // Don't put back in the free pool.
		} else {
			// If we reach here, this read is mapped, it is in the
			// padded shard, and it also has a mapped mate, so we
			// should be able to form a pair.
			var pair *readPair
			var ok bool
			completedPair := false

			// Get info by shard even if this record is not in
			// shard.  This is ok because we will correct for
			// records we see in the padding.
			info := m.shardInfo.GetInfoByShard(&shard)

			if mateInPaddedShard(&shard, record) {
				log.Debug.Printf("read %s should be within shard %v info %v", record.Name, shard, info)
				// Mate is in this shard including padding, so check if we saw it already
				pair, ok = pairsByName[record.Name]
				if ok {
					log.Debug.Printf("Found second read %s %v local readIdx %d", record.Name,
						record.Start(), readIdx)
					pair.addRead(record, readIdx+info.PaddingStartFileIdx)
					completedPair = true
					delete(pending, record.Name)
				} else {
					log.Debug.Printf("Found first read %s %v local readIdx %d", record.Name,
						record.Start(), readIdx)
					pairsByName[record.Name] = &readPair{record, nil, readIdx + info.PaddingStartFileIdx, 0}
					pending[record.Name] = true
				}
			} else {
				// Mate is in another ref or is outside this padded
				// shard, so its mate should be in distantMates.
				log.Debug.Printf("read %s has distant mate: different ref %v, distance %v",
					record.Name, record.Ref.ID() != record.MateRef.ID(), abs(record.Pos-record.MatePos))
				mate, mateFileIdx := m.distantMates.GetMate(shard.ShardIdx, record)
				if mate == nil {
					log.Fatalf("record %v, is missing distant mate, check that both reads are present and "+
						"bai index is valid", record)
				}

				if m.Opts.ClearExisting {
					clearDupFlagTags(mate)
				}

				// Make sure to clone the record below from
				// distantPairs because flagDuplicates() will
				// modify the record and make DistantMateTable
				// misbehave.
				clone := *mate
				log.Debug.Printf("adding distant mate as pair for %s", record.Name)
				pair = &readPair{record, nil, readIdx + info.PaddingStartFileIdx, 0}
				pair.addRead(&clone, mateFileIdx)

				completedPair = true
				pairsByName[record.Name] = pair
				log.Debug.Printf("pair is now %s", pair)
			}

			if completedPair {
				matcher.insertPair(pair.left, pair.right, pair.leftFileIdx, pair.rightFileIdx)
			}
		}
		readIdx++
	}
	for name := range pending {
		log.Error.Printf("Could not find mate for pending read: %v", name)
	}
	if len(pending) > 0 {
		log.Fatalf("Could not find mate for some reads")
	}
	t1 := time.Now()

	// Detect and mark duplicates.
	dupMetrics := flagDuplicates(m.Opts, &shard, m.readGroupLibrary, singlesByName, pairsByName, matcher)
	MetricsCollection.Merge(dupMetrics)
	t2 := time.Now()

	// Compress and write records.
	for _, r := range orderedReads {
		if r.Ref == nil {
			continue
		}
		if shard.RecordInShard(r) {
			if !m.Opts.RemoveDups || (r.Flags&sam.Duplicate) == 0 {
				writeCallback(r)
			}
		}
	}
	readCount += len(orderedReads)
	t3 := time.Now()

	// Update global metrics.
	m.globalMetrics.Merge(MetricsCollection)
	t4 := time.Now()

	log.Debug.Printf("worker %d finished shard %s, reads %d, process %v , mark %v, compress %v, metrics %v, total %v",
		worker, shard.String(), readCount, t1.Sub(t0), t2.Sub(t1), t3.Sub(t2), t4.Sub(t3), t4.Sub(t0))
}

func flagRead(opts *Opts, r *sam.Record, primary, optical bool, dupSetId uint64, dupSetSize, pcrDupSetSize int,
	corrected string) {
	if opts.TagDups && dupSetSize >= 0 {
		var tag sam.Aux
		var err error
		if dupSetSize >= 0 {
			if opts.IntDI {
				tag, err = sam.NewAux(diTag, int(dupSetId))
				if err != nil {
					log.Fatalf("error creating DI:i:%d tag: %v", dupSetId, err)
				}
			} else {
				tag, err = sam.NewAux(diTag, strconv.FormatUint(dupSetId, 10))
				if err != nil {
					log.Fatalf("error creating DI:Z:%d tag: %v", dupSetId, err)
				}
			}
			r.AuxFields = append(r.AuxFields, tag)
		}

		if dupSetSize >= 0 {
			tag, err = sam.NewAux(dsTag, dupSetSize)
			if err != nil {
				log.Fatalf("error creating DS:i:%d tag: %v", dupSetSize, err)
			}
			r.AuxFields = append(r.AuxFields, tag)
		}
		if pcrDupSetSize >= 0 {
			tag, err = sam.NewAux(dlTag, pcrDupSetSize)
			if err != nil {
				log.Fatalf("error creating DL:i:%d tag: %v", pcrDupSetSize, err)
			}
			r.AuxFields = append(r.AuxFields, tag)
		}

		if opts.TagDups && dupSetSize > 1 && len(corrected) > 0 {
			tag, err = sam.NewAux(duTag, corrected)
			if err != nil {
				log.Fatalf("error creating DU:Z:%s tag: %v", corrected, err)
			}
			r.AuxFields = append(r.AuxFields, tag)
		}
	}
	if !primary {
		r.Flags |= sam.Duplicate
		if opts.TagDups && opts.OpticalDetector != nil {
			if optical {
				tag, err := sam.NewAux(dtTag, "SQ")
				if err != nil {
					log.Fatalf("error creating DT:z:SQ tag: %v", err)
				}
				r.AuxFields = append(r.AuxFields, tag)
			} else {
				tag, err := sam.NewAux(dtTag, "LB")
				if err != nil {
					log.Fatalf("error creating DT:z:LB tag: %v", err)
				}
				r.AuxFields = append(r.AuxFields, tag)
			}
		}
	}
}

// SetupAndMark does some minimal setup for validating opts, and
// creating provider and then runs mark().
func SetupAndMark(ctx context.Context, provider bamprovider.Provider, opts *Opts) error {
	if err := validate(opts); err != nil {
		return err
	}

	// Prepare umi inputs.
	if len(opts.UmiFile) > 0 {
		var err error
		umiReader, err := file.Open(ctx, opts.UmiFile)
		if err != nil {
			log.Debug.Printf("Could not read umi file %s: %s", opts.UmiFile, err)
			return err
		}
		defer umiReader.Close(ctx) // nolint: errcheck
		opts.KnownUmis, err = ioutil.ReadAll(umiReader.Reader(ctx))
		if err != nil {
			log.Debug.Printf("Could not read umi file %s: %s", opts.UmiFile, err)
			return err
		}
		if len(opts.KnownUmis) == 0 {
			log.Debug.Printf("UMI list is empty: %s", opts.UmiFile)
			return err
		}
	}

	// Mark/remove those duplicates.
	markDuplicates := &MarkDuplicates{
		Provider: provider,
		Opts:     opts,
	}
	globalMetrics, err := markDuplicates.Mark(nil)
	if err != nil {
		log.Debug.Printf("Error marking duplicates: %v", err)
		return err
	}

	// Output metric and histogram files.
	if opts.MetricsFile != "" {
		if err := writeMetrics(ctx, opts, globalMetrics); err != nil {
			return err
		}
	}
	if opts.OpticalHistogram != "" {
		if err := writeOpticalHistogram(ctx, opts, globalMetrics); err != nil {
			return err
		}
	}
	return nil
}

func flagDuplicates(opts *Opts, shard *bam.Shard, readGroupLibrary map[string]string, singlesByName map[string]*readPair,
	pairsByName map[string]*readPair, matcher duplicateMatcher) *MetricsCollection {
	dupMetrics := newMetricsCollection()

	matcher.computeDupSets(dupMetrics)
	for {
		dupSet, ok := matcher.nextDupSet()
		if !ok {
			break
		}

		optDups := map[string]bool{}
		for _, name := range dupSet.opticals {
			optDups[name] = true
		}

		dupSetId := uint64(0)
		for i, qname := range dupSet.pairs {
			p := pairsByName[qname]
			if i == 0 {
				dupSetId = p.leftFileIdx
			}

			// The pair may contain a read from a different shard, so
			// verify the read is inShard before marking and counting.
			for _, r := range []*sam.Record{p.left, p.right} {
				if shard.RecordInShard(r) {
					if i == 0 {
						log.Debug.Printf("marking %s as primary of DI %d", r.Name, dupSetId)
						flagRead(opts, r, true, false, dupSetId, len(dupSet.pairs), len(dupSet.pairs)-len(optDups)-1,
							dupSet.corrected[r.Name])
					} else {
						log.Debug.Printf("marking %s as duplicate of DI %d optical %v", r.Name, dupSetId, optDups[qname])
						flagRead(opts, r, false, optDups[qname], dupSetId, len(dupSet.pairs), len(dupSet.pairs)-len(optDups)-1,
							dupSet.corrected[r.Name])
						metrics := dupMetrics.Get(GetLibrary(readGroupLibrary, r))
						metrics.ReadPairDups++
						if optDups[qname] {
							metrics.ReadPairOpticalDups++
						}
					}
				}
			}
		}
		for i, qname := range dupSet.singles {
			p := singlesByName[qname]
			if shard.RecordInShard(p.left) {
				// A mate-unmapped read cannot be an optical dup.  A
				// mate-unmapped read cannot be associated with a
				// particular dupSetId, or dupSetSize, even if the
				// only duplicates are also mate-unmapped (this
				// behavior is copied from picard).
				flagRead(opts, p.left, len(dupSet.pairs) == 0 && i == 0, false, 0, -1, -1, dupSet.corrected[p.left.Name])
				if len(dupSet.pairs) == 0 && i > 0 || len(dupSet.pairs) > 0 {
					metrics := dupMetrics.Get(GetLibrary(readGroupLibrary, p.left))
					metrics.UnpairedDups++
				}
			}
		}
	}
	return dupMetrics
}
