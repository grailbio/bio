package main

import (
	"fmt"
	"log"
	"runtime"

	gbam "github.com/grailbio/bio/encoding/bam"
	"github.com/grailbio/bio/encoding/bamprovider"
	"github.com/grailbio/hts/sam"
)

type aggrFlagstat struct {
	total         int
	mapped        int
	duplicate     int
	secondary     int
	supplementary int
	paired        int
	goodPair      int
	single        int
	pairMap       int
	diffChr       int
	diffHigh      int
	r1, r2        int
}

func (stat *aggrFlagstat) mergeFrom(src aggrFlagstat) {
	stat.total += src.total
	stat.mapped += src.mapped
	stat.duplicate += src.duplicate
	stat.secondary += src.secondary
	stat.supplementary += src.supplementary
	stat.paired += src.paired
	stat.goodPair += src.goodPair
	stat.single += src.single
	stat.pairMap += src.pairMap
	stat.diffChr += src.diffChr
	stat.diffHigh += src.diffHigh
	stat.r1 += src.r1
	stat.r2 += src.r2
}

func (stat *aggrFlagstat) record(r *sam.Record) {
	stat.total++
	f := r.Flags
	if (f & sam.Unmapped) == 0 {
		stat.mapped++
	}
	if (f & sam.Duplicate) != 0 {
		stat.duplicate++
	}
	if (f & sam.Secondary) != 0 {
		stat.secondary++
	} else if (f & sam.Supplementary) != 0 {
		stat.supplementary++
	} else if (f & sam.Paired) != 0 {
		stat.paired++
		if (f&sam.ProperPair) != 0 && (f&sam.Unmapped) == 0 {
			stat.goodPair++
		}
		if (f & sam.Read1) != 0 {
			stat.r1++
		}
		if (f & sam.Read2) != 0 {
			stat.r2++
		}
		if (f&sam.MateUnmapped) != 0 && (f&sam.Unmapped) == 0 {
			stat.single++
		}
		if (f&sam.Unmapped) == 0 && (f&sam.MateUnmapped) == 0 {
			stat.pairMap++
			if r.Ref.ID() != r.MateRef.ID() {
				stat.diffChr++
				if r.MapQ >= 5 {
					stat.diffHigh++
				}
			}
		}
	}
}

func percent(a int, b int) string {
	if b == 0 {
		return "N/A"
	}
	return fmt.Sprintf("%.2f%%", float64(a)*100/float64(b))
}

func flagstat(path string) error {
	provider := bamprovider.NewProvider(path, bamprovider.ProviderOpts{
		DropFields: []gbam.FieldType{
			gbam.FieldCigar,
			gbam.FieldMatePos,
			gbam.FieldTempLen,
			gbam.FieldName,
			gbam.FieldSeq,
			gbam.FieldQual,
			gbam.FieldAux,
		}})
	shards, err := provider.GenerateShards(bamprovider.GenerateShardsOpts{
		IncludeUnmapped:     true,
		SplitMappedCoords:   true,
		SplitUnmappedCoords: true,
	})
	if err != nil {
		return err
	}
	shardCh := gbam.NewShardChannel(shards)
	qcCh := make(chan aggrFlagstat, len(shards))
	failedCh := make(chan aggrFlagstat, len(shards))
	for i := 0; i < runtime.NumCPU(); i++ {
		go func() {
			for shard := range shardCh {
				qcStats := aggrFlagstat{}
				failedStats := aggrFlagstat{}
				iter := provider.NewIterator(shard)
				for iter.Scan() {
					rec := iter.Record()
					stat := &qcStats
					if (rec.Flags & sam.QCFail) != 0 {
						stat = &failedStats
					}
					stat.record(rec)
					sam.PutInFreePool(rec)
				}
				if err := iter.Close(); err != nil {
					log.Panicf("%v: close shard %v: %v", path, shard, err)
				}
				qcCh <- qcStats
				failedCh <- failedStats
			}
		}()
	}
	qc := aggrFlagstat{}
	failed := aggrFlagstat{}
	for range shards {
		qc.mergeFrom(<-qcCh)
		failed.mergeFrom(<-failedCh)
	}
	if err := provider.Close(); err != nil {
		return err
	}
	fmt.Printf("%d + %d in total (QC-passed reads + QC-failed reads)\n", qc.total, failed.total)
	fmt.Printf("%d + %d secondary\n", qc.secondary, failed.secondary)
	fmt.Printf("%d + %d supplementary\n", qc.supplementary, failed.supplementary)
	fmt.Printf("%d + %d duplicates\n", qc.duplicate, failed.duplicate)
	fmt.Printf("%d + %d mapped (%s:%s)\n", qc.mapped, failed.mapped,
		percent(qc.mapped, qc.total), percent(failed.mapped, failed.total))
	fmt.Printf("%d + %d paired in sequencing\n", qc.paired, failed.paired)
	fmt.Printf("%d + %d read1\n", qc.r1, failed.r1)
	fmt.Printf("%d + %d read2\n", qc.r2, failed.r1)
	fmt.Printf("%d + %d properly paired (%s:%s)\n", qc.goodPair, failed.goodPair,
		percent(qc.goodPair, qc.paired), percent(failed.goodPair, failed.paired))
	fmt.Printf("%d + %d with itself and mate mapped\n", qc.pairMap, failed.pairMap)
	fmt.Printf("%d + %d singletons (%s:%s)\n", qc.single, failed.single,
		percent(qc.single, qc.total), percent(failed.single, failed.total))
	fmt.Printf("%d + %d with mate mapped to a different chr\n", qc.diffChr, failed.diffChr)
	fmt.Printf("%d + %d with mate mapped to a different chr (mapQ>=5)\n", qc.diffHigh, failed.diffHigh)
	return nil
}
