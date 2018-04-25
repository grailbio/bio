package main

import (
	"encoding/binary"
	"encoding/json"
	"fmt"
	"hash"
	"runtime"

	"github.com/biogo/hts/sam"
	"github.com/blainsmith/seahash"
	"github.com/grailbio/base/errorreporter"
	"github.com/grailbio/base/log"
	"github.com/grailbio/base/unsafe"
	gbam "github.com/grailbio/bio/encoding/bam"
	"github.com/grailbio/bio/encoding/bamprovider"
)

type checksumOpts struct {
	// baiPath sets the name of the BAM index file. If empty, bampath+".bai" is used.
	baiPath string

	// all treats all the following bool fields to be true.  If all=true, then the
	// individual values of the following fields are ignored.
	all bool

	// name causes the record names to be added to the checksum.
	name bool
	// mapq causes the mapq tag values to be added to the checksum.
	mapQ bool
	// cigar causes the sequences to be added to the checksum
	cigar bool
	// mate causes the materef and matepos to be added to the checksum
	matePos bool
	// templen causes the templen to be added to the checksum
	tempLen bool
	// seq causes the sequences to be added to the checksum
	seq bool
	// qual causes the sequences to be added to the checksum
	qual bool
	// aux causes the aux tag values to be added to the checksum.
	aux bool
}

// RefChecksum is the checksum of reads for one chromosome.
type refChecksum struct {
	// Name is the name of the reference.
	Name string
	// NRecs is the # Of records found for this reference sequence.
	NRecs int64
	// SumPos is sum of all position values. A quick commutative hash.
	SumPos uint64
	// SumFlags is sum of all flag values.
	SumFlags uint64
	// SumTemplen is the sum of templen values.
	SumTempLen uint64
	// SumMapQ is the sum of mapq values.
	SumMapQ uint64
	// SumMatePos is the sum of mapref and mappos values.
	SumMatePos uint64
	// SumName is sum of all names.
	SumName uint64
	// SumSeq is sum of all seq strings.
	SumSeq uint64
	// SumCigar is sum of all cigar strings.
	SumCigar uint64
	// SumQual is the sum of seq quality values.
	SumQual uint64
	// SumAux is sum of all aux fields.
	SumAux uint64
}

func hashField(h hash.Hash64, pos [8]byte, value []byte) uint64 {
	h.Reset()
	h.Write(pos[:])
	h.Write(value)
	return h.Sum64()
}

func (c *refChecksum) add(r *sam.Record, h hash.Hash64, opts checksumOpts) {
	if c.Name == "" {
		r.Name = "*"
		if r.Ref != nil {
			c.Name = r.Ref.Name()
		}
	}
	c.NRecs++
	c.SumPos += uint64(r.Pos)

	pos := [8]byte{}
	binary.LittleEndian.PutUint32(pos[:], uint32(r.Ref.ID()))
	binary.LittleEndian.PutUint32(pos[4:], uint32(r.Pos))

	value := [16]byte{}
	binary.LittleEndian.PutUint32(value[:4], uint32(r.Flags))
	c.SumFlags += hashField(h, pos, value[:4])

	if opts.all || opts.tempLen {
		binary.LittleEndian.PutUint32(value[:4], uint32(r.TempLen))
		c.SumTempLen += hashField(h, pos, value[:4])
	}
	if opts.all || opts.mapQ {
		binary.LittleEndian.PutUint32(value[:4], uint32(r.MapQ))
		c.SumMapQ += hashField(h, pos, value[:4])
	}
	if opts.all || opts.matePos {
		binary.LittleEndian.PutUint32(value[:4], uint32(r.MateRef.ID()))
		binary.LittleEndian.PutUint32(value[4:], uint32(r.MatePos))
		c.SumMatePos += hashField(h, pos, value[:8])
	}
	if opts.all || opts.name {
		c.SumName += hashField(h, pos, unsafe.StringToBytes(r.Name))
	}
	if opts.all || opts.qual {
		c.SumQual += hashField(h, pos, r.Qual)
	}
	if opts.all || opts.seq {
		// TODO(saito) The below code drop the last base from the checksum.
		c.SumSeq += hashField(h, pos, gbam.UnsafeDoubletsToBytes(r.Seq.Seq[:r.Seq.Length/2]))
	}
	if opts.all || opts.qual {
		c.SumQual += hashField(h, pos, r.Qual)
	}
	if opts.all || opts.cigar {
		// TODO(saito) The below code drop the last base from the checksum.
		c.SumCigar += hashField(h, pos, gbam.UnsafeCigarToBytes(r.Cigar))
	}
	if opts.all || opts.aux {
		h.Reset()
		h.Write(pos[:])
		for _, aux := range r.AuxFields {
			h.Write(aux)
		}
		c.SumAux += h.Sum64()
	}
}

func (c *refChecksum) merge(other refChecksum) {
	if other.Name != "" {
		c.Name = other.Name
	}
	c.NRecs += other.NRecs
	c.SumName += other.SumName
	c.SumPos += other.SumPos
	c.SumMapQ += other.SumMapQ
	c.SumCigar += other.SumCigar
	c.SumFlags += other.SumFlags
	c.SumMatePos += other.SumMatePos
	c.SumTempLen += other.SumTempLen
	c.SumSeq += other.SumSeq
	c.SumQual += other.SumQual
	c.SumAux += other.SumAux
}

// fileChecksum represents the checksum of a file.
type fileChecksum struct {
	Refs     []refChecksum // One for each ref. Index is refid.
	Unmapped refChecksum   // For unmapped reads.
	err      errorreporter.T
}

func (csum *fileChecksum) resizeRefs(minSize int) {
	for len(csum.Refs) < minSize {
		csum.Refs = append(csum.Refs, refChecksum{})
	}
}

func (csum *fileChecksum) add(r *sam.Record, h hash.Hash64, opts checksumOpts) {
	if r.Ref == nil {
		csum.Unmapped.add(r, h, opts)
		return
	}
	csum.resizeRefs(r.Ref.ID() + 1)
	csum.Refs[r.Ref.ID()].add(r, h, opts)
}

func (csum *fileChecksum) merge(other fileChecksum) {
	csum.resizeRefs(len(other.Refs))
	for i := 0; i < len(other.Refs); i++ {
		csum.Refs[i].merge(other.Refs[i])
	}
	csum.Unmapped.merge(other.Unmapped)
	csum.err.Set(other.err.Err())
}

func checksumBAMShard(opts checksumOpts, ch chan gbam.Shard, provider bamprovider.Provider) fileChecksum {
	csum := fileChecksum{}
	h := seahash.New()
	for {
		shard, ok := <-ch
		if !ok {
			break
		}
		iter := provider.NewIterator(shard)
		for iter.Scan() {
			record := iter.Record()
			csum.add(record, h, opts)
			sam.PutInFreePool(record)
		}
		csum.err.Set(iter.Close())
	}
	return csum
}

func checksumFile(bamPath string, opts checksumOpts) fileChecksum {
	var csum fileChecksum
	bopts := bamprovider.ProviderOpts{Index: opts.baiPath}
	if !opts.all && !opts.name {
		bopts.DropFields = append(bopts.DropFields, gbam.FieldName)
	}
	if !opts.all && !opts.mapQ {
		bopts.DropFields = append(bopts.DropFields, gbam.FieldMapq)
	}
	if !opts.all && !opts.tempLen {
		bopts.DropFields = append(bopts.DropFields, gbam.FieldTempLen)
	}
	if !opts.all && !opts.matePos {
		bopts.DropFields = append(bopts.DropFields, gbam.FieldMateRefID, gbam.FieldMatePos)
	}
	if !opts.all && !opts.seq {
		bopts.DropFields = append(bopts.DropFields, gbam.FieldSeq)
	}
	if !opts.all && !opts.cigar {
		bopts.DropFields = append(bopts.DropFields, gbam.FieldCigar)
	}
	if !opts.all && !opts.qual {
		bopts.DropFields = append(bopts.DropFields, gbam.FieldQual)
	}
	if !opts.all && !opts.aux {
		bopts.DropFields = append(bopts.DropFields, gbam.FieldAux)
	}
	provider := bamprovider.NewProvider(bamPath, bopts)
	shardList, err := provider.GenerateShards(bamprovider.GenerateShardsOpts{
		Strategy:            bamprovider.ByteBased,
		IncludeUnmapped:     true,
		SplitUnmappedCoords: true,
		SplitMappedCoords:   true,
	})
	if err != nil {
		csum.err.Set(err)
		return csum
	}
	shardCh := gbam.NewShardChannel(shardList)
	parallelism := runtime.NumCPU()
	resultCh := make(chan fileChecksum, parallelism)
	for i := 0; i < parallelism; i++ {
		go func() {
			resultCh <- checksumBAMShard(opts, shardCh, provider)
		}()
	}
	for i := 0; i < parallelism; i++ {
		csum.merge(<-resultCh)
	}
	csum.err.Set(provider.Close())
	return csum
}

func checksum(path string, opts checksumOpts) error {
	csum := checksumFile(path, opts)
	if csum.err.Err() != nil {
		return csum.err.Err()
	}
	js, err := json.MarshalIndent(csum, "", "  ")
	if err != nil {
		log.Panic(err)
	}
	fmt.Println(string(js))
	return nil
}
