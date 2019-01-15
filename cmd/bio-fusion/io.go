package main

// This file defines fusionWriter and fusionReader. Type fusionWriter dumps
// GeneDB and Fragment into a recordio file, and fusionReader reads them back.
// The recordio file can be used to bypass the 1st phase of the af4 and run only
// the filtering phase.

import (
	"bytes"
	"context"
	"encoding/gob"
	"log"

	"github.com/grailbio/base/file"
	"github.com/grailbio/base/recordio"
	"github.com/grailbio/base/recordio/recordiozstd"
	"github.com/grailbio/bio/fusion"
)

// fusionWriter is for writing GeneDB and Fragments to a recordio file.
// The recordio file can be used to bypass the 1st phase of the af4 and run only
// the filtering phase.
type fusionWriter struct {
	out    file.File
	w      recordio.Writer
	geneDB *fusion.GeneDB
	opts   fusion.Opts
}

const (
	// <fileVersionHeader, fileVersion> is stored in a recordio header.B
	fileVersionHeader = "af4version"
	fileVersion       = "AF4_V1"
)

// fusionFileHeader is stored in the trailer section of the recordio file.
type fusionFileHeader struct {
	// Opts is the list of options used to generate the candidates.
	Opts fusion.Opts
	// Gene is the list of genes registered in geneDB. Kmers are dropped, since they aren't used in the[q 2nd phase.
	Gene []fusion.GeneInfo
}

func encodeGOB(gw *gob.Encoder, v interface{}) {
	if err := gw.Encode(v); err != nil {
		panic(err)
	}
}

func decodeGOB(gr *gob.Decoder, v interface{}) {
	if err := gr.Decode(v); err != nil {
		panic(err)
	}
}

func newFusionWriter(ctx context.Context, outPath string, geneDB *fusion.GeneDB, opts fusion.Opts) *fusionWriter {
	recordiozstd.Init()
	out, err := file.Create(ctx, outPath)
	if err != nil {
		log.Panicf("rio open %v: %v", outPath, err)
	}
	w := recordio.NewWriter(out.Writer(ctx), recordio.WriterOpts{
		Transformers: []string{recordiozstd.Name},
	})
	w.AddHeader(fileVersionHeader, fileVersion)
	w.AddHeader(recordio.KeyTrailer, true)
	return &fusionWriter{out: out, w: w, geneDB: geneDB, opts: opts}
}

// Write adds a fusion candidate. Any error will crash the process.
func (w *fusionWriter) Write(c fusion.Candidate) {
	b := bytes.NewBuffer(nil)
	gw := gob.NewEncoder(b)
	encodeGOB(gw, c)
	w.w.Append(b.Bytes())
}

// Close closes the writer. It must be called exactly once, after writing all
// the candidates.
func (w *fusionWriter) Close(ctx context.Context) {
	b := bytes.NewBuffer(nil)
	gw := gob.NewEncoder(b)
	h := fusionFileHeader{
		Opts: w.opts,
	}
	minGeneID, limitGeneID := w.geneDB.GeneIDRange()
	for gid := minGeneID; gid < limitGeneID; gid++ {
		h.Gene = append(h.Gene, *w.geneDB.GeneInfo(gid))
	}
	encodeGOB(gw, h)
	w.w.SetTrailer(b.Bytes())
	err := w.w.Finish()
	if err != nil {
		log.Panic("close", err)
	}
	if err := w.out.Close(ctx); err != nil {
		log.Panic("close", err)
	}
}

// fusionReader is for reading GeneDB and Fragments from a recordio file created
// by fusionWriter.  The recordio file can be used to bypass the 1st phase of
// the af4 and run only the filtering phase.
type fusionReader struct {
	in     file.File
	r      recordio.Scanner
	geneDB *fusion.GeneDB
	opts   fusion.Opts

	c fusion.Candidate // last candidate read by Scan.
}

func newFusionReader(ctx context.Context, inPath string) *fusionReader {
	in, err := file.Open(ctx, inPath)
	if err != nil {
		log.Panicf("open %s: %v", inPath, err)
	}
	recordiozstd.Init()
	r := recordio.NewScanner(in.Reader(ctx), recordio.ScannerOpts{})
	versionFound := false
	for _, kv := range r.Header() {
		if kv.Key == fileVersionHeader {
			if kv.Value.(string) != fileVersion {
				log.Panicf("TODO: fusion file version mismatch, got %v, expect %v",
					kv.Value.(string), fileVersion)
			}
			versionFound = true
			break
		}
	}
	if !versionFound {
		log.Panic(fileVersionHeader + " not found")
	}
	gr := gob.NewDecoder(bytes.NewReader(r.Trailer()))
	h := fusionFileHeader{}
	decodeGOB(gr, &h)

	geneDB := fusion.NewGeneDB(h.Opts)
	geneDB.PrepopulateGeneInfo(h.Gene)
	return &fusionReader{in: in, r: r, geneDB: geneDB, opts: h.Opts}
}

// Opts returns the fusion options written in the recordio file. This method can be called
// any time.
func (r *fusionReader) Opts() fusion.Opts { return r.opts }

// GeneDB returns the geneDB written in the recordio file. The DB doesn't
// contain kmer information. This method can be called any time.
func (r *fusionReader) GeneDB() *fusion.GeneDB { return r.geneDB }

// Scan reads the next fusion candidate.
//
// REQUIRES: Close hasn't been called.
func (r *fusionReader) Scan() bool {
	if !r.r.Scan() {
		return false
	}
	b := bytes.NewReader(r.r.Get().([]byte))
	gr := gob.NewDecoder(b)
	r.c = fusion.Candidate{}
	decodeGOB(gr, &r.c)
	return true
}

// Get yields the current candidate.
//
// REQUIRES: Last Scan call returned true.
func (r *fusionReader) Get() fusion.Candidate { return r.c }

// Close closes the reader. It must be called exactly once.  Any error will
// crash the writer.
func (r *fusionReader) Close(ctx context.Context) {
	if err := r.r.Err(); err != nil {
		log.Panic(err)
	}
	if err := r.in.Close(ctx); err != nil {
		log.Panic(err)
	}
}
