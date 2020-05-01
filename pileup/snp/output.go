// Copyright 2020 Grail Inc.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
package snp

import (
	"context"
	"os"
	"strconv"
	"strings"

	"github.com/grailbio/base/file"
	"github.com/grailbio/base/log"
	"github.com/grailbio/base/recordio"
	"github.com/grailbio/base/recordio/recordiozstd"
	"github.com/grailbio/base/tsv"
	"github.com/grailbio/bio/pileup"
	"github.com/grailbio/hts/bgzf"
)

// writeChromPosRef is a convenience function which appends the CHROM/POS/REF
// columns common to the TSV output formats.
// It converts pos from 0-based to 1-based, since for better or worse, our
// domain has settled on "0-based coordinates in binary files, 1-based in text"
// as a standard.
func writeChromPosRef(tsvw *tsv.Writer, refName string, pos PosType, refChar byte) {
	tsvw.WriteString(refName)         // CHROM
	tsvw.WriteUint32(uint32(pos + 1)) // POS (1-based in VCF text)
	tsvw.WriteByte(refChar)
}

func convertPileupRowsToTSV(ctx context.Context, tmpFiles []*os.File, mainPath string, colBitset int, bgzip bool, parallelism int, refNames []string, refSeqs [][]byte) (err error) {
	refPath := mainPath + ".ref.tsv"
	if bgzip {
		refPath = refPath + ".gz"
	}
	var dstRef file.File
	if dstRef, err = file.Create(ctx, refPath); err != nil {
		return
	}
	defer file.CloseAndReport(ctx, dstRef, &err)

	altPath := mainPath + ".alt.tsv"
	if bgzip {
		altPath = altPath + ".gz"
	}
	var dstAlt file.File
	if dstAlt, err = file.Create(ctx, altPath); err != nil {
		return
	}
	defer file.CloseAndReport(ctx, dstAlt, &err)

	// Write headers.
	var refTSV *tsv.Writer
	var altTSV *tsv.Writer
	if !bgzip {
		refTSV = tsv.NewWriter(dstRef.Writer(ctx))
		altTSV = tsv.NewWriter(dstAlt.Writer(ctx))
	} else {
		bgzfRefWriter := bgzf.NewWriter(dstRef.Writer(ctx), parallelism)
		bgzfAltWriter := bgzf.NewWriter(dstAlt.Writer(ctx), parallelism)
		refTSV = tsv.NewWriter(bgzfRefWriter)
		altTSV = tsv.NewWriter(bgzfAltWriter)
		defer func() {
			if e := bgzfRefWriter.Close(); e != nil && err == nil {
				err = e
			}
			if e := bgzfAltWriter.Close(); e != nil && err == nil {
				err = e
			}
		}()
	}
	refTSV.WriteString("#CHROM\tPOS\tREF")
	altTSV.WriteString("#CHROM\tPOS\tREF\tALT")
	if (colBitset & colBitDpRef) != 0 {
		refTSV.WriteString("DP")
	}
	if (colBitset & colBitDpAlt) != 0 {
		altTSV.WriteString("DP")
	}
	perReadStats := ((colBitset & colPerReadMask) != 0)
	var emptyPerReadStats []byte
	if perReadStats {
		if (colBitset & colBitEndDists) != 0 {
			refTSV.WriteString("5P_DISTS\t3P_DISTS")
			altTSV.WriteString("5P_DISTS\t3P_DISTS")
			emptyPerReadStats = append(emptyPerReadStats, ".\t.\t"...)
		}
		if (colBitset & colBitQuals) != 0 {
			refTSV.WriteString("QUALS")
			altTSV.WriteString("QUALS")
			emptyPerReadStats = append(emptyPerReadStats, ".\t"...)
		}
		if (colBitset & colBitFraglens) != 0 {
			refTSV.WriteString("FRAGLENS")
			altTSV.WriteString("FRAGLENS")
			emptyPerReadStats = append(emptyPerReadStats, ".\t"...)
		}
		if (colBitset & colBitStrands) != 0 {
			refTSV.WriteString("STRANDS")
			altTSV.WriteString("STRANDS")
			emptyPerReadStats = append(emptyPerReadStats, ".\t"...)
		}
	}
	// These two columns will be renamed once we've removed
	// targeted_to_tsv_snp2.py (used to create Conta-readable files) from the
	// pipeline.  The basestrand format should be *more* convenient for Conta...
	if (colBitset & colBitHighQ) != 0 {
		refTSV.WriteString("ref_depth_tier1")
		altTSV.WriteString("alt_depth_tier1")
	}
	if (colBitset & colBitLowQ) != 0 {
		refTSV.WriteString("ref_depth_tier2")
		altTSV.WriteString("alt_depth_tier2")
	}
	if err = refTSV.EndLine(); err != nil {
		return
	}
	if err = altTSV.EndLine(); err != nil {
		return
	}
	// Convert temporary-file bodies.
	lastRefID := uint32(0)
	curRefName := refNames[0]
	curRefSeq8 := refSeqs[0]
	for i, f := range tmpFiles {
		if _, err = f.Seek(0, 0); err != nil {
			return
		}
		scanner := recordio.NewScanner(f, recordio.ScannerOpts{
			Unmarshal: unmarshalPileupRow,
		})
		// Possible todo: parallelize pileupRow -> final-output-format rendering.
		// This intermediate-recordio design causes wall-clock time for the entire
		// run to increase by up to ~35% over the old
		// intermediate-TSVs-which-can-be-concatenated design.  (The design change
		// was still made because, if performance is an issue, you should be
		// requesting recordio final output instead of TSV anyway.)
		for scanner.Scan() {
			pr := scanner.Get().(*pileupRow)
			refID := pr.refID
			if refID != lastRefID {
				curRefName = refNames[refID]
				curRefSeq8 = refSeqs[refID]
				lastRefID = refID
			}
			pos := pr.pos
			refBase8 := curRefSeq8[pos]
			refChar := pileup.Seq8ToASCIITable[refBase8]
			writeChromPosRef(refTSV, curRefName, PosType(pos), refChar)
			refBase := PosType(pileup.Seq8ToEnumTable[refBase8])
			if (colBitset & colBitDpRef) != 0 {
				refTSV.WriteUint32(pr.payload.depth)
			}
			if perReadStats {
				if refBase == PosType(pileup.BaseX) {
					refTSV.WritePartialBytes(emptyPerReadStats)
				} else {
					refFeatures := pr.payload.perRead[refBase]
					if len(refFeatures) == 0 {
						refTSV.WritePartialBytes(emptyPerReadStats)
					} else {
						if (colBitset & colBitEndDists) != 0 {
							for _, f := range refFeatures {
								refTSV.WriteCsvUint32(uint32(f.dist5p))
							}
							refTSV.EndCsv()
							for _, f := range refFeatures {
								refTSV.WriteCsvUint32(uint32(f.fraglen - 1 - f.dist5p))
							}
							refTSV.EndCsv()
						}
						if (colBitset & colBitQuals) != 0 {
							for _, f := range refFeatures {
								refTSV.WriteCsvUint32(uint32(f.qual))
							}
							refTSV.EndCsv()
						}
						if (colBitset & colBitFraglens) != 0 {
							for _, f := range refFeatures {
								refTSV.WriteCsvUint32(uint32(f.fraglen))
							}
							refTSV.EndCsv()
						}
						if (colBitset & colBitStrands) != 0 {
							for _, f := range refFeatures {
								refTSV.WriteCsvByte(pileup.StrandTypeToASCIITable[f.strand])
							}
							refTSV.EndCsv()
						}
					}
				}
			}
			counts := &pr.payload.counts
			if (colBitset & colBitHighQ) != 0 {
				refTSV.WriteUint32(counts[refBase][0] + counts[refBase][1])
			}
			if (colBitset & colBitLowQ) != 0 {
				refTSV.WriteByte('0')
			}
			if err = refTSV.EndLine(); err != nil {
				return
			}
			// Do we want to report ALT=N?  Probably want to make this configurable,
			// since this was handled inconsistently in the past...
			// Current choice is to report counts, but nothing else, for Ns.
			// Note that, when stitch=true, mismatch between the two read-sides is
			// treated as N.
			for altBase := PosType(0); altBase < pileup.NBaseEnum; altBase++ {
				if altBase == refBase {
					continue
				}
				altCount := counts[altBase][0] + counts[altBase][1]
				if altCount != 0 {
					writeChromPosRef(altTSV, curRefName, PosType(pos), refChar)
					altTSV.WriteByte(pileup.EnumToASCIITable[altBase])
					if (colBitset & colBitDpAlt) != 0 {
						altTSV.WriteUint32(pr.payload.depth)
					}
					if perReadStats {
						if altBase == PosType(pileup.BaseX) {
							altTSV.WritePartialBytes(emptyPerReadStats)
						} else {
							altFeatures := pr.payload.perRead[altBase]
							if (colBitset & colBitEndDists) != 0 {
								for _, f := range altFeatures {
									altTSV.WriteCsvUint32(uint32(f.dist5p))
								}
								altTSV.EndCsv()
								for _, f := range altFeatures {
									altTSV.WriteCsvUint32(uint32(f.fraglen - 1 - f.dist5p))
								}
								altTSV.EndCsv()
							}
							if (colBitset & colBitQuals) != 0 {
								for _, f := range altFeatures {
									altTSV.WriteCsvUint32(uint32(f.qual))
								}
								altTSV.EndCsv()
							}
							if (colBitset & colBitFraglens) != 0 {
								for _, f := range altFeatures {
									altTSV.WriteCsvUint32(uint32(f.fraglen))
								}
								altTSV.EndCsv()
							}
							if (colBitset & colBitStrands) != 0 {
								for _, f := range altFeatures {
									altTSV.WriteCsvByte(pileup.StrandTypeToASCIITable[f.strand])
								}
								altTSV.EndCsv()
							}
						}
					}
					if (colBitset & colBitHighQ) != 0 {
						altTSV.WriteUint32(altCount)
					}
					if (colBitset & colBitLowQ) != 0 {
						altTSV.WriteByte('0')
					}
					if err = altTSV.EndLine(); err != nil {
						return
					}
				}
			}
		}
		if err = scanner.Err(); err != nil {
			return
		}
		curPath := f.Name()
		if err = f.Close(); err != nil {
			return
		}
		tmpFiles[i] = nil
		// os.Remove returns an error if we try to remove a file that isn't there.
		_ = os.Remove(curPath)
	}
	if err = refTSV.Flush(); err != nil {
		return
	}
	if err = altTSV.Flush(); err != nil {
		return
	}
	if bgzip {
		log.Printf("convertPileupRowsToTSV: done, final results written to %s.{ref,alt}.tsv.gz", mainPath)
	} else {
		log.Printf("convertPileupRowsToTSV: done, final results written to %s.{ref,alt}.tsv", mainPath)
	}
	return
}

func convertPileupRowsToBasestrandRio(ctx context.Context, tmpFiles []*os.File, mainPath string, refNames []string) (err error) {
	var dst file.File
	if dst, err = file.Create(ctx, mainPath+".basestrand.rio"); err != nil {
		return
	}
	defer file.CloseAndReport(ctx, dst, &err)

	out := dst.Writer(ctx)
	// WriteBaseStrandsRio doesn't quite have the interface we want.
	recordWriter := recordio.NewWriter(out, recordio.WriterOpts{
		Marshal:      marshalBaseStrand,
		Transformers: []string{recordiozstd.Name},
	})
	recordWriter.AddHeader(refNamesHeader, strings.Join(refNames, "\000"))
	recordWriter.AddHeader(recordio.KeyTrailer, true)
	var numPiles int
	for i, f := range tmpFiles {
		if _, err = f.Seek(0, 0); err != nil {
			return
		}
		scanner := recordio.NewScanner(f, recordio.ScannerOpts{
			Unmarshal: unmarshalPileupRow,
		})
		for scanner.Scan() {
			pr := scanner.Get().(*pileupRow)
			counts := &pr.payload.counts
			recordWriter.Append(&BaseStrandPile{
				RefID: pr.refID,
				Pos:   pr.pos,
				Counts: [4][2]uint32{
					counts[0],
					counts[1],
					counts[2],
					counts[3],
				},
			})
			numPiles++
		}
		if err = scanner.Err(); err != nil {
			return
		}
		curPath := f.Name()
		if err = f.Close(); err != nil {
			return
		}
		tmpFiles[i] = nil
		// os.Remove returns an error if we try to remove a file that isn't there.
		_ = os.Remove(curPath)
	}
	recordWriter.SetTrailer(baseStrandsRioTrailer(numPiles))
	if err = recordWriter.Finish(); err != nil {
		return
	}
	log.Printf("convertPileupRowsToBasestrandRio: done, final results written to %s.basestrand.rio", mainPath)
	return
}

func flushPlusAndMinusBuf(w *tsv.Writer, plusBufPtr, minusBufPtr *[]byte) {
	if len(*plusBufPtr) == 0 {
		w.WriteByte('.')
	} else {
		w.WritePartialBytes(*plusBufPtr)
		w.EndCsv()
		*plusBufPtr = (*plusBufPtr)[:0]
	}
	if len(*minusBufPtr) == 0 {
		w.WriteByte('.')
	} else {
		w.WritePartialBytes(*minusBufPtr)
		w.EndCsv()
		*minusBufPtr = (*minusBufPtr)[:0]
	}
}

func convertPileupRowsToBasestrandTSV(ctx context.Context, tmpFiles []*os.File, mainPath string, colBitset int, bgzip bool, parallelism int, refNames []string, refSeqs [][]byte) (err error) {
	fullPath := mainPath + ".basestrand.tsv"
	if bgzip {
		fullPath = fullPath + ".gz"
	}
	var dst file.File
	if dst, err = file.Create(ctx, fullPath); err != nil {
		return
	}
	defer file.CloseAndReport(ctx, dst, &err)

	var w *tsv.Writer
	if !bgzip {
		w = tsv.NewWriter(dst.Writer(ctx))
	} else {
		bgzfWriter := bgzf.NewWriter(dst.Writer(ctx), parallelism)
		w = tsv.NewWriter(bgzfWriter)
		defer func() {
			if e := bgzfWriter.Close(); e != nil && err == nil {
				err = e
			}
		}()
	}
	// Note that the recordio format does not include REF.
	w.WriteString("#CHROM\tPOS\tREF\tA+\tA-\tC+\tC-\tG+\tG-\tT+\tT-")
	perReadStats := ((colBitset & colPerReadMask) != 0)
	var emptyPerReadStats []byte
	if perReadStats {
		if (colBitset & colBitEndDists) != 0 {
			w.WriteString("\t5P_DISTS_A+\t5P_DISTS_A-\t5P_DISTS_C+\t5P_DISTS_C-\t5P_DISTS_G+\t5P_DISTS_G-\t5P_DISTS_T+\t5P_DISTS_T-\t3P_DISTS_A+\t3P_DISTS_A-\t3P_DISTS_C+\t3P_DISTS_C-\t3P_DISTS_G+\t3P_DISTS_G-\t3P_DISTS_T+\t3P_DISTS_T-")
			emptyPerReadStats = append(emptyPerReadStats, ".\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t"...)
		}
		if (colBitset & colBitQuals) != 0 {
			w.WriteString("\tQUALS_A+\tQUALS_A-\tQUALS_C+\tQUALS_C-\tQUALS_G+\tQUALS_G-\tQUALS_T+\tQUALS_T-")
			emptyPerReadStats = append(emptyPerReadStats, ".\t.\t.\t.\t.\t.\t.\t.\t"...)
		}
		if (colBitset & colBitFraglens) != 0 {
			w.WriteString("\tFRAGLENS_A+\tFRAGLENS_A-\tFRAGLENS_C+\tFRAGLENS_C-\tFRAGLENS_G+\tFRAGLENS_G-\tFRAGLENS_T+\tFRAGLENS_T-")
			emptyPerReadStats = append(emptyPerReadStats, ".\t.\t.\t.\t.\t.\t.\t.\t"...)
		}
		if (colBitset & colBitStrands) != 0 {
			w.WriteString("\tSTRANDS_A+\tSTRANDS_A-\tSTRANDS_C+\tSTRANDS_C-\tSTRANDS_G+\tSTRANDS_G-\tSTRANDS_T+\tSTRANDS_T-")
			emptyPerReadStats = append(emptyPerReadStats, ".\t.\t.\t.\t.\t.\t.\t.\t"...)
		}
	}
	if err = w.EndLine(); err != nil {
		return
	}
	lastRefID := uint32(0)
	curRefName := refNames[0]
	curRefSeq8 := refSeqs[0]
	plusBuf := make([]byte, 0, 256)
	minusBuf := make([]byte, 0, 256)
	for i := range tmpFiles {
		f := tmpFiles[i]
		if _, err = f.Seek(0, 0); err != nil {
			return
		}
		scanner := recordio.NewScanner(f, recordio.ScannerOpts{
			Unmarshal: unmarshalPileupRow,
		})
		for scanner.Scan() {
			pr := scanner.Get().(*pileupRow)
			refID := pr.refID
			if refID != lastRefID {
				curRefName = refNames[refID]
				curRefSeq8 = refSeqs[refID]
				lastRefID = refID
			}
			pos := pr.pos
			refBase8 := curRefSeq8[pos]
			refChar := pileup.Seq8ToASCIITable[refBase8]
			writeChromPosRef(w, curRefName, PosType(pos), refChar)
			for _, perStrandCounts := range pr.payload.counts[:4] {
				for _, c := range perStrandCounts {
					w.WriteUint32(c)
				}
			}
			if perReadStats {
				if pr.payload.depth == 0 {
					w.WritePartialBytes(emptyPerReadStats)
				} else {
					curPerRead := &pr.payload.perRead
					// Note that this code would be simpler if we added a strand
					// dimension to perRead.  But that has the drawback of significantly
					// increasing the base size of pileupPayload everywhere, just for a
					// currently-rare use case.
					if (colBitset & colBitEndDists) != 0 {
						for _, baseFeatures := range curPerRead {
							if len(baseFeatures) == 0 {
								w.WriteString(".\t.")
							} else {
								for _, f := range baseFeatures {
									curDist5p := uint64(f.dist5p)
									if f.strand == byte(pileup.StrandFwd) {
										plusBuf = strconv.AppendUint(plusBuf, curDist5p, 10)
										plusBuf = append(plusBuf, ',')
									} else {
										minusBuf = strconv.AppendUint(minusBuf, curDist5p, 10)
										minusBuf = append(minusBuf, ',')
									}
								}
								flushPlusAndMinusBuf(w, &plusBuf, &minusBuf)
							}
						}
						for _, baseFeatures := range curPerRead {
							if len(baseFeatures) == 0 {
								w.WriteString(".\t.")
							} else {
								for _, f := range baseFeatures {
									curDist3p := uint64(f.fraglen - 1 - f.dist5p)
									if f.strand == byte(pileup.StrandFwd) {
										plusBuf = strconv.AppendUint(plusBuf, curDist3p, 10)
										plusBuf = append(plusBuf, ',')
									} else {
										minusBuf = strconv.AppendUint(minusBuf, curDist3p, 10)
										minusBuf = append(minusBuf, ',')
									}
								}
								flushPlusAndMinusBuf(w, &plusBuf, &minusBuf)
							}
						}
					}
					if (colBitset & colBitQuals) != 0 {
						for _, baseFeatures := range curPerRead {
							if len(baseFeatures) == 0 {
								w.WriteString(".\t.")
							} else {
								for _, f := range baseFeatures {
									curQual := uint64(f.qual)
									if f.strand == byte(pileup.StrandFwd) {
										plusBuf = strconv.AppendUint(plusBuf, curQual, 10)
										plusBuf = append(plusBuf, ',')
									} else {
										minusBuf = strconv.AppendUint(minusBuf, curQual, 10)
										minusBuf = append(minusBuf, ',')
									}
								}
								flushPlusAndMinusBuf(w, &plusBuf, &minusBuf)
							}
						}
					}
					if (colBitset & colBitFraglens) != 0 {
						for _, baseFeatures := range curPerRead {
							if len(baseFeatures) == 0 {
								w.WriteString(".\t.")
							} else {
								for _, f := range baseFeatures {
									curFraglen := uint64(f.fraglen)
									if f.strand == byte(pileup.StrandFwd) {
										plusBuf = strconv.AppendUint(plusBuf, curFraglen, 10)
										plusBuf = append(plusBuf, ',')
									} else {
										minusBuf = strconv.AppendUint(minusBuf, curFraglen, 10)
										minusBuf = append(minusBuf, ',')
									}
								}
								flushPlusAndMinusBuf(w, &plusBuf, &minusBuf)
							}
						}
					}
					if (colBitset & colBitStrands) != 0 {
						for _, baseFeatures := range curPerRead {
							if len(baseFeatures) == 0 {
								w.WriteString(".\t.")
							} else {
								for _, f := range baseFeatures {
									if f.strand == byte(pileup.StrandFwd) {
										plusBuf = append(plusBuf, "+,"...)
									} else {
										minusBuf = append(minusBuf, "-,"...)
									}
								}
								flushPlusAndMinusBuf(w, &plusBuf, &minusBuf)
							}
						}
					}
				}
			}
			if err = w.EndLine(); err != nil {
				return
			}
		}
		if err = scanner.Err(); err != nil {
			return
		}
		curPath := f.Name()
		if err = f.Close(); err != nil {
			return
		}
		tmpFiles[i] = nil
		// os.Remove returns an error if we try to remove a file that isn't there.
		_ = os.Remove(curPath)
	}
	if err = w.Flush(); err != nil {
		return
	}
	log.Printf("convertPileupRowsToBasestrandTSV: done, final results written to %s", fullPath)
	return
}
