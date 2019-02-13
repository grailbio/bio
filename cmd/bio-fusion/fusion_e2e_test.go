package main

import (
	"context"
	"flag"
	"io"
	"os"
	"testing"

	"github.com/grailbio/base/file"
	"github.com/grailbio/base/grail"
	"github.com/grailbio/base/log"
	"github.com/grailbio/base/vcontext"
	"github.com/grailbio/bio/fusion"
	"github.com/grailbio/testutil"
	"github.com/grailbio/testutil/assert"
	"github.com/grailbio/testutil/benchmark"
)

var (
	manualFlag       = flag.Bool("run-manual-tests", false, "If true, run tests that access AWS.")
	updateGoldenFlag = flag.Bool("update-golden", false, "Update the golden files instead of comparing them.")
	cacheDirFlag     = flag.String("cache-dir", "/tmp/bio-target-rna-fusion-cache", "Temp dir used to store data files")
)

func copyFile(ctx context.Context, t *testing.T, dstPath, srcPath string) {
	log.Printf("Copying %s -> %s", srcPath, dstPath)
	out, err := file.Create(ctx, dstPath)
	assert.NoError(t, err)
	in, err := file.Open(ctx, srcPath)
	assert.NoError(t, err)
	_, err = io.Copy(out.Writer(ctx), in.Reader(ctx))
	assert.NoError(t, err)
	assert.NoError(t, in.Close(ctx))
	assert.NoError(t, out.Close(ctx))
}

func TestEndToEndSmall(t *testing.T) {
	if !*manualFlag {
		t.Skip("not enabled")
	}
	ctx := vcontext.Background()

	const s3Dir = "s3://grail-go-testing/bio-target-rna-fusion/small"
	cacheDir, err := benchmark.CacheDir(ctx, s3Dir, *cacheDirFlag)
	if err != nil {
		t.Fatal(err)
	}
	tmpDir := "/tmp/bio-target-rna-fusion-test"
	opts := fusion.DefaultOpts
	opts.KmerLength = 19
	opts.MaxGenesPerKmer = 2
	opts.MaxProximityDistance = 1000
	opts.MaxProximityGenes = 0
	log.Printf("Starting e2e test using testdata files in %s", cacheDir)
	DetectFusion(ctx, fusionFlags{
		fastaOutputPath:    tmpDir + "/all.fa",
		rioOutputPath:      tmpDir + "/all.rio",
		filteredOutputPath: tmpDir + "/filtered.fa",
		r1:                 cacheDir + "/smallr1.fastq.gz",
		r2:                 cacheDir + "/smallr2.fastq.gz",
		transcriptPath:     cacheDir + "/transcriptome.fa",
		cosmicFusionPath:   cacheDir + "/small_pairs.txt",
		geneListInputPath:  cacheDir + "/gene_names.txt"},
		opts)

	const (
		goldenAll      = "/all-2018-12-21.fa"
		goldenFiltered = "/filtered-2018-12-21.fa"
	)

	if *updateGoldenFlag {
		copyFile(ctx, t, s3Dir+goldenAll, tmpDir+"/all.fa")
		copyFile(ctx, t, s3Dir+goldenFiltered, tmpDir+"/filtered.fa")
	} else {
		testutil.CompareFiles(t, tmpDir+"/all.fa", cacheDir+goldenAll, nil)
		testutil.CompareFiles(t, tmpDir+"/filtered.fa", cacheDir+goldenFiltered, nil)
	}
}

func TestMain(m *testing.M) {
	shutdown := grail.Init()
	status := m.Run()
	shutdown()
	os.Exit(status)
}
