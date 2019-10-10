package cmd

import (
	"context"
	"flag"
	"io"
	"os"
	"testing"
	"time"

	"github.com/grailbio/base/file"
	"github.com/grailbio/base/grail"
	"github.com/grailbio/base/log"
	"github.com/grailbio/base/vcontext"
	"github.com/grailbio/bio/fusion"
	"github.com/grailbio/testutil"
	"github.com/grailbio/testutil/assert"
	"golang.org/x/sync/errgroup"
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
	if err := cacheDir(ctx, s3Dir, *cacheDirFlag); err != nil {
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
		r1:                 *cacheDirFlag + "/smallr1.fastq.gz",
		r2:                 *cacheDirFlag + "/smallr2.fastq.gz",
		transcriptPath:     *cacheDirFlag + "/transcriptome.fa",
		cosmicFusionPath:   *cacheDirFlag + "/small_pairs.txt",
		geneListInputPath:  *cacheDirFlag + "/gene_names.txt"},
		opts)

	const (
		goldenAll      = "/all-2018-12-21.fa"
		goldenFiltered = "/filtered-2018-12-21.fa"
	)

	if *updateGoldenFlag {
		copyFile(ctx, t, s3Dir+goldenAll, tmpDir+"/all.fa")
		copyFile(ctx, t, s3Dir+goldenFiltered, tmpDir+"/filtered.fa")
	} else {
		testutil.CompareFiles(t, tmpDir+"/all.fa", *cacheDirFlag+goldenAll, nil)
		testutil.CompareFiles(t, tmpDir+"/filtered.fa", *cacheDirFlag+goldenFiltered, nil)
	}
}

func TestMain(m *testing.M) {
	shutdown := grail.Init()
	status := m.Run()
	shutdown()
	os.Exit(status)
}

func cacheDir(ctx context.Context, srcDir, dstDir string) error {
	eg := errgroup.Group{}
	lister := file.List(ctx, srcDir, true /*recursive*/)
	for lister.Scan() {
		path := lister.Path()
		eg.Go(func() error {
			return cacheFile(ctx, path, dstDir+path[len(srcDir):])
		})
	}
	err := eg.Wait()
	if e := lister.Err(); e != nil && err == nil {
		err = e
	}
	return err
}

func cacheFile(ctx context.Context, srcPath string, dstPath string) (err error) {
	in, err := file.Open(ctx, srcPath)
	if err != nil {
		return err
	}
	defer file.CloseAndReport(ctx, in, &err)
	srcInfo, err := in.Stat(ctx)
	if err != nil {
		return err
	}
	dstInfo, err := file.Stat(ctx, dstPath)
	if err == nil {
		diff := dstInfo.ModTime().Sub(srcInfo.ModTime())
		if (diff >= -5*time.Second) && (diff <= 5*time.Second) {
			return nil
		}
	}

	out, err := file.Create(ctx, dstPath)
	if err != nil {
		return err
	}
	if _, err := io.Copy(out.Writer(ctx), in.Reader(ctx)); err != nil {
		file.CloseAndReport(ctx, out, &err)
		return err
	}
	file.CloseAndReport(ctx, out, &err)
	return os.Chtimes(dstPath, srcInfo.ModTime(), srcInfo.ModTime())
}
