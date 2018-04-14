package main_test

import (
	"path/filepath"
	"testing"

	"github.com/grailbio/internal/testutil"
	"github.com/stretchr/testify/assert"
	"v.io/x/lib/gosh"
)

func TestVerify(t *testing.T) {
	if !testutil.IsBazel() {
		t.Skip("not bazel")
	}
	bamPath := testutil.GetFilePath("@grailgo//bio/encoding/bam/testdata/170614_WGS_LOD_Pre_Library_B3_27961B_05.merged.10000.bam")
	sh := gosh.NewShell(nil)
	defer sh.Cleanup()
	pamtoolPath := testutil.GoExecutable(t, sh, "@grailgo//vendor/github.com/grailbio/bio/cmd/bio-pamtool/bio-pamtool")
	dir := sh.MakeTempDir()
	pamPath := filepath.Join(dir, "test.pam")
	sh.Cmd(pamtoolPath, "convert", bamPath, pamPath).Run()
	assert.NoError(t, sh.Err)

	bamCsum := sh.Cmd(pamtoolPath, "checksum", "-all", bamPath).Stdout()
	pamCsum := sh.Cmd(pamtoolPath, "checksum", "-all", pamPath).Stdout()
	assert.Equal(t, bamCsum, pamCsum)

	// Test error generation.
	bam2Path := testutil.GetFilePath("@grailgo//bio/encoding/bam/testdata/test.bam")
	bam2Csum := sh.Cmd(pamtoolPath, "checksum", bam2Path).Stdout()
	assert.NotEqual(t, pamCsum, bam2Csum)
}
