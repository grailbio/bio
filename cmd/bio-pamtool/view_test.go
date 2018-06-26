package main_test

import (
	"path/filepath"
	"testing"

	"github.com/grailbio/testutil"
	"github.com/stretchr/testify/assert"
	"v.io/x/lib/gosh"
)

func TestView(t *testing.T) {
	if !testutil.IsBazel() {
		t.Skip("not bazel")
	}
	bamPath := testutil.GetFilePath("//go/src/grail.com/bio/encoding/bam/testdata/test.bam")

	sh := gosh.NewShell(nil)
	defer sh.Cleanup()
	pamtoolPath := testutil.GoExecutable(t, "//go/src/github.com/grailbio/bio/cmd/bio-pamtool/bio-pamtool")
	dir := sh.MakeTempDir()
	pamPath := filepath.Join(dir, "test.pam")
	sh.Cmd(pamtoolPath, "convert", bamPath, pamPath).Run()
	assert.NoError(t, sh.Err)

	expected1 := "read1	0	chr1	123	60	10M	=	456	20	ACGTACGTAC	ABCDEFGHIJ	RG:Z:NA12878\n"
	assert.Equal(t, expected1, sh.Cmd(pamtoolPath, "view", "-filter", "rec_name==\"read1\"", bamPath).CombinedOutput())
	assert.Equal(t, expected1, sh.Cmd(pamtoolPath, "view", "-filter", "rec_name==\"read1\"", pamPath).CombinedOutput())
	expected2 := "read2	0	chr2	111	60	10M	=	222	20	ACGTACGTAC	ABCDEFGHIJ	RG:Z:NA12878\n"
	assert.Equal(t, expected2, sh.Cmd(pamtoolPath, "view", "-filter", "rec_name==\"read2\"", bamPath).CombinedOutput())
	assert.Equal(t, expected2, sh.Cmd(pamtoolPath, "view", "-filter", "rec_name==\"read2\"", pamPath).CombinedOutput())
}
