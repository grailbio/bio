package main_test

import (
	"path/filepath"
	"testing"

	"github.com/grailbio/testutil"
	"github.com/stretchr/testify/assert"
	"v.io/x/lib/gosh"
)

func TestFlagStat(t *testing.T) {
	if !testutil.IsBazel() {
		t.Skip("not bazel")
	}
	bamPath := testutil.GetFilePath("@grailgo//bio/encoding/bam/testdata/170614_WGS_LOD_Pre_Library_B3_27961B_05.merged.10000.bam")
	sh := gosh.NewShell(nil)
	defer sh.Cleanup()
	dir := sh.MakeTempDir()
	pamPath := filepath.Join(dir, "test.pam")

	pamtoolPath := testutil.GoExecutable(t, sh, "@grailgo//vendor/github.com/grailbio/bio/cmd/bio-pamtool/bio-pamtool")
	sh.Cmd(pamtoolPath, "convert", bamPath, pamPath).Run()
	assert.NoError(t, sh.Err)

	expected := `20042 + 0 in total (QC-passed reads + QC-failed reads)
0 + 0 secondary
42 + 0 supplementary
0 + 0 duplicates
18840 + 0 mapped (94.00%:N/A)
20000 + 0 paired in sequencing
10000 + 0 read1
10000 + 0 read2
17964 + 0 properly paired (89.82%:N/A)
18722 + 0 with itself and mate mapped
76 + 0 singletons (0.38%:N/A)
144 + 0 with mate mapped to a different chr
50 + 0 with mate mapped to a different chr (mapQ>=5)
`
	output := sh.Cmd(pamtoolPath, "flagstat", pamPath).Stdout()
	assert.NoError(t, sh.Err)
	assert.Equal(t, expected, output)
	output = sh.Cmd(pamtoolPath, "flagstat", bamPath).Stdout()
	assert.NoError(t, sh.Err)
	assert.Equal(t, expected, output)
}
