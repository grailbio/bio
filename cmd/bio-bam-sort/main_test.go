package main

import (
	"bytes"
	"fmt"
	"testing"

	"github.com/stretchr/testify/require"
	"grail.com/testutil"
	"v.io/x/lib/gosh"
)

func TestSort(t *testing.T) {
	samData := []string{
		`@HD	VN:1.3	SO:coordinate
@SQ	SN:chr1	LN:10000
@SQ	SN:chr2	LN:9000
read1	0	chr1	123	60	10M	=	456	20	ACGTACGTAC	ABCDEFGHIJ	RG:Z:NA12878
read2	0	chr2	111	60	10M	=	222	20	ACGTACGTAC	ABCDEFGHIJ	RG:Z:NA12878
`,
		`@HD	VN:1.3	SO:coordinate
@SQ	SN:chr1	LN:10000
@SQ	SN:chr2	LN:9000
@SQ	SN:chr3	LN:8000
read1	0	chr1	123	60	10M	=	456	20	ACGTACGTAC	ABCDEFGHIJ	RG:Z:NA12879
read2	0	chr2	111	60	10M	=	222	20	ACGTACGTAC	ABCDEFGHIJ	RG:Z:NA12879
read10	12	*	0	0	10M	*	0	0	ACGTACGTAC	ABCDEFGHIJ	RG:Z:NA12878
read10	12	*	0	0	10M	*	0	0	ACGTACGTAC	ABCDEFGHIJ	RG:Z:NA12879
`,
	}

	sh := gosh.NewShell(nil)
	sort := testutil.GoExecutable(t, sh, "@grailgo//vendor/github.com/grailbio/bio/cmd/bio-bam-sort/bio-bam-sort")
	pamtool := testutil.GoExecutable(t, sh, "@grailgo//vendor/github.com/grailbio/bio/cmd/bio-pamtool/bio-pamtool")
	tempDir := sh.MakeTempDir()
	defer sh.Cleanup()

	shardPaths := []string{}
	for i, sam := range samData {
		shardPath := fmt.Sprintf("%s/%d.sortshard", tempDir, i)
		shardPaths = append(shardPaths, shardPath)
		sortCmd := sh.Cmd(sort, "-shard-index", fmt.Sprintf("%d", i+1), "-", shardPath)
		sortCmd.SetStdinReader(bytes.NewReader([]byte(sam)))
		sortCmd.Run()
	}

	for _, ext := range []string{"bam", "pam"} {
		outPath := fmt.Sprintf("%s/out.%s", tempDir, ext)
		var args []string
		if ext == "bam" {
			args = append([]string{"--bam", outPath}, shardPaths...)
			sh.Cmd(sort, args...).Run()
			sh.Cmd("samtools", "index", outPath).Run()
		} else {
			args = append([]string{"--pam", outPath}, shardPaths...)
			sh.Cmd(sort, args...).Run()
		}
		require.Equal(t, `@HD	VN:1.3	SO:coordinate
@SQ	SN:chr1	LN:10000
@SQ	SN:chr2	LN:9000
@SQ	SN:chr3	LN:8000
read1	0	chr1	123	60	10M	=	456	20	ACGTACGTAC	ABCDEFGHIJ	RG:Z:NA12878
read1	0	chr1	123	60	10M	=	456	20	ACGTACGTAC	ABCDEFGHIJ	RG:Z:NA12879
read2	0	chr2	111	60	10M	=	222	20	ACGTACGTAC	ABCDEFGHIJ	RG:Z:NA12878
read2	0	chr2	111	60	10M	=	222	20	ACGTACGTAC	ABCDEFGHIJ	RG:Z:NA12879
read10	12	*	0	0	10M	*	0	0	ACGTACGTAC	ABCDEFGHIJ	RG:Z:NA12878
read10	12	*	0	0	10M	*	0	0	ACGTACGTAC	ABCDEFGHIJ	RG:Z:NA12879
`, sh.Cmd("sh", "-c", fmt.Sprintf("%s view -with-header %s | grep -v @PG", pamtool, outPath)).Stdout())
	}
}
