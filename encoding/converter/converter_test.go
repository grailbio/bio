package converter_test

import (
	"math"
	"path/filepath"
	"testing"

	gbam "github.com/grailbio/bio/encoding/bam"
	"github.com/grailbio/bio/encoding/bamprovider"
	"github.com/grailbio/bio/encoding/converter"
	"github.com/grailbio/bio/encoding/pam"
	"github.com/grailbio/hts/sam"
	"github.com/grailbio/testutil"
	"github.com/grailbio/testutil/assert"
	"v.io/x/lib/gosh"
	"v.io/x/lib/lookpath"
)

func TestPAM(t *testing.T) {
	tempDir, cleanup := testutil.TempDir(t, "", "")
	defer cleanup()

	bamPath := testutil.GetFilePath("//go/src/grail.com/bio/encoding/bam/testdata/170614_WGS_LOD_Pre_Library_B3_27961B_05.merged.10000.bam")
	pamPath := filepath.Join(tempDir, "test.pam")
	assert.NoError(t, converter.ConvertToPAM(pam.WriteOpts{}, pamPath, bamPath, "", math.MaxInt64))

	verifyFiles(t, bamPath, pamPath)
}

func TestBAM(t *testing.T) {
	sh := gosh.NewShell(t)
	defer sh.Cleanup()
	if !hasSamtools(t, sh) {
		return
	}
	tempDir, cleanup := testutil.TempDir(t, "", "")
	defer cleanup()

	bamPath := testutil.GetFilePath("//go/src/grail.com/bio/encoding/bam/testdata/170614_WGS_LOD_Pre_Library_B3_27961B_05.merged.10000.bam")
	pamPath := filepath.Join(tempDir, "test.pam")
	assert.NoError(t, converter.ConvertToPAM(pam.WriteOpts{}, pamPath, bamPath, "", math.MaxInt64))

	bam2Path := filepath.Join(tempDir, "test.bam")
	convertToBAM(t, sh, pamPath, bam2Path)
	verifyFiles(t, pamPath, bam2Path)
}

func hasSamtools(t *testing.T, sh *gosh.Shell) bool {
	if _, err := lookpath.Look(sh.Vars, "samtools"); err != nil {
		t.Skipf("samtools not found on the machine. Skipping the test")
		return false
	}
	return true
}

func convertToBAM(t *testing.T, sh *gosh.Shell, srcPath, destPath string) {
	p := bamprovider.NewProvider(srcPath)
	assert.NoError(t, converter.ConvertToBAM(destPath, p))
	assert.NoError(t, p.Close())

	cmd := sh.Cmd("samtools", "index", destPath)
	cmd.Run()
	assert.NoError(t, cmd.Err)
}

// verifyFiles verifes that files path0 and path1 store the same records in the
// same order.
func verifyFiles(t *testing.T, path0, path1 string) {
	p0 := bamprovider.NewProvider(path0)
	header0, err := p0.GetHeader()
	assert.NoError(t, err)
	p1 := bamprovider.NewProvider(path1)
	header1, err := p0.GetHeader()
	assert.NoError(t, err)

	iter0 := p0.NewIterator(gbam.UniversalShard(header0))
	iter1 := p1.NewIterator(gbam.UniversalShard(header1))
	n := 0
	for {
		ok0 := iter0.Scan()
		ok1 := iter1.Scan()
		assert.EQ(t, ok1, ok0, "path0=%s, path1=%s, n=%d", path0, path1, n)
		if !ok0 {
			break
		}
		r0 := iter0.Record()
		r1 := iter1.Record()
		assert.EQ(t, r1.String(), r0.String(), "path0=%s, path1=%s, n=%d", path0, path1, n)
		n++
		sam.PutInFreePool(r0)
		sam.PutInFreePool(r1)
	}
	assert.NoError(t, iter0.Close())
	assert.NoError(t, iter1.Close())
	assert.NoError(t, p0.Close())
	assert.NoError(t, p1.Close())
}
