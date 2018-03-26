package bam_test

import (
	"os"
	"testing"

	"bytes"
	"encoding/binary"
	"github.com/biogo/hts/bam"
	grailbam "github.com/grailbio/bio/encoding/bam"
	"github.com/grailbio/internal/testutil"
	"github.com/stretchr/testify/require"
	"io"
)

func TestMarshal(t *testing.T) {
	path := testutil.GetFilePath("@grailgo//bio/encoding/bam/testdata/170614_WGS_LOD_Pre_Library_B3_27961B_05.merged.10000.bam")
	in, err := os.Open(path)
	require.NoErrorf(t, err, "path: %s", path)
	r, err := bam.NewReader(in, 0)
	require.NoErrorf(t, err, "path: %s", path)
	for {
		rec, err := r.Read()
		if err == io.EOF {
			break
		}
		require.NoError(t, err)

		buf := bytes.NewBuffer(nil)
		require.NoError(t, grailbam.Marshal(rec, buf))
		serialized := buf.Bytes()
		serializedLen := int(binary.LittleEndian.Uint32(serialized[:4]))
		require.Equal(t, serializedLen, len(serialized)-4)

		rec2, err := grailbam.Unmarshal(serialized[4:], r.Header())
		require.NoError(t, err, "rec=", rec.String())
		require.Equal(t, rec.String(), rec2.String())
	}
	require.NoError(t, r.Close())
	require.NoError(t, in.Close())
}
