package bgzf

import (
	"bytes"
	"io/ioutil"
	"math/rand"
	"os"
	"testing"

	"github.com/grailbio/base/grail"
	"github.com/klauspost/compress/gzip"
	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
	"github.com/yasushi-saito/zlibng"
)

func TestWriter(t *testing.T) {
	// Create random bytes.
	for _, length := range []int{0, 1, 100, 65279, 65280, 65281, 500000} {
		t.Logf("length: %d", length)
		for _, useParams := range []bool{false, true} {
			input := make([]byte, length)
			n, err := rand.Read(input)
			require.Nil(t, err)
			assert.Equal(t, length, n)

			// Write bgzf
			var buf bytes.Buffer
			var w *Writer
			if useParams {
				w, err = NewWriterParams(&buf, 1, 0x0ff05, zlibng.RLEStrategy, 3)
			} else {
				w, err = NewWriter(&buf, 1)
			}
			require.Nil(t, err)
			n, err = w.Write(input)
			assert.Nil(t, err)
			assert.Equal(t, length, n)
			err = w.Close()
			assert.Nil(t, err)

			// Verify output
			if useParams && length > 0 {
				// The XFL field is set in all gzip headers, except
				// for the bgzf footer (which is a legal gzip block
				// containing zero compressed bytes).
				bufBytes := buf.Bytes()
				assert.Equal(t, byte(3), bufBytes[8], "length %d", len(bufBytes))
			}
			r, err := gzip.NewReader(&buf)
			require.Nil(t, err)
			actual, err := ioutil.ReadAll(r)
			require.Nil(t, err)
			assert.Equal(t, length, len(actual))
			assert.Equal(t, 0, bytes.Compare(input, actual))
		}
	}
}

func TestVOffset(t *testing.T) {
	// Set bgzf block size to 5.
	var buf bytes.Buffer
	w, err := NewWriterParams(&buf, 1, 5, zlibng.RLEStrategy, 0)
	require.Nil(t, err)

	// Write 4 bytes, should not cause block completion, so voffset should be (0, 4)
	_, err = w.Write([]byte("ABCD"))
	require.Nil(t, err)
	assert.Equal(t, uint64(4), w.VOffset())

	// Write 1 byte, should cause block completion, so voffset should be (non-zero, 0)
	_, err = w.Write([]byte("E"))
	require.Nil(t, err)
	voffset1 := w.VOffset()
	assert.Equal(t, uint64(0), voffset1&uint64(0xffff))
	assert.NotEqual(t, uint64(0), voffset1>>16)

	// Write 1 byte, should not cause block completion.  Coffset
	// should be the same, and uoffset should be 1.
	_, err = w.Write([]byte("F"))
	require.Nil(t, err)
	voffset2 := w.VOffset()
	assert.Equal(t, uint64(1), voffset2&uint64(0xffff))
	assert.Equal(t, voffset1>>16, voffset2>>16)
}

func TestMain(m *testing.M) {
	shutdown := grail.Init()
	defer shutdown()
	os.Exit(m.Run())
}
