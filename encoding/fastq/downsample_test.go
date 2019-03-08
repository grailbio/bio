package fastq_test

import (
	"bytes"
	"compress/gzip"
	"context"
	"fmt"
	"io/ioutil"
	"strings"
	"testing"

	"github.com/grailbio/bio/encoding/fastq"
	"github.com/grailbio/testutil"
	"github.com/grailbio/testutil/assert"
	"github.com/grailbio/testutil/expect"
)

func writeFile(t *testing.T, path string, data []string) {
	buf := bytes.Buffer{}
	gz := gzip.NewWriter(&buf)
	for _, line := range data {
		gz.Write([]byte(line + "\n"))
	}
	assert.NoError(t, gz.Close())
	assert.NoError(t, ioutil.WriteFile(path, buf.Bytes(), 0600))
}

func TestDownsample(t *testing.T) {
	tests := []struct {
		rate       float64
		r1InLines  []string
		r2InLines  []string
		r1OutLines []string
		r2OutLines []string
		err        error
	}{
		{
			1.0,
			[]string{"a", "b", "c", "d", "e", "f", "g", "h"},
			[]string{"i", "j", "k", "l", "m", "n", "o", "p"},
			[]string{"a", "b", "c", "d", "e", "f", "g", "h"},
			[]string{"i", "j", "k", "l", "m", "n", "o", "p"},
			nil,
		},
		{
			1.2,
			[]string{"a", "b", "c", "d", "e", "f", "g", "h"},
			[]string{"i", "j", "k", "l", "m", "n", "o", "p"},
			[]string{"a", "b", "c", "d", "e", "f", "g", "h"},
			[]string{"i", "j", "k", "l", "m", "n", "o", "p"},
			nil,
		},
		{
			0.0,
			[]string{"a", "b", "c", "d", "e", "f", "g", "h"},
			[]string{"i", "j", "k", "l", "m", "n", "o", "p"},
			[]string{},
			[]string{},
			nil,
		},
		{
			0.5,
			[]string{"a", "b", "c", "d", "e", "f", "g", "h"},
			[]string{"i", "j", "k", "l", "m", "n", "o", "p"},
			[]string{"e", "f", "g", "h"},
			[]string{"m", "n", "o", "p"},
			nil,
		},
		{
			1.0,
			[]string{"a", "b", "c", "d", "e", "f", "g", "h"},
			[]string{"i", "j", "k", "l"},
			nil,
			nil,
			fmt.Errorf("more reads in R1 input than in R2 input"),
		},
		{
			1.0,
			[]string{"a", "b", "c", "d"},
			[]string{"i", "j", "k", "l", "m", "n", "o", "p"},
			nil,
			nil,
			fmt.Errorf("more reads in R2 input than in R1 input"),
		},
		{
			1.0,
			[]string{"a", "b", "c", "d", "e"},
			[]string{"i", "j", "k", "l", "m", "n", "o", "p"},
			nil,
			nil,
			fmt.Errorf("error reading R1 input: too few lines in FASTQ record: want 4, got 1"),
		},
		{
			1.0,
			[]string{"a", "b", "c", "d", "e", "f", "g", "h"},
			[]string{"i", "j", "k", "l", "m", "n"},
			nil,
			nil,
			fmt.Errorf("error reading R2 input: too few lines in FASTQ record: want 4, got 2"),
		},
	}

	tempDir, cleanup := testutil.TempDir(t, "", "")
	defer cleanup()
	ctx := context.Background()
	for idx, test := range tests {
		t.Run(fmt.Sprint(idx), func(t *testing.T) {
			r1Path := fmt.Sprintf("%s/%dr1.fastq", tempDir, idx)
			r2Path := fmt.Sprintf("%s/%dr2.fastq", tempDir, idx)
			writeFile(t, r1Path, test.r1InLines)
			writeFile(t, r2Path, test.r2InLines)
			var r1Out, r2Out bytes.Buffer
			err := fastq.Downsample(ctx, test.rate, r1Path, r2Path, &r1Out, &r2Out)
			if err == nil && test.err != nil {
				t.Errorf("did not get expected error: %v", test.err)
				return
			}
			if err != nil && test.err == nil {
				t.Errorf("got unexpected error: %v", err)
				return
			}
			if test.err == nil {
				checkDownsampleOutput(t, test.r1OutLines, &r1Out)
				checkDownsampleOutput(t, test.r2OutLines, &r2Out)
			}
		})
	}
}

func TestDownsampleToCount(t *testing.T) {
	tests := []struct {
		count      int64
		r1InLines  []string
		r2InLines  []string
		r1OutLines []string
		r2OutLines []string
		err        error
	}{
		{
			2,
			[]string{"a", "b", "c", "d", "e", "f", "g", "h"},
			[]string{"i", "j", "k", "l", "m", "n", "o", "p"},
			[]string{"a", "b", "c", "d", "e", "f", "g", "h"},
			[]string{"i", "j", "k", "l", "m", "n", "o", "p"},
			nil,
		},
		{
			4,
			[]string{"a", "b", "c", "d", "e", "f", "g", "h"},
			[]string{"i", "j", "k", "l", "m", "n", "o", "p"},
			[]string{"a", "b", "c", "d", "e", "f", "g", "h"},
			[]string{"i", "j", "k", "l", "m", "n", "o", "p"},
			nil,
		},
		{
			1,
			[]string{"a", "b", "c", "d", "e", "f", "g", "h"},
			[]string{"i", "j", "k", "l", "m", "n", "o", "p"},
			[]string{"e", "f", "g", "h"},
			[]string{"m", "n", "o", "p"},
			nil,
		},
	}
	tempDir, cleanup := testutil.TempDir(t, "", "")
	defer cleanup()
	ctx := context.Background()
	for idx, test := range tests {
		t.Run(fmt.Sprint(idx), func(t *testing.T) {
			r1Path := fmt.Sprintf("%s/%dr1.fastq", tempDir, idx)
			r2Path := fmt.Sprintf("%s/%dr2.fastq", tempDir, idx)
			writeFile(t, r1Path, test.r1InLines)
			writeFile(t, r2Path, test.r2InLines)
			var r1Out, r2Out bytes.Buffer
			err := fastq.DownsampleToCount(ctx, test.count, r1Path, r2Path, &r1Out, &r2Out)
			if err == nil && test.err != nil {
				t.Errorf("did not get expected error: %v", test.err)
				return
			}
			if err != nil && test.err == nil {
				t.Errorf("got unexpected error: %v", err)
				return
			}
			if test.err == nil {
				checkDownsampleOutput(t, test.r1OutLines, &r1Out)
				checkDownsampleOutput(t, test.r2OutLines, &r2Out)
			}
		})
	}
}

func checkDownsampleOutput(t *testing.T, expected []string, actual *bytes.Buffer) {
	actualLines := strings.Split(strings.Trim(actual.String(), "\n"), "\n")
	if actual.String() == "" {
		// We need this special case due to the behavior of strings.Split().
		actualLines = []string{}
	}
	expect.EQ(t, actualLines, expected)
}

func TestDownsampleLarge(t *testing.T) {
	const nRead = 57209 // The below two files contain 57209 reads each.
	r1Path := testutil.GetFilePath("//reflow/modules/testdata/af4/test_1.fastq.gz")
	r2Path := testutil.GetFilePath("//reflow/modules/testdata/af4/test_2.fastq.gz")
	for _, count := range []int64{100, 1000, 10000} {
		t.Run("count-"+fmt.Sprint(count), func(t *testing.T) {
			t.Parallel()
			var r1Out, r2Out bytes.Buffer
			assert.NoError(t, fastq.DownsampleToCount(context.Background(), count, r1Path, r2Path, &r1Out, &r2Out))
			nLine1 := bytes.Count(r1Out.Bytes(), []byte("\n")) / 4
			nLine2 := bytes.Count(r2Out.Bytes(), []byte("\n")) / 4
			expect.EQ(t, nLine1, nLine2)
			expect.GE(t, nLine1, int(float64(count)*0.9))
			expect.LE(t, nLine1, int(float64(count)*1.1))
		})
	}

	for _, rate := range []float64{0.001, 0.01, 0.1} {
		t.Run("rate-"+fmt.Sprint(rate), func(t *testing.T) {
			t.Parallel()
			var r1Out, r2Out bytes.Buffer
			assert.NoError(t, fastq.Downsample(context.Background(), rate, r1Path, r2Path, &r1Out, &r2Out))
			nLine1 := bytes.Count(r1Out.Bytes(), []byte("\n")) / 4
			nLine2 := bytes.Count(r2Out.Bytes(), []byte("\n")) / 4
			expect.EQ(t, nLine1, nLine2)
			expect.GE(t, nLine1, int(nRead*rate*0.9))
			expect.LE(t, nLine1, int(nRead*rate*1.1))
		})
	}
}
