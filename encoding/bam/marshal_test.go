// Copyright 2018 GRAIL, Inc. All rights reserved.
// Use of this source code is governed by the Apache 2.0
// license that can be found in the LICENSE file.

package bam_test

import (
	"bytes"
	"encoding/binary"
	"io"
	"os"
	"testing"

	grailbam "github.com/grailbio/bio/encoding/bam"
	"github.com/grailbio/hts/bam"
	"github.com/grailbio/testutil"
	"github.com/grailbio/testutil/assert"
)

func TestMarshal(t *testing.T) {
	path := testutil.GetFilePath("//go/src/grail.com/bio/encoding/bam/testdata/170614_WGS_LOD_Pre_Library_B3_27961B_05.merged.10000.bam")
	in, err := os.Open(path)
	assert.NoError(t, err, "path: %s", path)
	r, err := bam.NewReader(in, 0)
	assert.NoError(t, err, "path: %s", path)
	for {
		rec, err := r.Read()
		if err == io.EOF {
			break
		}
		assert.NoError(t, err)

		buf := bytes.NewBuffer(nil)
		assert.NoError(t, bam.Marshal(rec, buf))
		serialized := buf.Bytes()
		serializedLen := int(binary.LittleEndian.Uint32(serialized[:4]))
		assert.EQ(t, len(serialized)-4, serializedLen)

		rec2, err := grailbam.Unmarshal(serialized[4:], r.Header())
		assert.NoError(t, err, "rec=", rec.String())
		assert.EQ(t, rec2.String(), rec.String())
	}
	assert.NoError(t, r.Close())
	assert.NoError(t, in.Close())
}

func BenchmarkRead(b *testing.B) {
	path := testutil.GetFilePath("//go/src/grail.com/bio/encoding/bam/testdata/170614_WGS_LOD_Pre_Library_B3_27961B_05.merged.10000.bam")
	for i := 0; i < b.N; i++ {
		in, err := os.Open(path)
		assert.NoError(b, err, "path: %s", path)
		r, err := bam.NewReader(in, 0)
		assert.NoError(b, err, "path: %s", path)
		for {
			if _, err := r.Read(); err == io.EOF {
				break
			}
		}
		assert.NoError(b, r.Close())
		assert.NoError(b, in.Close())
	}
}
