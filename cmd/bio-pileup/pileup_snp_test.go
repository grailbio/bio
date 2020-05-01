// Copyright 2020 Grail Inc.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
package main

import (
	"bytes"
	"fmt"
	"os/exec"
	"path/filepath"
	"strconv"
	"testing"

	gunsafe "github.com/grailbio/base/unsafe"
	"github.com/grailbio/testutil"
	"github.com/grailbio/testutil/assert"
	"github.com/grailbio/testutil/expect"
)

func run(t *testing.T, bin string, args ...string) []byte {
	cmd := exec.Command(bin, args...)
	stdout := bytes.NewBuffer(nil)
	stderr := bytes.NewBuffer(nil)
	cmd.Stdout = stdout
	cmd.Stderr = stderr
	assert.NoError(t, cmd.Run(), "Command '%s %s' failed: %s", bin, args, stderr.String())
	fmt.Printf("Ran %s, stderr %s", args, stderr.String())
	return stdout.Bytes()
}

func getLinesTrunc(b []byte) [][]byte {
	// Requires b to be nonempty.  Splits it into lines with the trailing
	// newlines removed, without an empty line at the end.
	lenMinus1 := len(b) - 1
	if b[lenMinus1] == '\n' {
		b = b[:lenMinus1]
	}
	return bytes.Split(b, []byte{'\n'})
}

func TestPileupSnp(t *testing.T) {
	cmpFile := func(out string, file string) {
		testutil.CompareFile(t, out, file, nil)
	}

	executable := testutil.GoExecutable(t, "//go/src/github.com/grailbio/bio/cmd/bio-pileup/bio-pileup")

	bamPath := filepath.Join("testdata", "wgs_chr2_subset.bam")
	bedPath := filepath.Join("testdata", "common_chr2_subset.bed")
	faPath := filepath.Join("testdata", "chr2_subset.fa")

	tmpdir, cleanup := testutil.TempDir(t, "", "")
	defer testutil.NoCleanupOnError(t, cleanup, tmpdir)

	outPrefix := filepath.Join(tmpdir, "generic")
	run(t, executable, "-bed", bedPath, "-mapq=0", "-min-base-qual=42", "-stitch=true", "-cols=-lowq", "-out="+outPrefix, bamPath, faPath)

	testutil.CompareFiles(t, outPrefix+".ref.tsv", filepath.Join("testdata", "generic.ref.tsv.expected"), nil)

	altNoN := run(t, "sh", "-c", "cat "+outPrefix+".alt.tsv"+" | awk '{if ($4 != \"N\") print $0}'")
	cmpFile(string(altNoN), filepath.Join("testdata", "generic.alt.tsv.no_n.expected"))

	// Check DPs for ordered-pileup.  Values should be between 1x and 2x the
	// tiered-pileup DPs, due to lack of stitching.
	dpTiered := run(t, "cut", "-f", "4", outPrefix+".ref.tsv")
	dpTieredLines := getLinesTrunc(dpTiered)
	outPrefix = filepath.Join(tmpdir, "ordered")
	run(t, executable, "-bed", bedPath, "-mapq=0", "-cols=dpref,quals,fraglens", "-out="+outPrefix, bamPath, faPath)
	dpOrdered := run(t, "cut", "-f", "4", outPrefix+".ref.tsv")
	dpOrderedLines := getLinesTrunc(dpOrdered)
	expect.EQ(t, len(dpTieredLines), len(dpOrderedLines))
	// Skip header lines.
	dpTieredLines = dpTieredLines[1:]
	dpOrderedLines = dpOrderedLines[1:]
	for i, tieredLine := range dpTieredLines {
		tieredDp, err := strconv.Atoi(gunsafe.BytesToString(tieredLine))
		assert.NoError(t, err)
		var orderedDp int
		orderedDp, err = strconv.Atoi(gunsafe.BytesToString(dpOrderedLines[i]))
		assert.NoError(t, err)
		expect.LE(t, tieredDp, orderedDp)
		expect.LE(t, orderedDp, tieredDp*2)
	}

	// Check 5', 3', and strand columns against a manually inspected truth set.
	outPrefix = filepath.Join(tmpdir, "ordered2")
	run(t, executable, "-bed", bedPath, "-cols=enddists,strands", "-mapq=60", "-out="+outPrefix, bamPath, faPath)
	expectedRefPath := filepath.Join("testdata", "ordered2.ref.tsv.expected")
	expectedAltPath := filepath.Join("testdata", "ordered2.alt.tsv.expected")
	actualRefPath := filepath.Join(tmpdir, "ordered2.ref.tsv")
	actualAltPath := filepath.Join(tmpdir, "ordered2.alt.tsv")
	testutil.CompareFiles(t, actualRefPath, expectedRefPath, nil)
	testutil.CompareFiles(t, actualAltPath, expectedAltPath, nil)
}
