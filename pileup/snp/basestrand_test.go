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
package snp_test

import (
	"bytes"
	"os"
	"path"
	"path/filepath"
	"reflect"
	"strings"
	"sync"
	"testing"

	"github.com/grailbio/base/tsv"
	"github.com/grailbio/bio/pileup/snp"
	"github.com/grailbio/testutil"
	"github.com/grailbio/testutil/assert"
)

const testDataPath = "//go/src/github.com/grailbio/bio/pileup/snp/testdata/"

// TestReadBaseStrandTsv_fromFile tests reading a raw base strand file with ReadBaseStrandTsv.
func TestReadBaseStrandTsv_fromFile(t *testing.T) {
	filePath := testutil.GetFilePath(filepath.Join(testDataPath, "sample1.basestrand.tsv"))
	file, err := os.Open(filePath)
	assert.NoError(t, err)

	rows, err := snp.ReadBaseStrandTsv(file)
	assert.NoError(t, err)

	expectedRows := []snp.BaseStrandTsvRow{
		{
			Chr:  "chr8",
			Pos:  50824362,
			Ref:  "C",
			FwdA: 0,
			RevA: 0,
			FwdC: 0,
			RevC: 2,
			FwdG: 0,
			RevG: 0,
			FwdT: 0,
			RevT: 0,
		},
		{
			Chr:  "chr8",
			Pos:  50824381,
			Ref:  "T",
			FwdA: 0,
			RevA: 0,
			FwdC: 0,
			RevC: 0,
			FwdG: 0,
			RevG: 0,
			FwdT: 0,
			RevT: 7,
		},
	}
	assert.EQ(t, expectedRows, rows)
}

// TestReadBaseStrandTsv tests ReadBaseStrandTsv.
func TestReadBaseStrandTsv(t *testing.T) {
	for _, test := range []struct {
		name          string
		input         [][]string
		expectedRows  []snp.BaseStrandTsvRow
		errorExpected bool
		expectedError string
	}{
		{
			"oneRow",
			[][]string{{"chr1", "569378", "C", "0", "0", "0", "23", "0", "0", "0", "0"}},
			[]snp.BaseStrandTsvRow{
				{
					Chr:  "chr1",
					Pos:  569378,
					Ref:  "C",
					FwdA: 0,
					RevA: 0,
					FwdC: 0,
					RevC: 23,
					FwdG: 0,
					RevG: 0,
					FwdT: 0,
					RevT: 0,
				},
			},
			false,
			"",
		}, {
			"twoRows",
			[][]string{
				{"chr1", "569378", "C", "0", "0", "0", "23", "0", "0", "0", "0"},
				{"chr1", "569381", "A", "0", "23", "0", "0", "0", "0", "0", "0"},
			},
			[]snp.BaseStrandTsvRow{
				{
					Chr:  "chr1",
					Pos:  569378,
					Ref:  "C",
					FwdA: 0,
					RevA: 0,
					FwdC: 0,
					RevC: 23,
					FwdG: 0,
					RevG: 0,
					FwdT: 0,
					RevT: 0,
				},
				{
					Chr:  "chr1",
					Pos:  569381,
					Ref:  "A",
					FwdA: 0,
					RevA: 23,
					FwdC: 0,
					RevC: 0,
					FwdG: 0,
					RevG: 0,
					FwdT: 0,
					RevT: 0,
				},
			},
			false,
			"",
		},
		{
			"invalidRow_missingColumn",
			[][]string{{"chr1", "569378", "C", "A", "0", "0", "23", "0", "0", "0", "0"}},
			[]snp.BaseStrandTsvRow{},
			true,
			"invalid syntax",
		},
	} {
		reader := bytes.NewBufferString(makeTestTsv(test.input))

		rows, err := snp.ReadBaseStrandTsv(reader)

		if test.errorExpected {
			assert.HasSubstr(t, err.Error(), test.expectedError)
		} else {
			assert.EQ(t, test.expectedRows, rows)
		}
	}
}

func TestReadSingleStrandBaseStrandTsv(t *testing.T) {
	for _, test := range []struct {
		name         string
		fwdInput     [][]string
		revInput     [][]string
		expectedRows []snp.BaseStrandTsvRow
	}{
		{
			"noRows",
			[][]string{{"CHROM", "POS", "DEPTH", "REF", "A", "C", "G", "T", "N", "INS", "DEL"}},
			[][]string{{"CHROM", "POS", "DEPTH", "REF", "A", "C", "G", "T", "N", "INS", "DEL"}},
			[]snp.BaseStrandTsvRow{},
		},
		{
			"multipleRows",
			[][]string{
				{"CHROM", "POS", "DEPTH", "REF", "A", "C", "G", "T", "N", "INS", "DEL"},
				{"chr1", "569378", "5", "A", "5", "0", "0", "0", "0", "0", "0"},
				{"chr1", "569381", "8", "G", "0", "0", "7", "1", "0", "0", "0"},
				{"chr1", "569384", "7", "G", "0", "0", "7", "0", "0", "0", "0"},
			},
			[][]string{
				{"CHROM", "POS", "DEPTH", "REF", "A", "C", "G", "T", "N", "INS", "DEL"},
				{"chr1", "569378", "8", "A", "8", "0", "0", "0", "0", "0", "0"},
				{"chr1", "569381", "7", "G", "0", "0", "6", "1", "0", "0", "0"},
				{"chr1", "569384", "7", "G", "0", "0", "7", "0", "0", "0", "0"},
			},
			[]snp.BaseStrandTsvRow{
				{
					Chr:  "chr1",
					Pos:  569378,
					Ref:  "A",
					FwdA: 5,
					RevA: 8,
					FwdC: 0,
					RevC: 0,
					FwdG: 0,
					RevG: 0,
					FwdT: 0,
					RevT: 0,
				},
				{
					Chr:  "chr1",
					Pos:  569381,
					Ref:  "G",
					FwdA: 0,
					RevA: 0,
					FwdC: 0,
					RevC: 0,
					FwdG: 7,
					RevG: 6,
					FwdT: 1,
					RevT: 1,
				},
				{
					Chr:  "chr1",
					Pos:  569384,
					Ref:  "G",
					FwdA: 0,
					RevA: 0,
					FwdC: 0,
					RevC: 0,
					FwdG: 7,
					RevG: 7,
					FwdT: 0,
					RevT: 0,
				},
			},
		},
	} {
		t.Run(test.name, func(t *testing.T) {
			fwdReader := bytes.NewBufferString(makeTestTsv(test.fwdInput))
			revReader := bytes.NewBufferString(makeTestTsv(test.revInput))

			rows, err := snp.ReadSingleStrandBaseStrandTsv(fwdReader, revReader)
			assert.NoError(t, err)
			assert.EQ(t, test.expectedRows, rows)
		})
	}
}

// makeTestInput makes a test file tsv by joining rows with "/t" and lines with "/n".
func makeTestTsv(data [][]string) string {
	testString := ""
	for _, row := range data {
		testString += strings.Join(row, "\t") + "\n"
	}
	return testString
}

func TestReadBaseStrandTsvIntoChannel(t *testing.T) {
	var pileupFile = testutil.GetFilePath(path.Join(testDataPath, "test_pileup_snp.tsv"))
	file, err := os.Open(pileupFile)
	assert.NoError(t, err)
	reader := tsv.NewReader(file)
	reader.HasHeaderRow = true
	reader.UseHeaderNames = true

	assert.NoError(t, err)

	channel := make(chan []snp.BaseStrandTsvRow)
	var waitGroup sync.WaitGroup
	go snp.ReadBaseStrandTsvIntoChannel(reader, channel, 100, pileupFile, &waitGroup)
	waitGroup.Wait()
	var pileupSlice []snp.BaseStrandTsvRow

	for s := range channel {
		pileupSlice = append(pileupSlice, s...)
	}

	assert.EQ(t,
		snp.BaseStrandTsvRow{Chr: "chr1", Pos: 713277, Ref: "C", FwdA: 0,
			FwdC: 0, FwdG: 0, FwdT: 8, RevA: 0, RevC: 12, RevG: 0, RevT: 0}, pileupSlice[30])
	assert.EQ(t, 999, len(pileupSlice))
	assert.NoError(t, file.Close())
}

func TestChrId(t *testing.T) {
	for _, test := range []struct {
		name           string
		input          string
		expectedResult int
	}{
		{
			"autosomal",
			"chr12",
			12,
		},
		{
			"sexChromosome",
			"chrX",
			23,
		},
	} {
		t.Run(test.name, func(t *testing.T) {
			result, err := snp.ChrId(test.input)
			assert.NoError(t, err)
			assert.EQ(t, test.expectedResult, result)
		})
	}
}

func TestWriteBaseStrandToTSV(t *testing.T) {
	tests := []struct {
		piles    []snp.BaseStrandPile
		refNames []string
		expected string
	}{
		{
			piles: []snp.BaseStrandPile{
				{
					RefID:  0,
					Pos:    123,
					Counts: [4][2]uint32{{1, 2}, {3, 4}, {5, 6}, {7, 8}},
				},
				{
					RefID:  1,
					Pos:    456,
					Counts: [4][2]uint32{{1, 1}, {2, 2}, {3, 3}, {4, 4}},
				},
			},
			refNames: []string{
				"chr1",
				"chr2",
			},
			expected: "CHROM\tPOS\tA+\tA-\tC+\tC-\tG+\tG-\tT+\tT-\n" +
				"chr1\t123\t1\t2\t3\t4\t5\t6\t7\t8\n" +
				"chr2\t456\t1\t1\t2\t2\t3\t3\t4\t4\n",
		},
	}

	for _, test := range tests {
		var buffer bytes.Buffer
		if err := snp.WriteBaseStrandToTSV(test.piles, test.refNames, &buffer); err != nil {
			t.Fatalf("Unexpected error: %v", err)
		}
		if test.expected != buffer.String() {
			t.Errorf("Wrong output: want %s, got %s", test.expected, buffer.String())
		}
	}
}

func TestWriteBaseStrandTsv(t *testing.T) {
	tests := []struct {
		basestrands []snp.BaseStrandTsvRow
		expected    string
	}{
		{
			basestrands: []snp.BaseStrandTsvRow{
				{"chr1", 1, "A", 1, 2, 3, 4, 5, 6, 7, 8},
			},
			expected: "#CHROM\tPOS\tREF\tA+\tA-\tC+\tC-\tG+\tG-\tT+\tT-\n" +
				"chr1\t1\tA\t1\t2\t3\t4\t5\t6\t7\t8\n",
		},
	}

	for _, test := range tests {
		var buffer bytes.Buffer
		if err := snp.WriteBaseStrandTsv(test.basestrands, &buffer); err != nil {
			t.Fatalf("Unexpected error: %v", err)
		}
		if test.expected != buffer.String() {
			t.Errorf("Wrong output: want %s, got %s", test.expected, buffer.String())
		}
	}
}

func TestReadWriteBaseStrandsRio(t *testing.T) {
	tests := []struct {
		piles    []snp.BaseStrandPile
		refNames []string
	}{
		{
			piles: []snp.BaseStrandPile{
				{
					RefID:  0,
					Pos:    123,
					Counts: [4][2]uint32{{0, 1}, {2, 3}, {4, 5}, {6, 7}},
				},
				{
					RefID:  1,
					Pos:    456,
					Counts: [4][2]uint32{{0, 1}, {2, 2}, {4, 4}, {8, 8}},
				},
				{
					RefID:  1,
					Pos:    777,
					Counts: [4][2]uint32{{9, 8}, {7, 6}, {5, 4}, {3, 2}},
				},
			},
			refNames: []string{"chr1", "chr2"},
		},
	}

	for _, test := range tests {
		var buffer bytes.Buffer
		err := snp.WriteBaseStrandsRio(test.piles, test.refNames, &buffer)
		assert.NoError(t, err)
		piles, refNames, err := snp.ReadBaseStrandsRio(bytes.NewReader(buffer.Bytes()))
		assert.NoError(t, err)
		if !reflect.DeepEqual(test.piles, piles) {
			t.Errorf("Wrong piles: want %v, got %v", test.piles, piles)
		}
		if !reflect.DeepEqual(test.refNames, refNames) {
			t.Errorf("Wrong ref names: want %v, got %v", test.refNames, refNames)
		}
	}
}
