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
package snp

import (
	"testing"
)

func TestQualPassTable(t *testing.T) {
	minBaseQual := byte(43)
	qpt, _ := newQualPassTable(minBaseQual)

	// Every row and column should be nonincreasing, and we spot-check a few
	// positions.
	for row := byte(0); row < nQual; row++ {
		prevRowPass := false
		for col := byte(0); col < nQual; col++ {
			curRowPass := qpt.lookup2(row, col)
			if (!curRowPass) && prevRowPass {
				t.Fatalf("qual-pass table.lookup2(%d, %d) is false, when previous value in same row was true", row, col)
			}
			prevRowPass = curRowPass
		}
	}
	for col := byte(0); col < nQual; col++ {
		prevColPass := false
		for row := byte(0); row < nQual; row++ {
			curColPass := qpt.lookup2(row, col)
			if (!curColPass) && prevColPass {
				t.Fatalf("qual-pass table.lookup2(%d, %d) is false, when previous value in same row was true", row, col)
			}
			prevColPass = curColPass
		}
	}
	if qpt.lookup2(15, 15) != false {
		t.Fatalf("Unexpected value %v for qpt.lookup2(15, 15) (expected false)", qpt.lookup2(15, 15))
	}
	if qpt.lookup2(21, 21) != false {
		t.Fatalf("Unexpected value %v for qpt.lookup2(21, 21) (expected false)", qpt.lookup2(21, 21))
	}
	if qpt.lookup2(25, 25) != true {
		t.Fatalf("Unexpected value %v for qpt.lookup2(25, 25) (expected true)", qpt.lookup2(25, 25))
	}
}
