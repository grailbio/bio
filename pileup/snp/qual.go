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
	"fmt"
	"math"
)

// This file contains qual phred-math routines.  Some of them may migrate to
// more central locations in the future.

// All functions here assume input qual scores are never larger than
// (nQual - 1), and never return qual scores larger than that.
const nQual = 96

// When stitching two identical bases, we can estimate the error probability as
//   (e1 * e2) / ((1 - e1) * (1 - e2) + e1 * e2).
// (This conservatively assumes that some sequencing errors can be far more
// likely than others; we can divide by ~3 if we instead assume errors are
// evenly distributed among the other three possibilities.)
// qualSumTable stores the phred scores corresponding to these precomputed
// values.
var qualSumTable [nQual][nQual]byte

func init() {
	// Fortunately, nQual is small enough that we don't have to worry about
	// floating-point underflow anywhere.
	var errProbs [nQual]float64
	for i := range errProbs {
		errProbs[i] = math.Exp(float64(i) * (-0.1 * math.Ln10))
	}
	for i := range qualSumTable {
		e1 := errProbs[i]
		for j := 0; j <= i; j++ {
			e2 := errProbs[j]
			errProduct := e1 * e2
			newErrProb := errProduct / ((1.0-e1)*(1.0-e2) + errProduct)
			curQual := byte(math.Round(math.Log(newErrProb) * (-10.0 * math.Log10E)))
			if curQual >= nQual {
				curQual = nQual - 1
			}
			qualSumTable[i][j] = curQual
			qualSumTable[j][i] = curQual
		}
	}
}

type qualPassTable [nQual][nQual]bool

func newQualPassTable(minBaseQual byte) (t qualPassTable, err error) {
	if minBaseQual >= nQual {
		err = fmt.Errorf("newQualPassTable: minBaseQual too large")
		return
	}
	for i := range t {
		for j := 0; j <= i; j++ {
			curQual := qualSumTable[i][j]
			if curQual >= minBaseQual {
				t[i][j] = true
				t[j][i] = true
			}
		}
	}
	return
}

func (t *qualPassTable) lookup2(q1, q2 byte) bool {
	return t[q1][q2]
}
