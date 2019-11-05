package markduplicates

/**
* MIT License
*
* Copyright (c) 2017 Broad Institute
*
* Permission is hereby granted, free of charge, to any person obtaining a copy
* of this software and associated documentation files (the "Software"), to deal
* in the Software without restriction, including without limitation the rights
* to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the Software is
* furnished to do so, subject to the following conditions:
*
* The above copyright notice and this permission notice shall be included in all
* copies or substantial portions of the Software.
*
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
* IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
* FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
* AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
* LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
* OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
* SOFTWARE.
 */

import (
	"errors"
	"fmt"
	"math"

	"github.com/grailbio/base/log"
)

/**
 * Estimates the size of a library based on the number of paired end molecules observed
 * and the number of unique pairs observed.
 * Based on the Lander-Waterman equation that states:
 *   C/X = 1 - exp( -N/X )
 * where
 *   X = number of distinct molecules in library
 *   N = number of read pairs
 *   C = number of distinct fragments observed in read pairs
 */
func estimateLibrarySize(readPairs, uniqueReadPairs uint64) (uint64, error) {
	f := func(x, c, n float64) float64 {
		return c/x + math.Expm1(-n/x)
	}

	readPairDuplicates := readPairs - uniqueReadPairs
	if readPairs > 0 && readPairDuplicates > 0 {
		n := float64(readPairs)
		c := float64(uniqueReadPairs)
		m := float64(1.0)
		M := float64(100.0)

		if c >= n || f(m*c, c, n) < 0 {
			log.Fatalf("Invalid values for pairs and unique pairs: %v, %v", n, c)
		}

		// If c and n are large and almost equal, M can go to +Inf
		// before f() becomes negative.  If that happens, break out,
		// and set M to +Inf-1 to avoid looping indefinitely.  The
		// result will be meaningless, but at least this won't hang forever.
		for f(M*c, c, n) >= 0 && !math.IsInf(M, 1) {
			M *= 10.0
			if math.IsInf(M, 1) {
				return 0, fmt.Errorf("could not find M to make f() negative with arguments (%v, %v)",
					readPairs, uniqueReadPairs)
			}
		}

		for i := 0; i < 40; i++ {
			r := (m + M) / 2.0
			u := f(r*c, c, n)
			if u == 0 {
				break
			} else if u > 0 {
				m = r
			} else if u < 0 {
				M = r
			}
		}
		return uint64(c * (m + M) / 2.0), nil
	}
	return 0, errors.New("no duplicates")
}
