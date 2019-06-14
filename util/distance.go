package util

import (
	"fmt"
	"strconv"
	"strings"
)

// matrix represents a 2 dimensional matrix.
type matrix struct {
	nRow, nCol int
	data       []int // row-major nRow*nCol array.
}

// matrix returns an n x m matrix.
func newMatrix(n, m int) (x matrix) {
	return matrix{
		nRow: n,
		nCol: m,
		data: make([]int, n*m),
	}
}

// String returns a string representation of a matrix.
// TODO(ayip): this could be implemented using text/tabwriter.
func (m matrix) String() (r string) {
	maxLength := 0
	for _, d := range m.data {
		if l := len(strconv.Itoa(d)); l > maxLength {
			maxLength = l
		}
	}

	lines := []string{"\n"}
	for i := 0; i < m.nRow; i++ {
		var parts []string
		for j := 0; j < m.nCol; j++ {
			parts = append(parts, fmt.Sprintf("%0*s", maxLength, strconv.Itoa(m.data[i*m.nCol+j])))
		}
		lines = append(lines, strings.Join(parts, " | "))
	}
	return strings.Join(lines, "\n")
}

// operation is a type that describes one of the three possible traversals in a
// Levenshtein edit distance matrix.
//
//   ___|___
//    1 | 3
//    2 | 4
//
// (1) diagonal (1 -> 4)
// (2) right (2 -> 4)
// (3) down (3 -> 4)
type operation uint8

// diagonal, right, and down refer to the three possible traversals allowed in
// the Levenshtein edit distance matrix.
const (
	diagonal operation = iota
	right
	down
)

// operations is a slice of operation types.
type operations []operation

// contains checks whether the slice contains any operations in a given operation slice.
func (o operations) contains(given operations) bool {
	for _, x := range given {
		for _, y := range o {
			if x == y {
				return true
			}
		}
	}
	return false
}

// computeRow computes cells in a Levenshtein matrix for a given row specified
// by i up to the column specified by 'col'.
func (m matrix) computeRow(i, col int, r1, r2 []byte) {
	for j := 0; j <= col; j++ {
		m.computeCell(i, j, r1, r2)
	}
}

// computeCol computes cells in a Levenshtein matrix for a given column
// specified by j up to the row specified by 'row'.
func (m matrix) computeCol(j, row int, r1, r2 []byte) {
	for i := 0; i <= row; i++ {
		m.computeCell(i, j, r1, r2)
	}
}

// computeCell computes the cell (i, j) in a Levenshtein matrix.
func (m matrix) computeCell(i, j int, r1, r2 []byte) operations {
	if i == 0 {
		m.data[i*m.nCol+j] = j
		return operations{}
	}
	if j == 0 {
		m.data[i*m.nCol+j] = i
		return operations{}
	}
	if r1[i-1] == r2[j-1] {
		m.data[i*m.nCol+j] = m.data[(i-1)*m.nCol+(j-1)]
		return operations{diagonal}
	}

	downValue := m.data[(i-1)*m.nCol+j] + 1
	diagonalValue := m.data[(i-1)*m.nCol+(j-1)] + 1
	rightValue := m.data[i*m.nCol+(j-1)] + 1

	minValue := downValue
	if diagonalValue < minValue {
		minValue = diagonalValue
	}
	if rightValue < minValue {
		minValue = rightValue
	}

	m.data[i*m.nCol+j] = minValue

	// Identify the operation(s) that led to the computed
	// minimum value.
	r := operations{}
	if downValue == minValue {
		r = append(r, down)
	}
	if diagonalValue == minValue {
		r = append(r, diagonal)
	}
	if rightValue == minValue {
		r = append(r, right)
	}
	return r
}

// Levenshtein computes the Levenshtein distance between two barcodes. The
// returned value - distance - is the number of insertions, deletions, and
// substitutions it takes to transform one barcode (s1) into another (s2). Each
// step in the transformation "costs" one distance point. Because a fixed
// number of barcode bases are always sequenced, bases downstream of the
// barcode will be read in the event of a deletion in the barcode sequence. To
// account for this situation, we take in the sequence downstream of both
// barcodes (a1 and a2).  Note that s1 and s2 must have the same length.
//
// TODO(ayip): we can optimize this for memory allocations by creating
// a reusable object that contains the working state for each
// invocation of Levenshtein().
func Levenshtein(s1, s2, a1, a2 string) (distance int) {
	if len(s1) != len(s2) {
		panic(fmt.Sprintf("s1 and s2 must have equal length: '%s', '%s'", s1, s2))
	}

	r1 := []byte(s1)
	r2 := []byte(s2)

	rows := len(r1)
	cols := len(r2)

	m := newMatrix(rows+len(a1)+1, cols+len(a2)+1)

	i := 1
	iEnd := rows
	j := 1
	jEnd := cols

	var cellOperations operations
	for {
		if i <= iEnd {
			m.computeRow(i, j-1, r1, r2)
		}

		if j <= jEnd {
			m.computeCol(j, i-1, r1, r2)
		}

		cellOperations = m.computeCell(i, j, r1, r2)

		if i < rows {
			i++
			j++
			continue
		}

		if i >= rows {
			done := true
			if cellOperations.contains(operations{down}) && len(a2) > 0 {
				r2 = append(r2, a2[0])
				a2 = a2[1:]
				done = false
				j++
				jEnd++
			}
			if cellOperations.contains(operations{right}) && len(a1) > 0 {
				r1 = append(r1, a1[0])
				a1 = a1[1:]
				done = false
				i++
				iEnd++
			}
			if done {

				if m.data[rows*m.nCol+cols] <= m.data[i*m.nCol+j] {
					return m.data[rows*m.nCol+cols]
				}
				return m.data[i*m.nCol+j]
			}
		}
	}
}
