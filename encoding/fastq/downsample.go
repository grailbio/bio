package fastq

import (
	"bufio"
	"bytes"
	"io"
	"math/rand"

	"github.com/pkg/errors"
)

const (
	linesPerRead = 4
)

// Downsample writes read pairs from r1In and r2In to r1Out and r2Out. Read pairs will be randomly
// selected for inclusion in the output at the given sampling rate.
func Downsample(rate float64, r1In, r2In io.Reader, r1Out, r2Out io.Writer) error {
	if rate < 0.0 || rate > 1.0 {
		return errors.New("rate must be between 0 and 1 (inclusive)")
	}
	random := rand.New(rand.NewSource(0))
	r1Scanner := bufio.NewScanner(r1In)
	r2Scanner := bufio.NewScanner(r2In)
	for {
		r1, r1Err := scanRead(r1Scanner)
		if r1Err != nil && r1Err != io.EOF {
			return errors.Wrap(r1Err, "error reading R1 input")
		}
		r2, r2Err := scanRead(r2Scanner)
		if r2Err != nil && r2Err != io.EOF {
			return errors.Wrap(r2Err, "error reading R2 input")
		}
		if r1Err == io.EOF && r2Err == io.EOF {
			// Both readers ended after the same number of reads, as expected.
			return nil
		} else if r1Err == io.EOF {
			return errors.New("more reads in R2 input than in R1 input")
		} else if r2Err == io.EOF {
			return errors.New("more reads in R1 input than in R2 input")
		}
		if random.Float64() < rate {
			r1Out.Write(r1)
			r2Out.Write(r2)
		}
	}
	return nil
}

func scanRead(scanner *bufio.Scanner) ([]byte, error) {
	var buffer bytes.Buffer
	for i := 0; i < linesPerRead; i++ {
		if !scanner.Scan() {
			if i == 0 && scanner.Err() == nil {
				// Reached end of input.
				return nil, io.EOF
			}
			// Something went wrong.
			if scanner.Err() != nil {
				return nil, scanner.Err()
			}
			return nil, errors.Errorf("too few lines in FASTQ record: want %d, got %d", linesPerRead, i)
		}
		buffer.WriteString(scanner.Text())
		buffer.WriteString("\n")
	}
	return buffer.Bytes(), nil
}
