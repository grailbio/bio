package fastq_test

import (
	"bytes"
	"fmt"
	"reflect"
	"strings"
	"testing"

	"github.com/grailbio/bio/encoding/fastq"
)

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
	for _, test := range tests {
		var r1In, r2In, r1Out, r2Out bytes.Buffer
		for _, line := range test.r1InLines {
			r1In.WriteString(line + "\n")
		}
		for _, line := range test.r2InLines {
			r2In.WriteString(line + "\n")
		}
		err := fastq.Downsample(test.rate, &r1In, &r2In, &r1Out, &r2Out)
		if err == nil && test.err != nil {
			t.Errorf("did not get expected error: %v", test.err)
			continue
		}
		if err != nil && test.err == nil {
			t.Errorf("got unexpected error: %v", err)
			continue
		}
		if test.err == nil {
			checkDownsampleOutput(t, test.r1OutLines, &r1Out)
			checkDownsampleOutput(t, test.r2OutLines, &r2Out)
		}
	}
}

func checkDownsampleOutput(t *testing.T, expected []string, actual *bytes.Buffer) {
	actualLines := strings.Split(strings.Trim(actual.String(), "\n"), "\n")
	if actual.String() == "" {
		// We need this special case due to the behavior of strings.Split().
		actualLines = []string{}
	}
	if len(actualLines) != len(expected) {
		t.Errorf("wrong number of output lines: want %d, got %d", len(expected), len(actualLines))
	}
	if !reflect.DeepEqual(expected, actualLines) {
		t.Errorf("wrong downsample output: want %v, got %v", expected, actualLines)
	}
}
