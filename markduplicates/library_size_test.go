package markduplicates

import (
	"testing"

	"github.com/stretchr/testify/assert"
)

func TestEstimateLibrarySize(t *testing.T) {
	tests := []struct {
		readPairs       uint64
		uniqueReadPairs uint64
		expected        uint64
	}{
		{1000000, 800000, 2154184},
		{171512300, 171512299, 14708234445116054},
	}

	for _, test := range tests {
		v, err := estimateLibrarySize(test.readPairs, test.uniqueReadPairs)
		assert.NoError(t, err)
		assert.InEpsilon(t, test.expected, v, 0.0000000001)
	}
}
