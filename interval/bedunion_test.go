package interval

import (
	"math"
	"reflect"
	"testing"

	"github.com/grailbio/testutil/expect"
	_ "grail.com/cloud/grailfile"
)

func TestLoadSortedBEDIntervals(t *testing.T) {
	tests := []struct {
		pathname              string
		invert, oneBasedInput bool
		want                  BEDUnion
	}{
		{"testdata/test1.bed",
			false,
			false,
			BEDUnion{
				nameMap: map[string]([]PosType){
					"chr1": []PosType{
						2488104, 2488172,
						2489165, 2489273,
						2489782, 2489907,
						2490320, 2490438,
						2491262, 2491417,
						2492063, 2492157,
						2493112, 2493254,
						2494304, 2494335,
						2494587, 2494712},
				},
				lastRefName: "",
				lastRefID:   -1,
			},
		},
		{"testdata/test2.bed",
			true,
			true,
			BEDUnion{
				nameMap: map[string]([]PosType){
					"chr1": []PosType{
						-1,
						2488103, 2488172,
						2489164, 2489273,
						2489781, 2489907,
						math.MaxInt32},
					"chr2": []PosType{
						-1,
						2490319, 2492157,
						2493111, 2494254,
						2494586, 2494712,
						math.MaxInt32},
				},
				lastRefName: "",
				lastRefID:   -1,
			},
		},
	}

	for _, tt := range tests {
		result, err := NewBEDUnionFromPath(
			tt.pathname,
			NewBEDOpts{
				Invert:        tt.invert,
				OneBasedInput: tt.oneBasedInput,
			},
		)
		expect.NoError(t, err)
		if !reflect.DeepEqual(result, tt.want) {
			t.Errorf("Wanted: %v  Got: %v", tt.want, result)
		}
	}
}

func TestParseRegionString(t *testing.T) {
	tests := []struct {
		region  string
		chrName string
		start0  PosType
		end     PosType
	}{
		{
			"chr1:1-1000",
			"chr1",
			0,
			1000,
		},
		{
			"chr1:1000",
			"chr1",
			999,
			1000,
		},
		{
			"chr1",
			"chr1",
			0,
			math.MaxInt32 - 1,
		},
	}

	for _, tt := range tests {
		result, err := ParseRegionString(tt.region)
		expect.NoError(t, err)
		expect.EQ(t, tt.chrName, result.RefName)
		expect.EQ(t, tt.start0, result.Start0)
		expect.EQ(t, tt.end, result.End)
	}
}
