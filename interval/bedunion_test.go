package interval

import (
	"math"
	"reflect"
	"testing"

	"github.com/grailbio/hts/sam"
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
		refName string
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
		expect.EQ(t, tt.refName, result.RefName)
		expect.EQ(t, tt.start0, result.Start0)
		expect.EQ(t, tt.end, result.End)
	}
}

type IntersectsByIDTestInstance struct {
	refID    int
	startPos PosType
	limitPos PosType
	want     bool
}

func TestIntersectsByID(t *testing.T) {
	tests := []struct {
		pathname string
		insts    []IntersectsByIDTestInstance
	}{
		{
			pathname: "testdata/test1.bed",
			insts: []IntersectsByIDTestInstance{
				{
					startPos: 2488103,
					limitPos: 2488104,
					want:     false,
				},
				{
					startPos: 2488103,
					limitPos: 2488105,
					want:     true,
				},
				{
					startPos: 2488103,
					limitPos: 2488171,
					want:     true,
				},
				{
					startPos: 2488103,
					limitPos: 2488172,
					want:     true,
				},
				{
					startPos: 2488103,
					limitPos: 2488173,
					want:     true,
				},
				{
					startPos: 2488104,
					limitPos: 2488105,
					want:     true,
				},
				{
					startPos: 2488104,
					limitPos: 2488171,
					want:     true,
				},
				{
					startPos: 2488104,
					limitPos: 2488172,
					want:     true,
				},
				{
					startPos: 2488104,
					limitPos: 2488173,
					want:     true,
				},
				{
					startPos: 2488172,
					limitPos: 2488173,
					want:     false,
				},
				{
					startPos: 2488105,
					limitPos: 2488171,
					want:     true,
				},
				{
					startPos: 2488105,
					limitPos: 2488172,
					want:     true,
				},
				{
					startPos: 2488105,
					limitPos: 2488173,
					want:     true,
				},
				{
					startPos: 2488171,
					limitPos: 2488172,
					want:     true,
				},
				{
					startPos: 2488171,
					limitPos: 2488173,
					want:     true,
				},
			},
		},
	}
	// Mock a *sam.Header and calling NewBEDUnionFromPath()
	ref1, _ := sam.NewReference("chr1", "", "", 249250621, nil, nil)
	ref2, _ := sam.NewReference("chr2", "", "", 243199373, nil, nil)
	samHeader, _ := sam.NewHeader(nil, []*sam.Reference{ref1, ref2})
	opts := NewBEDOpts{
		SAMHeader: samHeader,
	}
	for _, tt := range tests {
		bedUnion, err := NewBEDUnionFromPath(
			tt.pathname,
			opts,
		)
		expect.NoError(t, err)
		for _, ii := range tt.insts {
			result := bedUnion.IntersectsByID(ii.refID, ii.startPos, ii.limitPos)
			if result != ii.want {
				t.Errorf("Unexpected result for test case %v", ii.want, result)
			}
		}
	}
}

type SubsetTestInstance struct {
	startRefID int
	startPos   PosType
	limitRefID int
	limitPos   PosType
	want       BEDUnion
}

func TestSubset(t *testing.T) {
	tests := []struct {
		pathname string
		insts    []SubsetTestInstance
	}{
		{
			pathname: "testdata/test1.bed",
			insts: []SubsetTestInstance{
				{
					limitPos: 2488104,
					want: BEDUnion{
						nameMap:   make(map[string]([]PosType)),
						idMap:     make([][]PosType, 2),
						RefNames:  []string{"chr1", "chr2"},
						lastRefID: -1,
					},
				},
				{
					limitPos: 2488105,
					want: BEDUnion{
						nameMap: map[string]([]PosType){
							"chr1": []PosType{
								2488104, 2488105},
						},
						idMap: [][]PosType{
							[]PosType{
								2488104, 2488105},
							nil,
						},
						RefNames:  []string{"chr1", "chr2"},
						lastRefID: -1,
					},
				},
				{
					limitPos: 2488172,
					want: BEDUnion{
						nameMap: map[string]([]PosType){
							"chr1": []PosType{
								2488104, 2488172},
						},
						idMap: [][]PosType{
							[]PosType{
								2488104, 2488172},
							nil,
						},
						RefNames:  []string{"chr1", "chr2"},
						lastRefID: -1,
					},
				},
				{
					limitPos: 2489165,
					want: BEDUnion{
						nameMap: map[string]([]PosType){
							"chr1": []PosType{
								2488104, 2488172},
						},
						idMap: [][]PosType{
							[]PosType{
								2488104, 2488172},
							nil,
						},
						RefNames:  []string{"chr1", "chr2"},
						lastRefID: -1,
					},
				},
				{
					limitPos: 2489166,
					want: BEDUnion{
						nameMap: map[string]([]PosType){
							"chr1": []PosType{
								2488104, 2488172,
								2489165, 2489166,
							},
						},
						idMap: [][]PosType{
							[]PosType{
								2488104, 2488172,
								2489165, 2489166,
							},
							nil,
						},
						RefNames:  []string{"chr1", "chr2"},
						lastRefID: -1,
					},
				},
				{
					startPos: 2488104,
					limitPos: 2489166,
					want: BEDUnion{
						nameMap: map[string]([]PosType){
							"chr1": []PosType{
								2488104, 2488172,
								2489165, 2489166,
							},
						},
						idMap: [][]PosType{
							[]PosType{
								2488104, 2488172,
								2489165, 2489166,
							},
							nil,
						},
						RefNames:  []string{"chr1", "chr2"},
						lastRefID: -1,
					},
				},
				{
					startPos: 2488105,
					limitPos: 2489166,
					want: BEDUnion{
						nameMap: map[string]([]PosType){
							"chr1": []PosType{
								2488105, 2488172,
								2489165, 2489166,
							},
						},
						idMap: [][]PosType{
							[]PosType{
								2488105, 2488172,
								2489165, 2489166,
							},
							nil,
						},
						RefNames:  []string{"chr1", "chr2"},
						lastRefID: -1,
					},
				},
				{
					limitPos: 3000000,
					want: BEDUnion{
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
								2494587, 2494712,
							},
						},
						idMap: [][]PosType{
							[]PosType{
								2488104, 2488172,
								2489165, 2489273,
								2489782, 2489907,
								2490320, 2490438,
								2491262, 2491417,
								2492063, 2492157,
								2493112, 2493254,
								2494304, 2494335,
								2494587, 2494712,
							},
							nil,
						},
						RefNames:  []string{"chr1", "chr2"},
						lastRefID: -1,
					},
				},
				{
					startPos: 2494303,
					limitPos: 3000000,
					want: BEDUnion{
						nameMap: map[string]([]PosType){
							"chr1": []PosType{
								2494304, 2494335,
								2494587, 2494712,
							},
						},
						idMap: [][]PosType{
							[]PosType{
								2494304, 2494335,
								2494587, 2494712,
							},
							nil,
						},
						RefNames:  []string{"chr1", "chr2"},
						lastRefID: -1,
					},
				},
				{
					startPos: 2494304,
					limitPos: 3000000,
					want: BEDUnion{
						nameMap: map[string]([]PosType){
							"chr1": []PosType{
								2494304, 2494335,
								2494587, 2494712,
							},
						},
						idMap: [][]PosType{
							[]PosType{
								2494304, 2494335,
								2494587, 2494712,
							},
							nil,
						},
						RefNames:  []string{"chr1", "chr2"},
						lastRefID: -1,
					},
				},
				{
					startPos: 2494305,
					limitPos: 3000000,
					want: BEDUnion{
						nameMap: map[string]([]PosType){
							"chr1": []PosType{
								2494305, 2494335,
								2494587, 2494712,
							},
						},
						idMap: [][]PosType{
							[]PosType{
								2494305, 2494335,
								2494587, 2494712,
							},
							nil,
						},
						RefNames:  []string{"chr1", "chr2"},
						lastRefID: -1,
					},
				},
			},
		},
		{
			pathname: "testdata/test2.bed",
			insts: []SubsetTestInstance{
				{
					startPos:   2489166,
					limitRefID: 1,
					limitPos:   2490439,
					want: BEDUnion{
						nameMap: map[string]([]PosType){
							"chr1": []PosType{
								2489166, 2489273,
								2489782, 2489907,
							},
							"chr2": []PosType{
								2490320, 2490439,
							},
						},
						idMap: [][]PosType{
							[]PosType{
								2489166, 2489273,
								2489782, 2489907,
							},
							[]PosType{
								2490320, 2490439,
							},
						},
						RefNames:  []string{"chr1", "chr2"},
						lastRefID: -1,
					},
				},
			},
		},
	}

	// Mock a *sam.Header and calling NewBEDUnionFromPath()
	ref1, _ := sam.NewReference("chr1", "", "", 249250621, nil, nil)
	ref2, _ := sam.NewReference("chr2", "", "", 243199373, nil, nil)
	samHeader, _ := sam.NewHeader(nil, []*sam.Reference{ref1, ref2})
	opts := NewBEDOpts{
		SAMHeader: samHeader,
	}
	for _, tt := range tests {
		bedUnion, err := NewBEDUnionFromPath(
			tt.pathname,
			opts,
		)
		expect.NoError(t, err)
		for _, ii := range tt.insts {
			result := bedUnion.Subset(ii.startRefID, ii.startPos, ii.limitRefID, ii.limitPos)
			if !reflect.DeepEqual(result, ii.want) {
				t.Errorf("Wanted: %v  Got: %v", ii.want, result)
			}
		}
	}
}
