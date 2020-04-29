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
	"context"
)

type Opts struct {
	// Commandline options.
	BedPath      string
	Region       string
	BamIndexPath string
	Clip         int
	Cols         string
	FlagExclude  int
	Mapq         int
	MaxReadLen   int
	MaxReadSpan  int
	MinBagDepth  int
	MinBaseQual  int
	Parallelism  int
	PerStrand    bool
	RemoveSq     bool
	Stitch       bool
	TempDir      string
}

func Pileup(ctx context.Context, xampath, fapath, format, outPrefix string, opts *Opts) (err error) {
	return
}
