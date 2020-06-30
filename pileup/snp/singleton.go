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
	"github.com/grailbio/base/recordio/recordiozstd"
	"github.com/grailbio/bio/encoding/fasta"
	"github.com/grailbio/bio/pileup"
)

// PosType is the integer type used to represent genomic positions.
type PosType = pileup.PosType

// PosTypeMax is the maximum value that can be represented by a PosType.
const PosTypeMax = pileup.PosTypeMax

// FaEncoding is the fasta in-memory encoding expected by snp.Pileup().
// (Seq8 is actually worse than both ASCII and Base5 for this SNP-pileup, but
// it simplifies future extension to indels.)
const FaEncoding = fasta.Seq8

func init() {
	recordiozstd.Init()
}
