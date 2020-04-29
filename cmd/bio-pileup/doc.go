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

/*
Given a BAM or PAM, and a BED file describing genomic positions of interest,
bio-pileup reports the number of reads supporting each allele at each position.
This command is similar to "bcftools mpileup".

There are options for "collapsing" the two ends of a read-pair together, or
entire duplicate-sets ("bags") identified by doppelmark.

Only SNPs are currently reported, but indel support is very likely to be added
in the future.

Sample usage:
bio-pileup \
    --bed my-regions.bed \
    --out output-prefix \
    my.bam \
    ref.fa
*/
package main
