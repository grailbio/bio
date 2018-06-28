// Copyright 2018 GRAIL, Inc.  All rights reserved.
// Use of this source code is governed by the Apache-2.0
// license that can be found in the LICENSE file.

// Package biosimd provides access to SIMD-based implementations of several
// common .bam/.fa/etc.-specific operations on byte arrays which the compiler
// cannot be trusted to autovectorize within the next several years.
//
// See base/simd/doc.go for more comments on the overall design.
package biosimd
