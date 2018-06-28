// Copyright 2018 GRAIL, Inc.  All rights reserved.
// Use of this source code is governed by the Apache-2.0
// license that can be found in the LICENSE file.

// +build amd64,!appengine

package biosimd

import (
	"github.com/grailbio/base/simd"
)

// amd64 compile-time constants.  Private base/simd constants are recalculated
// here; probably want to change that.

// BytesPerWord is the number of bytes in a machine word.
const BytesPerWord = simd.BytesPerWord

// Log2BytesPerWord is log2(BytesPerWord).  This is relevant for manual
// bit-shifting when we know that's a safe way to divide and the compiler does
// not (e.g. dividend is of signed int type).
const Log2BytesPerWord = simd.Log2BytesPerWord

// const minPageSize = 4096  may be relevant for safe functions soon.

// These could be compile-time constants for now, but not after AVX2
// autodetection is added.

// bytesPerVec is the size of the maximum-width vector that may be used.  It is
// currently always 16, but it will be set to larger values at runtime in the
// future when AVX2/AVX-512/etc. is detected.
var bytesPerVec int

// log2BytesPerVec supports efficient division by bytesPerVec.
var log2BytesPerVec uint

// *** the following function is defined in biosimd_amd64.s

//go:noescape
func hasSSE42Asm() bool

// *** end assembly function signatures

func init() {
	if !hasSSE42Asm() {
		panic("SSE4.2 required.")
	}
	bytesPerVec = 16
	log2BytesPerVec = 4
}
