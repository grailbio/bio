// Copyright 2018 GRAIL, Inc. All rights reserved.
// Use of this source code is governed by the Apache 2.0
// license that can be found in the LICENSE file.

package pam

import (
	"v.io/x/lib/vlog"
)

// unsafeArena is an arena allocator. It supports allocating []bytes quickly.
type unsafeArena struct {
	buf []byte
	n   int // # bytes allocated so far
}

// align rounds "ub.n" up so that it is a multiple of 8. The next alloc() call
// returns a block aligned at a 8-byte boundary. Used when storing a pointer in
// []byte.  8-byte alignment is sufficient for all CPUs we care about.
func (ub *unsafeArena) align() {
	const pointerSize = 8
	ub.n = ((ub.n-1)/pointerSize + 1) * pointerSize
}

// alloc allocates a byte slice of "size" bytes.
//
// Requires: ub must have at least size bytes of free space.
func (ub *unsafeArena) alloc(size int) []byte {
	if ub.n+size > len(ub.buf) {
		vlog.Fatalf("Arena overflow, n=%d, size=%d, ub=%d", ub.n, size, len(ub.buf))
	}
	a := ub.buf[ub.n : ub.n+size]
	ub.n += size
	return a
}
