// Copyright 2018 GRAIL, Inc. All rights reserved.
// Use of this source code is governed by the Apache 2.0
// license that can be found in the LICENSE file.

package fieldio

import (
	"github.com/grailbio/base/log"
)

// UnsafeArena is an arena allocator. It supports allocating []bytes quickly.
type UnsafeArena struct {
	buf []byte
	n   int // # bytes allocated so far
}

// NewUnsafeArena creates an arena for filling the given buffer.
func NewUnsafeArena(buf []byte) UnsafeArena {
	return UnsafeArena{buf: buf, n: 0}
}

// Align rounds "ub.n" up so that it is a multiple of 8. The next alloc() call
// returns a block aligned at a 8-byte boundary. Used when storing a pointer in
// []byte.  8-byte alignment is sufficient for all CPUs we care about.
func (ub *UnsafeArena) Align() {
	const pointerSize = 8
	ub.n = ((ub.n-1)/pointerSize + 1) * pointerSize
}

// Alloc allocates a byte slice of "size" bytes.
//
// Requires: ub must have at least size bytes of free space.
func (ub *UnsafeArena) Alloc(size int) []byte {
	if ub.n+size > len(ub.buf) {
		log.Panicf("Arena overflow, n=%d, size=%d, ub=%d", ub.n, size, len(ub.buf))
	}
	a := ub.buf[ub.n : ub.n+size]
	ub.n += size
	return a
}
