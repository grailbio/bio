// Copyright 2018 GRAIL, Inc. All rights reserved.
// Use of this source code is governed by the Apache 2.0
// license that can be found in the LICENSE file.

package bam

import (
	"sync"
	"testing"

	"github.com/stretchr/testify/require"
)

// Test the case where each goroutine calls Get immediately followed by Put.
func TestIndependentGets(t *testing.T) {
	p := NewRecordFreePool(func() *Record { return &Record{Magic: Magic} }, -1)
	wg := sync.WaitGroup{}
	const numThreads = 100
	for i := 0; i < numThreads; i++ {
		wg.Add(1)
		go func() {
			defer wg.Done()
			for i := 0; i < 10000; i++ {
				v := p.Get()
				require.Equal(t, Magic, v.Magic)
				p.Put(v)
			}
		}()
	}
	wg.Wait()
	// Allow some slack per thread.,
	require.Truef(t, p.ApproxLen() <= numThreads*2, "Pool too large: %v", p.ApproxLen())
}

// Test the case where each goroutine calls Get, and lets another goroutine calls Put.
func TestPutsByAnotherThread(t *testing.T) {
	const numThreads = 100
	const getsPerThread = 1000
	ch := make(chan *Record, numThreads)
	p := NewRecordFreePool(func() *Record { return &Record{Magic: Magic} }, -1)
	// Getters
	getterWg := sync.WaitGroup{}
	for i := 0; i < numThreads; i++ {
		getterWg.Add(1)
		go func() {
			defer getterWg.Done()
			for i := 0; i < getsPerThread; i++ {
				v := p.Get()
				require.Equal(t, Magic, v.Magic)
				ch <- v
			}
		}()
	}

	// Putters
	putterWg := sync.WaitGroup{}
	for i := 0; i < numThreads/2; i++ {
		putterWg.Add(1)
		go func() {
			defer putterWg.Done()
			for v := range ch {
				require.Equal(t, Magic, v.Magic)
				p.Put(v)
			}
		}()
	}
	getterWg.Wait()
	close(ch)
	putterWg.Wait()
	// Allow some slack
	require.Truef(t, p.ApproxLen() <= numThreads*getsPerThread/20, "Pool too large: %v", p.ApproxLen())
}
