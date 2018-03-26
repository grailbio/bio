package main

import (
	"runtime"
	"testing"
)

func init() {
	failOnError = false
	verbose = false
}

/*
Initial benchmarks.

MacBook Pro (15-inch, Late 2016)
2.9 GHz Intel Core i7, 16 GB 2133 MHz LPDDR3

Benchmark_CountPairs1-8     	       3	 375762316 ns/op	211787285 B/op	 1758968 allocs/op
Benchmark_CountPairs2-8     	       5	 267890875 ns/op	218935889 B/op	 1848307 allocs/op
Benchmark_CountPairs3-8     	       5	 249959634 ns/op	225385998 B/op	 1934129 allocs/op
Benchmark_CountPairs4-8     	       5	 249393217 ns/op	235748201 B/op	 2027445 allocs/op
Benchmark_CountPairs8-8     	      10	 264018949 ns/op	269135950 B/op	 2382116 allocs/op
Benchmark_CountPairsMax-8   	       5	 268260828 ns/op	268969185 B/op	 2378565 allocs/op
PASS

MacBook (Retina, 12-inch, Early 2016)
1.1 GHz Intel Core m3, 8 GB 1867 MHz LPDDR3

Benchmark_CountPairs1-4     	       2	 568627733 ns/op	210097928 B/op	 1756363 allocs/op
Benchmark_CountPairs2-4     	       2	 523113283 ns/op	216598248 B/op	 1845816 allocs/op
Benchmark_CountPairs3-4     	       3	 458764289 ns/op	222152538 B/op	 1928867 allocs/op
Benchmark_CountPairs4-4     	       3	 474702404 ns/op	228799285 B/op	 2017212 allocs/op
Benchmark_CountPairsMax-4   	       2	 509178182 ns/op	229663360 B/op	 2021222 allocs/op
PASS

Mac Pro (Late 2013)
3.5 GHz 6-Core Intel Xeon E5, 32 GB 1866 MHz DDR3
Benchmark_CountPairs1-12      	       3	 340725526 ns/op	224403325 B/op	 1761717 allocs/op
Benchmark_CountPairs2-12      	       5	 271135703 ns/op	231939115 B/op	 1852050 allocs/op
Benchmark_CountPairs3-12      	       5	 273473505 ns/op	237919000 B/op	 1939311 allocs/op
Benchmark_CountPairs4-12      	       5	 251185258 ns/op	246385561 B/op	 2029610 allocs/op
Benchmark_CountPairs8-12      	       5	 257888054 ns/op	284316904 B/op	 2371840 allocs/op
Benchmark_CountPairs12-12     	       5	 282670759 ns/op	330229716 B/op	 2760856 allocs/op
Benchmark_CountPairsMax-12    	       5	 322025468 ns/op	333048316 B/op	 2767174 allocs/op
PASS

After D4897 and D4922:

Mac Pro (Late 2013)
3.5 GHz 6-Core Intel Xeon E5, 32 GB 1866 MHz DDR3
Benchmark_CountPairs1-12      	       5	 263996652 ns/op	158518982 B/op	  453966 allocs/op
Benchmark_CountPairs2-12      	       5	 207860727 ns/op	161997260 B/op	  544847 allocs/op
Benchmark_CountPairs3-12      	      10	 166112231 ns/op	172258724 B/op	  639017 allocs/op
Benchmark_CountPairs4-12      	      10	 157606547 ns/op	179560038 B/op	  732586 allocs/op
Benchmark_CountPairs8-12      	      10	 180858378 ns/op	213006017 B/op	 1106365 allocs/op
Benchmark_CountPairs12-12     	      10	 191116914 ns/op	243496251 B/op	 1478869 allocs/op
Benchmark_CountPairsMax-12    	      10	 183917619 ns/op	243799584 B/op	 1478957 allocs/op
*/

func benchmarkCountPairs(cpus int, b *testing.B) {
	if cpus > runtime.NumCPU() {
		b.Skipf("only have %v cpus", runtime.NumCPU())
	}
	b.ReportAllocs()
	for i := 0; i < b.N; i++ {
		countPairs(cpus, 100000)
	}
}

func Benchmark_CountPairs1(b *testing.B) {
	benchmarkCountPairs(1, b)
}

func Benchmark_CountPairs2(b *testing.B) {
	benchmarkCountPairs(2, b)
}

func Benchmark_CountPairs3(b *testing.B) {
	benchmarkCountPairs(3, b)
}

func Benchmark_CountPairs4(b *testing.B) {
	benchmarkCountPairs(4, b)
}

func Benchmark_CountPairs8(b *testing.B) {
	benchmarkCountPairs(8, b)
}

func Benchmark_CountPairs12(b *testing.B) {
	benchmarkCountPairs(12, b)
}

func Benchmark_CountPairs16(b *testing.B) {
	benchmarkCountPairs(16, b)
}

func Benchmark_CountPairs20(b *testing.B) {
	benchmarkCountPairs(20, b)
}

func Benchmark_CountPairsMax(b *testing.B) {
	benchmarkCountPairs(runtime.NumCPU(), b)
}
