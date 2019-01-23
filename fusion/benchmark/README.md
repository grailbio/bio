This directory contains scripts we used internally to run benchmarks presented
in the ISMB paper.  It may or may not run in your environment.

Usage:

1. Edit `util.py`, especially `DEFAULT_CACHE_DIR` and  `DEFAULT_RESULT_DIR`.

2. Run the benchmark

```
    ./simulated_benchmark.py # for simulated dataset
    ./titration_benchmark # for titration samples
    ./rna_benchmark.py # for RNA samples
```

See the toplevel README.md file for more details.
