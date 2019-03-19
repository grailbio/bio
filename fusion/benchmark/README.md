This directory contains scripts we used internally to run benchmarks presented
in the ISMB paper.  It may or may not run in your environment.

Usage:

1. Edit `util.py`, especially `DEFAULT_CACHE_DIR` and  `DEFAULT_RESULT_DIR`.

2. Run the benchmark

```
    ./simulated_benchmark.py # for simulated dataset
    ./titration_benchmark # for titration samples

    ./rna_benchmark.py # for AF4 and STAR-Fusion RNA cfRNA samples and AF4 for RNA-Seq brca samples
                       # all_art_lod_brca.txt file is used as the comsic fusion input file that included all fusions in those breast cancer samples. This file (s3://grail-publications/2019-ISMB/references/all_art_lod_brca.txt) is the s3://grail-publications/2019-ISMB/references/all_pair_art_lod_gpair_merged.txt with 27 added fusions from Edgren et al. publication table 1 (https://genomebiology.biomedcentral.com/articles/10.1186/gb-2011-12-1-r6)


    # for STAR-Fusion on simulated dataset
    ./benchmark_starfusion_simulated.py
    ./benchmark_starfusion_result_analysis.py

```

See the toplevel README.md file for more details.
