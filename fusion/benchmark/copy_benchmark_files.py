#!/usr/bin/env python3

"""Copy files from the original locations to under s3://grail-publications/2910-ISMB.

"""

import os
import logging
from typing import Set, List

import util

def copy_files(src_files: List[str], dest_dir: str) -> None:
    basenames: Set[str] = set()
    for path in src_files:
        basename = os.path.basename(path)
        if basename in basenames:
            raise Exception("Duplicate filename: " + path)
        basenames.add(basename)
    logging.info('%s -> %s', src_files, dest_dir)
    util.check_call([str(util.grail_file_path()), 'cp', '-v'] + src_files + [dest_dir])

def main() -> None:
    logging.basicConfig(level=logging.DEBUG)
    copy_files(util.REFERENCE_PATHS, util.REFERENCE_DIR)
    src_paths: List[str] = []
    for sim_sample in util.ORG_SIMULATED_SAMPLES:
        src_paths += [sim_sample.path.r1, sim_sample.path.r2]
    copy_files(src_paths, util.SIMULATED_BENCHMARK_DIR)

    src_paths = []
    for titration_sample in util.ORG_TITRATION_SAMPLES:
        src_paths += util.expand_fastq_files(titration_sample.paths)
    copy_files(src_paths, util.TITRATION_BENCHMARK_DIR)

    src_paths = []
    for rna_sample in util.ORG_RNA_SAMPLES:
        src_paths += util.expand_fastq_files(rna_sample.paths)
    copy_files(src_paths, util.RNA_BENCHMARK_DIR)

main()
