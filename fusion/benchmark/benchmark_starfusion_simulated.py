#!/usr/bin/env python3

"""
Run Starfusion on simulated datasets.
Results are stored in /scratch-nvme/xyang/result/
STAR-Fusion required files are in directory: /scratch-nvme/starfusion, which contain GRCh38_v27_CTAT_lib_Feb092018/ and
STAR-Fusion-v1.5.0/ (see starfusion wiki) 

Usage:
python3.6 sim_starfusion_benchmark.py --cache_dir /scratch-nvme/xyang/tmp --result_dir /scratch-nvme/xyang/result
"""

import re
import argparse
import logging
from pathlib import Path
import os
import sys
from typing import Any, List, Set
import util

from rna_benchmark import run_starfusion


def main() -> None:
    logging.basicConfig(
        level=logging.DEBUG,
        format="%(asctime)s:%(levelname)s: %(message)s")
    p = argparse.ArgumentParser()
    p.add_argument('--cache_dir', default=util.DEFAULT_CACHE_DIR,
                   help='Benchmark cache dir')
    p.add_argument('--result_dir', default=util.DEFAULT_RESULT_DIR,
                   help='Benchmark result dir')
    p.add_argument('--starfusion_data_dir', default='/scratch-nvme/starfusion',
                   help='Directory for expanding starfusion plug-n-play files')
    p.add_argument('--starfusion_plug_n_play_targz', default=os.environ['HOME'] + '/GRCh38_v27_CTAT_lib_Feb092018.plug-n-play.tar.gz',
                   help='Tar.gz file of starfusion plug-n-play file. https://github.com/STAR-Fusion/STAR-Fusion/wiki#data-resources-required')
    p.add_argument('--starfusion_targz', default=os.environ['HOME'] + '/STAR-Fusion-v1.5.0.FULL.tar.gz',
                   help='Tar.gz file of starfusion source package. https://github.com/STAR-Fusion/STAR-Fusion/wiki#data-resources-required')

    args = p.parse_args()


    for sample in util.SIMULATED_SAMPLES:
        util.s3_cache_files([sample.path.r1, sample.path.r2], args.cache_dir)
        fastq_files: List[str] = []
        cached_file_pairs: List[util.FASTQPair] = []

        fastq_files += [sample.path.r1, sample.path.r2]
        cached_file_pairs.append(util.FASTQPair(r1=args.cache_dir + '/' + os.path.basename(sample.path.r1),
                                                r2=args.cache_dir + '/' + os.path.basename(sample.path.r2)))
        print(cached_file_pairs)
        sample_name = str(sample.n) + "_" + str(sample.coverage)
        run_starfusion(sample_name, cached_file_pairs, args)

main()
