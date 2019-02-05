#!/usr/bin/env python3

import glob
import argparse
import collections
import logging
from pathlib import Path
import subprocess
import os
from typing import List, NamedTuple, Optional, Set, Counter, Dict

import util

def main() -> None:
    logging.basicConfig(
        level=logging.DEBUG,
        format="%(asctime)s:%(levelname)s: %(message)s")
    p = argparse.ArgumentParser()
    p.add_argument('--cache_dir', default=util.DEFAULT_CACHE_DIR,
                   help='Benchmark cache dir')
    p.add_argument('--result_dir', default=util.DEFAULT_RESULT_DIR,
                   help='Benchmark result dir')
    args = p.parse_args()

    util.s3_cache_files([util.REFERENCE_DIR + '/gencode.v26.whole_genes.fa',
                         util.REFERENCE_DIR + '/all_pair_art_lod_gpair_merged.txt'],
                        args.cache_dir)
    for sample in util.TITRATION_SAMPLES:
        logging.info('Start benchmark %s', sample.name)
        result_dir = args.result_dir + '/' + sample.name
        try:
            os.makedirs(result_dir, 0o755)
        except:
            logging.error("mkdir %s failed", result_dir)
        if os.path.exists(result_dir + "/filtered.fa"):
            logging.info("Skip %s", result_dir)
            continue
        util.s3_cache_files(util.expand_fastq_files(sample.paths), args.cache_dir)
        cached_r1 = ",".join([args.cache_dir + '/' + os.path.basename(fq.r1) for fq in sample.paths])
        cached_r2 = ",".join([args.cache_dir + '/' + os.path.basename(fq.r2) for fq in sample.paths])
        cached_ref = args.cache_dir + '/gencode.v26.whole_genes.fa'
        cached_cosmic_fusion = args.cache_dir + '/all_pair_art_lod_gpair_merged.txt'

        af4_args = [str(util.af4_path()),
                    f'-log_dir={result_dir}',
                    f'-pprof=:12345',
                    f'-mutex-profile-rate=1000',
                    f'-block-profile-rate=1000',
                    f'-r1={cached_r1}',
                    f'-r2={cached_r2}',
                    f'-fasta-output={result_dir}/all.fa',
                    f'-filtered-output={result_dir}/filtered.fa',
                    f'-transcript={cached_ref}',
                    f'-max-genes-per-kmer=2',
                    f'-max-proximity-distance=1000',
                    f'-max-proximity-genes=5',
                    f'-unstranded-prep',
                    f'-cosmic-fusion={cached_cosmic_fusion}']
        util.check_call(af4_args)
        logging.info('Finished benchmark %d: %s', sample.name)
        logging.info("Runtime stats: %s", util.run_stats(Path(result_dir)))
        for path in glob.glob(f'{args.cache_dir}/*rerun*'):
            try:
                os.remove(path)
            except:
                logging.error("failed to remove " + path)

if __name__ == "__main__":
    main()
