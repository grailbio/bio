#!/usr/bin/env python3

import glob
import argparse
import collections
import logging
from pathlib import Path
import subprocess
import util
import os
from typing import List, NamedTuple, Optional, Set, Counter, Dict

DIRS1 = [
    "s3://grail-clinical-fastq/MenloPark/170214_E00543_0033_AHCLLKALXX/2738376149/bcl2fastq-2.19.1.403-umi-1mismatch-noignore", # t1
    "s3://grail-clinical-fastq/MenloPark/170214_E00543_0034_BHCLJVALXX/2291556533/bcl2fastq-2.19.1.403-umi-1mismatch-noignore",# t2
    #"s3://grail-clinical-fastq/MenloPark/170214_E00564_0033_AHCWC7ALXX/3447259045/bcl2fastq-2.19.1.403-umi-1mismatch-noignore",# t3
    "s3://grail-clinical-fastq/MenloPark/170214_E00564_0034_BHCL2WALXX/2255989125/bcl2fastq-2.19.1.403-umi-1mismatch-noignore",# t4
    #"s3://grail-clinical-fastq/MenloPark/170214_E00543_0033_AHCLLKALXX/2738376149/bcl2fastq-2.19.1.403-umi-1mismatch-noignore",# t5
    #"s3://grail-clinical-fastq/MenloPark/170214_E00543_0034_BHCLJVALXX/2291556533/bcl2fastq-2.19.1.403-umi-1mismatch-noignore",# t6
    "s3://grail-clinical-fastq/MenloPark/170214_E00564_0033_AHCWC7ALXX/3447259045/bcl2fastq-2.19.1.403-umi-1mismatch-noignore",# t7
    #"s3://grail-clinical-fastq/MenloPark/170214_E00564_0034_BHCL2WALXX/2255989125/bcl2fastq-2.19.1.403-umi-1mismatch-noignore",# t8
    #"s3://grail-clinical-fastq/MenloPark/170214_E00543_0033_AHCLLKALXX/2738376149/bcl2fastq-2.19.1.403-umi-1mismatch-noignore",# t9
    "s3://grail-clinical-fastq/MenloPark/170214_E00543_0034_BHCLJVALXX/2291556533/bcl2fastq-2.19.1.403-umi-1mismatch-noignore",# t10
    #"s3://grail-clinical-fastq/MenloPark/170214_E00564_0033_AHCWC7ALXX/3447259045/bcl2fastq-2.19.1.403-umi-1mismatch-noignore",# t11
    #"s3://grail-clinical-fastq/MenloPark/170214_E00543_0033_AHCLLKALXX/2738376149/bcl2fastq-2.19.1.403-umi-1mismatch-noignore",# t12
    "s3://grail-clinical-fastq/MenloPark/170214_E00543_0034_BHCLJVALXX/2291556533/bcl2fastq-2.19.1.403-umi-1mismatch-noignore",# t13
    "s3://grail-clinical-fastq/MenloPark/170214_E00564_0033_AHCWC7ALXX/3447259045/bcl2fastq-2.19.1.403-umi-1mismatch-noignore"]# t14

DIRS0 = [
    "170206_ARTLoD_B1_01rerun",# t1
    "170206_ARTLoD_B1_02rerun",# t2
    #"170206_ARTLoD_B1_03rerun",# t3
    "170206_ARTLoD_B1_04rerun",# t4
    #"170206_ARTLoD_B1_05rerun",# t5
    #"170206_ARTLoD_B1_06rerun",# t6
    "170206_ARTLoD_B1_07rerun",# t7
    #"170206_ARTLoD_B1_08rerun",# t8
    #"170206_ARTLoD_B1_09rerun",# t9
    "170206_ARTLoD_B1_10rerun",# t10
    #"170206_ARTLoD_B1_11rerun",# t11
    #"170206_ARTLoD_B1_12rerun",# t12
    "170206_ARTLoD_B1_13rerun",# t13
    "170206_ARTLoD_B1_14rerun"]# t14

def main() -> None:
    logging.basicConfig(
        level=logging.DEBUG,
        format="%(asctime)s:%(levelname)s: %(message)s")
    p = argparse.ArgumentParser()
    p.add_argument('--cache_dir', default='/scratch-nvme/cache_tmp',
                   help='Benchmark cache dir')
    p.add_argument('--result_dir', default='/scratch-sata/af4_benchmark_results',
                   help='Benchmark result dir')
    args = p.parse_args()

    ref_path = "s3://grail-pecan/resources/annotation/gencode.v26.whole_genes.fa"
    cosmic_fusion_path = "s3://grail-jkim/tmp/all_pair_art_lod_gpair_merged.txt"
    for i in range(0, len(DIRS0)):
        logging.info('Start benchmark %d: %s', i, DIRS1[i])
        dir_path = DIRS1[i] + '/' + DIRS0[i]

        result_dir = args.result_dir + '/' + os.path.basename(DIRS0[i])
        try:
            os.makedirs(result_dir, 0o755)
        except:
            logging.error("mkdir %s failed", result_dir)
        if os.path.exists(result_dir + "/filtered.fa"):
            logging.info("Skip %s", result_dir)
            continue

        paths = util.s3_ls(dir_path)
        r1 = sorted([path for path in paths if "R1" in path])
        r2 = sorted([path for path in paths if "R2" in path])

        assert len(r1) == len(r2)
        for j in range(len(r1)):
            assert r1[j].replace("R1", "R2"), r2[j]


        util.s3_cache_files(r1 + r2 + [ref_path, cosmic_fusion_path], args.cache_dir)
        cached_r1 = ",".join([args.cache_dir + '/' + os.path.basename(path) for path in r1])
        cached_r2 = ",".join([args.cache_dir + '/' + os.path.basename(path) for path in r2])
        cached_ref = args.cache_dir + '/' + os.path.basename(ref_path)
        cached_cosmic_fusion = args.cache_dir + '/' + os.path.basename(cosmic_fusion_path)

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
                    f'-cosmic-fusion={cached_cosmic_fusion}']
        util.check_call(af4_args)
        logging.info('Finished benchmark %d: %s', i, DIRS1[i])
        logging.info("Runtime stats: %s", util.run_stats(Path(result_dir)))
        for path in glob.glob(f'{args.cache_dir}/*rerun*'):
            try:
                os.remove(path)
            except:
                logging.error("failed to remove " + path)
main()
