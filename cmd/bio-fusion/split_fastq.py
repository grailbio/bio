#!/usr/bin/env python3.6

"""This script splits one large fastq files into multiple smaller ones."""

import glob
import logging
import sys
import os
import re
import subprocess
from concurrent.futures import ThreadPoolExecutor
from typing import Any

# FASTQ files to split.
files = [
    's3://grail-cfrna-fastq/MissionBay/161206_E00501_0049_AHCYJYALXX/3747505862/bcl2fastq-2.19.0.316-umi-1mismatch-noignore/108CPREL315/108CPREL315_S4_L004_R1_001.fastq.gz',
    's3://grail-cfrna-fastq/MissionBay/161206_E00501_0049_AHCYJYALXX/3747505862/bcl2fastq-2.19.0.316-umi-1mismatch-noignore/108CPREL315/108CPREL315_S4_L004_R2_001.fastq.gz',
    's3://grail-cfrna-fastq/MissionBay/170105_K00215_0113_AHF2N2BBXX/3062013994/bcl2fastq-2.19.0.316-umi-1mismatch-noignore/74HPREL332/74HPREL332_S1_L001_R1_001.fastq.gz',
    's3://grail-cfrna-fastq/MissionBay/170105_K00215_0113_AHF2N2BBXX/3062013994/bcl2fastq-2.19.0.316-umi-1mismatch-noignore/74HPREL332/74HPREL332_S1_L001_R2_001.fastq.gz',
    's3://grail-cfrna-fastq/MissionBay/161206_E00501_0050_BHF2JCALXX/590215334/bcl2fastq-2.19.0.316-umi-1mismatch-noignore/118HPREL322/118HPREL322_S3_L003_R1_001.fastq.gz',
    's3://grail-cfrna-fastq/MissionBay/161206_E00501_0050_BHF2JCALXX/590215334/bcl2fastq-2.19.0.316-umi-1mismatch-noignore/118HPREL322/118HPREL322_S3_L003_R2_001.fastq.gz',
    's3://grail-cfrna-fastq/MissionBay/160927_K00122_0179_BHFWLKBBXX/3600619798/bcl2fastq-2.19.0.316-umi-1mismatch-noignore/53HPREL160/53HPREL160_S3_L003_R1_001.fastq.gz',
    's3://grail-cfrna-fastq/MissionBay/160927_K00122_0179_BHFWLKBBXX/3600619798/bcl2fastq-2.19.0.316-umi-1mismatch-noignore/53HPREL160/53HPREL160_S3_L003_R2_001.fastq.gz',
    's3://grail-cfrna-fastq/MissionBay/161206_E00501_0050_BHF2JCALXX/590215334/bcl2fastq-2.19.0.316-umi-1mismatch-noignore/117HPREL321/117HPREL321_S2_L002_R1_001.fastq.gz',
    's3://grail-cfrna-fastq/MissionBay/161206_E00501_0050_BHF2JCALXX/590215334/bcl2fastq-2.19.0.316-umi-1mismatch-noignore/117HPREL321/117HPREL321_S2_L002_R2_001.fastq.gz']

# Where results are stored. Paths are in form
#
# <dir>/<basename-of-original>-NN.fastq.gz
#
# where NN=00,01,...
dest_dir = 's3://grail-ysaito/af4_benchmark'

def main():
    pool = ThreadPoolExecutor(64)
    logging.basicConfig(level=logging.DEBUG)
    results = []
    for fastq_path in files:
        results.append(pool.submit(split_fastq, fastq_path))
    for r in results:
        print(r.result())

def gzip_and_upload(path: str):
    subprocess.check_call(['gzip', path])
    subprocess.check_call(['grail-file', 'cp', path, dest_dir + '/' + os.path.basename(path) + '.gz'])

def split_fastq(s3_fastq_path: str) -> None:
    fastq_path = os.path.basename(s3_fastq_path)
    if not os.path.exists(fastq_path):
        logging.info('cp %s -> %s', s3_fastq_path, fastq_path)
        subprocess.check_call(['grail-file', 'cp', s3_fastq_path, fastq_path])

    m = re.match('.*/([^/]+).fastq.gz', fastq_path)
    if not m:
        m = re.match('([^/]+).fastq.gz', fastq_path)
    assert m, fastq_path
    basename = m[1]
    cmdline = f'zcat {fastq_path} | split -d -l 209715200 - {basename}-'
    logging.info('Run: %s', cmdline)
    subprocess.check_call(cmdline, shell=True)
    pool = ThreadPoolExecutor(64)

    results: Any = []
    for split_shard_path in glob.glob(basename + "-[0-9][0-9]"):
        fastq_shard_path = split_shard_path + '.fastq'
        os.rename(split_shard_path, fastq_shard_path)
        results.append(pool.submit(gzip_and_upload, fastq_shard_path))
        #results.append(pool.submit(subprocess.check_call, ['gzip', fastq_shard_path]))

    for r in results:
        print(r.result())

main()
