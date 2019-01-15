#!/usr/bin/env python3

"""
Run cfRNA benchmarks. The default parameter settings are for ubuntu*.mpk machines.
Results are stored in /scratch-nvme/af4_benchmark_results

Usage:
rna_benchmark.py --run={af4,starfusion}

"""
import glob
import re
import argparse
import collections
import logging
from pathlib import Path
import subprocess
import util
import os
from typing import Any, List, NamedTuple, Optional, Set, Counter, Dict


FilePair = NamedTuple('FilePair', [('r1', str), ('r2', str)])

# Note: files in s3://grail-ysaito/af4_benchmark are created by running split_fastq.py on the original fastq pair.

SAMPLES = {
        '101CPREL277' :	[
            FilePair(r1='s3://grail-cfrna-fastq/MissionBay/161115_E00481_0058_AH53VWALXX/1151367959/bcl2fastq-2.19.0.316-umi-1mismatch-noignore/101CPREL277/101CPREL277_S1_L001_R1_001.fastq.gz',
                     r2='s3://grail-cfrna-fastq/MissionBay/161115_E00481_0058_AH53VWALXX/1151367959/bcl2fastq-2.19.0.316-umi-1mismatch-noignore/101CPREL277/101CPREL277_S1_L001_R2_001.fastq.gz'),
            FilePair(r1='s3://grail-cfrna-fastq/MissionBay/161115_E00481_0058_AH53VWALXX/1151367959/bcl2fastq-2.19.0.316-umi-1mismatch-noignore/101CPREL277/101CPREL277_S1_L002_R1_001.fastq.gz',
                     r2='s3://grail-cfrna-fastq/MissionBay/161115_E00481_0058_AH53VWALXX/1151367959/bcl2fastq-2.19.0.316-umi-1mismatch-noignore/101CPREL277/101CPREL277_S1_L002_R2_001.fastq.gz'),
            FilePair(r1='s3://grail-cfrna-fastq/MissionBay/161115_E00481_0058_AH53VWALXX/1151367959/bcl2fastq-2.19.0.316-umi-1mismatch-noignore/101CPREL277/101CPREL277_S1_L003_R1_001.fastq.gz',
                     r2='s3://grail-cfrna-fastq/MissionBay/161115_E00481_0058_AH53VWALXX/1151367959/bcl2fastq-2.19.0.316-umi-1mismatch-noignore/101CPREL277/101CPREL277_S1_L003_R2_001.fastq.gz'),
            FilePair(r1='s3://grail-cfrna-fastq/MissionBay/161115_E00481_0058_AH53VWALXX/1151367959/bcl2fastq-2.19.0.316-umi-1mismatch-noignore/101CPREL277/101CPREL277_S1_L004_R1_001.fastq.gz',
                     r2='s3://grail-cfrna-fastq/MissionBay/161115_E00481_0058_AH53VWALXX/1151367959/bcl2fastq-2.19.0.316-umi-1mismatch-noignore/101CPREL277/101CPREL277_S1_L004_R2_001.fastq.gz'),
            FilePair(r1='s3://grail-cfrna-fastq/MissionBay/161115_E00481_0058_AH53VWALXX/1151367959/bcl2fastq-2.19.0.316-umi-1mismatch-noignore/101CPREL277/101CPREL277_S1_L005_R1_001.fastq.gz',
                     r2='s3://grail-cfrna-fastq/MissionBay/161115_E00481_0058_AH53VWALXX/1151367959/bcl2fastq-2.19.0.316-umi-1mismatch-noignore/101CPREL277/101CPREL277_S1_L005_R2_001.fastq.gz'),
            FilePair(r1='s3://grail-cfrna-fastq/MissionBay/161115_E00481_0058_AH53VWALXX/1151367959/bcl2fastq-2.19.0.316-umi-1mismatch-noignore/101CPREL277/101CPREL277_S1_L006_R1_001.fastq.gz',
                     r2='s3://grail-cfrna-fastq/MissionBay/161115_E00481_0058_AH53VWALXX/1151367959/bcl2fastq-2.19.0.316-umi-1mismatch-noignore/101CPREL277/101CPREL277_S1_L006_R2_001.fastq.gz'),
            FilePair(r1='s3://grail-cfrna-fastq/MissionBay/161115_E00481_0058_AH53VWALXX/1151367959/bcl2fastq-2.19.0.316-umi-1mismatch-noignore/101CPREL277/101CPREL277_S1_L007_R1_001.fastq.gz',
                     r2='s3://grail-cfrna-fastq/MissionBay/161115_E00481_0058_AH53VWALXX/1151367959/bcl2fastq-2.19.0.316-umi-1mismatch-noignore/101CPREL277/101CPREL277_S1_L007_R2_001.fastq.gz'),
            FilePair(r1='s3://grail-cfrna-fastq/MissionBay/161115_E00481_0058_AH53VWALXX/1151367959/bcl2fastq-2.19.0.316-umi-1mismatch-noignore/101CPREL277/101CPREL277_S1_L008_R1_001.fastq.gz',
                     r2='s3://grail-cfrna-fastq/MissionBay/161115_E00481_0058_AH53VWALXX/1151367959/bcl2fastq-2.19.0.316-umi-1mismatch-noignore/101CPREL277/101CPREL277_S1_L008_R2_001.fastq.gz')],
        '108CPREL315': [
	    #FilePair(r1='s3://grail-cfrna-fastq/MissionBay/161206_E00501_0049_AHCYJYALXX/3747505862/bcl2fastq-2.19.0.316-umi-1mismatch-noignore/108CPREL315/108CPREL315_S4_L004_R1_001.fastq.gz',
            #         r2='s3://grail-cfrna-fastq/MissionBay/161206_E00501_0049_AHCYJYALXX/3747505862/bcl2fastq-2.19.0.316-umi-1mismatch-noignore/108CPREL315/108CPREL315_S4_L004_R2_001.fastq.gz')],

            FilePair(r1='s3://grail-ysaito/af4_benchmark/108CPREL315_S4_L004_R1_001-00.fastq.gz',
                     r2='s3://grail-ysaito/af4_benchmark/108CPREL315_S4_L004_R2_001-00.fastq.gz'),
            FilePair(r1='s3://grail-ysaito/af4_benchmark/108CPREL315_S4_L004_R1_001-01.fastq.gz',
                     r2='s3://grail-ysaito/af4_benchmark/108CPREL315_S4_L004_R2_001-01.fastq.gz'),
            FilePair(r1='s3://grail-ysaito/af4_benchmark/108CPREL315_S4_L004_R1_001-02.fastq.gz',
                     r2='s3://grail-ysaito/af4_benchmark/108CPREL315_S4_L004_R2_001-02.fastq.gz'),
            FilePair(r1='s3://grail-ysaito/af4_benchmark/108CPREL315_S4_L004_R1_001-03.fastq.gz',
                     r2='s3://grail-ysaito/af4_benchmark/108CPREL315_S4_L004_R2_001-03.fastq.gz'),
            FilePair(r1='s3://grail-ysaito/af4_benchmark/108CPREL315_S4_L004_R1_001-04.fastq.gz',
                     r2='s3://grail-ysaito/af4_benchmark/108CPREL315_S4_L004_R2_001-04.fastq.gz'),
            FilePair(r1='s3://grail-ysaito/af4_benchmark/108CPREL315_S4_L004_R1_001-05.fastq.gz',
                     r2='s3://grail-ysaito/af4_benchmark/108CPREL315_S4_L004_R2_001-05.fastq.gz'),
            FilePair(r1='s3://grail-ysaito/af4_benchmark/108CPREL315_S4_L004_R1_001-06.fastq.gz',
                     r2='s3://grail-ysaito/af4_benchmark/108CPREL315_S4_L004_R2_001-06.fastq.gz'),
            FilePair(r1='s3://grail-ysaito/af4_benchmark/108CPREL315_S4_L004_R1_001-07.fastq.gz',
                     r2='s3://grail-ysaito/af4_benchmark/108CPREL315_S4_L004_R2_001-07.fastq.gz'),
            FilePair(r1='s3://grail-ysaito/af4_benchmark/108CPREL315_S4_L004_R1_001-08.fastq.gz',
                     r2='s3://grail-ysaito/af4_benchmark/108CPREL315_S4_L004_R2_001-08.fastq.gz')],


        '74HPREL332': [
            #FilePair(r1='s3://grail-cfrna-fastq/MissionBay/170105_K00215_0113_AHF2N2BBXX/3062013994/bcl2fastq-2.19.0.316-umi-1mismatch-noignore/74HPREL332/74HPREL332_S1_L001_R1_001.fastq.gz',
            #         r2='s3://grail-cfrna-fastq/MissionBay/170105_K00215_0113_AHF2N2BBXX/3062013994/bcl2fastq-2.19.0.316-umi-1mismatch-noignore/74HPREL332/74HPREL332_S1_L001_R2_001.fastq.gz')],
            FilePair(r1='s3://grail-ysaito/af4_benchmark/74HPREL332_S1_L001_R1_001-00.fastq.gz',
                     r2='s3://grail-ysaito/af4_benchmark/74HPREL332_S1_L001_R2_001-00.fastq.gz'),
            FilePair(r1='s3://grail-ysaito/af4_benchmark/74HPREL332_S1_L001_R1_001-01.fastq.gz',
                     r2='s3://grail-ysaito/af4_benchmark/74HPREL332_S1_L001_R2_001-01.fastq.gz'),
            FilePair(r1='s3://grail-ysaito/af4_benchmark/74HPREL332_S1_L001_R1_001-02.fastq.gz',
                     r2='s3://grail-ysaito/af4_benchmark/74HPREL332_S1_L001_R2_001-02.fastq.gz'),
            FilePair(r1='s3://grail-ysaito/af4_benchmark/74HPREL332_S1_L001_R1_001-03.fastq.gz',
                     r2='s3://grail-ysaito/af4_benchmark/74HPREL332_S1_L001_R2_001-03.fastq.gz')],
        '118HPREL322': [
            #FilePair(r1='s3://grail-cfrna-fastq/MissionBay/161206_E00501_0050_BHF2JCALXX/590215334/bcl2fastq-2.19.0.316-umi-1mismatch-noignore/118HPREL322/118HPREL322_S3_L003_R1_001.fastq.gz',
            #        r2='s3://grail-cfrna-fastq/MissionBay/161206_E00501_0050_BHF2JCALXX/590215334/bcl2fastq-2.19.0.316-umi-1mismatch-noignore/118HPREL322/118HPREL322_S3_L003_R2_001.fastq.gz')],
            FilePair(r1='s3://grail-ysaito/af4_benchmark/118HPREL322_S3_L003_R1_001-00.fastq.gz',
                     r2='s3://grail-ysaito/af4_benchmark/118HPREL322_S3_L003_R2_001-00.fastq.gz'),
            FilePair(r1='s3://grail-ysaito/af4_benchmark/118HPREL322_S3_L003_R1_001-01.fastq.gz',
                     r2='s3://grail-ysaito/af4_benchmark/118HPREL322_S3_L003_R2_001-01.fastq.gz'),
            FilePair(r1='s3://grail-ysaito/af4_benchmark/118HPREL322_S3_L003_R1_001-02.fastq.gz',
                     r2='s3://grail-ysaito/af4_benchmark/118HPREL322_S3_L003_R2_001-02.fastq.gz'),
            FilePair(r1='s3://grail-ysaito/af4_benchmark/118HPREL322_S3_L003_R1_001-03.fastq.gz',
                     r2='s3://grail-ysaito/af4_benchmark/118HPREL322_S3_L003_R2_001-03.fastq.gz'),
            FilePair(r1='s3://grail-ysaito/af4_benchmark/118HPREL322_S3_L003_R1_001-04.fastq.gz',
                     r2='s3://grail-ysaito/af4_benchmark/118HPREL322_S3_L003_R2_001-04.fastq.gz'),
            FilePair(r1='s3://grail-ysaito/af4_benchmark/118HPREL322_S3_L003_R1_001-05.fastq.gz',
                     r2='s3://grail-ysaito/af4_benchmark/118HPREL322_S3_L003_R2_001-05.fastq.gz'),
            FilePair(r1='s3://grail-ysaito/af4_benchmark/118HPREL322_S3_L003_R1_001-06.fastq.gz',
                     r2='s3://grail-ysaito/af4_benchmark/118HPREL322_S3_L003_R2_001-06.fastq.gz'),
            FilePair(r1='s3://grail-ysaito/af4_benchmark/118HPREL322_S3_L003_R1_001-07.fastq.gz',
                     r2='s3://grail-ysaito/af4_benchmark/118HPREL322_S3_L003_R2_001-07.fastq.gz')],
        # '68HPREL273': [
        #     FilePair(r1='s3://grail-cfrna-fastq/MissionBay/161115_E00481_0058_AH53VWALXX/1151367959/bcl2fastq-2.19.0.316-umi-1mismatch-noignore/68HPREL273/68HPREL273_S5_L001_R1_001.fastq.gz',
        #              r2='s3://grail-cfrna-fastq/MissionBay/161115_E00481_0058_AH53VWALXX/1151367959/bcl2fastq-2.19.0.316-umi-1mismatch-noignore/68HPREL273/68HPREL273_S5_L001_R2_001.fastq.gz'),
        #     FilePair(r1='s3://grail-cfrna-fastq/MissionBay/161115_E00481_0058_AH53VWALXX/1151367959/bcl2fastq-2.19.0.316-umi-1mismatch-noignore/68HPREL273/68HPREL273_S5_L002_R1_001.fastq.gz',
        #              r2='s3://grail-cfrna-fastq/MissionBay/161115_E00481_0058_AH53VWALXX/1151367959/bcl2fastq-2.19.0.316-umi-1mismatch-noignore/68HPREL273/68HPREL273_S5_L002_R2_001.fastq.gz'),
        #     FilePair(r1='s3://grail-cfrna-fastq/MissionBay/161115_E00481_0058_AH53VWALXX/1151367959/bcl2fastq-2.19.0.316-umi-1mismatch-noignore/68HPREL273/68HPREL273_S5_L003_R1_001.fastq.gz',
        #              r2='s3://grail-cfrna-fastq/MissionBay/161115_E00481_0058_AH53VWALXX/1151367959/bcl2fastq-2.19.0.316-umi-1mismatch-noignore/68HPREL273/68HPREL273_S5_L003_R2_001.fastq.gz'),
        #     FilePair(r1='s3://grail-cfrna-fastq/MissionBay/161115_E00481_0058_AH53VWALXX/1151367959/bcl2fastq-2.19.0.316-umi-1mismatch-noignore/68HPREL273/68HPREL273_S5_L004_R1_001.fastq.gz',
        #              r2='s3://grail-cfrna-fastq/MissionBay/161115_E00481_0058_AH53VWALXX/1151367959/bcl2fastq-2.19.0.316-umi-1mismatch-noignore/68HPREL273/68HPREL273_S5_L004_R2_001.fastq.gz'),
        #     FilePair(r1='s3://grail-cfrna-fastq/MissionBay/161115_E00481_0058_AH53VWALXX/1151367959/bcl2fastq-2.19.0.316-umi-1mismatch-noignore/68HPREL273/68HPREL273_S5_L005_R1_001.fastq.gz',
        #              r2='s3://grail-cfrna-fastq/MissionBay/161115_E00481_0058_AH53VWALXX/1151367959/bcl2fastq-2.19.0.316-umi-1mismatch-noignore/68HPREL273/68HPREL273_S5_L005_R2_001.fastq.gz'),
        #     FilePair(r1='s3://grail-cfrna-fastq/MissionBay/161115_E00481_0058_AH53VWALXX/1151367959/bcl2fastq-2.19.0.316-umi-1mismatch-noignore/68HPREL273/68HPREL273_S5_L006_R1_001.fastq.gz',
        #              r2='s3://grail-cfrna-fastq/MissionBay/161115_E00481_0058_AH53VWALXX/1151367959/bcl2fastq-2.19.0.316-umi-1mismatch-noignore/68HPREL273/68HPREL273_S5_L006_R2_001.fastq.gz'),
        #     FilePair(r1='s3://grail-cfrna-fastq/MissionBay/161115_E00481_0058_AH53VWALXX/1151367959/bcl2fastq-2.19.0.316-umi-1mismatch-noignore/68HPREL273/68HPREL273_S5_L007_R1_001.fastq.gz',
        #              r2='s3://grail-cfrna-fastq/MissionBay/161115_E00481_0058_AH53VWALXX/1151367959/bcl2fastq-2.19.0.316-umi-1mismatch-noignore/68HPREL273/68HPREL273_S5_L007_R2_001.fastq.gz'),
        #     FilePair(r1='s3://grail-cfrna-fastq/MissionBay/161115_E00481_0058_AH53VWALXX/1151367959/bcl2fastq-2.19.0.316-umi-1mismatch-noignore/68HPREL273/68HPREL273_S5_L008_R1_001.fastq.gz',
        #              r2='s3://grail-cfrna-fastq/MissionBay/161115_E00481_0058_AH53VWALXX/1151367959/bcl2fastq-2.19.0.316-umi-1mismatch-noignore/68HPREL273/68HPREL273_S5_L008_R2_001.fastq.gz')],
        '53HPREL160': [
            #FilePair(r1='s3://grail-cfrna-fastq/MissionBay/160927_K00122_0179_BHFWLKBBXX/3600619798/bcl2fastq-2.19.0.316-umi-1mismatch-noignore/53HPREL160/53HPREL160_S3_L003_R1_001.fastq.gz',
             #        r2='s3://grail-cfrna-fastq/MissionBay/160927_K00122_0179_BHFWLKBBXX/3600619798/bcl2fastq-2.19.0.316-umi-1mismatch-noignore/53HPREL160/53HPREL160_S3_L003_R2_001.fastq.gz')],

            FilePair(r1='s3://grail-ysaito/af4_benchmark/53HPREL160_S3_L003_R1_001-00.fastq.gz',
                     r2='s3://grail-ysaito/af4_benchmark/53HPREL160_S3_L003_R2_001-00.fastq.gz'),
            FilePair(r1='s3://grail-ysaito/af4_benchmark/53HPREL160_S3_L003_R1_001-01.fastq.gz',
                     r2='s3://grail-ysaito/af4_benchmark/53HPREL160_S3_L003_R2_001-01.fastq.gz')],

        # '117HPREL321': [
        #     #FilePair(r1='s3://grail-cfrna-fastq/MissionBay/161206_E00501_0050_BHF2JCALXX/590215334/bcl2fastq-2.19.0.316-umi-1mismatch-noignore/117HPREL321/117HPREL321_S2_L002_R1_001.fastq.gz',
        #     #r2='s3://grail-cfrna-fastq/MissionBay/161206_E00501_0050_BHF2JCALXX/590215334/bcl2fastq-2.19.0.316-umi-1mismatch-noignore/117HPREL321/117HPREL321_S2_L002_R2_001.fastq.gz')]

        #     FilePair(r1='s3://grail-ysaito/af4_benchmark/117HPREL321_S2_L002_R1_001-00.fastq.gz',
        #             r2='s3://grail-ysaito/af4_benchmark/117HPREL321_S2_L002_R2_001-00.fastq.gz'),
        #     FilePair(r1='s3://grail-ysaito/af4_benchmark/117HPREL321_S2_L002_R1_001-01.fastq.gz',
        #              r2='s3://grail-ysaito/af4_benchmark/117HPREL321_S2_L002_R2_001-01.fastq.gz'),
        #     FilePair(r1='s3://grail-ysaito/af4_benchmark/117HPREL321_S2_L002_R1_001-02.fastq.gz',
        #              r2='s3://grail-ysaito/af4_benchmark/117HPREL321_S2_L002_R2_001-02.fastq.gz'),
        #     FilePair(r1='s3://grail-ysaito/af4_benchmark/117HPREL321_S2_L002_R1_001-03.fastq.gz',
        #              r2='s3://grail-ysaito/af4_benchmark/117HPREL321_S2_L002_R2_001-03.fastq.gz'),
        #     FilePair(r1='s3://grail-ysaito/af4_benchmark/117HPREL321_S2_L002_R1_001-04.fastq.gz',
        #              r2='s3://grail-ysaito/af4_benchmark/117HPREL321_S2_L002_R2_001-04.fastq.gz'),
        #     FilePair(r1='s3://grail-ysaito/af4_benchmark/117HPREL321_S2_L002_R1_001-05.fastq.gz',
        #              r2='s3://grail-ysaito/af4_benchmark/117HPREL321_S2_L002_R2_001-05.fastq.gz'),
        #     FilePair(r1='s3://grail-ysaito/af4_benchmark/117HPREL321_S2_L002_R1_001-06.fastq.gz',
        #              r2='s3://grail-ysaito/af4_benchmark/117HPREL321_S2_L002_R2_001-06.fastq.gz'),
        #     FilePair(r1='s3://grail-ysaito/af4_benchmark/117HPREL321_S2_L002_R1_001-07.fastq.gz',
        #              r2='s3://grail-ysaito/af4_benchmark/117HPREL321_S2_L002_R2_001-07.fastq.gz'),
        #     FilePair(r1='s3://grail-ysaito/af4_benchmark/117HPREL321_S2_L002_R1_001-08.fastq.gz',
        #              r2='s3://grail-ysaito/af4_benchmark/117HPREL321_S2_L002_R2_001-08.fastq.gz')],

}

def run_starfusion(sample_name: str, cached_file_pairs: List[FilePair], args: Any):
    match = re.match(r'.*/([^/]+)\.FULL\.tar\.gz$', args.starfusion_targz)
    assert match
    local_starfusion_dir = match[1]

    logging.info("LOCAL: %s", local_starfusion_dir)
    if not os.path.exists(os.path.join(args.starfusion_data_dir, local_starfusion_dir)):
        util.check_call(['tar', 'xzf', args.starfusion_targz, '-C', args.starfusion_data_dir])
        util.check_call(['make', '-C', os.path.join(args.starfusion_data_dir, local_starfusion_dir)])

    match = re.match(r'.*/([^/]+)\.tar\.gz$', args.starfusion_plug_n_play_targz)
    assert match
    local_plugnplay_dir = match[1]
    if not os.path.exists(args.starfusion_data_dir + local_plugnplay_dir):
        util.check_call(['tar', 'xzf', args.starfusion_plug_n_play_targz, '-C', args.starfusion_data_dir])

    result_dir = args.starfusion_result_dir + '/' + os.path.basename(sample_name + '-starfusion')
    logging.info('Start starfusion benchmark: %s', result_dir)
    try:
        os.makedirs(result_dir, 0o755)
    except:
        logging.error("mkdir %s failed", result_dir)

    cached_r1 = ",".join([args.cache_dir + '/' + os.path.basename(fp.r1) for fp in cached_file_pairs])
    cached_r2 = ",".join([args.cache_dir + '/' + os.path.basename(fp.r2) for fp in cached_file_pairs])

    starfusion_args = ['docker', 'run']
    mounted: Set[str] = set()
    for dir in [args.starfusion_data_dir, args.starfusion_result_dir, args.cache_dir]:
        if dir not in mounted:
            mounted.add(dir)
            starfusion_args += ['-v', f'{dir}:{dir}']
    starfusion_args += [
        '--rm',
        'trinityctat/ctatfusion',
        os.path.join(args.starfusion_data_dir, local_starfusion_dir, '/STAR-Fusion'),
        '--left_fq', cached_r1,
        '--right_fq', cached_r2,
        '--CPU', '56',
        '--genome_lib_dir', os.path.join(args.starfusion_data_dir, local_plugnplay_dir,'/ctat_genome_lib_build_dir'),
        '-O', result_dir,
        '--FusionInspector', 'validate']
    try:
        util.check_call(starfusion_args)
    except Exception as e:
        logging.error("Starfusion failed (ignoring): %s", e)
    logging.info('Finished starfusion benchmark: %s', result_dir)

def run_af4(sample_name: str, cached_file_pairs: List[FilePair], args: Any):
    ref_path = "s3://grail-pecan/resources/annotation/gencode.v26.whole_genes.fa"
    cosmic_fusion_path = "s3://grail-jkim/tmp/all_pair_art_lod_gpair_merged.txt"
    util.s3_cache_files([ref_path, cosmic_fusion_path], args.cache_dir)

    cached_r1 = ",".join([args.cache_dir + '/' + os.path.basename(fp.r1) for fp in cached_file_pairs])
    cached_r2 = ",".join([args.cache_dir + '/' + os.path.basename(fp.r2) for fp in cached_file_pairs])

    cached_ref = args.cache_dir + '/' + os.path.basename(ref_path)
    cached_cosmic_fusion = args.cache_dir + '/' + os.path.basename(cosmic_fusion_path)

    for mode in ['denovo', 'targeted']:
        result_dir = args.af4_result_dir + '/' + os.path.basename(sample_name + '-' + mode)
        if os.path.exists(result_dir + "/filtered.fa"):
            logging.info('Skipping benchmark: %s', result_dir)
            continue
        logging.info('Start af4 benchmark: %s', result_dir)
        try:
            os.makedirs(result_dir, 0o755)
        except:
            logging.error("mkdir %s failed", result_dir)
        af4_args = [str(util.af4_path()),
                    f'-log_dir={result_dir}',
                    f'-pprof=:12345',
                    f'-mutex-profile-rate=1000',
                    f'-block-profile-rate=1000',
                    f'-umi-in-read',
                    f'-r1={cached_r1}',
                    f'-r2={cached_r2}',
                    f'-fasta-output={result_dir}/all.fa',
                    f'-filtered-output={result_dir}/filtered.fa',
                    f'-transcript=' + args.cache_dir + '/gencode.v26.250padded_separate_jns_transcripts_parsed_no_mt_no_overlap_no_pary_no_versioned.fa']
        if mode == 'targeted':
            af4_args.append(f'-cosmic-fusion={cached_cosmic_fusion}')
        util.check_call(af4_args)
        logging.info('Finished benchmark: %s', result_dir)
        logging.info("Runtime stats: %s", util.run_stats(Path(result_dir)))

def main() -> None:
    logging.basicConfig(
        level=logging.DEBUG,
        format="%(asctime)s:%(levelname)s: %(message)s")
    p = argparse.ArgumentParser()
    p.add_argument('--cache_dir', default='/scratch-nvme/cache_tmp',
                   help='Benchmark cache dir')
    p.add_argument('--af4_result_dir', default='/scratch-sata/af4_benchmark_results',
                   help='Benchmark result dir')
    p.add_argument('--strafusion_result_dir', default='/scratch-nvme/af4_benchmark_results',
                   help='Benchmark result dir')
    p.add_argument('--starfusion_data_dir', default='/scratch-nvme/starfusion',
                   help='Directory for expanding starfusion plug-n-play files')
    p.add_argument('--run', action='append', choices=['af4', 'starfusion'],
                   help='List of systems to run. If unset, run all the configured systems')
    p.add_argument('--starfusion_plug_n_play_targz', default=os.environ['HOME'] + '/GRCh38_v27_CTAT_lib_Feb092018.plug-n-play.tar.gz',
                   help='Tar.gz file of starfusion plug-n-play file. https://github.com/STAR-Fusion/STAR-Fusion/wiki#data-resources-required')
    p.add_argument('--starfusion_targz', default=os.environ['HOME'] + '/STAR-Fusion-v1.5.0.FULL.tar.gz',
                   help='Tar.gz file of starfusion source package. https://github.com/STAR-Fusion/STAR-Fusion/wiki#data-resources-required')

    args = p.parse_args()
    if not args.run:
        args.run = ['af4', 'starfusion']

    for sample_name, file_pairs in SAMPLES.items():
        fastq_files: List[str] = []
        cached_file_pairs: List[FilePair] = []
        for fp in file_pairs:
            assert fp.r1.replace("R1", "R2"), fp.r2
            fastq_files += [fp.r1, fp.r2]
            cached_file_pairs.append(FilePair(r1=args.cache_dir + '/' + os.path.basename(fp.r1),
                                              r2=args.cache_dir + '/' + os.path.basename(fp.r2)))
        util.s3_cache_files(fastq_files, args.cache_dir)

        if 'af4' in args.run:
            run_af4(sample_name, cached_file_pairs, args)
        if 'starfusion' in args.run:
            run_starfusion(sample_name, cached_file_pairs, args)

        for fp in cached_file_pairs:
            try:
                #os.remove(fp.r1)
                #os.remove(fp.r2)
                pass
            except:
                logging.error("failed to remove %s", fp)
main()
