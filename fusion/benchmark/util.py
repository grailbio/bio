import logging
import glob
import os
import subprocess
import re
import itertools
from pathlib import Path
from typing import Iterable, List, NamedTuple, Optional

# The local directory to download the S3 files.
DEFAULT_CACHE_DIR = '/scratch-nvme/cache_tmp'

# The directory that stores benchmark results.
DEFAULT_RESULT_DIR = '/scratch-nvme/af4_results'

# FASTQPair is a pair of fastq pathnames, one for R1, the other for R2
FASTQPair = NamedTuple('FASTQPair', [('r1', str), ('r2', str)])

BENCHMARK_DATA_DIR = 's3://grail-publications/2019-ISMB'

ORG_REFERENCE_PATHS = [
    's3://grail-arao/fusion_references/gencode.v26.250padded_separate_jns_transcripts_parsed_no_mt_no_overlap_no_pary_no_versioned.fa',
    's3://grail-pecan/resources/annotation/gencode.v26.whole_genes.fa',
    's3://grail-jkim/tmp/all_pair_art_lod_gpair_merged.txt',
    's3://grail-arao/fusion_references/liu_gpair.txt']
REFERENCE_DIR = f'{BENCHMARK_DATA_DIR}/references'
REFERENCE_PATHS = [f'{REFERENCE_DIR}/{os.path.basename(path)}' for path in ORG_REFERENCE_PATHS]

# Describes a Liu et al simulated DNA fusion fastq filepair
SimulatedSample = NamedTuple('SimulatedSample',
                             [('n', str), # read length '50', '75', etc
                              ('coverage', str), # '5x', '20x', etc.
                              ('path', FASTQPair)]) # FASTQ files

# Original Grail-internal location.
ORG_SIMULATED_SAMPLES = [
    SimulatedSample(n=n,
                    coverage=x,
                    path=FASTQPair(
                        r1=f's3://grail-arao/liu_et_al_data/read{n}/{x}/{x}_FRG500_SD50_R{n}_1.fastq.gz',
                        r2=f's3://grail-arao/liu_et_al_data/read{n}/{x}/{x}_FRG500_SD50_R{n}_2.fastq.gz'))
    for n, x in itertools.product(['50', '75', '100'],
                                  ['5X', '50X', '20X', '100X', '200X'])]

# Public location.
SIMULATED_BENCHMARK_DIR = f'{BENCHMARK_DATA_DIR}/simulated_benchmark'
SIMULATED_SAMPLES = [
    SimulatedSample(n=fq.n,
                    coverage=fq.coverage,
                    path=FASTQPair(
                        r1=f'{SIMULATED_BENCHMARK_DIR}/{os.path.basename(fq.path.r1)}',
                        r2=f'{SIMULATED_BENCHMARK_DIR}/{os.path.basename(fq.path.r2)}'))
    for fq in ORG_SIMULATED_SAMPLES]

TitrationSample = NamedTuple('TitrationSample',
                             [('name', str),
                              ('paths', List[FASTQPair])])

# Files used in the titration benchmarks (grail-internal location)
ORG_TITRATION_SAMPLES = [
    TitrationSample(name='170206_ARTLoD_B1_01rerun', paths=[
        FASTQPair(r1='s3://grail-clinical-fastq/MenloPark/170214_E00543_0033_AHCLLKALXX/2738376149/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_01rerun/170206_ARTLoD_B1_01rerun_S1_L001_R1_001.fastq.gz',
                  r2='s3://grail-clinical-fastq/MenloPark/170214_E00543_0033_AHCLLKALXX/2738376149/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_01rerun/170206_ARTLoD_B1_01rerun_S1_L001_R2_001.fastq.gz'),
        FASTQPair(r1='s3://grail-clinical-fastq/MenloPark/170214_E00543_0033_AHCLLKALXX/2738376149/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_01rerun/170206_ARTLoD_B1_01rerun_S1_L002_R1_001.fastq.gz',
                  r2='s3://grail-clinical-fastq/MenloPark/170214_E00543_0033_AHCLLKALXX/2738376149/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_01rerun/170206_ARTLoD_B1_01rerun_S1_L002_R2_001.fastq.gz'),
        FASTQPair(r1='s3://grail-clinical-fastq/MenloPark/170214_E00543_0033_AHCLLKALXX/2738376149/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_01rerun/170206_ARTLoD_B1_01rerun_S1_L003_R1_001.fastq.gz',
                  r2='s3://grail-clinical-fastq/MenloPark/170214_E00543_0033_AHCLLKALXX/2738376149/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_01rerun/170206_ARTLoD_B1_01rerun_S1_L003_R2_001.fastq.gz'),
        FASTQPair(r1='s3://grail-clinical-fastq/MenloPark/170214_E00543_0033_AHCLLKALXX/2738376149/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_01rerun/170206_ARTLoD_B1_01rerun_S1_L004_R1_001.fastq.gz',
                  r2='s3://grail-clinical-fastq/MenloPark/170214_E00543_0033_AHCLLKALXX/2738376149/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_01rerun/170206_ARTLoD_B1_01rerun_S1_L004_R2_001.fastq.gz'),
        FASTQPair(r1='s3://grail-clinical-fastq/MenloPark/170214_E00543_0033_AHCLLKALXX/2738376149/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_01rerun/170206_ARTLoD_B1_01rerun_S1_L005_R1_001.fastq.gz',
                  r2='s3://grail-clinical-fastq/MenloPark/170214_E00543_0033_AHCLLKALXX/2738376149/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_01rerun/170206_ARTLoD_B1_01rerun_S1_L005_R2_001.fastq.gz'),
        FASTQPair(r1='s3://grail-clinical-fastq/MenloPark/170214_E00543_0033_AHCLLKALXX/2738376149/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_01rerun/170206_ARTLoD_B1_01rerun_S1_L006_R1_001.fastq.gz',
                  r2='s3://grail-clinical-fastq/MenloPark/170214_E00543_0033_AHCLLKALXX/2738376149/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_01rerun/170206_ARTLoD_B1_01rerun_S1_L006_R2_001.fastq.gz'),
        FASTQPair(r1='s3://grail-clinical-fastq/MenloPark/170214_E00543_0033_AHCLLKALXX/2738376149/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_01rerun/170206_ARTLoD_B1_01rerun_S1_L007_R1_001.fastq.gz',
                  r2='s3://grail-clinical-fastq/MenloPark/170214_E00543_0033_AHCLLKALXX/2738376149/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_01rerun/170206_ARTLoD_B1_01rerun_S1_L007_R2_001.fastq.gz'),
        FASTQPair(r1='s3://grail-clinical-fastq/MenloPark/170214_E00543_0033_AHCLLKALXX/2738376149/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_01rerun/170206_ARTLoD_B1_01rerun_S1_L008_R1_001.fastq.gz',
                  r2='s3://grail-clinical-fastq/MenloPark/170214_E00543_0033_AHCLLKALXX/2738376149/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_01rerun/170206_ARTLoD_B1_01rerun_S1_L008_R2_001.fastq.gz'),
    ]),
    TitrationSample(name='170206_ARTLoD_B1_02rerun', paths=[
        FASTQPair(r1='s3://grail-clinical-fastq/MenloPark/170214_E00543_0034_BHCLJVALXX/2291556533/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_02rerun/170206_ARTLoD_B1_02rerun_S1_L001_R1_001.fastq.gz',
                  r2='s3://grail-clinical-fastq/MenloPark/170214_E00543_0034_BHCLJVALXX/2291556533/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_02rerun/170206_ARTLoD_B1_02rerun_S1_L001_R2_001.fastq.gz'),
        FASTQPair(r1='s3://grail-clinical-fastq/MenloPark/170214_E00543_0034_BHCLJVALXX/2291556533/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_02rerun/170206_ARTLoD_B1_02rerun_S1_L002_R1_001.fastq.gz',
                  r2='s3://grail-clinical-fastq/MenloPark/170214_E00543_0034_BHCLJVALXX/2291556533/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_02rerun/170206_ARTLoD_B1_02rerun_S1_L002_R2_001.fastq.gz'),
        FASTQPair(r1='s3://grail-clinical-fastq/MenloPark/170214_E00543_0034_BHCLJVALXX/2291556533/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_02rerun/170206_ARTLoD_B1_02rerun_S1_L003_R1_001.fastq.gz',
                  r2='s3://grail-clinical-fastq/MenloPark/170214_E00543_0034_BHCLJVALXX/2291556533/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_02rerun/170206_ARTLoD_B1_02rerun_S1_L003_R2_001.fastq.gz'),
        FASTQPair(r1='s3://grail-clinical-fastq/MenloPark/170214_E00543_0034_BHCLJVALXX/2291556533/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_02rerun/170206_ARTLoD_B1_02rerun_S1_L004_R1_001.fastq.gz',
                  r2='s3://grail-clinical-fastq/MenloPark/170214_E00543_0034_BHCLJVALXX/2291556533/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_02rerun/170206_ARTLoD_B1_02rerun_S1_L004_R2_001.fastq.gz'),
        FASTQPair(r1='s3://grail-clinical-fastq/MenloPark/170214_E00543_0034_BHCLJVALXX/2291556533/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_02rerun/170206_ARTLoD_B1_02rerun_S1_L005_R1_001.fastq.gz',
                  r2='s3://grail-clinical-fastq/MenloPark/170214_E00543_0034_BHCLJVALXX/2291556533/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_02rerun/170206_ARTLoD_B1_02rerun_S1_L005_R2_001.fastq.gz'),
        FASTQPair(r1='s3://grail-clinical-fastq/MenloPark/170214_E00543_0034_BHCLJVALXX/2291556533/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_02rerun/170206_ARTLoD_B1_02rerun_S1_L006_R1_001.fastq.gz',
                  r2='s3://grail-clinical-fastq/MenloPark/170214_E00543_0034_BHCLJVALXX/2291556533/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_02rerun/170206_ARTLoD_B1_02rerun_S1_L006_R2_001.fastq.gz'),
        FASTQPair(r1='s3://grail-clinical-fastq/MenloPark/170214_E00543_0034_BHCLJVALXX/2291556533/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_02rerun/170206_ARTLoD_B1_02rerun_S1_L007_R1_001.fastq.gz',
                  r2='s3://grail-clinical-fastq/MenloPark/170214_E00543_0034_BHCLJVALXX/2291556533/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_02rerun/170206_ARTLoD_B1_02rerun_S1_L007_R2_001.fastq.gz'),
        FASTQPair(r1='s3://grail-clinical-fastq/MenloPark/170214_E00543_0034_BHCLJVALXX/2291556533/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_02rerun/170206_ARTLoD_B1_02rerun_S1_L008_R1_001.fastq.gz',
                  r2='s3://grail-clinical-fastq/MenloPark/170214_E00543_0034_BHCLJVALXX/2291556533/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_02rerun/170206_ARTLoD_B1_02rerun_S1_L008_R2_001.fastq.gz'),
    ]),
    # TitrationSample(name='170206_ARTLoD_B1_03rerun', paths=[
    #     FASTQPair(r1='s3://grail-clinical-fastq/MenloPark/170214_E00564_0033_AHCWC7ALXX/3447259045/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_03rerun/170206_ARTLoD_B1_03rerun_S1_L001_R1_001.fastq.gz',
    #               r2='s3://grail-clinical-fastq/MenloPark/170214_E00564_0033_AHCWC7ALXX/3447259045/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_03rerun/170206_ARTLoD_B1_03rerun_S1_L001_R2_001.fastq.gz'),
    #     FASTQPair(r1='s3://grail-clinical-fastq/MenloPark/170214_E00564_0033_AHCWC7ALXX/3447259045/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_03rerun/170206_ARTLoD_B1_03rerun_S1_L002_R1_001.fastq.gz',
    #               r2='s3://grail-clinical-fastq/MenloPark/170214_E00564_0033_AHCWC7ALXX/3447259045/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_03rerun/170206_ARTLoD_B1_03rerun_S1_L002_R2_001.fastq.gz'),
    #     FASTQPair(r1='s3://grail-clinical-fastq/MenloPark/170214_E00564_0033_AHCWC7ALXX/3447259045/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_03rerun/170206_ARTLoD_B1_03rerun_S1_L003_R1_001.fastq.gz',
    #               r2='s3://grail-clinical-fastq/MenloPark/170214_E00564_0033_AHCWC7ALXX/3447259045/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_03rerun/170206_ARTLoD_B1_03rerun_S1_L003_R2_001.fastq.gz'),
    #     FASTQPair(r1='s3://grail-clinical-fastq/MenloPark/170214_E00564_0033_AHCWC7ALXX/3447259045/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_03rerun/170206_ARTLoD_B1_03rerun_S1_L004_R1_001.fastq.gz',
    #               r2='s3://grail-clinical-fastq/MenloPark/170214_E00564_0033_AHCWC7ALXX/3447259045/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_03rerun/170206_ARTLoD_B1_03rerun_S1_L004_R2_001.fastq.gz'),
    #     FASTQPair(r1='s3://grail-clinical-fastq/MenloPark/170214_E00564_0033_AHCWC7ALXX/3447259045/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_03rerun/170206_ARTLoD_B1_03rerun_S1_L005_R1_001.fastq.gz',
    #               r2='s3://grail-clinical-fastq/MenloPark/170214_E00564_0033_AHCWC7ALXX/3447259045/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_03rerun/170206_ARTLoD_B1_03rerun_S1_L005_R2_001.fastq.gz'),
    #     FASTQPair(r1='s3://grail-clinical-fastq/MenloPark/170214_E00564_0033_AHCWC7ALXX/3447259045/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_03rerun/170206_ARTLoD_B1_03rerun_S1_L006_R1_001.fastq.gz',
    #               r2='s3://grail-clinical-fastq/MenloPark/170214_E00564_0033_AHCWC7ALXX/3447259045/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_03rerun/170206_ARTLoD_B1_03rerun_S1_L006_R2_001.fastq.gz'),
    #     FASTQPair(r1='s3://grail-clinical-fastq/MenloPark/170214_E00564_0033_AHCWC7ALXX/3447259045/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_03rerun/170206_ARTLoD_B1_03rerun_S1_L007_R1_001.fastq.gz',
    #               r2='s3://grail-clinical-fastq/MenloPark/170214_E00564_0033_AHCWC7ALXX/3447259045/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_03rerun/170206_ARTLoD_B1_03rerun_S1_L007_R2_001.fastq.gz'),
    #     FASTQPair(r1='s3://grail-clinical-fastq/MenloPark/170214_E00564_0033_AHCWC7ALXX/3447259045/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_03rerun/170206_ARTLoD_B1_03rerun_S1_L008_R1_001.fastq.gz',
    #               r2='s3://grail-clinical-fastq/MenloPark/170214_E00564_0033_AHCWC7ALXX/3447259045/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_03rerun/170206_ARTLoD_B1_03rerun_S1_L008_R2_001.fastq.gz'),
    # ]),
    TitrationSample(name='170206_ARTLoD_B1_04rerun', paths=[
        FASTQPair(r1='s3://grail-clinical-fastq/MenloPark/170214_E00564_0034_BHCL2WALXX/2255989125/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_04rerun/170206_ARTLoD_B1_04rerun_S1_L001_R1_001.fastq.gz',
                  r2='s3://grail-clinical-fastq/MenloPark/170214_E00564_0034_BHCL2WALXX/2255989125/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_04rerun/170206_ARTLoD_B1_04rerun_S1_L001_R2_001.fastq.gz'),
        FASTQPair(r1='s3://grail-clinical-fastq/MenloPark/170214_E00564_0034_BHCL2WALXX/2255989125/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_04rerun/170206_ARTLoD_B1_04rerun_S1_L002_R1_001.fastq.gz',
                  r2='s3://grail-clinical-fastq/MenloPark/170214_E00564_0034_BHCL2WALXX/2255989125/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_04rerun/170206_ARTLoD_B1_04rerun_S1_L002_R2_001.fastq.gz'),
        FASTQPair(r1='s3://grail-clinical-fastq/MenloPark/170214_E00564_0034_BHCL2WALXX/2255989125/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_04rerun/170206_ARTLoD_B1_04rerun_S1_L003_R1_001.fastq.gz',
                  r2='s3://grail-clinical-fastq/MenloPark/170214_E00564_0034_BHCL2WALXX/2255989125/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_04rerun/170206_ARTLoD_B1_04rerun_S1_L003_R2_001.fastq.gz'),
        FASTQPair(r1='s3://grail-clinical-fastq/MenloPark/170214_E00564_0034_BHCL2WALXX/2255989125/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_04rerun/170206_ARTLoD_B1_04rerun_S1_L004_R1_001.fastq.gz',
                  r2='s3://grail-clinical-fastq/MenloPark/170214_E00564_0034_BHCL2WALXX/2255989125/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_04rerun/170206_ARTLoD_B1_04rerun_S1_L004_R2_001.fastq.gz'),
        FASTQPair(r1='s3://grail-clinical-fastq/MenloPark/170214_E00564_0034_BHCL2WALXX/2255989125/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_04rerun/170206_ARTLoD_B1_04rerun_S1_L005_R1_001.fastq.gz',
                  r2='s3://grail-clinical-fastq/MenloPark/170214_E00564_0034_BHCL2WALXX/2255989125/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_04rerun/170206_ARTLoD_B1_04rerun_S1_L005_R2_001.fastq.gz'),
        FASTQPair(r1='s3://grail-clinical-fastq/MenloPark/170214_E00564_0034_BHCL2WALXX/2255989125/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_04rerun/170206_ARTLoD_B1_04rerun_S1_L006_R1_001.fastq.gz',
                  r2='s3://grail-clinical-fastq/MenloPark/170214_E00564_0034_BHCL2WALXX/2255989125/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_04rerun/170206_ARTLoD_B1_04rerun_S1_L006_R2_001.fastq.gz'),
        FASTQPair(r1='s3://grail-clinical-fastq/MenloPark/170214_E00564_0034_BHCL2WALXX/2255989125/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_04rerun/170206_ARTLoD_B1_04rerun_S1_L007_R1_001.fastq.gz',
                  r2='s3://grail-clinical-fastq/MenloPark/170214_E00564_0034_BHCL2WALXX/2255989125/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_04rerun/170206_ARTLoD_B1_04rerun_S1_L007_R2_001.fastq.gz'),
        FASTQPair(r1='s3://grail-clinical-fastq/MenloPark/170214_E00564_0034_BHCL2WALXX/2255989125/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_04rerun/170206_ARTLoD_B1_04rerun_S1_L008_R1_001.fastq.gz',
                  r2='s3://grail-clinical-fastq/MenloPark/170214_E00564_0034_BHCL2WALXX/2255989125/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_04rerun/170206_ARTLoD_B1_04rerun_S1_L008_R2_001.fastq.gz'),
    ]),
    # TitrationSample(name='170206_ARTLoD_B1_05rerun', paths=[
    #     FASTQPair(r1='s3://grail-clinical-fastq/MenloPark/170214_E00543_0033_AHCLLKALXX/2738376149/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_05rerun/170206_ARTLoD_B1_05rerun_S2_L001_R1_001.fastq.gz',
    #               r2='s3://grail-clinical-fastq/MenloPark/170214_E00543_0033_AHCLLKALXX/2738376149/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_05rerun/170206_ARTLoD_B1_05rerun_S2_L001_R2_001.fastq.gz'),
    #     FASTQPair(r1='s3://grail-clinical-fastq/MenloPark/170214_E00543_0033_AHCLLKALXX/2738376149/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_05rerun/170206_ARTLoD_B1_05rerun_S2_L002_R1_001.fastq.gz',
    #               r2='s3://grail-clinical-fastq/MenloPark/170214_E00543_0033_AHCLLKALXX/2738376149/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_05rerun/170206_ARTLoD_B1_05rerun_S2_L002_R2_001.fastq.gz'),
    #     FASTQPair(r1='s3://grail-clinical-fastq/MenloPark/170214_E00543_0033_AHCLLKALXX/2738376149/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_05rerun/170206_ARTLoD_B1_05rerun_S2_L003_R1_001.fastq.gz',
    #               r2='s3://grail-clinical-fastq/MenloPark/170214_E00543_0033_AHCLLKALXX/2738376149/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_05rerun/170206_ARTLoD_B1_05rerun_S2_L003_R2_001.fastq.gz'),
    #     FASTQPair(r1='s3://grail-clinical-fastq/MenloPark/170214_E00543_0033_AHCLLKALXX/2738376149/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_05rerun/170206_ARTLoD_B1_05rerun_S2_L004_R1_001.fastq.gz',
    #               r2='s3://grail-clinical-fastq/MenloPark/170214_E00543_0033_AHCLLKALXX/2738376149/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_05rerun/170206_ARTLoD_B1_05rerun_S2_L004_R2_001.fastq.gz'),
    #     FASTQPair(r1='s3://grail-clinical-fastq/MenloPark/170214_E00543_0033_AHCLLKALXX/2738376149/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_05rerun/170206_ARTLoD_B1_05rerun_S2_L005_R1_001.fastq.gz',
    #               r2='s3://grail-clinical-fastq/MenloPark/170214_E00543_0033_AHCLLKALXX/2738376149/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_05rerun/170206_ARTLoD_B1_05rerun_S2_L005_R2_001.fastq.gz'),
    #     FASTQPair(r1='s3://grail-clinical-fastq/MenloPark/170214_E00543_0033_AHCLLKALXX/2738376149/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_05rerun/170206_ARTLoD_B1_05rerun_S2_L006_R1_001.fastq.gz',
    #               r2='s3://grail-clinical-fastq/MenloPark/170214_E00543_0033_AHCLLKALXX/2738376149/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_05rerun/170206_ARTLoD_B1_05rerun_S2_L006_R2_001.fastq.gz'),
    #     FASTQPair(r1='s3://grail-clinical-fastq/MenloPark/170214_E00543_0033_AHCLLKALXX/2738376149/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_05rerun/170206_ARTLoD_B1_05rerun_S2_L007_R1_001.fastq.gz',
    #               r2='s3://grail-clinical-fastq/MenloPark/170214_E00543_0033_AHCLLKALXX/2738376149/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_05rerun/170206_ARTLoD_B1_05rerun_S2_L007_R2_001.fastq.gz'),
    #     FASTQPair(r1='s3://grail-clinical-fastq/MenloPark/170214_E00543_0033_AHCLLKALXX/2738376149/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_05rerun/170206_ARTLoD_B1_05rerun_S2_L008_R1_001.fastq.gz',
    #               r2='s3://grail-clinical-fastq/MenloPark/170214_E00543_0033_AHCLLKALXX/2738376149/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_05rerun/170206_ARTLoD_B1_05rerun_S2_L008_R2_001.fastq.gz'),
    # ]),
    # TitrationSample(name='170206_ARTLoD_B1_06rerun', paths=[
    #     FASTQPair(r1='s3://grail-clinical-fastq/MenloPark/170214_E00543_0034_BHCLJVALXX/2291556533/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_06rerun/170206_ARTLoD_B1_06rerun_S2_L001_R1_001.fastq.gz',
    #               r2='s3://grail-clinical-fastq/MenloPark/170214_E00543_0034_BHCLJVALXX/2291556533/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_06rerun/170206_ARTLoD_B1_06rerun_S2_L001_R2_001.fastq.gz'),
    #     FASTQPair(r1='s3://grail-clinical-fastq/MenloPark/170214_E00543_0034_BHCLJVALXX/2291556533/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_06rerun/170206_ARTLoD_B1_06rerun_S2_L002_R1_001.fastq.gz',
    #               r2='s3://grail-clinical-fastq/MenloPark/170214_E00543_0034_BHCLJVALXX/2291556533/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_06rerun/170206_ARTLoD_B1_06rerun_S2_L002_R2_001.fastq.gz'),
    #     FASTQPair(r1='s3://grail-clinical-fastq/MenloPark/170214_E00543_0034_BHCLJVALXX/2291556533/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_06rerun/170206_ARTLoD_B1_06rerun_S2_L003_R1_001.fastq.gz',
    #               r2='s3://grail-clinical-fastq/MenloPark/170214_E00543_0034_BHCLJVALXX/2291556533/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_06rerun/170206_ARTLoD_B1_06rerun_S2_L003_R2_001.fastq.gz'),
    #     FASTQPair(r1='s3://grail-clinical-fastq/MenloPark/170214_E00543_0034_BHCLJVALXX/2291556533/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_06rerun/170206_ARTLoD_B1_06rerun_S2_L004_R1_001.fastq.gz',
    #               r2='s3://grail-clinical-fastq/MenloPark/170214_E00543_0034_BHCLJVALXX/2291556533/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_06rerun/170206_ARTLoD_B1_06rerun_S2_L004_R2_001.fastq.gz'),
    #     FASTQPair(r1='s3://grail-clinical-fastq/MenloPark/170214_E00543_0034_BHCLJVALXX/2291556533/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_06rerun/170206_ARTLoD_B1_06rerun_S2_L005_R1_001.fastq.gz',
    #               r2='s3://grail-clinical-fastq/MenloPark/170214_E00543_0034_BHCLJVALXX/2291556533/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_06rerun/170206_ARTLoD_B1_06rerun_S2_L005_R2_001.fastq.gz'),
    #     FASTQPair(r1='s3://grail-clinical-fastq/MenloPark/170214_E00543_0034_BHCLJVALXX/2291556533/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_06rerun/170206_ARTLoD_B1_06rerun_S2_L006_R1_001.fastq.gz',
    #               r2='s3://grail-clinical-fastq/MenloPark/170214_E00543_0034_BHCLJVALXX/2291556533/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_06rerun/170206_ARTLoD_B1_06rerun_S2_L006_R2_001.fastq.gz'),
    #     FASTQPair(r1='s3://grail-clinical-fastq/MenloPark/170214_E00543_0034_BHCLJVALXX/2291556533/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_06rerun/170206_ARTLoD_B1_06rerun_S2_L007_R1_001.fastq.gz',
    #               r2='s3://grail-clinical-fastq/MenloPark/170214_E00543_0034_BHCLJVALXX/2291556533/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_06rerun/170206_ARTLoD_B1_06rerun_S2_L007_R2_001.fastq.gz'),
    #     FASTQPair(r1='s3://grail-clinical-fastq/MenloPark/170214_E00543_0034_BHCLJVALXX/2291556533/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_06rerun/170206_ARTLoD_B1_06rerun_S2_L008_R1_001.fastq.gz',
    #               r2='s3://grail-clinical-fastq/MenloPark/170214_E00543_0034_BHCLJVALXX/2291556533/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_06rerun/170206_ARTLoD_B1_06rerun_S2_L008_R2_001.fastq.gz'),
    # ]),
    TitrationSample(name='170206_ARTLoD_B1_07rerun', paths=[
        FASTQPair(r1='s3://grail-clinical-fastq/MenloPark/170214_E00564_0033_AHCWC7ALXX/3447259045/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_07rerun/170206_ARTLoD_B1_07rerun_S2_L001_R1_001.fastq.gz',
                  r2='s3://grail-clinical-fastq/MenloPark/170214_E00564_0033_AHCWC7ALXX/3447259045/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_07rerun/170206_ARTLoD_B1_07rerun_S2_L001_R2_001.fastq.gz'),
        FASTQPair(r1='s3://grail-clinical-fastq/MenloPark/170214_E00564_0033_AHCWC7ALXX/3447259045/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_07rerun/170206_ARTLoD_B1_07rerun_S2_L002_R1_001.fastq.gz',
                  r2='s3://grail-clinical-fastq/MenloPark/170214_E00564_0033_AHCWC7ALXX/3447259045/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_07rerun/170206_ARTLoD_B1_07rerun_S2_L002_R2_001.fastq.gz'),
        FASTQPair(r1='s3://grail-clinical-fastq/MenloPark/170214_E00564_0033_AHCWC7ALXX/3447259045/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_07rerun/170206_ARTLoD_B1_07rerun_S2_L003_R1_001.fastq.gz',
                  r2='s3://grail-clinical-fastq/MenloPark/170214_E00564_0033_AHCWC7ALXX/3447259045/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_07rerun/170206_ARTLoD_B1_07rerun_S2_L003_R2_001.fastq.gz'),
        FASTQPair(r1='s3://grail-clinical-fastq/MenloPark/170214_E00564_0033_AHCWC7ALXX/3447259045/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_07rerun/170206_ARTLoD_B1_07rerun_S2_L004_R1_001.fastq.gz',
                  r2='s3://grail-clinical-fastq/MenloPark/170214_E00564_0033_AHCWC7ALXX/3447259045/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_07rerun/170206_ARTLoD_B1_07rerun_S2_L004_R2_001.fastq.gz'),
        FASTQPair(r1='s3://grail-clinical-fastq/MenloPark/170214_E00564_0033_AHCWC7ALXX/3447259045/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_07rerun/170206_ARTLoD_B1_07rerun_S2_L005_R1_001.fastq.gz',
                  r2='s3://grail-clinical-fastq/MenloPark/170214_E00564_0033_AHCWC7ALXX/3447259045/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_07rerun/170206_ARTLoD_B1_07rerun_S2_L005_R2_001.fastq.gz'),
        FASTQPair(r1='s3://grail-clinical-fastq/MenloPark/170214_E00564_0033_AHCWC7ALXX/3447259045/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_07rerun/170206_ARTLoD_B1_07rerun_S2_L006_R1_001.fastq.gz',
                  r2='s3://grail-clinical-fastq/MenloPark/170214_E00564_0033_AHCWC7ALXX/3447259045/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_07rerun/170206_ARTLoD_B1_07rerun_S2_L006_R2_001.fastq.gz'),
        FASTQPair(r1='s3://grail-clinical-fastq/MenloPark/170214_E00564_0033_AHCWC7ALXX/3447259045/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_07rerun/170206_ARTLoD_B1_07rerun_S2_L007_R1_001.fastq.gz',
                  r2='s3://grail-clinical-fastq/MenloPark/170214_E00564_0033_AHCWC7ALXX/3447259045/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_07rerun/170206_ARTLoD_B1_07rerun_S2_L007_R2_001.fastq.gz'),
        FASTQPair(r1='s3://grail-clinical-fastq/MenloPark/170214_E00564_0033_AHCWC7ALXX/3447259045/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_07rerun/170206_ARTLoD_B1_07rerun_S2_L008_R1_001.fastq.gz',
                  r2='s3://grail-clinical-fastq/MenloPark/170214_E00564_0033_AHCWC7ALXX/3447259045/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_07rerun/170206_ARTLoD_B1_07rerun_S2_L008_R2_001.fastq.gz'),
    ]),
    # TitrationSample(name='170206_ARTLoD_B1_08rerun', paths=[
    #     FASTQPair(r1='s3://grail-clinical-fastq/MenloPark/170214_E00564_0034_BHCL2WALXX/2255989125/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_08rerun/170206_ARTLoD_B1_08rerun_S2_L001_R1_001.fastq.gz',
    #               r2='s3://grail-clinical-fastq/MenloPark/170214_E00564_0034_BHCL2WALXX/2255989125/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_08rerun/170206_ARTLoD_B1_08rerun_S2_L001_R2_001.fastq.gz'),
    #     FASTQPair(r1='s3://grail-clinical-fastq/MenloPark/170214_E00564_0034_BHCL2WALXX/2255989125/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_08rerun/170206_ARTLoD_B1_08rerun_S2_L002_R1_001.fastq.gz',
    #               r2='s3://grail-clinical-fastq/MenloPark/170214_E00564_0034_BHCL2WALXX/2255989125/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_08rerun/170206_ARTLoD_B1_08rerun_S2_L002_R2_001.fastq.gz'),
    #     FASTQPair(r1='s3://grail-clinical-fastq/MenloPark/170214_E00564_0034_BHCL2WALXX/2255989125/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_08rerun/170206_ARTLoD_B1_08rerun_S2_L003_R1_001.fastq.gz',
    #               r2='s3://grail-clinical-fastq/MenloPark/170214_E00564_0034_BHCL2WALXX/2255989125/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_08rerun/170206_ARTLoD_B1_08rerun_S2_L003_R2_001.fastq.gz'),
    #     FASTQPair(r1='s3://grail-clinical-fastq/MenloPark/170214_E00564_0034_BHCL2WALXX/2255989125/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_08rerun/170206_ARTLoD_B1_08rerun_S2_L004_R1_001.fastq.gz',
    #               r2='s3://grail-clinical-fastq/MenloPark/170214_E00564_0034_BHCL2WALXX/2255989125/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_08rerun/170206_ARTLoD_B1_08rerun_S2_L004_R2_001.fastq.gz'),
    #     FASTQPair(r1='s3://grail-clinical-fastq/MenloPark/170214_E00564_0034_BHCL2WALXX/2255989125/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_08rerun/170206_ARTLoD_B1_08rerun_S2_L005_R1_001.fastq.gz',
    #               r2='s3://grail-clinical-fastq/MenloPark/170214_E00564_0034_BHCL2WALXX/2255989125/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_08rerun/170206_ARTLoD_B1_08rerun_S2_L005_R2_001.fastq.gz'),
    #     FASTQPair(r1='s3://grail-clinical-fastq/MenloPark/170214_E00564_0034_BHCL2WALXX/2255989125/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_08rerun/170206_ARTLoD_B1_08rerun_S2_L006_R1_001.fastq.gz',
    #               r2='s3://grail-clinical-fastq/MenloPark/170214_E00564_0034_BHCL2WALXX/2255989125/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_08rerun/170206_ARTLoD_B1_08rerun_S2_L006_R2_001.fastq.gz'),
    #     FASTQPair(r1='s3://grail-clinical-fastq/MenloPark/170214_E00564_0034_BHCL2WALXX/2255989125/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_08rerun/170206_ARTLoD_B1_08rerun_S2_L007_R1_001.fastq.gz',
    #               r2='s3://grail-clinical-fastq/MenloPark/170214_E00564_0034_BHCL2WALXX/2255989125/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_08rerun/170206_ARTLoD_B1_08rerun_S2_L007_R2_001.fastq.gz'),
    #     FASTQPair(r1='s3://grail-clinical-fastq/MenloPark/170214_E00564_0034_BHCL2WALXX/2255989125/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_08rerun/170206_ARTLoD_B1_08rerun_S2_L008_R1_001.fastq.gz',
    #               r2='s3://grail-clinical-fastq/MenloPark/170214_E00564_0034_BHCL2WALXX/2255989125/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_08rerun/170206_ARTLoD_B1_08rerun_S2_L008_R2_001.fastq.gz'),
    # ]),

    # TitrationSample(name='170206_ARTLoD_B1_09rerun', paths=[
    #     FASTQPair(r1='s3://grail-clinical-fastq/MenloPark/170214_E00543_0033_AHCLLKALXX/2738376149/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_09rerun/170206_ARTLoD_B1_09rerun_S3_L001_R1_001.fastq.gz',
    #               r2='s3://grail-clinical-fastq/MenloPark/170214_E00543_0033_AHCLLKALXX/2738376149/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_09rerun/170206_ARTLoD_B1_09rerun_S3_L001_R2_001.fastq.gz'),
    #     FASTQPair(r1='s3://grail-clinical-fastq/MenloPark/170214_E00543_0033_AHCLLKALXX/2738376149/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_09rerun/170206_ARTLoD_B1_09rerun_S3_L002_R1_001.fastq.gz',
    #               r2='s3://grail-clinical-fastq/MenloPark/170214_E00543_0033_AHCLLKALXX/2738376149/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_09rerun/170206_ARTLoD_B1_09rerun_S3_L002_R2_001.fastq.gz'),
    #     FASTQPair(r1='s3://grail-clinical-fastq/MenloPark/170214_E00543_0033_AHCLLKALXX/2738376149/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_09rerun/170206_ARTLoD_B1_09rerun_S3_L003_R1_001.fastq.gz',
    #               r2='s3://grail-clinical-fastq/MenloPark/170214_E00543_0033_AHCLLKALXX/2738376149/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_09rerun/170206_ARTLoD_B1_09rerun_S3_L003_R2_001.fastq.gz'),
    #     FASTQPair(r1='s3://grail-clinical-fastq/MenloPark/170214_E00543_0033_AHCLLKALXX/2738376149/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_09rerun/170206_ARTLoD_B1_09rerun_S3_L004_R1_001.fastq.gz',
    #               r2='s3://grail-clinical-fastq/MenloPark/170214_E00543_0033_AHCLLKALXX/2738376149/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_09rerun/170206_ARTLoD_B1_09rerun_S3_L004_R2_001.fastq.gz'),
    #     FASTQPair(r1='s3://grail-clinical-fastq/MenloPark/170214_E00543_0033_AHCLLKALXX/2738376149/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_09rerun/170206_ARTLoD_B1_09rerun_S3_L005_R1_001.fastq.gz',
    #               r2='s3://grail-clinical-fastq/MenloPark/170214_E00543_0033_AHCLLKALXX/2738376149/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_09rerun/170206_ARTLoD_B1_09rerun_S3_L005_R2_001.fastq.gz'),
    #     FASTQPair(r1='s3://grail-clinical-fastq/MenloPark/170214_E00543_0033_AHCLLKALXX/2738376149/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_09rerun/170206_ARTLoD_B1_09rerun_S3_L006_R1_001.fastq.gz',
    #               r2='s3://grail-clinical-fastq/MenloPark/170214_E00543_0033_AHCLLKALXX/2738376149/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_09rerun/170206_ARTLoD_B1_09rerun_S3_L006_R2_001.fastq.gz'),
    #     FASTQPair(r1='s3://grail-clinical-fastq/MenloPark/170214_E00543_0033_AHCLLKALXX/2738376149/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_09rerun/170206_ARTLoD_B1_09rerun_S3_L007_R1_001.fastq.gz',
    #               r2='s3://grail-clinical-fastq/MenloPark/170214_E00543_0033_AHCLLKALXX/2738376149/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_09rerun/170206_ARTLoD_B1_09rerun_S3_L007_R2_001.fastq.gz'),
    #     FASTQPair(r1='s3://grail-clinical-fastq/MenloPark/170214_E00543_0033_AHCLLKALXX/2738376149/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_09rerun/170206_ARTLoD_B1_09rerun_S3_L008_R1_001.fastq.gz',
    #               r2='s3://grail-clinical-fastq/MenloPark/170214_E00543_0033_AHCLLKALXX/2738376149/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_09rerun/170206_ARTLoD_B1_09rerun_S3_L008_R2_001.fastq.gz'),
    # ]),
    TitrationSample(name='170206_ARTLoD_B1_10rerun', paths=[
        FASTQPair(r1='s3://grail-clinical-fastq/MenloPark/170214_E00543_0034_BHCLJVALXX/2291556533/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_10rerun/170206_ARTLoD_B1_10rerun_S3_L001_R1_001.fastq.gz',
                  r2='s3://grail-clinical-fastq/MenloPark/170214_E00543_0034_BHCLJVALXX/2291556533/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_10rerun/170206_ARTLoD_B1_10rerun_S3_L001_R2_001.fastq.gz'),
        FASTQPair(r1='s3://grail-clinical-fastq/MenloPark/170214_E00543_0034_BHCLJVALXX/2291556533/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_10rerun/170206_ARTLoD_B1_10rerun_S3_L002_R1_001.fastq.gz',
                  r2='s3://grail-clinical-fastq/MenloPark/170214_E00543_0034_BHCLJVALXX/2291556533/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_10rerun/170206_ARTLoD_B1_10rerun_S3_L002_R2_001.fastq.gz'),
        FASTQPair(r1='s3://grail-clinical-fastq/MenloPark/170214_E00543_0034_BHCLJVALXX/2291556533/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_10rerun/170206_ARTLoD_B1_10rerun_S3_L003_R1_001.fastq.gz',
                  r2='s3://grail-clinical-fastq/MenloPark/170214_E00543_0034_BHCLJVALXX/2291556533/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_10rerun/170206_ARTLoD_B1_10rerun_S3_L003_R2_001.fastq.gz'),
        FASTQPair(r1='s3://grail-clinical-fastq/MenloPark/170214_E00543_0034_BHCLJVALXX/2291556533/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_10rerun/170206_ARTLoD_B1_10rerun_S3_L004_R1_001.fastq.gz',
                  r2='s3://grail-clinical-fastq/MenloPark/170214_E00543_0034_BHCLJVALXX/2291556533/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_10rerun/170206_ARTLoD_B1_10rerun_S3_L004_R2_001.fastq.gz'),
        FASTQPair(r1='s3://grail-clinical-fastq/MenloPark/170214_E00543_0034_BHCLJVALXX/2291556533/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_10rerun/170206_ARTLoD_B1_10rerun_S3_L005_R1_001.fastq.gz',
                  r2='s3://grail-clinical-fastq/MenloPark/170214_E00543_0034_BHCLJVALXX/2291556533/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_10rerun/170206_ARTLoD_B1_10rerun_S3_L005_R2_001.fastq.gz'),
        FASTQPair(r1='s3://grail-clinical-fastq/MenloPark/170214_E00543_0034_BHCLJVALXX/2291556533/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_10rerun/170206_ARTLoD_B1_10rerun_S3_L006_R1_001.fastq.gz',
                  r2='s3://grail-clinical-fastq/MenloPark/170214_E00543_0034_BHCLJVALXX/2291556533/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_10rerun/170206_ARTLoD_B1_10rerun_S3_L006_R2_001.fastq.gz'),
        FASTQPair(r1='s3://grail-clinical-fastq/MenloPark/170214_E00543_0034_BHCLJVALXX/2291556533/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_10rerun/170206_ARTLoD_B1_10rerun_S3_L007_R1_001.fastq.gz',
                  r2='s3://grail-clinical-fastq/MenloPark/170214_E00543_0034_BHCLJVALXX/2291556533/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_10rerun/170206_ARTLoD_B1_10rerun_S3_L007_R2_001.fastq.gz'),
        FASTQPair(r1='s3://grail-clinical-fastq/MenloPark/170214_E00543_0034_BHCLJVALXX/2291556533/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_10rerun/170206_ARTLoD_B1_10rerun_S3_L008_R1_001.fastq.gz',
                  r2='s3://grail-clinical-fastq/MenloPark/170214_E00543_0034_BHCLJVALXX/2291556533/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_10rerun/170206_ARTLoD_B1_10rerun_S3_L008_R2_001.fastq.gz'),
    ]),
    # TitrationSample(name='170206_ARTLoD_B1_11rerun', paths=[
    #     FASTQPair(r1='s3://grail-clinical-fastq/MenloPark/170214_E00564_0033_AHCWC7ALXX/3447259045/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_11rerun/170206_ARTLoD_B1_11rerun_S3_L001_R1_001.fastq.gz',
    #               r2='s3://grail-clinical-fastq/MenloPark/170214_E00564_0033_AHCWC7ALXX/3447259045/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_11rerun/170206_ARTLoD_B1_11rerun_S3_L001_R2_001.fastq.gz'),
    #     FASTQPair(r1='s3://grail-clinical-fastq/MenloPark/170214_E00564_0033_AHCWC7ALXX/3447259045/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_11rerun/170206_ARTLoD_B1_11rerun_S3_L002_R1_001.fastq.gz',
    #               r2='s3://grail-clinical-fastq/MenloPark/170214_E00564_0033_AHCWC7ALXX/3447259045/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_11rerun/170206_ARTLoD_B1_11rerun_S3_L002_R2_001.fastq.gz'),
    #     FASTQPair(r1='s3://grail-clinical-fastq/MenloPark/170214_E00564_0033_AHCWC7ALXX/3447259045/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_11rerun/170206_ARTLoD_B1_11rerun_S3_L003_R1_001.fastq.gz',
    #               r2='s3://grail-clinical-fastq/MenloPark/170214_E00564_0033_AHCWC7ALXX/3447259045/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_11rerun/170206_ARTLoD_B1_11rerun_S3_L003_R2_001.fastq.gz'),
    #     FASTQPair(r1='s3://grail-clinical-fastq/MenloPark/170214_E00564_0033_AHCWC7ALXX/3447259045/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_11rerun/170206_ARTLoD_B1_11rerun_S3_L004_R1_001.fastq.gz',
    #               r2='s3://grail-clinical-fastq/MenloPark/170214_E00564_0033_AHCWC7ALXX/3447259045/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_11rerun/170206_ARTLoD_B1_11rerun_S3_L004_R2_001.fastq.gz'),
    #     FASTQPair(r1='s3://grail-clinical-fastq/MenloPark/170214_E00564_0033_AHCWC7ALXX/3447259045/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_11rerun/170206_ARTLoD_B1_11rerun_S3_L005_R1_001.fastq.gz',
    #               r2='s3://grail-clinical-fastq/MenloPark/170214_E00564_0033_AHCWC7ALXX/3447259045/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_11rerun/170206_ARTLoD_B1_11rerun_S3_L005_R2_001.fastq.gz'),
    #     FASTQPair(r1='s3://grail-clinical-fastq/MenloPark/170214_E00564_0033_AHCWC7ALXX/3447259045/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_11rerun/170206_ARTLoD_B1_11rerun_S3_L006_R1_001.fastq.gz',
    #               r2='s3://grail-clinical-fastq/MenloPark/170214_E00564_0033_AHCWC7ALXX/3447259045/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_11rerun/170206_ARTLoD_B1_11rerun_S3_L006_R2_001.fastq.gz'),
    #     FASTQPair(r1='s3://grail-clinical-fastq/MenloPark/170214_E00564_0033_AHCWC7ALXX/3447259045/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_11rerun/170206_ARTLoD_B1_11rerun_S3_L007_R1_001.fastq.gz',
    #               r2='s3://grail-clinical-fastq/MenloPark/170214_E00564_0033_AHCWC7ALXX/3447259045/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_11rerun/170206_ARTLoD_B1_11rerun_S3_L007_R2_001.fastq.gz'),
    #     FASTQPair(r1='s3://grail-clinical-fastq/MenloPark/170214_E00564_0033_AHCWC7ALXX/3447259045/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_11rerun/170206_ARTLoD_B1_11rerun_S3_L008_R1_001.fastq.gz',
    #               r2='s3://grail-clinical-fastq/MenloPark/170214_E00564_0033_AHCWC7ALXX/3447259045/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_11rerun/170206_ARTLoD_B1_11rerun_S3_L008_R2_001.fastq.gz'),
    # ]),
    # TitrationSample(name='170206_ARTLoD_B1_12rerun', paths=[
    #     FASTQPair(r1='s3://grail-clinical-fastq/MenloPark/170214_E00543_0033_AHCLLKALXX/2738376149/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_12rerun/170206_ARTLoD_B1_12rerun_S4_L001_R1_001.fastq.gz',
    #               r2='s3://grail-clinical-fastq/MenloPark/170214_E00543_0033_AHCLLKALXX/2738376149/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_12rerun/170206_ARTLoD_B1_12rerun_S4_L001_R2_001.fastq.gz'),
    #     FASTQPair(r1='s3://grail-clinical-fastq/MenloPark/170214_E00543_0033_AHCLLKALXX/2738376149/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_12rerun/170206_ARTLoD_B1_12rerun_S4_L002_R1_001.fastq.gz',
    #               r2='s3://grail-clinical-fastq/MenloPark/170214_E00543_0033_AHCLLKALXX/2738376149/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_12rerun/170206_ARTLoD_B1_12rerun_S4_L002_R2_001.fastq.gz'),
    #     FASTQPair(r1='s3://grail-clinical-fastq/MenloPark/170214_E00543_0033_AHCLLKALXX/2738376149/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_12rerun/170206_ARTLoD_B1_12rerun_S4_L003_R1_001.fastq.gz',
    #               r2='s3://grail-clinical-fastq/MenloPark/170214_E00543_0033_AHCLLKALXX/2738376149/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_12rerun/170206_ARTLoD_B1_12rerun_S4_L003_R2_001.fastq.gz'),
    #     FASTQPair(r1='s3://grail-clinical-fastq/MenloPark/170214_E00543_0033_AHCLLKALXX/2738376149/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_12rerun/170206_ARTLoD_B1_12rerun_S4_L004_R1_001.fastq.gz',
    #               r2='s3://grail-clinical-fastq/MenloPark/170214_E00543_0033_AHCLLKALXX/2738376149/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_12rerun/170206_ARTLoD_B1_12rerun_S4_L004_R2_001.fastq.gz'),
    #     FASTQPair(r1='s3://grail-clinical-fastq/MenloPark/170214_E00543_0033_AHCLLKALXX/2738376149/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_12rerun/170206_ARTLoD_B1_12rerun_S4_L005_R1_001.fastq.gz',
    #               r2='s3://grail-clinical-fastq/MenloPark/170214_E00543_0033_AHCLLKALXX/2738376149/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_12rerun/170206_ARTLoD_B1_12rerun_S4_L005_R2_001.fastq.gz'),
    #     FASTQPair(r1='s3://grail-clinical-fastq/MenloPark/170214_E00543_0033_AHCLLKALXX/2738376149/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_12rerun/170206_ARTLoD_B1_12rerun_S4_L006_R1_001.fastq.gz',
    #               r2='s3://grail-clinical-fastq/MenloPark/170214_E00543_0033_AHCLLKALXX/2738376149/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_12rerun/170206_ARTLoD_B1_12rerun_S4_L006_R2_001.fastq.gz'),
    #     FASTQPair(r1='s3://grail-clinical-fastq/MenloPark/170214_E00543_0033_AHCLLKALXX/2738376149/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_12rerun/170206_ARTLoD_B1_12rerun_S4_L007_R1_001.fastq.gz',
    #               r2='s3://grail-clinical-fastq/MenloPark/170214_E00543_0033_AHCLLKALXX/2738376149/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_12rerun/170206_ARTLoD_B1_12rerun_S4_L007_R2_001.fastq.gz'),
    #     FASTQPair(r1='s3://grail-clinical-fastq/MenloPark/170214_E00543_0033_AHCLLKALXX/2738376149/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_12rerun/170206_ARTLoD_B1_12rerun_S4_L008_R1_001.fastq.gz',
    #               r2='s3://grail-clinical-fastq/MenloPark/170214_E00543_0033_AHCLLKALXX/2738376149/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_12rerun/170206_ARTLoD_B1_12rerun_S4_L008_R2_001.fastq.gz'),
    # ]),
    TitrationSample(name='170206_ARTLoD_B1_13rerun', paths=[
        FASTQPair(r1='s3://grail-clinical-fastq/MenloPark/170214_E00543_0034_BHCLJVALXX/2291556533/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_13rerun/170206_ARTLoD_B1_13rerun_S4_L001_R1_001.fastq.gz',
                  r2='s3://grail-clinical-fastq/MenloPark/170214_E00543_0034_BHCLJVALXX/2291556533/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_13rerun/170206_ARTLoD_B1_13rerun_S4_L001_R2_001.fastq.gz'),
        FASTQPair(r1='s3://grail-clinical-fastq/MenloPark/170214_E00543_0034_BHCLJVALXX/2291556533/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_13rerun/170206_ARTLoD_B1_13rerun_S4_L002_R1_001.fastq.gz',
                  r2='s3://grail-clinical-fastq/MenloPark/170214_E00543_0034_BHCLJVALXX/2291556533/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_13rerun/170206_ARTLoD_B1_13rerun_S4_L002_R2_001.fastq.gz'),
        FASTQPair(r1='s3://grail-clinical-fastq/MenloPark/170214_E00543_0034_BHCLJVALXX/2291556533/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_13rerun/170206_ARTLoD_B1_13rerun_S4_L003_R1_001.fastq.gz',
                  r2='s3://grail-clinical-fastq/MenloPark/170214_E00543_0034_BHCLJVALXX/2291556533/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_13rerun/170206_ARTLoD_B1_13rerun_S4_L003_R2_001.fastq.gz'),
        FASTQPair(r1='s3://grail-clinical-fastq/MenloPark/170214_E00543_0034_BHCLJVALXX/2291556533/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_13rerun/170206_ARTLoD_B1_13rerun_S4_L004_R1_001.fastq.gz',
                  r2='s3://grail-clinical-fastq/MenloPark/170214_E00543_0034_BHCLJVALXX/2291556533/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_13rerun/170206_ARTLoD_B1_13rerun_S4_L004_R2_001.fastq.gz'),
        FASTQPair(r1='s3://grail-clinical-fastq/MenloPark/170214_E00543_0034_BHCLJVALXX/2291556533/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_13rerun/170206_ARTLoD_B1_13rerun_S4_L005_R1_001.fastq.gz',
                  r2='s3://grail-clinical-fastq/MenloPark/170214_E00543_0034_BHCLJVALXX/2291556533/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_13rerun/170206_ARTLoD_B1_13rerun_S4_L005_R2_001.fastq.gz'),
        FASTQPair(r1='s3://grail-clinical-fastq/MenloPark/170214_E00543_0034_BHCLJVALXX/2291556533/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_13rerun/170206_ARTLoD_B1_13rerun_S4_L006_R1_001.fastq.gz',
                  r2='s3://grail-clinical-fastq/MenloPark/170214_E00543_0034_BHCLJVALXX/2291556533/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_13rerun/170206_ARTLoD_B1_13rerun_S4_L006_R2_001.fastq.gz'),
        FASTQPair(r1='s3://grail-clinical-fastq/MenloPark/170214_E00543_0034_BHCLJVALXX/2291556533/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_13rerun/170206_ARTLoD_B1_13rerun_S4_L007_R1_001.fastq.gz',
                  r2='s3://grail-clinical-fastq/MenloPark/170214_E00543_0034_BHCLJVALXX/2291556533/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_13rerun/170206_ARTLoD_B1_13rerun_S4_L007_R2_001.fastq.gz'),
        FASTQPair(r1='s3://grail-clinical-fastq/MenloPark/170214_E00543_0034_BHCLJVALXX/2291556533/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_13rerun/170206_ARTLoD_B1_13rerun_S4_L008_R1_001.fastq.gz',
                  r2='s3://grail-clinical-fastq/MenloPark/170214_E00543_0034_BHCLJVALXX/2291556533/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_13rerun/170206_ARTLoD_B1_13rerun_S4_L008_R2_001.fastq.gz'),
    ]),
    TitrationSample(name='170206_ARTLoD_B1_14rerun', paths=[
        FASTQPair(r1='s3://grail-clinical-fastq/MenloPark/170214_E00564_0033_AHCWC7ALXX/3447259045/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_14rerun/170206_ARTLoD_B1_14rerun_S4_L001_R1_001.fastq.gz',
                  r2='s3://grail-clinical-fastq/MenloPark/170214_E00564_0033_AHCWC7ALXX/3447259045/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_14rerun/170206_ARTLoD_B1_14rerun_S4_L001_R2_001.fastq.gz'),
        FASTQPair(r1='s3://grail-clinical-fastq/MenloPark/170214_E00564_0033_AHCWC7ALXX/3447259045/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_14rerun/170206_ARTLoD_B1_14rerun_S4_L002_R1_001.fastq.gz',
                  r2='s3://grail-clinical-fastq/MenloPark/170214_E00564_0033_AHCWC7ALXX/3447259045/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_14rerun/170206_ARTLoD_B1_14rerun_S4_L002_R2_001.fastq.gz'),
        FASTQPair(r1='s3://grail-clinical-fastq/MenloPark/170214_E00564_0033_AHCWC7ALXX/3447259045/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_14rerun/170206_ARTLoD_B1_14rerun_S4_L003_R1_001.fastq.gz',
                  r2='s3://grail-clinical-fastq/MenloPark/170214_E00564_0033_AHCWC7ALXX/3447259045/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_14rerun/170206_ARTLoD_B1_14rerun_S4_L003_R2_001.fastq.gz'),
        FASTQPair(r1='s3://grail-clinical-fastq/MenloPark/170214_E00564_0033_AHCWC7ALXX/3447259045/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_14rerun/170206_ARTLoD_B1_14rerun_S4_L004_R1_001.fastq.gz',
                  r2='s3://grail-clinical-fastq/MenloPark/170214_E00564_0033_AHCWC7ALXX/3447259045/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_14rerun/170206_ARTLoD_B1_14rerun_S4_L004_R2_001.fastq.gz'),
        FASTQPair(r1='s3://grail-clinical-fastq/MenloPark/170214_E00564_0033_AHCWC7ALXX/3447259045/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_14rerun/170206_ARTLoD_B1_14rerun_S4_L005_R1_001.fastq.gz',
                  r2='s3://grail-clinical-fastq/MenloPark/170214_E00564_0033_AHCWC7ALXX/3447259045/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_14rerun/170206_ARTLoD_B1_14rerun_S4_L005_R2_001.fastq.gz'),
        FASTQPair(r1='s3://grail-clinical-fastq/MenloPark/170214_E00564_0033_AHCWC7ALXX/3447259045/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_14rerun/170206_ARTLoD_B1_14rerun_S4_L006_R1_001.fastq.gz',
                  r2='s3://grail-clinical-fastq/MenloPark/170214_E00564_0033_AHCWC7ALXX/3447259045/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_14rerun/170206_ARTLoD_B1_14rerun_S4_L006_R2_001.fastq.gz'),
        FASTQPair(r1='s3://grail-clinical-fastq/MenloPark/170214_E00564_0033_AHCWC7ALXX/3447259045/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_14rerun/170206_ARTLoD_B1_14rerun_S4_L007_R1_001.fastq.gz',
                  r2='s3://grail-clinical-fastq/MenloPark/170214_E00564_0033_AHCWC7ALXX/3447259045/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_14rerun/170206_ARTLoD_B1_14rerun_S4_L007_R2_001.fastq.gz'),
        FASTQPair(r1='s3://grail-clinical-fastq/MenloPark/170214_E00564_0033_AHCWC7ALXX/3447259045/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_14rerun/170206_ARTLoD_B1_14rerun_S4_L008_R1_001.fastq.gz',
                  r2='s3://grail-clinical-fastq/MenloPark/170214_E00564_0033_AHCWC7ALXX/3447259045/bcl2fastq-2.19.1.403-umi-1mismatch-noignore/170206_ARTLoD_B1_14rerun/170206_ARTLoD_B1_14rerun_S4_L008_R2_001.fastq.gz'),
    ]),
]

# Files used in the titration benchmarks (public location)
TITRATION_BENCHMARK_DIR = f'{BENCHMARK_DATA_DIR}/titration_benchmark'
TITRATION_SAMPLES = [
    TitrationSample(name=s.name,
                    paths=[FASTQPair(r1=f'{TITRATION_BENCHMARK_DIR}/{os.path.basename(fq.r1)}',
                                     r2=f'{TITRATION_BENCHMARK_DIR}/{os.path.basename(fq.r2)}') for fq in s.paths])
    for s in ORG_TITRATION_SAMPLES]

# Describes a RNA sample used in the paper
RNASample = NamedTuple('RNASample',
                       [('name', str), # sample name
                        ('paths', List[FASTQPair])]) # FASTQ files

# Files used in the RNA benchmarks (grail-internal location)
ORG_RNA_SAMPLES = [
    RNASample(name='101CPREL277', paths=[
        FASTQPair(r1='s3://grail-cfrna-fastq/MissionBay/161115_E00481_0058_AH53VWALXX/1151367959/bcl2fastq-2.19.0.316-umi-1mismatch-noignore/101CPREL277/101CPREL277_S1_L001_R1_001.fastq.gz',
                  r2='s3://grail-cfrna-fastq/MissionBay/161115_E00481_0058_AH53VWALXX/1151367959/bcl2fastq-2.19.0.316-umi-1mismatch-noignore/101CPREL277/101CPREL277_S1_L001_R2_001.fastq.gz'),
        FASTQPair(r1='s3://grail-cfrna-fastq/MissionBay/161115_E00481_0058_AH53VWALXX/1151367959/bcl2fastq-2.19.0.316-umi-1mismatch-noignore/101CPREL277/101CPREL277_S1_L002_R1_001.fastq.gz',
                  r2='s3://grail-cfrna-fastq/MissionBay/161115_E00481_0058_AH53VWALXX/1151367959/bcl2fastq-2.19.0.316-umi-1mismatch-noignore/101CPREL277/101CPREL277_S1_L002_R2_001.fastq.gz'),
        FASTQPair(r1='s3://grail-cfrna-fastq/MissionBay/161115_E00481_0058_AH53VWALXX/1151367959/bcl2fastq-2.19.0.316-umi-1mismatch-noignore/101CPREL277/101CPREL277_S1_L003_R1_001.fastq.gz',
                  r2='s3://grail-cfrna-fastq/MissionBay/161115_E00481_0058_AH53VWALXX/1151367959/bcl2fastq-2.19.0.316-umi-1mismatch-noignore/101CPREL277/101CPREL277_S1_L003_R2_001.fastq.gz'),
        FASTQPair(r1='s3://grail-cfrna-fastq/MissionBay/161115_E00481_0058_AH53VWALXX/1151367959/bcl2fastq-2.19.0.316-umi-1mismatch-noignore/101CPREL277/101CPREL277_S1_L004_R1_001.fastq.gz',
                  r2='s3://grail-cfrna-fastq/MissionBay/161115_E00481_0058_AH53VWALXX/1151367959/bcl2fastq-2.19.0.316-umi-1mismatch-noignore/101CPREL277/101CPREL277_S1_L004_R2_001.fastq.gz'),
        FASTQPair(r1='s3://grail-cfrna-fastq/MissionBay/161115_E00481_0058_AH53VWALXX/1151367959/bcl2fastq-2.19.0.316-umi-1mismatch-noignore/101CPREL277/101CPREL277_S1_L005_R1_001.fastq.gz',
                  r2='s3://grail-cfrna-fastq/MissionBay/161115_E00481_0058_AH53VWALXX/1151367959/bcl2fastq-2.19.0.316-umi-1mismatch-noignore/101CPREL277/101CPREL277_S1_L005_R2_001.fastq.gz'),
        FASTQPair(r1='s3://grail-cfrna-fastq/MissionBay/161115_E00481_0058_AH53VWALXX/1151367959/bcl2fastq-2.19.0.316-umi-1mismatch-noignore/101CPREL277/101CPREL277_S1_L006_R1_001.fastq.gz',
                  r2='s3://grail-cfrna-fastq/MissionBay/161115_E00481_0058_AH53VWALXX/1151367959/bcl2fastq-2.19.0.316-umi-1mismatch-noignore/101CPREL277/101CPREL277_S1_L006_R2_001.fastq.gz'),
        FASTQPair(r1='s3://grail-cfrna-fastq/MissionBay/161115_E00481_0058_AH53VWALXX/1151367959/bcl2fastq-2.19.0.316-umi-1mismatch-noignore/101CPREL277/101CPREL277_S1_L007_R1_001.fastq.gz',
                  r2='s3://grail-cfrna-fastq/MissionBay/161115_E00481_0058_AH53VWALXX/1151367959/bcl2fastq-2.19.0.316-umi-1mismatch-noignore/101CPREL277/101CPREL277_S1_L007_R2_001.fastq.gz'),
        FASTQPair(r1='s3://grail-cfrna-fastq/MissionBay/161115_E00481_0058_AH53VWALXX/1151367959/bcl2fastq-2.19.0.316-umi-1mismatch-noignore/101CPREL277/101CPREL277_S1_L008_R1_001.fastq.gz',
                  r2='s3://grail-cfrna-fastq/MissionBay/161115_E00481_0058_AH53VWALXX/1151367959/bcl2fastq-2.19.0.316-umi-1mismatch-noignore/101CPREL277/101CPREL277_S1_L008_R2_001.fastq.gz')]),
    # Note: the following filepair is the original. Files in s3://grail-ysaito/af4_benchmark are created by running split_fastq.py on the original fastq pair.
    #
    #FASTQPair(r1='s3://grail-cfrna-fastq/MissionBay/161206_E00501_0049_AHCYJYALXX/3747505862/bcl2fastq-2.19.0.316-umi-1mismatch-noignore/108CPREL315/108CPREL315_S4_L004_R1_001.fastq.gz',
    #         r2='s3://grail-cfrna-fastq/MissionBay/161206_E00501_0049_AHCYJYALXX/3747505862/bcl2fastq-2.19.0.316-umi-1mismatch-noignore/108CPREL315/108CPREL315_S4_L004_R2_001.fastq.gz')],
    RNASample(name='108CPREL315', paths=[
        FASTQPair(r1='s3://grail-ysaito/af4_benchmark/108CPREL315_S4_L004_R1_001-00.fastq.gz',
                  r2='s3://grail-ysaito/af4_benchmark/108CPREL315_S4_L004_R2_001-00.fastq.gz'),
        FASTQPair(r1='s3://grail-ysaito/af4_benchmark/108CPREL315_S4_L004_R1_001-01.fastq.gz',
                  r2='s3://grail-ysaito/af4_benchmark/108CPREL315_S4_L004_R2_001-01.fastq.gz'),
        FASTQPair(r1='s3://grail-ysaito/af4_benchmark/108CPREL315_S4_L004_R1_001-02.fastq.gz',
                  r2='s3://grail-ysaito/af4_benchmark/108CPREL315_S4_L004_R2_001-02.fastq.gz'),
        FASTQPair(r1='s3://grail-ysaito/af4_benchmark/108CPREL315_S4_L004_R1_001-03.fastq.gz',
                  r2='s3://grail-ysaito/af4_benchmark/108CPREL315_S4_L004_R2_001-03.fastq.gz'),
        FASTQPair(r1='s3://grail-ysaito/af4_benchmark/108CPREL315_S4_L004_R1_001-04.fastq.gz',
                  r2='s3://grail-ysaito/af4_benchmark/108CPREL315_S4_L004_R2_001-04.fastq.gz'),
        FASTQPair(r1='s3://grail-ysaito/af4_benchmark/108CPREL315_S4_L004_R1_001-05.fastq.gz',
                  r2='s3://grail-ysaito/af4_benchmark/108CPREL315_S4_L004_R2_001-05.fastq.gz'),
        FASTQPair(r1='s3://grail-ysaito/af4_benchmark/108CPREL315_S4_L004_R1_001-06.fastq.gz',
                  r2='s3://grail-ysaito/af4_benchmark/108CPREL315_S4_L004_R2_001-06.fastq.gz'),
        FASTQPair(r1='s3://grail-ysaito/af4_benchmark/108CPREL315_S4_L004_R1_001-07.fastq.gz',
                  r2='s3://grail-ysaito/af4_benchmark/108CPREL315_S4_L004_R2_001-07.fastq.gz'),
        FASTQPair(r1='s3://grail-ysaito/af4_benchmark/108CPREL315_S4_L004_R1_001-08.fastq.gz',
                  r2='s3://grail-ysaito/af4_benchmark/108CPREL315_S4_L004_R2_001-08.fastq.gz')]),

    #FASTQPair(r1='s3://grail-cfrna-fastq/MissionBay/170105_K00215_0113_AHF2N2BBXX/3062013994/bcl2fastq-2.19.0.316-umi-1mismatch-noignore/74HPREL332/74HPREL332_S1_L001_R1_001.fastq.gz',
    #         r2='s3://grail-cfrna-fastq/MissionBay/170105_K00215_0113_AHF2N2BBXX/3062013994/bcl2fastq-2.19.0.316-umi-1mismatch-noignore/74HPREL332/74HPREL332_S1_L001_R2_001.fastq.gz')],
    RNASample(name='74HPREL332', paths=[
        FASTQPair(r1='s3://grail-ysaito/af4_benchmark/74HPREL332_S1_L001_R1_001-00.fastq.gz',
                  r2='s3://grail-ysaito/af4_benchmark/74HPREL332_S1_L001_R2_001-00.fastq.gz'),
        FASTQPair(r1='s3://grail-ysaito/af4_benchmark/74HPREL332_S1_L001_R1_001-01.fastq.gz',
                  r2='s3://grail-ysaito/af4_benchmark/74HPREL332_S1_L001_R2_001-01.fastq.gz'),
        FASTQPair(r1='s3://grail-ysaito/af4_benchmark/74HPREL332_S1_L001_R1_001-02.fastq.gz',
                  r2='s3://grail-ysaito/af4_benchmark/74HPREL332_S1_L001_R2_001-02.fastq.gz'),
        FASTQPair(r1='s3://grail-ysaito/af4_benchmark/74HPREL332_S1_L001_R1_001-03.fastq.gz',
                  r2='s3://grail-ysaito/af4_benchmark/74HPREL332_S1_L001_R2_001-03.fastq.gz')]),
    RNASample(name='118HPREL322', paths=[
        #FASTQPair(r1='s3://grail-cfrna-fastq/MissionBay/161206_E00501_0050_BHF2JCALXX/590215334/bcl2fastq-2.19.0.316-umi-1mismatch-noignore/118HPREL322/118HPREL322_S3_L003_R1_001.fastq.gz',
        #        r2='s3://grail-cfrna-fastq/MissionBay/161206_E00501_0050_BHF2JCALXX/590215334/bcl2fastq-2.19.0.316-umi-1mismatch-noignore/118HPREL322/118HPREL322_S3_L003_R2_001.fastq.gz')],
        FASTQPair(r1='s3://grail-ysaito/af4_benchmark/118HPREL322_S3_L003_R1_001-00.fastq.gz',
                  r2='s3://grail-ysaito/af4_benchmark/118HPREL322_S3_L003_R2_001-00.fastq.gz'),
        FASTQPair(r1='s3://grail-ysaito/af4_benchmark/118HPREL322_S3_L003_R1_001-01.fastq.gz',
                  r2='s3://grail-ysaito/af4_benchmark/118HPREL322_S3_L003_R2_001-01.fastq.gz'),
        FASTQPair(r1='s3://grail-ysaito/af4_benchmark/118HPREL322_S3_L003_R1_001-02.fastq.gz',
                  r2='s3://grail-ysaito/af4_benchmark/118HPREL322_S3_L003_R2_001-02.fastq.gz'),
        FASTQPair(r1='s3://grail-ysaito/af4_benchmark/118HPREL322_S3_L003_R1_001-03.fastq.gz',
                  r2='s3://grail-ysaito/af4_benchmark/118HPREL322_S3_L003_R2_001-03.fastq.gz'),
        FASTQPair(r1='s3://grail-ysaito/af4_benchmark/118HPREL322_S3_L003_R1_001-04.fastq.gz',
                  r2='s3://grail-ysaito/af4_benchmark/118HPREL322_S3_L003_R2_001-04.fastq.gz'),
        FASTQPair(r1='s3://grail-ysaito/af4_benchmark/118HPREL322_S3_L003_R1_001-05.fastq.gz',
                  r2='s3://grail-ysaito/af4_benchmark/118HPREL322_S3_L003_R2_001-05.fastq.gz'),
        FASTQPair(r1='s3://grail-ysaito/af4_benchmark/118HPREL322_S3_L003_R1_001-06.fastq.gz',
                  r2='s3://grail-ysaito/af4_benchmark/118HPREL322_S3_L003_R2_001-06.fastq.gz'),
        FASTQPair(r1='s3://grail-ysaito/af4_benchmark/118HPREL322_S3_L003_R1_001-07.fastq.gz',
                  r2='s3://grail-ysaito/af4_benchmark/118HPREL322_S3_L003_R2_001-07.fastq.gz')]),
    # RNASample(name='68HPREL273', paths=[
    #     FASTQPair(r1='s3://grail-cfrna-fastq/MissionBay/161115_E00481_0058_AH53VWALXX/1151367959/bcl2fastq-2.19.0.316-umi-1mismatch-noignore/68HPREL273/68HPREL273_S5_L001_R1_001.fastq.gz',
    #              r2='s3://grail-cfrna-fastq/MissionBay/161115_E00481_0058_AH53VWALXX/1151367959/bcl2fastq-2.19.0.316-umi-1mismatch-noignore/68HPREL273/68HPREL273_S5_L001_R2_001.fastq.gz'),
    #     FASTQPair(r1='s3://grail-cfrna-fastq/MissionBay/161115_E00481_0058_AH53VWALXX/1151367959/bcl2fastq-2.19.0.316-umi-1mismatch-noignore/68HPREL273/68HPREL273_S5_L002_R1_001.fastq.gz',
    #              r2='s3://grail-cfrna-fastq/MissionBay/161115_E00481_0058_AH53VWALXX/1151367959/bcl2fastq-2.19.0.316-umi-1mismatch-noignore/68HPREL273/68HPREL273_S5_L002_R2_001.fastq.gz'),
    #     FASTQPair(r1='s3://grail-cfrna-fastq/MissionBay/161115_E00481_0058_AH53VWALXX/1151367959/bcl2fastq-2.19.0.316-umi-1mismatch-noignore/68HPREL273/68HPREL273_S5_L003_R1_001.fastq.gz',
    #              r2='s3://grail-cfrna-fastq/MissionBay/161115_E00481_0058_AH53VWALXX/1151367959/bcl2fastq-2.19.0.316-umi-1mismatch-noignore/68HPREL273/68HPREL273_S5_L003_R2_001.fastq.gz'),
    #     FASTQPair(r1='s3://grail-cfrna-fastq/MissionBay/161115_E00481_0058_AH53VWALXX/1151367959/bcl2fastq-2.19.0.316-umi-1mismatch-noignore/68HPREL273/68HPREL273_S5_L004_R1_001.fastq.gz',
    #              r2='s3://grail-cfrna-fastq/MissionBay/161115_E00481_0058_AH53VWALXX/1151367959/bcl2fastq-2.19.0.316-umi-1mismatch-noignore/68HPREL273/68HPREL273_S5_L004_R2_001.fastq.gz'),
    #     FASTQPair(r1='s3://grail-cfrna-fastq/MissionBay/161115_E00481_0058_AH53VWALXX/1151367959/bcl2fastq-2.19.0.316-umi-1mismatch-noignore/68HPREL273/68HPREL273_S5_L005_R1_001.fastq.gz',
    #              r2='s3://grail-cfrna-fastq/MissionBay/161115_E00481_0058_AH53VWALXX/1151367959/bcl2fastq-2.19.0.316-umi-1mismatch-noignore/68HPREL273/68HPREL273_S5_L005_R2_001.fastq.gz'),
    #     FASTQPair(r1='s3://grail-cfrna-fastq/MissionBay/161115_E00481_0058_AH53VWALXX/1151367959/bcl2fastq-2.19.0.316-umi-1mismatch-noignore/68HPREL273/68HPREL273_S5_L006_R1_001.fastq.gz',
    #              r2='s3://grail-cfrna-fastq/MissionBay/161115_E00481_0058_AH53VWALXX/1151367959/bcl2fastq-2.19.0.316-umi-1mismatch-noignore/68HPREL273/68HPREL273_S5_L006_R2_001.fastq.gz'),
    #     FASTQPair(r1='s3://grail-cfrna-fastq/MissionBay/161115_E00481_0058_AH53VWALXX/1151367959/bcl2fastq-2.19.0.316-umi-1mismatch-noignore/68HPREL273/68HPREL273_S5_L007_R1_001.fastq.gz',
    #              r2='s3://grail-cfrna-fastq/MissionBay/161115_E00481_0058_AH53VWALXX/1151367959/bcl2fastq-2.19.0.316-umi-1mismatch-noignore/68HPREL273/68HPREL273_S5_L007_R2_001.fastq.gz'),
    #     FASTQPair(r1='s3://grail-cfrna-fastq/MissionBay/161115_E00481_0058_AH53VWALXX/1151367959/bcl2fastq-2.19.0.316-umi-1mismatch-noignore/68HPREL273/68HPREL273_S5_L008_R1_001.fastq.gz',
    #              r2='s3://grail-cfrna-fastq/MissionBay/161115_E00481_0058_AH53VWALXX/1151367959/bcl2fastq-2.19.0.316-umi-1mismatch-noignore/68HPREL273/68HPREL273_S5_L008_R2_001.fastq.gz')]),

    #FASTQPair(r1='s3://grail-cfrna-fastq/MissionBay/160927_K00122_0179_BHFWLKBBXX/3600619798/bcl2fastq-2.19.0.316-umi-1mismatch-noignore/53HPREL160/53HPREL160_S3_L003_R1_001.fastq.gz',
    #        r2='s3://grail-cfrna-fastq/MissionBay/160927_K00122_0179_BHFWLKBBXX/3600619798/bcl2fastq-2.19.0.316-umi-1mismatch-noignore/53HPREL160/53HPREL160_S3_L003_R2_001.fastq.gz')],
    RNASample(name='53HPREL160', paths=[
        FASTQPair(r1='s3://grail-ysaito/af4_benchmark/53HPREL160_S3_L003_R1_001-00.fastq.gz',
                  r2='s3://grail-ysaito/af4_benchmark/53HPREL160_S3_L003_R2_001-00.fastq.gz'),
        FASTQPair(r1='s3://grail-ysaito/af4_benchmark/53HPREL160_S3_L003_R1_001-01.fastq.gz',
                  r2='s3://grail-ysaito/af4_benchmark/53HPREL160_S3_L003_R2_001-01.fastq.gz')]),
    # RNASample(name='117HPREL321', paths=[
    #     #FASTQPair(r1='s3://grail-cfrna-fastq/MissionBay/161206_E00501_0050_BHF2JCALXX/590215334/bcl2fastq-2.19.0.316-umi-1mismatch-noignore/117HPREL321/117HPREL321_S2_L002_R1_001.fastq.gz',
    #     #r2='s3://grail-cfrna-fastq/MissionBay/161206_E00501_0050_BHF2JCALXX/590215334/bcl2fastq-2.19.0.316-umi-1mismatch-noignore/117HPREL321/117HPREL321_S2_L002_R2_001.fastq.gz')])

    #     FASTQPair(r1='s3://grail-ysaito/af4_benchmark/117HPREL321_S2_L002_R1_001-00.fastq.gz',
    #             r2='s3://grail-ysaito/af4_benchmark/117HPREL321_S2_L002_R2_001-00.fastq.gz'),
    #     FASTQPair(r1='s3://grail-ysaito/af4_benchmark/117HPREL321_S2_L002_R1_001-01.fastq.gz',
    #              r2='s3://grail-ysaito/af4_benchmark/117HPREL321_S2_L002_R2_001-01.fastq.gz'),
    #     FASTQPair(r1='s3://grail-ysaito/af4_benchmark/117HPREL321_S2_L002_R1_001-02.fastq.gz',
    #              r2='s3://grail-ysaito/af4_benchmark/117HPREL321_S2_L002_R2_001-02.fastq.gz'),
    #     FASTQPair(r1='s3://grail-ysaito/af4_benchmark/117HPREL321_S2_L002_R1_001-03.fastq.gz',
    #              r2='s3://grail-ysaito/af4_benchmark/117HPREL321_S2_L002_R2_001-03.fastq.gz'),
    #     FASTQPair(r1='s3://grail-ysaito/af4_benchmark/117HPREL321_S2_L002_R1_001-04.fastq.gz',
    #              r2='s3://grail-ysaito/af4_benchmark/117HPREL321_S2_L002_R2_001-04.fastq.gz'),
    #     FASTQPair(r1='s3://grail-ysaito/af4_benchmark/117HPREL321_S2_L002_R1_001-05.fastq.gz',
    #              r2='s3://grail-ysaito/af4_benchmark/117HPREL321_S2_L002_R2_001-05.fastq.gz'),
    #     FASTQPair(r1='s3://grail-ysaito/af4_benchmark/117HPREL321_S2_L002_R1_001-06.fastq.gz',
    #              r2='s3://grail-ysaito/af4_benchmark/117HPREL321_S2_L002_R2_001-06.fastq.gz'),
    #     FASTQPair(r1='s3://grail-ysaito/af4_benchmark/117HPREL321_S2_L002_R1_001-07.fastq.gz',
    #              r2='s3://grail-ysaito/af4_benchmark/117HPREL321_S2_L002_R2_001-07.fastq.gz'),
    #     FASTQPair(r1='s3://grail-ysaito/af4_benchmark/117HPREL321_S2_L002_R1_001-08.fastq.gz',
    #              r2='s3://grail-ysaito/af4_benchmark/117HPREL321_S2_L002_R2_001-08.fastq.gz')],
]

RNA_BENCHMARK_DIR = f'{BENCHMARK_DATA_DIR}/rna_benchmark'
RNA_SAMPLES = [
    RNASample(name=s.name,
              paths=[FASTQPair(r1=f'{RNA_BENCHMARK_DIR}/{os.path.basename(fq.r1)}',
                               r2=f'{RNA_BENCHMARK_DIR}/{os.path.basename(fq.r2)}') for fq in s.paths])
    for s in ORG_RNA_SAMPLES]

BazelLabel = NamedTuple('BazelLabel',
                        [('label', str),
                         ('absolute', bool),
                         ('package', str),
                         ('target', str)])

def parse_bazel_label(label: str) -> BazelLabel:
    m = re.match('([^:]+):(.*)$', label)
    if m:
        package, target = m[1], m[2]
    else:
        m = re.match('([^:]+)/([^/]+)$', label)
        if not m:
            raise Exception(f'Failed to parse label {label}')
        package = m[0]
        target = m[2]
        label = f'{package}:{target}'

    return BazelLabel(label=label,
                      absolute=package.startswith('//'),
                      package=package,
                      target=target)

def check_call(args: List[str]) -> None:
    logging.info('Run: %s', ' '.join(args))
    subprocess.check_call(args)

def check_output(args: List[str]) -> str:
    logging.info('Run: %s', ' '.join(args))
    return subprocess.check_output(args, universal_newlines=True)

def repo_root() -> Path:
    """Get the root directory of the repository."""
    commit = check_output(['git', 'rev-parse', '--show-toplevel'])
    return Path(commit.strip())

def build(labels: List[str]) -> None:
    check_call(['bazel', 'build'] + labels)

def go_executable(label: str) -> Path:
    build([label])
    p = parse_bazel_label(label)
    if not p.absolute:
        raise Exception(f"Target {p} must be absolute")
    return repo_root() / 'bazel-bin' / Path(p.package[2:]) / 'linux_amd64_stripped' / Path(p.target)

def nongo_executable(label: str) -> Path:
    p = parse_bazel_label(label)
    if not p.absolute:
        raise Exception(f"Target {p} must be absolute")
    return repo_root() / 'bazel-bin' / Path(p.package[2:]) / Path(p.target)

GRAIL_FILE_PATH: Optional[Path] = None

def grail_file_path() -> Path:
    """Return the abspath of grail-file binary. Builds the binary if necessary"""
    global GRAIL_FILE_PATH
    if not GRAIL_FILE_PATH:
        target = '//go/src/github.com/grailbio/base/cmd/grail-file'
        build([target])
        GRAIL_FILE_PATH = go_executable(target)
    return GRAIL_FILE_PATH

AF4_PATH: Optional[Path] = None

def af4_path() -> Path:
    """Return the abspath of Go bio-target-rna-fusion binary. Builds the binary if necessary"""
    global AF4_PATH
    if not AF4_PATH:
        af4_label = '//go/src/github.com/grailbio/bio/cmd/bio-fusion'
        build([af4_label])
        AF4_PATH = go_executable(af4_label)
    return AF4_PATH

def expand_fastq_files(fq: Iterable[FASTQPair]) -> List[str]:
    return [x.r1 for x in fq] + [x.r2 for x in fq]

def s3_cache_files(src_paths: List[str], cache_dir: Path, force=False) -> None:
    """Copy src_paths in cache_dir if they haven't been copied already."""
    args = [str(grail_file_path()), 'cp', '-v']
    n = 0
    for src_path in src_paths:
        dest_path = str(cache_dir) + '/' + os.path.basename(src_path)
        if not os.path.exists(dest_path) or force:
            args.append(src_path)
            n += 1
    if n == 0:
        return
    args.append(str(cache_dir) + '/')
    check_call(args)

def s3_cache_dir(src_dir: str, cache_dir: Path) -> None:
    try:
        check_call([str(grail_file_path()), 'cp', '-R', '-v', src_dir, str(cache_dir)])
    except Exception as e:
        logging.error('Error(ignored): %s', e)

RunStats = NamedTuple('RunStats', [
    ('duration', float),  # duration of the run, in seconds
    ('n_fragments', int),  # total # of fragments
    ('all_candidates', int), # # candidates found in the 1st stage.
    ('low_complexity_substring', int), # # candidates found to have low-complexity substring.
    ('close_proximity', int), # # candidates found to have genepairs in close proximity.
    ('duplicates', int), # # candidates found to be UMI duplicates
    ('min_span', int), # # candidates filtered by minspan.
    ('abundant_partners', int), # # candidates to have genes w/ abundant partners.
    ('final_candidates', int), # # candidates found in the 1st stage.
    #('fusion_stats', str), # %+v dump of fusion.Stats
    ('n_genes', int),
    ('n_fragment_matches', List[int]),
])

def parse_time(m) -> float:
    """Parse the info log timestamp."""
    return int(m[1]) * 3600 + int(m[2]) * 60  + int(m[3]) + (int(m[4]) / 1000000.0)

def run_stats(dir_path: Path) -> RunStats:
    """Parse the INFO log file produced by bio-target-rna-fusion and extract high-level stats."""

    info_paths = glob.glob(str(dir_path / '*.INFO'))
    if len(info_paths) != 1:
        raise Exception(f'{dir_path}: No INFO file found ({info_paths})')
    ts_re = r'[IEW\d]+ (\d\d+):(\d\d):(\d\d).(\d\d\d\d\d\d)\s+\d+.*'
    start_time = 0.0
    end_time = 0.0
    n_fragments = 0

    all_candidates = -9999999
    low_complexity_substring = -9999999
    close_proximity = -9999999
    duplicates = -9999999
    n_remaining_after_close_proximity = -9999999
    n_remaining_after_duplicates = -9999999
    n_remaining_after_min_span = -9999999
    final_candidates = -9999999
    n_fragments2 = -999999
    n_genes = -999999
    n_fragment_matches: List[int] = []
    with open(info_paths[0]) as fd:
        for line in fd.readlines():
            m = re.match(ts_re + r'Start reading geneDB', line)
            if m:
                start_time = parse_time(m)
            m = re.match(ts_re + r'All done', line)
            if m:
                end_time = parse_time(m)
            m = re.match(r'.*Processed (\d+) reads in', line)
            if m:
                n_fragments += int(m[1])
            m = re.match(r'.*Starting filtering (\d+) candidates', line)
            if m:
                all_candidates = int(m[1])
            m = re.match(r'.*Stats: (\d+) candidates after stage 1', line)
            if m:
                all_candidates = int(m[1])
            m = re.match(r'.* (\d+) of \d+ remaining after removing (\d+) low-complex substring and (\d+) close proximity', line)
            if m:
                n_remaining_after_close_proximity = int(m[1])
                low_complexity_substring = int(m[2])
                close_proximity = int(m[3])
            m = re.match(r'.* (\d+) remaining after removing duplicates', line)
            if m:
                n_remaining_after_duplicates = int(m[1])
            m = re.match(r'.* (\d+) remaining after filtering by minspan', line)
            if m:
                n_remaining_after_min_span = int(m[1])
            m = re.match(r'.*Wrote (\d+) filtered candidates', line)  # old log format
            if m:
                final_candidates = int(m[1])
            m = re.match(r'.*Stats: (\d+) final candidates', line)
            if m:
                final_candidates = int(m[1])
            m = re.match(r'.*Stats:.*\{LowComplexity.* Genes:(\d+) Fragments:(\d+) FragmentsWithMatchingGenes:\[(\d+) (\d+) (\d+) (\d+) (\d+)\]', line)
            if m:
                n_genes = int(m[1])
                n_fragments2 = int(m[2])
                n_fragment_matches = [int(m[3]), int(m[4]), int(m[5]), int(m[6]), int(m[7])]

    return RunStats(duration=(end_time - start_time),
                    n_fragments=n_fragments,
                    all_candidates=all_candidates,
                    low_complexity_substring=low_complexity_substring,
                    close_proximity=close_proximity,
                    duplicates=(n_remaining_after_close_proximity - n_remaining_after_duplicates),
                    min_span=(n_remaining_after_duplicates - n_remaining_after_min_span),
                    abundant_partners=(n_remaining_after_min_span - final_candidates),
                    final_candidates=final_candidates,
                    n_genes=n_genes,
                    n_fragment_matches=n_fragment_matches)

def pretty_sample_name(path: str) -> str:
    """Given a path of a benchmark result dir,
    return the prettified sample name that's used in the paper.

    It returns the path basename by default.
    """

    name_map = {
        "170206_ARTLoD_B1_01rerun" : "T1 (0.0001)",
        "170206_ARTLoD_B1_02rerun" : "T2 (0.0002)",
        "170206_ARTLoD_B1_03rerun" : "T3 (0.0002)",
        "170206_ARTLoD_B1_04rerun" : "T4 (0.0004)",
        "170206_ARTLoD_B1_05rerun" : "T5 (0.0004)",
        "170206_ARTLoD_B1_06rerun" : "T6 (0.0004)",
        "170206_ARTLoD_B1_07rerun" : "T7 (0.0006)",
        "170206_ARTLoD_B1_08rerun" : "T8 (0.0006)",
        "170206_ARTLoD_B1_09rerun" : "T9 (0.0006)",
        "170206_ARTLoD_B1_10rerun" : "T10 (0.0008)",
        "170206_ARTLoD_B1_11rerun" : "T11 (0.0008)",
        "170206_ARTLoD_B1_12rerun" : "T12 (0.0008)",
        "170206_ARTLoD_B1_13rerun" : "T13 (0.01)",
        "170206_ARTLoD_B1_14rerun" : "T14 (0.01)",

        # RNA datasets
        "101CPREL277" : "Prc101",
        "108CPREL315" : "Prc108",
        "74HPREL332" : "HC332",
        "118HPREL322" : "HC118",
        "68HPREL273" : "HC273",  # not used
        "53HPREL160" : "HC160",
        "117HPREL321" : "HC117", # not used
    }
    if path[-1] == '/':
        name = os.path.basename(path[:-1])
    else:
        name = os.path.basename(path)
    for key, val in name_map.items():
        name = name.replace(key, val)
    return name
