#!/usr/bin/env python3

"""
A simple script to run the C++ and Go versions of AF4 in a controlled setting and compare results.

Example:
DIR=/scratch-nvme/fusion
run.py --cosmic_fusion=/scratch-nvme/fusion/all_pair.txt --r1=$DIR/32-1_S1_L001_R1_001.fastq.gz,$DIR/32-1_S1_L002_R1_001.fastq.gz,$DIR/32-1_S1_L003_R1_001.fastq.gz,$DIR/32-1_S1_L004_R1_001.fastq.gz,$DIR/32-1_S1_L005_R1_001.fastq.gz,$DIR/32-1_S1_L006_R1_001.fastq.gz,$DIR/32-1_S1_L007_R1_001.fastq.gz,$DIR/32-1_S1_L008_R1_001.fastq.gz --r2=$DIR/32-1_S1_L001_R2_001.fastq.gz,$DIR/32-1_S1_L002_R2_001.fastq.gz,$DIR/32-1_S1_L003_R2_001.fastq.gz,$DIR/32-1_S1_L004_R2_001.fastq.gz,$DIR/32-1_S1_L005_R2_001.fastq.gz,$DIR/32-1_S1_L006_R2_001.fastq.gz,$DIR/32-1_S1_L007_R2_001.fastq.gz,$DIR/32-1_S1_L008_R2_001.fastq.gz
"""

import argparse
import os
import logging
from pathlib import Path
import re
import subprocess
from typing import List, NamedTuple, Optional

import util

Config = NamedTuple('Config', [
    # Directory to download & cache benchmark files.
    ('cache_dir', Path),
    # If true, update the gene-list file so that the gene sort order is the same between Go and C++.
    ('update_gene_list', bool),
    # R1 and R2 fastq files
    ('r1', str),
    ('r2', str),
    # Transcriptome files
    ('transcript', Path),
    # Fusion event pair files. If none, run in denovo mode.
    ('cosmic_fusion', Optional[Path]), # None when denovo
    # 'c++' or 'go' or both.
    ('run', List[str])])

GO_WDIR = Path('/tmp/fusion-go')

def sort_file(in_path: Path, out_path: Optional[Path]) -> None:
    if out_path is None:
        out_path = Path(str(in_path) + '.sorted')
    lines = in_path.open().readlines()
    lines.sort()
    with out_path.open('w') as out:
        for line in lines:
            print(line, file=out, end='')

def sort_fasta_headers(in_path: Path, out_path: Path) -> None:
    lines: List[str] = []
    for line in in_path.open().readlines():
        if line.startswith('>'):
            lines.append(line)
    lines.sort()
    with out_path.open('w') as out:
        for line in lines:
            print(line, file=out, end='')

def minimal_go_args(config: Config) -> List[str]:
    args = [f'-r1={config.r1}',
            f'-r2={config.r2}',
            f'-transcript={config.transcript}',
            '--pprof=:12345']
    if config.cosmic_fusion:
        args.append(f'-cosmic-fusion={config.cosmic_fusion}')
    return args
#'-k=19',
#'-max-genes-per-kmer=2',
#'-max-proximity-distance=1000',
#'-max-proximity-genes=0',
#

def run_go(args: List[str]) -> None:
    bin_label = '//go/src/grail.com/cmd/bio-target-rna-fusion'
    util.build([bin_label])
    bin_path = util.go_executable(bin_label)
    logging.info("Start: go: %s %s", bin_path, ' '.join(args))
    subprocess.check_call([str(bin_path)] + args)


def gene_list_path(config: Config) -> Path:
    if not config.cosmic_fusion:
        return Path('/scratch-nvme/fusion/gene-names-denovo.txt')
    return Path('/scratch-nvme/fusion/gene-names.txt')

CPP_WDIR = Path('/tmp/fusion-cpp')
def default_cpp_args(config: Config) -> List[str]:
    return [f'--r1_fqgz={config.r1}',
            f'--r2_fqgz={config.r2}',
            f'--transcript_file={config.transcript}',
            f'--wdir={str(CPP_WDIR)}',
            '--p=64',
            '--stderrthreshold=0',
            '--v=1',
            '--k=19',
            '--cap=2',
            '--proximity_dist=1000',
            '--proximity_num=0',
            '--unstranded_library=1']

def run_cpp(args: List[str]) -> None:
    CPP_WDIR.mkdir(parents=True, exist_ok=True)
    bin_label = '//bio/rna/fusion:target_rna_fusion'
    util.build([bin_label])
    bin_path = util.nongo_executable(bin_label)
    logging.info("Start: c++: %s %s", bin_path, ' '.join(args))
    subprocess.check_call([str(bin_path)] + args)

def main() -> None:
    logging.basicConfig(
        level=logging.DEBUG,
        format="%(asctime)s:%(levelname)s: %(message)s")

    p = argparse.ArgumentParser()
    p.add_argument('--cache_dir', default='/scratch-nvme/cache_tmp',
                   help='Benchmark cache dir')
    p.add_argument('--cosmic_fusion', default='/scratch-nvme/fusion/all_pair.txt',
                   help='Fusion event file. If empty, af4 runs in the denovo mode')
    p.add_argument('--run', action='append', choices=['go', 'c++'])
    p.add_argument('--r1', default='/scratch-nvme/fusion/32-1_S1_L001_R1_001.fastq.gz')
    p.add_argument('--r2', default='/scratch-nvme/fusion/32-1_S1_L001_R2_001.fastq.gz')
    p.add_argument('--update_gene_list', default=True, action='store_true')
    args = p.parse_args()
    if not args.run:
        args.run = ['go', 'c++']
    config = Config(
        cache_dir=Path(args.cache_dir),
        run=args.run,
        cosmic_fusion=Path(args.cosmic_fusion) if args.cosmic_fusion else None,
        transcript=Path('/scratch-nvme/fusion/transcriptome.fa'),
        r1=args.r1,
        r2=args.r2,
        update_gene_list=args.update_gene_list)

    output_suffix = '-' + os.path.basename(os.path.splitext(config.r1)[0])
    if not config.cosmic_fusion:
        output_suffix += '-denovo'
    else:
        output_suffix += '-' + os.path.basename(os.path.splitext(config.cosmic_fusion)[0])

    if 'go' in config.run:
        gene_list = gene_list_path(config)
        if config.update_gene_list:
            run_go(minimal_go_args(config) + [f'-gene-list-output={gene_list}'])
            sort_file(gene_list, gene_list)
        run_go(minimal_go_args(config) +
               [f'-gene-list={gene_list}',
                f'-fasta-output={GO_WDIR}/all{output_suffix}.fa',
                f'-rio-output={GO_WDIR}/all{output_suffix}.rio',
                f'-filtered-output={GO_WDIR}/filtered{output_suffix}.fa'])
        sort_fasta_headers(GO_WDIR / f'all{output_suffix}.fa', GO_WDIR / f'all{output_suffix}.fa.sorted')
        sort_fasta_headers(GO_WDIR / f'filtered{output_suffix}.fa', GO_WDIR / f'filtered{output_suffix}.fa.sorted')

    if 'c++' in config.run:
        cpp_args = default_cpp_args(config)
        if not config.cosmic_fusion:
            cpp_args.append('--denovo')
        else:
            cpp_args.append(f'--cosmic_fusion_file={config.cosmic_fusion}')
        run_cpp(cpp_args)
        all_fa = list(CPP_WDIR.glob('fusion_[0-9]*.fa'))[0]
        sort_fasta_headers(all_fa, CPP_WDIR / f'all{output_suffix}.fa.sorted')
        filtered_fa = list(CPP_WDIR.glob('*final*[0-9].fa'))[0]
        sort_fasta_headers(filtered_fa, CPP_WDIR / f'filtered{output_suffix}.fa.sorted')
    logging.info("All done")

main()
