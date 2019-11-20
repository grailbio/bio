#!/usr/bin/env python3.6

"""
A simple script to run the C++ and Go versions of AF4 in a controlled setting and compare results.

Example:
DIR=/scratch-nvme/fusion
run.py --cosmic_fusion=/scratch-nvme/fusion/all_pair.txt --r1=$DIR/32-1_S1_L001_R1_001.fastq.gz,$DIR/32-1_S1_L002_R1_001.fastq.gz,$DIR/32-1_S1_L003_R1_001.fastq.gz,$DIR/32-1_S1_L004_R1_001.fastq.gz,$DIR/32-1_S1_L005_R1_001.fastq.gz,$DIR/32-1_S1_L006_R1_001.fastq.gz,$DIR/32-1_S1_L007_R1_001.fastq.gz,$DIR/32-1_S1_L008_R1_001.fastq.gz --r2=$DIR/32-1_S1_L001_R2_001.fastq.gz,$DIR/32-1_S1_L002_R2_001.fastq.gz,$DIR/32-1_S1_L003_R2_001.fastq.gz,$DIR/32-1_S1_L004_R2_001.fastq.gz,$DIR/32-1_S1_L005_R2_001.fastq.gz,$DIR/32-1_S1_L006_R2_001.fastq.gz,$DIR/32-1_S1_L007_R2_001.fastq.gz,$DIR/32-1_S1_L008_R2_001.fastq.gz
"""

import sys

sys.path.append("../../fusion/benchmark")

import argparse
import os
import glob
import logging
import enum
from pathlib import Path
import re
import subprocess
from typing import List, NamedTuple, Optional

import util


class RunMode(enum.Enum):
    DENOVO = 1
    TARGETED = 2
    TARGETED_WITH_READ_THROUGH = 3


Config = NamedTuple(
    "Config",
    [
        ("mode", RunMode),
        # Directory to download & cache benchmark files.
        ("cache_dir", Path),
        # If true, update the gene-list file so that the gene sort order is the same between Go and C++.
        ("update_gene_list", bool),
        # R1 and R2 fastq files
        ("r1", str),
        ("r2", str),
        # Transcriptome files
        ("transcript", Path),
        # Fusion event pair files. If none, run in denovo mode.
        ("cosmic_fusion", Optional[Path]),  # None when denovo
        # 'c++' or 'go' or both.
        ("run", List[str]),
    ],
)

GO_WDIR = Path("/tmp/fusion-go")


def sort_file(in_path: Path, out_path: Optional[Path]) -> None:
    if out_path is None:
        out_path = Path(str(in_path) + ".sorted")
    lines = in_path.open().readlines()
    lines.sort()
    with out_path.open("w") as out:
        for line in lines:
            print(line, file=out, end="")


def sort_fasta_headers(in_path: Path, out_path: Path) -> None:
    lines: List[str] = []
    for line in in_path.open().readlines():
        if line.startswith(">"):
            segments = line.split("|")
            name_segments = segments[0].split(":")
            lines.append(name_segments[0] + "|" + "|".join(segments[1:]))
    lines.sort()
    with out_path.open("w") as out:
        for line in lines:
            print(line, file=out, end="")


def go_args(config: Config, gene_list_path: Path, output_suffix: str) -> List[str]:
    args = [
        f"-r1={config.r1}",
        f"-r2={config.r2}",
        f"-transcript={config.transcript}",
        f"-gene-list={str(gene_list_path)}",
        f"-fasta-output={GO_WDIR}/all{output_suffix}.fa",
        f"-rio-output={GO_WDIR}/all{output_suffix}.rio",
        f"-filtered-output={GO_WDIR}/filtered{output_suffix}.fa",
        "--pprof=:12345",
    ]
    if config.mode == RunMode.DENOVO:
        args += ["-k=19", "-umi-in-name", "-max-genes-per-kmer=2"]
    elif config.mode == RunMode.TARGETED:
        args += ["--k=19", "--umi-in-name", f"--cosmic-fusion={config.cosmic_fusion}"]
    elif config.mode == RunMode.TARGETED_WITH_READ_THROUGH:
        args += [
            "--k=19",
            "--umi-in-name",
            "--max-proximity-distance=1000",
            "--max-proximity-genes=0",
            f"--cosmic-fusion={config.cosmic_fusion}",
        ]
    else:
        raise Exception(f"invalid mode: {config.mode}")
    return args


def run_go(args: List[str]) -> None:
    bin_label = "//go/src/github.com/grailbio/bio/cmd/bio-fusion"
    util.build([bin_label])
    bin_path = util.go_executable(bin_label)
    logging.info("Start: go: %s %s", bin_path, " ".join(args))
    subprocess.check_call([str(bin_path)] + args)


def gene_list_path(config: Config) -> Path:
    if not config.cosmic_fusion:
        return Path(f"{config.cache_dir}/gene-names-denovo.txt")
    return Path(f"{config.cache_dir}/gene-names-targeted.txt")


CPP_WDIR = Path("/tmp/fusion-cpp")


def cpp_args(config: Config) -> List[str]:
    args = [
        f"--r1_fqgz={config.r1}",
        f"--r2_fqgz={config.r2}",
        f"--transcript_file={config.transcript}",
        f"--wdir={str(CPP_WDIR)}",
        "--p=64",
        "--stderrthreshold=0",
        "--v=1",
    ]
    if config.mode == RunMode.DENOVO:
        args += [
            "--k=19",
            "--denovo",
            "--umi_in_name",
            "--cap=2",
            "--unstranded_library",
        ]
    elif config.mode == RunMode.TARGETED:
        args += [
            "--k=19",
            "--umi_in_name",
            f"--cosmic_fusion_file={config.cosmic_fusion}",
            "--unstranded_library",
        ]
    elif config.mode == RunMode.TARGETED_WITH_READ_THROUGH:
        args += [
            "--k=19",
            "--umi_in_name",
            "--proximity_dist=1000",
            "--proximity_num=0",
            f"--cosmic_fusion_file={config.cosmic_fusion}",
            "--unstranded_library",
        ]
    else:
        raise Exception(f"invalid mode: {config.mode}")
    return args


def run_cpp(args: List[str]) -> None:
    CPP_WDIR.mkdir(parents=True, exist_ok=True)
    bin_label = "//bio/rna/fusion:target_rna_fusion"
    util.build([bin_label])
    bin_path = util.nongo_executable(bin_label)
    logging.info("Start: c++: %s %s", bin_path, " ".join(args))
    subprocess.check_call([str(bin_path)] + args)


# For small tests on adhoc:
#
# ./run.py --cosmic_fusion=/tmp/bio-fusion-bench/small/small_pairs.txt --transcript=/tmp/bio-fusion-bench//small/transcriptome.fa --r1=/tmp/bio-fusion-bench//small/smallr1.fastq.gz --r2=/tmp/bio-fusion-bench//small/smallr2.fastq.gz --cache_dir=/tmp/fusion_cache
def main() -> None:
    logging.basicConfig(
        level=logging.DEBUG, format="%(asctime)s:%(levelname)s: %(message)s"
    )

    p = argparse.ArgumentParser()
    p.add_argument(
        "--cache_dir", default="/scratch-nvme/cache_tmp", help="Benchmark cache dir"
    )
    p.add_argument(
        "--cosmic_fusion",
        default="/scratch-nvme/fusion/all_pair.txt",
        help="Fusion event file. If empty, af4 runs in the denovo mode",
    )
    p.add_argument(
        "--transcript",
        default="/scratch-nvme/fusion/transcriptome.fa",
        help="Path to the transcript file",
    )
    p.add_argument("--run", action="append", choices=["go", "c++"])
    p.add_argument("--r1", default="/scratch-nvme/fusion/32-1_S1_L001_R1_001.fastq.gz")
    p.add_argument("--r2", default="/scratch-nvme/fusion/32-1_S1_L001_R2_001.fastq.gz")
    p.add_argument("--update_gene_list", default=True, action="store_true")
    args = p.parse_args()
    if not args.run:
        args.run = ["go", "c++"]
    config = Config(
        mode=RunMode.TARGETED_WITH_READ_THROUGH,
        cache_dir=Path(args.cache_dir),
        run=args.run,
        cosmic_fusion=Path(args.cosmic_fusion) if args.cosmic_fusion else None,
        transcript=Path(args.transcript),
        r1=args.r1,
        r2=args.r2,
        update_gene_list=args.update_gene_list,
    )

    output_suffix = "-" + os.path.basename(os.path.splitext(config.r1)[0])
    if not config.cosmic_fusion:
        output_suffix += "-denovo"
    else:
        output_suffix += "-" + os.path.basename(
            os.path.splitext(config.cosmic_fusion)[0]
        )

    if "go" in config.run:
        gene_list = gene_list_path(config)
        if config.update_gene_list:
            argv = [
                f"-transcript={config.transcript}",
                f"-gene-list-output={gene_list}",
            ]
            if config.mode != RunMode.DENOVO:
                argv.append(f"-cosmic-fusion={config.cosmic_fusion}")
            run_go(argv)
            sort_file(gene_list, gene_list)
        run_go(go_args(config, gene_list, output_suffix))
        sort_fasta_headers(
            GO_WDIR / f"all{output_suffix}.fa",
            GO_WDIR / f"all{output_suffix}.fa.sorted",
        )
        sort_fasta_headers(
            GO_WDIR / f"filtered{output_suffix}.fa",
            GO_WDIR / f"filtered{output_suffix}.fa.sorted",
        )

    if "c++" in config.run:
        for path in glob.glob(str(CPP_WDIR) + "/*"):
            logging.info("Remove " + path)
            os.remove(path)
        run_cpp(cpp_args(config))
        all_fa = list(CPP_WDIR.glob("fusion_[0-9]*.fa"))[0]
        sort_fasta_headers(all_fa, CPP_WDIR / f"all{output_suffix}.fa.sorted")
        filtered_fa = list(CPP_WDIR.glob("*final*[0-9].fa"))[0]
        sort_fasta_headers(filtered_fa, CPP_WDIR / f"filtered{output_suffix}.fa.sorted")

    try:
        logging.info("Diffing the 1st stage outputs")
        subprocess.check_call(
            [
                "diff",
                GO_WDIR / f"all{output_suffix}.fa.sorted",
                CPP_WDIR / f"all{output_suffix}.fa.sorted",
            ]
        )
    except Exception as e:
        logging.info(e)

    try:
        logging.info("Diffing the 2nd stage outputs")
        subprocess.check_call(
            [
                "diff",
                GO_WDIR / f"filtered{output_suffix}.fa.sorted",
                CPP_WDIR / f"filtered{output_suffix}.fa.sorted",
            ]
        )
    except Exception as e:
        logging.info(e)

    logging.info("All done")


main()
