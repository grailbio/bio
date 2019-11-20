#!/usr/bin/env python3

"""
Download brca data used in the af4 benchmark.

Usage:
download_brca_data.py --fastq-dump {binary_fastq_dump} --odir {output directory}

"""
import argparse
import logging
import subprocess
import os
from typing import List, NamedTuple

BrcaSample = NamedTuple("BrcaSample", [("sra_id", str), ("name", str)])

Samples = [
    BrcaSample(sra_id="SRR064286", name="MCF7"),
    BrcaSample(sra_id="SRR064287", name="KPL4"),
    BrcaSample(sra_id="SRR064437", name="NORMAL"),
    BrcaSample(sra_id="SRR064438", name="BT474"),
    BrcaSample(sra_id="SRR064439", name="BT474"),
    BrcaSample(sra_id="SRR064440", name="SKBR3"),
    BrcaSample(sra_id="SRR064441", name="SKBR3"),
]


def check_call(args: List[str]) -> None:
    logging.info("Run: %s", " ".join(args))
    subprocess.check_call(args)


def main() -> None:
    logging.basicConfig(
        level=logging.DEBUG, format="%(asctime)s:%(levelname)s: %(message)s"
    )
    p = argparse.ArgumentParser()
    p.add_argument(
        "--fastq_dump",
        default=os.environ["HOME"] + "/sratoolkit.2.9.4-1-ubuntu64/bin/fastq-dump",
        help="sra fastq-dump binary",
    )
    p.add_argument("--odir", default=os.environ["HOME"], help="output directory")
    args = p.parse_args()

    for sample in Samples:
        check_call(["mkdir", "-p", sample.name])
        check_call([args.fastq_dump, "--split-files", "--gzip", sample.sra_id])
        check_call(["mv", sample.sra_id + "*", sample.name])


main()
