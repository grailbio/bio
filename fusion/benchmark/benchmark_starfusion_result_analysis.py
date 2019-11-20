#!/usr/bin/evn python3

import logging
import argparse
import os
import sys
from pathlib import Path
from typing import NamedTuple, Any, List, Set

from simulated_benchmark import GenePair
from util import s3_cache_files, REFERENCE_DIR


def read_fusion_pair(p: Path) -> Set[GenePair]:
    sorted_targets: Set[GenePair] = set()
    with open(p) as i_f:
        for line in i_f:
            if line.startswith("Gene") or line.startswith("#"):
                continue
            sorted_targets.add(GenePair(line.strip().split()[0]))
    return sorted_targets


def main() -> None:
    logging.basicConfig(
        level=logging.DEBUG, format="%(asctime)s:%(levelname)s: %(messge)s"
    )
    p = argparse.ArgumentParser()

    p.add_argument(
        "--starfusion_dir",
        default="/scratch-nvme/xyang/result/",
        help="Starfusion result dir, which then contains individual sample result",
    )
    p.add_argument(
        "--cache_dir", default="/scratch-nvme/xyang/cache/", help="cache dir"
    )

    args = p.parse_args()
    local_truth_gpair = f"{args.cache_dir}/liu_gpair.txt"
    if not os.path.exists(local_truth_gpair):
        s3_cache_files([REFERENCE_DIR + "/liu_gpair.txt"], args.cache_dir)

    listOfFiles = list()
    for (dirpath, dirnames, filenames) in os.walk(args.starfusion_dir):
        listOfFiles += [os.path.join(dirpath, file) for file in filenames]

    # get all relevant starfusion result files
    filtered_results = [
        f
        for f in listOfFiles
        if f.find("star-fusion.fusion_predictions.abridged.tsv") != -1
    ]

    # compare with truth set
    true_gpairs = read_fusion_pair(local_truth_gpair)

    npos = len(true_gpairs)

    for f in filtered_results:
        print(f)
        results = read_fusion_pair(f)

        npredict = len(results)
        tp = len(true_gpairs.intersection(results))
        fn = npos - tp
        fp = npredict - tp

        precision = tp / (tp + fp)
        recall = tp / npos
        f1 = 2 * precision * recall / (precision + recall)
        f1 = "%.3f" % (f1)
        print(f"tp={tp}, fn={fn}, fp={fp}, f1={f1}")


main()
