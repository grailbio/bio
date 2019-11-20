#!/usr/bin/env python3

"""
Script that parses INFO log from bio-target-rna-fusion and shows key stats.

Usage:

"""

import os
import glob
import sys
from pathlib import Path

import util

if len(sys.argv) <= 1:
    raise Exception("Usage: runtime_stats dirglob...")

name_map = {
    "170206_ARTLoD_B1_01rerun": "T1 (0.0001)",
    "170206_ARTLoD_B1_02rerun": "T2 (0.0002)",
    "170206_ARTLoD_B1_03rerun": "T3 (0.0002)",
    "170206_ARTLoD_B1_04rerun": "T4 (0.0004)",
    "170206_ARTLoD_B1_05rerun": "T5 (0.0004)",
    "170206_ARTLoD_B1_06rerun": "T6 (0.0004)",
    "170206_ARTLoD_B1_07rerun": "T7 (0.0006)",
    "170206_ARTLoD_B1_08rerun": "T8 (0.0006)",
    "170206_ARTLoD_B1_09rerun": "T9 (0.0006)",
    "170206_ARTLoD_B1_10rerun": "T10 (0.0008)",
    "170206_ARTLoD_B1_11rerun": "T11 (0.0008)",
    "170206_ARTLoD_B1_12rerun": "T12 (0.0008)",
    "170206_ARTLoD_B1_13rerun": "T13 (0.01)",
    "170206_ARTLoD_B1_14rerun": "T14 (0.01)",
    # RNA datasets
    "101CPREL277": "Prc101",
    "108CPREL315": "Prc108",
    "74HPREL332": "HC332",
    "118HPREL322": "HC118",
    "68HPREL273": "HC273",  # not used
    "53HPREL160": "HC160",
    "117HPREL321": "HC117",  # not used
}


def pretty_name(path: str) -> str:
    if path[-1] == "/":
        name = os.path.basename(path[:-1])
    else:
        name = os.path.basename(path)
    for key, val in name_map.items():
        name = name.replace(key, val)
    return name


for pattern in sys.argv[1:]:
    for path in sorted(glob.glob(pattern)):
        s = util.run_stats(Path(path))
        name = util.pretty_sample_name(path)

        n_fragment_frac = []
        for v in s.n_fragment_matches:
            n_fragment_frac.append(100 * float(v) / float(s.n_fragments))

        print(name, s, "fragment_matches", n_fragment_frac)
        # print(f'{name} & {n_reads_str} & {s.all_candidates - s.min_span} & {s.low_complexity_substring} & {s.close_proximity} & {s.duplicates} & {s.abundant_partners} & {s.final_candidates} & {s.duration} & {s.fusion_stats}\\\\')
