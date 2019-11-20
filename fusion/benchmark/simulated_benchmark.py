#!/usr/bin/env python3

import argparse
import collections
import logging
from pathlib import Path
import os
from typing import List, NamedTuple, Optional, Set, Counter

import util


class GenePair:
    def __init__(self, pair: str) -> None:
        """Arg pair is a two gene names separated by '/' or '--'. For example "EXOSC7/KIAA1328" or EXOSC7--KIAA1328."""
        g1: str
        g2: str
        if "/" in pair:
            g1, g2 = pair.split("/")
        else:
            g1, g2 = pair.split("--")
        if g1 > g2:
            g1, g2 = g2, g1
        self.g1 = g1
        self.g2 = g2

    def __hash__(self):
        return hash(self.g1) + hash(self.g2)

    def __eq__(self, other) -> bool:
        return self.__dict__ == other.__dict__


Score = NamedTuple("Score", [("tp", int), ("fp", int), ("fn", int)])


class TargetedFusionStats:
    def __init__(self, input_targets: Path, caller_fa: Path) -> None:
        """
        Reads a target file to return tuples of pairs.

        :param str input_targets: Target file
        """
        assert os.path.exists(input_targets), input_targets
        assert os.path.exists(caller_fa), caller_fa
        self.sorted_targets: Set[GenePair] = set()
        with open(input_targets) as i_f:
            for line in i_f:
                if line.startswith("Gene"):
                    continue
                gp = GenePair(line.strip().split()[0])
                self.sorted_targets.add(GenePair(line.strip().split()[0]))

        self.sorted_fixed_targets = self.sorted_targets.copy()
        # logging.info('Adding reference gene pairs %s', self.sorted_fixed_targets)

        # KIAA1328/EXO7C manifests as KIAA1328/CLEC3B
        self.sorted_fixed_targets.remove(GenePair("EXOSC7/KIAA1328"))
        self.sorted_fixed_targets.add(GenePair("CLEC3B/KIAA1328"))
        # STAG3L1 is a pseudogene of STAG3
        self.sorted_fixed_targets.remove(GenePair("CCDC88C/STAG3L1"))
        self.sorted_fixed_targets.add(GenePair("CCDC88C/STAG3"))

        # fragments that have one fusion event called.
        unique_calls: Counter[GenePair] = collections.Counter()

        # Fragments with >1 fusion event called.  For such fragments, we pick
        # the the event with the most support and assign the fragment to it.
        multi_calls: List[List[GenePair]] = []

        with open(caller_fa) as i_f:
            for line in i_f:
                if line.startswith(">"):
                    gene_pairs = line.strip().split("|")[1].split(",")
                    if len(gene_pairs) == 1:
                        unique_calls[GenePair(gene_pairs[0])] += 1
                    else:
                        multi_calls.append([GenePair(x) for x in gene_pairs])

        # Fix the genepair -> frequency counts while assigning multicalls to
        # unique calls, to make the math order-independent.
        org_unique_calls = unique_calls.copy()
        for mc in multi_calls:
            best_count = -1
            best_call: Optional[GenePair] = None

            for gp in mc:
                freq = org_unique_calls.get(gp, 0)
                if freq > best_count:
                    best_count = freq
                    best_call = gp
            assert best_call
            unique_calls[best_call] += best_count

        self.calls = unique_calls

    def stats(self, threshold=2) -> Score:
        """
        Prints some simple stats about a caller's predictions wrt the targets
        """
        tp, fp, fn = 0, 0, 0
        calls_above_threshold = {x: y for x, y in self.calls.items() if y >= threshold}
        for call, count in calls_above_threshold.items():
            if call in self.sorted_fixed_targets:
                tp += 1
            else:
                fp += 1

        for call in self.sorted_fixed_targets:
            if call not in calls_above_threshold:
                fn += 1
        return Score(tp=tp, fn=fn, fp=fp)


def main() -> None:
    logging.basicConfig(
        level=logging.DEBUG, format="%(asctime)s:%(levelname)s: %(message)s"
    )

    p = argparse.ArgumentParser()
    p.add_argument(
        "--cache_dir", default=util.DEFAULT_CACHE_DIR, help="Benchmark cache dir"
    )
    p.add_argument(
        "--result_dir", default=util.DEFAULT_RESULT_DIR, help="Benchmark result dir"
    )
    p.add_argument(
        "--rerun_af4",
        action="store_true",
        help="Always run AF4 even if the result file already exists",
    )
    p.add_argument(
        "--recache_files",
        action="store_true",
        help="Always copy benchmark data files, even if they already exist locally.",
    )
    args = p.parse_args()
    util.s3_cache_files(
        [
            util.REFERENCE_DIR
            + "/gencode.v26.250padded_separate_jns_transcripts_parsed_no_mt_no_overlap_no_pary_no_versioned.fa",
            util.REFERENCE_DIR + "/all_pair_art_lod_gpair_merged.txt",
            util.REFERENCE_DIR + "/liu_gpair.txt",
        ],
        args.cache_dir,
    )
    for mode in ["denovo", "targeted"]:
        for sample in util.SIMULATED_SAMPLES:
            util.s3_cache_files([sample.path.r1, sample.path.r2], args.cache_dir)
            result_dir = (
                f"{args.result_dir}/synthetic-{mode}-{sample.n}-{sample.coverage}"
            )
            try:
                os.makedirs(result_dir, 0o755)
            except:
                logging.error("mkdir %s failed", result_dir)
            if not os.path.exists(f"{result_dir}/filtered.fa") or args.rerun_af4:
                logging.info("running benchmark in %s", result_dir)
                af4_args = [
                    str(util.af4_path()),
                    f"-log_dir={result_dir}",
                    f"-r1={args.cache_dir}/{sample.path.r1}",
                    f"-r2={args.cache_dir}/{sample.path.r2}",
                    f"-fasta-output={result_dir}/all.fa",
                    f"-filtered-output={result_dir}/filtered.fa",
                    f"-max-genes-per-kmer=2",
                    f"-max-proximity-distance=1000",
                    f"-max-proximity-genes=5",
                    "-transcript="
                    + args.cache_dir
                    + "/gencode.v26.250padded_separate_jns_transcripts_parsed_no_mt_no_overlap_no_pary_no_versioned.fa",
                ]
                if mode == "targeted":
                    af4_args.append(
                        "-cosmic-fusion="
                        + args.cache_dir
                        + "/all_pair_art_lod_gpair_merged.txt"
                    )
                util.check_call(af4_args)
                logging.info("Runtime stats: %s", util.run_stats(Path(result_dir)))

            stats = TargetedFusionStats(
                Path(f"{args.cache_dir}/liu_gpair.txt"),
                Path(f"{result_dir}/filtered.fa"),
            )

            s = stats.stats()
            tp = "%d" % (s.tp,)
            fp = "%d" % (s.fp,)
            fn = "%d" % (s.fn,)
            print(f"{mode} & {sample.n} & {sample.coverage} & {tp} & {fp} & {fn}\\\\")


if __name__ == "__main__":
    main()
