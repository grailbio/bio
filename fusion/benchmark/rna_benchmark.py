#!/usr/bin/env python3

"""
Run cfRNA benchmarks. The default parameter settings are for ubuntu*.mpk machines.
Results are stored in /scratch-nvme/af4_benchmark_results

Usage:
rna_benchmark.py --run={af4,starfusion}

"""
import re
import argparse
import logging
from pathlib import Path
import os
from typing import Any, List, Set

import util


def run_starfusion(
    sample_name: str, cached_file_pairs: List[util.FASTQPair], args: Any
):
    match = re.match(r".*/([^/]+)\.FULL\.tar\.gz$", args.starfusion_targz)
    assert match
    local_starfusion_dir = match[1]

    logging.info("LOCAL: %s", local_starfusion_dir)
    if not os.path.exists(os.path.join(args.starfusion_data_dir, local_starfusion_dir)):
        util.check_call(
            ["tar", "xzf", args.starfusion_targz, "-C", args.starfusion_data_dir]
        )
        util.check_call(
            ["make", "-C", os.path.join(args.starfusion_data_dir, local_starfusion_dir)]
        )

    match = re.match(r".*/([^/]+)\.tar\.gz$", args.starfusion_plug_n_play_targz)
    assert match
    local_plugnplay_dir = match[1]
    if not os.path.exists(args.starfusion_data_dir + local_plugnplay_dir):
        util.check_call(
            [
                "tar",
                "xzf",
                args.starfusion_plug_n_play_targz,
                "-C",
                args.starfusion_data_dir,
            ]
        )

    result_dir = args.result_dir + "/" + os.path.basename(sample_name + "-starfusion")
    logging.info("Start starfusion benchmark: %s", result_dir)
    try:
        os.makedirs(result_dir, 0o755)
    except:
        logging.error("mkdir %s failed", result_dir)

    cached_r1 = ",".join(
        [args.cache_dir + "/" + os.path.basename(fp.r1) for fp in cached_file_pairs]
    )
    cached_r2 = ",".join(
        [args.cache_dir + "/" + os.path.basename(fp.r2) for fp in cached_file_pairs]
    )

    starfusion_args = ["docker", "run"]
    mounted: Set[str] = set()
    for dir in [args.starfusion_data_dir, args.result_dir, args.cache_dir]:
        if dir not in mounted:
            mounted.add(dir)
            starfusion_args += ["-v", f"{dir}:{dir}"]
    starfusion_args += [
        "--rm",
        "trinityctat/ctatfusion",
        os.path.join(args.starfusion_data_dir, local_starfusion_dir, "STAR-Fusion"),
        "--left_fq",
        cached_r1,
        "--right_fq",
        cached_r2,
        "--CPU",
        "56",
        "--genome_lib_dir",
        os.path.join(
            args.starfusion_data_dir, local_plugnplay_dir, "ctat_genome_lib_build_dir"
        ),
        "-O",
        result_dir,
        "--FusionInspector",
        "validate",
    ]
    try:
        util.check_call(starfusion_args)
    except Exception as e:
        logging.error("Starfusion failed (ignoring): %s", e)
    logging.info("Finished starfusion benchmark: %s", result_dir)


def run_af4(
    sample_name: str,
    cached_file_pairs: List[util.FASTQPair],
    cosmic_fusion_path: str,
    args: Any,
):
    ref_path = "s3://grail-publications/resources/gencode.v26.whole_genes.fa"
    util.s3_cache_files([ref_path, cosmic_fusion_path], args.cache_dir)

    cached_r1 = ",".join(
        [args.cache_dir + "/" + os.path.basename(fp.r1) for fp in cached_file_pairs]
    )
    cached_r2 = ",".join(
        [args.cache_dir + "/" + os.path.basename(fp.r2) for fp in cached_file_pairs]
    )
    for mode in ["denovo", "targeted"]:
        result_dir = args.result_dir + "/" + os.path.basename(sample_name + "-" + mode)
        if os.path.exists(result_dir + "/filtered.fa"):
            logging.info("Skipping benchmark: %s", result_dir)
            continue
        logging.info("Start af4 benchmark: %s", result_dir)
        try:
            os.makedirs(result_dir, 0o755)
        except:
            logging.error("mkdir %s failed", result_dir)
        af4_args = [
            str(util.af4_path()),
            f"-log_dir={result_dir}",
            f"-pprof=:12345",
            f"-mutex-profile-rate=1000",
            f"-block-profile-rate=1000",
            f"-r1={cached_r1}",
            f"-r2={cached_r2}",
            f"-max-genes-per-kmer=2",
            f"-max-proximity-distance=1000",
            f"-max-proximity-genes=5",
            f"-fasta-output={result_dir}/all.fa",
            f"-filtered-output={result_dir}/filtered.fa",
            f"-transcript={args.cache_dir}/gencode.v26.250padded_separate_jns_transcripts_parsed_no_mt_no_overlap_no_pary_no_versioned.fa",
        ]
        if mode == "targeted":
            af4_args.append(
                f"-cosmic-fusion={args.cache_dir}/all_pair_art_lod_gpair_merged.txt"
            )
        util.check_call(af4_args)
        logging.info("Finished benchmark: %s", result_dir)
        logging.info("Runtime stats: %s", util.run_stats(Path(result_dir)))


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
        "--starfusion_data_dir",
        default="/scratch-nvme/starfusion",
        help="Directory for expanding starfusion plug-n-play files",
    )
    p.add_argument(
        "--run",
        action="append",
        choices=["af4", "starfusion"],
        help="List of systems to run. If unset, run all the configured systems",
    )
    p.add_argument(
        "--starfusion_plug_n_play_targz",
        default=os.environ["HOME"]
        + "/GRCh38_v27_CTAT_lib_Feb092018.plug-n-play.tar.gz",
        help="Tar.gz file of starfusion plug-n-play file. https://github.com/STAR-Fusion/STAR-Fusion/wiki#data-resources-required",
    )
    p.add_argument(
        "--starfusion_targz",
        default=os.environ["HOME"] + "/STAR-Fusion-v1.5.0.FULL.tar.gz",
        help="Tar.gz file of starfusion source package. https://github.com/STAR-Fusion/STAR-Fusion/wiki#data-resources-required",
    )
    p.add_argument(
        "--brca_data_dir",
        default="/scratch-nvme/xyang/brca_rnaseq_data",
        help="BT474, KPL4, MCF7, SKBR3 Breast cancer data directory",
    )

    args = p.parse_args()
    if not args.run:
        args.run = ["af4", "starfusion"]

    ## brca rna-seq for af4
    brca_samples = [
        os.path.join(args.brca_data_dir, s) for s in ["BT474", "KPL4", "MCF7", "SKBR3"]
    ]
    for s in brca_samples:
        if not os.path.exists(os.path.join(args.brca_data_dir, s)):
            util.check_call(
                [
                    "download_brca_data.py",
                    "--odir",
                    "/scratch-nvme/xyang/brca_rnaseq_data",
                ]
            )

    cosmic_fusion_path = (
        "s3://grail-publications/2019-ISMB/references/all_art_lod_brca.txt"
    )
    for sample in brca_samples:
        r1s: List[str] = []
        for fq in os.listdir(sample):
            if "_1" in fq:
                r1s.append(os.path.join(sample, fq))
        cached_file_pairs: List[util.FASTQPair] = []
        for r1 in r1s:
            assert os.path.exists(r1.replace("_1", "_2"))
            cached_file_pairs.append(util.FASTQPair(r1=r1, r2=r1.replace("_1", "_2")))
        print(os.path.basename(sample))
        print(cached_file_pairs)

        run_af4(os.path.basename(sample), cached_file_pairs, cosmic_fusion_path, args)

    ## cfrna for af4 and starfusion
    cosmic_fusion_path = (
        "s3://grail-publications/2019-ISMB/references/all_pair_art_lod_gpair_merged.txt"
    )
    for sample in util.RNA_SAMPLES:
        fastq_files: List[str] = []
        cached_file_pairs: List[util.FASTQPair] = []
        for fp in sample.paths:
            assert fp.r1.replace("R1", "R2") == fp.r2, fp.r2
            fastq_files += [fp.r1, fp.r2]
            cached_file_pairs.append(
                util.FASTQPair(
                    r1=args.cache_dir + "/" + os.path.basename(fp.r1),
                    r2=args.cache_dir + "/" + os.path.basename(fp.r2),
                )
            )
        util.s3_cache_files(fastq_files, args.cache_dir)

        if "af4" in args.run:
            run_af4(sample.name, cached_file_pairs, cosmic_fusion_path, args)
        if "starfusion" in args.run:
            run_starfusion(sample.name, cached_file_pairs, args)


if __name__ == "__main__":
    main()
