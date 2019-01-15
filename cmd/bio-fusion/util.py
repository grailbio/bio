import logging
import glob
import os
import subprocess
import re
from pathlib import Path
from typing import List, NamedTuple, Optional

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
        af4_label = '//go/src/grail.com/cmd/bio-target-rna-fusion'
        build([af4_label])
        AF4_PATH = go_executable(af4_label)
    return AF4_PATH

def s3_ls(path: str) -> List[str]:
    out = check_output([str(grail_file_path()), 'ls', '-R', path])
    return [x for x in out.split('\n') if x != '']

def s3_cache_files(src_paths: List[str], cache_dir: Path, force=False) -> None:
    args = [str(grail_file_path()), 'cp', '-v']
    n  = 0
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
