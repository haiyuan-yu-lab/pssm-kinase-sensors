from pathlib import Path
from typing import List, Dict
from peptide_search.tools.sequence import (greedy_pssm_generator,
                                           get_sequence_score,
                                           get_pssm_range)
from peptide_search.io.readers import read_pssm_file
import numpy as np
from rich.progress import track


def get_score_distribution(kinases: List[Dict], sequence: str) -> Dict:
    scores = {}
    for kinase, kd in kinases.items():
        scores[kinase] = get_sequence_score(kd["pssm"],
                                            kd["aminoacids"],
                                            sequence)
    return scores


def run(kinases_directory: Path,
        kinases: List[str],
        target_kinase: str,
        output_file: Path,
        separation_threshold=0.5,
        background_threshold=0.5,
        max_search=1000) -> None:
    assert kinases_directory.is_dir()
    assert 0.0 < separation_threshold <= 1.0
    assert 0.0 < background_threshold <= 1.0

    target_file = kinases_directory / f"{target_kinase}.csv"
    assert target_file.is_file()

    sub_targets_files = []
    for stk in kinases:
        stk_file = kinases_directory / f"{stk}.csv"
        assert stk_file.is_file()
        sub_targets_files.append(stk_file)

    background_files = []
    for bgk in kinases_directory.glob("*.csv"):
        if bgk.stem not in kinases and bgk.stem != target_kinase:
            background_files.append(bgk)

    subtargets = {}
    for st in sub_targets_files:
        pssm, aminoacids = read_pssm_file(st)
        subtargets[st.stem] = {
            "pssm": pssm,
            "aminoacids": aminoacids
        }
    del sub_targets_files

    background = {}
    for bf in background_files:
        pssm, aminoacids = read_pssm_file(bf)
        background[bf.stem] = {
            "pssm": pssm,
            "aminoacids": aminoacids
        }
    del background_files

    pssm, aminoacids = read_pssm_file(target_file)
    target_min, target_max = get_pssm_range(pssm, aminoacids)
    st_distance = separation_threshold * (target_max - target_min)
    bg_distance = background_threshold * (target_max - target_min)

    sequence_generator = greedy_pssm_generator(pssm, aminoacids)

    with output_file.open("w") as of:
        of.write("sequence\ttarget:pssm_score\tsubtarget:pssm_score\t"
                 "subtarget_range (min:max)\tbackground_range (min:max)\n")
        for i in track(range(max_search),
                       description="Generating sequences..."):
            seq = next(sequence_generator)
            score = get_sequence_score(pssm, aminoacids, seq)

            st_scores = get_score_distribution(subtargets, seq)
            st_scores_max = max(s for s in st_scores.values())
            st_scores_min = min(s for s in st_scores.values())

            if score - st_scores_max < st_distance:
                print("not valid")

            bg_scores = get_score_distribution(background, seq)
            bg_scores_max = max(s for s in bg_scores.values())
            bg_scores_min = min(s for s in bg_scores.values())

            if score - bg_scores_max < bg_distance:
                print("not valid (BG)")

            subtarget_string = "|".join(f"{k}:{s}"
                                        for k, s in st_scores.items())
            of.write(f"{seq}\t{target_kinase}:{score}\t{subtarget_string}\t"
                     f"{st_scores_min}:{st_scores_max}\t"
                     f"{bg_scores_min}:{bg_scores_max}\n")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="use multiple PSSMs to generate sensors in a greedy"
                    " fashion")
    parser.add_argument("-d", "--kinases-directory", required=True,
                        help="List to a directory of csv files to load kinases"
                             " from. File should be named <kinase-name>.csv")
    parser.add_argument("-k", "--kinases", required=True,
                        help="List to kinase names separated by spaces, should"
                             " not include the target kinase",
                        nargs="+")
    parser.add_argument("-t", "--target", required=True,
                        help="kinase name to use as target")
    parser.add_argument("-s", "--separation-threshold",
                        default=0.5, type=float,
                        help="The expected minimum separation threshold"
                             " between the target kinase and the sub-targets")
    parser.add_argument("-b", "--background-threshold",
                        default=0.5, type=float,
                        help="The expected minimum separation threshold"
                             " between the target kinase and the background")
    parser.add_argument("-m", "--max-search",
                        default=1000, type=int,
                        help="Maximum number of sequences to generate before"
                             " forcifully stopping the script.")
    parser.add_argument("-o", "--output-file", required=True,
                        help="Path to the output directory")
    args = parser.parse_args()
    run(Path(args.kinases_directory),
        args.kinases,
        args.target,
        Path(args.output_file),
        args.separation_threshold,
        args.background_threshold,
        args.max_search)
