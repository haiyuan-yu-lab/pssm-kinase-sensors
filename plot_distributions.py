from pathlib import Path
from typing import List, Dict
from peptide_search.tools.sequence import (get_sequence_score,
                                           get_pssm_range)
from peptide_search.io.readers import read_pssm_file
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib import colors
import seaborn as sns


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
        sequences: List[str]) -> None:
    assert kinases_directory.is_dir()

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

    assert all(len(s) == pssm.shape[1] for s in sequences)

    target_min, target_max = get_pssm_range(pssm, aminoacids)
    data = {k: [] for k in ["kinase", "seq", "score", "group"]}
    for seq in sequences:
        score = get_sequence_score(pssm, aminoacids, seq)

        data["kinase"].append(target_kinase)
        data["seq"].append(seq)
        data["score"].append(score)
        data["group"].append("target")

        st_scores = get_score_distribution(subtargets, seq)

        for st_kinase, st_score in st_scores.items():
            data["kinase"].append(st_kinase)
            data["seq"].append(seq)
            data["score"].append(st_score)
            data["group"].append("subtarget")

        bg_scores = get_score_distribution(background, seq)

        for bg_kinase, bg_score in bg_scores.items():
            data["kinase"].append(bg_kinase)
            data["seq"].append(seq)
            data["score"].append(bg_score)
            data["group"].append("background")

    df = pd.DataFrame(data)
    del data

    cols = {
        "target": colors.to_rgba("crimson"),
        "subtarget": colors.to_rgba("gold"),
        "background": colors.to_rgba("lightsteelblue"),
    }

    legend_elements = [Line2D([0], [0], marker='o', color='w',
                              label='Target',
                              markerfacecolor=cols["target"], markersize=9),
                       Line2D([0], [0], marker='o', color='w',
                              label='Sub-target',
                              markerfacecolor=cols["subtarget"], markersize=9),
                       Patch(facecolor=cols["background"], edgecolor='gray',
                             label='Background')]

    fig, ax = plt.subplots(figsize=(2.3 * len(sequences), 7))
    cond = df["group"] == "background"
    sns.violinplot(data=df[cond],
                   x="seq",
                   y="score",
                   color=cols["background"],
                   ax=ax)
    cond = df["group"] != "background"
    sns.swarmplot(data=df[cond],
                  x="seq",
                  y="score",
                  hue="group",
                  palette=cols,
                  size=9,
                  ax=ax)
    ax.set_xlabel("Sequence")
    ax.set_ylabel("Score")
    ax.legend(handles=legend_elements,
              loc=(1.01, 0.5))
    fig.suptitle(target_kinase)
    plt.savefig(output_file, bbox_inches="tight")
    plt.close("all")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="use multiple PSSMs to generate sensors in a greedy"
                    " fashion")
    parser.add_argument("-d", "--kinases-directory", required=True,
                        help="List to a directory of csv files to load kinases"
                             " from. File should be named <kinase-name>.csv")
    parser.add_argument("-k", "--kinases", required=True, nargs="+",
                        help="List to kinase names separated by spaces, should"
                             " not include the target kinase")
    parser.add_argument("-t", "--target", required=True,
                        help="kinase name to use as target")
    parser.add_argument("-s", "--sequences", required=True, nargs="+",
                        help="The sequences to use for calculating the scores")
    parser.add_argument("-o", "--output-file", required=True,
                        help="Path to the output file")
    args = parser.parse_args()
    run(Path(args.kinases_directory),
        args.kinases,
        args.target,
        Path(args.output_file),
        args.sequences)
