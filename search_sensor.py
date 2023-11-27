from pathlib import Path
from typing import List
from peptide_search.tools.sequence import (greedy_pssm_generator,
                                           get_sequence_score,
                                           get_pssm_range)
from peptide_search.io.readers import read_pssm_file
import numpy as np


def run(kinases_paths: List[Path], output_dir: Path) -> None:
    assert all(k.is_file() for k in kinases_paths)
    assert output_dir.is_dir()
        


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="use multiple PSSMs to generate sensors in a greedy"
                    " fashion")
    parser.add_argument("-k", "--kinases", required=True,
                        help="List to kinase files separated by comma")
    parser.add_argument("-t", "--target", required=True,
                        help="Path to the kinase file to use as target")
    parser.add_argument("-o", "--output-dir", required=True,
                        help="Path to the output directory")
    args = parser.parse_args()
    run(Path(args.kinases),
        Path(args.output_dir))
