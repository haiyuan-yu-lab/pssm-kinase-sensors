from pathlib import Path
from typing import Tuple
import pandas as pd
import numpy.typing as npt


def read_pssm_file(pssm_file: Path) -> Tuple[npt.NDArray]:
    assert pssm_file.is_file()
    df = pd.read_csv(pssm_file)
    df.columns = ["aa", '-5', '-4', '-3', '-2', '-1', '0', '1', '2', '3', '4']
    aas = df["aa"].values
    pssm = df[['-5', '-4', '-3', '-2', '-1', '0', '1', '2', '3', '4']].values
    return pssm, aas
