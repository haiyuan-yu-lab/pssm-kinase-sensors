import numpy as np
import numpy.typing as npt
from typing import Iterator, List


def is_valid_sequence(sequence: str) -> bool:
    """
    checks if the provided sequence complies with the following conditions:
        - aminoacid at position 0 is always S
        - No Cys(C) anywhere in the sequence
        - No K in positions -1 and 1
        - If a K/R appears in a negative position, it won't appear on a
          positive position, and vice-versa
    """
    if len(sequence) != 10 or sequence[5] != "S" or "K" in sequence[4:7]:
        return False
    first_half = sequence[:5]
    second_half = sequence[6:]
    if (("K" in first_half or "R" in first_half) and
            ("K" in second_half or "R" in second_half)):
        return False
    return True


def greedy_pssm_generator(pssm: npt.NDArray,
                          aminoacids: List[str]) -> Iterator[str]:
    """
    Generates amino acid sequences for the sensor, indexed from -5 to 4.
    Parameters
    ----------
    pssm : 2D matrix
        every row corresponds to a specific aminoacid, and each column of the
        array corresponds to the position it would occupy in the kinase.
    aminoacids : list of str
        the order of the aminoacids, corresponding to every row in `pssm`

    Returns
    -------
    Iterator of strings with the aminoacids in descending sequence score

    Notes
    -----

    The following conditions are applied to the generated sequence:

        - aminoacid at position 0 is always S
        - No Cys(C) anywhere in the sequence
        - No K in positions -1 and 1
        - If a K/R appears in a negative position, it won't appear on a
          positive position, and vice-versa
    """

    # sort sequences per position
    greedy_sort = np.flip(pssm.argsort(axis=0), axis=0)
    cur_positions = [0] * 9  # start with the maximum possible score
    cur_scores = np.zeros(9)
    while max(cur_positions) < len(aminoacids):
        sequence = []
        for col in range(len(cur_positions)):
            idx = greedy_sort[cur_positions[col], col]
            cur_scores[col] = pssm[cur_positions[col], col]
            sequence.append(aminoacids[idx])
        sequence.insert(5, "S")
        sequence = "".join(sequence)
        if is_valid_sequence(sequence):
            yield sequence
        # advance to the next position in a gredy fashion
        current_score = np.log2(np.prod(cur_scores))
        scores = []
        for i, position in enumerate(cur_positions):
            if position >= len(aminoacids) - 1:
                continue
            score = pssm[position+1, i]
            putative_scores = cur_scores[:i] + [score] + cur_scores[i+1:]
            putative_score = np.log2(np.prod(putative_scores))
            scores.append((i, current_score - putative_score))
        if len(scores) > 0:


            



def get_sequence_score(pssm: npt.NDArray, sequence: str):
