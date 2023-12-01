import heapq
import numpy as np
import numpy.typing as npt
from typing import Iterator, List, Tuple


def is_valid_sequence(sequence: List[str] | str) -> bool:
    """
    checks if the provided sequence complies with the following conditions:
        - aminoacid at position 0 is always S
        - No Cys(C) anywhere in the sequence
        - No K in positions -1 and 1
        - If a K/R appears in a negative position, it won't appear on a
          positive position, and vice-versa
    """
    if "C" in sequence:
        return False
    if len(sequence) != 10 or sequence[5] != "S" or "K" in sequence[4:7]:
        return False
    first_half = sequence[:5]
    second_half = sequence[6:]
    if (("K" in first_half or "R" in first_half) and
            ("K" in second_half or "R" in second_half)):
        return False
    return True


def greedy_pssm_generator(pssm: npt.NDArray,
                          aminoacids: List[str],
                          descending: bool = True) -> Iterator[str]:
    """
    Generates amino acid sequences for the sensor, indexed from -5 to 4 in
    score order.

    Parameters
    ----------
    pssm : 2D matrix
        every row corresponds to a specific aminoacid, and each column of the
        array corresponds to the position it would occupy in the kinase.
    aminoacids : list of str
        the order of the aminoacids, corresponding to every row in `pssm`
    descending: bool, default True
        if set to true, the greddy direction is ascending for the score.

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

    r_size, c_size = pssm.shape

    # sort sequences per position
    greedy_sort = pssm.argsort(axis=0)
    if descending:
        greedy_sort = np.flip(pssm.argsort(axis=0), axis=0)
    cols = list(range(c_size))
    initial_indices = [0] * c_size  # start with the maximum possible score
    pssm_indices = greedy_sort[initial_indices, cols]
    init_vals = pssm[pssm_indices, cols]
    init_seq = [aminoacids[i] for i in pssm_indices]
    init_seq[5] = "S"
    init_score = np.log2(np.prod(init_vals))
    if descending:
        init_score = -init_score
    init_state = (init_score, tuple(init_seq), tuple(initial_indices))

    pq = [init_state]
    visited = set([tuple(pssm_indices)])

    while pq:
        score, seq, indices = heapq.heappop(pq)
        if is_valid_sequence(seq):
            yield "".join(seq)
        for col in cols:
            if col == 5:
                continue
            if indices[col] + 1 < r_size:
                new_indices = list(indices)
                new_indices[col] += 1

                if tuple(new_indices) not in visited:
                    visited.add(tuple(new_indices))
                    pssm_indices = greedy_sort[new_indices, cols]
                    new_vals = pssm[pssm_indices, cols]
                    new_seq = [aminoacids[i] for i in pssm_indices]
                    new_seq[5] = "S"
                    new_score = np.log2(np.prod(new_vals))
                    if descending:
                        new_score = -new_score
                    heapq.heappush(pq, (new_score,
                                        tuple(new_seq),
                                        tuple(new_indices)))


def get_sequence_score(pssm: npt.NDArray,
                       aminoacids: List[str] | str,
                       sequence: List[str] | str) -> float:
    indices = []
    for c in sequence:
        indices.append(aminoacids.tolist().index(c))
    vals = pssm[indices, list(range(pssm.shape[1]))]
    return np.log2(np.prod(vals))


def get_pssm_range(pssm: npt.NDArray,
                   aminoacids: List[str] | str) -> Tuple[float]:
    max_seq = next(greedy_pssm_generator(pssm, aminoacids))
    min_seq = next(greedy_pssm_generator(pssm, aminoacids, descending=False))
    return (get_sequence_score(pssm, aminoacids, min_seq),
            get_sequence_score(pssm, aminoacids, max_seq))
