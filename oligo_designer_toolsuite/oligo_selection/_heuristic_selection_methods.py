############################################
# imports
############################################

import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix

############################################
# Heuristic Selection Methods
############################################


def heuristic_selection_independent_set(
    oligoset_init: list,
    oligos_scores: pd.Series,
    overlapping_matrix: csr_matrix,
    overlapping_matrix_ids: list,
    n_oligo: int,
    ascending: bool,
    n_trials: int = 100,
):
    """
    This method empirically finds an optimal set of non-overlapping oligos based on their scores. It iteratively
    selects the best non-overlapping oligos to form a set, evaluates their scores, and discards oligos unlikely
    to produce a better set than the best found set.

    :param oligos_scores: Series containing scores of oligos.
    :type oligos_scores: pd.Series
    :param overlapping_matrix: Dataframe indicating overlap between oligos.
    :type overlapping_matrix: pd.DataFrame
    :param n_oligo: The number of oligos to select in each set.
    :type n_oligo: int
    :param ascending: Determines if scores should be sorted in ascending order, defaults to True.
    :type ascending: bool
    :param n_trials: Number of top scoring oligos to consider for set formation, defaults to 100.
    :type n_trials: int
    :return: Updated scores, and the best set found.
    :rtype: tuple(dict, pd.Series, pd.Series)
    """

    # Sort the oligos by their score
    oligo_ids_sorted = oligos_scores.sort_values(ascending=ascending).index.to_list()

    # Sort overlap matrix by oligo scores, i.e. best performing oligo at first entry
    overlapping_matrix_indices = [overlapping_matrix_ids.index(oligo_id) for oligo_id in oligo_ids_sorted]
    overlapping_matrix_sorted = overlapping_matrix[overlapping_matrix_indices, :][
        :, overlapping_matrix_indices
    ]

    # Initialize max_score with score from initial oligoset
    best_oligoset = oligoset_init
    best_oligoset_max_score = oligos_scores.loc[oligoset_init].max()

    for first_idx in range(min(len(oligo_ids_sorted), n_trials)):
        # Use the integer index because the matrix is converted to np array
        oligoset_idxs = np.array([first_idx])
        for _ in range(n_oligo - 1):
            # Find first oligo in sorted array that is not overlapping with any selected oligo
            no_overlap = np.all(overlapping_matrix_sorted[oligoset_idxs].toarray(), axis=0)
            # Ensure not to select already selected oligos
            # no_overlap[oligoset_idxs] = False
            # Add the first entry withou overlap to set, i.e. the next oligo with the best score
            if np.any(no_overlap):
                oligoset_idxs = np.append(oligoset_idxs, np.where(no_overlap)[0][0])
            else:
                break

        # If enough non overlapping oligos are found and score is lower than existing set, replace existing with new set
        if len(oligoset_idxs) == n_oligo:
            oligoset = [oligo_ids_sorted[idx] for idx in oligoset_idxs]
            oligoset_max_score = oligos_scores.loc[oligoset].max()
            if oligoset_max_score < best_oligoset_max_score:
                best_oligoset_max_score = oligoset_max_score
                best_oligoset = oligoset

    oligos_scores = oligos_scores[oligos_scores <= best_oligoset_max_score]

    return best_oligoset, oligos_scores
