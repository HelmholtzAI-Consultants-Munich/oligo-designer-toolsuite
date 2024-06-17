############################################
# imports
############################################

import networkx as nx
import numpy as np
import pandas as pd

from oligo_designer_toolsuite.database import OligoDatabase

############################################
# Heuristic Selection Methods
############################################


def heuristic_selection_independent_set(
    oligos_scores: pd.Series,
    overlapping_matrix: pd.DataFrame,
    n_oligo: int,
    ascending: bool,
    n_trials=100,
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

    # sort the oligos by their score
    oligos_sorted = oligos_scores.sort_values(ascending=ascending)
    oligo_ids_sorted = oligos_sorted.index.tolist()

    # overlapping matrix must have a consistent order
    mat_sorted = overlapping_matrix.loc[oligo_ids_sorted, oligo_ids_sorted].values

    # Represent overlap matrix as graph
    G = nx.convert_matrix.from_numpy_array(mat_sorted)
    # First check if there are no cliques with n oligos
    cliques = nx.algorithms.clique.find_cliques(G)

    # initialize best_idx_set with arbitrary set of non-overlapping oligos with minimum n_oligo oligos
    for clique in cliques:
        if len(clique) >= n_oligo:
            best_idx_set = clique[:n_oligo]
            break

    # initialize max_score with score from set chosen above
    max_score = np.max(oligos_sorted.values[best_idx_set])

    for first_idx in range(min(len(oligo_ids_sorted), n_trials)):
        # use the integer index because the matric is a np array
        set_idxs = np.array([first_idx])
        for _ in range(n_oligo - 1):
            # find first oligo in sorted array that is not overlapping with any selected oligo
            no_overlap = np.all(mat_sorted[set_idxs], axis=0)
            if np.any(no_overlap):
                set_idxs = np.append(set_idxs, np.where(no_overlap)[0][0])
            else:
                break
        if len(set_idxs) == n_oligo:
            score = np.max(oligos_sorted.values[set_idxs])
            if score < max_score:
                max_score = score
                best_idx_set = set_idxs
    best_set = oligos_sorted.iloc[best_idx_set]

    for oligo_id in oligos_scores.index:
        if oligos_scores[oligo_id] > max_score:
            # delete the oligo
            oligos_scores.drop(oligo_id, inplace=True)

    return oligos_scores, best_set
