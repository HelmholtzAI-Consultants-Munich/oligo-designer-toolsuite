############################################
# imports
############################################

import numpy as np
import networkx as nx

############################################
# Heuristic Selection Methods
############################################

def padlock_heuristic_selection(
    database_region, oligos_scores, overlapping_matrix, n_oligo, n_trials=10000
):
    """This method tries to find empirically a good set of teh padlock design of oligos. For teh best n_trials oligos a set is cretaed
    by adding one by one the best non overlapping oligo untiil a set is created. The best obtained set saved and returned. Moreover, foor how this particular
    scoring method works we know that all the oligos with a higher score than the worse oligo of the best set are going to yield a worse set the one just found,
    hence are deleted form the dictionary and the series of scores.

    :param database_region: dictionary with all the oligos of the region
    :type database_region: dict
    :param oligos_scores: scores of the oligos
    :type oligos_scores: pandas.Series
    :param overlapping_matrix: matrix containig information about the overlapping of the oligos
    :type overlapping_matrix: pandas.DataFrame
    :param n_oligo: size of the set
    :type n_oligo: int
    :param n_trials: number of sets to be tried, defaults to 10000
    :type n_trials: int, optional
    :return: filtered dcit of the oligos, filtered scores fo teh oligos, bes set found
    :rtype: dict, pandas.Series, pandas.Series
    """

    # sort the oligos by their score
    oligos_sorted = oligos_scores.sort_values()
    oligo_ids_sorted = oligos_sorted.index.tolist()
    # overlapping matrix must have a consistent oreder
    mat_sorted = overlapping_matrix.loc[oligo_ids_sorted, oligo_ids_sorted].values

    # already sorted df, the max is the last entry
    # max_score = (oligos_sorted.iloc[-1] * 1.1)  

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

    best_idx_set = []
    for first_idx in range(
        min(len(oligo_ids_sorted), n_trials)
    ):  # use the integer index because the matric is a np array
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
    best_set = oligos_sorted.iloc[best_idx_set]  # adapt with the

    for oligo_id in oligos_scores.index:
        if oligos_scores[oligo_id] > max_score:
            # delete the oligo, both dictionary and dataframe are passed as a reference
            del database_region[oligo_id]
            oligos_scores.drop(oligo_id, inplace=True)
    return database_region, oligos_scores, best_set
