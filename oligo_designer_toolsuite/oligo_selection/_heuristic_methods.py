import numpy as np


def padlock_heuristic_selection(
    probes, probes_scores, overlapping_matrix, n_probe, n_trials=10000
):
    """This method tries to find empirically a good set of teh padlock design of probes. For teh best n_trials probes a set is cretaed
    by adding one by one the best non overlapping probe untiil a set is created. The best obtained set saved and returned. Moreover, foor how this particular
    scoring method works we know that all the probes with a higher score than the worse probe of the best set are going to yield a worse set the one just found,
    hence are deleted form the dictionary and the series of scores.

    :param probes: dictionary with all the probes
    :type probes: dict
    :param probes_scores: scores of the probes
    :type probes_scores: pandas.Series
    :param overlapping_matrix: matrix containig information about the overlapping of the probes
    :type overlapping_matrix: pandas.DataFrame
    :param n_probe: size of the set
    :type n_probe: int
    :param n_trials: number of sets to be tried, defaults to 10000
    :type n_trials: int, optional
    :return: filtered dcit of the probes, filtered scores fo teh probes, bes set found
    :rtype: dict, pandas.Series, pandas.Series
    """

    # sort the probes by their score
    probes_sorted = probes_scores.sort_values()
    probe_ids_sorted = probes_sorted.index.tolist()
    # overlapping matrix must have a consistent oreder
    mat_sorted = overlapping_matrix.loc[probe_ids_sorted, probe_ids_sorted].values

    max_score = (
        probes_sorted.iloc[-1] * 1.1
    )  # already sorted df, the max is the last entry
    best_idx_set = []
    for first_idx in range(
        min(len(probe_ids_sorted), n_trials)
    ):  # use the integer index because the matric is a np array
        set_idxs = np.array([first_idx])
        for _ in range(n_probe - 1):
            # find first probe in sorted array that is not overlapping with any selected probe
            no_overlap = np.all(mat_sorted[set_idxs], axis=0)
            if np.any(no_overlap):
                set_idxs = np.append(set_idxs, np.where(no_overlap)[0][0])
            else:
                break
        if len(set_idxs) == n_probe:
            score = np.max(probes_sorted.values[set_idxs])
            if score < max_score:
                max_score = score
                best_idx_set = set_idxs
    best_set = probes_sorted.iloc[best_idx_set]  # adapt with the

    for probe_id in probes_scores.index:
        if probes_scores[probe_id] > max_score:
            # delete the probe, both dictionary and dataframe are passed as a reference
            del probes[probe_id]
            probes_scores.drop(probe_id, inplace=True)
    return probes, probes_scores, best_set
