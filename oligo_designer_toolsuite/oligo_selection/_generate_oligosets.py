from typing import Callable

import networkx as nx
import numpy as np
import pandas as pd
from joblib import Parallel, delayed

from ..database import OligoDatabase
from ..oligo_efficiency import OligoScoringBase, SetScoringBase


class OligosetGenerator:
    """This class is used to generate ranked, non-overlapping oligo sets. Each oligo is scored according to the given ``oligos_scoring`` class
    and the sets according to the given ``set_scoring`` calss.

    :param oligoset_size: optimal number of oligos that each set should contain
    :type oligoset_size: int
    :param min_oligoset_size: minimum number of oligos that each set should contain
    :type min_oligoset_size: int
    :param oligos_scoring: class that scores the oligos and the sets of oligos
    :type oligos_scoring: OligoScoring class
    :param dir_oligo_sets: directory where the sets are written
    :type dir_oligo_sets: str
    :param heurustic_selection: functions that preselects the oligos making the seach of the best oligosets less demanding
    :type dir_oligo_sets: Callable
    """

    def __init__(
        self,
        oligoset_size: int,
        min_oligoset_size: int,
        oligos_scoring: OligoScoringBase,
        set_scoring: SetScoringBase,
        heurustic_selection: Callable = None,
    ) -> None:
        """Initialize the class."""

        self.oligoset_size = oligoset_size
        self.min_oligoset_size = min_oligoset_size
        self.heurustic_selection = heurustic_selection
        self.oligos_scoring = oligos_scoring
        self.set_scoring = set_scoring

    def apply(
        self, oligo_database: OligoDatabase, n_sets: int = 50, n_jobs: int = None
    ):
        """Generates in parallel the oligosets and selects the best ``n_sets`` according to the
        The database class is updated, in particular form the ``oligos_DB`` are filtered out all the oligos that don't belong to any oligoset and in the class attruibute ``oligosets`` are stored
        the computed oligosets. The latter is a dictionary having as keys the regions names and as values a pandas.DataFrame containinig the oligosets. The strucutre of the pandas.DataFrame is the following:

        +-------------+----------+----------+----------+-------+----------+-------------+-------------+-------+
        | oligoset_id | oligo_0  | oligo_1  | oligo_2  |  ...  | oligo_n  | set_score_1 | set_score_2 |  ...  |
        +-------------+----------+----------+----------+-------+----------+-------------+-------------+-------+
        | 0           | AGRN_184 | AGRN_133 | AGRN_832 |  ...  | AGRN_706 | 0.3445      | 1.2332      |  ...  |
        +-------------+----------+----------+-----+----+-------+----------+-------------+-------------+-------+

        :param oligo_database: class containg the oligo sequences and their features
        :type oligo_database: OligoDatabase class
        :param n_sets: maximal number of sets that will be generated, defaults to 50
        :type n_sets: int, optional
        :param n_jobs: nr of cores used, if None the value set in database class is used, defaults to None
        :type n_jobs: int
        :return: Updated database class
        :rtype: OligoDatabase class
        """

        # set the number of cores
        if n_jobs is None:
            n_jobs = oligo_database.n_jobs

        regions = list(oligo_database.database.keys())
        # get the oligo set for this region in parallel
        database_regions = Parallel(
            n_jobs=n_jobs
        )(  # there should be an explicit return
            delayed(self._get_oligo_set_for_gene)(
                region, oligo_database.database[region], n_sets
            )
            for region in regions
        )
        # restore the database
        for region, database_region in zip(regions, database_regions):
            if database_region is None:  # if no sets have been found
                oligo_database.database[region] = {}  # oligoset is not generated
            else:
                oligo_database.oligosets[region] = database_region["oligosets"]
                del database_region["oligosets"]
                oligo_database.database[region] = database_region
        oligo_database.remove_regions_with_insufficient_oligos(
            pipeline_step="oligoset generation"
        )

        return oligo_database

    def _get_oligo_set_for_gene(self, region: str, database_region: dict, n_sets: int):
        """Generate the oligosets for a region.

        :param region: region for which the oligosets are computed
        :type region: str
        :param database_region: dictionary with the information about the oligos in the region
        :type database_region: dict
        :param n_sets: number of sets to generate
        :type n_sets: int
        :return: updated_oligos
        :rtype: dict
        """
        # create the overlapping matrix
        overlapping_matrix = self._get_overlapping_matrix(database_region)
        # create the set
        n, oligosets, database_region = self._get_non_overlapping_sets(
            database_region, overlapping_matrix, n_sets
        )
        # delete the useless variable to free some memory(overlapping matrix)
        del overlapping_matrix  # free some memory
        # write the set
        if oligosets is None:  # value passed as a parameter
            return None  # no more oligos are left for this region
        else:
            # update the dictionary adding the key sets for the oligos in a set and  delleting the oligos not included in any set
            updated_database_region = {}
            for set_id, row in oligosets.iterrows():
                for i in range(1, n + 1):
                    if row[i] not in updated_database_region:
                        updated_database_region[row[i]] = database_region[row[i]]
            # add also the oligosets
            updated_database_region["oligosets"] = oligosets
            return updated_database_region

    def _get_overlapping_matrix(self, database_region):
        """Creates a matrix that encodes the overlapping of different oligos. the matrix has dimensions n_oligos * n_oligos where each
        row and column belong to a oligo. Each entry contains 1 if the the correspondent oligos don't overlap and 0 if they overlap, this
        is done because in the next this matrix will be used as an adjacency matrix and the sets of non overlapping oligos are cliques of this graph.

        :pram database_region: dictionary containing all the oligos of the region
        :type database_region: dict
        :return: overlapping matrix
        :rtype: pandas.DataFrame
        """

        def _get_overlap(seq1_intervals, seq2_intervals):
            for a in seq1_intervals:
                for b in seq2_intervals:
                    if min(a[1], b[1]) - max(a[0], b[0]) > -1:
                        return True
            return False

        intervals = []
        oligos_indices = []
        for oligo_id in database_region.keys():
            oligos_indices.append(oligo_id)  # keep track of the indices
            interval = []
            for start, end in zip(
                database_region[oligo_id]["start"], database_region[oligo_id]["end"]
            ):
                interval.append(
                    [start, end]
                )  # save a list of couples of [start,end] of the duplicates of that oligo
            intervals.append(interval)

        overlapping_matrix = np.zeros(
            (len(intervals), len(intervals)), dtype=int
        )  # on the diagonal we have overlaps
        for i in range(len(intervals)):
            for j in range(i + 1, len(intervals)):
                if _get_overlap(intervals[i], intervals[j]):
                    overlapping_matrix[i, j] = 1
        overlapping_matrix = (
            overlapping_matrix
            + np.transpose(overlapping_matrix)
            + np.eye(len(intervals))
        )
        overlapping_matrix = (
            np.ones((len(intervals), len(intervals)), dtype=int) - overlapping_matrix
        )
        overlapping_matrix = pd.DataFrame(
            data=overlapping_matrix,
            columns=oligos_indices,
            index=oligos_indices,
            dtype=int,
        )

        return overlapping_matrix

    def _get_non_overlapping_sets(self, database_region, overlapping_matrix, n_sets):
        """Generates the non overlapping sets and return the best n_sets. Firstly the oligos are scored and then, if it is available,
        an heuristic method is used to reduce the number of oligos to the more promising ones. Then all the possible combination of non overlapping
        sets are considered and the best n-sets are returned.

        :param database_region: dictionary containing all the oligos of the region
        :type database_region: dict
        :param overlapping_matrix: matrix containig information about which oligos overlap
        :type overlapping_matrix: pandas.DataFrame
        :param n_sets: number of sets to generate
        :type n_sets: int
        :return: actual size of the sets, DataFrame containnig the computed oligosets, dictionary of oligos
        :rtype: int, pandas.DataFrame, dict
        """
        oligo_indices = np.array(overlapping_matrix.index)
        # Score oligos and create a pd series with the same order of oligos as in overlapping matrix
        database_region, oligos_scores = self.oligos_scoring.apply(
            database_region, oligo_indices
        )  # add a entry score to the oligos
        # Represent overlap matrix as graph
        G = nx.convert_matrix.from_numpy_array(overlapping_matrix.values)
        # First check if there are no cliques with n oligos
        cliques = nx.algorithms.clique.find_cliques(G)
        n = self.oligoset_size
        n_max = 0
        for clique in cliques:
            n_max = max(len(clique), n_max)
            if n_max >= n:
                break
        if n_max < n:
            if (
                n_max <= self.min_oligoset_size
            ):  # in this case we don't need to compute the sets
                return n_max, None, None
            else:
                n = n_max
        # if we have an heuristic apply it
        heuristic_oligoset = None
        if self.heurustic_selection is not None and n == self.oligoset_size:
            # apply the heuristic
            database_region, oligos_scores, heuristic_set = self.heurustic_selection(
                database_region, oligos_scores, overlapping_matrix, n
            )
            heuristic_oligoset = self.set_scoring.apply(
                heuristic_set, n
            )  # make it a list as for all the other cliques for future use
            # recompute the cliques
            G = nx.convert_matrix.from_numpy_array(
                overlapping_matrix.loc[oligos_scores.index, oligos_scores.index].values
            )
        # recompute the cliques
        cliques = nx.algorithms.clique.find_cliques(G)

        # Search the best set
        oligosets = []
        # Note: Search could be further optimised by iteratively throwing out oligos with worse scores then current best set
        for count, clique in enumerate(cliques):
            # Limit the number of combinations we iterate through
            if count > 100000:  # set a meaningful value
                break
            if len(clique) >= n:
                # Get oligo_ids of clique, maybe create a function
                clique_oligos = oligos_scores.iloc[clique]
                oligoset = self.set_scoring.apply(clique_oligos, n)
                oligosets.append(oligoset)
        # put the sets in a dataframe
        oligosets = pd.DataFrame(
            columns=[f"oligo_{i}" for i in range(n)]
            + [f"set_score_{i}" for i in range(len(oligosets[0]) - n)],
            data=oligosets,
        )
        # add the heurustuc best set, if is in the best n-sets then it will be kept
        if heuristic_oligoset:
            oligosets.loc[len(oligosets)] = heuristic_oligoset
        # Sort oligosets by score
        oligosets.drop_duplicates(inplace=True, subset=oligosets.columns[:-1])
        oligosets.sort_values(list(oligosets.columns[n:]), ascending=True, inplace=True)
        oligosets = oligosets.head(n_sets)
        oligosets.reset_index(drop=True, inplace=True)
        oligosets.insert(0, "oligoset_id", oligosets.index)

        return n, oligosets, database_region
