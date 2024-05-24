############################################
# imports
############################################

import gc
from typing import Callable

import networkx as nx
import pandas as pd
from joblib import Parallel, delayed
from scipy.sparse import lil_matrix

from oligo_designer_toolsuite._constants import _TYPES_SEQ
from oligo_designer_toolsuite.database import OligoDatabase
from oligo_designer_toolsuite.oligo_efficiency_filter import (
    OligoScoringBase,
    SetScoringBase,
)

############################################
# Oligo set Generation Classes
############################################


class OligosetGeneratorIndependentSet:
    """
    Generator class for optimal sets of non-overlapping oligonucleotides, designed to manage and execute the scoring
    and selection of oligonucleotide sets based on given `oligo_scoring` and `set_scoring` classes.

    :param opt_oligoset_size: Optimal size of the oligo set to be generated.
    :type opt_oligoset_size: int
    :param min_oligoset_size: Minimum size of the oligo set considered for generation.
    :type min_oligoset_size: int
    :param oligos_scoring: An instance of OligoScoringBase used for scoring individual oligos.
    :type oligos_scoring: OligoScoringBase
    :param set_scoring: An instance of SetScoringBase used for scoring sets of oligos.
    :type set_scoring: SetScoringBase
    :param heurustic_selection: A callable for heuristic selection of oligo sets, default is None.
    :type heurustic_selection: Callable, optional
    :param distance_between_oligos: Distance between neighboring oligos, e.g. -x: oligos overlap x bases; 0: oligos can be next to each other; +x: oligos are x bases apart
    :param ascending: Determines if the scoring should sort in ascending order dependent on meaning of scores, defaults to True.
    :type ascending: bool
    :param max_oligos: Maximum number of oligos to consider in the generation process, defaults to 5000.
    :type max_oligos: int
    """

    def __init__(
        self,
        opt_oligoset_size: int,
        min_oligoset_size: int,
        oligos_scoring: OligoScoringBase,
        set_scoring: SetScoringBase,
        heuristic_selection: Callable = None,
        distance_between_oligos: int = 0,
        ascending: bool = True,
        max_oligos: int = 5000,
    ) -> None:
        """Constructor for the OligosetGenerator class."""

        self.opt_oligoset_size = opt_oligoset_size
        self.min_oligoset_size = min_oligoset_size
        self.heuristic_selection = heuristic_selection
        self.oligos_scoring = oligos_scoring
        self.set_scoring = set_scoring
        self.distance_between_oligos = distance_between_oligos
        self.ascending = ascending
        self.max_oligos = max_oligos

    def apply(
        self, oligo_database: OligoDatabase, sequence_type: _TYPES_SEQ, n_sets: int = 50, n_jobs: int = 1
    ):
        """
        Applies the oligo set generation process to an entire oligo database and returns updated database with selected best `n_sets` oligo sets.
        Oligosets are stored in the class attribute `oligosets`, which is a dictionary with region names as keys and oligoset dataframes as values.
        The structure of the pandas.DataFrame is the following:


        +-------------+----------+----------+----------+-------+----------+-------------+-------------+-------+
        | oligoset_id | oligo_0  | oligo_1  | oligo_2  |  ...  | oligo_n  | set_score_1 | set_score_2 |  ...  |
        +-------------+----------+----------+----------+-------+----------+-------------+-------------+-------+
        | 0           | AGRN_184 | AGRN_133 | AGRN_832 |  ...  | AGRN_706 | 0.3445      | 1.2332      |  ...  |
        +-------------+----------+----------+-----+----+-------+----------+-------------+-------------+-------+


        :param oligo_database: The database of oligos to process.
        :type oligo_database: OligoDatabase
        :param sequence_type: The type of sequence to be used in the oligos scoring.
        :type sequence_type: _TYPES_SEQ
        :param n_sets: Number of oligo sets to generate for each region, defaults to 50.
        :type n_sets: int
        :param n_jobs: Number of parallel jobs to use for processing, defaults to 1.
        :type n_jobs: int
        :return: The updated oligo database with selected oligo sets.
        :rtype: OligoDatabase
        """
        regions = list(oligo_database.database.keys())
        # get the oligo set for this region in parallel
        database_regions = Parallel(n_jobs=n_jobs)(  # there should be an explicit return
            delayed(self._get_oligo_set_for_gene)(oligo_database.database[region], sequence_type, n_sets)
            for region in regions
        )

        # restore the database
        for region, (database_region, oligosets_region) in zip(regions, database_regions):
            if database_region is None:  # if no sets have been found
                oligo_database.database[region] = {}  # oligoset is not generated
            else:
                oligo_database.oligosets[region] = oligosets_region
                oligo_database.database[region] = database_region

        oligo_database.remove_regions_with_insufficient_oligos(pipeline_step="oligoset generation")
        return oligo_database

    def _get_oligo_set_for_gene(self, database_region: dict, sequence_type: _TYPES_SEQ, n_sets: int):
        """Processes a single gene region from the oligo database to generate non-overlapping sets of
        oligos based on scoring and set selection criteria.

        :param database_region: A dictionary representing a specific gene region containing oligos and their attributes.
        :type database_region: dict
        :param sequence_type: The type of sequences being used, must match one of the predefined sequence types.
        :type sequence_type: _TYPES_SEQ
        :param n_sets: The number of oligo sets to attempt to generate for the region.
        :type n_sets: int
        :return: A tuple containing the updated database region with selected oligos and a DataFrame of the oligo sets.
        :rtype: tuple
        """

        # Score oligos and create a pd series
        database_region, oligos_scores = self.oligos_scoring.apply(
            oligos=database_region, sequence_type=sequence_type
        )

        # add a entry score to the oligos
        oligos_scores.sort_values(ascending=self.ascending, inplace=True)

        # hard limit on the number of oligos
        if len(oligos_scores) > self.max_oligos:
            # select the best oligos
            for oligo_id in oligos_scores.index[self.max_oligos :]:
                del database_region[oligo_id]
                oligos_scores.drop(oligo_id, inplace=True)

        # create the overlapping matrix
        overlapping_matrix = self._get_overlapping_matrix(database_region)

        # create the set
        n, oligosets, database_region = self._get_non_overlapping_sets(
            database_region, overlapping_matrix, oligos_scores, n_sets
        )

        # delete the useless variable to free some memory(overlapping matrix)
        del overlapping_matrix  # free some memory
        gc.collect()

        # write the set
        if oligosets is None:  # value passed as a parameter
            return None, None  # no more oligos are left for this region
        else:
            # update the dictionary adding the key sets for the oligos in a set and  delleting the oligos not included in any set
            updated_database_region = {}
            for _, row in oligosets.iterrows():
                for i in range(1, n + 1):
                    if row.iloc[i] not in updated_database_region:
                        updated_database_region[row.iloc[i]] = database_region[row.iloc[i]]
            return updated_database_region, oligosets

    def _get_overlapping_matrix(self, database_region: dict):
        """
        Generates a matrix indicating the overlap between oligos in a given database region.
        The matrix is computed based on the start and end intervals of each oligo,
        with an optional parameter to adjust for permissible overlap distance.
        The matrix has dimensions n_oligos * n_oligos. Each entry contains 1 if the
        correspondent oligos don't overlap and 0 if they overlap. This matrix will be used
        as an adjacency matrix, and the sets of non-overlapping oligos are cliques of this graph.


        :param database_region: Dictionary containing oligo IDs and their corresponding start and end intervals.
        :type database_region: dict
        :return: A DataFrame where each cell [i, j] is 0 if oligos i and j overlap and 1 otherwise, with oligo IDs as indices and columns.
        :rtype: pd.DataFrame
        """

        def _get_overlap(seq1_intervals, seq2_intervals):
            for a in seq1_intervals:
                for b in seq2_intervals:
                    if min(a[1], b[1]) - max(a[0], b[0]) >= -self.distance_between_oligos:
                        return True
            return False

        oligos_indices = list(database_region.keys())  # Keep track of the indices
        intervals = [
            [
                [start[0], end[0]]
                for start, end in zip(database_region[oligo_id]["start"], database_region[oligo_id]["end"])
            ]
            for oligo_id in oligos_indices
        ]

        # Create a sparse matrix
        n_oligos = len(intervals)
        overlapping_matrix = lil_matrix((n_oligos, n_oligos), dtype=int)

        for i in range(n_oligos):
            for j in range(i + 1, n_oligos):
                if _get_overlap(intervals[i], intervals[j]):
                    overlapping_matrix[i, j] = 1

        overlapping_matrix = overlapping_matrix.maximum(overlapping_matrix.transpose())
        overlapping_matrix.setdiag(1)  # Set diagonal elements to 1

        # Create a sparse matrix containing only ones
        ones_matrix = lil_matrix((n_oligos, n_oligos), dtype=int)
        ones_matrix[:, :] = 1

        # Invert the matrix by subtracting the overlapping matrix from the ones matrix
        overlapping_matrix = ones_matrix - overlapping_matrix

        overlapping_matrix = pd.DataFrame(
            data=overlapping_matrix.toarray(),
            columns=oligos_indices,
            index=oligos_indices,
            dtype=int,
        )

        return overlapping_matrix

    def _get_non_overlapping_sets(
        self, database_region: dict, overlapping_matrix: pd.DataFrame, oligos_scores: pd.Series, n_sets: int
    ):
        """
        Generates a list of non-overlapping oligo sets from the given region based on the specified criteria.
        This function finds cliques (non-overlapping sets) in an overlap graph of oligos and scores these sets
        to determine the optimal ones. If the "heuristic_selection" parameter is set, this heuristic is used to find
        non-overlapping oligo sets.

        :param database_region: Contains oligos data for a specific region.
        :type database_region: dict
        :param overlapping_matrix: A DataFrame indicating overlap between oligos.
        :type overlapping_matrix: pd.DataFrame
        :param oligos_scores: Scores for each oligo used to rank them within their sets.
        :type oligos_scores: pd.Series
        :param n_sets: Number of sets to return.
        :type n_sets: int
        :return: A tuple containing the number of oligos per set, DataFrame of oligo sets, and updated database region.
        :rtype: tuple
        """

        # Represent overlap matrix as graph
        G = nx.convert_matrix.from_numpy_array(overlapping_matrix.values)
        G = nx.relabel_nodes(G, {i: overlapping_matrix.index[i] for i in range(len(oligos_scores.index))})

        # First check if there are no cliques with n oligos
        cliques = nx.algorithms.clique.find_cliques(G)
        n = self.opt_oligoset_size
        n_max = 0

        for clique in cliques:
            n_max = max(len(clique), n_max)
            if n_max >= n:
                break

        if n_max < n:
            if n_max <= self.min_oligoset_size:  # in this case we don't need to compute the sets
                return n_max, None, None
            else:
                n = n_max

        # if we have an heuristic apply it
        heuristic_oligoset = None
        if self.heuristic_selection is not None and n == self.opt_oligoset_size:
            # apply the heuristic
            database_region, oligos_scores, heuristic_set = self.heuristic_selection(
                database_region, oligos_scores, overlapping_matrix, n, self.ascending
            )
            heuristic_oligoset = self.set_scoring.apply(
                heuristic_set, n
            )  # make it a list as for all the other cliques for future use
            # recompute the cliques
            overlapping_matrix = overlapping_matrix.loc[oligos_scores.index, oligos_scores.index]
            G = nx.convert_matrix.from_numpy_array(overlapping_matrix.values)
            G = nx.relabel_nodes(
                G,
                {i: overlapping_matrix.index[i] for i in range(len(oligos_scores.index))},
            )

        # recompute the cliques
        cliques = nx.algorithms.clique.find_cliques(G)
        # Search the best set

        if heuristic_oligoset:
            oligosets = [
                heuristic_oligoset
            ]  # add the heurustuc best set, if is in the best n-sets then it will be kept
        else:
            oligosets = []  # initialize the list of sets

        # Note: Search could be further optimised by iteratively throwing out oligos with worse scores then current best set
        for count, clique in enumerate(cliques):
            # Limit the number of combinations we iterate through
            if count > 100000:  # set a meaningful value
                break
            if len(clique) >= n:
                # Get oligo_ids of clique, maybe create a function
                clique_oligos = oligos_scores.loc[clique]
                oligoset = self.set_scoring.apply(clique_oligos, n)
                oligosets.append(oligoset)

        # put the sets in a dataframe
        if len(oligosets) > 0:
            oligosets = pd.DataFrame(
                columns=[f"oligo_{i}" for i in range(n)]
                + [f"set_score_{i}" for i in range(len(oligosets[0]) - n)],
                data=oligosets,
            )
        # Sort oligosets by score
        oligosets.drop_duplicates(inplace=True, subset=oligosets.columns[:-1])
        oligosets.sort_values(list(oligosets.columns[n:]), ascending=self.ascending, inplace=True)
        oligosets = oligosets.head(n_sets)
        oligosets.reset_index(drop=True, inplace=True)
        oligosets.insert(0, "oligoset_id", oligosets.index)

        return n, oligosets, database_region
