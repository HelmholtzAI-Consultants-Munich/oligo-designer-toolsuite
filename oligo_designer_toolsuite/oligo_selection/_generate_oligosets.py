############################################
# imports
############################################

import gc
import itertools
import random
from math import comb
from typing import Dict, List, Tuple

import pandas as pd
from joblib import Parallel, delayed
from joblib_progress import joblib_progress
from scipy.sparse import lil_matrix

from oligo_designer_toolsuite._constants import _TYPES_SEQ
from oligo_designer_toolsuite.database import OligoDatabase
from oligo_designer_toolsuite.oligo_efficiency_filter import (
    OligoScoringBase,
    SetScoringBase,
)

from ._heuristic_selection_methods import OligoSelectionPolicy

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
    :param max_oligos: Maximum number of oligos to consider in the generation process, defaults to 5000.
    :type max_oligos: int
    """

    def __init__(
        self,
        opt_oligoset_size: int,
        min_oligoset_size: int,
        oligos_scoring: OligoScoringBase,
        set_scoring: SetScoringBase,
        heuristic_selection: OligoSelectionPolicy,
        distance_between_oligos: int = 0,
        max_oligos: int = 5000,
    ) -> None:
        """Constructor for the OligosetGenerator class."""

        self.opt_oligoset_size = opt_oligoset_size
        self.min_oligoset_size = min_oligoset_size
        self.heuristic_selection = heuristic_selection
        self.oligos_scoring = oligos_scoring
        self.set_scoring = set_scoring
        self.distance_between_oligos = distance_between_oligos
        self.ascending = set_scoring.ascending
        self.max_oligos = max_oligos

    def apply(
        self,
        oligo_database: OligoDatabase,
        sequence_type: _TYPES_SEQ,
        pre_filter: bool,
        n_attempts: int = 10000,
        n_jobs: int = 1,
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
        :param n_attempts: Number of attempts to generate oligo sets, defaults to 10000.
        :type n_attempts: int
        :param n_jobs: Number of parallel jobs to use for processing, defaults to 1.
        :type n_jobs: int
        :return: The updated oligo database with selected oligo sets.
        :rtype: OligoDatabase
        """
        region_ids = list(oligo_database.database.keys())
        # get the oligo set for this region in parallel
        with joblib_progress(description="Find Oligosets", total=len(region_ids)):
            Parallel(
                n_jobs=n_jobs, prefer="threads", require="sharedmem"
            )(  # there should be an explicit return
                delayed(self._get_oligo_set_for_gene)(
                    oligo_database, region_id, sequence_type, pre_filter, n_attempts
                )
                for region_id in region_ids
            )

        oligo_database.remove_regions_with_insufficient_oligos(pipeline_step="oligoset generation")
        return oligo_database

    def _get_oligo_set_for_gene(
        self,
        oligo_database: OligoDatabase,
        region_id: str,
        sequence_type: _TYPES_SEQ,
        pre_filter: bool,
        n_attempts: int,
    ):
        """Processes a single gene region from the oligo database to generate non-overlapping sets of
        oligos based on scoring and set selection criteria.

        :param oligo_database: The database of oligos to process.
        :type oligo_database: OligoDatabase
        :param region_id: The ID of the region to process.
        :type region_id: str
        :param sequence_type: The type of sequences being used, must match one of the predefined sequence types.
        :type sequence_type: _TYPES_SEQ
        :param n_attempts: Number of attempts to generate oligo sets.
        :type n_attempts: int
        :return: None
        :rtype: None
        """

        # Score oligos and create a pd series
        oligo_database, oligos_scores = self.oligos_scoring.apply(
            oligo_database=oligo_database,
            region_id=region_id,
            sequence_type=sequence_type,
        )

        # sort oligos by score
        oligos_scores.sort_values(ascending=self.ascending, inplace=True)

        # hard limit on the number of oligos
        if len(oligos_scores) > self.max_oligos:
            # select the best oligos
            for oligo_id in oligos_scores.index[self.max_oligos :]:
                del oligo_database.database[region_id][oligo_id]
                oligos_scores.drop(oligo_id, inplace=True)

        # create the overlapping matrix
        overlapping_matrix = self._get_overlapping_matrix(oligo_database=oligo_database, region_id=region_id)

        # create the set
        oligosets = self.heuristic_selection.apply(
            oligos_scores=oligos_scores,
            overlapping_matrix=overlapping_matrix,
            pre_filter=pre_filter,
            n_attempts=n_attempts,
        )

        # Remove all oligos from database that are not part of oligosets
        oligos_keep = set()
        if oligosets is not None:
            oligosets_oligo_columns = [col for col in oligosets.columns if col.startswith("oligo_")]
            oligos_keep = set(oligosets[oligosets_oligo_columns].to_numpy().flatten())

        oligo_ids = list(oligo_database.database[region_id].keys())
        for oligo_id in oligo_ids:
            if not oligo_id in oligos_keep:
                del oligo_database.database[region_id][oligo_id]

        oligo_database.oligosets[region_id] = oligosets

        # delete the useless variable to free some memory(overlapping matrix)
        del overlapping_matrix  # free some memory
        gc.collect()

    def _get_overlapping_matrix(self, oligo_database: OligoDatabase, region_id: str):
        """
        Generates a matrix indicating the overlap between oligos in a given database region.
        The matrix is computed based on the start and end intervals of each oligo,
        with an optional parameter to adjust for permissible overlap distance.
        The matrix has dimensions n_oligos * n_oligos. Each entry contains 1 if the
        correspondent oligos don't overlap and 0 if they overlap. This matrix will be used
        as an adjacency matrix, and the sets of non-overlapping oligos are cliques of this graph.


        :param oligo_database: The database of oligos to process.
        :type oligo_database: OligoDatabase
        :return: A DataFrame where each cell [i, j] is 0 if oligos i and j overlap and 1 otherwise, with oligo IDs as indices and columns.
        :rtype: pd.DataFrame
        """

        def _get_overlap(seq1_intervals, seq2_intervals):
            # Determine if two ligos overlap based on a distance value
            return any(
                min(a[1], b[1]) - max(a[0], b[0]) >= -self.distance_between_oligos
                for a in seq1_intervals
                for b in seq2_intervals
            )

        # Keep track of the indices
        overlapping_matrix_ids = list(oligo_database.database[region_id].keys())

        # Get all intervals (start, end)
        intervals = [
            [
                [start[0], end[0]]
                for start, end in zip(
                    oligo_database.database[region_id][oligo_id]["start"],
                    oligo_database.database[region_id][oligo_id]["end"],
                )
            ]
            for oligo_id in overlapping_matrix_ids
        ]

        # Create a sparse overlap matrix
        n_oligos = len(overlapping_matrix_ids)
        overlapping_matrix = lil_matrix((n_oligos, n_oligos), dtype=int)

        # Calculate only upper triangle matrix since the matrix is symmetric
        for i in range(n_oligos):
            for j in range(i + 1, n_oligos):
                if _get_overlap(intervals[i], intervals[j]):
                    overlapping_matrix[i, j] = 1

        # Fill values of lower triangle
        overlapping_matrix = overlapping_matrix.maximum(overlapping_matrix.transpose())
        # Set diagonal elements to 1 as oligos always overlap with themselves
        overlapping_matrix.setdiag(1)

        # Create a sparse matrix containing only ones
        ones_matrix = lil_matrix((n_oligos, n_oligos), dtype=int)
        ones_matrix[:, :] = 1

        # Invert theoverlap matrix by subtracting the overlapping matrix from the ones matrix
        overlapping_matrix = ones_matrix - overlapping_matrix
        overlapping_matrix = overlapping_matrix.tocsr()

        # overlapping_matrix = pd.DataFrame(
        #    data=overlapping_matrix.toarray(),
        #    columns=oligos_indices,
        #    index=oligos_indices,
        #    dtype=int,
        # )

        return overlapping_matrix


class HomogeneousPropertyOligoSetGenerator:
    """
    A generator class for creating oligo sets that ensure homogeneity in specified properties.

    This class generates sets of oligos by selecting combinations with the lowest weighted sum
    of variances for specified properties, which ensures homogeneity within each set.

    :param set_size: The size of each oligo set.
    :type set_size: int
    :param properties: A dictionary where the key is the property name and the value is the weight for that property.
    :type properties: dict
    """

    def __init__(self, set_size: int, properties: Dict[str, float]) -> None:
        """Constructor for the HomogeneousPropertyOligoSetGenerator class."""
        self.set_size = set_size
        self.properties = properties

    def apply(
        self, oligo_database: OligoDatabase, n_attempts: int = 1, n_combinations: int = 1000, n_jobs: int = 1
    ) -> OligoDatabase:
        """
        Applies the oligo set generation process to an entire oligo database and returns an updated database with selected best `n_attempts` oligo sets.
        Oligosets are stored in the class attribute `oligosets`, which is a dictionary with region names as keys and oligoset dataframes as values.
        The structure of the pandas.DataFrame is the following:


        +-------------+----------+----------+----------+-------+----------+-------------+
        | oligoset_id | oligo_0  | oligo_1  | oligo_2  |  ...  | oligo_n  | set_score   |
        +-------------+----------+----------+----------+-------+----------+-------------+
        | 0           | AGRN_184 | AGRN_133 | AGRN_832 |  ...  | AGRN_706 | 0.3445      |
        +-------------+----------+----------+-----+----+-------+----------+-------------+

        :param oligo_database: The oligo database to generate sets from.
        :type oligo_database: OligoDatabase
        :param n_attempts: Number of sets to generate for each region, defaults to 1.
        :type n_attempts: int, optional
        :param n_combinations: Number of random combinations to generate, defaults to 1000.
        :type n_combinations: int, optional
        :param n_jobs: Number of parallel jobs to run, defaults to 1.
        :type n_jobs: int, optional
        :return: The updated oligo database with generated oligo sets.
        :rtype: OligoDatabase
        """

        region_ids = list(oligo_database.database.keys())
        with joblib_progress(description="Find Oligosets", total=len(region_ids)):
            Parallel(n_jobs=n_jobs, prefer="threads", require="sharedmem")(
                delayed(self._get_oligo_sets_for_region)(oligo_database, region_id, n_attempts, n_combinations)
                for region_id in region_ids
            )

        oligo_database.remove_regions_with_insufficient_oligos(pipeline_step="oligoset generation")
        return oligo_database

    def _get_oligo_sets_for_region(
        self, oligo_database: OligoDatabase, region_id: str, n_attempts: int, n_combinations: int
    ) -> None:
        """
        Generate oligo sets for a specific region.

        This method generates all possible combinations of oligos for a given region and scores them based on
        the specified properties. The top N sets with the lowest weighted sum of variances are selected.

        :param oligo_database: The oligo database to generate sets from.
        :type oligo_database: OligoDatabase
        :param region_id: The ID of the region to generate oligo sets for.
        :type region_id: str
        :param n_attempts: Number of sets to generate for the region.
        :type n_attempts: int
        :param n_combinations: Number of random combinations to generate.
        :type n_combinations: int
        """

        region_dict = oligo_database.database[region_id]
        oligo_df = pd.DataFrame.from_dict(region_dict, orient="index")

        # # check if all properties in self.properties are in oligo_df columns
        for property in self.properties:
            if property not in oligo_df.columns:
                raise ValueError(
                    f"Property '{property}' is not present in oligo database please calculate it first using oligo_designer_toolsuite.OligoAttributes()."
                )

        combinations = self._generate_random_combinations(oligo_df.index, self.set_size, n_combinations)

        scored_combinations = [
            self._score_combination(oligo_df, list(combination)) for combination in combinations
        ]
        sorted_combinations = sorted(scored_combinations, key=lambda x: x[1], reverse=False)
        best_combinations = [combination for combination in sorted_combinations[:n_attempts]]

        rows = [[idx] + oligos + [score] for idx, (oligos, score) in enumerate(best_combinations)]
        columns = ["oligoset_id"] + [f"oligo_{i}" for i in range(self.set_size)] + ["set_score"]

        oligo_database.oligosets[region_id] = pd.DataFrame(rows, columns=columns)

    def _score_combination(self, oligo_df: pd.DataFrame, combination: List[str]) -> Tuple[List[str], float]:
        """
        Score a combination of oligos based on the specified properties.

        This method calculates the score for a combination of oligos by computing the weighted sum of the variances
        of the specified properties. The lower the score, the more homogeneous the set is with respect to the
        specified properties.

        :param oligo_df: The DataFrame containing oligo information.
        :type oligo_df: pd.DataFrame
        :param combination: A list of oligo IDs representing a combination.
        :type combination: list
        :return: A tuple containing the combination and its score.
        :rtype: tuple
        """

        oligo_set = oligo_df.loc[combination]
        score = sum([oligo_set[property].var() * self.properties[property] for property in self.properties])
        return combination, score

    @staticmethod
    def _generate_random_combinations(arr, combination_size, number_of_combinations):
        total_combinations = comb(len(arr), combination_size)

        if total_combinations <= number_of_combinations:
            return list(itertools.combinations(arr, combination_size))

        seen_combinations = set()
        while len(seen_combinations) < number_of_combinations:
            combination = tuple(sorted(random.sample(list(arr), combination_size)))
            if combination not in seen_combinations:
                seen_combinations.add(combination)
        return list(seen_combinations)
