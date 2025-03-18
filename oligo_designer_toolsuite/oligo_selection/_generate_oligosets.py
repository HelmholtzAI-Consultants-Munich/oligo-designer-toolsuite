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
from scipy.sparse import csr_matrix, lil_matrix

from oligo_designer_toolsuite._constants import _TYPES_SEQ
from oligo_designer_toolsuite.database import OligoDatabase
from oligo_designer_toolsuite.oligo_efficiency_filter import (
    OligoScoringBase,
    SetScoringBase,
)

from ._selection_methods import OligoSelectionPolicy

############################################
# Oligo set Generation Classes
############################################


class OligosetGeneratorIndependentSet:
    """
    A class to generate multiple non-overlapping sets of oligos based on a selection policy and scoring methods.
    The class optimizes the sets of oligos according to a scoring function for individual oligos and the entire set,
    while also respecting constraints such as minimum and maximum set size and the distance between oligos.

    :param selection_policy: Policy used to select oligos for the set based on predefined criteria.
    :type selection_policy: OligoSelectionPolicy
    :param oligos_scoring: Scoring function used to evaluate individual oligos in the selection process.
    :type oligos_scoring: OligoScoringBase
    :param set_scoring: Scoring function used to evaluate the overall quality of each oligo set.
    :type set_scoring: SetScoringBase
    :param max_oligos: Maximum number of oligos to include in the set optimizatoin process. If None, there is no limit on the number of oligos.
    :type max_oligos: int, optional
    :param distance_between_oligos: Minimum allowed distance between oligos in the set to avoid clustering of oligos.
    :type distance_between_oligos: int, optional
    """

    def __init__(
        self,
        selection_policy: OligoSelectionPolicy,
        oligos_scoring: OligoScoringBase,
        set_scoring: SetScoringBase,
        max_oligos: int = None,
        distance_between_oligos: int = 0,
    ) -> None:
        """Constructor for the OligosetGeneratorIndependentSet class."""
        self.selection_policy = selection_policy
        self.oligos_scoring = oligos_scoring
        self.set_scoring = set_scoring
        self.max_oligos = max_oligos
        self.distance_between_oligos = distance_between_oligos

    def apply(
        self,
        oligo_database: OligoDatabase,
        sequence_type: _TYPES_SEQ,
        set_size_opt: int,
        set_size_min: int,
        n_sets: int = 1,
        n_jobs: int = 1,
    ) -> OligoDatabase:
        """
        Applies the oligo set generation process to an oligo database, finding the best oligos for each region and
        specified sequence type based on the selection policy and scoring schemes. The process is parallelized for efficiency.

        Oligosets are stored in the class attribute `oligosets`, which is a dictionary with region names as keys and oligoset
        dataframes as values. The structure of the pandas.DataFrame is the following:

        +-------------+----------+----------+----------+-------+----------+-------------+-------------+-------+
        | oligoset_id | oligo_0  | oligo_1  | oligo_2  |  ...  | oligo_n  | set_score_1 | set_score_2 |  ...  |
        +-------------+----------+----------+----------+-------+----------+-------------+-------------+-------+
        | 0           | AGRN_184 | AGRN_133 | AGRN_832 |  ...  | AGRN_706 | 0.3445      | 1.2332      |  ...  |
        +-------------+----------+----------+-----+----+-------+----------+-------------+-------------+-------+


        :param oligo_database: The OligoDatabase containing the oligonucleotides and their associated attributes.
        :type oligo_database: OligoDatabase
        :param sequence_type: The type of sequence to be used for the set calculations.
        :type sequence_type: _TYPES_SEQ["oligo", "target"]
        :param set_size_opt: The optimal size of each oligo set.
        :type set_size_opt: int
        :param set_size_min: The minimum allowed size of each oligo set.
        :type set_size_min: int
        :param n_sets: The number of oligo sets to generate.
        :type n_sets: int
        :param n_jobs: The number of parallel jobs to use for processing.
        :type n_jobs: int, optional
        :return: The updated oligo database with the generated oligo sets.
        :rtype: OligoDatabase
        """
        region_ids = list(oligo_database.database.keys())
        # get the oligo set for this region in parallel
        with joblib_progress(description="Find Oligosets", total=len(region_ids)):
            Parallel(
                n_jobs=n_jobs, prefer="threads", require="sharedmem"
            )(  # there should be an explicit return
                delayed(self._get_oligo_set_for_gene)(
                    oligo_database, region_id, sequence_type, set_size_opt, set_size_min, n_sets
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
        set_size_opt: int,
        set_size_min: int,
        n_sets: int,
    ) -> None:
        """
        Computes the oligo set for a specific gene region by scoring, filtering, and selecting oligos.
        This includes generating a proximity matrix and applying a selection policy to create the optimal oligo set.

        :param oligo_database: The OligoDatabase containing the oligonucleotides and their associated attributes.
        :type oligo_database: OligoDatabase
        :param region_id: Region ID to process.
        :type region_id: str
        :param sequence_type: The type of sequence to be used for the set calculations.
        :type sequence_type: _TYPES_SEQ["oligo", "target"]
        :param set_size_opt: The optimal size of each oligo set.
        :type set_size_opt: int
        :param set_size_min: The minimum allowed size of each oligo set.
        :type set_size_min: int
        :param n_sets: The number of oligo sets to generate.
        :type n_sets: int
        """
        # Score oligos and create a pd series
        oligo_database, oligos_scores = self.oligos_scoring.apply(
            oligo_database=oligo_database,
            region_id=region_id,
            sequence_type=sequence_type,
        )

        # sort oligos by score
        oligos_scores.sort_values(ascending=self.set_scoring.ascending, inplace=True)

        # hard limit on the number of oligos
        if self.max_oligos and len(oligos_scores) > self.max_oligos:
            # select the best oligos
            for oligo_id in oligos_scores.index[self.max_oligos :]:
                del oligo_database.database[region_id][oligo_id]
                oligos_scores.drop(oligo_id, inplace=True)

        # create the overlapping matrix
        non_overlap_matrix, non_overlap_matrix_ids = self._get_non_overlap_matrix(
            oligo_database=oligo_database, region_id=region_id
        )

        # create the set
        oligosets = self.selection_policy.apply(
            oligos_scores=oligos_scores,
            non_overlap_matrix=non_overlap_matrix,
            non_overlap_matrix_ids=non_overlap_matrix_ids,
            set_size_opt=set_size_opt,
            set_size_min=set_size_min,
            n_sets=n_sets,
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
        del non_overlap_matrix  # free some memory
        gc.collect()

    def _get_non_overlap_matrix(self, oligo_database: OligoDatabase, region_id: str) -> csr_matrix:
        """
        Generates a sparse matrix that represents the overlap between oligos in the specified region of the oligo database.
        The matrix is computed based on the intervals (start, end) of each oligo, with a distance threshold to determine overlap.
        The matrix has dimensions n_oligos * n_oligos. Each entry contains 1 if the correspondent oligos don't overlap and 0 if they overlap.

        :param oligo_database: The OligoDatabase containing the oligonucleotides and their associated attributes.
        :type oligo_database: OligoDatabase
        :param region_id: Region ID to process.
        :type region_id: str
        :return: A sparse matrix representing the non-overlap between oligos, and a list of oligo IDs.
        :rtype: tuple(csr_matrix, list)
        """

        def _get_overlap(seq1_intervals, seq2_intervals):
            # Determine if two ligos overlap based on a distance value
            return any(
                min(a[1], b[1]) - max(a[0], b[0]) >= -self.distance_between_oligos
                for a in seq1_intervals
                for b in seq2_intervals
            )

        # Keep track of the indices
        non_overlap_matrix_ids = list(oligo_database.database[region_id].keys())

        # Get all intervals (start, end)
        intervals = [
            [
                [start[0], end[0]]
                for start, end in zip(
                    oligo_database.database[region_id][oligo_id]["start"],
                    oligo_database.database[region_id][oligo_id]["end"],
                )
            ]
            for oligo_id in non_overlap_matrix_ids
        ]

        # Create a sparse overlap matrix
        n_oligos = len(non_overlap_matrix_ids)
        non_overlap_matrix = lil_matrix((n_oligos, n_oligos), dtype=int)

        # Calculate only upper triangle matrix since the matrix is symmetric
        for i in range(n_oligos):
            for j in range(i + 1, n_oligos):
                if _get_overlap(intervals[i], intervals[j]):
                    non_overlap_matrix[i, j] = 1

        # Fill values of lower triangle
        non_overlap_matrix = non_overlap_matrix.maximum(non_overlap_matrix.transpose())
        # Set diagonal elements to 1 as oligos always overlap with themselves
        non_overlap_matrix.setdiag(1)

        # Create a sparse matrix containing only ones
        ones_matrix = lil_matrix((n_oligos, n_oligos), dtype=int)
        ones_matrix[:, :] = 1

        # Invert theoverlap matrix by subtracting the overlapping matrix from the ones matrix
        non_overlap_matrix = ones_matrix - non_overlap_matrix
        non_overlap_matrix = non_overlap_matrix.tocsr()

        return non_overlap_matrix, non_overlap_matrix_ids


class HomogeneousPropertyOligoSetGenerator:
    """
    Generates oligo sets based on the homogeneity of specified oligo properties. The oligo sets are
    created by selecting combinations of oligos with the lowest weighted sum of variances for
    specified oligo properties, which ensures homogeneity within each set.

    :param set_size: The desired size of the oligo set to be generated.
    :type set_size: int
    :param properties: A dictionary of oligo properties (e.g., 'GC_content', 'length') and their respective weights.
    :type properties: Dict[str, float]
    """

    def __init__(self, set_size: int, properties: Dict[str, float]) -> None:
        """Constructor for the HomogeneousPropertyOligoSetGenerator class."""
        self.set_size = set_size
        self.properties = properties

    def apply(
        self, oligo_database: OligoDatabase, n_sets: int = 1, n_combinations: int = 1000, n_jobs: int = 1
    ) -> OligoDatabase:
        """
        Applies the oligo set generation process to the provided oligo database. For each region in the database,
        the method generates a set of oligos with lowest homogeneity of specified oligo properties. The process is
        parallelized across multiple regions using joblib.

        Oligosets are stored in the class attribute `oligosets`, which is a dictionary with region names as keys and oligoset
        dataframes as values. The structure of the pandas.DataFrame is the following:

        +-------------+----------+----------+----------+-------+----------+-------------+-------------+-------+
        | oligoset_id | oligo_0  | oligo_1  | oligo_2  |  ...  | oligo_n  | set_score_1 | set_score_2 |  ...  |
        +-------------+----------+----------+----------+-------+----------+-------------+-------------+-------+
        | 0           | AGRN_184 | AGRN_133 | AGRN_832 |  ...  | AGRN_706 | 0.3445      | 1.2332      |  ...  |
        +-------------+----------+----------+-----+----+-------+----------+-------------+-------------+-------+


        :param oligo_database: The OligoDatabase containing the oligonucleotides and their associated attributes.
        :type oligo_database: OligoDatabase
        :param n_sets: The number of oligo sets to generate.
        :type n_sets: int
        :param n_combinations: The number of random oligo combinations to generate per region (generating all combinations would be ideal but too costly), defaults to 1000.
        :type n_combinations: int, optional
        :param n_jobs: The number of parallel jobs to run, defaults to 1.
        :type n_jobs: int, optional

        :return: The updated oligo database with generated oligo sets for each region.
        :rtype: OligoDatabase
        """
        region_ids = list(oligo_database.database.keys())
        with joblib_progress(description="Find Oligosets", total=len(region_ids)):
            Parallel(n_jobs=n_jobs, prefer="threads", require="sharedmem")(
                delayed(self._get_oligo_sets_for_region)(oligo_database, region_id, n_sets, n_combinations)
                for region_id in region_ids
            )

        oligo_database.remove_regions_with_insufficient_oligos(pipeline_step="oligoset generation")
        return oligo_database

    def _get_oligo_sets_for_region(
        self, oligo_database: OligoDatabase, region_id: str, n_sets: int, n_combinations: int
    ) -> None:
        """
        Generates oligo sets for a specific region in the oligo database by scoring and sorting combinations
        of oligos based on the specified properties. The top n_sets sets with the lowest weighted sum of variances are selected.

        :param oligo_database: The OligoDatabase containing the oligonucleotides and their associated attributes.
        :type oligo_database: OligoDatabase
        :param region_id: Region ID to process.
        :type region_id: str
        :param n_sets: The number of oligo sets to generate.
        :type n_sets: int
        :param n_combinations: The number of random oligo combinations to generate per region (generating all combinations would be ideal but too costly), defaults to 1000.
        :type n_combinations: int, optional
        """
        oligo_df = pd.DataFrame({"oligo_id": oligo_database.database[region_id].keys()})
        oligo_df.set_index("oligo_id", inplace=True)

        # # check if all properties in self.properties are in oligo_df columns
        for property_name in self.properties:
            property_table = oligo_database.get_oligo_attribute_table(
                attribute=property_name, flatten=True, region_ids=region_id
            )
            property_table.set_index("oligo_id", inplace=True)

            if property_table[property_name].isnull().any():
                raise ValueError(
                    f"Property '{property_name}' is not present in oligo database please calculate it first using oligo_designer_toolsuite.OligoAttributes()."
                )
            else:
                if not (
                    pd.api.types.is_integer_dtype(property_table[property_name])
                    or pd.api.types.is_float_dtype(property_table[property_name])
                ):
                    raise ValueError(
                        f"Property '{property_name}' is not numeric. Cannot use for variance computation."
                    )
            oligo_df[property_name] = property_table[property_name]

        combinations = self._generate_random_combinations(oligo_df.index.to_list(), self.set_size, n_combinations)

        scored_combinations = [
            self._score_combination(oligo_df, list(combination)) for combination in combinations
        ]
        sorted_combinations = sorted(scored_combinations, key=lambda x: x[1], reverse=False)
        best_combinations = [combination for combination in sorted_combinations[:n_sets]]

        rows = [[idx] + oligos + [score] for idx, (oligos, score) in enumerate(best_combinations)]
        columns = ["oligoset_id"] + [f"oligo_{i}" for i in range(self.set_size)] + ["set_score"]

        oligo_database.oligosets[region_id] = pd.DataFrame(rows, columns=columns)

    def _score_combination(self, oligo_df: pd.DataFrame, combination: List[str]) -> Tuple[List[str], float]:
        """
        Scores a combination of oligos by calculating the variance of each oligo's properties in the set.

        :param oligo_df: The DataFrame containing the oligo information.
        :type oligo_df: pd.DataFrame
        :param combination: A list of oligo IDs for which the score is computed.
        :type combination: List[str]
        :return: A tuple containing the oligo combination and its score.
        :rtype: Tuple[List[str], float]
        """
        oligo_set = oligo_df.loc[combination]
        score = sum([oligo_set[property].var() * self.properties[property] for property in self.properties])
        return combination, score

    @staticmethod
    def _generate_random_combinations(arr, combination_size, number_of_combinations) -> list:
        """
        Generates oligo sets of specified size from random combinations of oligos.

        :param arr: The list of oligos to generate combinations from.
        :type arr: list
        :param combination_size: The size of each combination.
        :type combination_size: int
        :param number_of_combinations: The number of random combinations to generate (generating all combinations would be ideal but too costly).
        :type number_of_combinations: int
        :return: A list of random combinations.
        :rtype: list
        """
        total_combinations = comb(len(arr), combination_size)

        if total_combinations <= number_of_combinations:
            return list(itertools.combinations(arr, combination_size))

        seen_combinations = set()
        while len(seen_combinations) < number_of_combinations:
            combination = tuple(sorted(random.sample(list(arr), combination_size)))
            if combination not in seen_combinations:
                seen_combinations.add(combination)
        return list(seen_combinations)
