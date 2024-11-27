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

    def __init__(
        self,
        oligos_scoring: OligoScoringBase,
        set_scoring: SetScoringBase,
        selection_policy: OligoSelectionPolicy,
        distance_between_oligos: int = 0,
        max_oligos: int = None,
    ) -> None:

        self.selection_policy = selection_policy
        self.oligos_scoring = oligos_scoring
        self.set_scoring = set_scoring
        self.distance_between_oligos = distance_between_oligos
        self.ascending = set_scoring.ascending
        self.max_oligos = max_oligos

    def apply(
        self,
        oligo_database: OligoDatabase,
        sequence_type: _TYPES_SEQ,
        n_attempts: int = 10000,
        n_jobs: int = 1,
    ):

        region_ids = list(oligo_database.database.keys())
        # get the oligo set for this region in parallel
        with joblib_progress(description="Find Oligosets", total=len(region_ids)):
            Parallel(
                n_jobs=n_jobs, prefer="threads", require="sharedmem"
            )(  # there should be an explicit return
                delayed(self._get_oligo_set_for_gene)(oligo_database, region_id, sequence_type, n_attempts)
                for region_id in region_ids
            )

        oligo_database.remove_regions_with_insufficient_oligos(pipeline_step="oligoset generation")
        return oligo_database

    def _get_oligo_set_for_gene(
        self,
        oligo_database: OligoDatabase,
        region_id: str,
        sequence_type: _TYPES_SEQ,
        n_attempts: int,
    ):

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
        oligosets = self.selection_policy.apply(
            oligos_scores=oligos_scores,
            overlapping_matrix=overlapping_matrix,
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

        return overlapping_matrix


class HomogeneousPropertyOligoSetGenerator:

    def __init__(self, set_size: int, properties: Dict[str, float]) -> None:

        self.set_size = set_size
        self.properties = properties

    def apply(
        self, oligo_database: OligoDatabase, n_attempts: int = 1, n_combinations: int = 1000, n_jobs: int = 1
    ) -> OligoDatabase:

        region_ids = list(oligo_database.database.keys())
        with joblib_progress(description="Find Oligosets", total=len(region_ids)):
            Parallel(n_jobs=n_jobs, prefer="threads", require="sharedmem")(
                delayed(self._get_oligo_sets_for_region)(
                    oligo_database, region_id, n_attempts, n_combinations
                )
                for region_id in region_ids
            )

        oligo_database.remove_regions_with_insufficient_oligos(pipeline_step="oligoset generation")
        return oligo_database

    def _get_oligo_sets_for_region(
        self, oligo_database: OligoDatabase, region_id: str, n_attempts: int, n_combinations: int
    ) -> None:

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
