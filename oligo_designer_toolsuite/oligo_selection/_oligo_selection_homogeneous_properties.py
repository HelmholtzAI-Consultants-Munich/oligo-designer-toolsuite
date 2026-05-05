############################################
# imports
############################################

import itertools
import random
from math import comb
from typing import cast

import pandas as pd

from oligo_designer_toolsuite._exceptions import DatabaseError
from oligo_designer_toolsuite.database import OligoDatabase

from ._oligo_selection_base import BaseOligoSelection

############################################
# Oligoset Selection Classes
############################################


class HomogeneousPropertyOligoSelection(BaseOligoSelection):
    """
    Generates oligo sets based on the homogeneity of specified oligo properties. The oligo sets are
    created by selecting combinations of oligos with the lowest weighted sum of variances for
    specified oligo properties, which ensures homogeneity within each set.

    :param set_size: The desired size of the oligo set to be generated.
    :type set_size: int
    :param properties: A dictionary of oligo properties (e.g., 'GC_content', 'length') and their respective weights.
    :type properties: dict[str, float]
    :param n_combinations: The number of random oligo combinations to sample per region. Default is 1000.
    :type n_combinations: int, optional
    """

    def __init__(
        self,
        set_size: int,
        properties: dict[str, float],
        n_combinations: int = 1000,
    ) -> None:
        """Constructor for the HomogeneousPropertyOligoSetGenerator class."""
        self.set_size = set_size
        self.properties = properties
        self.n_combinations = n_combinations

    def _get_oligo_sets_for_region(
        self,
        oligo_database: OligoDatabase,
        sequence_type: str,
        region_id: str,
        n_sets: int,
    ) -> None:
        """
        Generates oligo sets for a specific region by scoring and sorting random combinations of oligos
        based on the specified properties. The top n_sets sets with the lowest weighted sum of
        variances are selected and stored in ``oligo_database.oligosets[region_id]``. Uses
        :attr:`n_combinations` from the constructor for the number of random combinations per region.

        :param oligo_database: The OligoDatabase instance containing oligonucleotide sequences and
            their associated properties, organized by genomic regions.
        :type oligo_database: OligoDatabase
        :param sequence_type: Type of sequence being processed.
        :type sequence_type: str
        :param region_id: Region ID to process.
        :type region_id: str
        :param n_sets: The number of oligo sets to generate.
        :type n_sets: int
        """
        oligo_df = oligo_database.get_oligo_property_table(
            properties=list(self.properties.keys()), flatten=True, region_ids=region_id
        )

        # # check if all properties in self.properties are in oligo_df columns
        for property in self.properties:
            if oligo_df[property].isnull().any():
                raise DatabaseError(
                    f"Property '{property}' is not present in oligo database. "
                    f"Please calculate it first using oligo_designer_toolsuite.oligo_property_calculator.PropertyCalculator()."
                )
            else:
                if not (
                    pd.api.types.is_integer_dtype(oligo_df[property])
                    or pd.api.types.is_float_dtype(oligo_df[property])
                ):
                    raise DatabaseError(
                        f"Property '{property}' is not numeric. Cannot use for variance computation. "
                        f"Properties used for variance computation must be numeric (integer or float)."
                    )

        combinations = self._generate_random_combinations(
            list(oligo_df.index), self.set_size, self.n_combinations
        )

        scored_combinations = [
            self._score_combination(oligo_df, list(combination)) for combination in combinations
        ]
        sorted_combinations = sorted(scored_combinations, key=lambda x: x[1], reverse=False)
        best_combinations = [combination for combination in sorted_combinations[:n_sets]]

        rows = [[idx] + oligos + [score] for idx, (oligos, score) in enumerate(best_combinations)]
        columns = ["oligoset_id"] + [f"oligo_{i}" for i in range(self.set_size)] + ["set_score"]

        oligo_database.oligosets[region_id] = pd.DataFrame(rows, columns=columns)

    def _score_combination(self, oligo_df: pd.DataFrame, combination: list[str]) -> tuple[list[str], float]:
        """
        Score a candidate oligo set by measuring how homogeneous its properties are.

        For the given list of oligo IDs, this method computes the variance of each configured
        property across the set, weights each variance by the weight of the property, and sums
        them to obtain a single scalar score. Lower scores correspond to more homogeneous
        (less variable) oligo sets.

        :param oligo_df: The DataFrame containing the oligo information.
        :type oligo_df: pd.DataFrame
        :param combination: A list of oligo IDs for which the score is computed.
        :type combination: list[str]
        :return: A tuple containing the oligo combination and its score.
        :rtype: tuple[list[str], float]
        """
        oligo_set = oligo_df.loc[combination]
        score = sum(
            [
                cast(float, oligo_set[property].var()) * self.properties[property]
                for property in self.properties
            ]
        )
        return combination, score

    @staticmethod
    def _generate_random_combinations(
        arr: list[str],
        combination_size: int,
        number_of_combinations: int,
        seed: int = 42,
    ) -> list[tuple[str, ...]]:
        """
        Generates oligo sets of specified size from random combinations of oligos.

        :param arr: The list of oligos to generate combinations from.
        :type arr: list[str]
        :param combination_size: The size of each combination.
        :type combination_size: int
        :param number_of_combinations: The number of random combinations to generate.
        :type number_of_combinations: int
        :param seed: The seed for the random number generator.
        :type seed: int, optional
        :return: A list of random combinations.
        :rtype: list[tuple[str, ...]]
        """
        random.seed(seed)
        total_combinations = comb(len(arr), combination_size)

        if total_combinations <= number_of_combinations:
            return list(itertools.combinations(arr, combination_size))

        seen_combinations = set[tuple[str, ...]]()
        while len(seen_combinations) < number_of_combinations:
            combination = tuple(sorted(random.sample(list(arr), combination_size)))
            if combination not in seen_combinations:
                seen_combinations.add(combination)
        return list(seen_combinations)
