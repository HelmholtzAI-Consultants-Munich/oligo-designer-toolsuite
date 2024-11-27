############################################
# imports
############################################

import pandas as pd

from abc import ABC, abstractmethod
from typing import Tuple

############################################
# Set Scoring Classes
############################################


class SetScoringBase(ABC):
    """
    Abstract base class for scoring oligonucleotide sets. Subclasses must implement the `apply` method.
    Implementations should define how to extract the best subset of oligonucleotides based on their scores.
    """

    @abstractmethod
    def apply(self, oligo_set: pd.Series, n: int) -> Tuple[list, dict]:
        """
        Abstract method to apply the scoring method to the oligo set and return the top `n` oligos along with their respective scores.
        This method should be implemented by subclasses to define the specific scoring logic.


        :param oligo_set: A pandas Series representing the set of oligos with their scores.
        :type oligo_set: pd.Series
        :param n: The number of top oligos to select.
        :type n: int
        :return: A tuple containing a list of selected oligo IDs and a dictionary of computed scores.
        :rtype: Tuple[list, dict]
        """


class LowestSetScoring(SetScoringBase):
    """
    Scores oligonucleotide sets based on the lowest score in the set (dependent on the meaning of the score).
    In case of ties, the sum of all oligo scores in the set is provided as well.

    :param ascending: Whether to sort the oligos in ascending or descending order (dependent on the meaning of the score) before selecting the top `n`.
    :type ascending: bool
    """

    def __init__(self, ascending: bool) -> None:
        """
        Constructor for the LowestSetScoring class.
        """
        self.ascending = ascending

    def apply(self, oligo_set: pd.Series, n: int) -> Tuple[list, dict]:
        """
        Apply the scoring method to the oligo set, selecting the top `n` oligos based on the lowest score.

        :param oligo_set: A pandas Series representing the set of oligos with their scores.
        :type oligo_set: pd.Series
        :param n: The number of top oligos to select.
        :type n: int
        :return: A tuple containing a list of selected oligo IDs and a dictionary with the lowest score and the sum of scores.
        :rtype: Tuple[list, dict]
        """
        best_n_oligos = oligo_set.sort_values(ascending=self.ascending).head(n)

        # Calculate performance of oligoset
        if self.ascending:
            set_score_lowest = best_n_oligos.max()
        else:
            set_score_lowest = best_n_oligos.min()

        set_score_sum = best_n_oligos.sum()
        oligoset = best_n_oligos.index.tolist()
        return oligoset, {
            "set_score_worst": round(set_score_lowest, 4),
            "set_score_sum": round(set_score_sum, 4),
        }


class AverageSetScoring(SetScoringBase):
    """
    Scores oligonucleotide sets based on the average score in the set.
    In case of ties, lowest oligo score (dependent on the meaning of the score) of the set is provided as well.

    :param ascending: Whether to sort the oligos in ascending or descending order (dependent on the meaning of the score) before selecting the top `n`.
    :type ascending: bool
    """

    def __init__(self, ascending: bool) -> None:
        """Constructor for the AverageSetScoring class."""
        self.ascending = ascending

    def apply(self, oligo_set: pd.Series, n: int) -> Tuple[list, dict]:
        """
        Apply the scoring method to the oligo set, selecting the top `n` oligos based on the average score.

        :param oligo_set: A pandas Series representing the set of oligos with their scores.
        :type oligo_set: pd.Series
        :param n: The number of top oligos to select.
        :type n: int
        :return: A tuple containing a list of selected oligo IDs and a dictionary with the average score and the lowest score.
        :rtype: Tuple[list, dict]
        """
        best_n_oligos = oligo_set.sort_values(ascending=self.ascending).head(n)

        # Calculate performance of oligoset
        if self.ascending:
            set_score_lowest = best_n_oligos.max()
        else:
            set_score_lowest = best_n_oligos.min()
        set_score_avg = best_n_oligos.mean()

        oligoset = best_n_oligos.index.tolist()
        return oligoset, {
            "set_score_average": round(set_score_avg, 4),
            "set_score_worst": round(set_score_lowest, 4),
        }
