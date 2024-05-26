############################################
# imports
############################################

from abc import ABC, abstractmethod

import pandas as pd

############################################
# Set Scoring Classes
############################################


class SetScoringBase(ABC):
    """Abstract base class for scoring sets of oligonucleotides.
    Implementations should define how to extract the best subset of oligonucleotides based on their scores.
    """

    @abstractmethod
    def apply(self, oligo_set: pd.Series, n: int):
        """
        Extracts the top 'n' oligonucleotides from a provided pandas Series based on their scores.
        The method returns the IDs and scores of the selected oligonucleotides.
        Scores can be sorted in ascending or descending order depending on their meaning.
        REMARK: the set can have different scores with increasing relevance, where the less relevant scores are used in case of
        ties. Therefore in case the fist scores are equal the sorting is done according to the second one and so on.

        :param oligo_set: Pandas Series with oligonucleotide IDs as index and their scores as values.
        :type oligo_set: pd.Series
        :param n: The number of top oligonucleotides to select from the set.
        :type n: int
        :return: List of tuples with selected oligonucleotide IDs and their respective scores.
        :rtype: list
        """


class LowestSetScoring(SetScoringBase):
    """Implements the SetScoringBase to score a set of oligonucletides by it's lowest oligo score (dependent on the
    meaning of the score). In case of ties, the sum of all oligo scores in the set is provided as well.

    :param ascending: If True, scores are sorted in ascending order; if False, in descending order. This depens on the meaning of the score.
    :type ascending: bool
    """

    def __init__(self, ascending: bool) -> None:
        """Constructor for the LowestSetScoring class."""
        self.ascending = ascending

    def apply(self, oligo_set: pd.Series, n: int):
        """Selects the top 'n' oligonucleotides based on their scores,
        The most relevant set score is the maximum/minimum oligo score in the set,
        and the subsequent set score is the sum of oligo scores in the set.

        :param oligo_set: Pandas Series with oligonucleotide IDs as indices and their scores as values.
        :type oligo_set: pd.Series
        :param n: Number of oligonucleotides to select.
        :type n: int
        :return: List containing the IDs of the selected oligonucleotides, followed by the maximum/minimum oligo score and the sum of oligo scores in the set.
        :rtype: list
        """
        best_n_oligos = oligo_set.sort_values(ascending=self.ascending).head(n)

        # Calculate performance of oligoset
        if self.ascending:
            set_score_lowest = best_n_oligos.max()
        else:
            set_score_lowest = best_n_oligos.min()

        set_score_sum = best_n_oligos.sum()
        oligoset = best_n_oligos.index.tolist()
        oligoset += [round(set_score_lowest, 4), round(set_score_sum, 4)]
        return oligoset


class AverageSetScoring(SetScoringBase):
    """Implements the SetScoringBase to score a set of oligonucletides by the average of all oligo scores in the set.
    In case of ties, lowest oligo score (dependent on the meaning of the score) of the set is provided as well.

    :param ascending: If True, scores are sorted in ascending order; if False, in descending order. This depens on the meaning of the score.
    :type ascending: bool
    """

    def __init__(self, ascending: bool) -> None:
        """Constructor for the AverageSetScoring class."""
        self.ascending = ascending

    def apply(self, oligo_set: pd.Series, n: int):
        """Selects the top 'n' oligonucleotides based on their scores,
        The most relevant set score is the average oligo score in the set,
        and the subsequent set score is the sum of oligo scores in the set.

        :param oligo_set: Pandas Series with oligonucleotide IDs as indices and their scores as values.
        :type oligo_set: pd.Series
        :param n: Number of oligonucleotides to select.
        :type n: int
        :return: List containing the IDs of the selected oligonucleotides, followed by the average oligo score and the sum of oligo scores in the set.
        :rtype: list
        """
        best_n_oligos = oligo_set.sort_values(ascending=self.ascending).head(n)

        # Calculate performance of oligoset
        if self.ascending:
            set_score_lowest = best_n_oligos.max()
        else:
            set_score_lowest = best_n_oligos.min()
        set_score_avg = best_n_oligos.mean()

        oligoset = best_n_oligos.index.tolist()
        oligoset += [round(set_score_avg, 4), round(set_score_lowest, 4)]
        return oligoset
