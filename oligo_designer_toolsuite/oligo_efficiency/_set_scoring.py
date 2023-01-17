from abc import ABC, abstractmethod

import pandas as pd


class SetScoringBase(ABC):
    """Template class for scoring the stest of oligos."""

    @abstractmethod
    def apply(self, clique_oligos: pd.Series, n: int):
        """From a set of non-overlapping oligos extracts the best subset of n oligos and its score. Recieves in input
        a pandas series with the oligos id as index and the score as value. The return is a list with the id of the selected
        oligos and the scores of the set.
        REMARK: the set can have different scores with increasing relevance, where the less relevan scores are used in case of
        ties. Therefore in case the fist scores are equal the sording is done according to the second one and so on.

        :param clique_oligos: Series of a set of non overlapping oligos
        :type clique_oligos: pandas.Series
        :param n: size of teh set
        :type n: int
        :return: optimal set and its score
        :rtype: list
        """


class PadlockSetScoring(SetScoringBase):
    """Scoring class for the sets of oligos used for teh Padlock designer pipeline."""

    def __init__(self) -> None:
        pass

    def apply(self, clique_oligos: pd.Series, n: int):
        """From a set of non-overlapping oligos extracts the best subset of n oligos and its scores. The scores are, in order of relevance,
         the maximal oligo score in the set and the avreage of the scores.

        :param clique_oligos: Series of a set of non overlapping oligos
        :type clique_oligos: pandas.Series
        :param n: size of teh set
        :type n: int
        :return: optimal set and its score
        :rtype: list
        """
        best_n_oligos = clique_oligos.sort_values(ascending=True).head(n)
        # Calculate performance of oligoset
        oligoset_error_max = best_n_oligos.max()
        oligoset_error_sum = best_n_oligos.sum()
        oligoset = best_n_oligos.index.tolist()
        oligoset += [oligoset_error_max, oligoset_error_sum]
        return oligoset


class AverageSetScoring(SetScoringBase):
    """Scoring class for the sets of oligos. It retunrs as score the average of hte  oligo scores of the set."""

    def __init__(self) -> None:
        pass

    def apply(self, clique_oligos: pd.Series, n: int):
        """From a set of non-overlapping oligos extracts the best subset of n oligos and its scores.
        The scores is the avreage of the scores of the oligoas in the set.

        :param clique_oligos: Series of a set of non overlapping oligos
        :type clique_oligos: pandas.Series
        :param n: size of teh set
        :type n: int
        :return: optimal set and its score
        :rtype: list
        """
        best_n_oligos = clique_oligos.sort_values(ascending=True).head(n)
        oligoset_error_avg = best_n_oligos.mean()
        oligoset = best_n_oligos.index.tolist()
        oligoset += [oligoset_error_avg]
        return oligoset


class MaxSetScoring(SetScoringBase):
    """Scoring class for the sets of oligos. It retunrs as score the max  of the oligo scores of the set."""

    def __init__(self) -> None:
        pass

    def apply(self, clique_oligos: pd.Series, n: int):
        """From a set of non-overlapping oligos extracts the best subset of n oligos and its scores.
        The scores is the max of the scores of the oligoa in the set.

        :param clique_oligos: Series of a set of non overlapping oligos
        :type clique_oligos: pandas.Series
        :param n: size of teh set
        :type n: int
        :return: optimal set and its score
        :rtype: list
        """
        best_n_oligos = clique_oligos.sort_values(ascending=True).head(n)
        oligoset_error_avg = best_n_oligos.max()
        oligoset = best_n_oligos.index.tolist()
        oligoset += [oligoset_error_avg]
        return oligoset
