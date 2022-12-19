from abc import ABC, abstractmethod

import pandas as pd


class SetScoringBase(ABC):
    """Template class for scoring the stest of probes."""

    @abstractmethod
    def apply(clique_probes: pd.Series, n: int):
        """From a set of non-overlapping probes extracts the best subset of n probes and its score. Recieves in input
        a pandas series with the probes id as index and the score as value. The return is a list with the id of the selected
        probes and the scores of the set.
        REMARK: the set can have different scores with increasing relevance, where the less relevan scores are used in case of
        ties. Therefore in case the fist scores are equal the sording is done according to the second one and so on.

        :param clique_probes: Series of a set of non overlapping probes
        :type clique_probes: pandas.Series
        :param n: size of teh set
        :type n: int
        :return: optimal set and its score
        :rtype: list
        """


class PadlockSetScoring(SetScoringBase):
    """Scoring class for the sets of probes used for teh Padlock designer pipeline."""

    def __init__(self) -> None:
        pass

    def apply(self, clique_probes: pd.Series, n: int):
        """From a set of non-overlapping probes extracts the best subset of n probes and its scores. The scores are, in order of relevance,
         the maximal probe score in the set and the avreage of the scores.

        :param clique_probes: Series of a set of non overlapping probes
        :type clique_probes: pandas.Series
        :param n: size of teh set
        :type n: int
        :return: optimal set and its score
        :rtype: list
        """
        best_n_probes = clique_probes.sort_values(ascending=True).head(n)
        # Calculate performance of probeset
        probeset_error_max = best_n_probes.max()
        probeset_error_sum = best_n_probes.sum()
        probeset = best_n_probes.index.tolist()
        probeset += [probeset_error_max, probeset_error_sum]
        return probeset


class AverageSetScoring(SetScoringBase):
    """Scoring class for the sets of probes. It retunrs as score the average of hte  probe scores of the set."""

    def __init__(self) -> None:
        pass

    def apply(clique_probes: pd.Series, n: int):
        """From a set of non-overlapping probes extracts the best subset of n probes and its scores.
        The scores is the avreage of the scores of the probeas in the set.

        :param clique_probes: Series of a set of non overlapping probes
        :type clique_probes: pandas.Series
        :param n: size of teh set
        :type n: int
        :return: optimal set and its score
        :rtype: list
        """
        best_n_probes = clique_probes.sort_values(ascending=True).head(n)
        probeset_error_avg = best_n_probes.mean()
        probeset = best_n_probes.index.tolist()
        probeset += [probeset_error_avg]


class MaxSetScoring(SetScoringBase):
    """Scoring class for the sets of probes. It retunrs as score the max  of the probe scores of the set."""

    def __init__(self) -> None:
        pass

    def apply(clique_probes: pd.Series, n: int):
        """From a set of non-overlapping probes extracts the best subset of n probes and its scores.
        The scores is the max of the scores of the probea in the set.

        :param clique_probes: Series of a set of non overlapping probes
        :type clique_probes: pandas.Series
        :param n: size of teh set
        :type n: int
        :return: optimal set and its score
        :rtype: list
        """
        best_n_probes = clique_probes.sort_values(ascending=True).head(n)
        probeset_error_avg = best_n_probes.max()
        probeset = best_n_probes.index.tolist()
        probeset += [probeset_error_avg]
