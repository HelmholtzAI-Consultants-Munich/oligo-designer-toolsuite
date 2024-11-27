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
    """ """

    @abstractmethod
    def apply(self, oligo_set: pd.Series, n: int) -> Tuple[list, dict]:
        """ """


class LowestSetScoring(SetScoringBase):
    """ """

    def __init__(self, ascending: bool) -> None:
        """
        Constructor for the LowestSetScoring class.
        """
        self.ascending = ascending

    def apply(self, oligo_set: pd.Series, n: int) -> Tuple[list, dict]:

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

    def __init__(self, ascending: bool) -> None:
        """Constructor for the AverageSetScoring class."""
        self.ascending = ascending

    def apply(self, oligo_set: pd.Series, n: int) -> Tuple[list, dict]:

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
