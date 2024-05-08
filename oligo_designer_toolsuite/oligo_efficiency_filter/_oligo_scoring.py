############################################
# imports
############################################

from abc import ABC, abstractmethod

import pandas as pd

############################################
# Oligo Scoring Classes
############################################


class OligoScoringBase(ABC):
    """Template class for scoring the oligos."""

    def apply(self, oligos: dict):
        """Scores all the oligos using the defiend scoring function. The scores are both saved in the dictionary
        and in a pandas.Series. The latter is generated because it is the fastest way to generate the sets.

        :param oligos: dictionary containing the oligos
        :type oligos: dict
        :return: updated dictionary of the oligos, series with the computed scores
        :rtype: dict, pandas.Series
        """
        oligos_indices = list(oligos.keys())
        oligos_scores = pd.Series(index=oligos_indices, dtype=float)
        for oligo_id in oligos_indices:
            score = self.scoring_function(oligos[oligo_id])
            oligos[oligo_id]["oligo_score"] = score
            oligos_scores[oligo_id] = score
        return oligos, oligos_scores

    @abstractmethod
    def scoring_function(self, oligo: dict):
        """Computes the score of the given oligo

        :param oligo: dictionary containing all the features of the given oligo
        :type oligo: dict
        :return: score of the oligo
        :rtype: float
        """


class PadlockOligoScoring(OligoScoringBase):
    """Oligos scoring class for the padlock experiment.

    :param Tm_min: minimal melting temperature
    :type Tm_min: float
    :param Tm_opt: minimal melting temperature
    :type Tm_opt: float
    :param Tm_max: maximal melting temperature
    :type Tm_max: float
    :param GC_min: minimal percentage of guanine and cytosine
    :type GC_min: float
    :param GC_opt: optimal percentage of guanine and cytosine
    :type GC_opt: float
    :param GC_max: maximal percentage of guanine and cytosine
    :type GC_max: float
    :param Tm_weight: relevance of the melting temperature in the scoring function, defaults to 1
    :type Tm_weight: int, optional
    :param GC_weight: relevance of the GC content in the scoring function, defaults to 1
    :type GC_weight: int, optional
    """

    def __init__(
        self,
        Tm_min: float,
        Tm_opt: float,
        Tm_max: float,
        GC_content_min: float,
        GC_content_opt: float,
        GC_content_max: float,
        Tm_weight: float = 1,
        GC_weight: float = 1,
    ):
        """Constructor method"""
        self.Tm_min = Tm_min
        self.Tm_opt = Tm_opt
        self.Tm_max = Tm_max
        self.GC_min = GC_content_min
        self.GC_opt = GC_content_opt
        self.GC_max = GC_content_max
        self.Tm_weight = Tm_weight
        self.GC_weight = GC_weight
        self.__generate_scoring_functions()

    def scoring_function(self, oligo):
        """Computes the score of the given oligo

        :param oligo: dictionary containing all the features of the given oligo
        :type oligo: dict
        :return: score of the oligo
        :rtype: float
        """
        # distance from the optimal melting temperature weightend by the how far is the optimum from the min/ max
        # the scoring is the lower the better
        Tm_dif = oligo["TmNN"] - self.Tm_opt  # check the names of the columns
        GC_dif = oligo["GC_content"] - self.GC_opt
        score = self.Tm_weight * self.Tm_error(Tm_dif) + self.GC_weight * self.GC_error(GC_dif)
        return score

    def __generate_scoring_functions(self):
        """Computes relevant parts of the scoring function."""
        # define the error function for the melting temperature
        Tm_dif_max = self.Tm_max - self.Tm_opt
        Tm_dif_min = self.Tm_opt - self.Tm_min
        if Tm_dif_max == Tm_dif_min:
            self.Tm_error = lambda Tm_dif: abs(Tm_dif) / Tm_dif_max
        else:
            self.Tm_error = lambda Tm_dif: abs(Tm_dif) / Tm_dif_max * (Tm_dif > 0) + abs(
                Tm_dif
            ) / Tm_dif_min * (Tm_dif < 0)
        # define the error function for the GC content
        GC_dif_max = self.GC_max - self.GC_opt
        GC_dif_min = self.GC_opt - self.GC_min
        if GC_dif_max == GC_dif_min:
            self.GC_error = lambda GC_dif: abs(GC_dif) / GC_dif_max
        else:
            self.GC_error = lambda GC_dif: abs(GC_dif) / GC_dif_max * (GC_dif > 0) + abs(
                GC_dif
            ) / GC_dif_min * (GC_dif < 0)


class SeqFISHOligoScoring(OligoScoringBase):
    """Oligos scoring class for the SeqFISH+ experiment.
    Scoring function has the following form: ((GC_content_of_sequence - GC_opt)/(GC_max-GC_min))^2

    :param GC_min: minimal percentage of guanine and cytosine
    :type GC_min: float
    :param GC_opt: optimal percentage of guanine and cytosine
    :type GC_opt: float
    :param GC_max: maximal percentage of guanine and cytosine
    :type GC_max: float

    """

    def __init__(
        self,
        GC_content_opt: float,
    ):
        """
        Initialize the class
        """
        self.GC_opt = GC_content_opt

    def scoring_function(self, oligo):
        """Computes the score of the given oligo

        :param oligo: dictionary containing all the features of the given oligo
        :type oligo: dict
        :return: score of the oligo
        :rtype: float
        """
        # distance from optimal GC (the lower the better)
        return abs(oligo["GC_content"] - self.GC_opt)
