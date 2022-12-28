from abc import ABC, abstractmethod

import numpy as np
import pandas as pd


class ProbeScoringBase(ABC):
    """Template class for scoring the probes."""

    def apply(self, probes: dict, probes_indices: np.array):
        """Scores all the probes using the defiend scoring function. The scores are both saved in the dictionary
        and in a pandas.Series. The latter is generated because it is the fastest way to generate the sets.

        :param probes: dictionary containing the probes
        :type probes: dict
        :param probes_indices: list of the indices of the probes used as a reference to keep the ordering of the probes consistent
        :type probes_indices: np.array
        :return: updated dictionary of the probes, series with the computed scores
        :rtype: dict, pandas.Series
        """

        probes_scores = pd.Series(index=probes_indices, dtype=float)
        for probe_id in probes_indices:
            score = self.scoring_function(probes[probe_id])
            probes[probe_id]["probe_score"] = score
            probes_scores[probe_id] = score
        return probes, probes_scores

    @abstractmethod
    def scoring_function(self, probe: dict):
        """Computes the score of the given probe

        :param probe: dictionary containing all the features of the given probe
        :type probe: dict
        :return: score of the probe
        :rtype: float
        """


class PadlockProbeScoring(ProbeScoringBase):
    """Probes scoring class for the padlock experiment.


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
        """
        Initialize the class
        """

        self.Tm_min = Tm_min
        self.Tm_opt = Tm_opt
        self.Tm_max = Tm_max
        self.GC_min = GC_content_min
        self.GC_opt = GC_content_opt
        self.GC_max = GC_content_max
        self.Tm_weight = Tm_weight
        self.GC_weight = GC_weight
        self.__generate_scoring_functions()

    def scoring_function(self, probe):
        """Computes the score of the given probe

        :param probe: dictionary containing all the features of the given probe
        :type probe: dict
        :return: score of the probe
        :rtype: float
        """
        # distance from the optimal melting temperature weightend by the how far is the optimum from the min/ max
        # teh scoring is the lower teh better
        Tm_dif = (
            probe["melting_temperature"] - self.Tm_opt
        )  # check the names of the columns
        GC_dif = probe["GC_content"] - self.GC_opt
        return self.Tm_weight * self.Tm_error(Tm_dif) + self.GC_weight * self.GC_error(
            GC_dif
        )

    def __generate_scoring_functions(self):
        """Computes relevant parts of the scoring function."""
        # define the error function for the melting temperature
        Tm_dif_max = self.Tm_max - self.Tm_opt
        Tm_dif_min = self.Tm_opt - self.Tm_min
        if Tm_dif_max == Tm_dif_min:
            self.Tm_error = lambda Tm_dif: abs(Tm_dif) / Tm_dif_max
        else:
            self.Tm_error = lambda Tm_dif: abs(Tm_dif) / Tm_dif_max * (
                Tm_dif > 0
            ) + abs(Tm_dif) / Tm_dif_min * (Tm_dif < 0)
        # define the error function for the GC content
        GC_dif_max = self.GC_max - self.GC_opt
        GC_dif_min = self.GC_opt - self.GC_min
        if GC_dif_max == GC_dif_min:
            self.GC_error = lambda GC_dif: abs(GC_dif) / GC_dif_max
        else:
            self.GC_error = lambda GC_dif: abs(GC_dif) / GC_dif_max * (
                GC_dif > 0
            ) + abs(GC_dif) / GC_dif_min * (GC_dif < 0)
