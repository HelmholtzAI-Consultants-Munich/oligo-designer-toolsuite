from abc import ABC, abstractmethod

import pandas as pd


class ProbeScoring(ABC):
    """Template class for scoring the probes and the sets of probes"""

    def __init__(self):
        """Initialize the class"""

    def apply(self, probes, probes_indices):
        """Scores all the probes using the defiend scoring function. The scores are both saved in the dictionary
        and in a pandas.Series, which will be used for convenience in the selection of the sets.

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
    def scoring_function(self, probe):
        """Computes the score of the given probe

        :param probe: dictionary containing all th efeatures of the given probe
        :type probe: dict
        :return: score of the probe
        :rtype: float
        """

    @abstractmethod
    def get_probeset(clique_probes, n):
        """From a set of non-overlapping probes extracts the best subset of n probes and its score. Recieven in input
        a pandas serues with the probes id as index and the score as value. The return is a list with the id of the selected
        probes and the score of the set.
        REMARK: the set can have different scores with increasing relevance, whenre the less relevan scores are used in case of
        ties, in this case the scores are saved in a sublist which is compared in lexographiocal order in the sorting procedure,
        therefore in case the fist scores are equal the sording is done according to the second one and so on.

        :param clique_probes: Series of a set of non overlapping probes
        :type clique_probes: pandas.Series
        :param n: size of teh set
        :type n: int
        :return: optimal set and its score
        :rtype: list
        """


class PadlockScoring(ProbeScoring):
    """Scoring class for padlock."""

    def __init__(
        self, Tm_min, Tm_opt, Tm_max, GC_min, GC_opt, GC_max, Tm_weight=1, GC_weight=1
    ):
        """Initialize the class

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

        super().__init__()
        self.Tm_min = Tm_min
        self.Tm_opt = Tm_opt
        self.Tm_max = Tm_max
        self.GC_min = GC_min
        self.GC_opt = GC_opt
        self.GC_max = GC_max
        self.Tm_weight = Tm_weight
        self.GC_weight = GC_weight
        self.generate_scoring_functions()

    def scoring_function(self, probe):
        """Computes the score of the given probe

        :param probe: dictionary containing all th efeatures of the given probe
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

    def generate_scoring_functions(self):
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

    def get_probeset(self, clique_probes, n):
        """From a set of non-overlapping probes extracts the best subset of n probes and its scores. The scores are, in order of relevance,
         the maximal probe score in the set and the avreage of the scores

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
        probeset.append([probeset_error_max, probeset_error_sum])
        return probeset
