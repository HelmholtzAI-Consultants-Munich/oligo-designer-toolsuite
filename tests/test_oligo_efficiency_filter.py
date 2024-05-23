############################################
# imports
############################################

import unittest
from pandas import Series

from Bio.SeqUtils import MeltingTemp as mt

from oligo_designer_toolsuite.oligo_efficiency_filter import (
    AverageSetScoring,
    LowestSetScoring,
    WeightedTmGCOligoScoring,
    GCOligoScoring,
)

############################################
# Global Parameters
############################################

TM_PARAMETERS = {
    "check": True,
    "strict": True,
    "c_seq": None,
    "shift": 0,
    "nn_table": getattr(mt, "DNA_NN3"),
    "tmm_table": getattr(mt, "DNA_TMM1"),
    "imm_table": getattr(mt, "DNA_IMM1"),
    "de_table": getattr(mt, "DNA_DE1"),
    "dnac1": 50,  # [nM]
    "dnac2": 0,
    "selfcomp": False,
    "saltcorr": 7,
    "Na": 50,  # [mM]
    "K": 75,  # [mM]
    "Tris": 20,  # [mM]
    "Mg": 10,  # [mM]
    "dNTPs": 0,
}

############################################
# Tests
############################################


class TestOligoScoring(unittest.TestCase):

    def setUp(self):
        self.score_gc = GCOligoScoring(GC_content_opt=43.75)
        self.score_weighted_gc_tm = WeightedTmGCOligoScoring(
            Tm_min=50,
            Tm_opt=65,
            Tm_max=70,
            GC_content_min=25,
            GC_content_opt=50,
            GC_content_max=75,
            Tm_parameters=TM_PARAMETERS,
        )

        self.oligo = "TATTTCGGGCATGCAT"

    def test_GC_score(self):
        oligo_score = self.score_gc.scoring_function(self.oligo)
        assert oligo_score == 0, "GC score failed!"

    def test_weighted_GC_Tm_score(self):
        oligo_score = self.score_weighted_gc_tm.scoring_function(self.oligo)
        assert abs(oligo_score - 0.74666) < 1e-5, "Weighted GC-Tm score failed!"


class TestSetScoring(unittest.TestCase):

    def setUp(self):
        self.score_max_sum = LowestSetScoring()
        self.score_ave_max = AverageSetScoring()
        self.oligo_set = Series(data=[0, 1, 8, 5, 2, 6, 7, 3])
        self.n_oligo_set = 5

    def test_max_sum(self):
        oligoset = self.score_max_sum.apply(self.oligo_set, self.n_oligo_set)
        assert oligoset == [0, 1, 4, 7, 3, 5, 11], "Max scoring failed"

    def test_ave_max(self):
        oligoset = self.score_ave_max.apply(self.oligo_set, self.n_oligo_set)
        assert oligoset == [0, 1, 4, 7, 3, 2.2, 5], "Average scoring failed"
