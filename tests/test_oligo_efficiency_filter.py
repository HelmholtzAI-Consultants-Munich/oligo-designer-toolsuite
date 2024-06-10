############################################
# imports
############################################

import unittest

from Bio.SeqUtils import MeltingTemp as mt
from pandas import Series

from oligo_designer_toolsuite.oligo_efficiency_filter import (
    AverageSetScoring,
    GCOligoScoring,
    LowestSetScoring,
    WeightedIsoformTmGCOligoScoring,
    WeightedTmGCOligoScoring,
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
        self.score_weighted_isoform_gc_tm = WeightedIsoformTmGCOligoScoring(
            Tm_min=50,
            Tm_opt=65,
            Tm_max=70,
            GC_content_min=25,
            GC_content_opt=50,
            GC_content_max=75,
            Tm_parameters=TM_PARAMETERS,
        )

        # in this test case we test for one oligos sequence that is found at two different sites in one genomic region
        # one site covers two transcripts and the other site only one transcript (which is already covered by the first site)
        # hence with 4 transcripts in total we should get an isoform consensus of 0.5
        self.oligo_attributes = {
            "oligo": "TATTTCGGGCATGCAT",
            "transcript_id": [["tr1"], ["tr1", "tr3"]],
            "number_transcripts": [[4], [4]],
        }
        self.sequence_type = "oligo"

    def test_GC_score(self):
        oligo_score = self.score_gc.get_score(
            oligo_attributes=self.oligo_attributes, sequence_type=self.sequence_type
        )
        assert oligo_score == 0, "GC score failed!"

    def test_weighted_GC_Tm_score(self):
        oligo_score = self.score_weighted_gc_tm.get_score(
            oligo_attributes=self.oligo_attributes, sequence_type=self.sequence_type
        )
        assert abs(oligo_score - 0.74666) < 1e-5, "Weighted GC-Tm score failed!"

    def test_weighted_isoform_GC_Tm_score(self):
        oligo_score = self.score_weighted_isoform_gc_tm.get_score(
            oligo_attributes=self.oligo_attributes, sequence_type=self.sequence_type
        )
        assert abs(oligo_score - 1.74666) < 1e-5, "Weighted GC-Tm score failed!"


class TestSetScoring(unittest.TestCase):
    def setUp(self):
        self.score_max_sum = LowestSetScoring(ascending=True)
        self.score_ave_max = AverageSetScoring(ascending=True)
        self.oligo_set = Series(data=[0, 1, 8, 5, 2, 6, 7, 3])
        self.n_oligo_set = 5

    def test_max_sum(self):
        oligoset = self.score_max_sum.apply(self.oligo_set, self.n_oligo_set)
        print(oligoset[1])
        assert oligoset[0] == [0, 1, 4, 7, 3], "Max scoring failed"
        assert oligoset[1] == {"set_score_lowest": 5, "set_score_sum": 11}, "Max scoring failed"

    def test_ave_max(self):
        oligoset = self.score_ave_max.apply(self.oligo_set, self.n_oligo_set)
        assert oligoset[0] == [0, 1, 4, 7, 3], "Average scoring failed"
        assert oligoset[1] == {"set_score_average": 2.2, "set_score_lowest": 5}, "Average scoring failed"
