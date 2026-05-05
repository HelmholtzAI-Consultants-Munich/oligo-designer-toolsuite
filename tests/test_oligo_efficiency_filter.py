############################################
# imports
############################################

import os
import shutil
import unittest

from Bio.SeqUtils import MeltingTemp as mt
from pandas import Series
from scipy.sparse import csr_matrix

from oligo_designer_toolsuite.database import OligoDatabase
from oligo_designer_toolsuite.oligo_efficiency_filter import (
    AverageSetScoring,
    DeviationFromOptimalGCContentScorer,
    DeviationFromOptimalTmScorer,
    IsoformConsensusScorer,
    LowestSetScoring,
    NormalizedDeviationFromOptimalGCContentScorer,
    NormalizedDeviationFromOptimalTmScorer,
    OligoScoring,
    OverlapTargetedExonsScorer,
    OverlapUTRScorer,
    UniformDistanceScorer,
)

############################################
# Global Parameters
############################################

FILE_DATABASE = "tests/data/databases/database_oligo_efficiency.tsv"

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
    def setUp(self) -> None:
        self.tmp_path = os.path.join(os.getcwd(), "tmp_oligo_efficiency")

        # in this test case we test for one oligos sequence that is found at two different sites in one genomic region
        # one site covers two transcripts and the other site only one transcript (which is already covered by the first site)
        # hence with 4 transcripts in total we should get an isoform consensus of 0.5
        self.oligo_database = OligoDatabase(
            min_oligos_per_region=2,
            write_regions_with_insufficient_oligos=True,
            max_entries_in_memory=10,
            n_jobs=4,
            database_name="test_oligo_database",
            dir_output=self.tmp_path,
        )
        self.oligo_database.set_database_sequence_types(["oligo", "target"])
        self.oligo_database.load_database_from_table(
            FILE_DATABASE, region_ids=None, database_overwrite=True, merge_databases_on_sequence_type="oligo"
        )
        self.sequence_type = "oligo"

        self.utr_scorer = OverlapUTRScorer(score_weight=10)
        self.exon_scorer = OverlapTargetedExonsScorer(targeted_exons=["21", "4"], score_weight=2)
        self.isoform_consensus_scorer = IsoformConsensusScorer(score_weight=2)
        self.Tm_scorer = DeviationFromOptimalTmScorer(
            Tm_opt=57.55,
            Tm_parameters=TM_PARAMETERS,
            Tm_salt_correction_parameters=None,
            Tm_chem_correction_parameters=None,
            score_weight=1,
        )
        self.GC_scorer = DeviationFromOptimalGCContentScorer(
            GC_content_opt=43.75,
            score_weight=1,
        )
        self.Tm_norm_scorer = NormalizedDeviationFromOptimalTmScorer(
            Tm_min=50,
            Tm_opt=65,
            Tm_max=70,
            Tm_parameters=TM_PARAMETERS,
            Tm_chem_correction_parameters=None,
            Tm_salt_correction_parameters=None,
            score_weight=1,
        )
        self.GC_norm_scorer = NormalizedDeviationFromOptimalGCContentScorer(
            GC_content_min=25,
            GC_content_opt=50,
            GC_content_max=80,
            score_weight=1,
        )
        self.uniform_distance_scorer = UniformDistanceScorer(average_oligo_length=5.0, score_weight=1.0)

    def tearDown(self) -> None:
        shutil.rmtree(self.tmp_path)

    def test_utr_scorer(self) -> None:
        oligo_score_wo_utr = self.utr_scorer.apply(
            oligo_database=self.oligo_database,
            region_id="region_1",
            oligo_id="region_1::1",
            sequence_type=self.sequence_type,
        )
        oligo_score_with_utr = self.utr_scorer.apply(
            oligo_database=self.oligo_database,
            region_id="region_1",
            oligo_id="region_1::2",
            sequence_type=self.sequence_type,
        )
        assert oligo_score_with_utr == 0, "error: scoring for region wo UTR incorrect."
        assert oligo_score_with_utr < oligo_score_wo_utr, "error: UTR score incorrect."

    def test_exon_scorer(self) -> None:
        oligo_score_wo_exon = self.exon_scorer.apply(
            oligo_database=self.oligo_database,
            region_id="region_1",
            oligo_id="region_1::1",
            sequence_type=self.sequence_type,
        )
        oligo_score_with_exon1 = self.exon_scorer.apply(
            oligo_database=self.oligo_database,
            region_id="region_1",
            oligo_id="region_1::2",
            sequence_type=self.sequence_type,
        )
        oligo_score_with_exon2 = self.exon_scorer.apply(
            oligo_database=self.oligo_database,
            region_id="region_1",
            oligo_id="region_1::3",
            sequence_type=self.sequence_type,
        )
        assert oligo_score_wo_exon == 2, "error: scoring for region wo exon incorrect."
        assert oligo_score_wo_exon > oligo_score_with_exon1, "error: exon score incorrect."
        assert oligo_score_with_exon2 == 0, "error: exon score incorrect."

    def test_isoform_consensus(self) -> None:
        oligo_score1 = self.isoform_consensus_scorer.apply(
            oligo_database=self.oligo_database,
            region_id="region_1",
            oligo_id="region_1::1",
            sequence_type=self.sequence_type,
        )
        assert oligo_score1 == 1, "error: scoring for isoform consensus incorrect."

    def test_Tm_scorer(self) -> None:
        oligo_score = self.Tm_scorer.apply(
            oligo_database=self.oligo_database,
            region_id="region_1",
            oligo_id="region_1::1",
            sequence_type=self.sequence_type,
        )
        assert oligo_score == 0, "error: scoring for Tm incorrect."

    def test_GC_scorer(self) -> None:
        oligo_score = self.GC_scorer.apply(
            oligo_database=self.oligo_database,
            region_id="region_1",
            oligo_id="region_1::1",
            sequence_type=self.sequence_type,
        )
        assert oligo_score == 0, "error: scoring for GC content incorrect."

    def test_Tm_norm_scorer(self) -> None:
        oligo_score = self.Tm_norm_scorer.apply(
            oligo_database=self.oligo_database,
            region_id="region_1",
            oligo_id="region_1::1",
            sequence_type=self.sequence_type,
        )
        assert abs(oligo_score - 0.496) < 1e-3, "error: scoring for Tm incorrect."

    def test_GC_norm_scorer(self) -> None:
        oligo_score = self.GC_norm_scorer.apply(
            oligo_database=self.oligo_database,
            region_id="region_1",
            oligo_id="region_1::1",
            sequence_type=self.sequence_type,
        )
        assert oligo_score == 0.25, "error: scoring for GC content incorrect."

    def test_uniform_distance_scorer(self) -> None:
        # Three oligos with pairwise distances 5, 15, 5 -> dist_opt=5; oligo 1 with set [2,3] has d_min=5 -> score 0
        matrix = csr_matrix([[0, 5, 15], [5, 0, 5], [15, 5, 0]])
        ids = ["region_1::1", "region_1::2", "region_1::3"]
        oligo_score = self.uniform_distance_scorer.apply(
            oligo_database=self.oligo_database,
            region_id="region_1",
            oligo_id="region_1::1",
            sequence_type=self.sequence_type,
            non_overlap_matrix=matrix,
            non_overlap_matrix_ids=ids,
            set_oligo_ids=["region_1::2", "region_1::3"],
            oligoset_size=3,
        )
        assert oligo_score == 0, "error: uniform distance scoring incorrect."

    def test_oligo_scoring(self) -> None:
        oligos_scoring = OligoScoring(
            scorers=[self.exon_scorer, self.isoform_consensus_scorer, self.Tm_scorer, self.GC_norm_scorer]
        )
        oligo_scores = oligos_scoring.apply(
            oligo_database=self.oligo_database,
            region_id="region_1",
            oligo_ids=["region_1::1", "region_1::2"],
            sequence_type=self.sequence_type,
        )

        assert oligo_scores.loc["region_1::1"] == 3.25, "error: wrong score computed for oligo region_1::1"
        assert oligo_scores.loc["region_1::2"] == 0.25, "error: wrong score computed for oligo region_1::2"


class TestSetScoring(unittest.TestCase):
    def setUp(self) -> None:
        self.score_max_sum = LowestSetScoring(ascending=True)
        self.score_ave_max = AverageSetScoring(ascending=True)
        self.oligo_set = Series(data=[0, 1, 8, 5, 2, 6, 7, 3])
        self.n_oligo_set = 5

    def test_max_sum(self) -> None:
        oligoset = self.score_max_sum.apply(self.oligo_set, self.n_oligo_set)
        assert oligoset[0] == [0, 1, 4, 7, 3], "Max scoring failed case"
        assert oligoset[1] == {"set_score_worst": 5, "set_score_sum": 11}, "Max scoring failed"

    def test_ave_max(self) -> None:
        oligoset = self.score_ave_max.apply(self.oligo_set, self.n_oligo_set)
        assert oligoset[0] == [0, 1, 4, 7, 3], "Average scoring failed"
        assert oligoset[1] == {"set_score_average": 2.2, "set_score_worst": 5}, "Average scoring failed"
