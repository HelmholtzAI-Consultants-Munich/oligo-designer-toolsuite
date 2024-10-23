############################################
# imports
############################################

import os
import shutil
import unittest

import pandas as pd

from scipy.sparse import csr_matrix
from Bio.SeqUtils import MeltingTemp as mt

from oligo_designer_toolsuite.database import OligoDatabase
from oligo_designer_toolsuite.oligo_efficiency_filter import (
    LowestSetScoring,
    WeightedTmGCOligoScoring,
)
from oligo_designer_toolsuite.oligo_selection import (
    OligosetGeneratorIndependentSet,
    heuristic_selection_independent_set,
)

############################################
# Global Parameters
############################################

FILE_DATABASE = "tests/data/oligo_selection/oligos_info.tsv"

TM_PARAMETERS = {
    "check": True,  # default
    "strict": True,  # default
    "c_seq": None,  # default
    "shift": 0,  # default
    "nn_table": getattr(mt, "DNA_NN3"),
    "tmm_table": getattr(mt, "DNA_TMM1"),
    "imm_table": getattr(mt, "DNA_IMM1"),
    "de_table": getattr(mt, "DNA_DE1"),
    "dnac1": 50,  # [nM]
    "dnac2": 0,  # [nM]
    "selfcomp": False,  # default
    "saltcorr": 7,  # Owczarzy et al. (2008)
    "Na": 1.25,  # [mM]
    "K": 75,  # [mM]
    "Tris": 20,  # [mM]
    "Mg": 10,  # [mM]
    "dNTPs": 0,  # [mM] default
}

TM_PARAMETERS_CHEM_CORR = {
    "DMSO": 0,  # default
    "fmd": 20,
    "DMSOfactor": 0.75,  # default
    "fmdfactor": 0.65,  # default
    "fmdmethod": 1,  # default
    "GC": None,  # default
}


############################################
# Tests
############################################


class TestOverlapMatrix(unittest.TestCase):
    def setUp(self):
        oligo_scoring = WeightedTmGCOligoScoring(
            Tm_min=52,
            Tm_opt=60,
            Tm_max=67,
            GC_content_min=40,
            GC_content_opt=50,
            GC_content_max=60,
            Tm_parameters=TM_PARAMETERS,
            Tm_chem_correction_parameters=TM_PARAMETERS_CHEM_CORR,
        )
        set_scoring = LowestSetScoring(ascending=True)
        self.oligoset_generator = OligosetGeneratorIndependentSet(
            opt_oligoset_size=5,
            min_oligoset_size=2,
            oligos_scoring=oligo_scoring,
            set_scoring=set_scoring,
            heuristic_selection=None,
        )

    def test_nonoverlapping_matrix_ovelapping_oligos(self):
        # check the overlapping matrix is created correctly with 2 oligos given in input
        # generate sinthetic oligos that overlap
        oligo_database = OligoDatabase()
        oligo_database.database = {
            "region_1": {
                "A_0": {"start": [[10], [50]], "end": [[15], [55]]},
                "A_1": {"start": [[20], [53]], "end": [[25], [58]]},
            }
        }
        computed_matrix, computed_matrix_ids = self.oligoset_generator._get_overlapping_matrix(
            oligo_database=oligo_database, region_id="region_1"
        )
        computed_matrix = pd.DataFrame(
            data=computed_matrix.toarray(),
            columns=computed_matrix_ids,
            index=computed_matrix_ids,
            dtype=int,
        )

        true_matrix = pd.DataFrame(data=[[0, 0], [0, 0]], columns=["A_0", "A_1"], index=["A_0", "A_1"])
        assert true_matrix.equals(
            computed_matrix
        ), "overlapping matrix for two overlapping oligos wrongly computed"

    def test_nonoverlapping_matrix_for_nonovelapping_oligos(self):
        # generate sinthetic oligos that overlap
        oligo_database = OligoDatabase()
        oligo_database.database = {
            "region_1": {
                "A_0": {"start": [[10], [50]], "end": [[15], [55]]},
                "A_1": {"start": [[20], [35]], "end": [[25], [40]]},
            }
        }
        computed_matrix, computed_matrix_ids = self.oligoset_generator._get_overlapping_matrix(
            oligo_database=oligo_database, region_id="region_1"
        )
        computed_matrix = pd.DataFrame(
            data=computed_matrix.toarray(),
            columns=computed_matrix_ids,
            index=computed_matrix_ids,
            dtype=int,
        )

        true_matrix = pd.DataFrame(data=[[0, 1], [1, 0]], columns=["A_0", "A_1"], index=["A_0", "A_1"])
        assert true_matrix.equals(
            computed_matrix
        ), "overlapping matrix for two non-overlapping oligos wrongly computed"


class TestOligoScoring(unittest.TestCase):
    def setUp(self):
        self.tmp_path = os.path.join(os.getcwd(), "tmp_oligo_selection")
        self.oligo_database = OligoDatabase(
            min_oligos_per_region=2,
            write_regions_with_insufficient_oligos=True,
            dir_output=self.tmp_path,
        )
        self.oligo_database.load_database_from_table(FILE_DATABASE, database_overwrite=True)

        self.oligo_scoring = WeightedTmGCOligoScoring(
            Tm_min=52,
            Tm_opt=60,
            Tm_max=67,
            GC_content_min=40,
            GC_content_opt=50,
            GC_content_max=60,
            Tm_parameters=TM_PARAMETERS,
            Tm_chem_correction_parameters=TM_PARAMETERS_CHEM_CORR,
        )
        set_scoring = LowestSetScoring(ascending=True)
        self.oligoset_generator = OligosetGeneratorIndependentSet(
            opt_oligoset_size=5,
            min_oligoset_size=2,
            oligos_scoring=self.oligo_scoring,
            set_scoring=set_scoring,
            heuristic_selection=heuristic_selection_independent_set,
        )

    def tearDown(self) -> None:
        shutil.rmtree(self.tmp_path)

    def test_oligoset_generation(self):
        oligos_database = self.oligoset_generator.apply(
            oligo_database=self.oligo_database, sequence_type="oligo", n_sets=100, n_jobs=1
        )
        for gene in oligos_database.oligosets.keys():
            computed_sets = oligos_database.oligosets[gene]
            computed_sets.drop(columns=["oligoset_id"], inplace=True)
            true_sets = pd.read_csv(
                filepath_or_buffer=f"tests/data/oligo_selection/ranked_oligosets_{gene}.txt",
                sep="\t",
                index_col=0,
                float_precision="round_trip",
            )

            # additional sorting to guarantee the order is the same
            true_sets.sort_values(by=list(true_sets.columns), inplace=True)
            true_sets.reset_index(inplace=True, drop=True)

            computed_sets.sort_values(by=list(computed_sets.columns), inplace=True)
            computed_sets.reset_index(inplace=True, drop=True)

            assert true_sets.equals(computed_sets), f"Sets for {gene} are not computed correctly!"

    def test_non_overlapping_sets(self):
        oligo_database = OligoDatabase()
        oligo_database.database = {
            "region_1": {
                "A_0": {"oligo": "GCATCTCACTAAGATGTCTGTATCTGCGTGTGCG"},
                "A_1": {"oligo": "AATTAGAAGCGTGTGCGCACATCCCGG"},
                "A_2": {"oligo": "GCATCTCACTAAGATGTCTGTATCTGCGTGTGCGCCCCCACATCC"},
                "A_3": {"oligo": "AAAGCGTGTGTTGTGTTGCGCCCCCACATCCCG"},
                "A_4": {"oligo": "AAAGCTGTTGCGCCCCCACATCC"},
            }
        }
        oligos, oligo_scores = self.oligo_scoring.apply(
            oligo_database, region_id="region_1", sequence_type="oligo"
        )
        index = ["A_0", "A_1", "A_2", "A_3", "A_4"]
        data = [
            [0, 1, 1, 0, 0],
            [1, 0, 1, 0, 0],
            [1, 1, 0, 1, 1],
            [0, 0, 1, 0, 1],
            [0, 0, 1, 1, 0],
        ]
        overlapping_matrix = csr_matrix(data)
        data = [[0, "A_0", "A_1", "A_2", 1.59, 2.36], [1, "A_4", "A_2", "A_3", 2.15, 4.93]]
        true_sets = pd.DataFrame(
            data=data,
            columns=[
                "oligoset_id",
                "oligo_0",
                "oligo_1",
                "oligo_2",
                "set_score_lowest",
                "set_score_sum",
            ],
        )

        computed_sets = self.oligoset_generator._get_non_overlapping_sets(
            overlapping_matrix=overlapping_matrix,
            overlapping_matrix_ids=index,
            oligos_scores=oligo_scores,
            n_sets=100,
        )
        computed_sets["set_score_lowest"] = computed_sets["set_score_lowest"].round(2)
        computed_sets["set_score_sum"] = computed_sets["set_score_sum"].round(2)

        assert true_sets.equals(computed_sets), "Sets are not computed correctly"
