############################################
# imports
############################################

import os
import shutil
import unittest

import pandas as pd
from Bio.SeqUtils import MeltingTemp as mt
from scipy.sparse import csr_matrix

from oligo_designer_toolsuite.database import OligoDatabase, OligoAttributes
from oligo_designer_toolsuite.oligo_efficiency_filter import (
    LowestSetScoring,
    WeightedTmGCOligoScoring,
)
from oligo_designer_toolsuite.oligo_selection import (
    GraphBasedSelectionPolicy,
    GreedySelectionPolicy,
    HomogeneousPropertyOligoSetGenerator,
    OligosetGeneratorIndependentSet,
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


class TestOligosetGeneratorIndependentSet(unittest.TestCase):
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
        self.set_scoring = LowestSetScoring(ascending=True)
        self.selection_policy = GraphBasedSelectionPolicy(
            set_scoring=self.set_scoring,
            pre_filter=False,
            n_attempts=100000,
            heuristic=True,
            heuristic_n_attempts=100,
        )

        self.oligoset_generator = OligosetGeneratorIndependentSet(
            selection_policy=self.selection_policy,
            oligos_scoring=self.oligo_scoring,
            set_scoring=self.set_scoring,
            max_oligos=5000,
        )

        self.set_size_opt = 5
        self.set_size_min = 2

    def tearDown(self) -> None:
        shutil.rmtree(self.tmp_path)

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
        computed_matrix, computed_matrix_ids = self.oligoset_generator._get_non_overlap_matrix(
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
        computed_matrix, computed_matrix_ids = self.oligoset_generator._get_non_overlap_matrix(
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

    def test_oligoset_generation(self):
        oligos_database = self.oligoset_generator.apply(
            oligo_database=self.oligo_database,
            sequence_type="oligo",
            set_size_opt=self.set_size_opt,
            set_size_min=self.set_size_min,
            n_sets=100,
            n_jobs=1,
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
                "set_score_worst",
                "set_score_sum",
            ],
        )

        computed_sets = self.selection_policy.apply(
            set_size_opt=self.set_size_opt,
            set_size_min=self.set_size_min,
            n_sets=10,
            oligos_scores=oligo_scores,
            non_overlap_matrix=overlapping_matrix,
            non_overlap_matrix_ids=index,
        )
        computed_sets["set_score_worst"] = computed_sets["set_score_worst"].round(2)
        computed_sets["set_score_sum"] = computed_sets["set_score_sum"].round(2)

        assert true_sets.equals(computed_sets), "Sets are not computed correctly"


class TestHomogeneousPropertyOligoSetGenerator(unittest.TestCase):
    def setUp(self):
        self.tmp_path = os.path.join(os.getcwd(), "tmp_oligo_selection")
        self.oligo_database = OligoDatabase(
            min_oligos_per_region=2,
            write_regions_with_insufficient_oligos=True,
            dir_output=self.tmp_path,
        )
        self.oligo_database.load_database_from_table(FILE_DATABASE, database_overwrite=True)

        oligo_attributes = OligoAttributes()
        oligo_database = oligo_attributes.calculate_GC_content(self.oligo_database)
        oligo_database = oligo_attributes.calculate_TmNN(oligo_database, TM_PARAMETERS)

        self.oligoset_generator = HomogeneousPropertyOligoSetGenerator(
            set_size=5,
            properties={"GC_content": 1, "TmNN": 1},
        )

    def tearDown(self) -> None:
        shutil.rmtree(self.tmp_path)

    def test_oligoset_generation(self):
        # ⚠️ The function is not stable and the results are not always the same (depends on n_combinations)
        # As a result, we only test the format
        oligos_database = self.oligoset_generator.apply(
            self.oligo_database, n_sets=2, n_combinations=1000, n_jobs=1
        )
        
        gene_ids = {"PLEKHN1", "MIB2", "UBE2J2", "DVL1", "AGRN", "LOC112268402_1"}
        assert gene_ids == set(oligos_database.oligosets.keys()), "The calculated oligosets regions are not correct!"
        for gene in oligos_database.oligosets.keys():
            assert len(oligos_database.oligosets[gene]) == 2, "The number of oligosets is not correct!"
            assert set(oligos_database.oligosets[gene].columns) == {
                "oligoset_id",
                "oligo_0",
                "oligo_1",
                "oligo_2",
                "oligo_3",
                "oligo_4",
                "set_score",
            }, f"The columns of the oligosets are not correct! got {set(oligos_database.oligosets[gene].columns)}"


class TestOligoSelectionPolicy(unittest.TestCase):
    def setUp(self):
        self.tmp_path = os.path.join(os.getcwd(), "tmp_oligo_selection")
        self.region_id = "AGRN" # We only test for one region (this region is small enough for the purposes of the tests)

        self.oligo_database = OligoDatabase(
            min_oligos_per_region=2,
            write_regions_with_insufficient_oligos=True,
            dir_output=self.tmp_path,
        )
        self.oligo_database.load_database_from_table(FILE_DATABASE, database_overwrite=True)
        self.set_scoring = LowestSetScoring(ascending=True)

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
        self.oligo_database, self.oligos_scores = self.oligo_scoring.apply(
            oligo_database=self.oligo_database,
            region_id=self.region_id,
            sequence_type="oligo",
        )

        # sort oligos by score
        self.oligos_scores.sort_values(ascending=True, inplace=True)

        oligo_generator = OligosetGeneratorIndependentSet(selection_policy=None, oligos_scoring = self.oligo_scoring, set_scoring=self.set_scoring)
        # create the overlapping matrix
        self.non_overlap_matrix, self.non_overlap_matrix_ids = oligo_generator._get_non_overlap_matrix(
            oligo_database=self.oligo_database, region_id=self.region_id
        )

    def tearDown(self) -> None:
        shutil.rmtree(self.tmp_path)


    def test_graph_based_selection_policy(self):
        selection_policy = GraphBasedSelectionPolicy(
            set_scoring=self.set_scoring, 
            pre_filter=True
        )
        oligosets = selection_policy.apply(
            oligos_scores=self.oligos_scores,
            non_overlap_matrix=self.non_overlap_matrix,
            non_overlap_matrix_ids=self.non_overlap_matrix_ids,
            set_size_opt=5,
            set_size_min=3,
            n_sets=2,
        ).round(3)

        true_oligosets =  pd.DataFrame(
            {
                "oligoset_id": [0, 1],
                "oligo_0": ["AGRN_pid258", "AGRN_pid258"],
                "oligo_1": ["AGRN_pid77", "AGRN_pid77"],
                "oligo_2": ["AGRN_pid285", "AGRN_pid288"],
                "oligo_3": ["AGRN_pid248", "AGRN_pid248"],
                "set_score_worst": [1.017, 1.017],
                "set_score_sum": [2.314, 2.314],
            }
        )

        assert true_oligosets.equals(oligosets), "The oligosets are not computed correctly!"

    def test_greedy_selection_policy(self):
        selection_policy = GreedySelectionPolicy(
            set_scoring=self.set_scoring, 
            score_criteria=self.set_scoring.score_1, 
            pre_filter=True
        )
        oligosets = selection_policy.apply(
            oligos_scores=self.oligos_scores,
            non_overlap_matrix=self.non_overlap_matrix,
            non_overlap_matrix_ids=self.non_overlap_matrix_ids,
            set_size_opt=5,
            set_size_min=3,
            n_sets=2,
        ).round(3)
        
        true_oligosets_1 =  pd.DataFrame(
            {
                "oligoset_id": [0, 1],
                "oligo_0": ["AGRN_pid258", "AGRN_pid261"],
                "oligo_1": ["AGRN_pid285", "AGRN_pid288"],
                "oligo_2": ["AGRN_pid77", "AGRN_pid77"],
                "oligo_3": ["AGRN_pid248", "AGRN_pid248"],
                "set_score_worst": [1.017, 1.017],
                "set_score_sum": [2.314, 2.428],
            }
        )
        
        true_oligosets_2 =  pd.DataFrame(
            {
                "oligoset_id": [0, 1],
                "oligo_0": ["AGRN_pid258", "AGRN_pid261"],
                "oligo_1": ["AGRN_pid288", "AGRN_pid285"],
                "oligo_2": ["AGRN_pid77", "AGRN_pid77"],
                "oligo_3": ["AGRN_pid248", "AGRN_pid248"],
                "set_score_worst": [1.017, 1.017],
                "set_score_sum": [2.314, 2.428],
            }
        )


        assert true_oligosets_1.equals(oligosets) or true_oligosets_2.equals(oligosets), "The oligosets are not computed correctly!"


