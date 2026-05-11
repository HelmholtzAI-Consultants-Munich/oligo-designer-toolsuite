############################################
# imports
############################################

import os
import shutil
import unittest

import pandas as pd
from Bio.SeqUtils import MeltingTemp as mt

from oligo_designer_toolsuite._exceptions import DatabaseError
from oligo_designer_toolsuite.database import OligoDatabase
from oligo_designer_toolsuite.oligo_efficiency_filter import (
    LowestSetScoring,
    NormalizedDeviationFromOptimalGCContentScorer,
    NormalizedDeviationFromOptimalTmScorer,
    OligoScoring,
)
from oligo_designer_toolsuite.oligo_property_calculator import (
    GCContentProperty,
    PropertyCalculator,
    TmNNProperty,
)
from oligo_designer_toolsuite.oligo_selection import (
    HomogeneousPropertyOligoSelection,
    IndependentSetsOligoSelection,
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


class TestIndependentSetsOligoSelection(unittest.TestCase):
    """Unified tests for IndependentSetsOligoSelection: validation, non-overlap matrix, apply (synthetic and file DB), edge cases."""

    def setUp(self) -> None:
        self.tmp_path = os.path.join(os.getcwd(), "tmp_oligo_selection")
        os.makedirs(self.tmp_path, exist_ok=True)
        self.oligo_database = OligoDatabase(
            min_oligos_per_region=2,
            write_regions_with_insufficient_oligos=True,
            dir_output=self.tmp_path,
        )
        self.oligo_database.load_database_from_table(
            FILE_DATABASE, database_overwrite=True, merge_databases_on_sequence_type="oligo"
        )

        Tm_scorer = NormalizedDeviationFromOptimalTmScorer(
            Tm_min=52,
            Tm_opt=60,
            Tm_max=67,
            Tm_parameters=TM_PARAMETERS,
            Tm_chem_correction_parameters=TM_PARAMETERS_CHEM_CORR,
            Tm_salt_correction_parameters=None,
            score_weight=1,
        )
        GC_scorer = NormalizedDeviationFromOptimalGCContentScorer(
            GC_content_min=40,
            GC_content_opt=50,
            GC_content_max=60,
            score_weight=1,
        )
        self.oligo_scoring = OligoScoring(scorers=[Tm_scorer, GC_scorer])
        self.set_scoring = LowestSetScoring(ascending=True)
        self.set_size_opt = 5
        self.set_size_min = 2
        self.oligoset_generator = IndependentSetsOligoSelection(
            oligos_scoring=self.oligo_scoring,
            set_scoring=self.set_scoring,
            set_size_opt=self.set_size_opt,
            set_size_min=self.set_size_min,
            distance_between_oligos=0,
            n_attempts_graph=100,
            n_attempts_clique_enum=100,
            diversification_fraction=0.3,
            jaccard_opt=1,
            jaccard_step=0,
        )

    def tearDown(self) -> None:
        if os.path.exists(self.tmp_path):
            shutil.rmtree(self.tmp_path)

    # --- Constructor validation ---

    def test_jaccard_opt_validation_rejects_invalid(self) -> None:
        """Constructor should raise ValueError for jaccard_opt outside (0, 1]."""
        with self.assertRaises(ValueError):
            IndependentSetsOligoSelection(
                oligos_scoring=self.oligo_scoring,
                set_scoring=self.set_scoring,
                set_size_opt=3,
                set_size_min=2,
                distance_between_oligos=0,
                n_attempts_graph=5,
                n_attempts_clique_enum=10,
                diversification_fraction=0.1,
                jaccard_opt=-0.1,
                jaccard_step=0.1,
            )
        with self.assertRaises(ValueError):
            IndependentSetsOligoSelection(
                oligos_scoring=self.oligo_scoring,
                set_scoring=self.set_scoring,
                set_size_opt=3,
                set_size_min=2,
                distance_between_oligos=0,
                n_attempts_graph=5,
                n_attempts_clique_enum=10,
                diversification_fraction=0.1,
                jaccard_opt=1.5,
                jaccard_step=0.1,
            )

    def test_jaccard_step_validation_rejects_invalid(self) -> None:
        """Constructor should raise ValueError for jaccard_step outside (0, 1]."""
        with self.assertRaises(ValueError):
            IndependentSetsOligoSelection(
                oligos_scoring=self.oligo_scoring,
                set_scoring=self.set_scoring,
                set_size_opt=3,
                set_size_min=2,
                distance_between_oligos=0,
                n_attempts_graph=5,
                n_attempts_clique_enum=10,
                diversification_fraction=0.1,
                jaccard_opt=0.5,
                jaccard_step=-0.1,
            )

    # --- Non-overlap matrix ---

    def test_nonoverlapping_matrix_overlapping_oligos(self) -> None:
        """Overlapping oligos should yield 0 (non-compatible) in the non-overlap matrix."""
        oligo_database = OligoDatabase(dir_output=self.tmp_path)
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
        self.assertTrue(
            true_matrix.equals(computed_matrix),
            "overlapping matrix for two overlapping oligos wrongly computed",
        )

    def test_nonoverlapping_matrix_for_nonoverlapping_oligos(self) -> None:
        """Non-overlapping oligos should yield gap distance in the non-overlap matrix."""
        oligo_database = OligoDatabase(dir_output=self.tmp_path)
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
        true_matrix = pd.DataFrame(data=[[0, 5], [5, 0]], columns=["A_0", "A_1"], index=["A_0", "A_1"])
        self.assertTrue(
            true_matrix.equals(computed_matrix),
            "overlapping matrix for two non-overlapping oligos wrongly computed",
        )

    def test_non_overlap_matrix_distance_boundary(self) -> None:
        """Gap > distance_between_oligos gives stored distance; gap <= gives 0."""
        oligo_database = OligoDatabase(dir_output=self.tmp_path)
        oligo_database.database = {
            "r1": {
                "a": {"oligo": "A" * 10, "start": [[1]], "end": [[10]]},
                "b": {"oligo": "T" * 10, "start": [[16]], "end": [[25]]},
            }
        }
        generator = IndependentSetsOligoSelection(
            oligos_scoring=self.oligo_scoring,
            set_scoring=self.set_scoring,
            set_size_opt=2,
            set_size_min=2,
            distance_between_oligos=5,
            n_attempts_graph=1,
            n_attempts_clique_enum=5,
            diversification_fraction=0.1,
            jaccard_opt=0.5,
            jaccard_step=0.1,
        )
        matrix, _ = generator._get_non_overlap_matrix(oligo_database=oligo_database, region_id="r1")
        self.assertEqual(matrix[0, 1], 6)
        self.assertEqual(matrix[1, 0], 6)

    # --- Apply on synthetic DB ---

    def test_non_overlapping_sets(self) -> None:
        """Apply on a small synthetic database and check output structure."""
        oligo_database = OligoDatabase(dir_output=self.tmp_path)
        oligo_database.database = {
            "region_1": {
                "A_0": {
                    "oligo": "GCATCTCACTAAGATGTCTGTATCTGCGTGTGCG",
                    "start": [[10], [50]],
                    "end": [[15], [55]],
                },
                "A_1": {"oligo": "AATTAGAAGCGTGTGCGCACATCCCGG", "start": [[20], [53]], "end": [[25], [58]]},
                "A_2": {
                    "oligo": "GCATCTCACTAAGATGTCTGTATCTGCGTGTGCGCCCCCACATCC",
                    "start": [[30], [60]],
                    "end": [[35], [65]],
                },
                "A_3": {
                    "oligo": "AAAGCGTGTGTTGTGTTGCGCCCCCACATCCCG",
                    "start": [[40], [70]],
                    "end": [[45], [75]],
                },
                "A_4": {"oligo": "AAAGCTGTTGCGCCCCCACATCC", "start": [[50], [80]], "end": [[55], [85]]},
            }
        }
        oligo_database.oligosets = {}

        result = self.oligoset_generator.apply(
            oligo_database=oligo_database,
            sequence_type="oligo",
            n_sets=10,
            n_jobs=1,
        )
        self.assertIn("region_1", result.oligosets, "No sets generated for region_1")
        computed_sets = result.oligosets["region_1"]
        self.assertGreaterEqual(len(computed_sets), 1, "Selection should return at least one set")
        self.assertIn("set_score_worst", computed_sets.columns)
        self.assertIn("set_score_sum", computed_sets.columns)
        oligo_cols = [c for c in computed_sets.columns if c.startswith("oligo_")]
        self.assertGreaterEqual(
            len(oligo_cols), self.set_size_min, "Each set should have at least set_size_min oligos"
        )

    def test_region_insufficient_oligos_does_not_crash(self) -> None:
        """Region with fewer oligos than set_size_min should not crash (empty oligosets)."""
        oligo_database = OligoDatabase(dir_output=self.tmp_path)
        oligo_database.database = {
            "small_region": {
                "o1": {"oligo": "A" * 30, "start": [[1]], "end": [[30]]},
                "o2": {"oligo": "T" * 30, "start": [[50]], "end": [[79]]},
            }
        }
        oligo_database.oligosets = {}
        generator = IndependentSetsOligoSelection(
            oligos_scoring=self.oligo_scoring,
            set_scoring=self.set_scoring,
            set_size_opt=5,
            set_size_min=3,
            distance_between_oligos=0,
            n_attempts_graph=2,
            n_attempts_clique_enum=5,
            diversification_fraction=0.1,
            jaccard_opt=0.5,
            jaccard_step=0.1,
        )
        result = generator.apply(
            oligo_database=oligo_database,
            sequence_type="oligo",
            n_sets=2,
            n_jobs=1,
        )
        self.assertIn("small_region", result.oligosets)
        self.assertTrue(result.oligosets["small_region"].empty or len(result.oligosets["small_region"]) == 0)

    def test_request_more_sets_than_candidates_returns_available(self) -> None:
        """Requesting n_sets larger than possible should return as many as found without error."""
        oligo_database = OligoDatabase(dir_output=self.tmp_path)
        oligo_database.database = {
            "tiny": {
                "x": {"oligo": "A" * 30, "start": [[1]], "end": [[30]]},
                "y": {"oligo": "T" * 30, "start": [[40]], "end": [[69]]},
                "z": {"oligo": "G" * 30, "start": [[80]], "end": [[109]]},
            }
        }
        oligo_database.oligosets = {}
        generator = IndependentSetsOligoSelection(
            oligos_scoring=self.oligo_scoring,
            set_scoring=self.set_scoring,
            set_size_opt=3,
            set_size_min=2,
            distance_between_oligos=0,
            n_attempts_graph=5,
            n_attempts_clique_enum=10,
            diversification_fraction=0.1,
            jaccard_opt=0.5,
            jaccard_step=0.1,
        )
        result = generator.apply(
            oligo_database=oligo_database,
            sequence_type="oligo",
            n_sets=100,
            n_jobs=1,
        )
        self.assertIn("tiny", result.oligosets)
        self.assertGreaterEqual(len(result.oligosets["tiny"]), 1)
        self.assertLessEqual(len(result.oligosets["tiny"]), 100)

    # --- Apply on file DB (full invariants) ---

    def test_oligoset_generation(self) -> None:
        """Run set generation on file DB and assert structural and semantic invariants (non-deterministic)."""
        oligos_database = self.oligoset_generator.apply(
            oligo_database=self.oligo_database,
            sequence_type="oligo",
            n_sets=100,
            n_jobs=1,
        )
        oligos_database.remove_regions_with_insufficient_oligos(pipeline_step="oligoset generation")

        # NOTE: in the new implementation the results are not deterministic,
        # so we cannot test based on "true sets" anymore

        expected_score_cols = ["set_score_worst", "set_score_sum"]
        for gene in oligos_database.oligosets.keys():
            computed = oligos_database.oligosets[gene]
            region_oligos = set(oligos_database.database[gene].keys())

            self.assertIn("oligoset_id", computed.columns, msg=f"{gene}: missing oligoset_id")
            for col in expected_score_cols:
                self.assertIn(col, computed.columns, msg=f"{gene}: missing {col}")

            oligo_cols = [c for c in computed.columns if c.startswith("oligo_")]
            self.assertGreaterEqual(len(oligo_cols), self.set_size_min, msg=f"{gene}: too few oligo columns")
            self.assertLessEqual(len(oligo_cols), self.set_size_opt, msg=f"{gene}: too many oligo columns")

            if len(computed) == 0:
                continue

            set_size = len(oligo_cols)
            for _, row in computed.iterrows():
                oligo_ids = [row[c] for c in oligo_cols]
                self.assertEqual(len(oligo_ids), len(set(oligo_ids)), msg=f"{gene}: duplicate oligo in set")
                for oid in oligo_ids:
                    self.assertIn(oid, region_oligos, msg=f"{gene}: oligo {oid} not in database")

            self.assertTrue(
                computed["set_score_worst"].between(0, 2).all(), msg=f"{gene}: set_score_worst out of range"
            )
            self.assertTrue(
                computed["set_score_sum"].between(0, set_size * 2).all(),
                msg=f"{gene}: set_score_sum out of range",
            )
            self.assertFalse(
                computed.duplicated(subset=oligo_cols).any(), msg=f"{gene}: duplicate oligo sets found"
            )


class TestHomogeneousPropertyOligoSelection(unittest.TestCase):
    """Unified tests for HomogeneousPropertyOligoSelection: edge cases and apply on file DB."""

    def setUp(self) -> None:
        self.tmp_path = os.path.join(os.getcwd(), "tmp_oligo_selection_homogeneous")
        os.makedirs(self.tmp_path, exist_ok=True)
        self.oligo_database = OligoDatabase(
            min_oligos_per_region=2,
            write_regions_with_insufficient_oligos=True,
            dir_output=self.tmp_path,
        )
        self.oligo_database.load_database_from_table(
            FILE_DATABASE, database_overwrite=True, merge_databases_on_sequence_type="oligo"
        )
        properties = [GCContentProperty(), TmNNProperty(Tm_parameters=TM_PARAMETERS)]
        calculator = PropertyCalculator(properties=properties)
        self.oligo_database = calculator.apply(
            oligo_database=self.oligo_database, sequence_type="oligo", n_jobs=1
        )
        self.oligoset_generator = HomogeneousPropertyOligoSelection(
            set_size=5,
            properties={"GC_content_oligo": 1, "TmNN_oligo": 1},
            n_combinations=1000,
        )

    def tearDown(self) -> None:
        if os.path.exists(self.tmp_path):
            shutil.rmtree(self.tmp_path)

    def test_missing_property_raises_database_error(self) -> None:
        """Region missing a required property should raise DatabaseError."""
        oligo_database = OligoDatabase(dir_output=self.tmp_path)
        oligo_database.database = {
            "r1": {
                "o1": {"oligo": "ACGT" * 10, "GC_content_oligo": 0.5},
            }
        }
        oligo_database.oligosets = {}
        generator = HomogeneousPropertyOligoSelection(
            set_size=1,
            properties={"GC_content_oligo": 1, "TmNN_oligo": 1},
            n_combinations=10,
        )
        with self.assertRaises(DatabaseError):
            generator.apply(
                oligo_database=oligo_database,
                sequence_type="oligo",
                n_sets=1,
                n_jobs=1,
            )

    def test_region_fewer_oligos_than_set_size_handled(self) -> None:
        """Region with fewer oligos than set_size should not crash (empty or no sets)."""
        oligo_database = OligoDatabase(dir_output=self.tmp_path)
        oligo_database.database = {
            "small": {
                "a": {"oligo": "A" * 20, "GC_content_oligo": 0.5, "TmNN_oligo": 60.0},
                "b": {"oligo": "T" * 20, "GC_content_oligo": 0.5, "TmNN_oligo": 60.0},
            }
        }
        oligo_database.oligosets = {}
        generator = HomogeneousPropertyOligoSelection(
            set_size=5,
            properties={"GC_content_oligo": 1, "TmNN_oligo": 1},
            n_combinations=10,
        )
        result = generator.apply(
            oligo_database=oligo_database,
            sequence_type="oligo",
            n_sets=1,
            n_jobs=1,
        )
        self.assertIn("small", result.oligosets)
        self.assertTrue(result.oligosets["small"].empty or len(result.oligosets["small"]) == 0)

    def test_oligoset_generation(self) -> None:
        """Apply on file DB and check regions, number of sets per region, and column format (non-deterministic)."""
        oligos_database = self.oligoset_generator.apply(
            oligo_database=self.oligo_database,
            sequence_type="oligo",
            n_sets=2,
            n_jobs=1,
        )
        gene_ids = {"PLEKHN1", "MIB2", "UBE2J2", "DVL1", "AGRN", "LOC112268402_1"}
        self.assertEqual(
            gene_ids,
            set(oligos_database.oligosets.keys()),
            "The calculated oligosets regions are not correct!",
        )
        for gene in oligos_database.oligosets.keys():
            self.assertEqual(
                len(oligos_database.oligosets[gene]), 2, f"The number of oligosets for {gene} is not correct!"
            )
            self.assertEqual(
                set(oligos_database.oligosets[gene].columns),
                {"oligoset_id", "oligo_0", "oligo_1", "oligo_2", "oligo_3", "oligo_4", "set_score"},
                f"The columns of the oligosets are not correct! got {set(oligos_database.oligosets[gene].columns)}",
            )
