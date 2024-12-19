############################################
# imports
############################################

import os
import shutil
import tempfile
import unittest
from abc import abstractmethod

import numpy as np
import pandas as pd
from Bio.Seq import Seq
from oligo_designer_toolsuite_ai_filters.api import APIBase

from oligo_designer_toolsuite.database import OligoDatabase, ReferenceDatabase
from oligo_designer_toolsuite.oligo_specificity_filter import (
    BlastNFilter,
    BlastNSeedregionLigationsiteFilter,
    Bowtie2Filter,
    BowtieFilter,
    CrossHybridizationFilter,
    ExactMatchFilter,
    HybridizationProbabilityFilter,
    RemoveAllPolicy,
    RemoveByDegreePolicy,
    RemoveByLargerRegionPolicy,
)

############################################
# Setup
############################################

# Global Parameters
FILE_DATABASE_OLIGOS_EXACT_MATCH = "tests/data/databases/database_oligos_tsv/database_oligos_exactmatch.tsv"
FILE_DATABASE_OLIGOS_MATCH = "tests/data/databases/database_oligos_tsv/database_oligos_match.tsv"
FILE_DATABASE_OLIGOS_NOMATCH = "tests/data/databases/database_oligos_tsv/database_oligos_nomatch.tsv"
FILE_DATABASE_REFERENCE = "tests/data/databases/database_reference/database_reference.fna"
FILE_DATABASE_OLIGOS_LIGATION_MATCH = (
    "tests/data/databases/database_oligos_tsv/database_oligos_ligation_match.tsv"
)
FILE_DATABASE_OLIGOS_LIGATION_NOMATCH = (
    "tests/data/databases/database_oligos_tsv/database_oligos_ligation_nomatch.tsv"
)

FILE_DATABASE_REFERENCE_LIGATION = "tests/data/databases/database_reference/database_reference_ligation.fna"

FILE_DATABASE_OLIGOS_CROSSHYB = (
    "tests/data/databases/database_oligos_tsv/database_oligos_crosshybridization.tsv"
)
SOLUTIONS_LARGER_REGION = [
    f"tests/data/databases/expected_results/solution_crosshyb_larger_region_{i}.tsv" for i in range(3)
]
SOLUTIONS_DEGREE = [
    f"tests/data/databases/expected_results/solution_crosshyb_degree_{i}.tsv" for i in range(8)
]

FILE_DATABASE_OLIGOS_AI = "tests/data/databases/database_oligos_tsv/database_oligos_ai.tsv"
FILE_TABLE_HITS_BLAST_AI = "tests/data/table_hits/table_hits_blast_ai.tsv"
FILE_TABLE_HITS_BOWTIE_AI = "tests/data/table_hits/table_hits_bowtie_ai.tsv"


############################################
# Tests
############################################
class TestExactMatchFilter(unittest.TestCase):
    def setUp(self):
        self.tmp_path = os.path.join(os.getcwd(), "tmp_output_exactmatch_filter")
        self.oligo_database = OligoDatabase(
            min_oligos_per_region=2,
            write_regions_with_insufficient_oligos=True,
            database_name="db_oligo_exactmatch_filters_match",
            dir_output=self.tmp_path,
        )
        self.oligo_database.load_database_from_table(
            FILE_DATABASE_OLIGOS_EXACT_MATCH, database_overwrite=True
        )

    def tearDown(self):
        shutil.rmtree(self.tmp_path)

    def test_exact_match_filter_no_policy(self):
        sequence_type = "oligo"
        policy = RemoveAllPolicy()
        filter = ExactMatchFilter(policy)
        res = filter.apply(
            oligo_database=self.oligo_database, reference_database=None, sequence_type=sequence_type, n_jobs=2
        )

        assert (
            "WASH7P::2" not in res.database["WASH7P"].keys()
        ), "A matching oligo has not been filtered from exact matches!"
        assert (
            "AGRN::1" not in res.database["AGRN"].keys()
        ), "A non-matching oligo has been filtered from exact mathces!"

    def test_exact_match_filter_policy(self):
        sequence_type = "oligo"
        policy = RemoveByLargerRegionPolicy()
        filter = ExactMatchFilter(policy)
        res = filter.apply(
            oligo_database=self.oligo_database, reference_database=None, sequence_type=sequence_type, n_jobs=2
        )

        assert (
            "WASH7P::2" not in res.database["WASH7P"].keys()
        ), "A matching oligo has not been filtered from exact matches!"
        assert (
            "AGRN::1" in res.database["AGRN"].keys()
        ), "A non-matching oligo has been filtered from exact mathces!"


class AlignmentFilterTestBase:
    def setUp(self):
        self.tmp_path = tempfile.mkdtemp()
        os.makedirs(self.tmp_path, exist_ok=True)
        self.filter = self.setup_filter()
        self._setup_databases(
            database_file_match=FILE_DATABASE_OLIGOS_MATCH,
            database_file_nomatch=FILE_DATABASE_OLIGOS_NOMATCH,
            database_reference=FILE_DATABASE_REFERENCE,
        )

    def tearDown(self):
        shutil.rmtree(self.tmp_path)

    @abstractmethod
    def setup_filter(self):
        pass

    def _setup_databases(self, database_file_match, database_file_nomatch, database_reference):
        self.oligo_database_match = OligoDatabase(
            min_oligos_per_region=2,
            write_regions_with_insufficient_oligos=True,
            database_name="db_oligo_alignment_filters_match",
            dir_output=self.tmp_path,
        )
        self.oligo_database_match.load_database_from_table(database_file_match, database_overwrite=True)

        self.oligo_database_nomatch = OligoDatabase(
            min_oligos_per_region=2,
            write_regions_with_insufficient_oligos=True,
            database_name="db_oligo_alignment_filters_nomatch",
            dir_output=self.tmp_path,
        )
        self.oligo_database_nomatch.load_database_from_table(database_file_nomatch, database_overwrite=True)

        self.reference_database = ReferenceDatabase(
            database_name="db_reference_alignment_filters", dir_output=self.tmp_path
        )

        self.reference_database.load_database_from_fasta(
            files_fasta=database_reference, database_overwrite=True
        )

    def test_filter_match(self):
        sequence_type = "target"

        res = self.filter.apply(
            sequence_type=sequence_type,
            oligo_database=self.oligo_database_match,
            reference_database=self.reference_database,
            n_jobs=2,
        )

        assert "WASH7P::1" not in res.database["WASH7P"].keys(), "A matching oligo has not been filtered!"

    def test_filter_nomatch(self):
        sequence_type = "target"
        res = self.filter.apply(
            sequence_type=sequence_type,
            oligo_database=self.oligo_database_nomatch,
            reference_database=self.reference_database,
            n_jobs=2,
        )

        assert "AGRN::1" in res.database["AGRN"].keys(), "A non matching oligo has been filtered by Blast!"


class TestBlastFilter(AlignmentFilterTestBase, unittest.TestCase):
    def setup_filter(self):
        blastn_search_parameters = {
            "perc_identity": 80,
            "strand": "plus",
            "word_size": 10,
        }
        hit_parameters = {"coverage": 50}

        return BlastNFilter(
            blastn_search_parameters, hit_parameters, filter_name="blast", dir_output=self.tmp_path
        )


