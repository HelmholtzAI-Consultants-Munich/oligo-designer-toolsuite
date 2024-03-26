import os
import shutil
import unittest
from abc import abstractmethod

from oligo_designer_toolsuite.database import OligoDatabase, ReferenceDatabase
from oligo_designer_toolsuite.oligo_specificity_filter import (
    BlastNFilter,
    BlastNSeedregionLigationsiteFilter,
    Bowtie2Filter,
    BowtieFilter,
    CrossHybridizationFilter,
    ExactMatchFilter,
    RemoveByDegreePolicy,
    RemoveByLargerRegionPolicy,
)

# Global Parameters
FILE_DATABASE_OLIGOS_EXACT_MATCH = "data/tests/databases/database_oligos_exactmatch.tsv"
FILE_DATABASE_OLIGOS_MATCH = "data/tests/databases/database_oligos_match.tsv"
FILE_DATABASE_OLIGOS_NOMATCH = "data/tests/databases/database_oligos_nomatch.tsv"
FILE_DATABASE_REFERENCE = "data/tests/databases/database_reference.fna"
FILE_DATABASE_OLIGOS_LIGATION_MATCH = "data/tests/databases/database_oligos_ligation_match.tsv"
FILE_DATABASE_OLIGOS_LIGATION_NOMATCH = "data/tests/databases/database_oligos_ligation_nomatch.tsv"

FILE_DATABASE_REFERENCE_LIGATION = "data/tests/databases/database_reference_ligation.fna"

FILE_DATABASE_OLIGOS_CROSSHYB = "data/tests/databases/database_oligos_crosshybridization.tsv"


class TestExactMatchFilter(unittest.TestCase):
    def setUp(self):
        self.tmp_path = os.path.join(os.getcwd(), "tmp_exact_match_outputs")
        self.filter = ExactMatchFilter()
        self.oligo_database = OligoDatabase(
            min_oligos_per_region=2,
            write_regions_with_insufficient_oligos=True,
            dir_output=os.path.join(self.tmp_path, "oligo_database_match"),
        )
        self.oligo_database.load_database(FILE_DATABASE_OLIGOS_EXACT_MATCH)

    def tearDown(self):
        shutil.rmtree(self.tmp_path)

    def test_exact_match_filter(self):
        sequence_type = "oligo"

        res = self.filter.apply(sequence_type, self.oligo_database, 2)

        assert (
            "WASH7P::2" not in res.database["WASH7P"].keys()
        ), "A matching oligo has not been filtered from exact matches!"
        assert (
            "AGRN::1" in res.database["AGRN"].keys()
        ), "A non-matching oligo has been filtered from exact mathces!"


class AlignmentFilterTestBase:
    def setUp(self):
        self.tmp_path = os.path.join(os.getcwd(), self.setup_tmp_path())
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
    def setup_tmp_path(self):
        pass

    @abstractmethod
    def setup_filter(self):
        pass

    def _setup_databases(self, database_file_match, database_file_nomatch, database_reference):
        self.oligo_database_match = OligoDatabase(
            min_oligos_per_region=2,
            write_regions_with_insufficient_oligos=True,
            dir_output=os.path.join(self.tmp_path, "oligo_database_match"),
        )
        self.oligo_database_match.load_database(database_file_match)

        self.oligo_database_nomatch = OligoDatabase(
            min_oligos_per_region=2,
            write_regions_with_insufficient_oligos=True,
            dir_output=os.path.join(self.tmp_path, "oligo_database_nomatch"),
        )
        self.oligo_database_nomatch.load_database(database_file_nomatch)

        self.reference_database = ReferenceDatabase(dir_output=self.tmp_path)

        self.reference_database.load_sequences_from_fasta(
            file_fasta=database_reference, database_overwrite=True
        )

    def test_filter_match(self):
        sequence_type = "target"

        res = self.filter.apply(sequence_type, self.oligo_database_match, 2, self.reference_database)

        assert "WASH7P::1" not in res.database["WASH7P"].keys(), "A matching oligo has not been filtered!"

    def test_filter_nomatch(self):
        sequence_type = "target"
        res = self.filter.apply(sequence_type, self.oligo_database_nomatch, 2, self.reference_database)

        assert "AGRN::1" in res.database["AGRN"].keys(), "A non matching oligo has been filtered by Blast!"


class TestBlastFilter(AlignmentFilterTestBase, unittest.TestCase):
    def setup_filter(self):
        blast_search_parameters = {
            "perc_identity": 80,
            "strand": "plus",
            "word_size": 10,
        }
        blast_hit_parameters = {"coverage": 50}

        return BlastNFilter(blast_search_parameters, blast_hit_parameters, dir_output=self.tmp_path)

    def setup_tmp_path(self):
        return "tmp_blast_outputs"


class TestBowtieFilter(AlignmentFilterTestBase, unittest.TestCase):
    def setup_filter(self):
        bowtie_search_parameters = {"-n": 3, "-l": 5}

        return BowtieFilter(bowtie_search_parameters, dir_output=self.tmp_path)

    def setup_tmp_path(self):
        return "__tmp_bowtie_outputs"


class TestBowtie2Filter(AlignmentFilterTestBase, unittest.TestCase):
    def setup_filter(self):
        bowtie2_search_parameters = {"-N": 0}

        return Bowtie2Filter(bowtie2_search_parameters, dir_output=self.tmp_path)

    def setup_tmp_path(self):
        return "tmp_bowtie2_outputs"


class TestBlastNSeedregionLigationsiteFilter(AlignmentFilterTestBase, unittest.TestCase):
    def setup_filter(self):
        blast_search_parameters = {
            "perc_identity": 80,
            "strand": "plus",
            "word_size": 10,
        }
        blast_hit_parameters = {"coverage": 50}
        seedregion_size = 10

        return BlastNSeedregionLigationsiteFilter(
            seedregion_size,
            blast_search_parameters,
            blast_hit_parameters,
            dir_output=self.tmp_path,
        )

    def setup_tmp_path(self):
        return "tests/tmp_blast_ligation_outputs"

    def setUp(self):
        super().setUp()
        self._setup_databases(
            database_file_match=FILE_DATABASE_OLIGOS_LIGATION_MATCH,
            database_file_nomatch=FILE_DATABASE_OLIGOS_LIGATION_NOMATCH,
            database_reference=FILE_DATABASE_REFERENCE_LIGATION,
        )


class TestCrossHybridizationFilter(unittest.TestCase):
    def setUp(self):
        self.tmp_path = os.path.join(os.getcwd(), "tmp_crosshybridization_outputs")
        os.makedirs(self.tmp_path, exist_ok=True)
        self.oligo_database_crosshyb = self._setup_database(FILE_DATABASE_OLIGOS_CROSSHYB)
        self.oligo_database_crosshyb_exactmatch = self._setup_database(FILE_DATABASE_OLIGOS_EXACT_MATCH)
        self.sequence_type = "oligo"

        # Blast parameters
        self.blast_search_parameters_crosshyb = {
            "perc_identity": 80,
            "strand": "minus",
            "word_size": 10,
        }
        self.blast_hit_parameters_crosshyb = {"coverage": 50}

        # Bowtie parameters
        self.bowtie_search_parameters_crosshyb = {"-n": 3, "-l": 5, "--nofw": ""}

        self.expected_oligos_larger_region = {
            "region_1": {
                "region_1::oligo_7",
                "region_1::oligo_5",
                "region_1::oligo_6",
                "region_1::oligo_8",
                "region_1::oligo_4",
            },
            "region_2": {
                "region_2::oligo_3",
                "region_2::oligo_2",
                "region_2::oligo_6",
                "region_2::oligo_5",
                "region_2::oligo_4",
            },
            "region_3": {
                "region_3::oligo_1",
                "region_3::oligo_4",
                "region_3::oligo_3",
                "region_3::oligo_2",
                "region_3::oligo_5",
            },
        }

        self.expected_oligos_degree = {
            "region_1": {
                "region_1::oligo_1",
                "region_1::oligo_4",
                "region_1::oligo_5",
                "region_1::oligo_6",
                "region_1::oligo_7",
            },
            "region_2": {
                "region_2::oligo_1",
                "region_2::oligo_2",
                "region_2::oligo_3",
                "region_2::oligo_4",
                "region_2::oligo_5",
                "region_2::oligo_6",
                "region_2::oligo_7",
            },
            "region_3": {
                "region_3::oligo_2",
                "region_3::oligo_3",
                "region_3::oligo_4",
                "region_3::oligo_5",
            },
        }

    def tearDown(self):
        shutil.rmtree(self.tmp_path)

    def _setup_database(self, file_database):
        oligos = OligoDatabase(
            min_oligos_per_region=2,
            write_regions_with_insufficient_oligos=True,
            dir_output=os.path.join(self.tmp_path, "oligo_database"),
        )
        oligos.load_database(file_database)
        return oligos

    def _apply_filter_and_assert(self, filter_instance, expected_oligos):
        res = filter_instance.apply(self.sequence_type, self.oligo_database_crosshyb, 2)
        filtered_oligos = {
            key: {key_2 for key_2 in list(res.database[key].keys())} for key in list(res.database.keys())
        }
        self.assertEqual(
            expected_oligos,
            filtered_oligos,
            f"The cross-hybridization filter didn't return the expected oligos. \n\nExpected:\n{expected_oligos}\n\nGot:\n{filtered_oligos}",
        )

    def test_crosshyb_filter_blast_larger_region_policy(self):
        filter_instance = BlastNFilter(
            self.blast_search_parameters_crosshyb,
            self.blast_hit_parameters_crosshyb,
            dir_output=os.path.join(self.tmp_path, "blast_larger_region"),
        )
        policy = RemoveByLargerRegionPolicy()
        cross_hyb_filter = CrossHybridizationFilter(policy, filter_instance, self.tmp_path)
        self._apply_filter_and_assert(cross_hyb_filter, self.expected_oligos_larger_region)

    def test_crosshyb_filter_blast_degree_policy(self):
        filter_instance = BlastNFilter(
            self.blast_search_parameters_crosshyb,
            self.blast_hit_parameters_crosshyb,
            dir_output=os.path.join(self.tmp_path, "blast_degree"),
        )
        policy = RemoveByDegreePolicy()
        cross_hyb_filter = CrossHybridizationFilter(policy, filter_instance, self.tmp_path)
        self._apply_filter_and_assert(cross_hyb_filter, self.expected_oligos_degree)

    def test_crosshyb_filter_bowtie_larger_region_policy(self):
        filter_instance = BowtieFilter(
            self.bowtie_search_parameters_crosshyb,
            dir_output=os.path.join(self.tmp_path, "bowtie_larger_region"),
        )
        policy = RemoveByLargerRegionPolicy()
        cross_hyb_filter = CrossHybridizationFilter(policy, filter_instance, self.tmp_path)
        self._apply_filter_and_assert(cross_hyb_filter, self.expected_oligos_larger_region)

    def test_crosshyb_filter_bowtie_degree_policy(self):
        filter_instance = BowtieFilter(
            self.bowtie_search_parameters_crosshyb,
            dir_output=os.path.join(self.tmp_path, "bowtie_degree"),
        )
        policy = RemoveByDegreePolicy()
        cross_hyb_filter = CrossHybridizationFilter(policy, filter_instance, self.tmp_path)
        self._apply_filter_and_assert(cross_hyb_filter, self.expected_oligos_degree)

    def test_crosshyb_filter_exactmatch_larger_region_policy(self):
        filter_instance = ExactMatchFilter()
        policy = RemoveByLargerRegionPolicy()

        cross_hyb_filter = CrossHybridizationFilter(policy, filter_instance, self.tmp_path)
        res = cross_hyb_filter.apply(self.sequence_type, self.oligo_database_crosshyb_exactmatch, 2)

        assert (
            "WASH7P::1" not in res.database["WASH7P"].keys()
        ), "A non matching oligo has been filtered by exact matches!"
        assert (
            "WASH7P::3" not in res.database["WASH7P"].keys()
        ), "A matching oligo has not been filtered by exact mathces!"
        assert (
            "AGRN::1" in res.database["AGRN"].keys()
        ), "A non matching oligo has been filtered by exact matches!"
