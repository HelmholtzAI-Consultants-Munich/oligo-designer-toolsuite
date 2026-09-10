############################################
# imports
############################################

import os
import shutil
import unittest
from abc import abstractmethod

from oligo_designer_toolsuite.database import OligoDatabase, ReferenceDatabase
from oligo_designer_toolsuite.oligo_specificity_filter import (
    AlignmentSpecificityFilter,
    BaseSpecificityFilter,
    BlastNFilter,
    BlastNSeedregionSiteFilter,
    Bowtie2Filter,
    BowtieFilter,
    CrossHybridizationFilter,
    ExactMatchFilter,
    RemoveAllFilterPolicy,
    RemoveByDegreeFilterPolicy,
    RemoveByLargerRegionFilterPolicy,
    VariantsFilter,
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

FILE_DATABASE_REFERENCE_VARIANTS = "tests/data/databases/database_reference/database_reference_variants.vcf"


############################################
# Tests
############################################
class TestExactMatchFilter(unittest.TestCase):
    def setUp(self) -> None:
        self.tmp_path = os.path.join(os.getcwd(), "tmp_output_exactmatch_filter")
        self.oligo_database = OligoDatabase(
            min_oligos_per_region=2,
            write_regions_with_insufficient_oligos=True,
            database_name="db_oligo_exactmatch_filters_match",
            dir_output=self.tmp_path,
        )
        self.oligo_database.set_database_sequence_types(["oligo", "target"])
        self.oligo_database.load_database_from_table(
            FILE_DATABASE_OLIGOS_EXACT_MATCH,
            database_overwrite=True,
            merge_databases_on_sequence_type="oligo",
        )

    def tearDown(self) -> None:
        shutil.rmtree(self.tmp_path)

    def test_exact_match_filter_no_policy(self) -> None:
        policy = RemoveAllFilterPolicy()
        filter = ExactMatchFilter(policy=policy)
        res = filter.apply(oligo_database=self.oligo_database, sequence_type="oligo", n_jobs=2)

        assert (
            "WASH7P::2" not in res.database["WASH7P"].keys()
        ), "A matching oligo has not been filtered from exact matches!"
        assert (
            "AGRN::1" not in res.database["AGRN"].keys()
        ), "A non-matching oligo has been filtered from exact mathces!"

    def test_exact_match_filter_policy(self) -> None:
        policy = RemoveByLargerRegionFilterPolicy()
        filter = ExactMatchFilter(policy=policy)
        res = filter.apply(oligo_database=self.oligo_database, sequence_type="oligo", n_jobs=2)

        assert (
            "WASH7P::2" not in res.database["WASH7P"].keys()
        ), "A matching oligo has not been filtered from exact matches!"
        assert (
            "AGRN::1" in res.database["AGRN"].keys()
        ), "A non-matching oligo has been filtered from exact mathces!"


class AlignmentFilterTestBase:
    def setUp(self) -> None:
        self.tmp_path = os.path.join(os.getcwd(), "tmp_output_alignment_filter")
        os.makedirs(self.tmp_path, exist_ok=True)
        self.filter = self._setup_filter()
        self._setup_databases(
            database_file_match=FILE_DATABASE_OLIGOS_MATCH,
            database_file_nomatch=FILE_DATABASE_OLIGOS_NOMATCH,
            database_reference=FILE_DATABASE_REFERENCE,
        )

    def tearDown(self) -> None:
        shutil.rmtree(self.tmp_path)

    @abstractmethod
    def _setup_filter(self) -> AlignmentSpecificityFilter:
        pass

    def _setup_databases(
        self, database_file_match: str, database_file_nomatch: str, database_reference: str
    ) -> None:
        self.oligo_database_match = OligoDatabase(
            min_oligos_per_region=2,
            write_regions_with_insufficient_oligos=True,
            database_name="db_oligo_alignment_filters_match",
            dir_output=self.tmp_path,
        )
        self.oligo_database_match.set_database_sequence_types(["oligo", "target"])
        self.oligo_database_match.load_database_from_table(
            database_file_match, database_overwrite=True, merge_databases_on_sequence_type="oligo"
        )

        self.oligo_database_nomatch = OligoDatabase(
            min_oligos_per_region=2,
            write_regions_with_insufficient_oligos=True,
            database_name="db_oligo_alignment_filters_nomatch",
            dir_output=self.tmp_path,
        )
        self.oligo_database_nomatch.set_database_sequence_types(["oligo", "target"])
        self.oligo_database_nomatch.load_database_from_table(
            database_file_nomatch, database_overwrite=True, merge_databases_on_sequence_type="oligo"
        )

        self.reference_database = ReferenceDatabase(
            database_name="db_reference_alignment_filters", dir_output=self.tmp_path
        )

        self.reference_database.load_database_from_file(
            files=database_reference, file_type="fasta", database_overwrite=True
        )

    def test_filter_match(self) -> None:
        self.filter.set_reference_database(self.reference_database)

        self.filter.remove_hits = False
        db_flag = self.filter.apply(
            oligo_database=self.oligo_database_match,
            sequence_type="target",
            n_jobs=2,
        )
        res_flag = db_flag.database["WASH7P"]["WASH7P::1"][self.filter.filter_name]
        assert res_flag is not None, "A matching oligo has not been flagged!"

        self.filter.remove_hits = True
        db_remove = self.filter.apply(
            oligo_database=self.oligo_database_match,
            sequence_type="target",
            n_jobs=2,
        )
        res_remove = db_remove.database["WASH7P"].keys()
        assert "WASH7P::1" not in res_remove, "A matching oligo has not been filtered!"

    def test_filter_nomatch(self) -> None:
        self.filter.set_reference_database(self.reference_database)

        self.filter.remove_hits = False
        db_flag = self.filter.apply(
            oligo_database=self.oligo_database_nomatch,
            sequence_type="target",
            n_jobs=2,
        )
        res_flag = db_flag.database["AGRN"]["AGRN::1"][self.filter.filter_name]
        assert res_flag is None, "A non matching oligo has been flagged!"

        self.filter.remove_hits = True
        db_remove = self.filter.apply(
            oligo_database=self.oligo_database_nomatch,
            sequence_type="target",
            n_jobs=2,
        )
        res_remove = db_remove.database["AGRN"].keys()
        assert "AGRN::1" in res_remove, "A non matching oligo has been filtered!"


class TestBlastFilter(AlignmentFilterTestBase, unittest.TestCase):
    def _setup_filter(self) -> BlastNFilter:
        blastn_search_parameters = {
            "-perc_identity": 80,
            "-strand": "plus",
            "-word_size": 10,
        }
        hit_parameters = {"coverage": 50}

        return BlastNFilter(
            search_parameters=blastn_search_parameters,
            hit_parameters=hit_parameters,
            filter_name="blast",
            dir_output=self.tmp_path,
        )


class TestBowtieFilter(AlignmentFilterTestBase, unittest.TestCase):
    def _setup_filter(self) -> BowtieFilter:
        bowtie_search_parameters = {"-n": 3, "-l": 5}

        return BowtieFilter(
            search_parameters=bowtie_search_parameters,
            filter_name="bowtie",
            dir_output=self.tmp_path,
        )


class TestBowtie2Filter(AlignmentFilterTestBase, unittest.TestCase):
    def _setup_filter(self) -> Bowtie2Filter:
        bowtie2_search_parameters = {"-N": 0}

        return Bowtie2Filter(
            search_parameters=bowtie2_search_parameters,
            filter_name="bowtie2",
            dir_output=self.tmp_path,
        )


class TestBlastNSeedregionSiteFilter(AlignmentFilterTestBase, unittest.TestCase):
    def _setup_filter(self) -> BlastNSeedregionSiteFilter:
        blastn_search_parameters = {
            "-perc_identity": 80,
            "-strand": "plus",
            "-word_size": 10,
        }
        hit_parameters = {"coverage": 50}
        seedregion_size = 10

        return BlastNSeedregionSiteFilter(
            seedregion_size=seedregion_size,
            seedregion_site_name="ligation_site",
            search_parameters=blastn_search_parameters,
            hit_parameters=hit_parameters,
            filter_name="blast_ligationsite",
            dir_output=self.tmp_path,
        )

    def setUp(self) -> None:
        super().setUp()
        self._setup_databases(
            database_file_match=FILE_DATABASE_OLIGOS_LIGATION_MATCH,
            database_file_nomatch=FILE_DATABASE_OLIGOS_LIGATION_NOMATCH,
            database_reference=FILE_DATABASE_REFERENCE_LIGATION,
        )


class TestCrossHybridizationFilter(unittest.TestCase):
    def setUp(self) -> None:
        self.tmp_path = os.path.join(os.getcwd(), "tmp_output_crosshybridization_filter")
        os.makedirs(self.tmp_path, exist_ok=True)
        self.oligo_database_crosshyb = self._setup_database(FILE_DATABASE_OLIGOS_CROSSHYB)
        self.oligo_database_crosshyb_exactmatch = self._setup_database(FILE_DATABASE_OLIGOS_EXACT_MATCH)

        # Blast parameters
        self.blastn_search_parameters_crosshyb = {
            "-perc_identity": 80,
            "-strand": "minus",
            "-word_size": 10,
        }
        self.hit_parameters_crosshyb = {"coverage": 50}

        # Bowtie parameters
        self.bowtie_search_parameters_crosshyb = {"-n": 3, "-l": 5, "--nofw": ""}

        self.expected_oligos_larger_region = []
        for i, solution_file in enumerate(SOLUTIONS_LARGER_REGION):
            solution = OligoDatabase(
                min_oligos_per_region=2,
                write_regions_with_insufficient_oligos=True,
                database_name=f"db_oligo_crosshybridization_filters_solution_larger_region_{i}",
                dir_output=self.tmp_path,
            )
            solution.set_database_sequence_types(["oligo", "target"])
            solution.load_database_from_table(
                solution_file, database_overwrite=True, merge_databases_on_sequence_type="oligo"
            )
            self.expected_oligos_larger_region.append(solution.database)

        self.expected_oligos_degree = []
        for i, solution_file in enumerate(SOLUTIONS_DEGREE):
            solution = OligoDatabase(
                min_oligos_per_region=2,
                write_regions_with_insufficient_oligos=True,
                database_name=f"db_oligo_crosshybridization_filters_solution_degree_{i}",
                dir_output=self.tmp_path,
            )
            solution.set_database_sequence_types(["oligo", "target"])
            solution.load_database_from_table(
                solution_file, database_overwrite=True, merge_databases_on_sequence_type="oligo"
            )
            self.expected_oligos_degree.append(solution.database)

    def tearDown(self) -> None:
        shutil.rmtree(self.tmp_path)

    def _setup_database(self, file_database: str) -> OligoDatabase:
        oligos = OligoDatabase(
            min_oligos_per_region=2,
            write_regions_with_insufficient_oligos=True,
            database_name="db_oligo_crosshybridization_filters",
            dir_output=self.tmp_path,
        )
        oligos.set_database_sequence_types(["oligo", "target"])
        oligos.load_database_from_table(
            file_database, database_overwrite=True, merge_databases_on_sequence_type="oligo"
        )
        return oligos

    def _apply_filter_and_assert(
        self, filter_instance: BaseSpecificityFilter, expected_oligos: list[OligoDatabase]
    ) -> None:
        res = filter_instance.apply(
            oligo_database=self.oligo_database_crosshyb,
            sequence_type="oligo",
            n_jobs=2,
        )
        assert (
            res.database in expected_oligos
        ), f"The cross-hybridization filter didn't return the expected oligos."

    def test_crosshyb_filter_blast_larger_region_policy(self) -> None:
        filter_instance = BlastNFilter(
            search_parameters=self.blastn_search_parameters_crosshyb,
            hit_parameters=self.hit_parameters_crosshyb,
            filter_name="crosshybridization_blast_larger_region",
            dir_output=self.tmp_path,
        )
        policy = RemoveByLargerRegionFilterPolicy()
        cross_hyb_filter = CrossHybridizationFilter(
            policy=policy,
            alignment_method=filter_instance,
            filter_name="crosshybridization_blast_larger_region",
            dir_output=self.tmp_path,
        )
        self._apply_filter_and_assert(cross_hyb_filter, self.expected_oligos_larger_region)

    def test_crosshyb_filter_blast_degree_policy(self) -> None:
        filter_instance = BlastNFilter(
            search_parameters=self.blastn_search_parameters_crosshyb,
            hit_parameters=self.hit_parameters_crosshyb,
            filter_name="crosshybridization_blast_degree",
            dir_output=self.tmp_path,
        )
        policy = RemoveByDegreeFilterPolicy()
        cross_hyb_filter = CrossHybridizationFilter(
            policy=policy,
            alignment_method=filter_instance,
            filter_name="crosshybridization_blast_degree",
            dir_output=self.tmp_path,
        )
        self._apply_filter_and_assert(cross_hyb_filter, self.expected_oligos_degree)

    def test_crosshyb_filter_bowtie_larger_region_policy(self) -> None:
        filter_instance = BowtieFilter(
            search_parameters=self.bowtie_search_parameters_crosshyb,
            filter_name="crosshybridization_bowtie_larger_region",
            dir_output=self.tmp_path,
        )
        policy = RemoveByLargerRegionFilterPolicy()
        cross_hyb_filter = CrossHybridizationFilter(
            policy=policy,
            alignment_method=filter_instance,
            filter_name="crosshybridization_bowtie_larger_region",
            dir_output=self.tmp_path,
        )
        self._apply_filter_and_assert(cross_hyb_filter, self.expected_oligos_larger_region)

    def test_crosshyb_filter_bowtie_degree_policy(self) -> None:
        filter_instance = BowtieFilter(
            search_parameters=self.bowtie_search_parameters_crosshyb,
            filter_name="crosshybridization_bowtie_degree",
            dir_output=self.tmp_path,
        )
        policy = RemoveByDegreeFilterPolicy()
        cross_hyb_filter = CrossHybridizationFilter(
            policy=policy,
            alignment_method=filter_instance,
            filter_name="crosshybridization_bowtie_degree",
            dir_output=self.tmp_path,
        )
        self._apply_filter_and_assert(cross_hyb_filter, self.expected_oligos_degree)


class TestVariantsFilter(unittest.TestCase):
    def setUp(self) -> None:
        self.tmp_path = os.path.join(os.getcwd(), "tmp_output_variants_filter")
        os.makedirs(self.tmp_path, exist_ok=True)
        self.filter = self._setup_filter()
        self._setup_databases(
            database_file_match=FILE_DATABASE_OLIGOS_MATCH,
            database_file_nomatch=FILE_DATABASE_OLIGOS_NOMATCH,
            database_reference=FILE_DATABASE_REFERENCE_VARIANTS,
        )

    def tearDown(self) -> None:
        shutil.rmtree(self.tmp_path)

    def _setup_filter(self) -> VariantsFilter:
        return VariantsFilter(filter_name="variants", dir_output=self.tmp_path)

    def _setup_databases(
        self, database_file_match: str, database_file_nomatch: str, database_reference: str
    ) -> None:
        self.oligo_database_match = OligoDatabase(
            min_oligos_per_region=2,
            write_regions_with_insufficient_oligos=True,
            database_name="db_oligo_variants_filters_match",
            dir_output=self.tmp_path,
        )
        self.oligo_database_match.set_database_sequence_types(["oligo", "target"])
        self.oligo_database_match.load_database_from_table(
            database_file_match, database_overwrite=True, merge_databases_on_sequence_type="oligo"
        )

        self.oligo_database_nomatch = OligoDatabase(
            min_oligos_per_region=2,
            write_regions_with_insufficient_oligos=True,
            database_name="db_oligo_variants_filters_nomatch",
            dir_output=self.tmp_path,
        )
        self.oligo_database_nomatch.set_database_sequence_types(["oligo", "target"])
        self.oligo_database_nomatch.load_database_from_table(
            database_file_nomatch, database_overwrite=True, merge_databases_on_sequence_type="oligo"
        )

        self.reference_database = ReferenceDatabase(
            database_name="db_reference_variants_filters", dir_output=self.tmp_path
        )

        self.reference_database.load_database_from_file(
            files=database_reference, file_type="vcf", database_overwrite=True
        )

    def test_filter_match(self) -> None:

        self.filter.set_reference_database(self.reference_database)

        self.filter.remove_hits = False
        db_flag = self.filter.apply(
            oligo_database=self.oligo_database_match,
            n_jobs=2,
        )
        res_flag = db_flag.database["WASH7P"]["WASH7P::1"][self.filter.filter_name]
        assert res_flag is not None, "A matching oligo has not been flagged!"

        self.filter.remove_hits = True
        db_remove = self.filter.apply(
            oligo_database=self.oligo_database_match,
            n_jobs=2,
        )
        res_remove = db_remove.database["WASH7P"].keys()
        assert "WASH7P::1" not in res_remove, "A matching oligo has not been filtered!"

    def test_filter_nomatch(self) -> None:
        self.filter.set_reference_database(self.reference_database)

        self.filter.remove_hits = False
        db_flag = self.filter.apply(
            oligo_database=self.oligo_database_nomatch,
            n_jobs=2,
        )
        res_flag = db_flag.database["AGRN"]["AGRN::1"][self.filter.filter_name]
        assert res_flag is None, "A non matching oligo has been flagged!"

        self.filter.remove_hits = True
        db_remove = self.filter.apply(
            oligo_database=self.oligo_database_nomatch,
            n_jobs=2,
        )
        res_remove = db_remove.database["AGRN"].keys()
        assert "AGRN::1" in res_remove, "A non matching oligo has been filtered!"
