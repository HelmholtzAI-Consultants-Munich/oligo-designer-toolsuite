############################################
# imports
############################################

import os
import shutil
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
        self.tmp_path = os.path.join(os.getcwd(), "tmp_output_alignment_filter")
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


class TestBowtieFilter(AlignmentFilterTestBase, unittest.TestCase):
    def setup_filter(self):
        bowtie_search_parameters = {"-n": 3, "-l": 5}

        return BowtieFilter(bowtie_search_parameters, filter_name="bowtie", dir_output=self.tmp_path)


class TestBowtie2Filter(AlignmentFilterTestBase, unittest.TestCase):
    def setup_filter(self):
        bowtie2_search_parameters = {"-N": 0}

        return Bowtie2Filter(bowtie2_search_parameters, filter_name="bowtie2", dir_output=self.tmp_path)


class TestBlastNSeedregionLigationsiteFilter(AlignmentFilterTestBase, unittest.TestCase):
    def setup_filter(self):
        blastn_search_parameters = {
            "perc_identity": 80,
            "strand": "plus",
            "word_size": 10,
        }
        hit_parameters = {"coverage": 50}
        seedregion_size = 10

        return BlastNSeedregionLigationsiteFilter(
            seedregion_size,
            blastn_search_parameters,
            hit_parameters,
            filter_name="blast_ligationsite",
            dir_output=self.tmp_path,
        )

    def setUp(self):
        super().setUp()
        self._setup_databases(
            database_file_match=FILE_DATABASE_OLIGOS_LIGATION_MATCH,
            database_file_nomatch=FILE_DATABASE_OLIGOS_LIGATION_NOMATCH,
            database_reference=FILE_DATABASE_REFERENCE_LIGATION,
        )


class TestCrossHybridizationFilter(unittest.TestCase):
    def setUp(self):
        self.tmp_path = os.path.join(os.getcwd(), "tmp_output_crosshybridization_filter")
        os.makedirs(self.tmp_path, exist_ok=True)
        self.oligo_database_crosshyb = self._setup_database(FILE_DATABASE_OLIGOS_CROSSHYB)
        self.oligo_database_crosshyb_exactmatch = self._setup_database(FILE_DATABASE_OLIGOS_EXACT_MATCH)
        self.sequence_type = "oligo"

        # Blast parameters
        self.blastn_search_parameters_crosshyb = {
            "perc_identity": 80,
            "strand": "minus",
            "word_size": 10,
        }
        self.hit_parameters_crosshyb = {"coverage": 50}

        # Bowtie parameters
        self.bowtie_search_parameters_crosshyb = {"-n": 3, "-l": 5, "--nofw": ""}

        self.expected_oligos_larger_region = []
        for i, solution_file in enumerate(SOLUTIONS_LARGER_REGION):
            solution = OligoDatabase(
                min_oligos_per_region=2,
                write_regions_with_insufficient_oligos=True,
                database_name="db_oligo_crosshybridization_filters_solution_larger_region_{i}",
                dir_output=self.tmp_path,
            )
            solution.load_database_from_table(solution_file, database_overwrite=True)
            self.expected_oligos_larger_region.append(solution.database)

        self.expected_oligos_degree = []
        for i, solution_file in enumerate(SOLUTIONS_DEGREE):
            solution = OligoDatabase(
                min_oligos_per_region=2,
                write_regions_with_insufficient_oligos=True,
                database_name="db_oligo_crosshybridization_filters_solution_degree_{i}",
                dir_output=self.tmp_path,
            )
            solution.load_database_from_table(solution_file, database_overwrite=True)
            self.expected_oligos_degree.append(solution.database)

    def tearDown(self):
        shutil.rmtree(self.tmp_path)

    def _setup_database(self, file_database):
        oligos = OligoDatabase(
            min_oligos_per_region=2,
            write_regions_with_insufficient_oligos=True,
            database_name="db_oligo_crosshybridization_filters",
            dir_output=self.tmp_path,
        )
        oligos.load_database_from_table(file_database, database_overwrite=True)
        return oligos

    def _apply_filter_and_assert(self, filter_instance, expected_oligos):
        res = filter_instance.apply(
            oligo_database=self.oligo_database_crosshyb,
            reference_database=None,
            sequence_type=self.sequence_type,
            n_jobs=2,
        )
        assert (
            res.database in expected_oligos
        ), f"The cross-hybridization filter didn't return the expected oligos."

    def test_crosshyb_filter_blast_larger_region_policy(self):
        filter_instance = BlastNFilter(
            self.blastn_search_parameters_crosshyb,
            self.hit_parameters_crosshyb,
            filter_name="blast_larger_region",
            dir_output=self.tmp_path,
        )
        policy = RemoveByLargerRegionPolicy()
        cross_hyb_filter = CrossHybridizationFilter(
            policy=policy,
            alignment_method=filter_instance,
            filter_name="crosshybridization_blast_larger_region",
            dir_output=self.tmp_path,
        )
        self._apply_filter_and_assert(cross_hyb_filter, self.expected_oligos_larger_region)

    def test_crosshyb_filter_blast_degree_policy(self):
        filter_instance = BlastNFilter(
            self.blastn_search_parameters_crosshyb,
            self.hit_parameters_crosshyb,
            filter_name="blast_degree",
            dir_output=self.tmp_path,
        )
        policy = RemoveByDegreePolicy()
        cross_hyb_filter = CrossHybridizationFilter(
            policy=policy,
            alignment_method=filter_instance,
            filter_name="crosshybridization_blast_degree",
            dir_output=self.tmp_path,
        )
        self._apply_filter_and_assert(cross_hyb_filter, self.expected_oligos_degree)

    def test_crosshyb_filter_bowtie_larger_region_policy(self):
        filter_instance = BowtieFilter(
            self.bowtie_search_parameters_crosshyb,
            filter_name="bowtie_larger_region",
            dir_output=self.tmp_path,
        )
        policy = RemoveByLargerRegionPolicy()
        cross_hyb_filter = CrossHybridizationFilter(
            policy=policy,
            alignment_method=filter_instance,
            filter_name="crosshybridization_bowtie_larger_region",
            dir_output=self.tmp_path,
        )
        self._apply_filter_and_assert(cross_hyb_filter, self.expected_oligos_larger_region)

    def test_crosshyb_filter_bowtie_degree_policy(self):
        filter_instance = BowtieFilter(
            self.bowtie_search_parameters_crosshyb,
            filter_name="bowtie_degree",
            dir_output=self.tmp_path,
        )
        policy = RemoveByDegreePolicy()
        cross_hyb_filter = CrossHybridizationFilter(
            policy=policy,
            alignment_method=filter_instance,
            filter_name="crosshybridization_bowtie_degree",
            dir_output=self.tmp_path,
        )
        self._apply_filter_and_assert(cross_hyb_filter, self.expected_oligos_degree)


class DummyAPI(APIBase):
    # Class that considers real hits all the hits that have a 100% match
    def predict(self, queries, gapped_queries, references, gapped_references):
        predictions = np.ndarray(shape=(len(queries),), dtype=np.float32)
        for i, (q, r) in enumerate(zip(gapped_queries, gapped_references)):
            if q == r:
                predictions[i] = 1
            else:
                predictions[i] = 0
        return predictions


class TestHybridizationProbabilityBalstn(unittest.TestCase):
    def setUp(self) -> None:
        self.tmp_path = os.path.join(os.getcwd(), "tmp_output_hybridization_probability_filter_blast")
        blastn_search_parameters = {
            "perc_identity": 80,
            "strand": "both",
            "word_size": 10,
        }
        hit_parameters = {"coverage": 50}
        self.alignment_filter = BlastNFilter(
            search_parameters=blastn_search_parameters,
            hit_parameters=hit_parameters,
            dir_output=self.tmp_path,
        )
        self.filter = HybridizationProbabilityFilter(
            alignment_method=self.alignment_filter,
            threshold=0.1,
            dir_output=self.tmp_path,
        )
        self.filter.model = DummyAPI()
        self.database = OligoDatabase(dir_output=self.tmp_path)
        self.database.load_database_from_table(FILE_DATABASE_OLIGOS_AI, database_overwrite=True)
        self.reference_database = ReferenceDatabase(dir_output=self.tmp_path)
        self.reference_database.load_database_from_fasta(
            files_fasta=FILE_DATABASE_REFERENCE, database_overwrite=True
        )
        self.file_reference = self.reference_database.write_database_to_fasta(filename="db_reference")
        self.table_hits = pd.read_csv(FILE_TABLE_HITS_BLAST_AI, sep="\t")
        self.sequence_type = "target"
        self.region_id = "region"

    def tearDown(self) -> None:
        shutil.rmtree(self.tmp_path)

    def test_ai_filter_blastn(self):
        filtered_database = self.filter.apply(
            sequence_type=self.sequence_type,
            oligo_database=self.database,
            reference_database=self.reference_database,
            n_jobs=2,
        )
        returned_oligos = set(filtered_database.database["region"].keys())
        expected_oligos = set(f"region::{i}" for i in range(2, 20))

        assert (
            returned_oligos == expected_oligos
        ), f"The Blast ai filter didn't return the expected oligos. \n\nExpected:\n{expected_oligos}\n\nGot:\n{returned_oligos}"

    def test_get_queries(self):
        returned_queries = self.alignment_filter._get_queries(
            oligo_database=self.database,
            table_hits=self.table_hits,
            sequence_type=self.sequence_type,
            region_id=self.region_id,
        )
        returned_queries = set(returned_queries)
        expected_queries = set(
            [
                Seq("GCTCGGGCTTGTCCACAGGATGGACCCAGCTGAGCAAGCT"),
                Seq("AGCTTGCTCAGCTGGGTCCATCCTGTGGACAAGCCCGAGC"),
                Seq("TACAGGCATGACCCACCATGCCTGGCCAACTTACATTTTT"),
                Seq("AAAAATGTAAGTTGGCCAGGCATGGTGGGTCATGCCTGTA"),
                Seq("AAGGCCAAGGTCTCTGGGGGGCTGGACAAGCCGCCCTCAT"),
                Seq("ATGAGGGCGGCTTGTCCAGCCCCCCAGAGACCTTGGCCTT"),
                Seq("TTTTGCACCAGCCCAGATCGCATCTTCTTTCACCTGTTTT"),
                Seq("AAAACAGGTGAAAGAAGATGCGATCTGGGCTGGTGCAAAA"),
                Seq("CCGCTCGGCTGCATGAAACCAAAACGGCTGTCCGGGGACA"),
                Seq("TGTCCCCGGACAGCCGTTTTGGTTTCATGCAGCCGAGCGG"),
                Seq("AACCCGGCATCACCAAGAGGAGGTTCAAGGGAACGCTGCA"),
                Seq("TGCAGCGTTCCCTTGAACCTCCTCTTGGTGATGCCGGGTT"),
                Seq("TGCCCGCGCCGGAGTTCTCCCCAGCCGGAGTCCGGCAGGG"),
                Seq("CCCTGCCGGACTCCGGCTGGGGAGAACTCCGGCGCGGGCA"),
                Seq("AACCTGGTTGCACCTCGGCCTGGTCCCAGCAGGTATGGTT"),
                Seq("AACCATACCTGCTGGGACCAGGCCGAGGTGCAACCAGGTT"),
                Seq("ACTGATTGCTGCAGACGCTCACCCCAGACACTCACTGCAC"),
                Seq("GTGCAGTGAGTGTCTGGGGTGAGCGTCTGCAGCAATCAGT"),
                Seq("TATATATTTTGCACACTTTAAAATATTGGGTTGTTTACCG"),
                Seq("CGGTAAACAACCCAATATTTTAAAGTGTGCAAAATATATA"),
            ]
        )
        assert (
            returned_queries == expected_queries
        ), f"The Blast ai filter didn't return the expected queries. \n\nExpected:\n{expected_queries}\n\nGot:\n{returned_queries}"

    def test_get_target_blastn(self):
        returned_references = self.alignment_filter._get_references(
            table_hits=self.table_hits, file_reference=self.file_reference, region_id=self.region_id
        )
        returned_references = set(returned_references)
        expected_references = set(
            [
                Seq("GCTCGGGCTTGTCCACAGGATGGACCCAGCTGAGCAAGCT"),
                Seq("AGCTTGCTCAGCTGGGTCCATCCTGTGGACAAGCCCGAGC"),
                Seq("TACAGGCATGAGCCACCATGCCTGGCCAACTCACATTTTT"),
                Seq("AAAAATGTGAGTTGGCCAGGCATGGTGGCTCATGCCTGTA"),
                Seq("AAGGCCGGGGTCTCTGGGGGGCTGGAGAAGCCTCCCTCAT"),
                Seq("ATGAGGGAGGCTTCTCCAGCCCCCCAGAGACCCCGGCCTT"),
                Seq("AGCAGCACCAGCCCAGATCGCATCTTCTTTCACCTGAACG"),
                Seq("CGTTCAGGTGAAAGAAGATGCGATCTGGGCTGGTGCTGCT"),
                Seq("CCGCTACCGGCTGCATGACAACCAAAACGGCTGGTCCGGGGACA"),
                Seq("TGTCCCCGGACCAGCCGTTTTGGTTGTCATGCAGCCGGTAGCGG"),
                Seq("AACCCCATCACCAAGAGGAGGTTCAGGGAAGCTGCA"),
                Seq("TGCAGCTTCCCTGAACCTCCTCTTGGTGATGGGGTT"),
                Seq("TGCCCGCGCCGGAGTTCTCCCCGGAGCCGGAGTCCGGCAGGG"),
                Seq("CCCTGCCGGACTCCGGCTCCGGGGAGAACTCCGGCGCGGGCA"),
                Seq("TCCCTGGGCACCTCGGCCTGGTCCCAGCAGGTATGGGC"),
                Seq("GCCCATACCTGCTGGGACCAGGCCGAGGTGCCCAGGGA"),
                Seq("----ATTGCTGCAGACGCTCACCCCAGACACTCACTGCAC"),
                Seq("GTGCAGTGAGTGTCTGGGGTGAGCGTCTGCAGCAAT----"),
                Seq("TATATATTTTGCACACTTTAAAATATTGGGTTGTTT----"),
                Seq("----AAACAACCCAATATTTTAAAGTGTGCAAAATATATA"),
            ]
        )
        assert (
            returned_references == expected_references
        ), f"The Blast ai filter didn't return the expected references. \n\nExpected:\n{expected_references}\n\nGot:\n{returned_references}"

    def test_add_alignment_gaps_queries(self):
        queries = self.alignment_filter._get_queries(
            oligo_database=self.database,
            table_hits=self.table_hits,
            sequence_type=self.sequence_type,
            region_id=self.region_id,
        )
        references = self.alignment_filter._get_references(
            table_hits=self.table_hits, file_reference=self.file_reference, region_id=self.region_id
        )
        gapped_queries, _ = self.alignment_filter._add_alignment_gaps(
            table_hits=self.table_hits,
            queries=queries,
            references=references,
        )
        gapped_queries = set(gapped_queries)
        expected_gapped_queries = set(
            [
                Seq("GCTCGGGCTTGTCCACAGGATGGACCCAGCTGAGCAAGCT"),
                Seq("AGCTTGCTCAGCTGGGTCCATCCTGTGGACAAGCCCGAGC"),
                Seq("TACAGGCATGACCCACCATGCCTGGCCAACTTACATTTTT"),
                Seq("AAAAATGTAAGTTGGCCAGGCATGGTGGGTCATGCCTGTA"),
                Seq("AAGGCCAAGGTCTCTGGGGGGCTGGACAAGCCGCCCTCAT"),
                Seq("ATGAGGGCGGCTTGTCCAGCCCCCCAGAGACCTTGGCCTT"),
                Seq("TTTTGCACCAGCCCAGATCGCATCTTCTTTCACCTGTTTT"),
                Seq("AAAACAGGTGAAAGAAGATGCGATCTGGGCTGGTGCAAAA"),
                Seq("CCGCT--CGGCTGCATGA-AACCAAAACGGCTG-TCCGGGGACA"),
                Seq("TGTCCCCGGA-CAGCCGTTTTGGTT-TCATGCAGCCG--AGCGG"),
                Seq("AACCCGGCATCACCAAGAGGAGGTTCAAGGGAACGCTGCA"),
                Seq("TGCAGCGTTCCCTTGAACCTCCTCTTGGTGATGCCGGGTT"),
                Seq("TGCCCGCGCCGGAGTTCTCCCC--AGCCGGAGTCCGGCAGGG"),
                Seq("CCCTGCCGGACTCCGGCT--GGGGAGAACTCCGGCGCGGGCA"),
                Seq("AACCTGGTTGCACCTCGGCCTGGTCCCAGCAGGTATGGTT"),
                Seq("AACCATACCTGCTGGGACCAGGCCGAGGTGCAACCAGGTT"),
                Seq("ACTGATTGCTGCAGACGCTCACCCCAGACACTCACTGCAC"),
                Seq("GTGCAGTGAGTGTCTGGGGTGAGCGTCTGCAGCAATCAGT"),
                Seq("TATATATTTTGCACACTTTAAAATATTGGGTTGTTTACCG"),
                Seq("CGGTAAACAACCCAATATTTTAAAGTGTGCAAAATATATA"),
            ]
        )
        assert (
            gapped_queries == expected_gapped_queries
        ), f"The Blast ai filter didn't return the expected gapped queries. \n\nExpected:\n{expected_gapped_queries}\n\nGot:\n{gapped_queries}"

    def test_add_alignment_gaps_references(self):
        queries = self.alignment_filter._get_queries(
            oligo_database=self.database,
            table_hits=self.table_hits,
            sequence_type=self.sequence_type,
            region_id=self.region_id,
        )
        references = self.alignment_filter._get_references(
            table_hits=self.table_hits, file_reference=self.file_reference, region_id=self.region_id
        )
        _, gapped_references = self.alignment_filter._add_alignment_gaps(
            table_hits=self.table_hits,
            queries=queries,
            references=references,
        )
        gapped_references = set(gapped_references)
        expected_gapped_references = set(
            [
                Seq("GCTCGGGCTTGTCCACAGGATGGACCCAGCTGAGCAAGCT"),
                Seq("AGCTTGCTCAGCTGGGTCCATCCTGTGGACAAGCCCGAGC"),
                Seq("TACAGGCATGAGCCACCATGCCTGGCCAACTCACATTTTT"),
                Seq("AAAAATGTGAGTTGGCCAGGCATGGTGGCTCATGCCTGTA"),
                Seq("AAGGCCGGGGTCTCTGGGGGGCTGGAGAAGCCTCCCTCAT"),
                Seq("ATGAGGGAGGCTTCTCCAGCCCCCCAGAGACCCCGGCCTT"),
                Seq("AGCAGCACCAGCCCAGATCGCATCTTCTTTCACCTGAACG"),
                Seq("CGTTCAGGTGAAAGAAGATGCGATCTGGGCTGGTGCTGCT"),
                Seq("CCGCTACCGGCTGCATGACAACCAAAACGGCTGGTCCGGGGACA"),
                Seq("TGTCCCCGGACCAGCCGTTTTGGTTGTCATGCAGCCGGTAGCGG"),
                Seq("AACCC--CATCACCAAGAGGAGGTTCA-GGGAA-GCTGCA"),
                Seq("TGCAGC-TTCCC-TGAACCTCCTCTTGGTGATG--GGGTT"),
                Seq("TGCCCGCGCCGGAGTTCTCCCCGGAGCCGGAGTCCGGCAGGG"),
                Seq("CCCTGCCGGACTCCGGCTCCGGGGAGAACTCCGGCGCGGGCA"),
                Seq("TCCCTGG--GCACCTCGGCCTGGTCCCAGCAGGTATGGGC"),
                Seq("GCCCATACCTGCTGGGACCAGGCCGAGGTGC--CCAGGGA"),
                Seq("----ATTGCTGCAGACGCTCACCCCAGACACTCACTGCAC"),
                Seq("GTGCAGTGAGTGTCTGGGGTGAGCGTCTGCAGCAAT----"),
                Seq("TATATATTTTGCACACTTTAAAATATTGGGTTGTTT----"),
                Seq("----AAACAACCCAATATTTTAAAGTGTGCAAAATATATA"),
            ]
        )
        assert (
            gapped_references == expected_gapped_references
        ), f"The Blast ai filter didn't return the expected gapped references. \n\nExpected:\n{expected_gapped_references}\n\nGot:\n{gapped_references}"


class TestHybridizationProbabilityBowtie(unittest.TestCase):
    def setUp(self) -> None:
        self.tmp_path = os.path.join(os.getcwd(), "tmp_output_hybridization_probability_filter_bowtie")
        bowtie_search_parameters = {"-n": 3, "-l": 5}
        self.alignment_filter = BowtieFilter(
            search_parameters=bowtie_search_parameters,
            dir_output=self.tmp_path,
        )
        self.filter = HybridizationProbabilityFilter(
            alignment_method=self.alignment_filter,
            threshold=0.1,
            dir_output=self.tmp_path,
        )
        self.filter.model = DummyAPI()
        self.database = OligoDatabase(dir_output=self.tmp_path)
        self.database.load_database_from_table(FILE_DATABASE_OLIGOS_AI, database_overwrite=True)
        self.reference_database = ReferenceDatabase(dir_output=self.tmp_path)
        self.reference_database.load_database_from_fasta(
            files_fasta=FILE_DATABASE_REFERENCE, database_overwrite=True
        )
        self.file_reference = self.reference_database.write_database_to_fasta(filename="db_reference_bowtie")
        self.table_hits = pd.read_csv(FILE_TABLE_HITS_BOWTIE_AI, sep="\t")
        self.sequence_type = "target"
        self.region_id = "region"

    def tearDown(self) -> None:
        shutil.rmtree(self.tmp_path)

    def test_ai_filter_bowtie(self):
        filtered_database = self.filter.apply(
            sequence_type=self.sequence_type,
            oligo_database=self.database,
            reference_database=self.reference_database,
            n_jobs=2,
        )
        returned_oligos = set(filtered_database.database["region"].keys())
        expected_oligos = set(f"region::{i}" for i in range(2, 20))

        assert (
            returned_oligos == expected_oligos
        ), f"The Bowtie ai filter didn't return the expected oligos. \n\nExpected:\n{expected_oligos}\n\nGot:\n{returned_oligos}"

    def test_get_queries(self):
        returned_queries = self.alignment_filter._get_queries(
            oligo_database=self.database,
            table_hits=self.table_hits,
            sequence_type=self.sequence_type,
            region_id=self.region_id,
        )
        returned_queries = set(returned_queries)
        expected_queries = set(
            [
                Seq("GCTCGGGCTTGTCCACAGGATGGACCCAGCTGAGCAAGCT"),
                Seq("AGCTTGCTCAGCTGGGTCCATCCTGTGGACAAGCCCGAGC"),
                Seq("TACAGGCATGACCCACCATGCCTGGCCAACTTACATTTTT"),
                Seq("AAAAATGTAAGTTGGCCAGGCATGGTGGGTCATGCCTGTA"),
            ]
        )
        assert (
            returned_queries == expected_queries
        ), f"The Bowtie ai filter didn't return the expected queries. \n\nExpected:\n{expected_queries}\n\nGot:\n{returned_queries}"

    def test_get_target_bowtie(self):
        returned_references = self.alignment_filter._get_references(
            table_hits=self.table_hits, file_reference=self.file_reference, region_id=self.region_id
        )
        returned_references = set(returned_references)
        expected_references = set(
            [
                Seq("GCTCGGGCTTGTCCACAGGATGGACCCAGCTGAGCAAGCT"),
                Seq("AGCTTGCTCAGCTGGGTCCATCCTGTGGACAAGCCCGAGC"),
                Seq("TACAGGCATGAGCCACCATGCCTGGCCAACTCACATTTTT"),
                Seq("AAAAATGTGAGTTGGCCAGGCATGGTGGCTCATGCCTGTA"),
            ]
        )
        assert (
            returned_references == expected_references
        ), f"The Bowtie ai filter didn't return the expected references. \n\nExpected:\n{expected_references}\n\nGot:\n{returned_references}"
