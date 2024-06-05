############################################
# imports
############################################

import os
import shutil
import unittest

from oligo_designer_toolsuite.database import (
    OligoAttributes,
    OligoDatabase,
    ReferenceDatabase,
)
from oligo_designer_toolsuite.sequence_generator import OligoSequenceGenerator
from oligo_designer_toolsuite.utils import FastaParser

############################################
# setup
############################################

# Global Parameters
FILE_NCBI_EXONS = "tests/data/genomic_regions/sequences_ncbi_exons.fna"
FILE_DATABASE_OLIGO_ATTRIBUTES = "tests/data/databases/database_oligo_attributes.fna"

REGION_IDS = [
    "AARS1",
    "DECR2",
    "FAM234A",
    "RHBDF1",
    "WASIR2",
    "this_gene_does_not_exist",
]

############################################
# tests
############################################


class TestReferenceDatabase(unittest.TestCase):
    def setUp(self):
        self.tmp_path = os.path.join(os.getcwd(), "tmp_reference_database")

        self.fasta_parser = FastaParser()

        self.reference = ReferenceDatabase(database_name="test_reference_database", dir_output=self.tmp_path)
        # self.reference.load_sequences_from_fasta(
        #     files_fasta=[FILE_NCBI_EXONS, FILE_NCBI_EXONS], database_overwrite=True
        # )
        self.reference.load_sequences_from_fasta(files_fasta=FILE_NCBI_EXONS, database_overwrite=False)

    def tearDown(self):
        try:
            shutil.rmtree(self.tmp_path)
        except:
            pass

    def test_filter_database(self):
        """Test creation of reference database as well as load, write and filter functionalities."""
        self.reference.filter_database("AARS1", remove_region=True)
        for entry in self.reference.database:
            (
                region,
                _,
                _,
            ) = self.fasta_parser.parse_fasta_header(entry.id)
            assert region != "AARS1", f"error: this region {region} should be filtered out."

    def test_write_database(self):
        print("#" * 100)
        print(self.reference.database)
        print("#" * 100)
        file_fasta_database = self.reference.write_database_to_fasta(filename="filtered_databse")
        assert (
            self.fasta_parser.check_fasta_format(file_fasta_database) == True
        ), f"error: wrong file format for database in {file_fasta_database}"


class TestOligoDatabase(unittest.TestCase):
    def setUp(self):
        self.tmp_path = os.path.join(os.getcwd(), "tmp_oligo_database")

        self.fasta_parser = FastaParser()

        self.oligo_sequence_generator = OligoSequenceGenerator(dir_output=self.tmp_path)
        self.oligo_database = OligoDatabase(
            min_oligos_per_region=2,
            write_regions_with_insufficient_oligos=True,
            database_name="test_oligo_database",
            dir_output=self.tmp_path,
        )

        self.file_random_seqs = self.oligo_sequence_generator.create_sequences_random(
            filename_out="random_sequences1",
            length_sequences=30,
            num_sequences=100,
            name_sequences="random_sequences1",
            base_alphabet_with_probability={"A": 0.1, "C": 0.3, "G": 0.4, "T": 0.2},
        )

        self.file_sliding_window = self.oligo_sequence_generator.create_sequences_sliding_window(
            filename_out="sliding_window_sequences",
            files_fasta_in=FILE_NCBI_EXONS,
            length_interval_sequences=(30, 31),
            n_jobs=2,
        )

    def tearDown(self):
        shutil.rmtree(self.tmp_path)

    def test_load_sequences_from_fasta(self):
        self.oligo_database.load_sequences_from_fasta(
            files_fasta=self.file_random_seqs,
            sequence_type="oligo",
            region_ids=["random_sequences1"],
            database_overwrite=True,
        )
        self.oligo_database.load_sequences_from_fasta(
            files_fasta=self.file_sliding_window,
            sequence_type="target",
            region_ids=REGION_IDS,
            database_overwrite=False,
        )

        assert len(self.oligo_database.database) > 0, "error: no sequences loaded into database"

    def test_save_load_database(self):
        self.oligo_database.load_from_tsv(FILE_DATABASE_OLIGO_ATTRIBUTES, database_overwrite=True)

        file_database = self.oligo_database.save_database(
            region_ids=["region_1", "region_2"], dir_database="database_region1_region2"
        )

        self.oligo_database.load_from_tsv(file_database, database_overwrite=True)

        assert len(self.oligo_database.database.keys()) == 2, "error: wrong number regions saved and loaded"

    def test_write_database_to_fasta(self):
        self.oligo_database.load_from_tsv(
            FILE_DATABASE_OLIGO_ATTRIBUTES, region_ids=["region_1", "region_2"], database_overwrite=True
        )
        file_fasta = self.oligo_database.write_database_to_fasta(filename="database_region1_region2")

        assert (
            self.fasta_parser.check_fasta_format(file_fasta) == True
        ), f"error: wrong file format for database in {file_fasta}"

    # TODO: add test
    def test_write_oligosets(self):
        pass

    def test_remove_regions_with_insufficient_oligos(self):
        self.oligo_database.load_sequences_from_fasta(
            files_fasta=self.file_sliding_window,
            sequence_type="target",
            region_ids=REGION_IDS,
            database_overwrite=True,
        )
        self.oligo_database.load_sequences_from_fasta(
            files_fasta=self.file_random_seqs,
            sequence_type="oligo",
            database_overwrite=False,
        )

        self.oligo_database.remove_regions_with_insufficient_oligos("database_generation")
        assert len(self.oligo_database.database.keys()) == (
            len(REGION_IDS) - 1 + 1  # one region removed but one added from random seqs
        ), "error: wrong number of regions in database"

    def test_get_sequence_list(self):
        self.oligo_database.load_sequences_from_fasta(
            files_fasta=self.file_random_seqs,
            sequence_type="oligo",
            database_overwrite=True,
        )

        list_sequences = self.oligo_database.get_sequence_list()
        assert len(list_sequences) == 100, "error: wrong number of sequences in database"

    def test_get_sequence_oligoid_mapping(self):
        self.oligo_database.load_from_tsv(
            file_database=FILE_DATABASE_OLIGO_ATTRIBUTES,
            region_ids=None,
            database_overwrite=True,
        )

        mapping = self.oligo_database.get_sequence_oligoid_mapping(sequence_type="oligo")
        assert len(mapping["CTCACTCGACTCTTACACAGTCATA"]) == 4, "error: wrong number of oligos for sequence"

    def test_get_oligo_attribute(self):
        self.oligo_database.load_from_tsv(
            file_database=FILE_DATABASE_OLIGO_ATTRIBUTES,
            region_ids=None,
            database_overwrite=True,
        )
        attribute = self.oligo_database.get_oligo_attribute(attribute="test_attribute")

        assert len(attribute["test_attribute"].unique()) == 2, "error: wrong attribute returned"

    def test_update_oligo_attribute(self):
        self.oligo_database.load_from_tsv(
            file_database=FILE_DATABASE_OLIGO_ATTRIBUTES,
            region_ids="region_3",
            database_overwrite=True,
        )
        new_attribute = {
            "region_3::1": {"GC_content": 63},
            "region_3::2": {"GC_content": 66},
            "region_3::3": {"GC_content": 80},
            "region_3::4": {"GC_content": 70},
            "region_3::5": {"GC_content": 40},
        }
        self.oligo_database.update_oligo_attributes(new_attribute)
        attribute = self.oligo_database.get_oligo_attribute(attribute="GC_content")

        assert len(attribute) == 5, "error: attribute not correctly updated"


class TestOligoAttributes(unittest.TestCase):
    def setUp(self):
        self.tmp_path = os.path.join(os.getcwd(), "tmp_oligo_attributes")

        self.oligo_database = OligoDatabase(
            min_oligos_per_region=2,
            write_regions_with_insufficient_oligos=True,
            database_name="test_oligo_attributes",
            dir_output=self.tmp_path,
        )
        self.oligo_database.load_from_tsv(FILE_DATABASE_OLIGO_ATTRIBUTES, database_overwrite=True)

        self.oligo_attributes = OligoAttributes()

    def tearDown(self):
        shutil.rmtree(self.tmp_path)

    def test_calculate_oligo_length(self):
        oligo_database = self.oligo_attributes.calculate_oligo_length(self.oligo_database)

        assert oligo_database.database["region_1"]["region_1::1"]["length"] == 20, "error: wrong oligo length"
        assert oligo_database.database["region_1"]["region_1::2"]["length"] == 29, "error: wrong oligo length"

    def test_calculate_num_targeted_transcripts(self):
        oligo_database = self.oligo_attributes.calculate_num_targeted_transcripts(self.oligo_database)

        assert (
            oligo_database.database["region_1"]["region_1::1"]["num_targeted_transcripts"] == 2
        ), "error: wrong number targeted transcripts"
        assert (
            oligo_database.database["region_2"]["region_2::1"]["num_targeted_transcripts"] == 1
        ), "error: wrong number targeted transcripts"
        assert (
            oligo_database.database["region_3"]["region_3::4"]["num_targeted_transcripts"] == 28
        ), "error: wrong number targeted transcripts"

    def test_calculate_isoform_consensus(self):
        oligo_database = self.oligo_attributes.calculate_isoform_consensus(self.oligo_database)

        assert (
            oligo_database.database["region_1"]["region_1::1"]["isoform_consensus"] == 100
        ), "error: wrong isoform consensus, should be 100%"
        assert (
            oligo_database.database["region_2"]["region_2::1"]["isoform_consensus"] == 50
        ), "error: wrong isoform consensus, should be 50%"

    def test_calculate_seedregion(self):
        oligo_database = self.oligo_attributes.calculate_seedregion(self.oligo_database, start=0.4, end=0.6)

        assert (oligo_database.database["region_1"]["region_1::1"]["seedregion_start"] == 8) and (
            oligo_database.database["region_1"]["region_1::1"]["seedregion_end"] == 12
        ), "error: wrong seedregion calculated"

    def test_calculate_seedregion_ligationsite(self):
        oligo_database = self.oligo_attributes.calculate_seedregion_ligationsite(
            self.oligo_database, seedregion_size=5
        )

        assert (oligo_database.database["region_1"]["region_1::1"]["seedregion_start"] == 6) and (
            oligo_database.database["region_1"]["region_1::1"]["seedregion_end"] == 15
        ), "error: wrong seedregion calculated"

    def test_calculate_GC_content(self):
        oligo_database = self.oligo_attributes.calculate_GC_content(
            self.oligo_database, sequence_type="oligo"
        )

        assert (
            oligo_database.database["region_1"]["region_1::1"]["GC_content"] == 50
        ), "error: wrong GC content calculated"

    def test_calculate_TmNN(self):
        oligo_database = self.oligo_attributes.calculate_TmNN(
            self.oligo_database, sequence_type="oligo", Tm_parameters={}
        )

        assert (
            oligo_database.database["region_1"]["region_1::1"]["TmNN"] == 53.57
        ), "error: wrong Tm calculated"

    def test_calculate_len_selfcomp(self):
        oligo_database = self.oligo_attributes.calculate_length_selfcomplement(
            self.oligo_database, sequence_type="oligo"
        )

        assert (
            oligo_database.database["region_3"]["region_3::3"]["length_selfcomplement"] == 18
        ), "error: wrong length of selfcomplement calculated"

    def test_calculate_secondary_structure_DG(self):
        oligo_database = self.oligo_attributes.calculate_DG_secondary_structure(
            self.oligo_database, sequence_type="oligo", T=37
        )

        assert (
            oligo_database.database["region_1"]["region_1::1"]["DG_secondary_structure"] == 0.8
        ), "error: wrong DG calculated"

    def test_calculate_padlock_arms(self):
        oligo_database = self.oligo_attributes.calculate_padlock_arms(
            self.oligo_database,
            arm_length_min=3,
            arm_Tm_dif_max=15,
            arm_Tm_min=30,
            arm_Tm_max=80,
            Tm_parameters={},
        )

        assert (
            oligo_database.database["region_1"]["region_1::2"]["ligation_site"] == 14
        ), "error: wrong padlock arms calculated"
