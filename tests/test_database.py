############################################
# imports
############################################

import os
import yaml
import shutil
import unittest
import pandas as pd

from oligo_designer_toolsuite.database import (
    OligoAttributes,
    OligoDatabase,
    ReferenceDatabase,
)
from oligo_designer_toolsuite.sequence_generator import OligoSequenceGenerator
from oligo_designer_toolsuite.utils import FastaParser, check_tsv_format

############################################
# setup
############################################

# Global Parameters
FILE_NCBI_EXONS = "tests/data/genomic_regions/sequences_ncbi_exons.fna"
FILE_DATABASE_OLIGO_ATTRIBUTES = "tests/data/databases/database_oligo_attributes.tsv"

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
        self.reference.load_database_from_fasta(files_fasta=[FILE_NCBI_EXONS], database_overwrite=False)
        self.reference.load_database_from_fasta(files_fasta=FILE_NCBI_EXONS, database_overwrite=True)

    def tearDown(self):
        try:
            shutil.rmtree(self.tmp_path)
        except:
            pass

    def test_write_database(self):
        file_fasta_database = self.reference.write_database_to_fasta(filename="filtered_databse")
        assert (
            self.fasta_parser.check_fasta_format(file_fasta_database) == True
        ), f"error: wrong file format for database in {file_fasta_database}"

    def test_filter_database_by_region(self):
        self.reference.filter_database_by_region(region_ids="AARS1", keep_region=False)
        for entry in self.reference.database:
            (
                region,
                _,
                _,
            ) = self.fasta_parser.parse_fasta_header(entry.id)
            assert region != "AARS1", f"error: this region {region} should be filtered out."

    def test_filter_database_by_attribute_category(self):
        self.reference.filter_database_by_attribute_category(
            attribute_name="gene_id", attribute_category="AARS1", keep_if_equals_category=False
        )
        for entry in self.reference.database:
            (
                region,
                _,
                _,
            ) = self.fasta_parser.parse_fasta_header(entry.id)
            assert region != "AARS1", f"error: this region {region} should be filtered out."


class TestOligoDatabase(unittest.TestCase):
    def setUp(self):
        self.tmp_path = os.path.join(os.getcwd(), "tmp_oligo_database")

        self.fasta_parser = FastaParser()

        self.oligo_sequence_generator = OligoSequenceGenerator(dir_output=self.tmp_path)
        self.oligo_database = OligoDatabase(
            min_oligos_per_region=2,
            write_regions_with_insufficient_oligos=True,
            lru_db_max_in_memory=10,
            n_jobs=4,
            database_name="test_oligo_database",
            dir_output=self.tmp_path,
        )

    def tearDown(self):
        shutil.rmtree(self.tmp_path)

    def load_database_from_fasta(self):
        file_random_seqs = self.oligo_sequence_generator.create_sequences_random(
            filename_out="random_sequences1",
            length_sequences=30,
            num_sequences=100,
            name_sequences="random_sequences1",
            base_alphabet_with_probability={"A": 0.1, "C": 0.3, "G": 0.4, "T": 0.2},
        )
        file_sliding_window = self.oligo_sequence_generator.create_sequences_sliding_window(
            files_fasta_in=FILE_NCBI_EXONS,
            region_ids=REGION_IDS,
            length_interval_sequences=(30, 31),
            n_jobs=4,
        )
        self.oligo_database.load_database_from_fasta(
            files_fasta=file_random_seqs,
            sequence_type="oligo",
            region_ids="random_sequences1",
            database_overwrite=True,
        )
        self.oligo_database.load_database_from_fasta(
            files_fasta=file_sliding_window,
            sequence_type="target",
            region_ids=REGION_IDS,
            database_overwrite=False,
        )

    def test_load_database_from_fasta(self):
        self.load_database_from_fasta()

        assert len(self.oligo_database.database) == 6, "error: wrong number of sequences loaded into database"

    def test_load_database_from_table(self):
        self.oligo_database.load_database_from_table(
            FILE_DATABASE_OLIGO_ATTRIBUTES, region_ids=["region_1", "region_2"], database_overwrite=True
        )

        assert len(self.oligo_database.database) == 2, "error: wrong number of sequences loaded into database"

    def test_load_save_database(self):
        self.oligo_database.load_database_from_table(FILE_DATABASE_OLIGO_ATTRIBUTES, database_overwrite=True)

        dir_database = self.oligo_database.save_database(
            region_ids=["region_1", "region_2"], dir_database="database_region1_region2"
        )
        self.oligo_database.load_database(dir_database, database_overwrite=True)

        assert len(self.oligo_database.database.keys()) == 2, "error: wrong number regions saved and loaded"

    def test_write_database_to_fasta(self):
        self.oligo_database.load_database_from_table(
            FILE_DATABASE_OLIGO_ATTRIBUTES, region_ids=["region_1", "region_2"], database_overwrite=True
        )
        file_fasta = self.oligo_database.write_database_to_fasta(
            sequence_type="oligo", save_description=True, filename="database_region1_region2"
        )

        assert (
            self.fasta_parser.check_fasta_format(file_fasta) == True
        ), f"error: wrong file format for database in {file_fasta}"

        assert (
            len(self.fasta_parser.get_fasta_regions(file_fasta_in=file_fasta)) == 2
        ), f"error: wrong number of regions stored in {file_fasta}"

    def test_write_database_to_table(self):
        self.oligo_database.load_database_from_table(FILE_DATABASE_OLIGO_ATTRIBUTES, database_overwrite=True)

        file_database = self.oligo_database.write_database_to_table(
            attributes=["test_attribute", "ligation_site", "chromosome", "start", "end", "strand"],
            flatten_attribute=True,
            filename="database_region1_region2_flattened",
            region_ids=["region_1", "region_2"],
        )

        self.oligo_database.load_database_from_table(file_database, database_overwrite=True)

        assert check_tsv_format(file_database) == True, f"error: wrong file format"
        assert len(self.oligo_database.database.keys()) == 2, "error: wrong number regions saved and loaded"
        assert (
            self.oligo_database.get_oligo_attribute_value(
                attribute="test_attribute", flatten=True, region_id="region_1", oligo_id="region_1::1"
            )
            == "red"
        ), f"error: wrong attribute stored in {file_database}"

        file_database = self.oligo_database.write_database_to_table(
            attributes=["test_attribute", "ligation_site", "chromosome", "start", "end", "strand"],
            flatten_attribute=False,
            filename="database_region1_region2_unflattened",
            region_ids=["region_1", "region_2"],
        )

        self.oligo_database.load_database_from_table(file_database, database_overwrite=True)

        assert check_tsv_format(file_database) == True, f"error: wrong file format"
        assert len(self.oligo_database.database.keys()) == 2, "error: wrong number regions saved and loaded"
        assert (
            self.oligo_database.get_oligo_attribute_value(
                attribute="test_attribute", flatten=True, region_id="region_1", oligo_id="region_1::1"
            )
            == "red"
        ), f"error: wrong attribute stored in {file_database}"

    def test_write_oligosets_to_yaml(self):
        self.oligo_database.load_database_from_table(
            FILE_DATABASE_OLIGO_ATTRIBUTES, database_overwrite=True, region_ids="region_1"
        )

        oligoset = pd.DataFrame(
            data=[
                [0, "region_1::1", "region_1::2", "region_1::5", 1.59, 2.36],
                [1, "region_1::6", "region_1::4", "region_1::9", 2.15, 4.93],
            ],
            columns=[
                "oligoset_id",
                "oligo_0",
                "oligo_1",
                "oligo_2",
                "set_score_lowest",
                "set_score_sum",
            ],
        )

        self.oligo_database.oligosets["region_1"] = oligoset

        file_yaml = self.oligo_database.write_oligosets_to_yaml(
            attributes=[
                "test_attribute",
                "ligation_site",
                "chromosome",
                "start",
                "end",
                "strand",
                "transcript_id",
            ],
            top_n_sets=2,
            ascending=True,
        )

        with open(file_yaml, "r") as handle:
            yaml_oligosets = yaml.safe_load(handle)

        assert yaml_oligosets["region_1"]["Oligoset 1"]["Oligoset Score"] == {
            "set_score_lowest": 1.59,
            "set_score_sum": 2.36,
        }, f"error: wrong oligoset loaded"
        assert yaml_oligosets["region_1"]["Oligoset 1"]["Oligo 1"]["test_attribute"] == [
            ["red"]
        ], f"error: wrong oligoset loaded"

    def test_write_oligosets_to_table(self):
        self.oligo_database.load_database_from_table(
            FILE_DATABASE_OLIGO_ATTRIBUTES, database_overwrite=True, region_ids="region_1"
        )

        oligoset = pd.DataFrame(
            data=[
                [0, "region_1::1", "region_1::2", "region_1::5", 1.59, 2.36],
                [1, "region_1::6", "region_1::4", "region_1::9", 2.15, 4.93],
            ],
            columns=[
                "oligoset_id",
                "oligo_0",
                "oligo_1",
                "oligo_2",
                "set_score_lowest",
                "set_score_sum",
            ],
        )

        self.oligo_database.oligosets["region_1"] = oligoset

        folder_oligosets = self.oligo_database.write_oligosets_to_table()
        file_oligosets = os.path.join(folder_oligosets, "oligosets_region_1.tsv")
        assert (
            check_tsv_format(file=file_oligosets) == True
        ), f"error: incorrect file format of {file_oligosets}"

    def test_remove_regions_with_insufficient_oligos(self):
        self.load_database_from_fasta()
        self.oligo_database.remove_regions_with_insufficient_oligos("database_generation")
        assert len(self.oligo_database.database.keys()) == (
            len(REGION_IDS) - 1 + 1  # one region removed but one added from random seqs
        ), "error: wrong number of regions in database"

    def test_get_attribute_list(self):
        self.oligo_database.load_database_from_table(
            file_database=FILE_DATABASE_OLIGO_ATTRIBUTES,
            region_ids=None,
            database_overwrite=True,
        )

        list_attributes = self.oligo_database.get_attribute_list()
        print(list_attributes)
        assert len(list_attributes) == 13, "error: wrong number of attributes in database"
        assert "oligo" in list_attributes, "error: missing attribute"

    def test_get_oligoid_list(self):
        self.oligo_database.load_database_from_table(
            file_database=FILE_DATABASE_OLIGO_ATTRIBUTES,
            region_ids=None,
            database_overwrite=True,
        )

        list_oligoids = self.oligo_database.get_oligoid_list()
        assert len(list_oligoids) == 21, "error: wrong number of oligoids in database"
        assert "region_3::1" in list_oligoids, "error: missing oligoid"

    def test_get_sequence_list(self):
        self.oligo_database.load_database_from_table(
            file_database=FILE_DATABASE_OLIGO_ATTRIBUTES,
            region_ids=None,
            database_overwrite=True,
        )

        list_sequences = self.oligo_database.get_sequence_list(sequence_type="oligo")
        assert len(list_sequences) == 21, "error: wrong number of sequences in database"
        assert "TATAACCCTGAGGAGGTATACCTAG" in list_sequences, "error: missing sequence"

    def test_get_oligoid_sequence_mapping(self):
        self.oligo_database.load_database_from_table(
            file_database=FILE_DATABASE_OLIGO_ATTRIBUTES,
            region_ids=None,
            database_overwrite=True,
        )

        mapping = self.oligo_database.get_oligoid_sequence_mapping(sequence_type="oligo")
        assert mapping["region_1::1"] == "ATGCCCCAATGGATGACGAT", "error: wrong sequence for oligoid"
        assert mapping["region_3::5"] == "CTCACTCGACTCTTACACAGTCATA", "error: wrong sequence for oligoid"

        mapping = self.oligo_database.get_oligoid_sequence_mapping(sequence_type="target")
        assert mapping["region_1::1"] == "ATCGTCATCCATTGGGGCAT", "error: wrong sequence for oligoid"
        assert mapping["region_3::5"] == "TATGACTGTGTAAGAGTCGAGTGAG", "error: wrong sequence for oligoid"

    def test_get_sequence_oligoid_mapping(self):
        self.oligo_database.load_database_from_table(
            file_database=FILE_DATABASE_OLIGO_ATTRIBUTES,
            region_ids=None,
            database_overwrite=True,
        )

        mapping = self.oligo_database.get_sequence_oligoid_mapping(sequence_type="oligo")
        assert len(mapping["CTCACTCGACTCTTACACAGTCATA"]) == 4, "error: wrong number of oligos for sequence"

    def test_get_oligo_attribute_table(self):
        self.oligo_database.load_database_from_table(
            file_database=FILE_DATABASE_OLIGO_ATTRIBUTES,
            region_ids=None,
            database_overwrite=True,
        )
        attribute = self.oligo_database.get_oligo_attribute_table(attribute="test_attribute", flatten=True)

        assert len(attribute["test_attribute"].unique()) == 2, "error: wrong attribute returned"

    def test_get_oligo_attribute_value(self):
        self.oligo_database.load_database_from_table(
            file_database=FILE_DATABASE_OLIGO_ATTRIBUTES,
            region_ids=None,
            database_overwrite=True,
        )
        attribute1 = self.oligo_database.get_oligo_attribute_value(
            attribute="test_attribute", flatten=True, region_id="region_1", oligo_id="region_1::5"
        )
        attribute2 = self.oligo_database.get_oligo_attribute_value(
            attribute="test_attribute", flatten=False, region_id="region_3", oligo_id="region_3::3"
        )

        assert attribute1 == "red", "error: wrong attribute value returned"
        assert attribute2 == [["blue"]], "error: wrong attribute value returned"

    def test_update_oligo_attribute(self):
        self.oligo_database.load_database_from_table(
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
        attribute = self.oligo_database.get_oligo_attribute_table(attribute="GC_content", flatten=True)

        assert len(attribute) == 5, "error: attribute not correctly updated"

    def test_filter_database_by_region(self):
        self.oligo_database.load_database_from_table(
            file_database=FILE_DATABASE_OLIGO_ATTRIBUTES,
            region_ids=None,
            database_overwrite=True,
        )

        self.oligo_database.filter_database_by_region(remove_region=True, region_ids="region_3")

        assert len(self.oligo_database.database.keys()) == 2, "error: remove region was kept"

        self.oligo_database.filter_database_by_region(
            remove_region=False, region_ids=["region_1", "region_2"]
        )

        assert len(self.oligo_database.database.keys()) == 2, "error: keep regions were removed"

    def test_filter_database_by_attribute_threshold(self):
        self.oligo_database.load_database_from_table(
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

        self.oligo_database.filter_database_by_attribute_threshold(
            attribute_name="GC_content", attribute_thr=65, remove_if_smaller_threshold=True
        )
        assert len(self.oligo_database.get_oligoid_list()) == 3, "error: wrong number of oligos filtered"

        self.oligo_database.filter_database_by_attribute_threshold(
            attribute_name="GC_content", attribute_thr=75, remove_if_smaller_threshold=False
        )
        assert len(self.oligo_database.get_oligoid_list()) == 2, "error: wrong number of oligos filtered"

    def test_filter_database_by_attribute_category(self):
        self.oligo_database.load_database_from_table(
            file_database=FILE_DATABASE_OLIGO_ATTRIBUTES,
            region_ids=None,
            database_overwrite=True,
        )

        self.oligo_database.filter_database_by_attribute_category(
            attribute_name="test_attribute", attribute_category="red", remove_if_equals_category=True
        )
        assert len(self.oligo_database.get_oligoid_list()) == 9, "error: wrong number of oligos filtered"


class TestOligoAttributes(unittest.TestCase):
    def setUp(self):
        self.tmp_path = os.path.join(os.getcwd(), "tmp_oligo_attributes")

        self.oligo_database = OligoDatabase(
            min_oligos_per_region=2,
            write_regions_with_insufficient_oligos=True,
            database_name="test_oligo_attributes",
            dir_output=self.tmp_path,
        )
        self.oligo_database.load_database_from_table(FILE_DATABASE_OLIGO_ATTRIBUTES, database_overwrite=True)

        self.oligo_attributes = OligoAttributes()

    def tearDown(self):
        shutil.rmtree(self.tmp_path)

    def test_calculate_oligo_length(self):
        oligo_database = self.oligo_attributes.calculate_oligo_length(self.oligo_database)

        length1 = oligo_database.get_oligo_attribute_value(
            attribute="length", flatten=True, region_id="region_1", oligo_id="region_1::1"
        )
        length2 = oligo_database.get_oligo_attribute_value(
            attribute="length", flatten=True, region_id="region_1", oligo_id="region_1::2"
        )

        assert length1 == 20, "error: wrong oligo length"
        assert length2 == 29, "error: wrong oligo length"

    def test_calculate_num_targeted_transcripts(self):
        oligo_database = self.oligo_attributes.calculate_num_targeted_transcripts(self.oligo_database)

        num_targeted_transcripts1 = oligo_database.get_oligo_attribute_value(
            attribute="num_targeted_transcripts", flatten=True, region_id="region_1", oligo_id="region_1::1"
        )
        num_targeted_transcripts2 = oligo_database.get_oligo_attribute_value(
            attribute="num_targeted_transcripts", flatten=True, region_id="region_2", oligo_id="region_2::1"
        )
        num_targeted_transcripts3 = oligo_database.get_oligo_attribute_value(
            attribute="num_targeted_transcripts", flatten=True, region_id="region_3", oligo_id="region_3::4"
        )

        assert num_targeted_transcripts1 == 2, "error: wrong number targeted transcripts"
        assert num_targeted_transcripts2 == 1, "error: wrong number targeted transcripts"
        assert num_targeted_transcripts3 == 28, "error: wrong number targeted transcripts"

    def test_calculate_isoform_consensus(self):
        oligo_database = self.oligo_attributes.calculate_num_targeted_transcripts(self.oligo_database)
        oligo_database = self.oligo_attributes.calculate_isoform_consensus(oligo_database)

        isoform_consensus1 = oligo_database.get_oligo_attribute_value(
            attribute="isoform_consensus", flatten=True, region_id="region_1", oligo_id="region_1::1"
        )
        isoform_consensus2 = oligo_database.get_oligo_attribute_value(
            attribute="isoform_consensus", flatten=True, region_id="region_2", oligo_id="region_2::1"
        )

        assert isoform_consensus1 == 100, "error: wrong isoform consensus, should be 100%"
        assert isoform_consensus2 == 50, "error: wrong isoform consensus, should be 50%"

    def test_calculate_seedregion(self):
        oligo_database = self.oligo_attributes.calculate_seedregion(self.oligo_database, start=0.4, end=0.6)

        seedregion_start = oligo_database.get_oligo_attribute_value(
            attribute="seedregion_start", flatten=True, region_id="region_1", oligo_id="region_1::1"
        )
        seedregion_end = oligo_database.get_oligo_attribute_value(
            attribute="seedregion_end", flatten=True, region_id="region_1", oligo_id="region_1::1"
        )

        assert (seedregion_start == 8) and (seedregion_end == 12), "error: wrong seedregion calculated"

    def test_calculate_seedregion_ligationsite(self):
        oligo_database = self.oligo_attributes.calculate_seedregion_ligationsite(
            self.oligo_database, seedregion_size=5
        )

        seedregion_start = oligo_database.get_oligo_attribute_value(
            attribute="seedregion_start", flatten=True, region_id="region_1", oligo_id="region_1::1"
        )
        seedregion_end = oligo_database.get_oligo_attribute_value(
            attribute="seedregion_end", flatten=True, region_id="region_1", oligo_id="region_1::1"
        )

        assert (seedregion_start == 6) and (seedregion_end == 15), "error: wrong seedregion calculated"

    def test_calculate_GC_content(self):
        oligo_database = self.oligo_attributes.calculate_GC_content(
            self.oligo_database, sequence_type="oligo"
        )

        GC_content = oligo_database.get_oligo_attribute_value(
            attribute="GC_content", flatten=True, region_id="region_1", oligo_id="region_1::1"
        )

        assert GC_content == 50, "error: wrong GC content calculated"

    def test_calculate_TmNN(self):
        oligo_database = self.oligo_attributes.calculate_TmNN(
            self.oligo_database, sequence_type="oligo", Tm_parameters={}
        )

        TmNN = oligo_database.get_oligo_attribute_value(
            attribute="TmNN", flatten=True, region_id="region_1", oligo_id="region_1::1"
        )

        assert TmNN == 53.57, "error: wrong Tm calculated"

    def test_calculate_length_selfcomplement(self):
        oligo_database = self.oligo_attributes.calculate_length_selfcomplement(
            self.oligo_database, sequence_type="oligo"
        )

        length_selfcomplement = oligo_database.get_oligo_attribute_value(
            attribute="length_selfcomplement", flatten=True, region_id="region_3", oligo_id="region_3::3"
        )

        assert length_selfcomplement == 18, "error: wrong length of selfcomplement calculated"

    def test_calculate_length_complement(self):
        oligo_database = self.oligo_attributes.calculate_length_complement(
            self.oligo_database, comparison_sequence="AGTC", sequence_type="oligo"
        )

        length_complement = oligo_database.get_oligo_attribute_value(
            attribute="length_complement_AGTC", flatten=True, region_id="region_1", oligo_id="region_1::3"
        )

        assert length_complement == 3, "error: wrong length of complement calculated"

    def test_calculate_secondary_structure_DG(self):
        oligo_database = self.oligo_attributes.calculate_DG_secondary_structure(
            self.oligo_database, sequence_type="oligo", T=37
        )

        DG_secondary_structure = oligo_database.get_oligo_attribute_value(
            attribute="DG_secondary_structure", flatten=True, region_id="region_1", oligo_id="region_1::1"
        )

        assert DG_secondary_structure == 0.8, "error: wrong DG calculated"

    def test_calculate_padlock_arms(self):
        oligo_database = self.oligo_attributes.calculate_padlock_arms(
            self.oligo_database,
            arm_length_min=3,
            arm_Tm_dif_max=15,
            arm_Tm_min=30,
            arm_Tm_max=80,
            Tm_parameters={},
        )

        ligation_site = oligo_database.get_oligo_attribute_value(
            attribute="ligation_site", flatten=True, region_id="region_1", oligo_id="region_1::2"
        )

        assert ligation_site == 14, "error: wrong padlock arms calculated"

    def test_calculate_detect_oligo(self):
        oligo_database = self.oligo_attributes.calculate_detect_oligo(
            self.oligo_database,
            detect_oligo_length_min=8,
            detect_oligo_length_max=12,
            min_thymines=2,
        )

        detect_oligo_even = oligo_database.get_oligo_attribute_value(
            attribute="detect_oligo_even", flatten=True, region_id="region_1", oligo_id="region_1::2"
        )
        detect_oligo_long_left = oligo_database.get_oligo_attribute_value(
            attribute="detect_oligo_long_left", flatten=True, region_id="region_1", oligo_id="region_1::2"
        )
        detect_oligo_long_right = oligo_database.get_oligo_attribute_value(
            attribute="detect_oligo_long_right", flatten=True, region_id="region_1", oligo_id="region_1::2"
        )

        assert detect_oligo_even == "AGGGAATCGAAT", "error: wrong detection oligo even calculated"
        assert detect_oligo_long_left == None, "error: wrong detection oligo left calculated"
        assert detect_oligo_long_right == None, "error: wrong detection oligo right calculated"
