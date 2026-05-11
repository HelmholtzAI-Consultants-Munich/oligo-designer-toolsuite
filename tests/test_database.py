############################################
# imports
############################################

import os
import shutil
import unittest

import pandas as pd
import yaml

from oligo_designer_toolsuite.database import OligoDatabase, ReferenceDatabase
from oligo_designer_toolsuite.sequence_generator import OligoSequenceGenerator
from oligo_designer_toolsuite.utils import FastaParser, VCFParser, check_tsv_format

############################################
# setup
############################################

# Global Parameters
FILE_NCBI_EXONS = "tests/data/genomic_regions/sequences_ncbi_exons.fna"
FILE_VARIANTS = "tests/data/annotations/custom_GCF_000001405.40.chr16.vcf"
FILE_DATABASE_OLIGO_PROPERTIES = "tests/data/databases/database_oligo_properties.tsv"

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
    def setUp(self) -> None:
        self.tmp_path = os.path.join(os.getcwd(), "tmp_reference_database")

        self.fasta_parser = FastaParser()
        self.vcf_parser = VCFParser()

        self.reference_fasta = ReferenceDatabase(
            database_name="test_reference_database_fasta", dir_output=self.tmp_path
        )
        self.reference_fasta.load_database_from_file(
            files=[FILE_NCBI_EXONS], file_type="fasta", database_overwrite=False
        )
        self.reference_fasta.load_database_from_file(
            files=FILE_NCBI_EXONS, file_type="fasta", database_overwrite=True
        )

        self.reference_vcf = ReferenceDatabase(
            database_name="test_reference_database_vcf", dir_output=self.tmp_path
        )
        self.reference_vcf.load_database_from_file(
            files=FILE_VARIANTS, file_type="vcf", database_overwrite=True
        )

    def tearDown(self) -> None:
        try:
            shutil.rmtree(self.tmp_path)
        except:
            pass

    def test_write_database(self) -> None:
        file_fasta_database = self.reference_fasta.write_database_to_file(filename="ref_db_filtered_fasta")
        assert (
            self.fasta_parser.check_fasta_format(file_fasta_database) == True
        ), f"error: wrong file format for database in {file_fasta_database}"

        # Verify that .fai index file is created when writing FASTA database
        index_file = f"{file_fasta_database}.fai"
        assert os.path.exists(index_file), f"error: .fai index file should be created at {index_file}"

        file_vcf_database = self.reference_vcf.write_database_to_file(filename="ref_db_filtered_vcf")
        assert (
            self.vcf_parser.check_vcf_format(file_vcf_database) == True
        ), f"error: wrong file format for database in {file_vcf_database}"

    def test_load_database_from_file(self) -> None:
        """Test that index files are removed when database_overwrite=True."""
        # Create a temporary FASTA database with index file
        file_fasta_new = os.path.join(self.tmp_path, "tmp_new_fasta.fna")
        file_index_old = f"{self.reference_fasta.database_file}.fai"

        with open(file_fasta_new, "w") as f:
            f.write(">test_seq::test_info::chr1:100-200(+)\nATCGATCGATCG\n")
        with open(file_index_old, "w") as f:
            f.write("dummy index content\n")

        self.reference_fasta.load_database_from_file(
            files=file_fasta_new, file_type="fasta", database_overwrite=True
        )

        # Verify that the old index file was removed
        assert not os.path.exists(
            file_index_old
        ), "error: old .fai index file should have been removed when database_overwrite=True"

    def test_filter_database_by_region(self) -> None:
        self.reference_fasta.filter_database_by_region(region_ids="AARS1", keep_region=False)
        assert self.reference_fasta.database_file is not None, "error: database file is not set"
        fasta_sequences = self.fasta_parser.read_fasta_sequences(
            file_fasta_in=self.reference_fasta.database_file
        )
        for entry in fasta_sequences:
            (
                region,
                _,
                _,
            ) = self.fasta_parser.parse_fasta_header(entry.id)
            assert region != "AARS1", f"error: this region {region} should be filtered out."

    def test_filter_database_by_property_category(self) -> None:
        self.reference_fasta.filter_database_by_property_category(
            property_name="gene_id",
            property_category="AARS1",
            keep_if_equals_category=False,
        )
        assert self.reference_fasta.database_file is not None, "error: database file is not set"
        fasta_sequences = self.fasta_parser.read_fasta_sequences(
            file_fasta_in=self.reference_fasta.database_file
        )
        for entry in fasta_sequences:
            (
                region,
                _,
                _,
            ) = self.fasta_parser.parse_fasta_header(entry.id)
            assert region != "AARS1", f"error: this region {region} should be filtered out."


class TestOligoDatabase(unittest.TestCase):
    def setUp(self) -> None:
        self.tmp_path = os.path.join(os.getcwd(), "tmp_oligo_database")

        self.fasta_parser = FastaParser()

        self.oligo_sequence_generator = OligoSequenceGenerator(dir_output=self.tmp_path)
        self.oligo_database = OligoDatabase(
            min_oligos_per_region=2,
            write_regions_with_insufficient_oligos=True,
            max_entries_in_memory=10,
            n_jobs=4,
            database_name="test_oligo_database",
            dir_output=self.tmp_path,
        )
        self.oligo_database.set_database_sequence_types(["oligo", "target"])

    def tearDown(self) -> None:
        shutil.rmtree(self.tmp_path)

    def load_database_from_fasta(self) -> None:
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

    def test_load_database_from_fasta(self) -> None:
        self.load_database_from_fasta()

        assert len(self.oligo_database.database) == 6, "error: wrong number of sequences loaded into database"
        # Verify that sequence types are set correctly
        assert (
            "oligo" in self.oligo_database.database_sequence_types
        ), "error: 'oligo' should be in database_sequence_types"
        assert (
            "target" in self.oligo_database.database_sequence_types
        ), "error: 'target' should be in database_sequence_types"

    def test_load_database_from_table(self) -> None:
        self.oligo_database.load_database_from_table(
            FILE_DATABASE_OLIGO_PROPERTIES,
            region_ids=["region_1", "region_2"],
            database_overwrite=True,
            merge_databases_on_sequence_type="oligo",
        )

        assert len(self.oligo_database.database) == 2, "error: wrong number of sequences loaded into database"
        # Verify that sequence type is set when loading from table
        assert (
            "oligo" in self.oligo_database.database_sequence_types
        ), "error: 'oligo' should be in database_sequence_types after loading from table"

    def test_load_save_database(self) -> None:
        self.oligo_database.load_database_from_table(
            FILE_DATABASE_OLIGO_PROPERTIES, database_overwrite=True, merge_databases_on_sequence_type="oligo"
        )

        dir_database = self.oligo_database.save_database(
            region_ids=["region_1", "region_2"], name_database="database_region1_region2"
        )
        self.oligo_database.load_database(
            dir_database, database_overwrite=True, merge_databases_on_sequence_type="oligo"
        )

        assert len(self.oligo_database.database.keys()) == 2, "error: wrong number regions saved and loaded"

    def test_write_database_to_fasta(self) -> None:
        self.oligo_database.load_database_from_table(
            FILE_DATABASE_OLIGO_PROPERTIES,
            region_ids=["region_1", "region_2"],
            database_overwrite=True,
            merge_databases_on_sequence_type="oligo",
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

    def test_write_database_to_table(self) -> None:
        self.oligo_database.load_database_from_table(
            FILE_DATABASE_OLIGO_PROPERTIES, database_overwrite=True, merge_databases_on_sequence_type="oligo"
        )

        file_database = self.oligo_database.write_database_to_table(
            properties=["test_property", "ligation_site", "chromosome", "start", "end", "strand"],
            flatten_property=True,
            filename="database_region1_region2_flattened",
            region_ids=["region_1", "region_2"],
        )

        self.oligo_database.load_database_from_table(
            file_database, database_overwrite=True, merge_databases_on_sequence_type="oligo"
        )

        assert check_tsv_format(file_database) == True, f"error: wrong file format"
        assert len(self.oligo_database.database.keys()) == 2, "error: wrong number regions saved and loaded"
        assert (
            self.oligo_database.get_oligo_property_value(
                property="test_property", flatten=True, region_id="region_1", oligo_id="region_1::1"
            )
            == "red"
        ), f"error: wrong property stored in {file_database}"

        file_database = self.oligo_database.write_database_to_table(
            properties=["test_property", "ligation_site", "chromosome", "start", "end", "strand"],
            flatten_property=False,
            filename="database_region1_region2_unflattened",
            region_ids=["region_1", "region_2"],
        )

        self.oligo_database.load_database_from_table(
            file_database, database_overwrite=True, merge_databases_on_sequence_type="oligo"
        )

        assert check_tsv_format(file_database) == True, f"error: wrong file format"
        assert len(self.oligo_database.database.keys()) == 2, "error: wrong number regions saved and loaded"
        assert (
            self.oligo_database.get_oligo_property_value(
                property="test_property", flatten=True, region_id="region_1", oligo_id="region_1::1"
            )
            == "red"
        ), f"error: wrong property stored in {file_database}"

    def test_write_database_to_bed(self) -> None:
        self.oligo_database.load_database_from_table(
            FILE_DATABASE_OLIGO_PROPERTIES,
            region_ids=["region_1", "region_2"],
            database_overwrite=True,
            merge_databases_on_sequence_type="oligo",
        )

        file_bed = self.oligo_database.write_database_to_bed(
            filename="database_region1_region2_bed", region_ids=["region_1", "region_2"]
        )

        assert check_tsv_format(file_bed) == True, f"error: wrong file format"

        bed_table = pd.read_csv(
            file_bed, sep="\t", names=["chromosome", "start", "end", "name", "score", "strand"]
        )

        assert bed_table.loc[0, "start"] == 70289456
        assert bed_table.loc[0, "end"] == 70289485

    def test_write_oligosets_to_yaml(self) -> None:
        self.oligo_database.load_database_from_table(
            FILE_DATABASE_OLIGO_PROPERTIES,
            database_overwrite=True,
            region_ids="region_1",
            merge_databases_on_sequence_type="oligo",
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

        self.oligo_database.write_oligosets_to_yaml(
            properties=[
                "test_property",
                "ligation_site",
                "chromosome",
                "start",
                "end",
                "strand",
                "transcript_id",
            ],
            ascending=True,
        )

        # Construct file path since function returns None
        file_yaml = os.path.join(os.path.dirname(self.oligo_database.dir_output), "oligosets.yml")
        with open(file_yaml, "r") as handle:
            yaml_oligosets = yaml.safe_load(handle)

        assert yaml_oligosets["region_1"]["Oligoset 1"]["Oligoset Score"] == {
            "set_score_lowest": 1.59,
            "set_score_sum": 2.36,
        }, f"error: wrong oligoset loaded"

        assert yaml_oligosets["region_1"]["Oligoset 1"]["Oligo 1"]["test_property"] == [
            "red"
        ], f"error: wrong oligoset loaded"

        assert yaml_oligosets["region_1"]["Oligoset 1"]["Oligo 1"]["transcript_id"] == [
            ["NM_001605.3"],
            ["XM_047433666.1"],
        ], f"error: wrong oligoset loaded"

    def test_write_ready_to_order_yaml(self) -> None:
        self.oligo_database.load_database_from_table(
            FILE_DATABASE_OLIGO_PROPERTIES,
            database_overwrite=True,
            region_ids="region_1",
            merge_databases_on_sequence_type="oligo",
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

        self.oligo_database.write_ready_to_order_yaml(
            properties=["oligo", "target", "test_property"],
            ascending=True,
            filename="test_ready_to_order",
        )

        # Verify YAML file exists
        file_yaml = os.path.join(os.path.dirname(self.oligo_database.dir_output), "test_ready_to_order.yml")
        assert os.path.exists(file_yaml), f"error: YAML file {file_yaml} was not created"

        # Load and verify YAML structure
        with open(file_yaml, "r") as handle:
            yaml_order = yaml.safe_load(handle)

        # Verify structure: region_id -> oligoset_id -> oligo_id -> properties
        assert "region_1" in yaml_order, "error: region_1 should be in YAML"
        assert "oligoset_1" in yaml_order["region_1"], "error: oligoset_1 should be in region_1"
        assert "oligoset_2" in yaml_order["region_1"], "error: oligoset_2 should be in region_1"

        # Verify first oligoset contains expected oligos
        oligoset_1 = yaml_order["region_1"]["oligoset_1"]
        assert "region_1::1" in oligoset_1, "error: region_1::1 should be in oligoset_1"
        assert "region_1::2" in oligoset_1, "error: region_1::2 should be in oligoset_1"
        assert "region_1::5" in oligoset_1, "error: region_1::5 should be in oligoset_1"

    def test_write_oligosets_to_table(self) -> None:
        self.oligo_database.load_database_from_table(
            FILE_DATABASE_OLIGO_PROPERTIES,
            database_overwrite=True,
            region_ids="region_1",
            merge_databases_on_sequence_type="oligo",
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

        self.oligo_database.write_oligosets_to_table(
            properties=["oligo", "target", "test_property"],
            ascending=True,
        )

        # Construct file path since function returns None
        file_oligosets_tsv = os.path.join(os.path.dirname(self.oligo_database.dir_output), "oligosets.tsv")

        # Verify TSV file exists and has correct format
        assert os.path.exists(file_oligosets_tsv), f"error: TSV file {file_oligosets_tsv} was not created"
        assert (
            check_tsv_format(file=file_oligosets_tsv) == True
        ), f"error: incorrect file format of {file_oligosets_tsv}"

        # Verify Excel file exists (same directory, .xlsx extension)
        file_oligosets_excel = file_oligosets_tsv.replace(".tsv", ".xlsx")
        assert os.path.exists(
            file_oligosets_excel
        ), f"error: Excel file {file_oligosets_excel} was not created"

        # Verify Excel file structure (one sheet per region_id)
        try:
            excel_data = pd.read_excel(file_oligosets_excel, sheet_name=None, engine="openpyxl")
            # Check that at least one sheet exists (should be one for region_1)
            assert len(excel_data) > 0, "error: Excel file should contain at least one sheet"
            # Verify the sheet doesn't have region_id column (it should be removed)
            region_sheet = list(excel_data.values())[0]
            assert (
                "region_id" not in region_sheet.columns
            ), "error: region_id column should not be in individual Excel sheets"
            # Verify the sheet contains expected columns
            assert (
                "oligoset_id" in region_sheet.columns
            ), "error: Excel sheet should contain oligoset_id column"
            assert "oligo_id" in region_sheet.columns, "error: Excel sheet should contain oligo_id column"
            # Verify that specified properties are included
            assert "oligo" in region_sheet.columns, "error: Excel sheet should contain oligo property"
            assert "target" in region_sheet.columns, "error: Excel sheet should contain target property"
            assert "test_property" in region_sheet.columns, "error: Excel sheet should contain test_property"
        except ImportError:
            # Skip Excel verification if openpyxl is not installed
            pass

    def test_remove_regions_with_insufficient_oligos(self) -> None:
        self.load_database_from_fasta()
        self.oligo_database.remove_regions_with_insufficient_oligos("database_generation")
        assert len(self.oligo_database.database.keys()) == (
            len(REGION_IDS) - 1 + 1  # one region removed but one added from random seqs
        ), "error: wrong number of regions in database"

    def test_get_property_list(self) -> None:
        self.oligo_database.load_database_from_table(
            file_database=FILE_DATABASE_OLIGO_PROPERTIES,
            region_ids=None,
            database_overwrite=True,
            merge_databases_on_sequence_type="oligo",
        )

        list_properties = self.oligo_database.get_property_list()

        assert len(list_properties) == 13, "error: wrong number of properties in database"
        assert "oligo" in list_properties, "error: missing property"

    def test_get_regionid_list(self) -> None:
        self.oligo_database.load_database_from_table(
            file_database=FILE_DATABASE_OLIGO_PROPERTIES,
            region_ids=None,
            database_overwrite=True,
            merge_databases_on_sequence_type="oligo",
        )

        list_regionid = self.oligo_database.get_regionid_list()
        assert len(list_regionid) == 3, "error: wrong number of regionids in database"
        assert "region_3" in list_regionid, "error: missing regionid"

    def test_get_oligoid_list(self) -> None:
        self.oligo_database.load_database_from_table(
            file_database=FILE_DATABASE_OLIGO_PROPERTIES,
            region_ids=None,
            database_overwrite=True,
            merge_databases_on_sequence_type="oligo",
        )

        list_oligoids = self.oligo_database.get_oligoid_list()
        assert len(list_oligoids) == 21, "error: wrong number of oligoids in database"
        assert "region_3::1" in list_oligoids, "error: missing oligoid"

    def test_get_sequence_list(self) -> None:
        self.oligo_database.load_database_from_table(
            file_database=FILE_DATABASE_OLIGO_PROPERTIES,
            region_ids=None,
            database_overwrite=True,
            merge_databases_on_sequence_type="oligo",
        )

        list_sequences = self.oligo_database.get_sequence_list(sequence_type="oligo")
        assert len(list_sequences) == 21, "error: wrong number of sequences in database"
        assert "TATAACCCTGAGGAGGTATACCTAG" in list_sequences, "error: missing sequence"

    def test_get_oligoid_sequence_mapping(self) -> None:
        self.oligo_database.load_database_from_table(
            file_database=FILE_DATABASE_OLIGO_PROPERTIES,
            region_ids=None,
            database_overwrite=True,
            merge_databases_on_sequence_type="oligo",
        )

        mapping = self.oligo_database.get_oligoid_sequence_mapping(sequence_type="oligo")
        assert mapping["region_1::1"] == "ATGCCCCAATGGATGACGAT", "error: wrong sequence for oligoid"
        assert mapping["region_3::5"] == "CTCACTCGACTCTTACACAGTCATA", "error: wrong sequence for oligoid"

        mapping = self.oligo_database.get_oligoid_sequence_mapping(sequence_type="target")
        assert mapping["region_1::1"] == "ATCGTCATCCATTGGGGCAT", "error: wrong sequence for oligoid"
        assert mapping["region_3::5"] == "TATGACTGTGTAAGAGTCGAGTGAG", "error: wrong sequence for oligoid"

    def test_get_sequence_oligoid_mapping(self) -> None:
        self.oligo_database.load_database_from_table(
            file_database=FILE_DATABASE_OLIGO_PROPERTIES,
            region_ids=None,
            database_overwrite=True,
            merge_databases_on_sequence_type="oligo",
        )

        mapping = self.oligo_database.get_sequence_oligoid_mapping(sequence_type="oligo")
        assert len(mapping["CTCACTCGACTCTTACACAGTCATA"]) == 4, "error: wrong number of oligos for sequence"

    def test_get_oligo_property_table(self) -> None:
        self.oligo_database.load_database_from_table(
            file_database=FILE_DATABASE_OLIGO_PROPERTIES,
            region_ids=None,
            database_overwrite=True,
            merge_databases_on_sequence_type="oligo",
        )
        property_table = self.oligo_database.get_oligo_property_table(
            properties="test_property", flatten=True
        )

        assert (
            len(property_table.explode("test_property")["test_property"].unique()) == 2
        ), "error: wrong property returned"

    def test_get_oligo_property_value(self) -> None:
        self.oligo_database.load_database_from_table(
            file_database=FILE_DATABASE_OLIGO_PROPERTIES,
            region_ids=None,
            database_overwrite=True,
            merge_databases_on_sequence_type="oligo",
        )
        property1 = self.oligo_database.get_oligo_property_value(
            property="test_property", flatten=True, region_id="region_1", oligo_id="region_1::5"
        )
        property2 = self.oligo_database.get_oligo_property_value(
            property="test_property", flatten=False, region_id="region_3", oligo_id="region_3::3"
        )

        assert property1 == "red", "error: wrong property value returned"
        assert property2 == [["blue"]], "error: wrong property value returned"

    def test_update_oligo_properties(self) -> None:
        self.oligo_database.load_database_from_table(
            file_database=FILE_DATABASE_OLIGO_PROPERTIES,
            region_ids="region_3",
            database_overwrite=True,
            merge_databases_on_sequence_type="oligo",
        )
        new_property = {
            "region_3::1": {"GC_content": 63},
            "region_3::2": {"GC_content": 66},
            "region_3::3": {"GC_content": 80},
            "region_3::4": {"GC_content": 70},
            "region_3::5": {"GC_content": 40},
        }
        self.oligo_database.update_oligo_properties(new_oligo_property=new_property)
        property_table = self.oligo_database.get_oligo_property_table(properties="GC_content", flatten=True)

        assert len(property_table) == 5, "error: property not correctly updated"

    def test_filter_database_by_region(self) -> None:
        self.oligo_database.load_database_from_table(
            file_database=FILE_DATABASE_OLIGO_PROPERTIES,
            region_ids=None,
            database_overwrite=True,
            merge_databases_on_sequence_type="oligo",
        )

        self.oligo_database.filter_database_by_region(remove_region=True, region_ids="region_3")

        assert len(self.oligo_database.database.keys()) == 2, "error: remove region was kept"

        self.oligo_database.filter_database_by_region(
            remove_region=False, region_ids=["region_1", "region_2"]
        )

        assert len(self.oligo_database.database.keys()) == 2, "error: keep regions were removed"

    def test_filter_database_by_oligo(self) -> None:
        self.oligo_database.load_database_from_table(
            file_database=FILE_DATABASE_OLIGO_PROPERTIES,
            region_ids=None,
            database_overwrite=True,
            merge_databases_on_sequence_type="oligo",
        )

        self.oligo_database.filter_database_by_oligo(remove_region=True, oligo_ids="region_3::5")

        assert len(self.oligo_database.database["region_3"].keys()) == 4, "error: remove oligo was kept"

        self.oligo_database.filter_database_by_oligo(
            remove_region=False, oligo_ids=["region_3::3", "region_3::4"]
        )

        assert len(self.oligo_database.database["region_3"].keys()) == 2, "error: keep oligo were removed"

    def test_filter_database_by_property_threshold(self) -> None:
        self.oligo_database.load_database_from_table(
            file_database=FILE_DATABASE_OLIGO_PROPERTIES,
            region_ids="region_3",
            database_overwrite=True,
            merge_databases_on_sequence_type="oligo",
        )
        new_property = {
            "region_3::1": {"GC_content": 63},
            "region_3::2": {"GC_content": 66},
            "region_3::3": {"GC_content": 80},
            "region_3::4": {"GC_content": 70},
            "region_3::5": {"GC_content": 40},
        }
        self.oligo_database.update_oligo_properties(new_oligo_property=new_property)

        self.oligo_database.filter_database_by_property_threshold(
            property_name="GC_content", property_thr=65, remove_if_smaller_threshold=True
        )
        assert len(self.oligo_database.get_oligoid_list()) == 3, "error: wrong number of oligos filtered"

        self.oligo_database.filter_database_by_property_threshold(
            property_name="GC_content", property_thr=75, remove_if_smaller_threshold=False
        )
        assert len(self.oligo_database.get_oligoid_list()) == 2, "error: wrong number of oligos filtered"

    def test_filter_database_by_property_category(self) -> None:
        self.oligo_database.load_database_from_table(
            file_database=FILE_DATABASE_OLIGO_PROPERTIES,
            region_ids=None,
            database_overwrite=True,
            merge_databases_on_sequence_type="oligo",
        )

        self.oligo_database.filter_database_by_property_category(
            property_name="test_property", property_category="red", remove_if_equals_category=True
        )
        assert len(self.oligo_database.get_oligoid_list()) == 9, "error: wrong number of oligos filtered"

    def test_set_database_sequence_types(self) -> None:
        """Test that set_database_sequence_types correctly adds sequence types."""
        self.oligo_database.set_database_sequence_types("oligo")
        assert (
            "oligo" in self.oligo_database.database_sequence_types
        ), "error: 'oligo' should be added to database_sequence_types"

        self.oligo_database.set_database_sequence_types(["target", "oligo_pair_L"])
        assert (
            "target" in self.oligo_database.database_sequence_types
        ), "error: 'target' should be added to database_sequence_types"
        assert (
            "oligo_pair_L" in self.oligo_database.database_sequence_types
        ), "error: 'oligo_pair_L' should be added to database_sequence_types"
        assert (
            "oligo" in self.oligo_database.database_sequence_types
        ), "error: 'oligo' should still be in database_sequence_types"

        # Test that duplicates are not added
        initial_length = len(self.oligo_database.database_sequence_types)
        self.oligo_database.set_database_sequence_types("oligo")
        assert (
            len(self.oligo_database.database_sequence_types) == initial_length
        ), "error: duplicate sequence types should not be added"
