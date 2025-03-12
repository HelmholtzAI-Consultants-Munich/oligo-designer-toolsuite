############################################
# imports
############################################

import os
import shutil
import unittest
import warnings
import pandas as pd

from pathlib import Path
from effidict import LRUPickleDict

from oligo_designer_toolsuite.database import OligoDatabase
from oligo_designer_toolsuite.sequence_generator import OligoSequenceGenerator
from oligo_designer_toolsuite.utils import FastaParser, GffParser
from oligo_designer_toolsuite.utils import (
    check_if_dna_sequence,
    check_if_key_exists,
    check_if_list,
    check_if_list_of_lists,
    check_tsv_format,
    generate_unique_filename,
)
from oligo_designer_toolsuite.utils import (
    merge_databases,
    collapse_attributes_for_duplicated_sequences,
    format_oligo_attributes,
    check_if_region_in_database,
    flatten_attribute_list,
)

############################################
# Global Parameters
############################################

FILE_GFF = "tests/data/annotations/custom_GCF_000001405.40_GRCh38.p14_genomic_chr16.gff"
FILE_GTF = "tests/data/annotations/custom_GCF_000001405.40_GRCh38.p14_genomic_chr16.gtf"
FILE_FASTA = "tests/data/annotations/custom_GCF_000001405.40_GRCh38.p14_genomic_chr16.fna"
FILE_TSV = "tests/data/annotations/custom_GCF_000001405.40_GRCh38.p14_genomic_chr16.gtf.tsv"
FILE_PICKLE = "tests/data/annotations/custom_GCF_000001405.40_GRCh38.p14_genomic_chr16_gtf.pickle"
FILE_NCBI_EXONS = "tests/data/genomic_regions/sequences_ncbi_exons.fna"


############################################
# Tests
############################################
class TestCheckers(unittest.TestCase):
    def test_check_if_dna_sequence_valid(self):
        """Test if check_if_dna_sequence works correctly for a valid DNA sequence."""
        seq = "GGctAAgTTCCaGTttGCA"
        valid_characters = ["A", "C", "T", "G"]
        assert check_if_dna_sequence(seq, valid_characters), "error: check_if_dna_sequence failed"

    def test_check_if_dna_sequence_invalid(self):
        """Test if check_if_dna_sequence works correctly for an invalid DNA sequence."""
        seq = "GGctAAgTuuTCCaGTttGCA"
        valid_characters = ["A", "C", "T", "G", "W", "X"]
        assert not check_if_dna_sequence(
            seq, valid_characters
        ), "error: check_if_dna_sequence succeeded when it should have failed"

    def test_check_if_key_exists_empty(self):
        """Test the check_if_key_exists function with an empty cache."""
        empty_dict = LRUPickleDict()
        assert not check_if_key_exists(empty_dict, "a"), "Failed: Should return False for empty dictionary"

    def test_check_if_key_exists_flat(self):
        """Test the check_if_key_exists function with a flat cache."""
        flat_database = LRUPickleDict()
        flat_database.load_from_dict({"a": 1, "b": 2})

        assert check_if_key_exists(flat_database, "a"), "Failed: Key 'a' should exist in flat_database"
        assert not check_if_key_exists(
            flat_database, "z"
        ), "Failed: Key 'z' should not exist in flat_database"

    def test_check_if_key_exists_nested(self):
        """Test the check_if_key_exists function with a nested cache."""

        nested_database = LRUPickleDict()
        nested_database.load_from_dict({"a": {"b": {"c": 1}}, "d": 2, "e": {"f": {"g": {"h": 3}}}})

        assert check_if_key_exists(nested_database, "c"), "Failed: Key 'c' should exist in nested_database"
        assert not check_if_key_exists(
            nested_database, "z"
        ), "Failed: Key 'z' should not exist in nested_database"
        assert check_if_key_exists(
            nested_database, "h"
        ), "Failed: Key 'h' should exist deep within nested_database"

    def test_check_if_list_str(self):
        """Test if check_if_list works correctly for a string."""
        value = "test"
        result = check_if_list(value)
        assert result == [value], f"error: check_if_list failed. Expected: [{value}], got: {result}"

    def test_check_if_list_list(self):
        """Test if check_if_list works correctly for a list."""
        value = ["test", ["test2"]]
        result = check_if_list(value)
        assert result == value, f"error: check_if_list failed. Expected: {value}, got: {result}"

    def test_check_if_list_of_lists_str(self):
        """Test if check_if_list works correctly for a string."""
        value = "test"
        result = check_if_list_of_lists(value)
        assert result == [[value]], f"error: check_if_list failed. Expected: [[{value}]], got: {result}"

    def test_check_if_list_of_lists_list(self):
        """Test if check_if_list works correctly for a list."""
        value = ["test", ["test2"]]
        result = check_if_list_of_lists(value)
        assert result == [value], f"error: check_if_list failed. Expected: [{value}], got: {result}"

    def test_check_if_list_of_lists_list_of_lists(self):
        """Test if check_if_list works correctly for a list."""
        value = [["test", ["test2"]]]
        result = check_if_list_of_lists(value)
        assert result == value, f"error: check_if_list failed. Expected: {value}, got: {result}"

    def test_check_tsv_format(self):
        """Test if the parser extracts fasta header correctly."""
        try:
            check_tsv_format(FILE_TSV)
        except Exception as e:
            assert False, f"error: checker: check_tsv_format raised an exception: {e}, with file {FILE_TSV}"

    def test_generate_unique_filename(self):
        dir_output = "dir_test"
        base_name = "testfile"
        extension = ".txt"
        file = generate_unique_filename(dir_output=dir_output, base_name=base_name, extension=extension)
        assert os.path.splitext(file)[1] == extension, "error: wrong file extension"
        assert base_name in os.path.splitext(os.path.basename(file))[0], "error: basename not in file"
        assert os.path.dirname(file) == dir_output, "error: wrong directory in file"


class TestDatabaseProcessor(unittest.TestCase):
    def setUp(self):
        self.tmp_path = os.path.join(os.getcwd(), "tmp_oligo_sequence_generator")
        Path(self.tmp_path).mkdir(parents=True, exist_ok=True)

        self.oligo_sequence_generator = OligoSequenceGenerator()
        file_fasta_exons = self.oligo_sequence_generator.create_sequences_sliding_window(
            files_fasta_in=FILE_NCBI_EXONS,
            length_interval_sequences=(100, 100),
            region_ids="AARS1",
        )

        # create two database with identical entries
        self.oligo_database1 = OligoDatabase()
        self.oligo_database1.load_database_from_fasta(
            files_fasta=file_fasta_exons,
            database_overwrite=True,
            sequence_type="oligo",
            region_ids="AARS1",
        )

        self.oligo_database1.filter_database_by_oligo(
            remove_region=False, oligo_ids=["AARS1::1", "AARS1::2", "AARS1::3"]
        )

        self.oligo_database2 = OligoDatabase()
        self.oligo_database2.load_database_from_fasta(
            files_fasta=file_fasta_exons,
            database_overwrite=True,
            sequence_type="oligo",
            region_ids="AARS1",
        )

        self.oligo_database2.filter_database_by_oligo(
            remove_region=False, oligo_ids=["AARS1::1", "AARS1::2", "AARS1::5"]
        )
        self.oligo_database2.update_oligo_attributes({"AARS1::2": {"start": 70265560, "end": 70265660}})

    def tearDown(self):
        shutil.rmtree(self.tmp_path)

    def test_merge_databases(self):
        oligo_database_merged = merge_databases(
            self.oligo_database1.database,
            self.oligo_database2.database,
            sequence_type="oligo",
            dir_cache_files=self.oligo_database1._dir_cache_files,
            lru_db_max_in_memory=self.oligo_database1.lru_db_max_in_memory,
        )

        assert len(oligo_database_merged["AARS1"]) == 4, "error: region not succesfully merged"
        assert oligo_database_merged["AARS1"]["AARS1::1"]["start"] == [
            [70265563]
        ], "error: attributes incorrectly merged"
        assert oligo_database_merged["AARS1"]["AARS1::2"]["start"] == [
            [70265562],
            [70265560],
        ], "error: attributes incorrectly merged"

    def test_collapse_attributes_for_duplicated_sequences_dict_identical(self):
        dict1 = {"chromosome": [["10"]], "start": [[1000]], "end": [[2000]], "strand": [["+"]]}
        dict2 = {"chromosome": [["10"]], "start": [[1000]], "end": [[2000]], "strand": [["+"]]}

        dict_merged = collapse_attributes_for_duplicated_sequences(dict1, dict2)

        assert dict_merged == dict1, "error: identical dict should not have duplicated elements"

    def test_collapse_attributes_for_duplicated_sequences_dict_different(self):
        dict1 = {"chromosome": [["10"]], "start": [[1000]], "end": [[2000]], "strand": [["+"]]}
        dict2 = {"chromosome": [["11"]], "start": [[1020]], "end": [[2020]], "strand": [["-"]]}

        dict_merged = collapse_attributes_for_duplicated_sequences(dict1, dict2)

        assert dict_merged["chromosome"] == [["10"], ["11"]], "error: different dicts should have been merged"
        assert dict_merged["start"] == [[1000], [1020]], "error: different dicts should have been merged"
        assert dict_merged["end"] == [[2000], [2020]], "error: different dicts should have been merged"
        assert dict_merged["strand"] == [["+"], ["-"]], "error: different dicts should have been merged"

    def test_format_oligo_attributes(self):
        oligo_attributes = {"chromosome": "10", "start": [1000], "end": [[2000]], "strand": [["+"], ["-"]]}

        oligo_attributes = format_oligo_attributes(oligo_attributes=oligo_attributes)

        assert oligo_attributes["chromosome"] == [["10"]], "error: oligo attribute not correctly formatted"
        assert oligo_attributes["start"] == [[1000]], "error: oligo attribute not correctly formatted"
        assert oligo_attributes["end"] == [[2000]], "error: oligo attribute not correctly formatted"
        assert oligo_attributes["strand"] == [["+"], ["-"]], "error: oligo attribute not correctly formatted"

    def test_check_if_region_in_database(self):
        file_removed_regions = os.path.join(self.tmp_path, "removed_regions.tsv")
        check_if_region_in_database(
            database=self.oligo_database1.database,
            region_ids=["no_region1", "no_region2"],
            write_regions_with_insufficient_oligos=True,
            file_removed_regions=file_removed_regions,
        )

        removed_regions = pd.read_csv(
            filepath_or_buffer=file_removed_regions, sep="\t", header=None, names=["region", "step"]
        )
        assert removed_regions.region[0] == "no_region1", "error: region was not removed"

    def test_flatten_attribute_list(self):
        oligo_attributes = {
            "chromosome": [["10"]],
            "start": [[1000], [1020]],
            "end": [[2000]],
            "strand": [["+", "+"], ["-"]],
        }

        assert flatten_attribute_list(oligo_attributes["chromosome"]) == [
            "10"
        ], "error: attribute not flattened"
        assert flatten_attribute_list(oligo_attributes["strand"]) == [
            "+",
            "+",
            "-",
        ], "error: attribute not flattened"
        assert flatten_attribute_list(oligo_attributes["start"]) == [
            1000,
            1020,
        ], "error: attribute not flattened"


class TestGffParser(unittest.TestCase):
    def setUp(self):
        self.parser = GffParser()

    def test_check_gff_format(self):
        """Test parsing GFF annotation data."""
        try:
            self.parser.check_gff_format(FILE_GFF)
        except Exception as e:
            assert False, f"error: checker: check_gff_format raised an exception: {e}, with file {FILE_GFF}"

    def test_check_gtf_format(self):
        """Test parsing GTF annotation data."""
        try:
            self.parser.check_gff_format(FILE_GTF)
        except Exception as e:
            assert False, f"error: checker: check_gff_format raised an exception: {e}, with file {FILE_GTF}"

    def test_parse_annotation_from_gff(self):
        """Test parsing GFF annotation."""
        result = self.parser.parse_annotation_from_gff(FILE_GFF, target_lines=10)
        assert result.shape[1] == 23, "error: GFF3 dataframe not correctly loaded"

    def test_parse_annotation_from_gtf(self):
        """Test parsing GTF annotation."""
        result = self.parser.parse_annotation_from_gff(FILE_GTF, target_lines=10)
        assert result.shape[1] == 20, "error: GTF dataframe not correctly loaded"

    def test_load_annotation_from_pickle_file(self):
        """Test loading annotation from a pickle file."""
        result = self.parser.load_annotation_from_pickle(FILE_PICKLE)
        assert type(result) == pd.DataFrame, f"error: GTF dataframe not correctly loaded from pickle file"


class TestFastaParser(unittest.TestCase):
    def setUp(self):
        self.parser = FastaParser()

    def test_check_fasta_format(self):
        """Test parsing fasta file."""
        try:
            self.parser.check_fasta_format(FILE_FASTA)
        except Exception as e:
            assert (
                False
            ), f"error: checker: check_fasta_format raised an exception: {e}, with file {FILE_FASTA}"

    def test_read_fasta_sequences_existing_regions(self):
        """Test parsing fasta file."""
        ids = ["16"]
        result = self.parser.read_fasta_sequences(FILE_FASTA, region_ids=ids)

        assert len(result) == 1, f"error: the function loaded {len(result)} entries instead of 1"
        assert result[0].name == "16", f"error: the name should be '16' instead of {result[0].name}"
        assert (
            result[0].description == "16 Homo sapiens chromosome 16, GRCh38.p14 Primary Assembly"
        ), f"error: the description should be '16 Homo sapiens chromosome 16, GRCh38.p14 Primary Assembly' instead of {result[0].description}"
        assert (
            result[0].dbxrefs == []
        ), f"error: the dbxrefs should be an empty list instead of {result[0].dbxrefs}"

    def test_get_fasta_regions(self):
        """Test if the parser extracts fasta regions correctly."""
        expected_result = ["16"]
        result = self.parser.get_fasta_regions(FILE_FASTA)
        assert (
            result == expected_result
        ), f"error: fasta regions not correctly extracted. Expected ['16'] got {result}"

    def test_parse_fasta_header(self):
        """Test if the parser extracts fasta header correctly."""
        header = "ARPG3::transcript_id=XM4581;exon_id=XM4581_exon1::16:70265537-70265662(-)"
        region, additional_information, coordinates = self.parser.parse_fasta_header(header)
        assert region == "ARPG3", "error: wrong region parsed"
        assert coordinates["chromosome"] == ["16"], "error: wrong chrom parsed"
        assert coordinates["start"] == [70265537], "error: wrong start parsed"
        assert coordinates["end"] == [70265662], "error: wrong end parsed"
        assert coordinates["strand"] == ["-"], "error: wrong strand parsed"
        assert additional_information == {
            "transcript_id": ["XM4581"],
            "exon_id": ["XM4581_exon1"],
        }, f"error: wrong additional information parsed: {additional_information}"
