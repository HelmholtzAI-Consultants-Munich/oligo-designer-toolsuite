############################################
# imports
############################################

import unittest
import warnings

import pandas as pd
from effidict import LRUPickleDict

from oligo_designer_toolsuite.utils import FastaParser, GffParser
from oligo_designer_toolsuite.utils._checkers_and_helpers import (
    check_if_dna_sequence,
    check_if_key_exists,
    check_if_list,
    check_tsv_format,
)

############################################
# Global Parameters
############################################

FILE_GFF = "tests/data/annotations/custom_GCF_000001405.40_GRCh38.p14_genomic_chr16.gff"
FILE_GTF = "tests/data/annotations/custom_GCF_000001405.40_GRCh38.p14_genomic_chr16.gtf"
FILE_FASTA = "tests/data/annotations/custom_GCF_000001405.40_GRCh38.p14_genomic_chr16.fna"
FILE_TSV = "tests/data/annotations/custom_GCF_000001405.40_GRCh38.p14_genomic_chr16.gtf.tsv"
FILE_PICKLE = "tests/data/annotations/custom_GCF_000001405.40_GRCh38.p14_genomic_chr16_gtf.pickle"

############################################
# Tests
############################################


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

    def test_read_fasta_sequences_non_existing_region(self):
        """Test parsing fasta file."""
        ids = ["1"]
        result = self.parser.read_fasta_sequences(FILE_FASTA, region_ids=ids)
        # check if a warning was raised
        with warnings.catch_warnings(record=True) as w:
            result = self.parser.read_fasta_sequences(FILE_FASTA, region_ids=ids)
            assert len(w) > 0, "error: no warning was raised"
        assert len(result) == 0, f"error: the function loaded {len(result)} entries instead of an empty list"

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


class TestCheckers(unittest.TestCase):
    def test_check_tsv_format(self):
        """Test if the parser extracts fasta header correctly."""
        try:
            check_tsv_format(FILE_TSV)
        except Exception as e:
            assert False, f"error: checker: check_tsv_format raised an exception: {e}, with file {FILE_TSV}"

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
        s = "test"
        result = check_if_list(s)
        assert result == [s], f"error: check_if_list failed. Expected: [{s}], got: {result}"

    def test_check_if_list_list(self):
        """Test if check_if_list works correctly for a list."""
        l = ["test", ["test2"]]
        result = check_if_list(l)
        assert result == l, f"error: check_if_list failed. Expected: {l}, got: {result}"

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
