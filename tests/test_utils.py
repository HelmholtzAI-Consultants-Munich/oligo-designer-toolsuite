############################################
# imports
############################################

import os
import shutil
import unittest
from pathlib import Path
from typing import Any

import pandas as pd
from Bio.SeqUtils import MeltingTemp as mt
from effidict import EffiDict, LRUReplacement, PickleBackend

from oligo_designer_toolsuite.database import OligoDatabase
from oligo_designer_toolsuite.pipelines._utils import (
    get_highly_abundant_kmer_sequences,
    preprocess_tm_parameters,
)
from oligo_designer_toolsuite.sequence_generator import OligoSequenceGenerator
from oligo_designer_toolsuite.utils import (
    FastaParser,
    GffParser,
    VCFParser,
    cast_to_list,
    cast_to_list_of_lists,
    check_if_dna_sequence,
    check_if_key_exists,
    check_if_region_in_database,
    check_tsv_format,
    collapse_properties_for_duplicated_sequences,
    count_kmer_abundance,
    flatten_property_list,
    format_oligo_properties,
    generate_unique_filename,
    merge_databases,
    remove_index_files,
)

############################################
# Global Parameters
############################################

FILE_GFF = "tests/data/annotations/custom_GCF_000001405.40_GRCh38.p14_genomic_chr16.gff"
FILE_GTF = "tests/data/annotations/custom_GCF_000001405.40_GRCh38.p14_genomic_chr16.gtf"
FILE_GTF_COMPLEX = "tests/data/annotations/custom_GCF_000001405.40_GRCh38.p14_genomic_chr16_KI270728-1.gtf"
FILE_FASTA = "tests/data/annotations/custom_GCF_000001405.40_GRCh38.p14_genomic_chr16.fna"
FILE_VCF = "tests/data/annotations/custom_GCF_000001405.40.chr16.vcf"
FILE_TSV = "tests/data/annotations/custom_GCF_000001405.40_GRCh38.p14_genomic_chr16.gtf.tsv"
FILE_PICKLE = "tests/data/annotations/custom_GCF_000001405.40_GRCh38.p14_genomic_chr16_gtf.pickle"
FILE_NCBI_EXONS = "tests/data/genomic_regions/sequences_ncbi_exons.fna"


############################################
# Tests
############################################
class TestPreprocessTmParameters(unittest.TestCase):
    def test_converts_tables_from_strings(self) -> None:
        tm_parameters = {
            "nn_table": "DNA_NN3",
            "tmm_table": "DNA_TMM1",
            "imm_table": "DNA_IMM1",
            "de_table": "DNA_DE1",
            "Na": 1000,
        }

        processed = preprocess_tm_parameters(tm_parameters.copy())

        assert processed["nn_table"] == mt.DNA_NN3, "error: nn_table not converted to MeltingTemp table"
        assert processed["tmm_table"] == mt.DNA_TMM1, "error: tmm_table not converted to MeltingTemp table"
        assert processed["imm_table"] == mt.DNA_IMM1, "error: imm_table not converted to MeltingTemp table"
        assert processed["de_table"] == mt.DNA_DE1, "error: de_table not converted to MeltingTemp table"
        assert processed["Na"] == 1000, "error: non-table parameters should remain unchanged"

    def test_modifies_dictionary_in_place(self) -> None:
        tm_parameters = {
            "nn_table": "DNA_NN3",
            "tmm_table": "DNA_TMM1",
            "imm_table": "DNA_IMM1",
            "de_table": "DNA_DE1",
        }

        processed = preprocess_tm_parameters(tm_parameters)

        assert processed is tm_parameters, "error: tm_parameters should be modified in place"
        assert processed["nn_table"] == mt.DNA_NN3, "error: nn_table not converted in place"
        assert processed["tmm_table"] == mt.DNA_TMM1, "error: tmm_table not converted in place"
        assert processed["imm_table"] == mt.DNA_IMM1, "error: imm_table not converted in place"
        assert processed["de_table"] == mt.DNA_DE1, "error: de_table not converted in place"


class TestCheckers(unittest.TestCase):
    def setUp(self) -> None:
        self.tmp_path = os.path.join(os.getcwd(), "tmp_utils")
        os.makedirs(self.tmp_path, exist_ok=True)

        self._dir_cache_files = os.path.join(self.tmp_path, "cache_files")
        self.backend = PickleBackend(storage_path=self._dir_cache_files)
        self.strategy = LRUReplacement(disk_backend=self.backend, max_in_memory=100)

    def tearDown(self) -> None:
        shutil.rmtree(self.tmp_path)

    def test_check_if_dna_sequence_valid(self) -> None:
        """Test if check_if_dna_sequence works correctly for a valid DNA sequence."""
        seq = "GGctAAgTTCCaGTttGCA"
        valid_characters = ["A", "C", "T", "G"]
        assert check_if_dna_sequence(seq, valid_characters), "error: check_if_dna_sequence failed"

    def test_check_if_dna_sequence_invalid(self) -> None:
        """Test if check_if_dna_sequence works correctly for an invalid DNA sequence."""
        seq = "GGctAAgTuuTCCaGTttGCA"
        valid_characters = ["A", "C", "T", "G", "W", "X"]
        assert not check_if_dna_sequence(
            seq, valid_characters
        ), "error: check_if_dna_sequence succeeded when it should have failed"

    def test_check_if_key_exists_empty(self) -> None:
        """Test the check_if_key_exists function with an empty cache."""
        empty_dict = EffiDict(disk_backend=self.backend, replacement_strategy=self.strategy)
        assert not check_if_key_exists(empty_dict, "a"), "Failed: Should return False for empty dictionary"

    def test_check_if_key_exists_flat(self) -> None:
        """Test the check_if_key_exists function with a flat cache."""
        flat_database = EffiDict(disk_backend=self.backend, replacement_strategy=self.strategy)
        flat_database.load_from_dict({"a": 1, "b": 2})

        assert check_if_key_exists(flat_database, "a"), "Failed: Key 'a' should exist in flat_database"
        assert not check_if_key_exists(
            flat_database, "z"
        ), "Failed: Key 'z' should not exist in flat_database"

    def test_check_if_key_exists_nested(self) -> None:
        """Test the check_if_key_exists function with a nested cache."""

        nested_database = EffiDict(disk_backend=self.backend, replacement_strategy=self.strategy)
        nested_database.load_from_dict({"a": {"b": {"c": 1}}, "d": 2, "e": {"f": {"g": {"h": 3}}}})

        assert check_if_key_exists(nested_database, "c"), "Failed: Key 'c' should exist in nested_database"
        assert not check_if_key_exists(
            nested_database, "z"
        ), "Failed: Key 'z' should not exist in nested_database"
        assert check_if_key_exists(
            nested_database, "h"
        ), "Failed: Key 'h' should exist deep within nested_database"

    def test_cast_to_list_str(self) -> None:
        """Test if cast_to_list works correctly for a string."""
        value = "test"
        result = cast_to_list(value)
        assert result == [value], f"error: cast_to_list failed. Expected: [{value}], got: {result}"

    def test_cast_to_list_list(self) -> None:
        """Test if cast_to_list works correctly for a list."""
        value = ["test", ["test2"]]
        result = cast_to_list(value)
        assert result == value, f"error: cast_to_list failed. Expected: {value}, got: {result}"

    def test_cast_to_list_of_lists_str(self) -> None:
        """Test if cast_to_list_of_lists works correctly for a string."""
        value = "test"
        result = cast_to_list_of_lists(value)
        assert result == [
            [value]
        ], f"error: cast_to_list_of_lists failed. Expected: [[{value}]], got: {result}"

    def test_cast_to_list_of_lists_list(self) -> None:
        """Test if cast_to_list_of_lists works correctly for a list."""
        value = ["test", ["test2"]]
        result = cast_to_list_of_lists(value)
        assert result == [value], f"error: cast_to_list_of_lists failed. Expected: [{value}], got: {result}"

    def test_cast_to_list_of_lists_list_of_lists(self) -> None:
        """Test if cast_to_list_of_lists works correctly for a list."""
        value = [["test", ["test2"]]]
        result = cast_to_list_of_lists(value)
        assert result == value, f"error: cast_to_list_of_lists failed. Expected: {value}, got: {result}"

    def test_check_tsv_format(self) -> None:
        """Test if the parser extracts fasta header correctly."""
        try:
            check_tsv_format(FILE_TSV)
        except Exception as e:
            assert False, f"error: checker: check_tsv_format raised an exception: {e}, with file {FILE_TSV}"

    def test_generate_unique_filename(self) -> None:
        dir_output = "dir_test"
        base_name = "testfile"
        extension = ".txt"
        file = generate_unique_filename(dir_output=dir_output, base_name=base_name, extension=extension)
        assert os.path.splitext(file)[1] == extension, "error: wrong file extension"
        assert base_name in os.path.splitext(os.path.basename(file))[0], "error: basename not in file"
        assert os.path.dirname(file) == dir_output, "error: wrong directory in file"


class TestDatabaseProcessor(unittest.TestCase):
    def setUp(self) -> None:
        self.tmp_path = os.path.join(os.getcwd(), "tmp_oligo_sequence_generator")
        Path(self.tmp_path).mkdir(parents=True, exist_ok=True)

        self.oligo_sequence_generator = OligoSequenceGenerator(dir_output=self.tmp_path)
        file_fasta_exons = self.oligo_sequence_generator.create_sequences_sliding_window(
            files_fasta_in=FILE_NCBI_EXONS,
            length_interval_sequences=(100, 100),
            region_ids="AARS1",
        )

        # create two database with identical entries
        self.oligo_database1 = OligoDatabase(dir_output=self.tmp_path)
        self.oligo_database1.load_database_from_fasta(
            files_fasta=file_fasta_exons,
            database_overwrite=True,
            sequence_type="oligo",
            region_ids="AARS1",
        )

        self.oligo_database1.filter_database_by_oligo(
            remove_region=False, oligo_ids=["AARS1::1", "AARS1::2", "AARS1::3"]
        )

        self.oligo_database2 = OligoDatabase(dir_output=self.tmp_path)
        self.oligo_database2.load_database_from_fasta(
            files_fasta=file_fasta_exons,
            database_overwrite=True,
            sequence_type="oligo",
            region_ids="AARS1",
        )

        self.oligo_database2.filter_database_by_oligo(
            remove_region=False, oligo_ids=["AARS1::1", "AARS1::2", "AARS1::5"]
        )
        self.oligo_database2.update_oligo_properties({"AARS1::2": {"start": 70265560, "end": 70265660}})

    def tearDown(self) -> None:
        shutil.rmtree(self.tmp_path)

    def test_merge_databases(self) -> None:
        oligo_database_merged = merge_databases(
            self.oligo_database1.database,
            self.oligo_database2.database,
            sequence_type="oligo",
            database_sequence_types=["oligo"],
            dir_cache_files=self.oligo_database1._dir_cache_files,
            max_entries_in_memory=self.oligo_database1._max_entries_in_memory,
        )

        assert len(oligo_database_merged["AARS1"]) == 4, "error: region not succesfully merged"
        assert oligo_database_merged["AARS1"]["AARS1::1"]["start"] == [
            [70265563]
        ], "error: properties incorrectly merged"
        assert oligo_database_merged["AARS1"]["AARS1::2"]["start"] == [
            [70265562],
            [70265560],
        ], "error: properties incorrectly merged"

    def test_collapse_properties_for_duplicated_sequences_dict_identical(self) -> None:
        dict1 = {"chromosome": [["10"]], "start": [[1000]], "end": [[2000]], "strand": [["+"]]}
        dict2 = {"chromosome": [["10"]], "start": [[1000]], "end": [[2000]], "strand": [["+"]]}

        dict_merged = collapse_properties_for_duplicated_sequences(dict1, dict2, database_sequence_types=[])

        assert dict_merged == dict1, "error: identical dict should not have duplicated elements"

    def test_collapse_properties_for_duplicated_sequences_dict_different(self) -> None:
        dict1 = {"chromosome": [["10"]], "start": [[1000]], "end": [[2000]], "strand": [["+"]]}
        dict2 = {"chromosome": [["11"]], "start": [[1020]], "end": [[2020]], "strand": [["-"]]}

        dict_merged = collapse_properties_for_duplicated_sequences(dict1, dict2, database_sequence_types=[])

        assert dict_merged["chromosome"] == [["10"], ["11"]], "error: different dicts should have been merged"
        assert dict_merged["start"] == [[1000], [1020]], "error: different dicts should have been merged"
        assert dict_merged["end"] == [[2000], [2020]], "error: different dicts should have been merged"
        assert dict_merged["strand"] == [["+"], ["-"]], "error: different dicts should have been merged"

    def test_collapse_properties_for_duplicated_sequences_warns_on_different_sequence_values(self) -> None:
        """Test that a warning is issued when sequence type values differ between dictionaries."""
        dict1 = {"oligo": "ATCGATCGATCG", "chromosome": [["10"]], "start": [[1000]]}
        dict2 = {"oligo": "GCTAGCTAGCTA", "chromosome": [["10"]], "start": [[1000]]}

        with self.assertWarns(UserWarning) as warning_context:
            collapse_properties_for_duplicated_sequences(dict1, dict2, database_sequence_types=["oligo"])

        self.assertIn(
            "oligo",
            str(warning_context.warning),
            "error: warning should mention the key with different values",
        )

    def test_format_oligo_properties(self) -> None:
        oligo_properties = {"chromosome": "10", "start": [1000], "end": [[2000]], "strand": [["+"], ["-"]]}

        oligo_properties = format_oligo_properties(
            oligo_properties=oligo_properties, database_sequence_types=[]
        )

        assert oligo_properties["chromosome"] == [["10"]], "error: oligo property not correctly formatted"
        assert oligo_properties["start"] == [[1000]], "error: oligo property not correctly formatted"
        assert oligo_properties["end"] == [[2000]], "error: oligo property not correctly formatted"
        assert oligo_properties["strand"] == [["+"], ["-"]], "error: oligo property not correctly formatted"

    def test_collapse_properties_for_duplicated_sequences_with_sequence_type(self) -> None:
        """Test that sequence types are handled correctly in collapse_properties_for_duplicated_sequences."""
        dict1 = {"oligo": "ATCG", "chromosome": [["10"]], "start": [[1000]]}
        dict2 = {"oligo": "ATCG", "chromosome": [["11"]], "start": [[1020]]}

        dict_merged = collapse_properties_for_duplicated_sequences(
            dict1, dict2, database_sequence_types=["oligo"]
        )

        # Sequence type should warn if different, but here they're the same
        assert dict_merged["oligo"] == "ATCG", "error: sequence type should be preserved"
        assert dict_merged["chromosome"] == [
            ["10"],
            ["11"],
        ], "error: non-sequence properties should be merged"

    def test_format_oligo_properties_with_sequence_type(self) -> None:
        """Test that sequence types are not formatted in format_oligo_properties."""
        oligo_properties = {"oligo": "ATCG", "chromosome": "10", "start": [1000]}

        oligo_properties = format_oligo_properties(
            oligo_properties=oligo_properties, database_sequence_types=["oligo"]
        )

        assert oligo_properties["oligo"] == "ATCG", "error: sequence type should not be formatted"
        assert oligo_properties["chromosome"] == [["10"]], "error: non-sequence property should be formatted"
        assert oligo_properties["start"] == [[1000]], "error: non-sequence property should be formatted"

    def test_merge_databases_sequence_type_validation(self) -> None:
        """Test that merge_databases validates sequence_type is in database_sequence_types."""
        with self.assertRaises(ValueError):
            merge_databases(
                self.oligo_database1.database,
                self.oligo_database2.database,
                sequence_type="invalid_type",
                database_sequence_types=["oligo"],
                dir_cache_files=self.oligo_database1._dir_cache_files,
                max_entries_in_memory=self.oligo_database1._max_entries_in_memory,
            )

    def test_check_if_region_in_database(self) -> None:
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

    def test_flatten_property_list(self) -> None:
        oligo_properties: dict[str, list[list[Any]]] = {
            "chromosome": [["10"]],
            "start": [[1000], [1020]],
            "end": [[2000]],
            "strand": [["+", "+"], ["-"]],
        }

        assert flatten_property_list(oligo_properties["chromosome"]) == [
            "10"
        ], "error: property not flattened"
        assert flatten_property_list(oligo_properties["strand"]) == [
            "+",
            "+",
            "-",
        ], "error: property not flattened"
        assert flatten_property_list(oligo_properties["start"]) == [
            1000,
            1020,
        ], "error: property not flattened"


class TestGffParser(unittest.TestCase):
    def setUp(self) -> None:
        self.parser = GffParser()

    def test_check_gff_format(self) -> None:
        """Test parsing GFF annotation data."""
        try:
            self.parser.check_gff_format(FILE_GFF)
        except Exception as e:
            assert False, f"error: checker: check_gff_format raised an exception: {e}, with file {FILE_GFF}"

    def test_check_gtf_format(self) -> None:
        """Test parsing GTF annotation data."""
        try:
            self.parser.check_gff_format(FILE_GTF)
        except Exception as e:
            assert False, f"error: checker: check_gff_format raised an exception: {e}, with file {FILE_GTF}"

    def test_parse_annotation_from_gff(self) -> None:
        """Test parsing GFF annotation."""
        result: pd.DataFrame = self.parser.parse_annotation_from_gff(FILE_GFF, target_lines=10)
        assert result.shape[1] == 23, "error: GFF3 dataframe not correctly loaded"

    def test_parse_annotation_from_gtf(self) -> None:
        """Test parsing GTF annotation."""
        result: pd.DataFrame = self.parser.parse_annotation_from_gff(FILE_GTF, target_lines=10)
        assert result.shape[1] == 20, "error: GTF dataframe not correctly loaded"

    def test_parse_annotation_from_gtf_no_duplicates(self) -> None:
        """Test when parsing GTF annotation chromosomes are not read in as both integers and strings."""
        result: pd.DataFrame = self.parser.parse_annotation_from_gff(FILE_GTF_COMPLEX)
        self.assertListEqual(
            result["seqid"].unique().tolist(),
            ["16", "KI270728.1"],
            "error: GTF parsing does not generate unique chromosome values",
        )

    def test_load_annotation_from_pickle_file(self) -> None:
        """Test loading annotation from a pickle file."""
        result = self.parser.load_annotation_from_pickle(FILE_PICKLE)
        assert type(result) == pd.DataFrame, f"error: GTF dataframe not correctly loaded from pickle file"


class TestFastaParser(unittest.TestCase):
    def setUp(self) -> None:
        self.tmp_path = os.path.join(os.getcwd(), "tmp_fasta_parser")
        Path(self.tmp_path).mkdir(parents=True, exist_ok=True)

        self.parser = FastaParser()

    def tearDown(self) -> None:
        shutil.rmtree(self.tmp_path)

    def test_check_fasta_format(self) -> None:
        """Test parsing fasta file."""
        try:
            out = self.parser.check_fasta_format(FILE_FASTA)
            assert (
                out == True
            ), f"error: checker: check_fasta_format should have passed with file {FILE_FASTA}"
        except Exception as e:
            assert (
                False
            ), f"error: checker: check_fasta_format raised an exception: {e}, with file {FILE_FASTA}"

        try:
            out = self.parser.check_fasta_format(FILE_GFF)
            assert (
                out == False
            ), f"error: checker: check_fasta_format did not raise an exception with file {FILE_GFF}"
        except Exception:
            pass  # should go into this case

    def test_is_coordinate(self) -> None:
        """Test coordinate check with regular expression."""
        entry_true = "17:15-20(+)"
        entry_false1 = "17:15(-)"
        entry_false2 = "17:15-20"
        entry_false3 = "15-20(-)"

        assert (
            self.parser.is_coordinate(entry_true) == True
        ), f"error: {entry_true} should be recognized as coordinate."
        assert (
            self.parser.is_coordinate(entry_false1) == False
        ), f"error: {entry_false1} should not be recognized as coordinate."
        assert (
            self.parser.is_coordinate(entry_false2) == False
        ), f"error: {entry_false2} should not be recognized as coordinate."
        assert (
            self.parser.is_coordinate(entry_false3) == False
        ), f"error: {entry_false3} should not be recognized as coordinate."

    def test_get_fasta_regions(self) -> None:
        """Test if the parser extracts fasta regions correctly."""
        expected_result = ["16"]
        result = self.parser.get_fasta_regions(FILE_FASTA)
        assert (
            result == expected_result
        ), f"error: fasta regions not correctly extracted. Expected ['16'] got {result}"

    def test_read_fasta_sequences_existing_regions(self) -> None:
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

    def test_parse_fasta_header(self) -> None:
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

    def test_write_fasta_sequences(self) -> None:
        """Test if sequences are correctly written to fasta file."""
        file_out = os.path.join(self.tmp_path, "output.fna")

        fasta_sequences = self.parser.read_fasta_sequences(FILE_FASTA)
        self.parser.write_fasta_sequences(fasta_sequences=fasta_sequences, file_out=file_out)

        try:
            self.parser.check_fasta_format(file=file_out)
        except Exception as e:
            assert False, f"error: raised an exception: {e}, with written file."

    def test_merge_fasta_files(self) -> None:
        """Test if fasta files are merged correctly."""
        file_out = os.path.join(self.tmp_path, "output_merged.fna")

        self.parser.merge_fasta_files(
            files_in=[FILE_FASTA, FILE_NCBI_EXONS], file_out=file_out, overwrite=True
        )

        try:
            self.parser.check_fasta_format(file=file_out)
        except Exception as e:
            assert False, f"error: raised an exception: {e}, with merged files."

    def test_index_fasta_file_creates_index(self) -> None:
        """Test that .fai file is created when indexing a FASTA file."""
        file_fasta = os.path.join(self.tmp_path, "test_index.fna")
        index_file = f"{file_fasta}.fai"

        # Copy a test FASTA file
        shutil.copy2(FILE_FASTA, file_fasta)

        # Ensure index doesn't exist initially
        if os.path.exists(index_file):
            os.remove(index_file)

        # Index the FASTA file
        self.parser.index_fasta_file(file_fasta=file_fasta)

        # Verify that the index file was created
        assert os.path.exists(index_file), f"error: index file {index_file} should have been created"

    def test_index_fasta_file_removes_stale_index(self) -> None:
        """Test that existing stale index is removed before creating new one."""
        file_fasta = os.path.join(self.tmp_path, "test_index_stale.fna")
        index_file = f"{file_fasta}.fai"

        # Copy a test FASTA file
        shutil.copy2(FILE_FASTA, file_fasta)

        # Create a stale index file (with dummy content)
        with open(index_file, "w") as f:
            f.write("stale_index_content\n")

        # Verify stale index exists
        assert os.path.exists(index_file), "error: stale index file should exist before test"

        # Index the FASTA file (should remove stale index and create new one)
        self.parser.index_fasta_file(file_fasta=file_fasta)

        # Verify that a new index file was created (not the stale one)
        assert os.path.exists(index_file), f"error: index file {index_file} should exist after indexing"
        # Verify the index file is not the stale content
        with open(index_file, "r") as f:
            content = f.read()
            assert "stale_index_content" not in content, "error: stale index content should have been removed"


class TestVCFParser(unittest.TestCase):
    def setUp(self) -> None:
        self.tmp_path = os.path.join(os.getcwd(), "tmp_vcf_parser")
        Path(self.tmp_path).mkdir(parents=True, exist_ok=True)

        self.parser = VCFParser()

    def tearDown(self) -> None:
        shutil.rmtree(self.tmp_path)

    def test_check_vcf_format(self) -> None:
        """Test parsing fasta file."""
        try:
            out = self.parser.check_vcf_format(FILE_VCF)
            assert out == True, f"error: checker: check_fasta_format should have passed with file {FILE_VCF}"
        except Exception as e:
            assert False, f"error: checker: check_fasta_format raised an exception: {e}, with file {FILE_VCF}"

        try:
            out = self.parser.check_vcf_format(FILE_GFF)
            assert (
                out == False
            ), f"error: checker: check_fasta_format did not raise an exception with file {FILE_GFF}"
        except Exception:
            pass  # should go into this case

    def test_read_vcf_variants(self) -> None:
        variants, vcf_in = self.parser.read_vcf_variants(FILE_VCF)

        variant_type = variants[0].INFO.get("VC")
        variant_id = variants[0].ID

        assert variant_type == "SNV", f"error: wrong variant {variant_type} loaded."
        assert variant_id == "rs931559949", f"error: wrong variant {variant_id} loaded."

    def test_write_vcf_variants(self) -> None:
        file_out = os.path.join(self.tmp_path, "variants.vcf")

        variants, vcf_in = self.parser.read_vcf_variants(FILE_VCF)
        self.parser.write_vcf_variants(vcf_variants=variants, vcf_in=vcf_in, file_out=file_out)

        assert self.parser.check_vcf_format(file_out) == True, "error: vcf file stored in wrong format"

    def test_merge_vcf_files(self) -> None:
        file_out = os.path.join(self.tmp_path, "variants_merged.vcf")
        self.parser.merge_vcf_files(files_in=[FILE_VCF, FILE_VCF], file_out=file_out)

        assert self.parser.check_vcf_format(file_out) == True, "error: vcf file stored in wrong format"


class TestCountKmerAbundance(unittest.TestCase):
    def setUp(self) -> None:
        self.tmp_path = os.path.join(os.getcwd(), "tmp_count_kmer_abundance")
        Path(self.tmp_path).mkdir(parents=True, exist_ok=True)

        self.fasta_file_1 = os.path.join(self.tmp_path, "sequences1.fna")
        with open(self.fasta_file_1, "w") as f:
            f.write(">seq1\nAAGGGAAAAA\n>seq2\nACTGGATAATGCATGC\n>seq3\nTTTTTCCCTT\n")

        self.fasta_file_2 = os.path.join(self.tmp_path, "sequences2.fna")
        with open(self.fasta_file_2, "w") as f:
            f.write(">seq1\nAAAAAAAAAA\n>seq2\nATGCATGCATGC\n>seq3\nTTTTTTTTTT\n")

        self.file_empty = os.path.join(self.tmp_path, "empty.fna")
        with open(self.file_empty, "w") as f:
            f.write("")

        self.file_lowercase = os.path.join(self.tmp_path, "test_lowercase.fna")
        with open(self.file_lowercase, "w") as f:
            f.write(">test_seq\n")
            f.write("atcgatcg\n")

    def tearDown(self) -> None:
        shutil.rmtree(self.tmp_path)

    def test_count_kmer_abundance_single_k_value(self) -> None:
        """Test k-mer counting with single k value (int)."""
        result = count_kmer_abundance(files_fasta=self.fasta_file_1, k=4)

        assert 4 in result, "error: k=4 should be in result"
        assert isinstance(result[4], dict), "error: result[4] should be a dictionary"
        assert len(result) == 1, "error: should only have one k value in result"

    def test_count_kmer_abundance_k_range(self) -> None:
        """Test k-mer counting with k range (tuple)."""
        result = count_kmer_abundance(files_fasta=self.fasta_file_1, k=(3, 5))

        assert 3 in result, "error: k=3 should be in result"
        assert 4 in result, "error: k=4 should be in result"
        assert 5 in result, "error: k=5 should be in result"
        assert len(result) == 3, "error: should have 3 k values in result (3, 4, 5)"

    def test_count_kmer_abundance_k_list(self) -> None:
        """Test k-mer counting with k list."""
        result = count_kmer_abundance(files_fasta=self.fasta_file_1, k=[3, 5, 7])

        assert 3 in result, "error: k=3 should be in result"
        assert 5 in result, "error: k=5 should be in result"
        assert 7 in result, "error: k=7 should be in result"
        assert 4 not in result, "error: k=4 should not be in result"
        assert len(result) == 3, "error: should have 3 k values in result"

    def test_count_kmer_abundance_single_fasta_file(self) -> None:
        """Test k-mer counting with single FASTA file."""
        result = count_kmer_abundance(files_fasta=self.fasta_file_1, k=4)

        assert 4 in result, "error: k=4 should be in result"
        assert isinstance(result[4], dict), "error: result[4] should be a dictionary"
        assert len(result[4]) > 0, "error: should have k-mers in result"

    def test_count_kmer_abundance_multiple_fasta_files(self) -> None:
        """Test k-mer counting with multiple FASTA files."""
        result = count_kmer_abundance(files_fasta=[self.fasta_file_1, self.fasta_file_2], k=4)

        assert 4 in result, "error: k=4 should be in result"
        assert isinstance(result[4], dict), "error: result[4] should be a dictionary"
        assert len(result[4]) > 0, "error: should have k-mers in result"

    def test_count_kmer_abundance_fractional_abundance(self) -> None:
        """Test that fractional abundances sum to 1.0."""
        result = count_kmer_abundance(files_fasta=self.fasta_file_1, k=4)

        assert 4 in result, "error: k=4 should be in result"
        total_abundance = sum(result[4].values())
        assert (
            abs(total_abundance - 1.0) < 1e-10
        ), f"error: fractional abundances should sum to 1.0, got {total_abundance}"

    def test_count_kmer_abundance_sorting(self) -> None:
        """Test that k-mers are sorted by abundance in descending order."""
        result = count_kmer_abundance(files_fasta=self.fasta_file_1, k=4)

        assert 4 in result, "error: k=4 should be in result"
        abundances = list(result[4].values())
        # Check that abundances are in descending order
        for i in range(len(abundances) - 1):
            assert (
                abundances[i] >= abundances[i + 1]
            ), f"error: abundances should be sorted in descending order, got {abundances[i]} < {abundances[i + 1]}"

    def test_count_kmer_abundance_empty_fasta_file(self) -> None:
        """Test k-mer counting with empty FASTA file."""
        result = count_kmer_abundance(files_fasta=self.file_empty, k=4)

        assert 4 in result, "error: k=4 should be in result"
        assert len(result[4]) == 0, "error: empty file should result in empty k-mer dictionary"

    def test_count_kmer_abundance_invalid_k_negative(self) -> None:
        """Test that negative k values raise ValueError."""
        try:
            count_kmer_abundance(files_fasta=self.fasta_file_1, k=-1)
            assert False, "error: should have raised ValueError for negative k"
        except ValueError:
            pass  # Expected

    def test_count_kmer_abundance_invalid_k_tuple(self) -> None:
        """Test that invalid tuple raises ValueError."""
        try:
            count_kmer_abundance(files_fasta=self.fasta_file_1, k=(5, 3))  # k_min > k_max
            assert False, "error: should have raised ValueError for k_min > k_max"
        except ValueError:
            pass  # Expected

        try:
            count_kmer_abundance(files_fasta=self.fasta_file_1, k=(1, 2, 3))  # type: ignore
            assert False, "error: should have raised ValueError for tuple with wrong length"
        except ValueError:
            pass  # Expected

    def test_count_kmer_abundance_sequences_uppercase(self) -> None:
        """Test that sequences are converted to uppercase."""
        result = count_kmer_abundance(files_fasta=self.file_lowercase, k=4)

        # Check that all k-mers are uppercase
        for kmer in result[4].keys():
            assert kmer.isupper(), f"error: k-mer {kmer} should be uppercase"
        # Check that lowercase k-mer is not present, but uppercase version is
        assert "atcg" not in result[4], "error: lowercase k-mer should not be present"
        assert "ATCG" in result[4], "error: uppercase k-mer should be present"


class TestGetHighlyAbundantKmerSequences(unittest.TestCase):
    def setUp(self) -> None:
        self.tmp_path = os.path.join(os.getcwd(), "tmp_utils")
        os.makedirs(self.tmp_path, exist_ok=True)

        self.fasta_file = os.path.join(self.tmp_path, "kmer_example.fna")
        with open(self.fasta_file, "w") as f:
            f.write(">seq1\nAAAAAAAAAA\n>seq2\nATGCATGCATGC\n")

        self.thresholds = {3: 0.4, 4: 0.3, 5: 0.3}

    def tearDown(self) -> None:
        shutil.rmtree(self.tmp_path)

    def test_identifies_high_abundance_kmers(self) -> None:
        kmers = get_highly_abundant_kmer_sequences(
            files_fasta=self.fasta_file, kmer_abundance_threshold=self.thresholds
        )
        kmer_set = set(kmers)

        assert {"AAA", "AAAA", "AAAAA"}.issubset(
            kmer_set
        ), "error: high-abundance k-mers should be identified"

    def test_excludes_low_abundance_kmers(self) -> None:
        kmers = get_highly_abundant_kmer_sequences(
            files_fasta=str(self.fasta_file), kmer_abundance_threshold=self.thresholds
        )
        kmer_set = set(kmers)

        assert "ATG" not in kmer_set, "error: low-abundance 3-mer should not be included"
        assert "ATGC" not in kmer_set, "error: low-abundance 4-mer should not be included"
        assert "ATGCA" not in kmer_set, "error: low-abundance 5-mer should not be included"


class TestRemoveIndexFiles(unittest.TestCase):
    def setUp(self) -> None:
        self.tmp_path = os.path.join(os.getcwd(), "tmp_remove_index_files")
        Path(self.tmp_path).mkdir(parents=True, exist_ok=True)

    def tearDown(self) -> None:
        shutil.rmtree(self.tmp_path)

    def test_remove_fai_files_fasta(self) -> None:
        """Test removal of .fai files for FASTA."""
        file_reference = os.path.join(self.tmp_path, "reference.fna")
        index_file = os.path.join(self.tmp_path, "reference.fna.fai")

        # Create reference file and index file
        with open(file_reference, "w") as f:
            f.write(">seq1\nATCGATCG\n")
        with open(index_file, "w") as f:
            f.write("dummy index content\n")

        assert os.path.exists(index_file), "error: index file should exist before removal"
        assert os.path.exists(file_reference), "error: reference file should exist"

        remove_index_files(file_reference=file_reference, dir_output=self.tmp_path)

        assert not os.path.exists(index_file), "error: .fai index file should have been removed"
        assert os.path.exists(file_reference), "error: reference file should NOT be removed"

    def test_remove_csi_tbi_files_vcf(self) -> None:
        """Test removal of .csi and .tbi files for VCF."""
        file_reference_1 = os.path.join(self.tmp_path, "reference.vcf.gz")
        index_file_csi = os.path.join(self.tmp_path, "reference.vcf.gz.csi")

        # Create reference file and index files
        with open(file_reference_1, "w") as f:
            f.write("##fileformat=VCFv4.3\n")

        with open(index_file_csi, "w") as f:
            f.write("dummy csi content\n")

        assert os.path.exists(file_reference_1), "error: reference file should exist"
        assert os.path.exists(index_file_csi), "error: .csi index file should exist before removal"

        remove_index_files(file_reference=file_reference_1, dir_output=self.tmp_path)

        assert os.path.exists(file_reference_1), "error: reference file should NOT be removed"
        assert not os.path.exists(index_file_csi), "error: .csi index file should have been removed"

        file_reference_2 = os.path.join(self.tmp_path, "reference.vcf")
        index_file_tbi = os.path.join(self.tmp_path, "reference.vcf.tbi")

        with open(file_reference_2, "w") as f:
            f.write("##fileformat=VCFv4.3\n")
        with open(index_file_tbi, "w") as f:
            f.write("dummy tbi content\n")

        assert os.path.exists(file_reference_2), "error: reference file should exist"
        assert os.path.exists(index_file_tbi), "error: .tbi index file should exist before removal"

        remove_index_files(file_reference=file_reference_2, dir_output=self.tmp_path)

        assert os.path.exists(file_reference_2), "error: reference file should NOT be removed"
        assert not os.path.exists(index_file_tbi), "error: .tbi index file should have been removed"

    def test_remove_nhr_nin_nsq_files_blast(self) -> None:
        """Test removal of .nhr, .nin, .nsq files for BLAST."""
        file_reference = os.path.join(self.tmp_path, "reference.fna")
        index_file_nhr = os.path.join(self.tmp_path, "reference.fna.nhr")
        index_file_nin = os.path.join(self.tmp_path, "reference.fna.nin")
        index_file_nsq = os.path.join(self.tmp_path, "reference.fna.nsq")

        # Create reference file and index files
        with open(file_reference, "w") as f:
            f.write(">seq1\nATCGATCG\n")
        with open(index_file_nhr, "w") as f:
            f.write("dummy nhr content\n")
        with open(index_file_nin, "w") as f:
            f.write("dummy nin content\n")
        with open(index_file_nsq, "w") as f:
            f.write("dummy nsq content\n")

        assert os.path.exists(index_file_nhr), "error: .nhr index file should exist before removal"
        assert os.path.exists(index_file_nin), "error: .nin index file should exist before removal"
        assert os.path.exists(index_file_nsq), "error: .nsq index file should exist before removal"
        assert os.path.exists(file_reference), "error: reference file should exist"

        remove_index_files(file_reference=file_reference, dir_output=self.tmp_path)

        assert not os.path.exists(index_file_nhr), "error: .nhr index file should have been removed"
        assert not os.path.exists(index_file_nin), "error: .nin index file should have been removed"
        assert not os.path.exists(index_file_nsq), "error: .nsq index file should have been removed"
        assert os.path.exists(file_reference), "error: reference file should NOT be removed"

    def test_remove_ebwt_files_bowtie(self) -> None:
        """Test removal of .ebwt files for Bowtie."""
        file_reference = os.path.join(self.tmp_path, "reference.fna")
        index_file_1 = os.path.join(self.tmp_path, "reference.fna.1.ebwt")
        index_file_2 = os.path.join(self.tmp_path, "reference.fna.2.ebwt")
        index_file_rev = os.path.join(self.tmp_path, "reference.fna.rev.1.ebwt")

        # Create reference file and index files
        with open(file_reference, "w") as f:
            f.write(">seq1\nATCGATCG\n")
        with open(index_file_1, "w") as f:
            f.write("dummy ebwt content 1\n")
        with open(index_file_2, "w") as f:
            f.write("dummy ebwt content 2\n")
        with open(index_file_rev, "w") as f:
            f.write("dummy ebwt rev content\n")

        assert os.path.exists(index_file_1), "error: .1.ebwt index file should exist before removal"
        assert os.path.exists(index_file_2), "error: .2.ebwt index file should exist before removal"
        assert os.path.exists(index_file_rev), "error: .rev.1.ebwt index file should exist before removal"
        assert os.path.exists(file_reference), "error: reference file should exist"

        remove_index_files(file_reference=file_reference, dir_output=self.tmp_path)

        assert not os.path.exists(index_file_1), "error: .1.ebwt index file should have been removed"
        assert not os.path.exists(index_file_2), "error: .2.ebwt index file should have been removed"
        assert not os.path.exists(index_file_rev), "error: .rev.1.ebwt index file should have been removed"
        assert os.path.exists(file_reference), "error: reference file should NOT be removed"

    def test_remove_bt2_files_bowtie2(self) -> None:
        """Test removal of .bt2 files for Bowtie2."""
        file_reference = os.path.join(self.tmp_path, "reference.fna")
        index_file_1 = os.path.join(self.tmp_path, "reference.fna.1.bt2")
        index_file_2 = os.path.join(self.tmp_path, "reference.fna.2.bt2")
        index_file_rev = os.path.join(self.tmp_path, "reference.fna.rev.1.bt2")

        # Create reference file and index files
        with open(file_reference, "w") as f:
            f.write(">seq1\nATCGATCG\n")
        with open(index_file_1, "w") as f:
            f.write("dummy bt2 content 1\n")
        with open(index_file_2, "w") as f:
            f.write("dummy bt2 content 2\n")
        with open(index_file_rev, "w") as f:
            f.write("dummy bt2 rev content\n")

        assert os.path.exists(index_file_1), "error: .1.bt2 index file should exist before removal"
        assert os.path.exists(index_file_2), "error: .2.bt2 index file should exist before removal"
        assert os.path.exists(index_file_rev), "error: .rev.1.bt2 index file should exist before removal"
        assert os.path.exists(file_reference), "error: reference file should exist"

        remove_index_files(file_reference=file_reference, dir_output=self.tmp_path)

        assert not os.path.exists(index_file_1), "error: .1.bt2 index file should have been removed"
        assert not os.path.exists(index_file_2), "error: .2.bt2 index file should have been removed"
        assert not os.path.exists(index_file_rev), "error: .rev.1.bt2 index file should have been removed"
        assert os.path.exists(file_reference), "error: reference file should NOT be removed"

    def test_reference_file_not_removed(self) -> None:
        """Test that reference file itself is NOT removed."""
        file_reference = os.path.join(self.tmp_path, "reference.fna")
        index_file = os.path.join(self.tmp_path, "reference.fna.fai")
        other_file = os.path.join(self.tmp_path, "other_file.fna")

        # Create reference file, index file, and another unrelated file
        with open(file_reference, "w") as f:
            f.write(">seq1\nATCGATCG\n")
        with open(index_file, "w") as f:
            f.write("dummy index content\n")
        with open(other_file, "w") as f:
            f.write(">seq2\nGCTAGCTA\n")

        assert os.path.exists(file_reference), "error: reference file should exist"
        assert os.path.exists(index_file), "error: index file should exist"
        assert os.path.exists(other_file), "error: other file should exist"

        remove_index_files(file_reference=file_reference, dir_output=self.tmp_path)

        assert os.path.exists(file_reference), "error: reference file should NOT be removed"
        assert not os.path.exists(index_file), "error: index file should have been removed"
        assert os.path.exists(other_file), "error: unrelated file should NOT be removed"
