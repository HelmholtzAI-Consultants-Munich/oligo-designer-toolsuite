############################################
# imports
############################################

import os
import shutil
import unittest
import warnings

from oligo_designer_toolsuite.utils import (
    FastaParser,
    GffParser,
    check_if_dna_sequence,
    check_if_key_exists,
    check_if_list,
    check_tsv_format,
    get_sequence_from_annotation,
)

############################################
# Global Parameters
############################################


############################################
# Tests
############################################


class TestGffParser(unittest.TestCase):
    def setUp(self):
        self.parser = GffParser()
        self.gff_file = "data/annotations/custom_GCF_000001405.40_GRCh38.p14_genomic_chr16.gff"
        self.gtf_file = "data/annotations/custom_GCF_000001405.40_GRCh38.p14_genomic_chr16.gtf"
        self.pickle_file = ...

    def test_check_gff_format(self):
        """Test parsing GFF annotation data."""
        try:
            self.parser.check_gff_format(self.gff_file)
        except Exception as e:
            assert (
                False
            ), f"error: checker: check_gff_format raised an exception: {e}, with file {self.gff_file}"

    def test_check_gtf_format(self):
        """Test parsing GTF annotation data."""
        try:
            self.parser.check_gff_format(self.gtf_file)
        except Exception as e:
            assert (
                False
            ), f"error: checker: check_gff_format raised an exception: {e}, with file {self.gtf_file}"

    # TODO
    def test_load_annotation_from_pickle_file(self):
        """Test loading annotation from a pickle file."""
        expected_result = ...
        result = self.parser.load_annotation_from_pickle(self.pickle_file)
        assert result == expected_result, f"error: loading annotation from pickle file failed"

    def test_parse_annotation_from_gff(self):
        """Test parsing GFF annotation."""
        result = self.parser.parse_annotation_from_gff(self.gff_file, target_lines=10)
        assert result.shape[1] == 23, "error: GFF3 dataframe not correctly loaded"

    def test_parse_annotation_from_gtf(self):
        """Test parsing GTF annotation."""
        result = self.parser.parse_annotation_from_gff(self.gtf_file, target_lines=10)
        assert result.shape[1] == 20, "error: GTF dataframe not correctly loaded"


class TestFastaParser(unittest.TestCase):
    def setUp(self):
        self.parser = FastaParser()
        self.fasta_file = "data/annotations/custom_GCF_000001405.40_GRCh38.p14_genomic_chr16.fna"

    def test_check_fasta_format(self):
        """Test parsing fasta file."""
        try:
            self.parser.check_fasta_format(self.fasta_file)
        except Exception as e:
            assert (
                False
            ), f"error: checker: check_fasta_format raised an exception: {e}, with file {self.fasta_file}"

    # TODO
    def test_read_fasta_sequences_existing_regions(self):
        """Test parsing fasta file."""
        expected_result = ...
        ids = ...

        result = self.parser.read_fasta_sequences(self.fasta_file, region_ids=ids)
        assert result == expected_result, "error: fasta dataframe not correctly loaded"

    # TODO
    def test_read_fasta_sequences_non_existing_region(self):
        """Test parsing fasta file."""
        expected_result = ...
        ids = ...

        result = self.parser.read_fasta_sequences(self.fasta_file, region_ids=ids)
        # check if a warning was raised
        with warnings.catch_warnings(record=True) as w:
            result = self.parser.read_fasta_sequences(self.fasta_file, region_ids=ids)
            assert len(w) > 0, "error: no warning was raised"
        assert result == expected_result, "error: fasta dataframe not correctly loaded"

    def test_get_fasta_regions(self):
        """Test if the parser extracts fasta regions correctly."""
        expected_result = ...
        result = self.parser.get_fasta_regions(self.fasta_file)
        assert result == expected_result, "error: fasta regions not correctly extracted"

    def test_parse_fasta_header(self):
        """Test if the parser extracts fasta header correctly."""
        header = "ARPG3::transcript_id=XM4581;exon_id=XM4581_exon1::16:70265537-70265662(-)"
        region, additional_information, coordinates = self.parser.parse_fasta_header(header)
        assert region == "ARPG3", "error: wrong region parsed"
        assert coordinates["chromosome"] == ["16"], "error: wrong chrom parsed"
        assert coordinates["start"] == [70265537], "error: wrong start parsed"
        assert coordinates["end"] == [70265662], "error: wrong end parsed"
        assert coordinates["strand"] == ["-"], "error: wrong strand parsed"
        assert (
            additional_information == "transcript_id=XM4581;exon_id=XM4581_exon1"
        ), "error: wrong additional information parsed"


class TestCheckers(unittest.TestCase):
    def test_check_tsv_format(self):
        """Test if the parser extracts fasta header correctly."""
        file_tsv = "data/annotations/custom_GCF_000001405.40_GRCh38.p14_genomic_chr16.gtf.tsv"
        try:
            check_tsv_format(file_tsv)
        except Exception as e:
            assert False, f"error: checker: check_tsv_format raised an exception: {e}, with file {file_tsv}"

    # TODO
    def test_check_if_key_exists(self):
        """Test if the parser extracts fasta header correctly."""
        nested_dict = ...
        key = ...
        expected_result = ...
        result = check_if_key_exists(nested_dict, key)
        assert result == expected_result, "error: check_if_key_exists failed"

    def test_check_if_list_str(self):
        """Test if the parser extracts fasta header correctly."""
        s = "test"
        result = check_if_list(s)
        assert result == [s], f"error: check_if_list failed. Expected: [{s}], got: {result}"

    def test_check_if_list_list(self):
        """Test if the parser extracts fasta header correctly."""
        l = ["test", ["test2"]]
        result = check_if_list(l)
        assert result == l, f"error: check_if_list failed. Expected: {l}, got: {result}"

    def test_check_if_dna_sequence_valid(self):
        """Test if the parser extracts fasta header correctly."""
        seq = "GGctAAgTTCCaGTttGCA"
        valid_characters = ["A", "C", "T", "G"]
        assert check_if_dna_sequence(seq, valid_characters), "error: check_if_dna_sequence failed"

    def test_check_if_dna_sequence_invalid(self):
        """Test if the parser extracts fasta header correctly."""
        seq = "GGctAAgTuuTCCaGTttGCA"
        valid_characters = ["A", "C", "T", "G", "W", "X"]
        assert not check_if_dna_sequence(
            seq, valid_characters
        ), "error: check_if_dna_sequence succeeded when it should have failed"


class TestSequenceProcessor(unittest.TestCase):
    def setUp(self):
        self.tmp_path = "tests/tmp"
        os.makedirs(self.tmp_path, exist_ok=True)
        self.file_bed = "data/annotations/custom_GCF_000001405.40_GRCh38.p14_genomic_chr16.bed"

    def tearDown(self):
        shutil.rmtree(self.tmp_path)

    def test_get_sequence_from_annotation(self):
        """Test if the parser extracts fasta header correctly."""
        parser = FastaParser()

        file_reference_fasta = "data/annotations/custom_GCF_000001405.40_GRCh38.p14_genomic_chr16.fna"
        file_fasta = os.path.join(self.tmp_path, "test_ann2seq_function.fna")
        get_sequence_from_annotation(
            self.file_bed,
            file_reference_fasta,
            file_fasta,
            split=False,
            strand=True,
            name=True,
        )
        res = parser.check_fasta_format(file_fasta)
        assert res == True, f"error: the created sequence file is not a fasta file"

    # TODO
    def test_get_complement_regions(self):
        ...

    # """Test if the parser extracts fasta header correctly."""
    # file_bed_in = ...
    # file_chromosome_length = ...
    # file_bed_out = ...
    # get_complement_regions(file_bed_in, file_chromosome_length, file_bed_out)
    # res = ...
    # assert res == ..., f"error: get_complement_regions failed"
