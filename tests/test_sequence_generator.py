############################################
# imports
############################################

import os
import shutil
import unittest
from abc import ABC, abstractmethod
from typing import Any, Callable, Generic, cast
from unittest.mock import patch

from Bio import SeqIO

from oligo_designer_toolsuite._exceptions import ConfigurationError
from oligo_designer_toolsuite.database import OligoDatabase
from oligo_designer_toolsuite.genomic_database_downloader import (
    GenomicDatabaseDownloaderEnsembl,
    GenomicDatabaseDownloaderNcbi,
    GenomicEntity,
    NcbiGenomicEntity,
)
from oligo_designer_toolsuite.genomic_database_downloader._genomic_database_downloader_base import (
    GenomicDatabaseDownloaderBase,
    TGenomicEntity,
)
from oligo_designer_toolsuite.genomic_database_downloader._genomic_database_downloader_ncbi import (
    AssemblyMode,
    SpeciesMode,
)
from oligo_designer_toolsuite.oligo_property_calculator import (
    LengthProperty,
    PropertyCalculator,
)
from oligo_designer_toolsuite.sequence_generator import (
    CustomGenomicRegionGenerator,
    OligoSequenceGenerator,
)
from oligo_designer_toolsuite.utils import FastaParser, check_if_dna_sequence, get_complement_regions

from .expected_values_region_generator import (
    EXPECTED_HEADER_VALUES_BACTERIA_NCBI,
    EXPECTED_HEADER_VALUES_HUMAN_ENSEMBL,
    EXPECTED_HEADER_VALUES_HUMAN_NCBI,
    EXPECTED_HEADER_VALUES_MOUSE_NCBI,
    RegionHeaderSpec,
)

############################################
# Setup
############################################

FILE_ANNOTATION_ENSEMBL = "tests/data/annotations/custom_Homo_sapiens.GRCh38.108.chr16.gtf"
FILE_SEQUENCE_ENSEMBL = "tests/data/annotations/custom_Homo_sapiens.GRCh38.dna_sm.chromosome.16.fa"

FILE_ANNOTATION_NCBI = "tests/data/annotations/custom_GCF_000001405.40_GRCh38.p14_genomic_chr16.gtf"
FILE_SEQUENCE_NCBI = "tests/data/annotations/custom_GCF_000001405.40_GRCh38.p14_genomic_chr16.fna"

METADATA_NCBI = {
    "files_source": "NCBI",
    "species": "Homo_sapiens",
    "annotation_release": "110",
    "genome_assembly": "GRCh38",
}

METADATA_ENSEMBL = {
    "files_source": "Ensembl",
    "species": "Homo_sapiens",
    "annotation_release": "108",
    "genome_assembly": "GRCh38",
}

FILE_ANNOTATION_MOUSE_NCBI = (
    "tests/data/annotations/custom_GCF_000001635.26_GRCm38.p6_genomic_NC_000085.6.gtf"
)
FILE_SEQUENCE_MOUSE_NCBI = "tests/data/annotations/custom_GCF_000001635.26_GRCm38.p6_genomic_NC_000085.6.fna"

METADATA_MOUSE_NCBI = {
    "files_source": "NCBI",
    "species": "Mus_musculus",
    "annotation_release": "108.20200622",
    "genome_assembly": "GRCm38.p6",
}

FILE_ANNOTATION_BACTERIA_NCBI = "tests/data/annotations/GCF_000068585.1_ASM6858v1_genomic.gtf"
FILE_SEQUENCE_BACTERIA_NCBI = "tests/data/annotations/GCF_000068585.1_ASM6858v1_genomic.fna"

METADATA_BACTERIA_NCBI = {
    "files_source": "NCBI",
    "species": "Chlamydia_trachomatis",
    "annotation_release": "no_annotation_info",
    "genome_assembly": "ASM6858v1",
}

FILE_NCBI_EXONS = "tests/data/genomic_regions/sequences_ncbi_exons.fna"
FILE_NCBI_EXON_EXON_JUNCTIONS = "tests/data/genomic_regions/sequences_ABAT_ncbi_exon_exon_junctions.fna"
FILE_NCBI_EXON_EXON_JUNCTIONS_SHORT = (
    "tests/data/genomic_regions/sequences_AARS1_ncbi_exon_exon_junctions_short.fna"
)


############################################
# Tests
############################################
class GenomicDatabaseDownloaderTestBase(Generic[TGenomicEntity]):
    def setUp(self) -> None:
        self.tmp_path = os.path.join(os.getcwd(), "tmp_genomic_database_downloader")
        os.makedirs(self.tmp_path, exist_ok=True)

        self.loader = self.setup_genomic_database_downloader()
        self.entity = self.setup_genomic_entity()

    def tearDown(self) -> None:
        shutil.rmtree(self.tmp_path, ignore_errors=True)

    @abstractmethod
    def setup_genomic_database_downloader(self) -> GenomicDatabaseDownloaderBase[TGenomicEntity]:
        pass

    @abstractmethod
    def setup_genomic_entity(self) -> TGenomicEntity:
        pass


class GenomicDatabaseDownloaderDownloadBase(GenomicDatabaseDownloaderTestBase):
    def setUp(self) -> None:
        self.tmp_path = os.path.join(os.getcwd(), "tmp_genomic_database_downloader")
        os.makedirs(self.tmp_path, exist_ok=True)
        self.loader = self.setup_genomic_database_downloader()
        self.genomic_entity = self.setup_genomic_entity()

    def tearDown(self) -> None:
        shutil.rmtree(self.tmp_path)

    def test_download(self) -> None:
        _ = self.loader.fetch_genomic_entity(self.genomic_entity, ["fasta"])


class TestGenomicDatabaseDownloaderNCBICurrent(GenomicDatabaseDownloaderDownloadBase, unittest.TestCase):
    def setup_genomic_database_downloader(self) -> GenomicDatabaseDownloaderNcbi:
        return GenomicDatabaseDownloaderNcbi(self.tmp_path)

    def setup_genomic_entity(self) -> NcbiGenomicEntity:
        return NcbiGenomicEntity(
            # Parameters
            mode="species",
            taxon="vertebrate_mammalian",  # taxon the species belongs to
            species="Homo_sapiens",
            annotation_release="current",
            assembly_name=None,
            assembly_accession=None,
            assembly_source=None,
        )


class TestGenomicDatabaseDownloaderNCBICurrentProkaryotes(
    GenomicDatabaseDownloaderDownloadBase, unittest.TestCase
):
    def setup_genomic_database_downloader(self) -> GenomicDatabaseDownloaderNcbi:
        return GenomicDatabaseDownloaderNcbi(self.tmp_path)

    def setup_genomic_entity(self) -> NcbiGenomicEntity:
        return NcbiGenomicEntity(
            # Parameters
            mode="species",
            taxon="bacteria",  # taxon the species belongs to
            species="Actinomadura_yumaensis",
            annotation_release="current",
            assembly_name=None,
            assembly_accession=None,
            assembly_source=None,
        )


class TestGenomicDatabaseDownloaderEnsemblCurrent(GenomicDatabaseDownloaderDownloadBase, unittest.TestCase):
    def setup_genomic_database_downloader(self) -> GenomicDatabaseDownloaderEnsembl:
        return GenomicDatabaseDownloaderEnsembl(self.tmp_path)

    def setup_genomic_entity(self) -> GenomicEntity:
        return GenomicEntity(taxon=None, species="homo_sapiens", release="current")


class TestGenomicDatabaseDownloaderNCBIModeValidation(unittest.TestCase):
    def setUp(self) -> None:
        self.tmp_path = os.path.join(os.getcwd(), "tmp_genomic_database_downloader_mode_validation")
        os.makedirs(self.tmp_path, exist_ok=True)
        self.loader = GenomicDatabaseDownloaderNcbi(self.tmp_path)

    def tearDown(self) -> None:
        shutil.rmtree(self.tmp_path, ignore_errors=True)

    def test_assembly_mode_requires_both_fields(self) -> None:
        entity = NcbiGenomicEntity(
            mode="assembly",
            taxon=None,
            species=None,
            annotation_release=None,
            assembly_name=None,
            assembly_accession="GCF_000001405.38",
            assembly_source=None,
        )

        with self.assertRaises(ConfigurationError):
            self.loader._get_mode(entity)

    def test_assembly_mode_rejects_species_fields(self) -> None:
        entity = NcbiGenomicEntity(
            mode="assembly",
            taxon="vertebrate_mammalian",
            species="Homo_sapiens",
            annotation_release="current",
            assembly_name="GRCh38.p12",
            assembly_accession="GCF_000001405.38",
            assembly_source=None,
        )

        with self.assertRaises(ConfigurationError):
            self.loader._get_mode(entity)

    def test_assembly_mode_accepts_valid_configuration(self) -> None:
        entity = NcbiGenomicEntity(
            mode="assembly",
            taxon=None,
            species=None,
            annotation_release=None,
            assembly_name="GRCh38.p12",
            assembly_accession="GCF_000001405.38",
            assembly_source=None,
        )

        mode = self.loader._get_mode(entity)

        self.assertIsInstance(mode, AssemblyMode)

    def test_requires_explicit_mode(self) -> None:
        entity = NcbiGenomicEntity(
            mode=None,
            taxon="vertebrate_mammalian",
            species="Homo_sapiens",
            annotation_release="current",
            assembly_name=None,
            assembly_accession=None,
            assembly_source=None,
        )

        with self.assertRaises(ConfigurationError):
            self.loader._get_mode(entity)

    def test_rejects_unknown_mode(self) -> None:
        entity = NcbiGenomicEntity(
            mode="unknown",
            taxon="vertebrate_mammalian",
            species="Homo_sapiens",
            annotation_release="current",
            assembly_name=None,
            assembly_accession=None,
            assembly_source=None,
        )

        with self.assertRaises(ConfigurationError):
            self.loader._get_mode(entity)


class TestGenomicDatabaseDownloaderNCBIUnitBehavior(unittest.TestCase):
    def setUp(self) -> None:
        self.tmp_path = os.path.join(os.getcwd(), "tmp_genomic_database_downloader_unit_behavior")
        os.makedirs(self.tmp_path, exist_ok=True)
        self.loader = GenomicDatabaseDownloaderNcbi(self.tmp_path)

    def tearDown(self) -> None:
        shutil.rmtree(self.tmp_path)

    def getAssemblyMode(self, entity: NcbiGenomicEntity) -> AssemblyMode:
        mode = self.loader._get_mode(entity)
        assert isinstance(mode, AssemblyMode)
        return mode

    def getSpeciesMode(self, entity: NcbiGenomicEntity) -> SpeciesMode:
        mode = self.loader._get_mode(entity)
        assert isinstance(mode, SpeciesMode)
        return mode

    def test_direct_mode_builds_expected_ncbi_all_path(self) -> None:
        entity = NcbiGenomicEntity(
            mode="assembly",
            taxon=None,
            species=None,
            annotation_release=None,
            assembly_name="GRCh38.p12",
            assembly_accession="GCF_000001405.38",
            assembly_source=None,
        )

        mode = self.getAssemblyMode(entity)

        assert (
            mode._resolve_directory_from_direct_assembly()
            == "genomes/all/GCF/000/001/405/GCF_000001405.38_GRCh38.p12/"
        )

    def test_direct_mode_rejects_invalid_accession_format(self) -> None:
        entity = NcbiGenomicEntity(
            mode="assembly",
            taxon=None,
            species=None,
            annotation_release=None,
            assembly_name="GRCh38.p12",
            assembly_accession="GCF_123",
            assembly_source=None,
        )

        mode = self.getAssemblyMode(entity)

        with self.assertRaises(ConfigurationError):
            mode._resolve_directory_from_direct_assembly()

    def test_taxon_mode_requires_required_fields(self) -> None:
        entity = NcbiGenomicEntity(
            mode="species",
            taxon=None,
            species="vertebrate_mammalian",
            annotation_release=None,
            assembly_name=None,
            assembly_accession=None,
            assembly_source=None,
        )

        with self.assertRaises(ConfigurationError):
            self.loader.fetch_genomic_entity(entity)

    def test_unsupported_taxon_is_rejected(self) -> None:
        entity = NcbiGenomicEntity(
            mode="species",
            taxon="mitochondrion",
            species="Homo_sapiens",
            annotation_release="current",
            assembly_name=None,
            assembly_accession=None,
            assembly_source=None,
        )

        with self.assertRaises(ConfigurationError):
            self.loader.fetch_genomic_entity(entity)

    def test_unknown_taxon_is_rejected(self) -> None:
        entity = NcbiGenomicEntity(
            mode="species",
            taxon="foo_taxon",
            species="Homo_sapiens",
            annotation_release="current",
            assembly_name=None,
            assembly_accession=None,
            assembly_source=None,
        )

        with self.assertRaises(ConfigurationError):
            self.loader.fetch_genomic_entity(entity)

    def test_assembly_source_incompatible_with_taxon_is_rejected(self) -> None:
        entity = NcbiGenomicEntity(
            mode="species",
            taxon="viral",
            species="SARS-CoV-2",
            annotation_release="current",
            assembly_name=None,
            assembly_accession=None,
            assembly_source="annotation_releases",
        )

        with self.assertRaises(ConfigurationError):
            self.loader.fetch_genomic_entity(entity)

    def test_numeric_release_is_rejected_for_latest_assembly_versions(self) -> None:
        entity = NcbiGenomicEntity(
            mode="species",
            taxon="bacteria",
            species="Actinomadura_yumaensis",
            annotation_release="110",
            assembly_source="latest_assembly_versions",
            assembly_name=None,
            assembly_accession=None,
        )

        with self.assertRaises(ConfigurationError):
            self.loader.fetch_genomic_entity(entity)

    def test_current_annotation_release_uses_first_entry(self) -> None:
        entity = NcbiGenomicEntity(
            mode="species",
            taxon="vertebrate_mammalian",
            species="Homo_sapiens",
            annotation_release="current",
            assembly_source="annotation_releases",
            assembly_name=None,
            assembly_accession=None,
        )

        species_mode = self.getSpeciesMode(entity)

        with patch.object(
            species_mode,
            "_list_ftp_entries",
            return_value=["GCF_000001405.40-RS_2025_08", "GCF_009914755.1-RS_2025_08"],
        ):
            ftp_directory = species_mode._resolve_base_directory("annotation_releases")
        assert species_mode.annotation_release == "GCF_000001405.40-RS_2025_08"
        assert (
            ftp_directory == "genomes/refseq/vertebrate_mammalian/Homo_sapiens/annotation_releases/"
            "GCF_000001405.40-RS_2025_08/"
        )

    def test_non_annotation_source_uses_first_gcf_entry(self) -> None:
        entity = NcbiGenomicEntity(
            mode="species",
            taxon="bacteria",
            species="Actinomadura_yumaensis",
            annotation_release="current",
            assembly_source="latest_assembly_versions",
            assembly_name=None,
            assembly_accession=None,
        )

        species_mode = self.getSpeciesMode(entity)

        with patch.object(
            species_mode,
            "_list_ftp_entries",
            return_value=["README.txt", "GCF_003054545.1_ASM305454v1", "GCA_111111111.1_OTHER"],
        ):
            ftp_directory = species_mode._resolve_base_directory("latest_assembly_versions")
        assert (
            ftp_directory == "genomes/refseq/bacteria/Actinomadura_yumaensis/latest_assembly_versions/"
            "GCF_003054545.1_ASM305454v1/"
        )

    def test_non_annotation_source_without_gcf_entry_raises(self) -> None:
        entity = NcbiGenomicEntity(
            mode="species",
            taxon="bacteria",
            species="Actinomadura_yumaensis",
            annotation_release="current",
            assembly_source="latest_assembly_versions",
            assembly_name=None,
            assembly_accession=None,
        )

        species_mode = self.getSpeciesMode(entity)

        with patch.object(
            species_mode, "_list_ftp_entries", return_value=["README.txt", "GCA_123456789.1_OTHER"]
        ):
            with self.assertRaises(ConfigurationError):
                species_mode._resolve_base_directory("latest_assembly_versions")

    def test_resolve_assembly_metadata_falls_back_to_assembly_report(self) -> None:
        entity = NcbiGenomicEntity(
            mode="species",
            taxon="bacteria",
            species="Xanthobacter_sp._SG618",
            annotation_release="current",
            assembly_source="reference",
            assembly_name=None,
            assembly_accession=None,
        )

        species_mode = self.getSpeciesMode(entity)

        assembly_report_file_name = "GCF_012932745.1_ASM1293274v1_assembly_report.txt"
        assembly_report = os.path.join(self.tmp_path, assembly_report_file_name)
        assembly_report_copy = os.path.join(self.tmp_path, f"{assembly_report_file_name}-copy")

        with open(assembly_report, "w") as handle:
            handle.write("# Assembly name: ASM1293274v1\n")
            handle.write("# RefSeq assembly accession: GCF_012932745.1\n")

        def _download_side_effect(ftp_directory: str, file_name: str) -> str:
            if "README_" in file_name:
                raise FileNotFoundError("README missing")
            if "_assembly_report" in file_name:
                shutil.copy(assembly_report, assembly_report_copy)
                return assembly_report_copy
            raise FileNotFoundError("No file")

        def _get_exact_filename_side_effect(ftp_directory: str, pattern: str) -> str:
            return assembly_report_file_name

        with (
            patch.object(species_mode.ftp_util, "download", side_effect=_download_side_effect),
            patch.object(
                species_mode.ftp_util, "get_exact_filename", side_effect=_get_exact_filename_side_effect
            ),
        ):
            assembly_accession, assembly_name = species_mode._resolve_assembly_metadata("dummy_dir/")

        assert assembly_name == "ASM1293274v1"
        assert assembly_accession == "GCF_012932745.1"


class GenomicDatabaseDownloaderFilesBase(GenomicDatabaseDownloaderTestBase):
    @abstractmethod
    def get_correct_metadata(self) -> tuple[str, str]:
        pass

    @abstractmethod
    def get_correct_gtf(self) -> str:
        pass

    @abstractmethod
    def get_correct_fasta(self) -> str:
        pass

    def assert_file_name_is_included(self, file_path: str | None, expected_file_name: str) -> None:
        assert file_path is not None
        assert file_path.find(expected_file_name) > -1, "error: wrong file downloaded"

    def test_metadata_loader(self) -> None:
        result = self.loader.fetch_genomic_entity(
            self.entity,
            files=["fasta"],
        )

        annotation_release, assembly_name = self.get_correct_metadata()

        assert result["annotation_release"] == annotation_release, "error: wrong annotation release retrieved"

        assert result["genome_assembly"] == assembly_name, "error: wrong assembly retrieved"

    def test_gtf_loader(self) -> None:
        result = self.loader.fetch_genomic_entity(
            self.entity,
            files=["gtf"],
        )

        self.assert_file_name_is_included(result["fetched_files"]["gtf"], self.get_correct_gtf())

    def test_fasta_loader(self) -> None:
        result = self.loader.fetch_genomic_entity(
            self.entity,
            files=["fasta"],
        )

        self.assert_file_name_is_included(result["fetched_files"]["fasta"], self.get_correct_fasta())


class TestGenomicDatabaseDownloaderNCBIOldAnnotations(
    GenomicDatabaseDownloaderFilesBase,
    unittest.TestCase,
):
    def setup_genomic_database_downloader(self) -> GenomicDatabaseDownloaderNcbi:
        return GenomicDatabaseDownloaderNcbi(self.tmp_path)

    def setup_genomic_entity(self) -> NcbiGenomicEntity:
        return NcbiGenomicEntity(
            mode="species",
            taxon="vertebrate_mammalian",
            species="Homo_sapiens",
            annotation_release="110",
            assembly_source=None,
            assembly_accession=None,
            assembly_name=None,
        )

    def get_correct_metadata(self) -> tuple[str, str]:
        return (
            "NCBI_Homo_sapiens_Annotation_Release_110",
            "GRCh38.p14",
        )

    def get_correct_gtf(self) -> str:
        return "GCF_000001405.40_GRCh38.p14_genomic.gtf"

    def get_correct_fasta(self) -> str:
        return "GCF_000001405.40_GRCh38.p14_genomic.fna"


class TestGenomicDatabaseDownloaderNCBIAssemblyNumber(
    GenomicDatabaseDownloaderFilesBase,
    unittest.TestCase,
):
    def setup_genomic_database_downloader(self) -> GenomicDatabaseDownloaderNcbi:
        return GenomicDatabaseDownloaderNcbi(self.tmp_path)

    def setup_genomic_entity(self) -> NcbiGenomicEntity:
        return NcbiGenomicEntity(
            mode="assembly",
            taxon=None,
            species=None,
            annotation_release=None,
            assembly_source=None,
            assembly_accession="GCF_000068585.1",
            assembly_name="ASM6858v1",
        )

    def get_correct_metadata(self) -> tuple[str, str]:
        annotation_release = "no_annotation_info"
        assembly_name = "ASM6858v1"

        return annotation_release, assembly_name

    def get_correct_gff(self) -> str:
        return "GCF_000068585.1_ASM6858v1_genomic.gff"

    def get_correct_gtf(self) -> str:
        return "GCF_000068585.1_ASM6858v1_genomic.gtf"

    def get_correct_fasta(self) -> str:
        return "GCF_000068585.1_ASM6858v1_genomic.fna"


class TestGenomicDatabaseDownloaderEnsemblOldAnnotations(
    GenomicDatabaseDownloaderFilesBase, unittest.TestCase
):
    def setup_genomic_database_downloader(self) -> GenomicDatabaseDownloaderEnsembl:
        return GenomicDatabaseDownloaderEnsembl(self.tmp_path)

    def setup_genomic_entity(self) -> GenomicEntity:
        return GenomicEntity(taxon=None, species="homo_sapiens", release="108")

    def get_correct_metadata(self) -> tuple[str, str]:
        annotation_release = "108"
        assembly_name = "GRCh38"

        return annotation_release, assembly_name

    def get_correct_gff(self) -> str:
        return "Homo_sapiens.GRCh38.108.gff3"

    def get_correct_gtf(self) -> str:
        return "Homo_sapiens.GRCh38.108.gtf"

    def get_correct_fasta(self) -> str:
        return "Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa"

    def test_download_ncrna_file(self) -> None:
        ensembl_loader = cast(GenomicDatabaseDownloaderEnsembl, self.loader)
        result = ensembl_loader.fetch_genomic_entity(self.entity, ["fasta"], sequence_nature="ncrna")

        self.assert_file_name_is_included(result["fetched_files"]["fasta"], "Homo_sapiens.GRCh38.ncrna.fa")


class GenomicRegionGeneratorBase(unittest.TestCase, ABC):
    expected_generation_behavior: dict[str, str] = {}
    expected_header_values: dict[str, RegionHeaderSpec] = {}

    def setUp(self) -> None:
        self.tmp_path = os.path.join(os.getcwd(), "tmp_genomic_region_generator")
        self.fasta_parser = FastaParser()
        self.region_generator = self.setup_region_generator()

    def tearDown(self) -> None:
        try:
            shutil.rmtree(self.tmp_path)
        except:
            pass

    @abstractmethod
    def setup_region_generator(self) -> CustomGenomicRegionGenerator:
        pass

    def _run_generation_test(self, region_type: str, generator_func: Callable) -> None:
        expected_generation_behavior = self.expected_generation_behavior[region_type]

        if expected_generation_behavior == "error":
            with self.assertRaises(ConfigurationError):
                generator_func()
            return

        if expected_generation_behavior == "warning":
            with self.assertWarns(Warning):
                result = generator_func()
        else:
            result = generator_func()

        expected_header_values = self.expected_header_values[region_type]

        self.assertTrue(
            self.fasta_parser.check_fasta_format(result), f"error: wrong file format for file: {result}"
        )

        if expected_header_values:
            fasta_sequences = SeqIO.index(result, "fasta")
            additional_info: dict[str, Any] = {}
            for idx in fasta_sequences:
                region, ai, coordinates = self.fasta_parser.parse_fasta_header(idx)
                # in case parse_additional_info in parse_fasta_header is False, additional_info
                # is returned as str and mypy complains about that, therefore cast here
                additional_info = cast(dict[str, Any], ai)
                if region == expected_header_values["region"]:
                    break
            for key, expected in expected_header_values["additional_info"].items():
                with self.subTest(header_key=key):
                    self.assertEqual(additional_info.get(key), expected)
            for key, expected in expected_header_values["coordinates"].items():
                with self.subTest(header_key=key):
                    self.assertEqual(coordinates.get(key), expected)

    def test_gene(self) -> None:
        self._run_generation_test("gene", self.region_generator.get_sequence_gene)

    def test_exon(self) -> None:
        self._run_generation_test("exon", self.region_generator.get_sequence_exon)

    def test_exon_exon_junction(self) -> None:
        self._run_generation_test(
            "exon_exon_junction", lambda: self.region_generator.get_sequence_exon_exon_junction(block_size=50)
        )

    def test_CDS(self) -> None:
        self._run_generation_test("CDS", self.region_generator.get_sequence_CDS)

    def test_UTR(self) -> None:
        self._run_generation_test(
            "UTR", lambda: self.region_generator.get_sequence_UTR(five_prime=True, three_prime=True)
        )

    def test_intergenic(self) -> None:
        self._run_generation_test("intergenic", self.region_generator.get_sequence_intergenic)

    def test_intron(self) -> None:
        self._run_generation_test("intron", self.region_generator.get_sequence_intron)

    def test_chromosome_length_file_keeps_first_row_with_header_none(self) -> None:
        # chrom.sizes / bedtools -g files are headerless; header=0 would drop the first chromosome.
        file_fasta = os.path.join(self.tmp_path, "genome.fna")
        with open(file_fasta, "w") as handle:
            handle.write(">chr1\n")
            handle.write("A" * 1000 + "\n")
            handle.write(">chr2\n")
            handle.write("C" * 2000 + "\n")

        self.region_generator.sequence_file = file_fasta
        file_chromosome_length = self.region_generator._get_chromosome_length()

        # Same load path as get_sequence_intergenic when a chromosome has no genes.
        chromosome_length = self.region_generator._read_chromosome_length_file(file_chromosome_length)
        assert list(chromosome_length["seqid"]) == ["chr1", "chr2"]
        assert list(chromosome_length["length"]) == [1000, 2000]

        # Same file is the -g argument to bedtools complement via get_complement_regions.
        file_bed_in = os.path.join(self.tmp_path, "annotation_in.bed")
        file_bed_out = os.path.join(self.tmp_path, "annotation_out.bed")
        with open(file_bed_in, "w") as handle:
            handle.write("chr1\t100\t200\n")

        get_complement_regions(file_bed_in, file_chromosome_length, file_bed_out)

        complement = self.region_generator.bed_parser.read_bed(file_bed_out, names=["seqid", "start", "end"])
        assert list(complement["seqid"]) == ["chr1", "chr1"]
        assert list(complement["start"]) == [0, 200]
        assert list(complement["end"]) == [100, 1000]


class TestGenomicRegionGeneratorNCBI(GenomicRegionGeneratorBase):
    expected_generation_behavior = {
        "gene": "pass",
        "exon": "pass",
        "exon_exon_junction": "pass",
        "CDS": "pass",
        "UTR": "pass",
        "intergenic": "pass",
        "intron": "pass",
    }
    expected_header_values = EXPECTED_HEADER_VALUES_HUMAN_NCBI

    def setup_region_generator(self) -> CustomGenomicRegionGenerator:

        return CustomGenomicRegionGenerator(
            FILE_ANNOTATION_NCBI,
            FILE_SEQUENCE_NCBI,
            files_source=METADATA_NCBI["files_source"],
            species=METADATA_NCBI["species"],
            annotation_release=METADATA_NCBI["annotation_release"],
            genome_assembly=METADATA_NCBI["genome_assembly"],
            dir_output=self.tmp_path,
        )


class TestGenomicRegionGeneratorEnsembl(GenomicRegionGeneratorBase):
    expected_generation_behavior = {
        "gene": "pass",
        "exon": "pass",
        "exon_exon_junction": "pass",
        "CDS": "pass",
        "UTR": "pass",
        "intergenic": "pass",
        "intron": "pass",
    }
    expected_header_values = EXPECTED_HEADER_VALUES_HUMAN_ENSEMBL

    def setup_region_generator(self) -> CustomGenomicRegionGenerator:

        return CustomGenomicRegionGenerator(
            FILE_ANNOTATION_ENSEMBL,
            FILE_SEQUENCE_ENSEMBL,
            files_source=METADATA_ENSEMBL["files_source"],
            species=METADATA_ENSEMBL["species"],
            annotation_release=METADATA_ENSEMBL["annotation_release"],
            genome_assembly=METADATA_ENSEMBL["genome_assembly"],
            dir_output=self.tmp_path,
        )


class TestGenomicRegionGeneratorMouseNCBI(GenomicRegionGeneratorBase):
    expected_generation_behavior = {
        "gene": "pass",
        "exon": "warning",
        "exon_exon_junction": "warning",
        "CDS": "warning",
        "UTR": "warning",
        "intergenic": "pass",
        "intron": "pass",
    }
    expected_header_values = EXPECTED_HEADER_VALUES_MOUSE_NCBI

    def setup_region_generator(self) -> CustomGenomicRegionGenerator:

        return CustomGenomicRegionGenerator(
            FILE_ANNOTATION_MOUSE_NCBI,
            FILE_SEQUENCE_MOUSE_NCBI,
            files_source=METADATA_MOUSE_NCBI["files_source"],
            species=METADATA_MOUSE_NCBI["species"],
            annotation_release=METADATA_MOUSE_NCBI["annotation_release"],
            genome_assembly=METADATA_MOUSE_NCBI["genome_assembly"],
            dir_output=self.tmp_path,
        )


class TestGenomicRegionGeneratorBacteriaNCBI(GenomicRegionGeneratorBase):
    expected_generation_behavior = {
        "gene": "pass",
        "exon": "warning",
        "exon_exon_junction": "error",
        "CDS": "warning",
        "UTR": "error",
        "intergenic": "pass",
        "intron": "error",
    }
    expected_header_values = EXPECTED_HEADER_VALUES_BACTERIA_NCBI

    def setup_region_generator(self) -> CustomGenomicRegionGenerator:

        return CustomGenomicRegionGenerator(
            FILE_ANNOTATION_BACTERIA_NCBI,
            FILE_SEQUENCE_BACTERIA_NCBI,
            files_source=METADATA_BACTERIA_NCBI["files_source"],
            species=METADATA_BACTERIA_NCBI["species"],
            annotation_release=METADATA_BACTERIA_NCBI["annotation_release"],
            genome_assembly=METADATA_BACTERIA_NCBI["genome_assembly"],
            dir_output=self.tmp_path,
        )


class TestOligoSequenceGenerator(unittest.TestCase):
    def setUp(self) -> None:
        self.tmp_path = os.path.join(os.getcwd(), "tmp_oligo_sequence_generator")

        self.oligo_database_1 = OligoDatabase(dir_output=self.tmp_path, database_name="db1")
        self.oligo_database_2 = OligoDatabase(dir_output=self.tmp_path, database_name="db2")
        self.oligo_sequence_generator = OligoSequenceGenerator(dir_output=self.tmp_path)
        self.fasta_parser = FastaParser()

        self.sequence_type = "oligo"

    def tearDown(self) -> None:
        shutil.rmtree(self.tmp_path)

    def test_create_sequences_random(self) -> None:
        file_fasta_random_seqs1 = self.oligo_sequence_generator.create_sequences_random(
            filename_out="random_sequences1",
            length_sequences=30,
            num_sequences=100,
            name_sequences="random_sequences1",
            base_alphabet_with_probability={"A": 0.1, "C": 0.3, "G": 0.4, "T": 0.2},
        )
        assert (
            self.fasta_parser.check_fasta_format(file_fasta_random_seqs1) == True
        ), f"error: wrong file format for file: {file_fasta_random_seqs1}"

        self.oligo_database_1.load_database_from_fasta(
            files_fasta=file_fasta_random_seqs1,
            database_overwrite=True,
            sequence_type=self.sequence_type,
            region_ids=None,
        )
        length_property = LengthProperty()
        calculator = PropertyCalculator(properties=[length_property])
        self.oligo_database_1 = calculator.apply(
            oligo_database=self.oligo_database_1, sequence_type=self.sequence_type, n_jobs=1
        )

        num_sequences = self.oligo_database_1.get_oligoid_list(region_ids="random_sequences1")
        length_sequence = self.oligo_database_1.get_oligo_property_value(
            property=f"length_{self.sequence_type}",
            flatten=True,
            region_id="random_sequences1",
            oligo_id="random_sequences1::1",
        )

        assert len(num_sequences) == 100, "error: wrong number sequences created"

        assert length_sequence == 30, "error: wrong sequence length"

        assert check_if_dna_sequence(
            self.oligo_database_1.database["random_sequences1"]["random_sequences1::50"]["oligo"]
        ), "error: the craeted sequence is not a DNA seuqnece"

    def test_create_sequences_sliding_window(self) -> None:

        # test if warning is raised if no oligos can be created because of too short
        # exon-exon-junction sequences
        with self.assertWarns(Warning):
            file_fasta_exon_exon_junctions_short = (
                self.oligo_sequence_generator.create_sequences_sliding_window(
                    files_fasta_in=FILE_NCBI_EXON_EXON_JUNCTIONS_SHORT,
                    length_interval_sequences=(30, 31),
                    region_ids=[
                        "AARS1",
                    ],
                )
            )

        # test sliding window without strides
        file_fasta_exons = self.oligo_sequence_generator.create_sequences_sliding_window(
            files_fasta_in=FILE_NCBI_EXONS,
            length_interval_sequences=(30, 31),
            region_ids=[
                "AARS1",
                "DECR2",
                "FAM234A",
                "RHBDF1",
                "WASIR2",
            ],
        )
        file_fasta_exon_exon_junctions = self.oligo_sequence_generator.create_sequences_sliding_window(
            files_fasta_in=FILE_NCBI_EXON_EXON_JUNCTIONS,
            length_interval_sequences=(30, 31),
            region_ids=[
                "ABAT",
            ],
        )

        self.oligo_database_1.load_database_from_fasta(
            files_fasta=file_fasta_exons + file_fasta_exon_exon_junctions,
            database_overwrite=True,
            sequence_type="oligo",
            region_ids=["AARS1", "ABAT"],
        )
        length_property = LengthProperty()
        calculator = PropertyCalculator(properties=[length_property])
        self.oligo_database_1 = calculator.apply(
            oligo_database=self.oligo_database_1, sequence_type=self.sequence_type, n_jobs=1
        )

        length_sequence = self.oligo_database_1.get_oligo_property_value(
            property=f"length_{self.sequence_type}",
            flatten=True,
            region_id="AARS1",
            oligo_id="AARS1::1",
        )

        sequence = self.oligo_database_1.get_oligo_property_value(
            property="oligo",
            flatten=True,
            region_id="ABAT",
            oligo_id="ABAT::1",
        )

        num_start_value = self.oligo_database_1.get_oligo_property_value(
            property="start",
            flatten=True,
            region_id="ABAT",
            oligo_id="ABAT::1",
        )

        assert isinstance(num_start_value, list), "error: start property is not a list"
        num_start = len(num_start_value)

        assert "AARS1" in self.oligo_database_1.database.keys(), "error: region missing"

        assert length_sequence == 30, "error: wrong sequence length"

        assert sequence == "ATGGGCGCGTTCCATGGGAGGACCATGGGT", "error: wrong sequence"

        assert num_start == 2

        assert check_if_dna_sequence(
            self.oligo_database_1.database["AARS1"]["AARS1::50"]["oligo"]
        ), "error: the craeted sequence is not a DNA seuqnece"

        # test sliding window with strides
        file_fasta_exons_stride = self.oligo_sequence_generator.create_sequences_sliding_window(
            files_fasta_in=FILE_NCBI_EXONS,
            length_interval_sequences=(30, 31),
            stride=2,
            region_ids=[
                "AARS1",
                "DECR2",
                "FAM234A",
                "RHBDF1",
                "WASIR2",
            ],
        )

        self.oligo_database_2.load_database_from_fasta(
            files_fasta=file_fasta_exons_stride,
            database_overwrite=True,
            sequence_type="oligo",
            region_ids="AARS1",
        )

        # the first sequence is always the same no matter the stride,
        # and with a stride of 1, the third sequence should be equal to the second sequence with a stride of 2
        sequence_no_stride = self.oligo_database_1.get_oligo_property_value(
            property="oligo",
            flatten=True,
            region_id="AARS1",
            oligo_id="AARS1::3",
        )
        sequence_stride = self.oligo_database_2.get_oligo_property_value(
            property="oligo",
            flatten=True,
            region_id="AARS1",
            oligo_id="AARS1::2",
        )

        assert sequence_no_stride == sequence_stride, "error: sequences are not equal"
