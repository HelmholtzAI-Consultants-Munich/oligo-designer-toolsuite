############################################
# imports
############################################

import os
import shutil
import unittest
from abc import abstractmethod
from pathlib import Path
from typing import cast

from oligo_designer_toolsuite.database import OligoDatabase
from oligo_designer_toolsuite.oligo_property_calculator import (
    LengthProperty,
    PropertyCalculator,
)
from oligo_designer_toolsuite.sequence_generator import (
    CustomGenomicRegionGenerator,
    FtpLoaderEnsembl,
    FtpLoaderNCBI,
    OligoSequenceGenerator,
)
from oligo_designer_toolsuite.utils import FastaParser, check_if_dna_sequence
from oligo_designer_toolsuite.validation.models._general import BaseProbabilities

############################################
# Setup
############################################

FILE_ANNOTATION_ENSEMBL = "tests/data/annotations/custom_Homo_sapiens.GRCh38.108.chr16.gtf"
FILE_SEQUENCE_ENSEMBL = "tests/data/annotations/custom_Homo_sapiens.GRCh38.dna_sm.chromosome.16.fa"

FILE_ANNOTATION_NCBI = "tests/data/annotations/custom_GCF_000001405.40_GRCh38.p14_genomic_chr16.gtf"
FILE_SEQUENCE_NCBI = "tests/data/annotations/custom_GCF_000001405.40_GRCh38.p14_genomic_chr16.fna"

METDATA_NCBI = {
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

FILE_NCBI_EXONS = "tests/data/genomic_regions/sequences_ncbi_exons.fna"
FILE_NCBI_EXON_EXON_JUNCTIONS = "tests/data/genomic_regions/sequences_ABAT_ncbi_exon_exon_junctions.fna"
FILE_NCBI_EXON_EXON_JUNCTIONS_SHORT = (
    "tests/data/genomic_regions/sequences_AARS1_ncbi_exon_exon_junctions_short.fna"
)

############################################
# Tests
############################################


class FTPLoaderDownloadBase:
    def setUp(self) -> None:
        self.tmp_path = os.path.join(os.getcwd(), "tmp_ftp_loader")
        os.makedirs(self.tmp_path, exist_ok=True)
        self.loader = self.setup_ftp_loader()

    def tearDown(self) -> None:
        shutil.rmtree(self.tmp_path)

    @abstractmethod
    def setup_ftp_loader(self) -> FtpLoaderNCBI | FtpLoaderEnsembl:
        pass

    def test_download(self) -> None:
        _, _, _ = self.loader.download_files("fasta")


class TestFTPLoaderNCBICurrent(FTPLoaderDownloadBase, unittest.TestCase):
    def setup_ftp_loader(self) -> FtpLoaderNCBI:
        # Parameters
        taxon = "vertebrate_mammalian"  # taxon the species belongs to
        species = "Homo_sapiens"
        annotation_release = "current"

        return FtpLoaderNCBI(self.tmp_path, taxon, species, annotation_release)


class TestFTPLoaderEnsemblCurrent(FTPLoaderDownloadBase, unittest.TestCase):
    def setup_ftp_loader(self) -> FtpLoaderEnsembl:
        # Parameters
        species = "homo_sapiens"
        annotation_release = "current"

        return FtpLoaderEnsembl(self.tmp_path, species, annotation_release)


class FTPLoaderFilesBase:
    def setUp(self) -> None:
        self.tmp_path = os.path.join(os.getcwd(), "tmp_ftp_loader")
        self.loader = self.setup_ftp_loader()

    def tearDown(self) -> None:
        try:
            shutil.rmtree(self.tmp_path)
        except:
            pass

    @abstractmethod
    def setup_ftp_loader(self) -> FtpLoaderNCBI | FtpLoaderEnsembl:
        pass

    @abstractmethod
    def get_correct_metadata(self) -> tuple[str, str]:
        pass

    @abstractmethod
    def get_correct_gff(self) -> str:
        pass

    @abstractmethod
    def get_correct_gtf(self) -> str:
        pass

    @abstractmethod
    def get_correct_fasta(self) -> str:
        pass

    @abstractmethod
    def get_true_asserts(self) -> bool:
        pass

    def test_metadata_loader(self) -> None:
        _, annotation_release, assembly_name = self.loader.download_files("fasta")
        self.annotation_release, self.assembly_name = self.get_correct_metadata()
        assert annotation_release == self.annotation_release, "error: wrong annotation release retrieved"
        assert assembly_name == self.assembly_name, "error: wrong assembly name retrieved"

    def test_gff_loader(self) -> None:
        file_gff, _, _ = self.loader.download_files("gff")
        assert Path(file_gff).name == self.get_correct_gff(), "error: wrong file downloaded"

    def test_gtf_loader(self) -> None:
        file_gtf, _, _ = self.loader.download_files("gtf")
        assert Path(file_gtf).name == self.get_correct_gtf(), "error: wrong file downloaded"

    def test_fasta_loader(self) -> None:
        file_fasta, _, _ = self.loader.download_files("fasta")
        assert Path(file_fasta).name == self.get_correct_fasta(), "error: wrong file downloaded"


class TestFTPLoaderNCBIOldAnnotations(FTPLoaderFilesBase, unittest.TestCase):
    def setup_ftp_loader(self) -> FtpLoaderNCBI:
        # Parameters
        taxon = "vertebrate_mammalian"  # taxon the species belongs to
        species = "Homo_sapiens"
        annotation_release = "110"

        return FtpLoaderNCBI(self.tmp_path, taxon, species, annotation_release)

    def get_correct_metadata(self) -> tuple[str, str]:
        annotation_release = "110"
        assembly_name = "GRCh38.p14"

        return annotation_release, assembly_name

    def get_correct_gff(self) -> str:
        return "GCF_000001405.40_GRCh38.p14_genomic.gff"

    def get_correct_gtf(self) -> str:
        return "GCF_000001405.40_GRCh38.p14_genomic.gtf"

    def get_correct_fasta(self) -> str:
        return "GCF_000001405.40_GRCh38.p14_genomic.fna"


class TestFTPLoaderEnsemblOldAnnotations(FTPLoaderFilesBase, unittest.TestCase):
    def setup_ftp_loader(self) -> FtpLoaderEnsembl:
        # Parameters
        species = "homo_sapiens"
        annotation_release = "108"

        return FtpLoaderEnsembl(self.tmp_path, species, annotation_release)

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
        ensembl_loader = cast(FtpLoaderEnsembl, self.loader)
        file_fasta, _, _ = ensembl_loader.download_files("fasta", sequence_nature="ncrna")
        assert Path(file_fasta).name == "Homo_sapiens.GRCh38.ncrna.fa", "error: wrong file downloaded"


class GenomicRegionGeneratorBase:
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

    def test_gene(self) -> None:
        genes = self.region_generator.get_sequence_gene()
        assert (
            self.fasta_parser.check_fasta_format(genes) == True
        ), f"error: wrong file format for file: {genes}"

    def test_exon(self) -> None:
        exon = self.region_generator.get_sequence_exon()
        assert (
            self.fasta_parser.check_fasta_format(exon) == True
        ), f"error: wrong file format for file: {exon}"

    def test_exon_exon_junction(self) -> None:
        exon_exon_junction = self.region_generator.get_sequence_exon_exon_junction(block_size=50)
        assert (
            self.fasta_parser.check_fasta_format(exon_exon_junction) == True
        ), f"error: wrong file format for file: {exon_exon_junction}"

    def test_CDS(self) -> None:
        cds = self.region_generator.get_sequence_CDS()
        assert self.fasta_parser.check_fasta_format(cds) == True, f"error: wrong file format for file: {cds}"

    def test_UTR(self) -> None:
        utr = self.region_generator.get_sequence_UTR(five_prime=True, three_prime=True)
        assert self.fasta_parser.check_fasta_format(utr) == True, f"error: wrong file format for file: {utr}"

    def test_intergenic(self) -> None:
        intergenic = self.region_generator.get_sequence_intergenic()
        assert (
            self.fasta_parser.check_fasta_format(intergenic) == True
        ), f"error: wrong file format for file: {intergenic}"

    def test_introns(self) -> None:
        introns = self.region_generator.get_sequence_intron()
        assert (
            self.fasta_parser.check_fasta_format(introns) == True
        ), f"error: wrong file format for file: {introns}"


class TestGenomicRegionGeneratorNCBI(GenomicRegionGeneratorBase, unittest.TestCase):
    def setup_region_generator(self) -> CustomGenomicRegionGenerator:

        return CustomGenomicRegionGenerator(
            FILE_ANNOTATION_NCBI,
            FILE_SEQUENCE_NCBI,
            files_source=METDATA_NCBI["files_source"],
            species=METDATA_NCBI["species"],
            annotation_release=METDATA_NCBI["annotation_release"],
            genome_assembly=METDATA_NCBI["genome_assembly"],
            dir_output=self.tmp_path,
        )


class TestGenomicRegionGeneratorEnsembl(GenomicRegionGeneratorBase, unittest.TestCase):
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
            base_alphabet_with_probability=BaseProbabilities(A=0.1, C=0.3, G=0.4, T=0.2),
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
