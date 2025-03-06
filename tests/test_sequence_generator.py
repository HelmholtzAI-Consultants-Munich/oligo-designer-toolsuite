############################################
# imports
############################################

import os
import shutil
import unittest
from abc import abstractmethod
from pathlib import Path

from oligo_designer_toolsuite.database import OligoAttributes, OligoDatabase
from oligo_designer_toolsuite.sequence_generator import (
    CustomGenomicRegionGenerator,
    FtpLoaderEnsembl,
    FtpLoaderNCBI,
    OligoSequenceGenerator,
)
from oligo_designer_toolsuite.utils import FastaParser, check_if_dna_sequence

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

############################################
# Tests
############################################


class FTPLoaderDownloadBase:
    def setUp(self):
        self.tmp_path = os.path.join(os.getcwd(), "tmp_ftp_loader")
        os.makedirs(self.tmp_path, exist_ok=True)
        self.loader = self.setup_ftp_loader()

    def tearDown(self):
        shutil.rmtree(self.tmp_path)

    @abstractmethod
    def setup_ftp_loader(self):
        pass

    def test_download(self):
        _, _, _ = self.loader.download_files("fasta")


class TestFTPLoaderNCBICurrent(FTPLoaderDownloadBase, unittest.TestCase):
    def setup_ftp_loader(self):
        # Parameters
        taxon = "vertebrate_mammalian"  # taxon the species belongs to
        species = "Homo_sapiens"
        annotation_release = "current"

        return FtpLoaderNCBI(self.tmp_path, taxon, species, annotation_release)


class TestFTPLoaderEnsemblCurrent(FTPLoaderDownloadBase, unittest.TestCase):
    def setup_ftp_loader(self):
        # Parameters
        species = "homo_sapiens"
        annotation_release = "current"

        return FtpLoaderEnsembl(self.tmp_path, species, annotation_release)


class FTPLoaderFilesBase:
    def setUp(self):
        self.tmp_path = os.path.join(os.getcwd(), "tmp_ftp_loader")
        self.loader = self.setup_ftp_loader()

    def tearDown(self):
        try:
            shutil.rmtree(self.tmp_path)
        except:
            pass

    @abstractmethod
    def setup_ftp_loader(self):
        pass

    @abstractmethod
    def get_true_asserts(self):
        pass

    def test_metadata_loader(self):
        _, annotation_release, assembly_name = self.loader.download_files("fasta")
        self.annotation_release, self.assembly_name = self.get_correct_metadata()
        assert annotation_release == self.annotation_release, "error: wrong annotation release retrieved"
        assert assembly_name == self.assembly_name, "error: wrong assembly name retrieved"

    def test_gff_loader(self):
        file_gff, _, _ = self.loader.download_files("gff")
        assert Path(file_gff).name == self.get_correct_gff(), "error: wrong file downloaded"

    def test_gtf_loader(self):
        file_gtf, _, _ = self.loader.download_files("gtf")
        assert Path(file_gtf).name == self.get_correct_gtf(), "error: wrong file downloaded"

    def test_fasta_loader(self):
        file_fasta, _, _ = self.loader.download_files("fasta")
        assert Path(file_fasta).name == self.get_correct_fasta(), "error: wrong file downloaded"


class TestFTPLoaderNCBIOldAnnotations(FTPLoaderFilesBase, unittest.TestCase):
    def setup_ftp_loader(self):
        # Parameters
        taxon = "vertebrate_mammalian"  # taxon the species belongs to
        species = "Homo_sapiens"
        annotation_release = "110"

        return FtpLoaderNCBI(self.tmp_path, taxon, species, annotation_release)

    def get_correct_metadata(self):
        annotation_release = "110"
        assembly_name = "GRCh38.p14"

        return annotation_release, assembly_name

    def get_correct_gff(self):
        return "GCF_000001405.40_GRCh38.p14_genomic.gff"

    def get_correct_gtf(self):
        return "GCF_000001405.40_GRCh38.p14_genomic.gtf"

    def get_correct_fasta(self):
        return "GCF_000001405.40_GRCh38.p14_genomic.fna"


class TestFTPLoaderEnsemblOldAnnotations(FTPLoaderFilesBase, unittest.TestCase):
    def setup_ftp_loader(self):
        # Parameters
        species = "homo_sapiens"
        annotation_release = "108"

        return FtpLoaderEnsembl(self.tmp_path, species, annotation_release)

    def get_correct_metadata(self):
        annotation_release = "108"
        assembly_name = "GRCh38"

        return annotation_release, assembly_name

    def get_correct_gff(self):
        return "Homo_sapiens.GRCh38.108.gff3"

    def get_correct_gtf(self):
        return "Homo_sapiens.GRCh38.108.gtf"

    def get_correct_fasta(self):
        return "Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa"

    def test_download_ncrna_file(self):
        file_fasta, _, _ = self.loader.download_files("fasta", sequence_nature="ncrna")
        assert Path(file_fasta).name == "Homo_sapiens.GRCh38.ncrna.fa", "error: wrong file downloaded"


class GenomicRegionGeneratorBase:
    def setUp(self):
        self.tmp_path = os.path.join(os.getcwd(), "tmp_genomic_region_generator")
        self.fasta_parser = FastaParser()
        self.region_generator = self.setup_region_generator()

    def tearDown(self):
        try:
            shutil.rmtree(self.tmp_path)
        except:
            pass

    @abstractmethod
    def setup_region_generator(self):
        pass

    def test_gene(self):
        genes = self.region_generator.get_sequence_gene()
        assert (
            self.fasta_parser.check_fasta_format(genes) == True
        ), f"error: wrong file format for file: {genes}"

    def test_exon(self):
        exon = self.region_generator.get_sequence_exon()
        assert (
            self.fasta_parser.check_fasta_format(exon) == True
        ), f"error: wrong file format for file: {exon}"

    def test_exon_exon_junction(self):
        exon_exon_junction = self.region_generator.get_sequence_exon_exon_junction(block_size=50)
        assert (
            self.fasta_parser.check_fasta_format(exon_exon_junction) == True
        ), f"error: wrong file format for file: {exon_exon_junction}"

    def test_CDS(self):
        cds = self.region_generator.get_sequence_CDS()
        assert self.fasta_parser.check_fasta_format(cds) == True, f"error: wrong file format for file: {cds}"

    def test_UTR(self):
        utr = self.region_generator.get_sequence_UTR(five_prime=True, three_prime=True)
        assert self.fasta_parser.check_fasta_format(utr) == True, f"error: wrong file format for file: {utr}"

    def test_intergenic(self):
        intergenic = self.region_generator.get_sequence_intergenic()
        assert (
            self.fasta_parser.check_fasta_format(intergenic) == True
        ), f"error: wrong file format for file: {intergenic}"

    def test_introns(self):
        introns = self.region_generator.get_sequence_intron()
        assert (
            self.fasta_parser.check_fasta_format(introns) == True
        ), f"error: wrong file format for file: {introns}"


class TestGenomicRegionGeneratorNCBI(GenomicRegionGeneratorBase, unittest.TestCase):
    def setup_region_generator(self):

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
    def setup_region_generator(self):

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
    def setUp(self):
        self.tmp_path = os.path.join(os.getcwd(), "tmp_oligo_sequence_generator")

        self.oligo_database = OligoDatabase()
        self.oligo_attributes = OligoAttributes()
        self.oligo_sequence_generator = OligoSequenceGenerator(dir_output=self.tmp_path)
        self.fasta_parser = FastaParser()

    def tearDown(self):
        shutil.rmtree(self.tmp_path)

    def test_create_sequences_random(self):
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

        self.oligo_database.load_database_from_fasta(
            files_fasta=file_fasta_random_seqs1,
            database_overwrite=True,
            sequence_type="oligo",
            region_ids=None,
        )

        assert (
            len(self.oligo_database.database["random_sequences1"].keys()) == 100
        ), "error: wrong number sequences created"
        self.oligo_database = self.oligo_attributes.calculate_oligo_length(oligo_database=self.oligo_database)
        assert (
            self.oligo_database.get_oligo_attribute_value(
                attribute="length",
                flatten=True,
                region_id="random_sequences1",
                oligo_id="random_sequences1::1",
            )
            == 30
        ), "error: wrong sequence length"
        assert check_if_dna_sequence(
            self.oligo_database.database["random_sequences1"]["random_sequences1::50"]["oligo"]
        ), "error: the craeted sequence is not a DNA seuqnece"

    def test_create_sequences_sliding_window(self):
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

        self.oligo_database.load_database_from_fasta(
            files_fasta=file_fasta_exons,
            database_overwrite=True,
            sequence_type="oligo",
            region_ids="AARS1",
        )

        assert "AARS1" in self.oligo_database.database.keys(), "error: region missing"
        self.oligo_database = self.oligo_attributes.calculate_oligo_length(oligo_database=self.oligo_database)
        assert (
            self.oligo_database.get_oligo_attribute_value(
                attribute="length",
                flatten=True,
                region_id="AARS1",
                oligo_id="AARS1::1",
            )
            == 30
        ), "error: wrong sequence length"
        assert check_if_dna_sequence(
            self.oligo_database.database["AARS1"]["AARS1::50"]["oligo"]
        ), "error: the craeted sequence is not a DNA seuqnece"
