############################################
# imports
############################################

import os
import shutil
import unittest

from pathlib import Path
from abc import abstractmethod

from oligo_designer_toolsuite.sequence_generator import (
    CustomGenomicRegionGenerator,
    FtpLoaderEnsembl,
    FtpLoaderNCBI,
)
from oligo_designer_toolsuite.utils import FastaParser

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


# returns error: ftplib.error_perm: 550 Failed to change directory.
# class TestFTPLoaderEnsemblCurrent(FTPLoaderDownloadBase, unittest.TestCase):
#     def setup_ftp_loader(self):
#         # Parameters
#         species = "homo_sapiens"
#         annotation_release = "current"

#         return FtpLoaderEnsembl(self.tmp_path, species, annotation_release)


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


class TestFTPLoaderEnsembl(FTPLoaderFilesBase, unittest.TestCase):
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
