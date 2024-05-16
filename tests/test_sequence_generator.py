############################################
# imports
############################################

from pathlib import Path

import pytest

from oligo_designer_toolsuite.sequence_generator import (
    CustomGenomicRegionGenerator,
    FtpLoaderEnsembl,
    FtpLoaderNCBI,
)
from oligo_designer_toolsuite.utils import FastaParser

############################################
# Global Parameters
############################################

region_ids = [
    "AARS1",
    "DECR2",
    "FAM234A",
    "RHBDF1",
    "WASIR2",
    "this_gene_does_not_exist",
]
annotation_file_ensemble = "data/tests/annotations/custom_Homo_sapiens.GRCh38.108.chr16.gtf"
sequence_file_ensemble = "data/tests/annotations/custom_Homo_sapiens.GRCh38.dna_sm.chromosome.16.fa"

annotation_file_ncbi = "data/tests/annotations/custom_GCF_000001405.40_GRCh38.p14_genomic_chr16.gtf"
sequence_file_ncbi = "data/tests/annotations/custom_GCF_000001405.40_GRCh38.p14_genomic_chr16.fna"

metadata_ncbi = {
    "files_source": "NCBI",
    "species": "Homo_sapiens",
    "annotation_release": "110",
    "genome_assembly": "GRCh38",
}

metadata_ensemble = {
    "files_source": "Ensembl",
    "species": "Homo_sapiens",
    "annotation_release": "108",
    "genome_assembly": "GRCh38",
}

file_oligo_attributes = "/data/tests/sequences_oligo_attributes.fna"

############################################
# Tests
############################################


# tmp_path generates a temporary directory unique to the test invocation, created in the base temporary directory
# for more information see here: https://docs.pytest.org/en/6.2.x/tmpdir.html
@pytest.fixture(scope="session")
def file_ncbi_exons(tmpdir_factory):
    base_temp = tmpdir_factory.getbasetemp()
    region_generator_ncbi = CustomGenomicRegionGenerator(
        annotation_file_ncbi,
        sequence_file_ncbi,
        files_source=metadata_ncbi["files_source"],
        species=metadata_ncbi["species"],
        annotation_release=metadata_ncbi["annotation_release"],
        genome_assembly=metadata_ncbi["genome_assembly"],
        dir_output=base_temp,
    )

    ncbi_exons = region_generator_ncbi.get_sequence_exon()
    return ncbi_exons


def test_ftp_loader_ncbi(tmp_path):
    """Test if ftp download for NCBI works correctly."""
    # Parameters
    taxon = "vertebrate_mammalian"  # taxon the species belongs to
    species = "Homo_sapiens"
    annotation_release = "110"

    # initialize
    loader_ncbi = FtpLoaderNCBI(tmp_path, taxon, species, annotation_release)

    # retrieve files
    file_gff, annotation_release, assembly_name = loader_ncbi.download_files("gff")

    assert Path(file_gff).name == "GCF_000001405.40_GRCh38.p14_genomic.gff", "error: wrong file downloaded"
    assert annotation_release == "110", "error: wrong annotation release retrieved"
    assert assembly_name == "GRCh38.p14", "error: wrong assembly name retrieved"

    file_gff, annotation_release, assembly_name = loader_ncbi.download_files("gtf")

    assert Path(file_gff).name == "GCF_000001405.40_GRCh38.p14_genomic.gtf", "error: wrong file downloaded"
    assert annotation_release == "110", "error: wrong annotation release retrieved"
    assert assembly_name == "GRCh38.p14", "error: wrong assembly name retrieved"

    file_fasta, annotation_release, assembly_name = loader_ncbi.download_files("fasta")

    assert Path(file_fasta).name == "GCF_000001405.40_GRCh38.p14_genomic.fna", "error: wrong file downloaded"
    assert annotation_release == "110", "error: wrong annotation release retrieved"
    assert assembly_name == "GRCh38.p14", "error: wrong assembly name retrieved"

    ## Test if download for releases > 110 works -> they changed the folder structure there
    # Parameters
    taxon = "vertebrate_mammalian"  # taxon the species belongs to
    species = "Homo_sapiens"
    annotation_release = "current"

    # initialize and retrieve files
    loader_ncbi = FtpLoaderNCBI(tmp_path, taxon, species, annotation_release)
    file_gff, annotation_release, assembly_name = loader_ncbi.download_files("gff")


def test_ftp_loader_ensemble(tmp_path):
    """Test if ftp download for Ensemble works correctly."""
    # Parameters
    species = "homo_sapiens"
    annotation_release = "108"

    # initialize
    loader_ensemble = FtpLoaderEnsembl(tmp_path, species, annotation_release)

    # retrieve files
    file_gff, annotation_release, assembly_name = loader_ensemble.download_files("gff")

    assert Path(file_gff).name == "Homo_sapiens.GRCh38.108.gff3", "error: wrong file downloaded"
    assert annotation_release == "108", "error: wrong annotation release retrieved"
    assert assembly_name == "GRCh38", "error: wrong assembly name retrieved"

    file_gff, annotation_release, assembly_name = loader_ensemble.download_files("gtf")

    assert Path(file_gff).name == "Homo_sapiens.GRCh38.108.gtf", "error: wrong file downloaded"
    assert annotation_release == "108", "error: wrong annotation release retrieved"
    assert assembly_name == "GRCh38", "error: wrong assembly name retrieved"

    file_fasta, annotation_release, assembly_name = loader_ensemble.download_files("fasta")

    assert (
        Path(file_fasta).name == "Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa"
    ), f"error: wrong file: {file_fasta} downloaded"
    assert annotation_release == "108", "error: wrong annotation release retrieved"
    assert assembly_name == "GRCh38", "error: wrong assembly name retrieved"

    file_fasta, annotation_release, assembly_name = loader_ensemble.download_files(
        "fasta", sequence_nature="ncrna"
    )

    assert Path(file_fasta).name == "Homo_sapiens.GRCh38.ncrna.fa", "error: wrong file downloaded"
    assert annotation_release == "108", "error: wrong annotation release retrieved"
    assert assembly_name == "GRCh38", "error: wrong assembly name retrieved"


def test_region_generator_ncbi(tmpdir_factory):
    """Test region generation of genome, transcriptome and CDS with NCBI files."""
    base_temp = tmpdir_factory.getbasetemp()
    fasta_parser = FastaParser()
    region_generator_ncbi = CustomGenomicRegionGenerator(
        annotation_file_ncbi,
        sequence_file_ncbi,
        files_source=metadata_ncbi["files_source"],
        species=metadata_ncbi["species"],
        annotation_release=metadata_ncbi["annotation_release"],
        genome_assembly=metadata_ncbi["genome_assembly"],
        dir_output=base_temp,
    )

    ncbi_genes = region_generator_ncbi.get_sequence_gene()
    assert (
        fasta_parser.check_fasta_format(ncbi_genes) == True
    ), f"error: wrong file format for file: {ncbi_genes}"

    ncbi_exons = region_generator_ncbi.get_sequence_exon()
    assert (
        fasta_parser.check_fasta_format(ncbi_exons) == True
    ), f"error: wrong file format for file: {ncbi_exons}"

    ncbi_CDS = region_generator_ncbi.get_sequence_CDS()
    assert fasta_parser.check_fasta_format(ncbi_CDS) == True, f"error: wrong file format for file: {ncbi_CDS}"

    ncbi_UTR = region_generator_ncbi.get_sequence_UTR(five_prime=True, three_prime=True)
    assert fasta_parser.check_fasta_format(ncbi_UTR) == True, f"error: wrong file format for file: {ncbi_UTR}"

    ncbi_junction = region_generator_ncbi.get_sequence_exon_exon_junction(block_size=50)
    assert (
        fasta_parser.check_fasta_format(ncbi_junction) == True
    ), f"error: wrong file format for file: {ncbi_junction}"

    ncbi_intergenic = region_generator_ncbi.get_sequence_intergenic()
    assert (
        fasta_parser.check_fasta_format(ncbi_intergenic) == True
    ), f"error: wrong file format for file: {ncbi_intergenic}"

    ncbi_introns = region_generator_ncbi.get_sequence_intron()
    assert (
        fasta_parser.check_fasta_format(ncbi_introns) == True
    ), f"error: wrong file format for file: {ncbi_introns}"


def test_region_generator_ensemble(tmpdir_factory):
    """Test region generation of genome, transcriptome and CDS with ensemble files."""
    base_temp = tmpdir_factory.getbasetemp()
    fasta_parser = FastaParser()
    region_generator_ensembl = CustomGenomicRegionGenerator(
        annotation_file_ensemble,
        sequence_file_ensemble,
        files_source=metadata_ensemble["files_source"],
        species=metadata_ensemble["species"],
        annotation_release=metadata_ensemble["annotation_release"],
        genome_assembly=metadata_ensemble["genome_assembly"],
        dir_output=base_temp,
    )

    ensembl_gene = region_generator_ensembl.get_sequence_gene()
    assert (
        fasta_parser.check_fasta_format(ensembl_gene) == True
    ), f"error: wrong file format for file: {ensembl_gene}"

    ensembl_exon = region_generator_ensembl.get_sequence_exon()
    assert (
        fasta_parser.check_fasta_format(ensembl_exon) == True
    ), f"error: wrong file format for file: {ensembl_exon}"

    ensembl_CDS = region_generator_ensembl.get_sequence_CDS()
    assert (
        fasta_parser.check_fasta_format(ensembl_CDS) == True
    ), f"error: wrong file format for file: {ensembl_CDS}"

    ensembl_UTR = region_generator_ensembl.get_sequence_UTR(five_prime=True, three_prime=True)
    assert (
        fasta_parser.check_fasta_format(ensembl_UTR) == True
    ), f"error: wrong file format for file: {ensembl_UTR}"

    ensembl_junction = region_generator_ensembl.get_sequence_exon_exon_junction(block_size=50)
    assert (
        fasta_parser.check_fasta_format(ensembl_junction) == True
    ), f"error: wrong file format for file: {ensembl_junction}"

    ensembl_intergenic = region_generator_ensembl.get_sequence_intergenic()
    assert (
        fasta_parser.check_fasta_format(ensembl_intergenic) == True
    ), f"error: wrong file format for file: {ensembl_intergenic}"

    ensembl_intron = region_generator_ensembl.get_sequence_intron()
    assert (
        fasta_parser.check_fasta_format(ensembl_intron) == True
    ), f"error: wrong file format for file: {ensembl_intron}"
