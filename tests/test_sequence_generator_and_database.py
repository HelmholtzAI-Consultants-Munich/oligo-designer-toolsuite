############################################
# imports
############################################

from pathlib import Path

import pytest

from oligo_designer_toolsuite.database import OligoDatabase, ReferenceDatabase
from oligo_designer_toolsuite.sequence_generator import (
    CustomGenomicRegionGenerator,
    FtpLoaderEnsembl,
    FtpLoaderNCBI,
    OligoSequenceGenerator,
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
annotation_file_ensemble = "data/annotations/custom_Homo_sapiens.GRCh38.108.chr16.gtf"
sequence_file_ensemble = (
    "data/annotations/custom_Homo_sapiens.GRCh38.dna_rm.primary_assembly_chr16.fa"
)
annotation_file_ncbi = (
    "data/annotations/custom_GCF_000001405.40_GRCh38.p14_genomic_chr16.gtf"
)
sequence_file_ncbi = (
    "data/annotations/custom_GCF_000001405.40_GRCh38.p14_genomic_chr16.fna"
)

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

    assert (
        Path(file_gff).name == "GCF_000001405.40_GRCh38.p14_genomic.gff"
    ), "error: wrong file downloaded"
    assert annotation_release == "110", "error: wrong annotation release retrieved"
    assert assembly_name == "GRCh38.p14", "error: wrong assembly name retrieved"

    file_gff, annotation_release, assembly_name = loader_ncbi.download_files("gtf")

    assert (
        Path(file_gff).name == "GCF_000001405.40_GRCh38.p14_genomic.gtf"
    ), "error: wrong file downloaded"
    assert annotation_release == "110", "error: wrong annotation release retrieved"
    assert assembly_name == "GRCh38.p14", "error: wrong assembly name retrieved"

    file_fasta, annotation_release, assembly_name = loader_ncbi.download_files("fasta")

    assert (
        Path(file_fasta).name == "GCF_000001405.40_GRCh38.p14_genomic.fna"
    ), "error: wrong file downloaded"
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

    assert (
        Path(file_gff).name == "Homo_sapiens.GRCh38.108.gff3"
    ), "error: wrong file downloaded"
    assert annotation_release == "108", "error: wrong annotation release retrieved"
    assert assembly_name == "GRCh38", "error: wrong assembly name retrieved"

    file_gff, annotation_release, assembly_name = loader_ensemble.download_files("gtf")

    assert (
        Path(file_gff).name == "Homo_sapiens.GRCh38.108.gtf"
    ), "error: wrong file downloaded"
    assert annotation_release == "108", "error: wrong annotation release retrieved"
    assert assembly_name == "GRCh38", "error: wrong assembly name retrieved"

    file_fasta, annotation_release, assembly_name = loader_ensemble.download_files(
        "fasta"
    )

    assert (
        Path(file_fasta).name == "Homo_sapiens.GRCh38.dna_rm.primary_assembly.fa"
    ), "error: wrong file downloaded"
    assert annotation_release == "108", "error: wrong annotation release retrieved"
    assert assembly_name == "GRCh38", "error: wrong assembly name retrieved"

    file_fasta, annotation_release, assembly_name = loader_ensemble.download_files(
        "fasta", sequence_nature="ncrna"
    )

    assert (
        Path(file_fasta).name == "Homo_sapiens.GRCh38.ncrna.fa"
    ), "error: wrong file downloaded"
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
    ), "error: wrong file format"

    ncbi_exons = region_generator_ncbi.get_sequence_exon()
    assert (
        fasta_parser.check_fasta_format(ncbi_exons) == True
    ), "error: wrong file format"

    ncbi_CDS = region_generator_ncbi.get_sequence_CDS()
    assert fasta_parser.check_fasta_format(ncbi_CDS) == True, "error: wrong file format"

    ncbi_UTR = region_generator_ncbi.get_sequence_UTR(five_prime=True, three_prime=True)
    assert fasta_parser.check_fasta_format(ncbi_UTR) == True, "error: wrong file format"

    ncbi_junction = region_generator_ncbi.get_sequence_exon_exon_junction(block_size=50)
    assert (
        fasta_parser.check_fasta_format(ncbi_junction) == True
    ), "error: wrong file format"

    ncbi_intergenic = region_generator_ncbi.get_sequence_intergenic()
    assert (
        fasta_parser.check_fasta_format(ncbi_intergenic) == True
    ), "error: wrong file format"

    ncbi_introns = region_generator_ncbi.get_sequence_intron()
    assert (
        fasta_parser.check_fasta_format(ncbi_introns) == True
    ), "error: wrong file format"


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
    ), "error: wrong file format"

    ensembl_exon = region_generator_ensembl.get_sequence_exon()
    assert (
        fasta_parser.check_fasta_format(ensembl_exon) == True
    ), "error: wrong file format"

    ensembl_CDS = region_generator_ensembl.get_sequence_CDS()
    assert (
        fasta_parser.check_fasta_format(ensembl_CDS) == True
    ), "error: wrong file format"

    ensembl_UTR = region_generator_ensembl.get_sequence_UTR(
        five_prime=True, three_prime=True
    )
    assert (
        fasta_parser.check_fasta_format(ensembl_UTR) == True
    ), "error: wrong file format"

    ensembl_junction = region_generator_ensembl.get_sequence_exon_exon_junction(
        block_size=50
    )
    assert (
        fasta_parser.check_fasta_format(ensembl_junction) == True
    ), "error: wrong file format"

    ensembl_intergenic = region_generator_ensembl.get_sequence_intergenic()
    assert (
        fasta_parser.check_fasta_format(ensembl_intergenic) == True
    ), "error: wrong file format"

    ensembl_intron = region_generator_ensembl.get_sequence_intron()
    assert (
        fasta_parser.check_fasta_format(ensembl_intron) == True
    ), "error: wrong file format"


def test_reference_database(file_ncbi_exons):
    """Test creation of reference database as well as load, write and filter functionalities."""
    fasta_parser = FastaParser()
    reference = ReferenceDatabase()
    reference.load_metadata(metadata=metadata_ncbi)
    reference.load_sequences_fom_fasta(
        file_fasta=file_ncbi_exons, database_overwrite=True
    )
    reference.load_sequences_fom_fasta(
        file_fasta=file_ncbi_exons, database_overwrite=False
    )

    reference.filter_database("AARS1")
    for entry in reference.database:
        (
            region,
            _,
            _,
        ) = fasta_parser.parse_fasta_header(entry.id)
        assert region == "AARS1", f"error: this region {region} should be filtered out."

    file_fasta_database = reference.write_database_to_fasta(filename="filtered_databse")
    file_metadata_database = reference.write_metadata_to_yaml(
        filename="filtered_databse"
    )
    assert (
        fasta_parser.check_fasta_format(file_fasta_database) == True
    ), "error: wrong file format"


def test_oligo_database(file_ncbi_exons):
    """Test creation of oligo database as well as save, load and write to fasta functionalities."""
    oligo_sequence_generator = OligoSequenceGenerator()

    file_fasta_random_seqs1 = oligo_sequence_generator.create_sequences_random(
        filename_out="random_sequences1",
        length_sequences=30,
        num_sequences=100,
        name_sequences="random_sequences1",
        base_alphabet_with_probability={"A": 0.1, "C": 0.3, "G": 0.4, "T": 0.2},
    )
    file_fasta_random_seqs2 = oligo_sequence_generator.create_sequences_random(
        filename_out="random_sequences2",
        length_sequences=15,
        num_sequences=3,
        name_sequences="random_sequences2",
    )
    file_fasta_exons = oligo_sequence_generator.create_sequences_sliding_window(
        filename_out="sliding_window_sequences",
        file_fasta_in=file_ncbi_exons,
        length_interval_sequences=(30, 31),
    )

    oligos = OligoDatabase(
        min_oligos_per_region=2, write_regions_with_insufficient_oligos=True
    )
    oligos2 = OligoDatabase(
        min_oligos_per_region=4, write_regions_with_insufficient_oligos=True
    )

    # check if we can load sequences into the database for oligo or target and if we can merge databases
    oligos.load_metadata(metadata_ncbi)
    oligos.load_sequences_from_fasta(
        file_fasta_in=file_fasta_random_seqs1,
        sequence_type="oligo",
        region_ids=["random_sequences1"],
        database_overwrite=True,
    )
    oligos.load_sequences_from_fasta(
        file_fasta_in=file_fasta_random_seqs2,
        sequence_type="oligo",
        database_overwrite=False,
    )
    oligos.load_sequences_from_fasta(
        file_fasta_in=file_fasta_exons,
        sequence_type="target",
        region_ids=region_ids,
        database_overwrite=False,
    )

    # check if calculation of number of targeted transcripts and isoform consensus works
    oligos.calculate_num_targeted_transcripts()
    oligos.calculate_isoform_consensus()

    # check if removale of regions works
    oligos.remove_regions_with_insufficient_oligos("database_generation")
    assert len(oligos.database.keys()) == (
        len(region_ids) - 1 + 2
    ), "error: wrong number of regions in database"

    # check if save and load works
    file_database, file_metadata = oligos.save_database(
        region_ids="random_sequences2", filename_out="database_random_sequences2"
    )

    oligos2.load_metadata(file_metadata)
    oligos2.load_database(file_database, database_overwrite=True)
    oligos2.load_sequences_from_fasta(
        file_fasta_in=file_fasta_random_seqs1,
        sequence_type="oligo",
        database_overwrite=False,
    )

    # check if removale of regions works
    oligos2.remove_regions_with_insufficient_oligos("database_generation")
    assert (
        len(oligos2.database.keys()) == 1
    ), "error: wrong number of regions in database"

    # check if we get the correct number of sequences returned
    list_sequences = oligos2.get_sequence_list()
    assert len(list_sequences) == 100, "error: wrong number of sequences in database"

    oligos2.write_database_to_fasta()
