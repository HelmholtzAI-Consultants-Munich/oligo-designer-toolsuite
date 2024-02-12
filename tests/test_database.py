############################################
# imports
############################################


import pytest
from Bio.SeqUtils import MeltingTemp as mt

from oligo_designer_toolsuite.database import (
    CustomGenomicRegionGenerator,
    OligoDatabase,
    ReferenceDatabase,
)
from oligo_designer_toolsuite.oligo_property_filter import (
    GCContentFilter,
    HardMaskedSequenceFilter,
    MeltingTemperatureNNFilter,
    PropertyFilter,
)
from oligo_designer_toolsuite.utils import check_fasta_format, parse_fasta_header

############################################
# Global Parameters
############################################

genes = ["AARS1", "DECR2", "FAM234A", "RHBDF1", "WASIR2"]
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

Tm_parameters = {
    "check": True,
    "strict": True,
    "c_seq": None,
    "shift": 0,
    "nn_table": getattr(mt, "DNA_NN3"),
    "tmm_table": getattr(mt, "DNA_TMM1"),
    "imm_table": getattr(mt, "DNA_IMM1"),
    "de_table": getattr(mt, "DNA_DE1"),
    "dnac1": 50,  # [nM]
    "dnac2": 0,
    "selfcomp": False,
    "dNTPs": 0,
    "saltcorr": 7,
    "Na": 1.25,  # [mM]
    "K": 75,  # [mM]
    "Tris": 20,  # [mM]
    "Mg": 10,  # [mM]
}

Tm_correction_parameters = {
    "DMSO": 0,
    "DMSOfactor": 0.75,
    "fmdfactor": 0.65,
    "fmdmethod": 1,
    "GC": None,
    "fmd": 20,
}

############################################
# Tests
############################################


# tmp_path generates a temporary directory unique to the test invocation, created in the base temporary directory
# for more information see here: https://docs.pytest.org/en/6.2.x/tmpdir.html
@pytest.fixture(scope="session")
def file_ncbi_transcriptome(tmpdir_factory):
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
    ncbi_transcriptome = (
        region_generator_ncbi.generate_transcript_reduced_representation()
    )
    return ncbi_transcriptome


def test_region_generator_ncbi(tmpdir_factory):
    """Test region generation of genome, transcriptome and CDS with NCBI files."""
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

    ncbi_genome = region_generator_ncbi.generate_genome()
    assert check_fasta_format(ncbi_genome) == True, "error: wrong file format"

    ncbi_transcriptome = (
        region_generator_ncbi.generate_transcript_reduced_representation()
    )
    assert check_fasta_format(ncbi_transcriptome) == True, "error: wrong file format"

    ncbi_CDS = region_generator_ncbi.generate_CDS_reduced_representation()
    assert check_fasta_format(ncbi_CDS) == True, "error: wrong file format"


def test_region_generator_ensemble(tmpdir_factory):
    base_temp = tmpdir_factory.getbasetemp()
    """Test region generation of genome, transcriptome and CDS with ensemble files."""
    region_generator_ensembl = CustomGenomicRegionGenerator(
        annotation_file_ensemble,
        sequence_file_ensemble,
        files_source=metadata_ensemble["files_source"],
        species=metadata_ensemble["species"],
        annotation_release=metadata_ensemble["annotation_release"],
        genome_assembly=metadata_ensemble["genome_assembly"],
        dir_output=base_temp,
    )

    ensembl_genome = region_generator_ensembl.generate_genome()
    assert check_fasta_format(ensembl_genome) == True, "error: wrong file format"

    ensembl_transcriptome = (
        region_generator_ensembl.generate_transcript_reduced_representation()
    )
    assert check_fasta_format(ensembl_transcriptome) == True, "error: wrong file format"

    ensembl_CDS = region_generator_ensembl.generate_CDS_reduced_representation()
    assert check_fasta_format(ensembl_CDS) == True, "error: wrong file format"


def test_reference_database(file_ncbi_transcriptome):
    """Test creation of reference database as well as load, write and filter functionalities."""
    reference = ReferenceDatabase(
        file_ncbi_transcriptome,
        metadata_ncbi,
    )

    reference.load_fasta_into_database()
    reference.filter_database(["AARS1"])
    for entry in reference.database:
        (
            region,
            _,
            _,
        ) = parse_fasta_header(entry.id)
        assert region == "AARS1", f"error: this region {region} should be filtered out."
    file_fasta_database = reference.write_fasta_from_database(
        filename="filtered_databse"
    )
    assert check_fasta_format(file_fasta_database) == True, "error: wrong file format"


def test_oligo_database(file_ncbi_transcriptome):
    """Test creation of oligo database as well as save, load and write to fasta functionalities."""
    oligos = OligoDatabase(
        min_oligos_per_region=0,
        metadata=metadata_ncbi,
        n_jobs=2,
    )
    oligos.create_database(
        file_fasta=file_ncbi_transcriptome,
        oligo_length_min=90,
        oligo_length_max=90,
        region_ids=genes,
    )
    database = oligos.database

    # check if database changes when saved and loaded
    file_database, file_metadata = oligos.write_database()
    oligos.load_database(file_database)
    assert oligos.metadata == {}, "error: metadata not cleared"
    for oligo_id in database[genes[0]].keys():
        assert (
            database[genes[0]][oligo_id] == oligos.database[genes[0]][oligo_id]
        ), f"error: the database changes when it is saved and loaded again. "

    # check if function write a correct fasta file
    file_fasta = oligos.write_fasta_from_database()
    assert check_fasta_format(file_fasta) == True, "error: wrong file format"


def test_oligo_database_filters(file_ncbi_transcriptome):
    """Test base filter functions on oligo database."""

    def _get_sequences_from_database(database):
        sequences = []
        for region_id, oligo in database.items():
            for oligo_id, oligo_attributes in oligo.items():
                sequences.append(oligo_attributes["sequence"])
        sequences.sort()  # needed to compare
        return sequences

    oligos = OligoDatabase(
        min_oligos_per_region=0,
        metadata=metadata_ncbi,
        n_jobs=2,
    )

    oligos.create_database(
        file_fasta=file_ncbi_transcriptome,
        oligo_length_min=90,
        oligo_length_max=90,
        region_ids=genes,
    )

    masked_sequences = HardMaskedSequenceFilter(mask="N")
    GC_content = GCContentFilter(GC_content_min=40, GC_content_max=60)
    melting_temperature = MeltingTemperatureNNFilter(
        Tm_min=52,
        Tm_max=67,
        Tm_parameters=Tm_parameters,
        Tm_chem_correction_parameters=Tm_correction_parameters,
    )

    filters = [masked_sequences, GC_content, melting_temperature]
    property_filter = PropertyFilter(filters=filters)
    property_filter.apply(oligos)

    sequences = _get_sequences_from_database(oligos.database)
    assert (
        len(sequences) == 7
    ), "error: filtering did not return the correct number of sequences"
