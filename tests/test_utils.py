############################################
# imports
############################################

import os
from pathlib import Path

from Bio import SeqIO

from oligo_designer_toolsuite.utils import (
    FtpLoaderEnsembl,
    FtpLoaderNCBI,
    GffParser,
    check_gff_format,
    check_fasta_format,
    check_tsv_format,
    get_sequence_from_annotation,
    merge_fasta,
    parse_fasta_header,
)

############################################
# Global Parameters
############################################


############################################
# Tests
############################################


def test_data_parser(tmp_path):
    """Test if data parser functionalities work correctly."""
    # test file format checkers
    file_gff = (
        "tests/data/annotations/custom_GCF_000001405.40_GRCh38.p14_genomic_chr16.gff"
    )
    res = check_gff_format(file_gff)
    assert (
        res == True
    ), f"error: gff file format checker did not recognize gff file {file_gff}"

    file_gtf = (
        "tests/data/annotations/custom_GCF_000001405.40_GRCh38.p14_genomic_chr16.gtf"
    )
    res = check_gff_format(file_gtf)
    assert (
        res == True
    ), f"error: gff file format checker did not recognize gtf file {file_gtf}"

    file_fasta = (
        "tests/data/annotations/custom_GCF_000001405.40_GRCh38.p14_genomic_chr16.fna"
    )
    res = check_fasta_format(file_fasta)
    assert (
        res == True
    ), f"error: fasta file format checker did not recognize fasta file {file_fasta}"

    file_tsv = "tests/data/annotations/custom_GCF_000001405.40_GRCh38.p14_genomic_chr16.gtf.tsv"
    res = check_tsv_format(file_tsv)
    assert (
        res == True
    ), f"error: tsv file format checker did not recognize tsv file {file_tsv}"

    # test sequence extraction from annotation and fasta file
    file_bed = (
        "tests/data/annotations/custom_GCF_000001405.40_GRCh38.p14_genomic_chr16.bed"
    )
    file_reference_fasta = (
        "tests/data/annotations/custom_GCF_000001405.40_GRCh38.p14_genomic_chr16.fna"
    )
    file_fasta = os.path.join(tmp_path, "test_ann2seq_function.fna")
    get_sequence_from_annotation(
        file_bed, file_reference_fasta, file_fasta, split=False, strand=True, name=True
    )
    res = check_fasta_format(file_fasta)
    assert res == True, f"error: the created sequence file is not a fasta file"

    files_fasta = [file_reference_fasta, file_reference_fasta]
    file_merged_fasta = os.path.join(tmp_path, "test_merge_fasta_function.fna")
    merge_fasta(files_fasta, file_merged_fasta)
    res = len(list(SeqIO.parse(file_merged_fasta, "fasta")))
    assert res == 2, f"error: the fasta files were not merged correctly"

    header = "ARPG3::transcript_id=XM4581;exon_id=XM4581_exon1::16:70265537-70265662(-)"
    region, additional_information, coordinates = parse_fasta_header(header)
    assert region == "ARPG3", "error: wrong region parsed"
    assert coordinates["chromosome"] == ["16"], "error: wrong chrom parsed"
    assert coordinates["start"] == [70265537], "error: wrong start parsed"
    assert coordinates["end"] == [70265662], "error: wrong end parsed"
    assert coordinates["strand"] == ["-"], "error: wrong strand parsed"
    assert (
        additional_information == "transcript_id=XM4581;exon_id=XM4581_exon1"
    ), "error: wrong additional information parsed"


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


def test_GFF_parser():
    """Test of GFF/GTF parser parses file correctly."""
    ##### Test GFF3 parsing
    parser = GffParser()
    file_gff = (
        "tests/data/annotations/custom_GCF_000001405.40_GRCh38.p14_genomic_chr16.gff"
    )

    dataframe_gff = parser.read_gff(file_gff, target_lines=10)
    assert dataframe_gff.shape[1] == 23, "error: GFF3 dataframe not correctly loaded"

    ##### Test GTF parsing
    parser = GffParser()
    file_gtf = (
        "tests/data/annotations/custom_GCF_000001405.40_GRCh38.p14_genomic_chr16.gtf"
    )

    dataframe_gtf = parser.read_gff(file_gtf, target_lines=10)
    assert dataframe_gtf.shape[1] == 20, "error: GTF dataframe not correctly loaded"
