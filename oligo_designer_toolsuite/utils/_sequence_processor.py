############################################
# imports
############################################

import os
from subprocess import Popen

############################################
# Collection of utility functions
############################################


def get_sequence_from_annotation(
    file_bed,
    file_reference_fasta,
    file_fasta,
    split=False,
    strand=False,
    nameOnly=False,
    name=False,
):
    """Generate a FASTA file with sequences based on a BED file and a reference FASTA file.

    :param file_bed: The input BED file containing genomic coordinates.
    :type file_bed: str
    :param file_reference_fasta: The reference FASTA file containing genomic sequences.
    :type file_reference_fasta: str
    :param file_fasta: The output FASTA file to store the generated sequences.
    :type file_fasta: str
    :param split: If True, split exons by introns in the BED file.
    :type split: bool, optional
    :param strand: If True, strand-specific information is considered.
    :type strand: bool, optional
    :param nameOnly: If True, only the name field is used as the FASTA header.
    :type nameOnly: bool, optional
    :param name: If True, include the name field in the FASTA header.
    :type name: bool, optional
    """
    cmd = "bedtools getfasta"
    cmd += " -fi " + file_reference_fasta
    cmd += " -bed " + file_bed
    cmd += " -fo " + file_fasta
    if split:
        cmd += " -split"
    if strand:
        cmd += " -s"
    if nameOnly:
        cmd += " -nameOnly"
    if name:
        cmd += " -name"

    process = Popen(cmd, shell=True).wait()


def get_complement_regions(file_bed_in, file_chromosome_length, file_bed_out):
    """Generate the complement regions of the input BED file with respect to the given chromosome lengths.

    :param file_bed_in: The input BED file containing genomic coordinates.
    :type file_bed_in: str
    :param file_chromosome_length: The file containing chromosome lengths.
    :type file_chromosome_length: str
    :param file_bed_out: The output BED file to store the complement regions.
    :type file_bed_out: str
    """
    cmd = "bedtools complement"
    cmd += " -i " + file_bed_in
    cmd += " -g " + file_chromosome_length
    cmd += " -L "
    cmd += " > " + file_bed_out

    process = Popen(cmd, shell=True).wait()
