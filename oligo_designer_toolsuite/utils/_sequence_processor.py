############################################
# imports
############################################

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
) -> None:
    """
    Extracts sequences from a reference FASTA file based on regions specified in a BED file using `bedtools getfasta`.

    :param file_bed: Path to the BED file containing regions of interest.
    :type file_bed: str
    :param file_reference_fasta: Path to the reference FASTA file.
    :type file_reference_fasta: str
    :param file_fasta: Output FASTA file path to store extracted sequences.
    :type file_fasta: str
    :param split: Whether to split the sequences. Given BED12 input, extract and concatenate the sequences from the BED “blocks” (e.g., exons), defaults to False.
    :type split: bool
    :param strand: Whether to force strandedness. If the feature occupies the antisense strand, the sequence will be reverse complemented, defaults to False.
    :type strand: bool
    :param nameOnly: Whether to use the name field for the FASTA header, defaults to False.
    :type nameOnly: bool
    :param name: Whether to use the name field and coordinates for the FASTA header, defaults to False.
    :type name: bool
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


def get_complement_regions(file_bed_in: str, file_chromosome_length: str, file_bed_out: str) -> None:
    """
    Generates the complement of genomic regions specified in a BED file using `bedtools complement`.

    :param file_bed_in: Path to the input BED file containing regions of interest.
    :type file_bed_in: str
    :param file_chromosome_length: Path to the file specifying chromosome lengths.
    :type file_chromosome_length: str
    :param file_bed_out: Path to the output BED file where complement regions will be saved.
    :type file_bed_out: str
    """
    cmd = "bedtools complement"
    cmd += " -i " + file_bed_in
    cmd += " -g " + file_chromosome_length
    cmd += " -L "
    cmd += " > " + file_bed_out

    process = Popen(cmd, shell=True).wait()
