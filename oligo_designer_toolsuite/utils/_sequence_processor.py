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
):
    """
    Get sequence for regions annotated in file_bed using file_reference_fasta and save the output as a fasta file in file_fasta.

    :param file_bed: Path to bed file with annotated genomic regions.
    :type file_bed: string
    :param file_reference_fasta: Path to fasta file with reference sequence, e.g. transcriptome.
    :type file_reference_fasta: string
    :param file_fasta: Path to fasta file where retrieved sequences are written to.
    :type file_fasta: string
    :param split: Given BED12 fmt., extract and concatenate the sequences from the BED "blocks" (e.g., exons), defaults to False.
    :type split: bool
    :param strand: Force strandedness. If the feature occupies the antisense strand, the sequence will be reverse complemented, defaults to False.
    :type strand: bool
    :param nameOnly: Use the name field for the FASTA header, defaults to False.
    :type nameOnly: bool
    :param name: Use the name field and coordinates for the FASTA header, defaults to False.
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


def get_complement_regions(file_bed_in, file_chromosome_length, file_bed_out):
    """Get all regions in a genome that are not covered by at least one interval in the input file.
    Chromosomes that are in the chromosome length file but not in the input file will be suppressed.

    :param file_bed_in: Input bed file.
    :type file_bed_in: str
    :param file_chromosome_length: Length of each chromosome.
    :type file_chromosome_length: str
    :param file_bed_out: Output bed file with complementary regions to input bed file.
    :type file_bed_out: str
    """
    cmd = "bedtools complement"
    cmd += " -i " + file_bed_in
    cmd += " -g " + file_chromosome_length
    cmd += " -L "
    cmd += " > " + file_bed_out

    process = Popen(cmd, shell=True).wait()
