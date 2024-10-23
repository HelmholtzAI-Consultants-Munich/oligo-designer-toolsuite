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
    cmd = "bedtools complement"
    cmd += " -i " + file_bed_in
    cmd += " -g " + file_chromosome_length
    cmd += " -L "
    cmd += " > " + file_bed_out

    process = Popen(cmd, shell=True).wait()
