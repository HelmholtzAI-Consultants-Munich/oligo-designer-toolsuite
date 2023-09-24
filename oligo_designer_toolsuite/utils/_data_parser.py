############################################
# imports
############################################

import os
import csv
import shutil

# import pybedtools

from subprocess import Popen
from Bio import SeqIO

from ..utils._gff_parser import GffParser

############################################
# Collection of utility functions
############################################


def check_gff_format(file):
    """File check. Is the file a GFF3/ or GTF file?

    :param file: Path to file.
    :type file: str
    :return: Returns True if file has GFF3 or GTF format.
    :rtype: bool
    """
    parser = GffParser()
    gtf = parser.read_gff(file, target_lines=100)
    return any(gtf)


def check_fasta_format(file):
    """File check. Is the file a fasta file?

    :param file: Path to file.
    :type file: str
    :return: Returns True if file has fasta format.
    :rtype: bool
    """
    # taken from https://stackoverflow.com/questions/44293407/how-can-i-check-whether-a-given-file-is-fasta
    with open(file, "r") as handle:
        fasta = SeqIO.parse(handle, "fasta")
        return any(fasta)  # False when `fasta` is empty, i.e. wasn't a FASTA file


def check_tsv_format(file):
    """File check. Is the file a tsv file?

    :param file: Path to file.
    :type file: str
    :return: Returns True if file has tsv format.
    :rtype: bool
    """
    with open(file, "r") as tsv:
        read_tsv = csv.reader(tsv, delimiter="\t")
        return any(read_tsv)


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


def merge_fasta(files_fasta, file_merged_fasta):
    """Merge two fast files by concatenating file content.

    :param files_fasta: List of fasta files to be concatenated.
    :type files_fasta: list of str
    :param file_merged_fasta: File name of merged fasta file.
    :type file_merged_fasta: str
    """

    if files_fasta == []:
        raise ValueError("No fasta files provided for merge.")

    with open(file_merged_fasta, "wb") as handle_DB:
        for file in files_fasta:
            if os.path.exists(file):
                shutil.copyfileobj(open(file, "rb"), handle_DB)
                # shutil.rmtree(file)
                # os.remove(file)
            else:
                raise ValueError(f"Fasta file {file} does not exist!")


def parse_fasta_header(header):
    """Helper function to parse the header of a sequence entry in a fasta file.

    :param header: Header of sequence entry, starting with '>'.
    :type header: string
    :return: Parsed header information, i.e. region, coordinates (chromosome, start, end, strand), and (optional) additional information
    :rtype: str, dict, str
    """
    header = header.split("::")
    region = header[0]

    if len(header) == 1:
        coordinates = {
            "chromosome": [None],
            "start": [None],
            "end": [None],
            "strand": [None],
        }

    else:
        header_coordinates = header[-1].split(";")
        coordinates = {}
        for header_coordinate in header_coordinates:
            coordinates.setdefault("chromosome", []).append(header_coordinate.split(":")[0])
            coordinates.setdefault("start", []).append(int(header_coordinate.split(":")[1].split("-")[0]))
            coordinates.setdefault("end", []).append(
                int(header_coordinate.split(":")[1].split("-")[1].split("(")[0])
            )
            coordinates.setdefault("strand", []).append(header_coordinate.split("(")[1].split(")")[0])

    additional_information = {}

    if len(header) > 2:
        info_list = header[1]

        if ";" in info_list:
            info_list = info_list.split(";")

            for infos in info_list:
                key_values = infos.split(",")

                for key_value in key_values:
                    key, value = key_value.split("=")

                    if key in additional_information:
                        additional_information[key].append(value)
                    else:
                        additional_information[key] = [value]

        else:
            key, value = info_list.split("=")
            additional_information[key] = [value]

    return region, additional_information, coordinates
