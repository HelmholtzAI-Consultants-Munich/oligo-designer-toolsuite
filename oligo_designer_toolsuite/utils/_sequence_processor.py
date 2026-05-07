############################################
# imports
############################################

import os
import subprocess
from collections import Counter

from Bio import SeqIO

from ._checkers_and_helpers import cast_to_list
from ._sequence_parser import FastaParser

############################################
# Collection of utility functions
############################################


def get_sequence_from_annotation(
    file_bed: str,
    file_reference_fasta: str,
    file_fasta: str,
    split: bool = False,
    strand: bool = False,
    nameOnly: bool = False,
    name: bool = False,
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
    args = ["bedtools", "getfasta", "-fi", file_reference_fasta, "-bed", file_bed, "-fo", file_fasta]
    if split:
        args.append("-split")
    if strand:
        args.append("-s")
    if nameOnly:
        args.append("-nameOnly")
    if name:
        args.append("-name")

    subprocess.run(args, check=True, stdout=subprocess.DEVNULL)


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
    args = ["bedtools", "complement", "-i", file_bed_in, "-g", file_chromosome_length, "-L"]

    # redirect stdout to output file
    with open(file_bed_out, "w") as f:
        subprocess.run(args, stdout=f, check=True)


def get_intersection(file_A: str, file_B: list[str] | str, file_bed_out: str) -> None:
    """
    Compute the intersection between genomic regions in two BED files.

    This function uses `bedtools intersect` to find overlapping regions between the input BED files.
    The output is saved in a specified BED file.

    :param file_A: Path to the first BED file.
    :type file_A: str
    :param file_B: list of paths to the second BED file(s).
    :type file_B: list[str] | str
    :param file_bed_out: Path to the output BED file where the intersection results will be saved.
    :type file_bed_out: str
    """
    file_B_list: list[str] = cast_to_list(file_B)

    args = [
        "bedtools",
        "intersect",
        "-wa",
        "-wb",
        "-bed",
        "-a",
        file_A,
        "-b",
        *file_B_list,
    ]

    # redirect stdout to output file
    with open(file_bed_out, "w") as f:
        subprocess.run(args, stdout=f, check=True)


def append_nucleotide_to_sequences(input_fasta: str, nucleotide: str) -> str:
    """
    Appends a specific nucleotide to each sequence in a FASTA file.

    :param input_fasta: Path to the input FASTA file.
    :type input_fasta: str
    :param nucleotide: Nucleotide to append to each sequence (e.g., 'A', 'T', 'C', 'G').
    :type nucleotide: str
    :return: Path to the output FASTA file with modified sequences.
    :rtype: str
    """
    base, ext = os.path.splitext(input_fasta)
    output_fasta = f"{base}_modified{ext}"
    # Open the input and output FASTA files
    with open(input_fasta, "r") as infile, open(output_fasta, "w") as outfile:
        for record in SeqIO.parse(infile, "fasta"):
            # Append the nucleotide to the sequence
            record.seq = record.seq + nucleotide

            # Write the modified record to the output file
            SeqIO.write(record, outfile, "fasta")

    return output_fasta


def remove_index_files(file_reference: str, dir_output: str) -> None:
    """
    Removes the index files for a reference file.

    This function removes all files in the output directory that start with
    `file_reference.` (the reference file basename followed by a dot). This
    approach handles various index file naming patterns:

    - FASTA files: removes .fai index files (e.g., reference.fna.fai)
    - VCF files: removes .csi and .tbi index files (e.g., reference.vcf.gz.csi)
    - BLAST databases: removes .nhr, .nin, .nsq index files (e.g., reference.fna.nhr)
    - Bowtie/Bowtie2 indexes: removes .ebwt and .bt2 index files, including
      multi-part indexes (e.g., reference.fna.1.ebwt, reference.fna.rev.1.ebwt)

    **Note**: This function removes ALL files starting with `file_reference.`, not
    just known index file extensions. This ensures compatibility with complex index
    naming patterns (e.g., Bowtie files with multiple dots) but means any file
    matching this prefix pattern will be removed. The reference file itself is
    never removed.

    :param file_reference: The base name of the reference file.
    :type file_reference: str
    :param dir_output: The directory where the reference file is located.
    :type dir_output: str
    """
    file_reference_basename = os.path.basename(file_reference)

    # Index files are extensions of the input file: input_file.extension
    # So we remove any file that starts with input_file. (with a dot)
    prefix = file_reference_basename + "."

    for root, _, files in os.walk(dir_output):
        for file in files:
            # Remove files that start with the prefix (index files)
            # but not the reference file itself
            if file.startswith(prefix):
                file_path = os.path.join(root, file)
                os.remove(file_path)


def count_kmer_abundance(
    files_fasta: str | list[str],
    k: int | tuple[int, int] | list[int],
) -> dict[int, dict[str, float]]:
    """
    Counts k-mer abundances across multiple FASTA files and returns fractional abundances.

    This function iterates through one or more FASTA files, extracts k-mers of specified length(s),
    counts their occurrences, and calculates fractional abundances (each k-mer count divided by
    the total k-mer count for that k value). The reported fractional abundances are reported across
    all FASTA files.The function supports processing a single k value, a range of k values, or a
    list of specific k values.

    :param files_fasta: Path(s) to FASTA file(s) to process. Can be a single file path (str) or
        a list of file paths (list[str]).
    :type files_fasta: str | list[str]
    :param k: K-mer length(s) to process. Can be:
        - A single integer (e.g., 4) to process one k value
        - A tuple of two integers (k_min, k_max) to process an inclusive range of k values
        - A list of integers to process specific k values (e.g., [3, 5, 7])
    :type k: int | tuple[int, int] | list[int]
    :return: A nested dictionary mapping k value to a dictionary of k-mer sequences and their
        fractional abundances. The inner dictionaries are sorted by abundance in descending order.
        Format: {k: {kmer: fraction, ...}, ...}
    :rtype: dict[int, dict[str, float]]
    :raises ValueError: If k is invalid (e.g., negative values, empty range, or invalid tuple/list format).

    Examples:
        >>> # Count k-mers for k=4
        >>> kmer_fractions = count_kmer_abundance(files_fasta=["file1.fna", "file2.fna"], k=4)
        >>> # Returns: {4: {"ATCG": 0.001, "GCTA": 0.002, ...}}

        >>> # Count k-mers for range k=3 to k=5
        >>> kmer_fractions = count_kmer_abundance(files_fasta=["file1.fna"], k=(3, 5))
        >>> # Returns: {3: {...}, 4: {...}, 5: {...}}

        >>> # Count k-mers for specific k values
        >>> kmer_fractions = count_kmer_abundance(files_fasta=["file1.fna"], k=[3, 5, 7])
        >>> # Returns: {3: {...}, 5: {...}, 7: {...}}

        >>> # Filter top 10 most abundant for k=4
        >>> k4_fractions = kmer_fractions[4]
        >>> top_10 = dict(list(sorted(k4_fractions.items(), key=lambda x: x[1], reverse=True))[:10])

        >>> # Filter k-mers with fraction > 0.01 for k=4
        >>> high_abundance = {k: v for k, v in kmer_fractions[4].items() if v > 0.01}
    """
    # Normalize k input to a list of k values
    if isinstance(k, int):
        k_values = [k]
    elif isinstance(k, tuple):
        if len(k) != 2:
            raise ValueError(f"Tuple for k must have exactly 2 elements (k_min, k_max), got {len(k)}")
        k_min, k_max = k
        if k_min > k_max:
            raise ValueError(f"k_min ({k_min}) must be <= k_max ({k_max})")
        k_values = list(range(k_min, k_max + 1))
    elif isinstance(k, list):
        if len(k) == 0:
            raise ValueError("List of k values cannot be empty")
        k_values = k
    else:
        raise ValueError(f"k must be int, tuple[int, int], or list[int], got {type(k)}")

    if any(not isinstance(ki, int) or ki < 1 for ki in k_values):
        raise ValueError("All k values must be positive integers")

    # Normalize files_fasta input to a list
    files_fasta = cast_to_list(files_fasta)

    # Initialize FastaParser
    fasta_parser = FastaParser()

    # Initialize result dictionary
    result: dict[int, dict[str, float]] = {}

    # Process each k value
    for k_val in k_values:
        kmer_counts: Counter[str] = Counter()

        # Process each FASTA file
        for file_fasta in files_fasta:
            # Read sequences from FASTA file using FastaParser
            fasta_sequences = fasta_parser.read_fasta_sequences(file_fasta)
            for record in fasta_sequences:
                sequence = str(record.seq).upper()

                # Extract and count k-mers
                kmers = (sequence[i : i + k_val] for i in range(len(sequence) - k_val + 1))
                kmer_counts.update(kmers)

        # Calculate fractional abundances
        total_kmers = sum(kmer_counts.values())
        if total_kmers > 0:
            kmer_fractions = {kmer: count / total_kmers for kmer, count in kmer_counts.items()}
            # Sort by abundance (descending) for convenience
            result[k_val] = dict(sorted(kmer_fractions.items(), key=lambda x: x[1], reverse=True))
        else:
            # No k-mers found for this k value
            result[k_val] = {}

    return result
