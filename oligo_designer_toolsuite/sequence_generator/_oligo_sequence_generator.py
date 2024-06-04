############################################
# imports
############################################

import os
import random
import subprocess

from pathlib import Path
from joblib import Parallel, delayed
from joblib_progress import joblib_progress

from oligo_designer_toolsuite.utils import FastaParser

from ..utils._checkers import check_if_list

############################################
# Oligo Database Class
############################################


class OligoSequenceGenerator:
    """A class for generating oligo sequences.

    This class provides functionality for generating oligo sequences and managing output directories.

    The generated sequences are saved as fasta file with region id, additional information and coordinates in header.
    The header of each sequence must start with '>' and contain the following information:
    region_id, additional_information (optional) and coordinates (chrom, start, end, strand),
    where the region_id is compulsory and the other fields are optional.

    Input Format (per sequence):
    >region_id::additional information::chromosome:start-end(strand)
    sequence

    Example:
    >ASR1::transcrip_id=XM456,exon_number=5::16:54552-54786(+)
    AGTTGACAGACCCCAGATTAAAGTGTGTCGCGCAACAC

    :param dir_output: The directory path for output files.
    :type dir_output: str
    """

    def __init__(
        self,
        dir_output: str = "output",
    ):
        """Constructor for the OligoSequenceGenerator class."""
        self.dir_output = os.path.abspath(os.path.join(dir_output, "annotation"))
        Path(self.dir_output).mkdir(parents=True, exist_ok=True)

        self.fasta_parser = FastaParser()

    def create_sequences_random(
        self,
        filename_out: str,
        length_sequences: int,
        num_sequences: int,
        name_sequences: str = "randomsequence",
        base_alphabet_with_probability: dict = {
            "A": 0.25,
            "C": 0.25,
            "G": 0.25,
            "T": 0.25,
        },
    ):
        """Create a FASTA file containing random DNA sequences.

        :param filename_out: The name of the output FASTA file.
        :type filename_out: str
        :param length_sequences: The length of each random DNA sequence.
        :type length_sequences: int
        :param num_sequences: The number of random DNA sequences to generate.
        :type num_sequences: int
        :param name_sequences: The prefix for sequence names in the output FASTA file.
        :type name_sequences: str, optional
        :param base_alphabet_with_probability: A dictionary mapping DNA bases to their probabilities.
        :type base_alphabet_with_probability: dict, optional
        :return: The path to the generated FASTA file.
        :rtype: str
        """

        def get_sequence_random(sequence_length, base_alphabet_with_probability):
            """Generate a random DNA sequence of a given length based on a specified base alphabet with probabilities.

            :param sequence_length: The length of the generated DNA sequence.
            :type sequence_length: int
            :param base_alphabet_with_probability: A dictionary mapping DNA bases to their probabilities.
            :type base_alphabet_with_probability: dict
            :return: The generated random DNA sequence.
            :rtype: str
            """
            bases = list(base_alphabet_with_probability.keys())
            sequence = "".join(
                random.choices(
                    bases,
                    weights=[base_alphabet_with_probability[n] for n in bases],
                    k=sequence_length,
                )
            )
            return sequence

        sequences_set = set()
        while len(sequences_set) < num_sequences:
            num_missing_sequences = num_sequences - len(sequences_set)
            new_sequences = {
                get_sequence_random(length_sequences, base_alphabet_with_probability)
                for _ in range(num_missing_sequences)
            }
            sequences_set.update(new_sequences)

        file_fasta_out = os.path.join(self.dir_output, f"{filename_out}.fna")

        with open(file_fasta_out, "w") as handle_fasta:
            for seq in sequences_set:
                handle_fasta.write(f">{name_sequences}::regiontype=random_sequence\n{seq}\n")
        return file_fasta_out

    def create_sequences_sliding_window(
        self,
        filename_out: str,
        files_fasta_in: list[str],
        length_interval_sequences: tuple,
        region_ids: list[str] = None,
        n_jobs: int = 1,
    ):
        """Generate sequences with sliding windows from input FASTA files and write them to an output FASTA file.

        :param filename_out: The name of the output FASTA file.
        :type filename_out: str
        :param files_fasta_in: List of input FASTA files.
        :type files_fasta_in: list[str]
        :param length_interval_sequences: A tuple representing the interval of sliding window lengths to generate.
        :type length_interval_sequences: tuple
        :param region_ids: List of region IDs to consider. If None, all regions are considered.
        :type region_ids: list[str], optional
        :param n_jobs: The number of jobs to use for parallel processing.
        :type n_jobs: int
        :return: The path to the generated output FASTA file.
        :rtype: str
        """

        def generate_unique_filename(region):
            while True:
                random_number = random.randint(1, 1000000)
                file_fasta_region = os.path.join(self.dir_output, f"{region}_{random_number}.fna")
                if not os.path.exists(file_fasta_region):
                    return file_fasta_region

        def get_sliding_window_sequence(entry, length_interval_sequences):
            """Extract sliding window sequences from a given entry and write them to a FASTA file.

            :param entry: Bio.SeqRecord.SeqRecord object representing the sequence entry.
            :type entry: Bio.SeqRecord.SeqRecord
            :param length_interval_sequences: The max and min length of each sliding window sequence.
            :type length_interval_sequences: List[int]
            :param handle_fasta: The file handle for writing the output FASTA file.
            :type handle_fasta: _io.TextIOWrapper
            """

            entry_sequence = entry.seq
            region, additional_info, coordinates = self.fasta_parser.parse_fasta_header(
                header=entry.id, parse_additional_info=False
            )
            # chromosome and strand information is the same for an entry but parsed as
            # list in the parse_fasta_header, hence, we take only the first element for each
            # if no information provided in fasta file, those entries are None
            chromosome = coordinates["chromosome"][0]
            strand = coordinates["strand"][0]
            list_of_coordinates = []
            # if the fasta header DOES NOT contain coordinates information
            if chromosome is None:
                list_of_coordinates = [None for i in range(len(entry_sequence))]
            # if the fasta header DOES contain coordinates information
            else:
                # coordinates in fasta file use 1-base indixing, which go from 1 (for base 1) to n (for base n)
                # range produces values until end-1 (e.g. range(10) goes until 9) -> add +1
                # the start and end coordinates can be a list for regions spanning split sequences (e.g. exon junctions, UTRs)
                for start, end in zip(coordinates["start"], coordinates["end"]):
                    list_of_coordinates.extend(range(start, end + 1))
            # sort reverse on minus strand becaus the sequence is translated into the reverse complement by fasta -strand option
            if strand == "-":
                list_of_coordinates.reverse()

            file_fasta_region = generate_unique_filename(region=region)
            with open(file_fasta_region, "w") as handle_fasta:
                for sequence_length in range(length_interval_sequences[0], length_interval_sequences[1] + 1):
                    # generate sequences with sliding window and write to fasta file (use lock to ensure that hat only one process can write to the file at any given time)
                    if len(entry_sequence) < sequence_length:
                        continue
                    num_sequences = len(entry_sequence) - (sequence_length - 1)
                    sequences = [entry_sequence[i : i + sequence_length] for i in range(num_sequences)]

                    for i in range(num_sequences):
                        seq = sequences[i]
                        seq_start_end = [
                            list_of_coordinates[i],
                            list_of_coordinates[(i + sequence_length - 1)],
                        ]
                        # sort coordinates for oligos on minus strand
                        start_seq = min(seq_start_end)  # 1-base index
                        end_seq = max(seq_start_end)

                        header = f"{region}::{additional_info}::{chromosome}:{start_seq}-{end_seq}({strand})"
                        handle_fasta.write(f">{header}\n{seq}\n")
            return file_fasta_region

        files_fasta_in = check_if_list(files_fasta_in)
        for file_fasta in files_fasta_in:
            self.fasta_parser.check_fasta_format(file_fasta)

        if region_ids:
            region_ids = check_if_list(region_ids)
        else:
            region_ids = [
                self.fasta_parser.get_fasta_regions(file_fasta_in=file_fasta) for file_fasta in files_fasta_in
            ]
            # make keys unique
            region_ids = list(set(region_ids))

        # delete previous content
        for region_id in region_ids:
            file_fasta_region = os.path.join(self.dir_output, f"{region_id}.fna")
            if os.path.isfile(file_fasta_region):
                os.remove(file_fasta_region)

        # create oligos and store all oligos of one region in seperate files
        file_fasta_out = set()
        for file_fasta in files_fasta_in:
            self.fasta_parser.check_fasta_format(file_fasta)
            fasta_sequences = self.fasta_parser.read_fasta_sequences(file_fasta, region_ids)
            with joblib_progress(
                description=f"Oligo Generation from {os.path.basename(file_fasta)}",
                total=len(fasta_sequences),
            ):
                files_fasta_oligos = Parallel(n_jobs=n_jobs)(
                    delayed(get_sliding_window_sequence)(entry, length_interval_sequences)
                    for entry in fasta_sequences
                )
            for region_id in region_ids:
                files_fasta_oligos_region = [
                    file for file in files_fasta_oligos if os.path.basename(file).startswith(f"{region_id}_")
                ]
                if len(files_fasta_oligos_region) > 0:
                    file_fasta_region = os.path.join(self.dir_output, f"{region_id}.fna")
                    file_fasta_out.add(file_fasta_region)
                    # merge the fasta files into one single file per region
                    self.fasta_parser.merge_fasta_files(
                        files_in=files_fasta_oligos_region, file_out=file_fasta_region, overwrite=False
                    )

            # remove the temporary fasta files
            for file_fasta_oligos in files_fasta_oligos:
                if os.path.isfile(file_fasta_oligos):
                    os.remove(file_fasta_oligos)

        return sorted(list(file_fasta_out))
