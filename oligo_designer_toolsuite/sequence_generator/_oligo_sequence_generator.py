############################################
# imports
############################################

import os
import random
from pathlib import Path

from joblib import Parallel, delayed

from ..utils._checkers import check_if_list
from ..utils._sequence_parser import FastaParser, merge_fasta_files

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
        base_alphabet_with_probability: list[float] = {
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

        # make sure to have n unique sequences
        sequences_list = []
        while len(sequences_list) < num_sequences:
            num_missing_sequences = num_sequences - len(sequences_list)
            missing_sequences = list(
                set(
                    [
                        get_sequence_random(length_sequences, base_alphabet_with_probability)
                        for i in range(num_missing_sequences)
                    ]
                )
            )
            sequences_list = list(set(sequences_list + missing_sequences))

        file_fasta_out = os.path.join(self.dir_output, f"{filename_out}.fna")

        with open(file_fasta_out, "w") as handle_fasta:
            for seq in sequences_list:
                handle_fasta.write(f">{name_sequences}::regiontype=random_sequence\n{seq}\n")
        return file_fasta_out

    def create_sequences_sliding_window(
        self,
        filename_out: str,
        file_fasta_in: list[str],
        length_interval_sequences: tuple,
        region_ids: list[str] = None,
        n_jobs: int = 1
    ):
        """Generate sequences with sliding windows from input FASTA files and write them to an output FASTA file.

        :param filename_out: The name of the output FASTA file.
        :type filename_out: str
        :param file_fasta_in: List of input FASTA files.
        :type file_fasta_in: list[str]
        :param length_interval_sequences: A tuple representing the interval of sliding window lengths to generate.
        :type length_interval_sequences: tuple
        :param region_ids: List of region IDs to consider. If None, all regions are considered.
        :type region_ids: list[str], optional
        :param n_jobs: The number of jobs to use for parallel processing. 
        :type n_jobs: int
        :return: The path to the generated output FASTA file.
        :rtype: str
        """

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
                list_of_coordinates = [None for i in range(len(seq))]
            # if the fasta header DOES contain coordinates information
            else:
                # coordinates in fasta file use 1-base indixing, which go from 1 (for base 1) to n (for base n)
                # range produces values until end-1 (e.g. range(10) goes until 9) -> add +1
                # the start and end coordinates can be a list for regions spanning split sequences (e.g. exon junctions, UTRs)
                for i in range(len(coordinates["start"])):
                    list_of_coordinates.extend(
                        list(
                            range(
                                coordinates["start"][i],
                                coordinates["end"][i] + 1,
                            )
                        )
                    )

            # sort reverse on minus strand becaus the sequence is translated into the reverse complement by fasta -strand option
            if strand == "-":
                list_of_coordinates.reverse()
            file_fasta_region = os.path.join(self.dir_output, f"{region}_{random.randint(0, 1e5)}.fna")
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

        region_ids = check_if_list(region_ids)
        file_fasta_in = check_if_list(file_fasta_in)

        # file_fasta_out = os.path.join(self.dir_output, f"{filename_out}.fna")
        # with open(file_fasta_out, "w") as handle_fasta:
        #     for file_fasta in file_fasta_in:
        #         self.fasta_parser.check_fasta_format(file_fasta)
        #         fasta_sequences = self.fasta_parser.read_fasta_sequences(file_fasta, region_ids)
        #         # did not parallize this function because workers would be writing simultaneously to the same file
        #         # which could lead to data corruption. I tried using the file_lock=Lock() function from multiprocessing
        #         # as input to get_sliding_window_sequence(file_lock) and then use with file_lock: before writing to the file
        #         # but this lead to an pickle error from joblib which I could not solve, so I left it like that for now
        #         for entry in fasta_sequences:
        #             for length_sequences in range(
        #                 length_interval_sequences[0], length_interval_sequences[1] + 1
        #             ):
        #                 get_sliding_window_sequence(entry, length_sequences, handle_fasta)
        file_fasta_out = os.path.join(self.dir_output, f"{filename_out}.fna")
        for file_fasta in file_fasta_in:
            self.fasta_parser.check_fasta_format(file_fasta)
            fasta_sequences = self.fasta_parser.read_fasta_sequences(file_fasta, region_ids)
            files_fasta_region = Parallel(n_jobs=n_jobs)(
            delayed(get_sliding_window_sequence)(entry, length_interval_sequences) for entry in fasta_sequences)
            # merge the fasta files into one single file
            merge_fasta_files(files_in=files_fasta_region, file_out=file_fasta_out, overwrite=False)
            # remove the temporary fasta files
            for file_fasta_region in files_fasta_region:
                if os.path.isfile(file_fasta_region):
                    os.remove(file_fasta_region)
            
        return file_fasta_out
