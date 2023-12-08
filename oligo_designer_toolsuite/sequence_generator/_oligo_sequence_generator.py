############################################
# imports
############################################

import os
import random

from pathlib import Path

from ..utils._utils import check_if_list
from ..utils._sequence_parser import FastaParser


############################################
# Oligo Database Class
############################################


class OligoSequenceGenerator:
    """This class generates all possible oligos that can be designed for a given list of regions (e.g. genes),
    based on the transcriptome or the gene CDS annotation or the whole genome provided as fasta file.

    The header of each sequence must start with '>' and contain the following information:
    region_id, additional_information (optional) and coordinates (chrom, start, end, strand),
    where the region_id is compulsory and the other fileds are opional.

    Input Format (per sequence):
    >region_id::additional information::chromosome:start-end(strand)
    sequence

    Example:
    >ASR1::transcrip_id=XM456,exon_number=5::16:54552-54786(+)
    AGTTGACAGACCCCAGATTAAAGTGTGTCGCGCAACAC

    Moreover, the database can be saved and loaded to/from a tsv file.

    :param min_oligos_per_region: Minimum number of oligos per region, if lower, region is removed from database, defaults to 0.
    :type min_oligos_per_region: int, optional
    :param write_regions_with_insufficient_oligos: write removed regions to file ``regions_with_insufficient_oligos.txt``, defaults to True
    :type write_regions_with_insufficient_oligos: bool, optional
    :param n_jobs: Number of parallel processes, if set to None all available CPUs are use, defaults to None.
    :type n_jobs: int, optional
    :param dir_output: Output directory, defaults to 'output'.
    :type dir_output: str, optional
    """

    def __init__(
        self,
        dir_output: str = "output",
    ):
        """Constructor"""
        self.dir_output = os.path.abspath(os.path.join(dir_output, "annotation"))
        Path(self.dir_output).mkdir(parents=True, exist_ok=True)

        self.fasta_parser = FastaParser()

    def create_sequences_random(
        self,
        filename_out: str,
        length_sequences: int,
        num_sequences: int,
        name_sequences: str = "randomsequence",
        base_alphabet_with_probability: list[float] = {"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25},
    ):
        """create sequences from random sequences with desired length and sequence composition.
        save as fasta file, with appropriate header, e.g. region id and type=random_sequence."""

        def get_sequence_random(sequence_length, base_alphabet_with_probability):
            bases = list(base_alphabet_with_probability.keys())
            sequence = "".join(
                random.choices(
                    bases, weights=[base_alphabet_with_probability[n] for n in bases], k=sequence_length
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
    ):
        """create new fasta files with tiled sequences from input sequences via sliding window.
        Load sequences from input fasta file, and run a sliding window over them to create
        tiled sequences with desired length.
        save the generated sequences in fasta file with header information taken from input fasta
        file but adjust coordinates."""

        def get_sliding_window_sequence(entry, sequence_length, handle_fasta):
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

            # generate sequences with sliding window and write to fasta file (use lock to ensure that hat only one process can write to the file at any given time)
            if len(entry_sequence) > sequence_length:
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

        region_ids = check_if_list(region_ids)
        file_fasta_in = check_if_list(file_fasta_in)

        file_fasta_out = os.path.join(self.dir_output, f"{filename_out}.fna")
        with open(file_fasta_out, "w") as handle_fasta:
            for file_fasta in file_fasta_in:
                fasta_sequences = self.fasta_parser.read_fasta_sequences(file_fasta, region_ids)
                # did not parallize this function because workers would be writing simultaneously to the same file
                # which could lead to data corruption. I tried using the file_lock=Lock() function from multiprocessing
                # as input to get_sliding_window_sequence(file_lock) and then use with file_lock: before writing to the file
                # but this lead to an pickle error from joblib which I could not solve, so I left it like that for now
                for entry in fasta_sequences:
                    for length_sequences in range(
                        length_interval_sequences[0], length_interval_sequences[1] + 1
                    ):
                        get_sliding_window_sequence(entry, length_sequences, handle_fasta)

        return file_fasta_out
