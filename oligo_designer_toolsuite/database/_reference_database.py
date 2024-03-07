############################################
# imports
############################################

import os
import warnings
from pathlib import Path
from typing import Union

import yaml
from Bio import SeqIO

from ..utils._checkers import check_if_list
from ..utils._sequence_parser import FastaParser

############################################
# Reference Database Class
############################################


class ReferenceDatabase:
    """The ReferenceDatabase class provides functionality for managing a reference sequence database.

    The header of each sequence must start with '>' and contain the following information:
    region_id, additional_information (optional) and coordinates (chrom, start, end, strand),
    where the region_id is compulsory and the other fileds are opional.

    Input Format (per sequence):
    >region_id::additional information::chromosome:start-end(strand)
    sequence

    Example:
    >ASR1::transcrip_id=XM456,exon_number=5::16:54552-54786(+)
    AGTTGACAGACCCCAGATTAAAGTGTGTCGCGCAACAC

    :param dir_output: The directory path for storing the reference database files. Defaults to "output".
    :type dir_output: str
    """

    def __init__(
        self,
        dir_output: str = "output",
    ):
        """Constructor for the ReferenceDatabase class."""
        self.fasta_parser = FastaParser()

        self.metadata = {}
        self.database = []

        self.dir_output = os.path.abspath(
            os.path.join(dir_output, "reference_database")
        )

    def load_metadata(self, metadata: Union[str, dict]):
        """Load metadata into the ReferenceDatabase object.

        If metadata already exists, a warning is issued about overwriting the existing metadata.
        The new metadata can be provided either as a path to a YAML file or directly as a dictionary.

        :param metadata: Path to a YAML file or a dictionary containing metadata information.
        :type metadata: Union[str, dict]

        :raises ValueError: If metadata has an incorrect format.
        """
        if self.metadata:
            warnings.warn(
                "Metadata not empty! Overwriting metadata with new metadata from file!"
            )

        if type(metadata) is str and os.path.exists(metadata):
            with open(metadata) as handle:
                self.metadata = yaml.safe_load(handle)
        elif type(metadata) is dict:
            self.metadata = metadata
        else:
            raise ValueError("Metadat has icorrect format!")

    def load_sequences_fom_fasta(
        self, file_fasta: str, database_overwrite: bool = False
    ) -> None:
        """Load sequences from a FASTA file into the ReferenceDatabase object.

        This function reads sequences from a FASTA file and adds them to the ReferenceDatabase object.
        If database_overwrite is True, the existing database is replaced with the new sequences;
        otherwise, the new sequences are appended to the existing database.

        :param file_fasta: Path to the FASTA file containing sequences.
        :type file_fasta: str
        :param database_overwrite: Flag indicating whether to overwrite the existing database, defaults to False.
        :type database_overwrite: bool, optional
        """
        self.fasta_parser.check_fasta_format(file_fasta)
        fasta_sequences = self.fasta_parser.read_fasta_sequences(file_fasta)
        if database_overwrite:
            warnings.warn("Overwriting database!")
            self.database = fasta_sequences
        else:
            self.database.extend(fasta_sequences)

    def write_database_to_fasta(self, filename: str) -> str:
        """Write sequences from the database to a FASTA file.

        This function writes the sequences from the database to a FASTA file. The file is saved in the specified
        directory with the given filename.

        :param filename: The name of the output FASTA file (without extension).
        :type filename: str
        :return: Path to the generated FASTA file.
        :rtype: str

        :raises ValueError: If the database is empty.
        """
        Path(self.dir_output).mkdir(parents=True, exist_ok=True)
        file_database = os.path.join(self.dir_output, f"{filename}.fna")
        if self.database:
            with open(file_database, "w") as handle_fasta:
                SeqIO.write(self.database, handle_fasta, "fasta")
            self.file_fasta = file_database
        else:
            raise ValueError("Database is empty! Nothing to be written to fasta file.")

        return file_database

    def write_metadata_to_yaml(self, filename: str) -> str:
        """Write metadata to a YAML file.

        This function writes the metadata of the OligoDatabase object to a YAML file.
        The file is saved in the specified directory with the given filename.

        :param filename: The name of the output YAML file (without extension).
        :type filename: str
        :return: Path to the generated YAML file.
        :rtype: str
        """
        Path(self.dir_output).mkdir(parents=True, exist_ok=True)
        file_metadata = os.path.join(self.dir_output, f"{filename}.yaml")

        with open(file_metadata, "w") as handle:
            yaml.safe_dump(
                self.metadata, handle, sort_keys=True, default_flow_style=False
            )

        return file_metadata

    def filter_database(self, region_ids: list[str] = None) -> None:
        """Filter the database entries based on specified region IDs.

        This function filters the database entries, keeping only those corresponding to the specified region IDs.
        If the database is empty, it raises a ValueError.

        :param region_ids: A list of region IDs to filter the database entries. If None, no filtering is applied.
        :type region_ids: list[str], optional

        :raises ValueError: If the database is empty.
        """
        region_ids = check_if_list(region_ids)

        if self.database:
            database_filtered = []
            for entry in self.database:
                region, _, _ = self.fasta_parser.parse_fasta_header(
                    entry.id, parse_additional_info=False
                )
                if region in region_ids:
                    database_filtered.append(entry)
            self.database = database_filtered
        else:
            raise ValueError(
                "Can not filter. Database is empty! Call the method load_fasta_into_database() first."
            )
