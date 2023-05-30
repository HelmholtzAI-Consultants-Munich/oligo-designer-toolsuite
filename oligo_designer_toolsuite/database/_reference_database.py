############################################
# imports
############################################

import os
import yaml
import warnings
from Bio import SeqIO
from pathlib import Path

from ..utils._data_parser import check_fasta_format, parse_fasta_header

############################################
# Reference Database Class
############################################


class ReferenceDatabase:
    """This class stores a refernce database, which is created from
    a sequence fasta file.
    The header of each sequence must start with '>' and contain the following information:
    region_id, additional_information (optional) and coordinates (chrom, start, end, strand),
    where the region_id is compulsory and the other fileds are opional.


    Input Format (per sequence):
    >region_id::additional information::chromosome:start-end(strand)
    sequence

    Example:
    >ASR1::transcrip_id=XM456,exon_number=5::16:54552-54786(+)
    AGTTGACAGACCCCAGATTAAAGTGTGTCGCGCAACAC

    :param file_fasta: Path to the fasta file.
    :type file_fasta: str
    :param metadata: database metadata like species, annotation release, genome assembly, ect.
    :type metadata: dict, optional
    :param dir_output: Output directory, defaults to 'output'.
    :type dir_output: str, optional
    """

    def __init__(
        self,
        file_fasta: str,
        metadata: dict = {},
        dir_output: str = "output",
    ):
        """Constructor"""
        self.metadata = metadata

        self.file_fasta = file_fasta
        if os.path.exists(self.file_fasta):
            if not check_fasta_format(self.file_fasta):
                raise ValueError("Fasta file has incorrect format!")
        else:
            raise ValueError("Fasta file does not exist!")

        self.database = []

        self.dir_output = os.path.abspath(
            os.path.join(dir_output, "reference_database")
        )

    def load_fasta_into_database(self):
        """Load sequences from fasta file and stored in databse."""
        with open(self.file_fasta, "r") as handle:
            self.database = list(SeqIO.parse(handle, "fasta"))

    def write_fasta_from_database(
        self,
        filename: str,
    ):
        """Write sequences stored in database to fasta file.

        :param filename: Fasta file name prefix.
        :type filename: str
        :return: Path to fasta file.
        :rtype: str
        """
        Path(self.dir_output).mkdir(parents=True, exist_ok=True)
        file_database = os.path.join(self.dir_output, f"{filename}.fna")
        if self.database:
            with open(file_database, "w") as handle_fasta:
                SeqIO.write(self.database, handle_fasta, "fasta")
            self.file_fasta = file_database
        else:
            raise ValueError("Database is empty! Nothing written to fasta file.")

        return file_database

    def write_metadata_from_database(
        self,
        filename: str,
    ):
        Path(self.dir_output).mkdir(parents=True, exist_ok=True)
        file_metadata = os.path.join(self.dir_output, f"{filename}.yaml")

        with open(file_metadata, "w") as handle:
            yaml.safe_dump(
                self.metadata, handle, sort_keys=True, default_flow_style=False
            )

        return file_metadata

    def filter_database(
        self,
        region_ids: list[str] = None,
    ):
        """Filter database based on region id and overwrite database with filtered database.

        :param region_ids: List of region ids.
        :type region_ids: list
        """
        if self.database:
            database_filtered = []
            for entry in self.database:
                region, _, _ = parse_fasta_header(entry.id)
                if region in region_ids:
                    database_filtered.append(entry)
            self.database = database_filtered
        else:
            raise ValueError(
                "Can not filter. Database is empty! Call the method load_fasta_into_database() first."
            )
