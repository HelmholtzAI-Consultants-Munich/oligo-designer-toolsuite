############################################
# imports
############################################

import os
import warnings
from pathlib import Path

from Bio import SeqIO

from oligo_designer_toolsuite.utils import FastaParser, check_if_list

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

    :param database_name: Subdirectory path for the output, i.e. <dir_output>/<database_name>, defaults to "db_reference".
    :type database_name: str, optional
    :param dir_output: Directory path for the output, defaults to "output".
    :type dir_output: str, optional
    """

    def __init__(self, database_name: str = "db_reference", dir_output: str = "output"):
        """Constructor for the ReferenceDatabase class."""
        self.database_name = database_name
        self.dir_output = os.path.abspath(os.path.join(dir_output, database_name))
        Path(self.dir_output).mkdir(parents=True, exist_ok=True)

        self.fasta_parser = FastaParser()

        self.database = []

    def load_database_from_fasta(self, files_fasta: list[str], database_overwrite: bool = False) -> None:
        """Load sequences from one or more FASTA files into the ReferenceDatabase object.

        This function reads sequences from FASTA file(s) and adds them to the ReferenceDatabase object.
        If 'database_overwrite' is True, it clears the existing database before loading new sequences;
        otherwise, the new sequences are appended to the existing database.

        :param files_fasta: Paths to the FASTA file(s) containing sequences.
        :type files_fasta: list[str]
        :param database_overwrite: Flag indicating whether to overwrite the existing database, defaults to False.
        :type database_overwrite: bool, optional
        """
        if database_overwrite:
            warnings.warn("Overwriting database!")
            self.database = []

        files_fasta = check_if_list(files_fasta)
        for file in files_fasta:
            self.fasta_parser.check_fasta_format(file)
            fasta_sequences = self.fasta_parser.read_fasta_sequences(file)
            self.database.extend(fasta_sequences)

    def write_database_to_fasta(self, filename: str) -> str:
        """Write sequences from the database to a FASTA file.

        This function writes the sequences from the database to a FASTA file. The file is saved in the specified
        directory with the given filename.

        :param filename: The name of the output FASTA file (without extension).
        :type filename: str
        :param dir_output: Directory path for the output.
        :type dir_output: str
        :return: Path to the generated FASTA file.
        :rtype: str

        :raises ValueError: If the database is empty.
        """
        file_database = os.path.join(self.dir_output, f"{filename}.fna")

        if self.database:
            with open(file_database, "w") as handle_fasta:
                SeqIO.write(self.database, handle_fasta, "fasta")
            self.file_fasta = file_database
        else:
            raise ValueError("Database is empty! Nothing to be written to fasta file.")

        return file_database

    def filter_database(self, region_ids: list[str] = None, remove_region: bool = True) -> None:
        """Filters the database entries based on specified region IDs, either keeping or removing entries from those regions.

        This method modifies the database in-place, selectively keeping or excluding entries based on their region ID. It's
        essential to have the database loaded before invoking this function; otherwise, a ValueError is raised to indicate an
        empty database.

        :param region_ids: A list of region IDs to filter the database entries by. If None, no filtering is applied.
        :type region_ids: list[str], optional
        :param remove_region: A flag indicating whether to remove (True) or keep (False) the entries from the specified regions.
        :type remove_region: bool
        :raises ValueError: If the database is empty.
        """
        region_ids = check_if_list(region_ids)

        if self.database:
            database_filtered = []
            for entry in self.database:
                region, _, _ = self.fasta_parser.parse_fasta_header(entry.id, parse_additional_info=False)
                if remove_region:
                    if region not in region_ids:
                        database_filtered.append(entry)
                else:
                    if region in region_ids:
                        database_filtered.append(entry)
            self.database = database_filtered
        else:
            raise ValueError(
                "Can not filter. Database is empty! Call the method load_fasta_into_database() first."
            )
