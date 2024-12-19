############################################
# imports
############################################

import os
from pathlib import Path
from typing import Union, List

from Bio import SeqIO

from oligo_designer_toolsuite.utils import FastaParser, check_if_list

############################################
# Reference Database Class
############################################


class ReferenceDatabase:
    """
    The `ReferenceDatabase` class manages a reference sequence database used for oligonucleotide design.
    It handles the initialization, storage, and management of sequence data in a specified output directory.

    :param database_name: The name of the ReferenceDatabase, defaults to "db_reference".
    :type database_name: str
    :param dir_output: The directory where the database will be stored, defaults to "output".
    :type dir_output: str
    """

    def __init__(self, database_name: str = "db_reference", dir_output: str = "output") -> None:
        """Constructor for the ReferenceDatabase class."""
        self.database_name = database_name
        self.dir_output = os.path.abspath(os.path.join(dir_output, database_name))
        Path(self.dir_output).mkdir(parents=True, exist_ok=True)

        self.fasta_parser = FastaParser()

        self.database = []

    def load_database_from_fasta(self, files_fasta: Union[str, List[str]], database_overwrite: bool) -> None:
        """
        Loads a list of sequences from one or more FASTA file(s) into the database.
        If the `database_overwrite` flag is set to True, the existing database will be cleared before loading the new sequences.

        The header of each sequence must start with '>' and contain the following information:
        region_id, additional_information (optional) and coordinates (chrom, start, end, strand),
        where the region_id is compulsory and the other fileds are opional.

        Input Format (per sequence):
        >region_id::additional information::chromosome:start-end(strand)
        sequence

        Example:
        >ASR1::transcrip_id=XM456,exon_number=5::16:54552-54786(+)
        AGTTGACAGACCCCAGATTAAAGTGTGTCGCGCAACAC


        :param files_fasta: A single FASTA file or a list of FASTA files to be loaded into the database.
        :type files_fasta: Union[str, List[str]]
        :param database_overwrite: If True, the existing database content will be cleared before loading the new sequences.
        :type database_overwrite: bool
        """
        if database_overwrite:
            self.database = []

        files_fasta = check_if_list(files_fasta)
        for file in files_fasta:
            self.fasta_parser.check_fasta_format(file)
            fasta_sequences = self.fasta_parser.read_fasta_sequences(file)
            self.database.extend(fasta_sequences)

    def write_database_to_fasta(self, filename: str) -> str:
        """
        Writes the current database sequences to a FASTA file with the specified filename.
        If the database is empty, an error is raised.

        :param filename: The name of the output FASTA file (without extension) to save the database sequences.
        :type filename: str
        :return: The path to the written FASTA file.
        :rtype: str
        """
        file_database = os.path.join(self.dir_output, f"{filename}.fna")

        if self.database:
            with open(file_database, "w") as handle_fasta:
                SeqIO.write(self.database, handle_fasta, "fasta")
            self.file_fasta = file_database
        else:
            raise ValueError("Database is empty! Nothing to be written to fasta file.")

        return file_database

    def filter_database_by_region(self, region_ids: Union[str, List[str]], keep_region: bool) -> None:
        """
        Filters the database to include or exclude regions based on the provided list of region IDs.
        The regions to keep or remove are determined by the `keep_region` flag.

        :param region_ids: A single region ID or a list of region IDs to filter by.
        :type region_ids: Union[str, List[str]]
        :param keep_region: If True, only the specified regions will be kept in the database. If False, the specified regions will be removed from the database.
        :type keep_region: bool
        """
        region_ids = check_if_list(region_ids)

        if self.database:
            database_filtered = []
            for entry in self.database:
                region, _, _ = self.fasta_parser.parse_fasta_header(entry.id, parse_additional_info=False)
                if keep_region and (region in region_ids):
                    database_filtered.append(entry)
                elif not keep_region and (region not in region_ids):
                    database_filtered.append(entry)
            self.database = database_filtered
        else:
            raise ValueError(
                "Can not filter. Database is empty! Call the method load_database_from_fasta() first."
            )

    def filter_database_by_attribute_category(
        self, attribute_name: str, attribute_category: Union[str, List[str]], keep_if_equals_category: bool
    ) -> None:
        """
        Filters the database based on a specific attribute's category.
        Sequences are either kept or removed based on whether their attribute values match the specified category.

        :param attribute_name: The name of the attribute to filter by.
        :type attribute_name: str
        :param attribute_category: A single category or a list of categories to filter by.
        :type attribute_category: Union[str, List[str]]
        :param keep_if_equals_category: If True, sequences with matching categories are kept. If False, sequences with matching categories are removed.
        :type keep_if_equals_category: bool
        """
        attribute_category = check_if_list(attribute_category)

        if self.database:
            database_filtered = []
            for entry in self.database:
                region_id, attributes, _ = self.fasta_parser.parse_fasta_header(
                    entry.id, parse_additional_info=True
                )
                if attribute_name in attributes:
                    attribute_values = check_if_list(attributes[attribute_name])
                    if keep_if_equals_category and any(
                        item in attribute_category for item in attribute_values
                    ):
                        database_filtered.append(entry)
                    elif not keep_if_equals_category and all(
                        item not in attribute_category for item in attribute_values
                    ):
                        database_filtered.append(entry)
            self.database = database_filtered
        else:
            raise ValueError(
                "Can not filter. Database is empty! Call the method load_fasta_into_database() first."
            )
