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

    def __init__(self, database_name: str = "db_reference", dir_output: str = "output") -> None:
        """Constructor for the ReferenceDatabase class."""
        self.database_name = database_name
        self.dir_output = os.path.abspath(os.path.join(dir_output, database_name))
        Path(self.dir_output).mkdir(parents=True, exist_ok=True)

        self.fasta_parser = FastaParser()

        self.database = []

    def load_database_from_fasta(self, files_fasta: Union[str, List[str]], database_overwrite: bool) -> None:

        if database_overwrite:
            self.database = []

        files_fasta = check_if_list(files_fasta)
        for file in files_fasta:
            self.fasta_parser.check_fasta_format(file)
            fasta_sequences = self.fasta_parser.read_fasta_sequences(file)
            self.database.extend(fasta_sequences)

    def write_database_to_fasta(self, filename: str) -> str:

        file_database = os.path.join(self.dir_output, f"{filename}.fna")

        if self.database:
            with open(file_database, "w") as handle_fasta:
                SeqIO.write(self.database, handle_fasta, "fasta")
            self.file_fasta = file_database
        else:
            raise ValueError("Database is empty! Nothing to be written to fasta file.")

        return file_database

    def filter_database_by_region(self, region_ids: Union[str, List[str]], keep_region: bool) -> None:

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
