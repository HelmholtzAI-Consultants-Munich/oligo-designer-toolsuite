############################################
# imports
############################################

import os
import shutil
from pathlib import Path
from typing import Any, get_args

from oligo_designer_toolsuite._constants import _TYPES_REF
from oligo_designer_toolsuite._exceptions import DatabaseError
from oligo_designer_toolsuite.utils import FastaParser, VCFParser, cast_to_list, remove_index_files

############################################
# Reference Database Class
############################################


class ReferenceDatabase:
    """
    The `ReferenceDatabase` class manages a reference sequence database used for oligonucleotide design.
    It handles the initialization, storage, and management of sequence data in a specified output directory.

    :param database_name: The name of the ReferenceDatabase, defaults to "db_reference".
    :type database_name: str
    :param dir_output: Directory path where output files will be saved. Defaults to "output".
    :type dir_output: str
    """

    def __init__(self, database_name: str = "db_reference", dir_output: str = "output") -> None:
        """Constructor for the ReferenceDatabase class."""
        self.database_name = database_name
        self.dir_output = os.path.abspath(dir_output)

        self.database_file: str | None = None
        self.database_type: str | None = None

        self.fasta_parser = FastaParser()
        self.vcf_parser = VCFParser()

    def __del__(self) -> None:
        """Destructor for the ReferenceDatabase class."""
        if self.database_file is not None and os.path.exists(self.database_file):
            os.remove(self.database_file)

    def load_database_from_file(
        self,
        files: str | list[str],
        file_type: _TYPES_REF,
        database_overwrite: bool,
    ) -> None:
        """
        Load a database from one or more files and set the database type.
        This function only lazy loads the database and stores the file path.
        The file content is only loaded when the database is manipulated with any of the filter functions.
        If the `database_overwrite` flag is set to True, the existing database will be cleared before loading.

        If the database is laoded form a fasta file, the header of each sequence must start with '>' and contain the following information:
        region_id, additional_information (optional) and coordinates (chrom, start, end, strand) (optional),
        where the region_id is compulsory and the other fileds are opional.

        Input Format (per sequence):
        >region_id::additional information::chromosome:start-end(strand)
        sequence

        Example:
        >ASR1::transcrip_id=XM456,exon_number=5::16:54552-54786(+)
        AGTTGACAGACCCCAGATTAAAGTGTGTCGCGCAACAC

        :param files: Path(s) to the file(s) containing the database sequences.
        :type files: str | list[str]
        :param file_type: Type of the reference sequences (must be a valid type).
        :type file_type: _TYPES_REF["fasta", "vcf"]
        :param database_overwrite: If True, the existing database content will be cleared before loading the new sequences.
        :type database_overwrite: bool
        """
        # Check if output folder exists
        Path(self.dir_output).mkdir(parents=True, exist_ok=True)

        # Check if file type is correct
        options = get_args(_TYPES_REF)
        assert file_type in options, f"Sequence type not supported! '{file_type}' is not in {options}."
        files = cast_to_list(files)

        # remove all files if database should be overwritten
        if self.database_file is not None and database_overwrite:
            # Remove the main database file
            remove_index_files(file_reference=self.database_file, dir_output=self.dir_output)
            if os.path.exists(self.database_file):
                os.remove(self.database_file)

            self.database_file = None
            self.database_type = None

        # set retrieve input file list
        if self.database_file is None:
            files_in = files
            self.database_type = file_type
        elif self.database_file is not None and self.database_type == file_type:
            files_in = files + [self.database_file]
        else:
            raise DatabaseError(
                f"Cannot mix {file_type} and {self.database_type} databases. "
                f"Use database_overwrite=True to replace the existing database."
            )

        if self.database_type == "fasta":
            for file in files_in:
                self.fasta_parser.check_fasta_format(file)
            self.database_file = os.path.join(self.dir_output, f"tmp_{self.database_name}.fna")
            self.fasta_parser.merge_fasta_files(
                files_in=files_in, file_out=self.database_file, overwrite=True
            )
        elif self.database_type == "vcf":
            self.database_type = file_type
            self.database_file = os.path.join(self.dir_output, f"tmp_{self.database_name}.vcf.gz")
            self.vcf_parser.merge_vcf_files(files_in=files_in, file_out=self.database_file)
        else:
            raise DatabaseError(f"Database type {self.database_type} not supported.")

    def write_database_to_file(self, filename: str, dir_output: str | None = None) -> str:
        """
        Write the loaded database to a file based on its type.

        :param filename: Name of the output file (without extension).
        :type filename: str
        :param dir_output: Directory path where output files will be saved. If None, the directory defined in the init function will be used.
        :type dir_output: str | None
        :return: Path to the written database file.
        :rtype: str
        :raises ValueError: If the database type is not supported or if the database is empty.
        """
        if dir_output:
            Path(dir_output).mkdir(parents=True, exist_ok=True)
        else:
            dir_output = self.dir_output

        if self.database_file:
            file_ending = "fna" if self.database_type == "fasta" else "vcf.gz"
            file_database = os.path.join(dir_output, f"{filename}.{file_ending}")
            shutil.copy2(self.database_file, file_database)
            if self.database_type == "fasta":
                self.fasta_parser.index_fasta_file(file_fasta=file_database)
            return file_database
        else:
            raise DatabaseError("Database is empty. Nothing to be written to file.")

    def filter_database_by_region(self, region_ids: str | list[str], keep_region: bool) -> None:
        """
        Filter the database to retain or remove specific regions.

        :param region_ids: Region identifier(s) to process. Can be a single region ID (str) or a list of region IDs (list[str]). If None, all regions in the database are processed.
        :type region_ids: str | list[str]
        :param keep_region: Whether to keep (True) or remove (False) the specified regions.
        :type keep_region: bool
        :return: Path to the filtered database file.
        :rtype: str
        :raises ValueError: If the database is empty or filtering is attempted on a non-FASTA database.
        """
        region_ids = cast_to_list(region_ids)

        if self.database_file:
            if self.database_type == "fasta":
                file_database_filtered = self._filter_fasta_database_by_region(
                    region_ids=region_ids,
                    keep_region=keep_region,
                )
                os.remove(self.database_file)
                self.database_file = file_database_filtered
            else:
                raise DatabaseError("Filter only available for database in fasta format.")
        else:
            raise DatabaseError(
                "Cannot filter database: database is empty. Call load_database_from_file() first."
            )

    def filter_database_by_property_category(
        self,
        property_name: str,
        property_category: str | list[str],
        keep_if_equals_category: bool,
    ) -> Any:
        """
        Filter the database to retain or remove specific property categories.

        :param property_name: The property used for filtering.
        :type property_name: str
        :param property_category: list of property values to filter by.
        :type property_category: str | list[str]
        :param keep_if_equals_category: Whether to keep (True) or remove (False) records with the specified property values.
        :type keep_if_equals_category: bool
        :return: Path to the filtered database file.
        :rtype: str
        :raises ValueError: If the database is empty or filtering is attempted on a non-FASTA database.
        """
        property_category = cast_to_list(property_category)

        if self.database_file:
            if self.database_type == "fasta":
                file_database_filtered = self._filter_fasta_database_by_property_category(
                    property_name=property_name,
                    property_category=property_category,
                    keep_if_equals_category=keep_if_equals_category,
                )
                os.remove(self.database_file)
                self.database_file = file_database_filtered
            else:
                raise DatabaseError("Filter only available for database in fasta format.")
        else:
            raise DatabaseError(
                "Cannot filter database: database is empty. Call load_database_from_file() first."
            )

        return file_database_filtered

    def _filter_fasta_database_by_region(
        self,
        region_ids: str | list[str],
        keep_region: bool,
    ) -> str:
        """
        Filter a FASTA database to retain or exclude specific regions.
        Therefore, merge fasta files and load fasta content for filtering.

        :param region_ids: Region identifier(s) to process. Can be a single region ID (str) or a list of region IDs (list[str]). If None, all regions in the database are processed.
        :type region_ids: str | list[str]
        :param keep_region: Whether to keep (True) or remove (False) the specified regions.
        :type keep_region: bool
        """
        if self.database_file:
            file: str = self.database_file
        else:
            raise DatabaseError("Database file is not set. Call load_database_from_file() first.")

        if keep_region:
            regions_to_keep = region_ids
        else:
            regions_fasta = self.fasta_parser.get_fasta_regions(file_fasta_in=file)
            regions_to_keep = [region for region in regions_fasta if region not in region_ids]

        fasta_sequences_filtered = self.fasta_parser.read_fasta_sequences(
            file_fasta_in=file, region_ids=regions_to_keep
        )

        file_database_filtered = os.path.join(self.dir_output, f"{self.database_name}_filtered.fna")
        self.fasta_parser.write_fasta_sequences(
            fasta_sequences=fasta_sequences_filtered, file_out=file_database_filtered
        )

        return file_database_filtered

    def _filter_fasta_database_by_property_category(
        self,
        property_name: str,
        property_category: str | list[str],
        keep_if_equals_category: bool,
    ) -> str:
        """
        Filter a FASTA database based on property categories in sequence headers.
        Therefore, merge fasta files, load fasta content and process header for filtering.

        :param property_name: The property in the FASTA header used for filtering.
        :type property_name: str
        :param property_category: list of property categories used to filter sequences.
        :type property_category: str | list[str]
        :param keep_if_equals_category: Whether to keep (True) or remove (False) sequences matching the property category.
        :type keep_if_equals_category: bool
        """
        if self.database_file:
            file: str = self.database_file
        else:
            raise DatabaseError("Database file is not set. Call load_database_from_file() first.")

        fasta_sequences = self.fasta_parser.read_fasta_sequences(file_fasta_in=file)

        fasta_sequences_filtered = []

        for entry in fasta_sequences:
            _, properties, _ = self.fasta_parser.parse_fasta_header(entry.id, parse_additional_info=True)
            if isinstance(properties, str):
                continue
            if property_name in properties:
                property_values = cast_to_list(properties[property_name])
                if keep_if_equals_category and any(item in property_category for item in property_values):
                    fasta_sequences_filtered.append(entry)
                elif not keep_if_equals_category and all(
                    item not in property_category for item in property_values
                ):
                    fasta_sequences_filtered.append(entry)

        file_database_filtered = os.path.join(self.dir_output, f"{self.database_name}_filtered.fna")
        self.fasta_parser.write_fasta_sequences(
            fasta_sequences=fasta_sequences_filtered, file_out=file_database_filtered
        )

        return file_database_filtered
