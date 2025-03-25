############################################
# imports
############################################

import os
import shutil

from pathlib import Path
from typing import Union, List, get_args

from oligo_designer_toolsuite.utils import FastaParser, VCFParser, check_if_list
from oligo_designer_toolsuite._constants import _TYPES_REF

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
        self.dir_output = os.path.abspath(dir_output)

        self.database_file = None
        self.database_type = None

        self.fasta_parser = FastaParser()
        self.vcf_parser = VCFParser()

    def __del__(self):
        """Destructor for the ReferenceDatabase class."""
        if os.path.exists(self.database_file):
            os.remove(self.database_file)

    def load_database_from_file(
        self,
        files: Union[str, List[str]],
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
        :type files: Union[str, List[str]]
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
        files = check_if_list(files)

        # remove all files if database should be overwritten
        if self.database_file is not None and database_overwrite:
            os.remove(self.database_file)
            self.database_file = None
            self.database_type = None

        # set retrieve input file list
        if self.database_file is None:
            files_in = files
            self.database_type = file_type
        elif self.database_file is not None and self.database_type == file_type:
            files_in = files.extend(self.database_file)
        else:
            raise ValueError(f"Cannot mix {file_type} and {self.database_type} databases!")

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
            ValueError(f"Database type {self.database_type} not supported.")

    def write_database_to_file(self, filename: str, dir_output: str = None) -> str:
        """
        Write the loaded database to a file based on its type.

        :param filename: Name of the output file (without extension).
        :type filename: str
        :param dir_output: The directory where the database will be stored, if None, the directory defined in the init function will be used.
        :type dir_output: str
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
            return file_database
        else:
            raise ValueError("Database is empty! Nothing to be written to file.")

    def filter_database_by_region(self, region_ids: Union[str, List[str]], keep_region: bool) -> None:
        """
        Filter the database to retain or remove specific regions.

        :param region_ids: List of region identifiers to keep or remove.
        :type region_ids: Union[str, List[str]]
        :param keep_region: Whether to keep (True) or remove (False) the specified regions.
        :type keep_region: bool
        :return: Path to the filtered database file.
        :rtype: str
        :raises ValueError: If the database is empty or filtering is attempted on a non-FASTA database.
        """
        region_ids = check_if_list(region_ids)

        if self.database_file:
            if self.database_type == "fasta":
                file_database_filtered = self._filter_fasta_database_by_region(
                    region_ids=region_ids,
                    keep_region=keep_region,
                )
                os.remove(self.database_file)
                self.database_file = file_database_filtered
            else:
                raise ValueError("Filter only available for database in fasta format.")
        else:
            raise ValueError(
                "Can not filter. Database is empty! Call the method load_database_from_file() first."
            )

    def filter_database_by_attribute_category(
        self,
        attribute_name: str,
        attribute_category: Union[str, List[str]],
        keep_if_equals_category: bool,
    ) -> None:
        """
        Filter the database to retain or remove specific attribute categories.

        :param attribute_name: The attribute used for filtering.
        :type attribute_name: str
        :param attribute_category: List of attribute values to filter by.
        :type attribute_category: Union[str, List[str]]
        :param keep_if_equals_category: Whether to keep (True) or remove (False) records with the specified attribute values.
        :type keep_if_equals_category: bool
        :return: Path to the filtered database file.
        :rtype: str
        :raises ValueError: If the database is empty or filtering is attempted on a non-FASTA database.
        """
        attribute_category = check_if_list(attribute_category)

        if self.database_file:
            if self.database_type == "fasta":
                file_database_filtered = self._filter_fasta_database_by_attribute_category(
                    attribute_name=attribute_name,
                    attribute_category=attribute_category,
                    keep_if_equals_category=keep_if_equals_category,
                )
                os.remove(self.database_file)
                self.database_file = file_database_filtered
            else:
                raise ValueError("Filter only available for database in fasta format.")
        else:
            raise ValueError(
                "Can not filter. Database is empty! Call the method load_database_from_file() first."
            )

        return file_database_filtered

    def _filter_fasta_database_by_region(
        self,
        region_ids: Union[str, List[str]],
        keep_region: bool,
    ):
        """
        Filter a FASTA database to retain or exclude specific regions.
        Therefore, merge fasta files and load fasta content for filtering.

        :param region_ids: List of region identifiers to retain or remove.
        :type region_ids: Union[str, List[str]]
        :param keep_region: Whether to keep (True) or remove (False) the specified regions.
        :type keep_region: bool
        """
        if keep_region:
            regions_to_keep = region_ids
        else:
            regions_fasta = self.fasta_parser.get_fasta_regions(file_fasta_in=self.database_file)
            regions_to_keep = [region for region in regions_fasta if region not in region_ids]

        fasta_sequences_filtered = self.fasta_parser.read_fasta_sequences(
            file_fasta_in=self.database_file, region_ids=regions_to_keep
        )

        file_database_filtered = os.path.join(self.dir_output, f"{self.database_name}_filtered.fna")
        self.fasta_parser.write_fasta_sequences(
            fasta_sequences=fasta_sequences_filtered, file_out=file_database_filtered
        )

        return file_database_filtered

    def _filter_fasta_database_by_attribute_category(
        self,
        attribute_name: str,
        attribute_category: Union[str, List[str]],
        keep_if_equals_category: bool,
    ):
        """
        Filter a FASTA database based on attribute categories in sequence headers.
        Therefore, merge fasta files, load fasta content and process header for filtering.

        :param attribute_name: The attribute in the FASTA header used for filtering.
        :type attribute_name: str
        :param attribute_category: List of attribute values used to filter sequences.
        :type attribute_category: Union[str, List[str]]
        :param keep_if_equals_category: Whether to keep (True) or remove (False) sequences matching the attribute category.
        :type keep_if_equals_category: bool
        """
        fasta_sequences = self.fasta_parser.read_fasta_sequences(file_fasta_in=self.database_file)

        fasta_sequences_filtered = []

        for entry in fasta_sequences:
            _, attributes, _ = self.fasta_parser.parse_fasta_header(entry.id, parse_additional_info=True)
            if attribute_name in attributes:
                attribute_values = check_if_list(attributes[attribute_name])
                if keep_if_equals_category and any(item in attribute_category for item in attribute_values):
                    fasta_sequences_filtered.append(entry)
                elif not keep_if_equals_category and all(
                    item not in attribute_category for item in attribute_values
                ):
                    fasta_sequences_filtered.append(entry)

        file_database_filtered = os.path.join(self.dir_output, f"{self.database_name}_filtered.fna")
        self.fasta_parser.write_fasta_sequences(
            fasta_sequences=fasta_sequences_filtered, file_out=file_database_filtered
        )

        return file_database_filtered
