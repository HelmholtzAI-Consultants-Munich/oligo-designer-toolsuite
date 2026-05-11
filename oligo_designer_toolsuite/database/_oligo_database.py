############################################
# imports
############################################

import os
import pickle
import warnings
from pathlib import Path
from typing import Any

import pandas as pd
import yaml
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from effidict import EffiDict, LRUReplacement, PickleBackend
from joblib import Parallel, delayed
from joblib_progress import joblib_progress

from oligo_designer_toolsuite._constants import SEPARATOR_OLIGO_ID
from oligo_designer_toolsuite._exceptions import DatabaseError, FileFormatError
from oligo_designer_toolsuite.utils import (
    CustomYamlDumper,
    FastaParser,
    cast_to_list,
    cast_to_list_of_lists,
    check_if_key_in_database,
    check_if_region_in_database,
    check_tsv_format,
    collapse_properties_for_duplicated_sequences,
    flatten_property_list,
    format_oligo_properties,
    merge_databases,
)

CustomYamlDumper.add_representer(list, CustomYamlDumper.represent_list)
CustomYamlDumper.add_representer(dict, CustomYamlDumper.represent_dict)

############################################
# Oligo Database Class
############################################


class OligoDatabase:
    """
    The `OligoDatabase` class provides a comprehensive framework for managing oligonucleotide data along with their associated properties
    across various genomic regions. It supports loading, saving, filtering, and updating oligo data, ensuring efficient handling of large
    datasets through the use of an LRU (Least Recently Used) cache system. It provides functionalities to
    load data from different sources (such as FASTA files and tables), save data in various formats, and apply
    filtering based on specific criteria like regions, oligo IDs, or property thresholds. Additionally, the class
    includes methods to update, retrieve, and analyze oligo properties, facilitating the selection and evaluation
    of oligonucleotides for research and practical applications.

    :param min_oligos_per_region: Minimum number of oligos required per region to retain the region in the database and oligosets.
    :type min_oligos_per_region: int
    :param write_regions_with_insufficient_oligos: Flag to log regions with insufficient oligos.
    :type write_regions_with_insufficient_oligos: bool
    :param max_entries_in_memory: Maximum number of database entries to keep in memory.
    :type max_entries_in_memory: int
    :param n_jobs: Number of parallel jobs to use for processing.
    :type n_jobs: int
    :param database_name: Name of the database for storing oligo data.
    :type database_name: str
    :param dir_output: Directory path where output files will be saved.
    :type dir_output: str
    """

    def __init__(
        self,
        min_oligos_per_region: int = 0,
        write_regions_with_insufficient_oligos: bool = True,
        max_entries_in_memory: int = 10,
        n_jobs: int = 1,
        database_name: str = "db_oligo",
        dir_output: str = "output",
    ) -> None:
        """Constructor for the OligoDatabase class."""

        self.min_oligos_per_region = min_oligos_per_region
        self.write_regions_with_insufficient_oligos = write_regions_with_insufficient_oligos
        self._max_entries_in_memory = max_entries_in_memory
        self.n_jobs = n_jobs

        self.database_name = database_name
        self.dir_output = os.path.abspath(os.path.join(dir_output, database_name))
        Path(self.dir_output).mkdir(parents=True, exist_ok=True)

        self._dir_cache_files = os.path.join(self.dir_output, "cache_files")

        self.fasta_parser = FastaParser()

        # Initialize databse object
        backend = PickleBackend(storage_path=self._dir_cache_files)
        strategy = LRUReplacement(disk_backend=backend, max_in_memory=self._max_entries_in_memory)
        self.database = EffiDict(disk_backend=backend, replacement_strategy=strategy)

        # will be used later in the generation of oligo sets
        backend = PickleBackend(storage_path=self._dir_cache_files)
        strategy = LRUReplacement(disk_backend=backend, max_in_memory=self._max_entries_in_memory)
        self.oligosets = EffiDict(disk_backend=backend, replacement_strategy=strategy)

        # Database metadata
        self.database_sequence_types: list[str] = []

        # Initialize the file for regions with insufficient oligos
        if self.write_regions_with_insufficient_oligos:
            self.file_removed_regions = os.path.join(
                self.dir_output,
                f"regions_with_insufficient_oligos_for_{self.database_name}.txt",
            )
            with open(self.file_removed_regions, "a") as handle:
                handle.write(f"Region\tPipeline step\n")

    ############################################
    # Load Functions
    ############################################

    def load_database_from_fasta(
        self,
        files_fasta: str | list[str],
        database_overwrite: bool,
        sequence_type: str,
        region_ids: str | list[str] | None = None,
    ) -> None:
        """
        Loads oligonucleotide data from one or more FASTA files into the database, optionally overwriting the existing database.

        This function reads sequences from FASTA file(s) and adds them to the OligoDatabase, either as 'oligo' or
        'target' sequence. It parses the headers of the FASTA entries to extract oligo property information,
        and assigns unique IDs to the oligos within each region.

        The header of each sequence must start with '>' and contain the following information:
        region_id, additional_information (optional) and coordinates (chrom, start, end, strand),
        where the region_id is compulsory and the other fileds are opional.

        Input Format (per sequence):
        >region_id::additional information::chromosome:start-end(strand)
        sequence

        Example:
        >ASR1::transcrip_id=XM456,exon_number=5::16:54552-54786(+)
        AGTTGACAGACCCCAGATTAAAGTGTGTCGCGCAACAC

        :param files_fasta: Path(s) to FASTA file(s) to load. Can be a single file path (str) or a list of file paths (list[str]).
        :type files_fasta: str | list[str]
        :param database_overwrite: If True, the existing database will be overwritten.
        :type database_overwrite: bool
        :param sequence_type: Type of sequence being processed.
        :type sequence_type: str
        :param region_ids: Region identifier(s) to process. Can be a single region ID (str) or a list of region IDs (list[str]). If None, all regions in the database are processed, defaults to None.
        :type region_ids: str | list[str] | None, optional
        """

        def _load_fasta_file(file: str) -> None:
            """
            Load a FASTA file and integrate its sequences into the existing database.

            This function checks the format of a provided FASTA file and, if valid, reads the sequences.
            It then processes each sequence to extract regions, additional information, and coordinates from the headers.
            The sequences are stored in a structured database, ensuring that any duplicated sequences within the same region are appropriately merged.

            :param file: The path to the FASTA file to be loaded.
            :type file: str
            :return: None
            """
            if self.fasta_parser.check_fasta_format(file):
                fasta_sequences = self.fasta_parser.read_fasta_sequences(file, region_ids)
                sequences: dict[str, dict[str, Any]] = {}
                for entry in fasta_sequences:
                    region, additional_info, coordinates = self.fasta_parser.parse_fasta_header(entry.id)
                    if isinstance(additional_info, str):
                        oligo_properties = coordinates
                    else:
                        oligo_properties = coordinates | additional_info
                    oligo_properties = format_oligo_properties(oligo_properties, self.database_sequence_types)
                    if region in sequences:
                        if entry.seq in sequences[region]:
                            oligo_properties = collapse_properties_for_duplicated_sequences(
                                oligo_properties1=sequences[region][entry.seq],
                                oligo_properties2=oligo_properties,
                                database_sequence_types=self.database_sequence_types,
                            )
                        sequences[region][str(entry.seq)] = oligo_properties
                    else:
                        sequences[region] = {str(entry.seq): oligo_properties}

                database_region: dict[str, dict[str, Any]] = {region: {} for region in sequences.keys()}
                for region, sequences_region in sequences.items():
                    i = 1
                    for oligo_sequence, oligo_properties in sequences_region.items():
                        oligo_id = f"{region}{SEPARATOR_OLIGO_ID}{i}"
                        oligo_seq_info = {sequence_type: oligo_sequence} | oligo_properties
                        database_region[region][oligo_id] = oligo_seq_info
                        i += 1

                # only merge if there are common keys
                if len(set(self.database) & set(database_region)) > 0:
                    self.database = merge_databases(
                        database1=self.database,
                        database2=database_region,
                        sequence_type=sequence_type,
                        database_sequence_types=self.database_sequence_types,
                        dir_cache_files=self._dir_cache_files,
                        max_entries_in_memory=self._max_entries_in_memory,
                    )
                else:
                    for region in database_region.keys():
                        self.database[region] = database_region[region]

        # Check formatting
        region_ids = cast_to_list(region_ids) if region_ids else None
        files_fasta = cast_to_list(files_fasta)

        # Set database sequence types
        self.set_database_sequence_types(sequence_type)

        # Clear database if it should be overwritten
        if database_overwrite:
            backend = PickleBackend(storage_path=self._dir_cache_files)
            strategy = LRUReplacement(disk_backend=backend, max_in_memory=self._max_entries_in_memory)
            self.database = EffiDict(disk_backend=backend, replacement_strategy=strategy)

        # Load files parallel into database
        with joblib_progress(description=f"Database Loading", total=len(files_fasta)):
            Parallel(n_jobs=self.n_jobs, prefer="threads", require="sharedmem")(
                delayed(_load_fasta_file)(file_fasta) for file_fasta in files_fasta
            )

        # add this step to log regions which are not available in database
        if region_ids:
            check_if_region_in_database(
                database=self.database,
                region_ids=region_ids,
                write_regions_with_insufficient_oligos=self.write_regions_with_insufficient_oligos,
                file_removed_regions=self.file_removed_regions,
            )

    def load_database_from_table(
        self,
        file_database: str,
        database_overwrite: bool,
        merge_databases_on_sequence_type: str,
        region_ids: str | list[str] | None = None,
    ) -> None:
        """
        Loads oligonucleotide data from a tab-delimited table (TSV) file into the database, optionally overwriting the existing database.

        This function loads the OligoDatabase from a tab-separated values (TSV) file. The file must contain
        columns such as 'region_id', 'oligo_id' and optionally the oligo or target sequence as well as oligo properties.
        The database can be optionally filtered by specifying a list of region IDs.

        ⚠️ For big databases, it is not recommended to load the whole TSV file at once. Instead, the database should be
        split into smaller files. The function will merge the databases if there are common keys (regions).

        :param file_database: Path to the TSV file containing the database.
        :type file_database: str
        :param database_overwrite: If True, the existing database will be overwritten.
        :type database_overwrite: bool
        :param merge_databases_on_sequence_type: The sequence type on which two databases should be merged on if database_overwrite = False,
                must be one of the predefined sequence types, i.e. "oligo" or "target".
        :type sequence_type: str
        :param region_ids: Region identifier(s) to process. Can be a single region ID (str) or a list of region IDs (list[str]). If None, all regions in the database are processed, defaults to None.
        :type region_ids: str | list[str] | None, optional
        """
        # Check formatting
        region_ids = cast_to_list(region_ids) if region_ids else None

        # Check if file exists and has correct format
        if os.path.exists(file_database):
            if not check_tsv_format(file_database):
                raise FileFormatError(
                    f"Database file '{file_database}' has incorrect format. Expected TSV format."
                )
        else:
            raise FileFormatError(f"Database file '{file_database}' does not exist.")

        # Clear database if it should be overwritten
        if database_overwrite:
            backend = PickleBackend(storage_path=self._dir_cache_files)
            strategy = LRUReplacement(disk_backend=backend, max_in_memory=self._max_entries_in_memory)
            self.database = EffiDict(disk_backend=backend, replacement_strategy=strategy)

        # Load file and process content
        file_tsv_content = pd.read_table(file_database, sep="\t")
        file_tsv_content = file_tsv_content.apply(
            lambda col: col.apply(lambda x: x if pd.notna(x) else [None])
        )

        # convert lists represented as string to proper list format in the table with the eval function
        file_tsv_content = file_tsv_content.apply(
            lambda col: col.apply(
                lambda x: (
                    eval(x)
                    if isinstance(x, str) and x.startswith("[") and x.endswith("]")
                    else ([int(x)] if isinstance(x, str) and x.isdigit() else x)
                )
            )
        )

        # Merge loaded database with existing one
        database_tmp1 = file_tsv_content.to_dict(orient="records")

        backend = PickleBackend(storage_path=self._dir_cache_files)
        strategy = LRUReplacement(disk_backend=backend, max_in_memory=self._max_entries_in_memory)
        database_tmp2 = EffiDict(disk_backend=backend, replacement_strategy=strategy)

        for entry in database_tmp1:
            region_id, oligo_id = entry.pop("region_id"), entry.pop("oligo_id")
            if (not region_ids) or (region_ids and region_id in region_ids):
                if region_id not in database_tmp2:
                    database_tmp2[region_id] = {}
                database_tmp2[region_id][oligo_id] = format_oligo_properties(
                    entry, self.database_sequence_types
                )

        if not database_overwrite and self.database:
            database_tmp2 = merge_databases(
                database1=self.database,
                database2=database_tmp2,
                sequence_type=merge_databases_on_sequence_type,
                database_sequence_types=self.database_sequence_types,
                dir_cache_files=self._dir_cache_files,
                max_entries_in_memory=self._max_entries_in_memory,
            )

        # Filter for region ids
        if region_ids:
            check_if_region_in_database(
                database=database_tmp2,
                region_ids=region_ids,
                write_regions_with_insufficient_oligos=self.write_regions_with_insufficient_oligos,
                file_removed_regions=self.file_removed_regions,
            )

        self.database = database_tmp2

    def load_database(
        self,
        dir_database: str,
        database_overwrite: bool,
        merge_databases_on_sequence_type: str,
        region_ids: str | list[str] | None = None,
    ) -> None:
        """
        Loads oligonucleotide data from a directory of pickled files into the database, optionally overwriting the existing database.
        Each file in the folder represents a region in the database and contains the region ID,
        the oligo sequence and property information as well as the oligosets for this region.

        :param dir_database: Path to the directory containing the pickled database files.
        :type dir_database: str
        :param database_overwrite: If True, the existing database will be overwritten.
        :type database_overwrite: bool
        :param merge_databases_on_sequence_type: The sequence type on which two databases should be merged on (if database_overwrite = False),
                must be one of the predefined sequence types, i.e. "oligo" or "target".
        :type sequence_type: str
        :param region_ids: Region identifier(s) to process. Can be a single region ID (str) or a list of region IDs (list[str]). If None, all regions in the database are processed, defaults to None.
        :type region_ids: str | list[str] | None, optional
        """

        def _load_database_file(file: str) -> None:
            """
            Loads a database file and integrates its content into the existing database and oligosets.

            This function opens a pickled database file, extracts the region ID, and the associated database and oligoset regions.
            If the region ID is already present in the current database, it merges the new data with the existing one.
            Otherwise, it adds the new region data to the database and oligosets.

            :param file: The path to the database file to be loaded.
            :type file: str
            :return: None
            """
            # extract region ID and content from the file
            with open(file, "rb") as handle:
                content = pickle.load(handle)
                region_id = content["region_id"]
                database_region = content["database_region"]
                oligoset_region = content["oligoset_region"]

            # Only process selected regions
            if not region_ids or region_id in region_ids:
                # only merge if there are common keys
                if region_id in self.database.keys():
                    self.database = merge_databases(
                        database1=self.database,
                        database2={region_id: database_region},
                        sequence_type=merge_databases_on_sequence_type,
                        database_sequence_types=self.database_sequence_types,
                        dir_cache_files=self._dir_cache_files,
                        max_entries_in_memory=self._max_entries_in_memory,
                    )
                    self.oligosets[region_id] = pd.concat([self.oligosets[region_id], oligoset_region])
                else:
                    self.database[region_id] = database_region
                    self.oligosets[region_id] = oligoset_region

        # Check formatting
        region_ids = cast_to_list(region_ids) if region_ids else None

        if not os.path.isdir(dir_database):
            raise DatabaseError(f"Database directory '{dir_database}' does not exist.")

        if database_overwrite:
            backend = PickleBackend(storage_path=self._dir_cache_files)
            strategy = LRUReplacement(disk_backend=backend, max_in_memory=self._max_entries_in_memory)
            self.database = EffiDict(disk_backend=backend, replacement_strategy=strategy)

        # retrieve all files in the directory
        path = os.path.abspath(dir_database)
        files_database = [entry.path for entry in os.scandir(path) if entry.is_file()]

        # Load files parallel into database
        with joblib_progress(description=f"Database Loading", total=len(files_database)):
            Parallel(n_jobs=self.n_jobs, prefer="threads", require="sharedmem")(
                delayed(_load_database_file)(file_database) for file_database in files_database
            )

        # add this step to log regions which are not available in database
        if region_ids:
            check_if_region_in_database(
                database=self.database,
                region_ids=region_ids,
                write_regions_with_insufficient_oligos=self.write_regions_with_insufficient_oligos,
                file_removed_regions=self.file_removed_regions,
            )

    ############################################
    # Save Functions
    ############################################

    def save_database(
        self,
        name_database: str = "db_oligo",
        dir_output: str | None = None,
        region_ids: str | list[str] | None = None,
    ) -> str:
        """
        Saves the current database and associated oligosets as pickl files into a specified directory.

        :param name_database: Directory path where the database files should be saved. Default is "db_oligo".
        :type name_database: str
        :param dir_output: Directory path where output files will be saved.
        :type dir_output: str
        :param region_ids: Region identifier(s) to process. Can be a single region ID (str) or a list of region IDs (list[str]). If None, all regions in the database are processed, defaults to None.
        :type region_ids: str | list[str] | None, optional
        :return: The directory path where the database was saved.
        :rtype: str
        """
        # Check formatting
        region_ids = cast_to_list(region_ids) if region_ids else self.database.keys()

        if dir_output:
            dir_database = os.path.join(dir_output, name_database)
        else:
            dir_database = os.path.join(self.dir_output, name_database)
        Path(dir_database).mkdir(parents=True, exist_ok=True)

        for region_id in region_ids:
            database_region = self.database[region_id]
            if self.oligosets and region_id in self.oligosets:
                oligoset_region = self.oligosets[region_id]
            else:
                oligoset_region = None
            file_output = os.path.join(dir_database, region_id)
            with open(file_output, "wb") as file:
                pickle.dump(
                    {
                        "region_id": region_id,
                        "database_region": database_region,
                        "oligoset_region": oligoset_region,
                    },
                    file,
                )

        return dir_database

    def write_database_to_fasta(
        self,
        sequence_type: str,
        save_description: bool,
        filename: str = "db_oligo",
        dir_output: str | None = None,
        region_ids: str | list[str] | None = None,
    ) -> str:
        """
        Writes the current database to a FASTA file. Associated sequence properties can optionally be included in the sequence header.

        :param sequence_type: Type of sequence being processed.
        :type sequence_type: str
        :param save_description: Whether to include the sequence properties in the sequence header.
        :type save_description: bool
        :param filename: The base name of the output FASTA file, defaults to "db_oligo".
        :type filename: str
        :param dir_output: Directory path where output files will be saved.
        :type dir_output: str
        :param region_ids: Region identifier(s) to process. Can be a single region ID (str) or a list of region IDs (list[str]). If None, all regions in the database are processed, defaults to None.
        :type region_ids: str | list[str] | None, optional
        :return: The path to the saved FASTA file.
        :rtype: str
        """
        # Check if sequence type exists in database
        assert check_if_key_in_database(
            self.database, sequence_type
        ), f"Sequence type '{sequence_type}' not found in database."

        # Check formatting
        region_ids = cast_to_list(region_ids) if region_ids else self.database.keys()

        dir_output = dir_output if dir_output else self.dir_output
        file_fasta = os.path.join(dir_output, f"{filename}.fna")
        output_fasta = []

        with open(file_fasta, "w") as handle_fasta:
            for region_id in region_ids:
                database_region = self.database[region_id]
                for oligo_id, oligo_properties in database_region.items():
                    description = sequence_type if save_description else ""
                    seq_record = SeqRecord(
                        Seq(oligo_properties[sequence_type]),
                        id=oligo_id,
                        name=oligo_id.split(SEPARATOR_OLIGO_ID)[0],
                        description=description,
                    )
                    output_fasta.append(seq_record)

            SeqIO.write(output_fasta, handle_fasta, "fasta")

        return file_fasta

    def write_database_to_bed(
        self,
        filename: str = "db_oligo",
        dir_output: str | None = None,
        region_ids: str | list[str] | None = None,
    ) -> str:
        """
        Write the oligo database to a BED file format.

        This function exports oligo sequences with their genomic coordinates to a BED file,
        which can be used for downstream analysis. Entries which contain None for either
        chromosome, start, end or strand will be removed before the content is written to BED file.

        :param filename: Name of the output BED file (without extension). Default is "db_oligo".
        :type filename: str
        :param dir_output: Directory path where output files will be saved.
        :type dir_output: str
        :param region_ids: Region identifier(s) to process. Can be a single region ID (str) or a list of region IDs (list[str]). If None, all regions in the database are processed, defaults to None.
        :type region_ids: str | list[str] | None, optional
        :return: The path to the saved BED file.
        :rtype: str
        """
        # Check formatting
        region_ids = cast_to_list(region_ids) if region_ids else self.database.keys()

        dir_output = dir_output if dir_output else self.dir_output
        file_bed = os.path.join(dir_output, f"{filename}.bed")

        # retrieve relevant information from database
        property_table = self.get_oligo_property_table(
            properties=["chromosome", "start", "end", "strand"], flatten=True, region_ids=region_ids
        )

        # remove rows that contain None
        mask = property_table.isnull().any(axis=1)
        property_table = property_table[~mask]

        if mask.sum() > 0:
            warnings.warn(f"Removing {mask.sum()} row(s) containing None/NaN values.", UserWarning)

        # expand rows which contain lists for chr, start, end, strand columns into seperate rows
        property_table_extended = property_table.explode(
            ["chromosome", "start", "end", "strand"], ignore_index=True
        )
        property_table_extended["score"] = "."

        # save tabel content as BED file
        property_table_extended[["chromosome", "start", "end", "oligo_id", "score", "strand"]].to_csv(
            file_bed, sep="\t", header=False, index=False
        )

        return file_bed

    def write_database_to_table(
        self,
        properties: str | list[str],
        flatten_property: bool,
        filename: str = "oligo_database_table",
        dir_output: str | None = None,
        region_ids: str | list[str] | None = None,
    ) -> str:
        """
        Writes the current database and selected properties to a TSV table file,
        optionally flattening nested or list properties to a unique set of values.


        :param properties: A list of properties to include in the table.
        :type properties: str | list[str]
        :param flatten_property: Whether to flatten list properties to a unique set of values in the table.
        :type flatten_property: bool
        :param filename: The base name of the output TSV file, defaults to "oligo_database_table".
        :type filename: str
        :param dir_output: Directory path where output files will be saved.
        :type dir_output: str
        :param region_ids: Region identifier(s) to process. Can be a single region ID (str) or a list of region IDs (list[str]). If None, all regions in the database are processed, defaults to None.
        :type region_ids: str | list[str] | None, optional
        :return: The path to the saved TSV file.
        :rtype: str
        """
        # Check formatting
        region_ids = cast_to_list(region_ids) if region_ids else self.database.keys()
        properties = cast_to_list(properties)

        dir_output = dir_output if dir_output else self.dir_output
        file_table = os.path.join(os.path.dirname(dir_output), f"{filename}.tsv")

        first_entry = True
        for region_id in region_ids:
            file_tsv_content = []
            for oligo_id in self.database[region_id].keys():
                entry = {"region_id": region_id, "oligo_id": oligo_id}
                for property in properties:
                    if property in self.database[region_id][oligo_id]:
                        oligo_property = self.database[region_id][oligo_id][property]
                        if flatten_property:
                            oligo_property = flatten_property_list(oligo_property)
                            if oligo_property:
                                oligo_property = list(set(oligo_property))
                            entry[property] = (
                                str(oligo_property).replace("'", "").replace("[", "").replace("]", "")
                            )
                        else:
                            entry[property] = str(oligo_property).replace("[[", "[").replace("]]", "]")
                file_tsv_content.append(entry)
            file_tsv_content_as_df = pd.DataFrame(data=file_tsv_content)
            file_tsv_content_as_df.to_csv(file_table, sep="\t", index=False, mode="a", header=first_entry)
            first_entry = False

        return file_table

    def write_oligosets_to_yaml(
        self,
        properties: str | list[str],
        ascending: bool,
        filename: str = "oligosets",
        dir_output: str | None = None,
        region_ids: str | list[str] | None = None,
    ) -> None:
        """
        Writes the current top n oligosets and selected properties to a YAML file.
        The oligosets are sorted based on their scores in ascending or descending order.

        :param properties: A list of properties to include in the YAML file.
        :type properties: str | list[str]
        :param ascending: If True, sort oligosets by score in ascending order (the smaller the score the better the oligo set);
            otherwise, in descending order (the higher the score the better the oligo set).
        :type ascending: bool
        :param filename: Base name for the output YAML file, defaults to "oligosets".
        :type filename: str
        :param dir_output: Directory path where output files will be saved.
        :type dir_output: str
        :param region_ids: Region identifier(s) to process. Can be a single region ID (str) or a list of region IDs (list[str]). If None, all regions in the database are processed, defaults to None.
        :type region_ids: str | list[str] | None, optional
        :return: None
        """
        # Check formatting
        properties = cast_to_list(properties)
        region_ids = cast_to_list(region_ids) if region_ids else self.database.keys()

        yaml_dict: dict[str, dict[str, Any]] = {region_id: {} for region_id in region_ids}

        for region_id in region_ids:
            oligosets_region = self.oligosets[region_id]
            oligosets_oligo_columns = [col for col in oligosets_region.columns if col.startswith("oligo_")]
            oligosets_score_columns = [
                col for col in oligosets_region.columns if col.startswith("set_score_")
            ]

            oligosets_region = oligosets_region.sort_values(by=oligosets_score_columns, ascending=ascending)
            oligosets_region_oligos = oligosets_region[oligosets_oligo_columns]
            oligosets_region_scores = oligosets_region[oligosets_score_columns]

            for idx, oligoset in oligosets_region_oligos.iterrows():
                oligoset_id = f"Oligoset {idx + 1}"
                yaml_dict[region_id][oligoset_id] = {
                    "Oligoset Score": oligosets_region_scores.loc[idx].to_dict(),
                }

                for oligo_idx, oligo_id in enumerate(oligoset):
                    yaml_dict_oligo_entry = {"oligo_id": oligo_id}

                    # iterate through all properties that should be written
                    for property in properties:
                        if property in self.database[region_id][oligo_id]:
                            oligo_property = self.database[region_id][oligo_id][property]
                            # format oligo properties: flatten lists of lists, join string lists with comma, keep strings as-is, None -> empty list
                            if oligo_property:
                                if (
                                    sum(len(sublist) for sublist in cast_to_list_of_lists(oligo_property))
                                    == 1
                                ):
                                    oligo_property = flatten_property_list(oligo_property)
                            yaml_dict_oligo_entry[property] = oligo_property

                    oligo_id_yaml = f"Oligo {oligo_idx + 1}"
                    yaml_dict[region_id][oligoset_id][oligo_id_yaml] = yaml_dict_oligo_entry

        dir_output = dir_output if dir_output else self.dir_output
        file_yaml = os.path.join(os.path.dirname(dir_output), f"{filename}.yml")

        with open(file_yaml, "w") as handle:
            yaml.dump(yaml_dict, handle, Dumper=CustomYamlDumper, default_flow_style=False, sort_keys=False)

    def write_oligosets_to_table(
        self,
        properties: str | list[str],
        ascending: bool,
        filename: str = "oligosets",
        dir_output: str | None = None,
        region_ids: str | list[str] | None = None,
    ) -> None:
        # Check formatting
        properties = cast_to_list(properties)
        region_ids = cast_to_list(region_ids) if region_ids else self.database.keys()

        csv_table = list()

        for region_id in region_ids:
            oligosets_region = self.oligosets[region_id]
            oligosets_oligo_columns = [col for col in oligosets_region.columns if col.startswith("oligo_")]
            oligosets_score_columns = [
                col for col in oligosets_region.columns if col.startswith("set_score_")
            ]

            oligosets_region = oligosets_region.sort_values(by=oligosets_score_columns, ascending=ascending)
            oligosets_region = oligosets_region[oligosets_oligo_columns]
            oligosets_region = oligosets_region.reset_index(drop=True)

            # iterate through all oligo sets
            for oligoset_idx, oligoset in oligosets_region.iterrows():
                oligoset_id = f"oligoset_{oligoset_idx + 1}"
                for oligo_id in oligoset:
                    entry = {
                        "region_id": region_id,
                        "oligoset_id": oligoset_id,
                        "oligo_id": oligo_id,
                    }

                    for property in properties:
                        if property in self.database[region_id][oligo_id]:
                            oligo_property = self.database[region_id][oligo_id][property]
                            # format oligo properties: flatten lists of lists, join string lists with comma, keep strings as-is, None -> empty list
                            if oligo_property:
                                if (
                                    sum(len(sublist) for sublist in cast_to_list_of_lists(oligo_property))
                                    == 1
                                ):
                                    oligo_property = flatten_property_list(oligo_property)
                                    if len(oligo_property) == 1:
                                        oligo_property = oligo_property[0]
                            entry[property] = oligo_property

                    csv_table.append(entry)

        dir_output = dir_output if dir_output else self.dir_output
        file_table = os.path.join(os.path.dirname(dir_output), f"{filename}.tsv")

        csv_table_df = pd.DataFrame(csv_table)
        csv_table_df.to_csv(file_table, sep="\t", index=False)

        # Also write Excel file with one sheet per region_id
        file_excel = os.path.join(os.path.dirname(dir_output), f"{filename}.xlsx")
        try:
            with pd.ExcelWriter(file_excel, engine="openpyxl") as writer:
                for region_id in region_ids:
                    # Filter data for this region
                    region_data = csv_table_df[csv_table_df["region_id"] == region_id].copy()
                    if not region_data.empty:
                        # Remove region_id column from individual sheets since it's redundant
                        region_data = region_data.drop(columns=["region_id"])
                        # Write to sheet (Excel sheet names are limited to 31 characters and cannot contain certain characters)
                        sheet_name = str(region_id)[:31]
                        # Replace invalid characters for Excel sheet names
                        invalid_chars = ["\\", "/", "*", "[", "]", ":", "?"]
                        for char in invalid_chars:
                            sheet_name = sheet_name.replace(char, "_")
                        region_data.to_excel(writer, sheet_name=sheet_name, index=False)
        except ImportError:
            warnings.warn(
                "openpyxl is not installed. Excel file generation skipped. "
                "Install openpyxl to enable Excel export: pip install openpyxl",
                UserWarning,
            )
        except Exception as e:
            warnings.warn(
                f"Failed to write Excel file: {e}. TSV file was written successfully.",
                UserWarning,
            )

    def write_ready_to_order_yaml(
        self,
        properties: str | list[str],
        ascending: bool,
        filename: str = "ready_to_order",
        dir_output: str | None = None,
        region_ids: str | list[str] | None = None,
    ) -> None:
        """
        Writes a YAML file that only contains order information for the oligosets.

        :param filename: Base name for the output YAML file, defaults to "ready_to_order".
        :type filename: str
        :param dir_output: Directory path where output files will be saved.
        :type dir_output: str
        :return: None
        """
        properties = cast_to_list(properties)
        region_ids = cast_to_list(region_ids) if region_ids else self.database.keys()

        yaml_dict: dict[str, dict] = {}

        for region_id in region_ids:
            yaml_dict[region_id] = {}
            oligosets_region = self.oligosets[region_id]
            oligosets_oligo_columns = [col for col in oligosets_region.columns if col.startswith("oligo_")]
            oligosets_score_columns = [
                col for col in oligosets_region.columns if col.startswith("set_score_")
            ]

            oligosets_region = oligosets_region.sort_values(by=oligosets_score_columns, ascending=ascending)
            oligosets_region = oligosets_region[oligosets_oligo_columns]
            oligosets_region = oligosets_region.reset_index(drop=True)

            # iterate through all oligo sets
            for oligoset_idx, oligoset in oligosets_region.iterrows():
                oligoset_id = f"oligoset_{oligoset_idx + 1}"
                yaml_dict[region_id][oligoset_id] = {}
                for oligo_id in oligoset:
                    entry = {}
                    for property in properties:
                        if property in self.database[region_id][oligo_id]:
                            oligo_property = self.database[region_id][oligo_id][property]
                            # format oligo properties: flatten lists of lists, join string lists with comma, keep strings as-is, None -> empty list
                            if oligo_property:
                                if (
                                    sum(len(sublist) for sublist in cast_to_list_of_lists(oligo_property))
                                    == 1
                                ):
                                    oligo_property = flatten_property_list(oligo_property)
                            entry[property] = oligo_property
                    yaml_dict[region_id][oligoset_id][oligo_id] = entry

        dir_output = dir_output if dir_output else self.dir_output
        file_yaml = os.path.join(os.path.dirname(dir_output), f"{filename}.yml")

        with open(file_yaml, "w") as outfile:
            yaml.dump(yaml_dict, outfile, Dumper=CustomYamlDumper, default_flow_style=False, sort_keys=False)

    def remove_regions_with_insufficient_oligos(self, pipeline_step: str) -> None:
        """
        Removes regions from the database and oligoset collections if the number of oligos is less or equal to
        the specified minimum threshold (`min_oligos_per_region`). If the option `write_regions_with_insufficient_oligos`
        was enabled for the OligoDatabase class, it logs the removed regions along with the pipeline step in the `file_removed_regions`.

        :param pipeline_step: The name of the pipeline step during which the removal is performed, used for logging.
        :type pipeline_step: str
        """

        region_ids = list(self.database.keys())
        regions_to_remove = [
            region_id
            for region_id in region_ids
            if len(self.database[region_id]) < self.min_oligos_per_region
        ]

        for region in regions_to_remove:
            self.database[region] = None
            del self.database[region]

            self.oligosets[region] = None
            del self.oligosets[region]

        if self.write_regions_with_insufficient_oligos and regions_to_remove:
            with open(self.file_removed_regions, "a") as handle:
                handle.write("\n".join(f"{region}\t{pipeline_step}" for region in regions_to_remove) + "\n")

    ############################################
    # Getter Functions
    ############################################

    def get_property_list(self) -> list[str]:
        """Retrieves a list of property names stored in the database.

        :return: A list of property names stored in the database.
        :rtype: list[str]
        """
        region_id = next(iter(self.database.values()))
        oligo_id = next(iter(region_id.values()))
        properties = list(oligo_id.keys())

        return properties

    def get_regionid_list(self) -> list[str]:
        """
        Retrieves a list of all region IDs present in the database.

        :return: A list of region IDs in the database.
        :rtype: list[str]
        """
        region_ids = list(self.database.keys())

        return region_ids

    def get_oligoid_list(self, region_ids: str | list[str] | None = None) -> list[str]:
        """
        Retrieves a list of all oligo IDs present in the database for a specific region or all regions in the database.

        :param region_ids: Region identifier(s) to process. Can be a single region ID (str) or a list of region IDs (list[str]). If None, all regions in the database are processed, defaults to None.
        :type region_ids: str | list[str] | None, optional
        :return: A list of oligo IDs from all regions in the database.
        :rtype: list[str]
        """
        region_ids = cast_to_list(region_ids) if region_ids else list(self.database.keys())
        oligo_ids = [
            oligo_id
            for region_id in region_ids
            if self.database[region_id]  # skip regions without any oligos
            for oligo_id in self.database[region_id].keys()
        ]

        return oligo_ids

    def get_sequence_list(self, sequence_type: str) -> list[str]:
        """
        Retrieves a list of sequences of the specified type from the database.

        :param sequence_type: Type of sequence being processed.
        :type sequence_type: str
        :return: A list of sequences corresponding to the specified sequence type from all regions in the database.
        :rtype: list[str]
        """
        assert check_if_key_in_database(
            self.database, sequence_type
        ), f"Sequence type '{sequence_type}' not found in database."
        sequences = [
            str(oligo_properties[sequence_type])
            for region_id, database_region in self.database.items()
            if database_region  # skip regions without any oligos
            for oligo_id, oligo_properties in database_region.items()
        ]

        return sequences

    def get_oligoid_sequence_mapping(self, sequence_type: str, sequence_to_upper: bool = False) -> dict:
        """
        Generates a mapping of oligo IDs to their corresponding sequences of the specified type.

        :param sequence_type: Type of sequence being processed.
        :type sequence_type: str
        :param sequence_to_upper: Whether to convert sequences to uppercase, defaults to False.
        :type sequence_to_upper: bool
        :return: A dictionary mapping oligo IDs to their corresponding sequences.
        :rtype: dict
        """
        assert check_if_key_in_database(
            self.database, sequence_type
        ), f"Sequence type '{sequence_type}' not found in database."

        oligoid_sequence_mapping = {}

        for region_id, database_region in self.database.items():
            if not database_region:
                # Skip regions without any oligos
                continue
            for oligo_id, oligo_properties in database_region.items():
                seq = oligo_properties[sequence_type]
                if sequence_to_upper:
                    seq = seq.upper()
                oligoid_sequence_mapping[oligo_id] = seq

        return oligoid_sequence_mapping

    def get_sequence_oligoid_mapping(self, sequence_type: str, sequence_to_upper: bool = False) -> dict:
        """
        Generates a mapping of sequences to their corresponding oligo IDs for the specified sequence type.

        :param sequence_type: Type of sequence being processed.
        :type sequence_type: str
        :param sequence_to_upper: Whether to convert sequences to uppercase, defaults to False.
        :type sequence_to_upper: bool
        :return: A dictionary mapping sequences to their corresponding oligo IDs.
        :rtype: dict
        """
        assert check_if_key_in_database(
            self.database, sequence_type
        ), f"Sequence type '{sequence_type}' not found in database."

        sequence_oligoids_mapping = {}

        for region_id, database_region in self.database.items():
            if not database_region:
                # Skip regions without any oligos
                continue
            for oligo_id, oligo_properties in database_region.items():
                seq = oligo_properties[sequence_type]
                if sequence_to_upper:
                    seq = seq.upper()
                if seq not in sequence_oligoids_mapping:
                    # If the sequence key doesn't exist, create a new entry with the oligo ID in a list
                    sequence_oligoids_mapping[seq] = [oligo_id]
                else:
                    # If the sequence key already exists, append the oligo ID to the existing list
                    sequence_oligoids_mapping[seq].append(oligo_id)

        return sequence_oligoids_mapping

    def get_oligo_property_table(
        self,
        properties: str | list[str],
        flatten: bool,
        region_ids: str | list[str] | None = None,
    ) -> pd.DataFrame:
        """
        Generates a DataFrame containing oligo IDs and the specified property for each oligo,
        optionally flattening nested or list properties to a unique set of values.

        :param property: The name of the property to retrieve.
        :type property: str
        :param flatten: Whether to flatten list properties to a unique set of values in the table.
        :type flatten: bool
        :param region_ids: Region identifier(s) to process. Can be a single region ID (str) or a list of region IDs (list[str]). If None, all regions in the database are processed, defaults to None.
        :type region_ids: str | list[str] | None, optional
        :return: A DataFrame with oligo IDs and the corresponding property values.
        :rtype: pd.DataFrame
        """

        def _flatten_if_one(x: Any) -> Any:
            """
            Flatten lists with only one element, i.e. if x is a list of length 1, return that single element.
            If x is an empty list, return None.
            """
            if isinstance(x, list):
                if len(x) == 1:
                    return x[0]
                elif len(x) == 0:
                    return None
            return x

        # Check formatting
        region_ids = cast_to_list(region_ids) if region_ids else self.database.keys()
        properties = [properties] if isinstance(properties, str) else properties

        properties_dict = {}

        for region_id in region_ids:
            if region_id in self.database.keys():
                region_db = self.database[region_id]
                for oligo_id, oligo_properties in region_db.items():
                    key = (region_id, oligo_id)
                    if key not in properties_dict:
                        properties_dict[key] = {
                            "region_id": region_id,
                            "oligo_id": oligo_id,
                        }

                    # Retrieve each requested property
                    for property in properties:
                        if property not in oligo_properties:
                            val = None
                        elif flatten:
                            val = flatten_property_list(oligo_properties[property])
                        else:
                            val = oligo_properties[property]
                        properties_dict[key][property] = val

        # Convert the dict-of-dicts to a DataFrame
        properties_table = pd.DataFrame.from_dict(properties_dict, orient="index")
        properties_table = properties_table.reset_index(drop=True)

        # If all lists in the dataframe only contain one element, flatten all list entries
        properties_table = properties_table.map(_flatten_if_one)

        return properties_table

    def get_oligo_property_value(
        self, property: str, flatten: bool, region_id: str, oligo_id: str
    ) -> Any | list[Any] | None:
        """
        Retrieve the value of a specified property for a given oligo and region ID,
        optionally flattening nested or list properties to a unique set of values.

        :param property: The name of the property to retrieve.
        :type property: str
        :param flatten: Whether to flatten list properties to a unique set of values in the table.
        :type flatten: bool
        :param region_id: The ID of the region where the oligo is located.
        :type region_id: str
        :param oligo_id: The ID of the oligo for which the property value is retrieved.
        :type oligo_id: str
        :return: The value of the specified property, possibly flattened. Can be any type (int, float, str, bool, etc.) or a list of values, or None if the property doesn't exist.
        :rtype: Any | list[Any] | None
        :raises ValueError: If the specified region or oligo does not exist in the database.
        """
        if not region_id in self.database:
            raise DatabaseError(f"Region '{region_id}' does not exist in the database.")

        if not oligo_id in self.database[region_id]:
            raise DatabaseError(f"Oligo '{oligo_id}' does not exist in region '{region_id}'.")

        oligo_properties = self.database[region_id][oligo_id]
        if property not in oligo_properties:
            property_value = None
        elif flatten:
            property_value = flatten_property_list(self.database[region_id][oligo_id][property])
            if property_value and len(property_value) == 1:
                property_value = property_value[0]
        else:
            property_value = self.database[region_id][oligo_id][property]

        return property_value

    ############################################
    # Manipulation Functions
    ############################################

    def set_database_sequence_types(self, sequence_types: str | list[str]) -> None:
        """
        Sets or adds sequence types to the database metadata.

        :param sequence_types: Sequence type(s) to add to the database metadata. Can be a single sequence type (str) or a list of sequence types (list[str]).
        :type sequence_types: str | list[str]
        """
        sequence_types_list = cast_to_list(sequence_types)
        for seq_type in sequence_types_list:
            if seq_type not in self.database_sequence_types:
                self.database_sequence_types.append(seq_type)

    def update_oligo_properties(
        self, new_oligo_property: dict, region_ids: str | list[str] | None = None
    ) -> None:
        """
        Updates the properties and properties of oligos in the database with the values provided in `new_oligo_property`.
        This function iterates over the current database and updates each oligo's properties if a corresponding
        entry exists in the `new_oligo_property` dictionary. The update is done in place, modifying the existing
        properties of the oligonucleotides.

        Input format of new property:
        new_property = {<oligo_id>: {<property_name>: <property_value>}}

        Example:
        new_property = {"AARS1::1": {"length": 110}}

        :param new_oligo_property: A dictionary containing new properties for oligos, where keys are oligo IDs
                                    and values are the properties to be updated.
        :type new_oligo_property: dict
        :param region_ids: If provided, only oligos in these regions are updated. Can be a single region ID (str)
                           or a list of region IDs (list[str]). Otherwise all regions are scanned.
        :type region_ids: str | list[str] | None
        """
        region_ids = cast_to_list(region_ids) if region_ids is not None else list(self.database.keys())
        for region_id in region_ids:
            if region_id not in self.database:
                continue
            database_region = self.database[region_id]
            for oligo_id, oligo_properties in database_region.items():
                if oligo_id in new_oligo_property:
                    oligo_properties.update(
                        format_oligo_properties(new_oligo_property[oligo_id], self.database_sequence_types)
                    )

    def filter_database_by_region(self, remove_region: bool, region_ids: str | list[str]) -> None:
        """
        Filters the OligoDatabase based on the specified region IDs. Depending on the `remove_region` flag,
        this function either removes the specified regions from the database or retains only the specified regions
        and removes the others.

        :param remove_region: Flag indicating whether to remove the specified regions (True) or keep only the specified regions (False).
        :type remove_region: bool
        :param region_ids: Region identifier(s) to process. Can be a single region ID (str) or a list of region IDs (list[str]). If None, all regions in the database are processed.
        :type region_ids: str | list[str]
        """
        # Check formatting
        region_ids = cast_to_list(region_ids)
        if self.database:
            for region_id in self.database.keys():
                if remove_region and (region_id in region_ids):
                    del self.database[region_id]
                elif not remove_region and (region_id not in region_ids):
                    del self.database[region_id]
        else:
            raise DatabaseError(
                "Cannot filter database: database is empty. Call load_database() or load_database_from_fasta() first."
            )

    def filter_database_by_oligo(self, remove_region: bool, oligo_ids: str | list[str]) -> None:
        """
        Filters the OligoDatabase based on the specified oligo IDs. Depending on the `remove_region` flag,
        this function either removes the specified oligos from the database or retains only the specified oligos
        and removes the others.

        :param remove_region: Flag indicating whether to remove the specified oligos (True) or keep only the specified oligos (False).
        :type remove_region: bool
        :param oligo_ids: A single oligo ID or a list of oligo IDs to be removed or retained in the database.
        :type oligo_ids: str | list[str]
        """
        # Check formatting
        oligo_ids = cast_to_list(oligo_ids)
        if self.database:
            for region_id in self.database.keys():
                oligo_ids_region = list(self.database[region_id].keys())
                for oligo_id in oligo_ids_region:
                    if remove_region and (oligo_id in oligo_ids):
                        del self.database[region_id][oligo_id]
                    elif not remove_region and (oligo_id not in oligo_ids):
                        del self.database[region_id][oligo_id]
        else:
            raise DatabaseError(
                "Cannot filter database: database is empty. Call load_database() or load_database_from_fasta() first."
            )

    def filter_database_by_property_threshold(
        self, property_name: str, property_thr: float, remove_if_smaller_threshold: bool
    ) -> None:
        """
        Filters the OligoDatabase based on the specified property threshold. The function iterates through all
        oligos in the database and removes those that do not meet the given threshold criteria.

        :param property_name: The name of the property and its threshold to be evaluated.
        :type property_name: str
        :param property_thr: The threshold value for the property.
        :type property_thr: float
        :param remove_if_smaller_threshold: If True, removes oligos with property values smaller than the threshold;
                                            if False, removes oligos with property values larger than the threshold.
        :type remove_if_smaller_threshold: bool
        """
        oligos_to_delete = []
        for region_id in self.database.keys():
            for oligo_id in self.database[region_id].keys():
                property_values = self.get_oligo_property_value(
                    property=property_name, region_id=region_id, oligo_id=oligo_id, flatten=True
                )
                if property_values is not None:
                    property_values = cast_to_list(property_values)
                    if (
                        remove_if_smaller_threshold and any(item < property_thr for item in property_values)
                    ) or (
                        not remove_if_smaller_threshold
                        and all(item > property_thr for item in property_values)
                    ):
                        oligos_to_delete.append((region_id, oligo_id))

        for region_id, oligo_id in oligos_to_delete:
            del self.database[region_id][oligo_id]

    def filter_database_by_property_category(
        self, property_name: str, property_category: str | list[str], remove_if_equals_category: bool
    ) -> None:
        """
        Filters the OligoDatabase by the specified property category. The function removes oligos based on whether
        their property values match or do not match the given category/categories.

        :param property_name: The name of the property to evaluate.
        :type property_name: str
        :param property_category: The category or categories to compare against the property values.
        :type property_category: str | list[str]
        :param remove_if_equals_category: If True, removes oligos with property values that match the specified category;
                                        if False, removes oligos with property values that do not match the specified category.
        :type remove_if_equals_category: bool
        """
        # Check formatting
        property_category = cast_to_list(property_category)
        oligos_to_delete = []

        for region_id in self.database.keys():
            for oligo_id in self.database[region_id].keys():
                property_values = cast_to_list(
                    self.get_oligo_property_value(
                        property=property_name, region_id=region_id, oligo_id=oligo_id, flatten=True
                    )
                )
                if property_values:
                    # remove if any of the items match category
                    if remove_if_equals_category and any(
                        item in property_category for item in property_values
                    ):
                        oligos_to_delete.append((region_id, oligo_id))
                    # remove if all of the items don't match the category
                    elif not remove_if_equals_category and all(
                        item not in property_category for item in property_values
                    ):
                        oligos_to_delete.append((region_id, oligo_id))

        for region_id, oligo_id in oligos_to_delete:
            del self.database[region_id][oligo_id]
