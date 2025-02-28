############################################
# imports
############################################

import os
import pickle
from pathlib import Path
from typing import List, Union, get_args

import pandas as pd
import yaml
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from effidict import LRUPickleDict
from joblib import Parallel, delayed
from joblib_progress import joblib_progress

from oligo_designer_toolsuite._constants import _TYPES_SEQ, SEPARATOR_OLIGO_ID
from oligo_designer_toolsuite.utils import (
    FastaParser,
    CustomYamlDumper,
    check_if_list,
    check_if_region_in_database,
    check_tsv_format,
    collapse_attributes_for_duplicated_sequences,
    format_oligo_attributes,
    merge_databases,
    flatten_attribute_list,
)

CustomYamlDumper.add_representer(list, CustomYamlDumper.represent_list)
CustomYamlDumper.add_representer(dict, CustomYamlDumper.represent_dict)

############################################
# Oligo Database Class
############################################


class OligoDatabase:
    """
    The `OligoDatabase` class provides a comprehensive framework for managing oligonucleotide data along with their associated attributes
    across various genomic regions. It supports loading, saving, filtering, and updating oligo data, ensuring efficient handling of large
    datasets through the use of an LRU (Least Recently Used) cache system. It provides functionalities to
    load data from different sources (such as FASTA files and tables), save data in various formats, and apply
    filtering based on specific criteria like regions, oligo IDs, or attribute thresholds. Additionally, the class
    includes methods to update, retrieve, and analyze oligo attributes, facilitating the selection and evaluation
    of oligonucleotides for research and practical applications.

    :param min_oligos_per_region: Minimum number of oligos required per region to retain the region in the database and oligosets.
    :type min_oligos_per_region: int
    :param write_regions_with_insufficient_oligos: Flag to log regions with insufficient oligos.
    :type write_regions_with_insufficient_oligos: bool
    :param lru_db_max_in_memory: Maximum number of database entries to keep in memory.
    :type lru_db_max_in_memory: int
    :param n_jobs: Number of parallel jobs to use for processing.
    :type n_jobs: int
    :param database_name: Name of the database for storing oligo data.
    :type database_name: str
    :param dir_output: Directory for saving outputs and cached files.
    :type dir_output: str
    """

    def __init__(
        self,
        min_oligos_per_region: int = 0,
        write_regions_with_insufficient_oligos: bool = True,
        lru_db_max_in_memory: int = 10,
        n_jobs: int = 1,
        database_name: str = "db_oligo",
        dir_output: str = "output",
    ) -> None:
        """Constructor for the OligoDatabase class."""

        self.min_oligos_per_region = min_oligos_per_region
        self.write_regions_with_insufficient_oligos = write_regions_with_insufficient_oligos
        self.lru_db_max_in_memory = lru_db_max_in_memory
        self.n_jobs = n_jobs

        self.database_name = database_name
        self.dir_output = os.path.abspath(os.path.join(dir_output, database_name))
        Path(self.dir_output).mkdir(parents=True, exist_ok=True)

        self._dir_cache_files = os.path.join(self.dir_output, "cache_files")

        self.fasta_parser = FastaParser()

        # Initialize databse object
        self.database = LRUPickleDict(
            max_in_memory=self.lru_db_max_in_memory,
            storage_path=self._dir_cache_files,
        )

        self.oligosets = LRUPickleDict(
            max_in_memory=self.lru_db_max_in_memory,
            storage_path=self._dir_cache_files,
        )  # will be used later in the gereration of oligo sets

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
        files_fasta: Union[str, List[str]],
        database_overwrite: bool,
        sequence_type: _TYPES_SEQ,
        region_ids: Union[str, List[str]] = None,
    ) -> None:
        """
        Loads oligonucleotide data from one or more FASTA files into the database, optionally overwriting the existing database.

        This function reads sequences from FASTA file(s) and adds them to the OligoDatabase, either as 'oligo' or
        'target' sequence. It parses the headers of the FASTA entries to extract oligo attribute information,
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

        :param files_fasta: Path(s) to the FASTA file(s) to load.
        :type files_fasta: Union[str, List[str]]
        :param database_overwrite: If True, the existing database will be overwritten.
        :type database_overwrite: bool
        :param sequence_type: The type of sequence being loaded, must be one of the predefined sequence types, i.e. "oligo" or "target".
        :type sequence_type: _TYPES_SEQ["oligo", "target"]
        :param region_ids: List of region IDs to load. If None, all regions in the FASTA files are loaded, defaults to None.
        :type region_ids: Union[str, List[str]], optional
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

                sequences = {}
                for entry in fasta_sequences:
                    region, additional_info, coordinates = self.fasta_parser.parse_fasta_header(entry.id)
                    oligo_attributes = coordinates | additional_info
                    oligo_attributes = format_oligo_attributes(oligo_attributes)
                    if region in sequences:
                        if entry.seq in sequences[region]:
                            oligo_attributes = collapse_attributes_for_duplicated_sequences(
                                oligo_attributes1=sequences[region][entry.seq],
                                oligo_attributes2=oligo_attributes,
                            )
                        sequences[region][str(entry.seq)] = oligo_attributes
                    else:
                        sequences[region] = {str(entry.seq): oligo_attributes}

                database_region = {region: {} for region in sequences.keys()}
                for region, sequences_region in sequences.items():
                    i = 1
                    for oligo_sequence, oligo_attributes in sequences_region.items():
                        oligo_id = f"{region}{SEPARATOR_OLIGO_ID}{i}"
                        oligo_seq_info = {sequence_type: oligo_sequence} | oligo_attributes
                        database_region[region][oligo_id] = oligo_seq_info
                        i += 1

                # only merge if there are common keys
                if len(set(self.database) & set(database_region)) > 0:
                    self.database = merge_databases(
                        database1=self.database,
                        database2=database_region,
                        dir_cache_files=self._dir_cache_files,
                        lru_db_max_in_memory=self.lru_db_max_in_memory,
                    )
                else:
                    for region in database_region.keys():
                        self.database[region] = database_region[region]

        # Check formatting
        region_ids = check_if_list(region_ids)
        files_fasta = check_if_list(files_fasta)

        # Check if sequence type is correct
        options = get_args(_TYPES_SEQ)
        assert (
            sequence_type in options
        ), f"Sequence type not supported! '{sequence_type}' is not in {options}."

        # Clear database if it should be overwritten
        if database_overwrite:
            self.database = LRUPickleDict(
                max_in_memory=self.lru_db_max_in_memory,
                storage_path=self._dir_cache_files,
            )

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
        region_ids: Union[str, List[str]] = None,
    ) -> None:
        """
        Loads oligonucleotide data from a tab-delimited table (TSV) file into the database, optionally overwriting the existing database.

        This function loads the OligoDatabase from a tab-separated values (TSV) file. The file must contain
        columns such as 'region_id', 'oligo_id' and optionally the oligo or target sequence as well as oligo attributes.
        The database can be optionally filtered by specifying a list of region IDs.

        ⚠️ For big databases, it is not recommended to load the whole TSV file at once. Instead, the database should be
        split into smaller files. The function will merge the databases if there are common keys (regions).

        :param file_database: Path to the TSV file containing the database.
        :type file_database: str
        :param database_overwrite: If True, the existing database will be overwritten.
        :type database_overwrite: bool
        :param region_ids: List of region IDs to load. If None, all regions in the FASTA files are loaded, defaults to None.
        :type region_ids: Union[str, List[str]], optional
        """
        # Check formatting
        region_ids = check_if_list(region_ids)

        # Check if file exists and has correct format
        if os.path.exists(file_database):
            if not check_tsv_format(file_database):
                raise ValueError("Database has incorrect format!")
        else:
            raise ValueError("Database file does not exist!")

        # Clear database if it should be overwritten
        if database_overwrite:
            self.database = LRUPickleDict(
                max_in_memory=self.lru_db_max_in_memory,
                storage_path=self._dir_cache_files,
            )

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
        database_tmp2 = LRUPickleDict(
            max_in_memory=self.lru_db_max_in_memory,
            storage_path=self._dir_cache_files,
        )
        for entry in database_tmp1:
            region_id, oligo_id = entry.pop("region_id"), entry.pop("oligo_id")
            if (not region_ids) or (region_ids and region_id in region_ids):
                if region_id not in database_tmp2:
                    database_tmp2[region_id] = {}
                database_tmp2[region_id][oligo_id] = format_oligo_attributes(entry)

        if not database_overwrite and self.database:
            database_tmp2 = merge_databases(
                database1=self.database,
                database2=database_tmp2,
                dir_cache_files=self._dir_cache_files,
                lru_db_max_in_memory=self.lru_db_max_in_memory,
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
        region_ids: Union[str, List[str]] = None,
    ) -> None:
        """
        Loads oligonucleotide data from a directory of pickled files into the database, optionally overwriting the existing database.
        Each file in the folder represents a region in the database and contains the region ID,
        the oligo sequence and attribute information as well as the oligosets for this region.

        :param dir_database: Path to the directory containing the pickled database files.
        :type dir_database: str
        :param database_overwrite: If True, the existing database will be overwritten.
        :type database_overwrite: bool
        :param region_ids: List of region IDs to load. If None, all regions in the FASTA files are loaded, defaults to None.
        :type region_ids: Union[str, List[str]], optional
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
            # extract region ID from the file name and remove the extension
            with open(file, "rb") as handle:
                content = pickle.load(handle)
                if (not region_ids) or (region_ids and region_id in region_ids):
                    region_id = content["region_id"]
                    database_region = content["database_region"]
                    oligoset_region = content["oligoset_region"]

            # only merge if there are common keys
            if region_id in self.database.keys():
                self.database = merge_databases(
                    database1=self.database,
                    database2={region_id: database_region},
                    dir_cache_files=self._dir_cache_files,
                    lru_db_max_in_memory=self.lru_db_max_in_memory,
                )
                self.oligosets[region_id] = pd.concat([self.oligosets[region_id], oligoset_region])
            else:
                self.database[region_id] = database_region
                self.oligosets[region_id] = oligoset_region

        # Check formatting
        region_ids = check_if_list(region_ids)

        if not os.path.isdir(dir_database):
            raise ValueError("Database directory does not exist!")

        if database_overwrite:
            self.database = LRUPickleDict(
                max_in_memory=self.lru_db_max_in_memory,
                storage_path=self._dir_cache_files,
            )

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
        dir_database: str = "db_oligo",
        region_ids: Union[str, List[str]] = None,
    ) -> str:
        """
        Saves the current database and associated oligosets as pickl files into a specified directory.

        :param dir_database: Directory path where the database files should be saved. Default is "db_oligo".
        :type dir_database: str
        :param region_ids: List of region IDs to save. If None, all regions in the database are saved, defaults to None.
        :type region_ids: Union[str, List[str]], optional
        :return: The directory path where the database was saved.
        :rtype: str
        """
        # Check formatting
        region_ids = check_if_list(region_ids) if region_ids else self.database.keys()

        dir_database = os.path.join(self.dir_output, dir_database)
        Path(dir_database).mkdir(parents=True, exist_ok=True)

        for region_id in region_ids:
            database_region = self.database[region_id]
            oligoset_region = self.oligosets[region_id]
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
        sequence_type: _TYPES_SEQ,
        save_description: bool,
        filename: str = "db_oligo",
        region_ids: Union[str, List[str]] = None,
    ) -> str:
        """
        Writes the current database to a FASTA file. Associated sequence attributes can optionally be included in the sequence header.

        :param sequence_type: The type of sequence being written, must be one of the predefined sequence types, i.e. "oligo" or "target".
        :type sequence_type: _TYPES_SEQ["oligo", "target"]
        :param save_description: Whether to include the sequence attributes in the sequence header.
        :type save_description: bool
        :param filename: The base name of the output FASTA file, defaults to "db_oligo".
        :type filename: str
        :param region_ids: List of region IDs to write. If None, all regions in the database are written, defaults to None.
        :type region_ids: Union[str, List[str]], optional
        :return: The path to the saved FASTA file.
        :rtype: str
        """
        # Check if sequence type is correct
        options = get_args(_TYPES_SEQ)
        assert (
            sequence_type in options
        ), f"Sequence type not supported! '{sequence_type}' is not in {options}."

        # Check formatting
        region_ids = check_if_list(region_ids) if region_ids else self.database.keys()

        file_fasta = os.path.join(self.dir_output, f"{filename}.fna")
        output_fasta = []

        with open(file_fasta, "w") as handle_fasta:
            for region_id in region_ids:
                database_region = self.database[region_id]
                for oligo_id, oligo_attributes in database_region.items():
                    description = sequence_type if save_description else ""
                    seq_record = SeqRecord(
                        Seq(oligo_attributes[sequence_type]),
                        id=oligo_id,
                        name=oligo_id.split(SEPARATOR_OLIGO_ID)[0],
                        description=description,
                    )
                    output_fasta.append(seq_record)

            SeqIO.write(output_fasta, handle_fasta, "fasta")

        return file_fasta

    def write_database_to_table(
        self,
        attributes: Union[str, List[str]],
        flatten_attribute: bool,
        filename: str = "oligo_database_table",
        region_ids: list[str] = None,
    ) -> str:
        """
        Writes the current database and selected attributes to a TSV table file,
        optionally flattening nested or list attributes to a unique set of values.


        :param attributes: A list of attributes to include in the table.
        :type attributes: Union[str, List[str]]
        :param flatten_attribute: Whether to flatten list attributes to a unique set of values in the table.
        :type flatten_attribute: bool
        :param filename: The base name of the output TSV file, defaults to "oligo_database_table".
        :type filename: str
        :param region_ids: List of region IDs to write. If None, all regions in the database are written, defaults to None.
        :type region_ids: Union[str, List[str]], optional
        :return: The path to the saved TSV file.
        :rtype: str
        """
        # Check formatting
        region_ids = check_if_list(region_ids) if region_ids else self.database.keys()
        attributes = check_if_list(attributes)

        file_table = os.path.join(os.path.dirname(self.dir_output), f"{filename}.tsv")

        first_entry = True
        for region_id in region_ids:
            file_tsv_content = []
            for oligo_id in self.database[region_id].keys():
                entry = {"region_id": region_id, "oligo_id": oligo_id}
                for attribute in attributes:
                    if attribute in self.database[region_id][oligo_id]:
                        oligo_attribute = self.database[region_id][oligo_id][attribute]
                        if flatten_attribute:
                            oligo_attribute = flatten_attribute_list(oligo_attribute)
                            oligo_attribute = list(
                                set(
                                    oligo_attribute
                                    if isinstance(oligo_attribute, list)
                                    else [oligo_attribute]
                                )
                            )
                            entry[attribute] = (
                                str(oligo_attribute).replace("'", "").replace("[", "").replace("]", "")
                            )
                        else:
                            entry[attribute] = str(oligo_attribute).replace("[[", "[").replace("]]", "]")
                file_tsv_content.append(entry)
            file_tsv_content = pd.DataFrame(data=file_tsv_content)
            file_tsv_content.to_csv(file_table, sep="\t", index=False, mode="a", header=first_entry)
            first_entry = False

        return file_table

    def write_oligosets_to_yaml(
        self,
        attributes: Union[str, List[str]],
        top_n_sets: int,
        ascending: bool,
        filename: str = "oligosets",
        region_ids: list[str] = None,
    ) -> str:
        """
        Writes the current top n oligosets and selected attributes to a YAML file.
        The oligosets are sorted based on their scores in ascending or descending order.

        :param attributes: A list of attributes to include in the YAML file.
        :type attributes: Union[str, List[str]]
        :param top_n_sets: Number of top oligosets to include based on their score.
        :type top_n_sets: int
        :param ascending: If True, sort oligosets by score in ascending order (the smaller the score the better the oligo set);
            otherwise, in descending order (the higher the score the better the oligo set).
        :type ascending: bool
        :param filename: Base name for the output YAML file, defaults to "oligosets".
        :type filename: str
        :param region_ids: List of region IDs to write. If None, all regions in the database are written, defaults to None.
        :type region_ids: Union[str, List[str]], optional
        :return: Path to the saved YAML file.
        :rtype: str
        """
        # Check formatting
        attributes = check_if_list(attributes)
        region_ids = check_if_list(region_ids) if region_ids else self.database.keys()

        yaml_dict = {region_id: {} for region_id in region_ids}

        for region_id in region_ids:
            oligosets_region = self.oligosets[region_id]
            oligosets_oligo_columns = [col for col in oligosets_region.columns if col.startswith("oligo_")]
            oligosets_score_columns = [
                col for col in oligosets_region.columns if col.startswith("set_score_")
            ]

            oligosets_region = oligosets_region.sort_values(by=oligosets_score_columns, ascending=ascending)
            oligosets_region_oligos = oligosets_region.head(top_n_sets)[oligosets_oligo_columns]
            oligosets_region_scores = oligosets_region.head(top_n_sets)[oligosets_score_columns]

            for idx, oligoset in oligosets_region_oligos.iterrows():
                oligoset_id = f"Oligoset {idx + 1}"
                yaml_dict[region_id][oligoset_id] = {
                    "Oligoset Score": oligosets_region_scores.iloc[idx].to_dict(),
                }

                for oligo_idx, oligo_id in enumerate(oligoset):
                    yaml_dict_oligo_entry = {"oligo_id": oligo_id}

                    # iterate through all attributes that should be written
                    for attribute in attributes:
                        if attribute in self.database[region_id][oligo_id]:
                            oligo_attribute = self.database[region_id][oligo_id][attribute]
                            yaml_dict_oligo_entry[attribute] = oligo_attribute

                    oligo_id_yaml = f"Oligo {oligo_idx + 1}"
                    yaml_dict[region_id][oligoset_id][oligo_id_yaml] = yaml_dict_oligo_entry

        file_yaml = os.path.join(os.path.dirname(self.dir_output), f"{filename}.yml")

        with open(file_yaml, "w") as handle:
            yaml.dump(yaml_dict, handle, Dumper=CustomYamlDumper, default_flow_style=False, sort_keys=False)

        return file_yaml

    def write_oligosets_to_table(
        self, dir_output: str = "oligosets", region_ids: Union[str, List[str]] = None
    ) -> str:
        """
        Writes the current oligosets for each specified region to a tab-separated value (TSV) file, saved in the specified folder.

        :param dir_output: Name of the output folder where the oligoset TSV files will be saved, defaults to "oligosets".
        :type dir_output: str
        :param region_ids: List of region IDs to write. If None, all regions in the database are written, defaults to None.
        :type region_ids: Union[str, List[str]], optional
        :return: Path to the directory where the oligoset TSV files were saved.
        :rtype: str
        """
        # Check formatting
        region_ids = check_if_list(region_ids) if region_ids else self.database.keys()

        dir_oligosets = os.path.join(os.path.dirname(self.dir_output), dir_output)
        Path(dir_oligosets).mkdir(parents=True, exist_ok=True)

        for region_id in self.oligosets.keys():
            if region_id in region_ids:
                file_oligosets = os.path.join(dir_oligosets, f"oligosets_{region_id}.tsv")
                self.oligosets[region_id].to_csv(file_oligosets, sep="\t", index=False)

        return dir_oligosets

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
            if len(self.database[region_id]) <= self.min_oligos_per_region
        ]

        for region in regions_to_remove:
            # TODO: this is a workaround due to a bug fix in EffiDict which needs to be fixed
            self.database[region] = None
            del self.database[region]
            del self.database[region]

            self.oligosets[region] = None
            del self.oligosets[region]
            del self.oligosets[region]

        if self.write_regions_with_insufficient_oligos and regions_to_remove:
            with open(self.file_removed_regions, "a") as handle:
                handle.write("\n".join(f"{region}\t{pipeline_step}" for region in regions_to_remove) + "\n")

    ############################################
    # Getter Functions
    ############################################

    def get_attribute_list(self) -> list[str]:
        """Retrieves a list of attribute names stored in the database.

        :return: A list of attribute names stored in the database.
        :rtype: list[str]
        """
        region_id = next(iter(self.database.values()))
        oligo_id = next(iter(region_id.values()))
        attributes = list(oligo_id.keys())

        return attributes

    def get_regionid_list(self) -> list[str]:
        """
        Retrieves a list of all region IDs present in the database.

        :return: A list of region IDs in the database.
        :rtype: list[str]
        """
        region_ids = list(self.database.keys())

        return region_ids

    def get_oligoid_list(self) -> list[str]:
        """
        Retrieves a list of all oligo IDs present in the database.

        :return: A list of oligo IDs from all regions in the database.
        :rtype: list[str]
        """
        oligo_ids = [
            oligo_id for database_region in self.database.values() for oligo_id in database_region.keys()
        ]

        return oligo_ids

    def get_sequence_list(self, sequence_type: _TYPES_SEQ) -> list[str]:
        """
        Retrieves a list of sequences of the specified type from the database.

        :param sequence_type: The type of sequence to retrieve, must be one of the predefined sequence types, i.e. "oligo" or "target".
        :type sequence_type: _TYPES_SEQ["oligo", "target"]
        :return: A list of sequences corresponding to the specified sequence type from all regions in the database.
        :rtype: list[str]
        """
        options = get_args(_TYPES_SEQ)
        assert (
            sequence_type in options
        ), f"Sequence type not supported! '{sequence_type}' is not in {options}."
        sequences = [
            str(oligo_attributes[sequence_type])
            for region_id, database_region in self.database.items()
            for oligo_id, oligo_attributes in database_region.items()
        ]

        return sequences

    def get_oligoid_sequence_mapping(
        self, sequence_type: _TYPES_SEQ, sequence_to_upper: bool = False
    ) -> dict:
        """
        Generates a mapping of oligo IDs to their corresponding sequences of the specified type.

        :param sequence_type: The type of sequence to map, must be one of the predefined sequence types, i.e. "oligo" or "target".
        :type sequence_type: _TYPES_SEQ["oligo", "target"]
        :param sequence_to_upper: Whether to convert sequences to uppercase, defaults to False.
        :type sequence_to_upper: bool
        :return: A dictionary mapping oligo IDs to their corresponding sequences.
        :rtype: dict
        """
        options = get_args(_TYPES_SEQ)
        assert (
            sequence_type in options
        ), f"Sequence type not supported! '{sequence_type}' is not in {options}."

        oligoid_sequence_mapping = {}

        for region_id, database_region in self.database.items():
            for oligo_id, oligo_attributes in database_region.items():
                seq = oligo_attributes[sequence_type]
                if sequence_to_upper:
                    seq = seq.upper()
                oligoid_sequence_mapping[oligo_id] = seq

        return oligoid_sequence_mapping

    def get_sequence_oligoid_mapping(
        self, sequence_type: _TYPES_SEQ, sequence_to_upper: bool = False
    ) -> dict:
        """
        Generates a mapping of sequences to their corresponding oligo IDs for the specified sequence type.

        :param sequence_type: The type of sequence to map, must be one of the predefined sequence types, i.e. "oligo" or "target".
        :type sequence_type: _TYPES_SEQ["oligo", "target"]
        :param sequence_to_upper: Whether to convert sequences to uppercase, defaults to False.
        :type sequence_to_upper: bool
        :return: A dictionary mapping sequences to their corresponding oligo IDs.
        :rtype: dict
        """
        options = get_args(_TYPES_SEQ)
        assert (
            sequence_type in options
        ), f"Sequence type not supported! '{sequence_type}' is not in {options}."

        sequence_oligoids_mapping = {}

        for region_id, database_region in self.database.items():
            for oligo_id, oligo_attributes in database_region.items():
                seq = oligo_attributes[sequence_type]
                if sequence_to_upper:
                    seq = seq.upper()
                if seq not in sequence_oligoids_mapping:
                    # If the sequence key doesn't exist, create a new entry with the oligo ID in a list
                    sequence_oligoids_mapping[seq] = [oligo_id]
                else:
                    # If the sequence key already exists, append the oligo ID to the existing list
                    sequence_oligoids_mapping[seq].append(oligo_id)

        return sequence_oligoids_mapping

    def get_oligo_attribute_table(
        self, attribute: str, flatten: bool, region_ids: Union[str, List[str]] = None
    ) -> pd.DataFrame:
        """
        Generates a DataFrame containing oligo IDs and the specified attribute for each oligo,
        optionally flattening nested or list attributes to a unique set of values.

        :param attribute: The name of the attribute to retrieve.
        :type attribute: str
        :param flatten: Whether to flatten list attributes to a unique set of values in the table.
        :type flatten: bool
        :param region_ids: List of region IDs to retrieve. If None, all regions in the database are retrieved, defaults to None.
        :type region_ids: Union[str, List[str]], optional
        :return: A DataFrame with oligo IDs and the corresponding attribute values.
        :rtype: pd.DataFrame
        """
        # Check formatting
        region_ids = check_if_list(region_ids) if region_ids else self.database.keys()

        oligo_ids = []
        attributes = []

        for region_id in region_ids:
            database_region = self.database[region_id]
            for oligo_id, oligo_attributes in database_region.items():
                oligo_ids.append(oligo_id)
                if attribute not in oligo_attributes:
                    attributes.append(None)
                elif flatten:
                    attributes.append(flatten_attribute_list(oligo_attributes[attribute]))
                else:
                    attributes.append(oligo_attributes[attribute])

        return pd.DataFrame({"oligo_id": oligo_ids, attribute: attributes})

    def get_oligo_attribute_value(
        self, attribute: str, flatten: bool, region_id: str, oligo_id: str
    ) -> Union[str, List[str]]:
        """
        Retrieve the value of a specified attribute for a given oligo and region ID,
        optionally flattening nested or list attributes to a unique set of values.

        :param attribute: The name of the attribute to retrieve.
        :type attribute: str
        :param flatten: Whether to flatten list attributes to a unique set of values in the table.
        :type flatten: bool
        :param region_id: The ID of the region where the oligo is located.
        :type region_id: str
        :param oligo_id: The ID of the oligo for which the attribute value is retrieved.
        :type oligo_id: str
        :return: The value of the specified attribute, possibly flattened.
        :rtype: Union[str, List[str]]
        :raises ValueError: If the specified region or oligo does not exist.
        """
        if not region_id in self.database:
            raise ValueError(f"Region {region_id} does not exist.")

        if not oligo_id in self.database[region_id]:
            raise ValueError(f"Region {oligo_id} does not exist.")

        oligo_attributes = self.database[region_id][oligo_id]
        if attribute not in oligo_attributes:
            attribute_value = None
        elif flatten:
            attribute_value = flatten_attribute_list(self.database[region_id][oligo_id][attribute])
        else:
            attribute_value = self.database[region_id][oligo_id][attribute]

        return attribute_value

    ############################################
    # Manipulation Functions
    ############################################

    def update_oligo_attributes(self, new_oligo_attribute: dict) -> None:
        """
        Updates the attributes of oligos in the database with the values provided in `new_oligo_attribute`.
        This function iterates over the current database and updates each oligo's attributes if a corresponding
        entry exists in the `new_oligo_attribute` dictionary. The update is done in place, modifying the existing
        attributes of the oligonucleotides.

        Input format of new attribute:
        new_attribute = {<oligo_id>: {<attribute_name>: <attribute_value>}}

        Example:
        new_attribute = {"AARS1::1": {"length": 110}}

        :param new_oligo_attribute: A dictionary containing new attributes for oligos, where keys are oligo IDs
                                    and values are the attributes to be updated.
        :type new_oligo_attribute: dict
        """
        for region_id, database_region in self.database.items():
            for oligo_id, oligo_attributes in database_region.items():
                if oligo_id in new_oligo_attribute:
                    oligo_attributes.update(format_oligo_attributes(new_oligo_attribute[oligo_id]))

    def filter_database_by_region(self, remove_region: bool, region_ids: Union[str, List[str]]) -> None:
        """
        Filters the OligoDatabase based on the specified region IDs. Depending on the `remove_region` flag,
        this function either removes the specified regions from the database or retains only the specified regions
        and removes the others.

        :param remove_region: Flag indicating whether to remove the specified regions (True) or keep only the specified regions (False).
        :type remove_region: bool
        :param region_ids: A single region ID or a list of region IDs to be removed or retained in the database.
        :type region_ids: Union[str, List[str]]
        """
        # Check formatting
        region_ids = check_if_list(region_ids)
        if self.database:
            for region_id in self.database.keys():
                if remove_region and (region_id in region_ids):
                    del self.database[region_id]
                elif not remove_region and (region_id not in region_ids):
                    del self.database[region_id]
        else:
            raise ValueError(
                "Can not filter. Database is empty! Call the method load_database() or load_database_from_fasta() first."
            )

    def filter_database_by_oligo(self, remove_region: bool, oligo_ids: Union[str, List[str]]) -> None:
        """
        Filters the OligoDatabase based on the specified oligo IDs. Depending on the `remove_region` flag,
        this function either removes the specified oligos from the database or retains only the specified oligos
        and removes the others.

        :param remove_region: Flag indicating whether to remove the specified oligos (True) or keep only the specified oligos (False).
        :type remove_region: bool
        :param oligo_ids: A single oligo ID or a list of oligo IDs to be removed or retained in the database.
        :type oligo_ids: Union[str, List[str]]
        """
        # Check formatting
        oligo_ids = check_if_list(oligo_ids)
        if self.database:
            for region_id in self.database.keys():
                oligo_ids_region = list(self.database[region_id].keys())
                for oligo_id in oligo_ids_region:
                    if remove_region and (oligo_id in oligo_ids):
                        del self.database[region_id][oligo_id]
                    elif not remove_region and (oligo_id not in oligo_ids):
                        del self.database[region_id][oligo_id]
        else:
            raise ValueError(
                "Can not filter. Database is empty! Call the method load_database() or load_database_from_fasta() first."
            )

    def filter_database_by_attribute_threshold(
        self, attribute_name: str, attribute_thr: float, remove_if_smaller_threshold: bool
    ) -> None:
        """
        Filters the OligoDatabase based on the specified attribute threshold. The function iterates through all
        oligos in the database and removes those that do not meet the given threshold criteria.

        :param attribute_name: The name of the attribute to be evaluated.
        :type attribute_name: str
        :param attribute_thr: The threshold value for the attribute.
        :type attribute_thr: float
        :param remove_if_smaller_threshold: If True, removes oligos with attribute values smaller than the threshold;
                                            if False, removes oligos with attribute values larger than the threshold.
        :type remove_if_smaller_threshold: bool
        """
        oligos_to_delete = []
        for region_id in self.database.keys():
            for oligo_id in self.database[region_id].keys():
                attribute_values = check_if_list(
                    self.get_oligo_attribute_value(
                        attribute=attribute_name, region_id=region_id, oligo_id=oligo_id, flatten=True
                    )
                )
                if attribute_values:
                    if (
                        remove_if_smaller_threshold and any(item < attribute_thr for item in attribute_values)
                    ) or (
                        not remove_if_smaller_threshold
                        and all(item > attribute_thr for item in attribute_values)
                    ):
                        oligos_to_delete.append((region_id, oligo_id))

        for region_id, oligo_id in oligos_to_delete:
            del self.database[region_id][oligo_id]

    def filter_database_by_attribute_category(
        self, attribute_name: str, attribute_category: Union[str, List[str]], remove_if_equals_category: bool
    ) -> None:
        """
        Filters the OligoDatabase by the specified attribute category. The function removes oligos based on whether
        their attribute values match or do not match the given category/categories.

        :param attribute_name: The name of the attribute to evaluate.
        :type attribute_name: str
        :param attribute_category: The category or categories to compare against the attribute values.
        :type attribute_category: Union[str, List[str]]
        :param remove_if_equals_category: If True, removes oligos with attribute values that match the specified category;
                                        if False, removes oligos with attribute values that do not match the specified category.
        :type remove_if_equals_category: bool
        """
        # Check formatting
        attribute_category = check_if_list(attribute_category)
        oligos_to_delete = []

        for region_id in self.database.keys():
            for oligo_id in self.database[region_id].keys():
                attribute_values = check_if_list(
                    self.get_oligo_attribute_value(
                        attribute=attribute_name, region_id=region_id, oligo_id=oligo_id, flatten=True
                    )
                )
                if attribute_values:
                    # remove if any of the items match category
                    if remove_if_equals_category and any(
                        item in attribute_category for item in attribute_values
                    ):
                        oligos_to_delete.append((region_id, oligo_id))
                    # remove if all of the items don't match the category
                    elif not remove_if_equals_category and all(
                        item not in attribute_category for item in attribute_values
                    ):
                        oligos_to_delete.append((region_id, oligo_id))

        for region_id, oligo_id in oligos_to_delete:
            del self.database[region_id][oligo_id]
