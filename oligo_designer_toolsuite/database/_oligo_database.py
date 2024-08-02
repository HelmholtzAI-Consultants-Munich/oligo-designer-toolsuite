############################################
# imports
############################################

import os
import pickle
import warnings
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
    check_if_key_exists,
    check_if_list,
    check_tsv_format,
    check_if_region_in_database,
    check_tsv_format,
    collapse_attributes_for_duplicated_sequences,
    filter_dabase_for_region,
    format_oligo_attributes,
    merge_databases,
    filter_dabase_for_region,
)

############################################
# Oligo Database Class
############################################


class OligoDatabase:
    """Class for managing the oligo databases. This class provides functionality for handling oligo databases,
    allowing users to load, manipulate, and save oligo information efficiently. It includes methods for
    loading and managing oligo databases from various sources (e.g. fasta file or saved database),
    and performing operations such as removing regions with insufficient oligos and calculating the number of targeted
    transcripts for each oligo.

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

    :param min_oligos_per_region: Minimum number of oligos required per region (default is 0).
    :type min_oligos_per_region: int, optional
    :param write_regions_with_insufficient_oligos: Flag to enable writing regions with insufficient oligos to a file (default is True).
    :type write_regions_with_insufficient_oligos: bool, optional
    :param lru_db_max_in_memory: Maximum number of dictionary entries stored in RAM, defaults to 1.
    :type lru_db_max_in_memory: int, optional
    :param database_name: Subdirectory path for the output, i.e. <dir_output>/<database_name>, defaults to "db_oligo".
    :type database_name: str, optional
    :param dir_output: Directory path for the output, defaults to "output".
    :type dir_output: str, optional
    :param n_jobs: The number of parallel jobs to run. Default is 1.
    :type n_jobs: int
    """

    def __init__(
        self,
        min_oligos_per_region: int = 0,
        write_regions_with_insufficient_oligos: bool = True,
        lru_db_max_in_memory: int = 10,
        database_name: str = "db_oligo",
        dir_output: str = "output",
        n_jobs: int = 1,
    ):
        """Constructor for the OligoDatabase class."""
        self.min_oligos_per_region = min_oligos_per_region
        self.write_regions_with_insufficient_oligos = write_regions_with_insufficient_oligos

        self.database_name = database_name
        self.dir_output = os.path.abspath(os.path.join(dir_output, database_name))
        Path(self.dir_output).mkdir(parents=True, exist_ok=True)

        self._dir_cache_files = os.path.join(self.dir_output, "cache_files")
        self.lru_db_max_in_memory = lru_db_max_in_memory

        self.n_jobs = n_jobs

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

    def load_database(
        self,
        dir_database: str,
        region_ids: Union[str, List[str]] = None,
        database_overwrite: bool = False,
    ) -> None:
        """Load a previously saved oligo database from a folder containing pickled files.

        This function loads the oligo database from a folder containing pickled files. Each file in the folder
        represents a region in the database. The file must contain a dictionary with oligo IDs as keys and
        oligo sequences and additional metadata as values.

        The database can be optionally filtered by specifying a list of region IDs.

        :param dir_database: Path to the folder containing pickled files of the oligo database.
        :type dir_database: str
        :param region_ids: List of region IDs to filter the database. Defaults to None (all regions are loaded).
        :type region_ids: Union[str, List[str]], optional
        :param database_overwrite: If True, overwrite the existing database. Defaults to False.
        :type database_overwrite: bool, optional

        :raises ValueError: If the database file has an incorrect format or does not exist.
        """

        def load(file):
            # extract region ID from the file name and remove the extension
            region_id = os.path.basename(file).split(".")[0]
            with open(file, "rb") as handle:
                content = pickle.load(handle)
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
                delayed(load)(file_database) for file_database in files_database
            )

        # add this step to log regions which are not available in database
        if region_ids:
            check_if_region_in_database(
                database=self.database,
                region_ids=region_ids,
                write_regions_with_insufficient_oligos=self.write_regions_with_insufficient_oligos,
                file_removed_regions=self.file_removed_regions,
            )

    def load_database_from_fasta(
        self,
        files_fasta: list[str],
        sequence_type: _TYPES_SEQ,
        region_ids: list[str] = None,
        database_overwrite: bool = False,
    ) -> None:
        """Load "oligo" or "target" sequences from one or more FASTA files into the oligo database.

        This function reads sequences from FASTA file(s) and adds them to the oligo database, either as 'oligo' or
        'target' sequence and computes the reverse compliment sequence for the other sequence type. It parses the
        headers of the FASTA entries to extract region information, and assigns unique IDs to the oligos within
        each region.

        :param files_fasta: Path to the FASTA file(s) containing the sequences.
        :type files_fasta: list[str]
        :param sequence_type: Type of sequence to load, either 'target' or 'oligo'.
        :type sequence_type: _TYPES_SEQ
        :param region_ids: List of region IDs to filter the database. Defaults to None.
        :type region_ids: list[str], optional
        :param database_overwrite: If True, overwrite the existing database. Defaults to False.
        :type database_overwrite: bool, optional

        :raises AssertionError: If the provided sequence type is not supported.
        """

        def load(file):
            """Load sequences from a FASTA file into the oligo database.

            This function reads sequences from a provided FASTA file, processes them to extract relevant information,
            and stores them in the oligo database. If there are common keys between the existing database and the new
            data, the databases are merged.

            :param file: Path to the FASTA file to be loaded.
            :type file: str
            """
            self.fasta_parser.check_fasta_format(file)
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
                    oligo_sequence_reverse_complement = str(Seq(oligo_sequence).reverse_complement())
                    oligo_seq_info = {
                        sequence_type: oligo_sequence,
                        sequence_type_reverse_complement: oligo_sequence_reverse_complement,
                    } | oligo_attributes
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

        # Check if sequence type is correct
        options = get_args(_TYPES_SEQ)
        assert (
            sequence_type in options
        ), f"Sequence type not supported! '{sequence_type}' is not in {options}."
        sequence_type_reverse_complement = options[0] if options[0] != sequence_type else options[1]

        # check if region ids and fasta files are given as lists
        region_ids = check_if_list(region_ids)
        files_fasta = check_if_list(files_fasta)

        # Clear database if it should be overwritten
        if database_overwrite:
            self.database = LRUPickleDict(
                max_in_memory=self.lru_db_max_in_memory,
                storage_path=self._dir_cache_files,
            )

        # Load files parallel into database
        with joblib_progress(description=f"Database Loading", total=len(files_fasta)):
            Parallel(n_jobs=self.n_jobs, prefer="threads", require="sharedmem")(
                delayed(load)(file_fasta) for file_fasta in files_fasta
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
        region_ids: Union[str, List[str]] = None,
        database_overwrite: bool = False,
    ) -> None:
        """Load a previously written oligo database from a TSV file.

        This function loads the oligo database from a tab-separated values (TSV) file. The file must contain
        columns such as 'region_id', 'oligo_id', 'sequence', and additional attributes, like information from the
        fasta headers or information computed by the filtering classes.
        The order of columns in the database file is:

        +-----------+----------+----------+------------+-------+-----+--------+--------+------------------+
        | region_id | oligo_id | sequence | chromosome | start | end | strand | length | additional feat. |
        +-----------+----------+----------+------------+-------+-----+--------+--------+------------------+

        The database can be optionally filtered by specifying a list of region IDs.

        ⚠️ For big databases, it is not recommended to load the whole TSV file at once. Instead, the database should be
        split into smaller files. The function will merge the databases if there are common keys (regions).

        :param file_database: Path to the TSV file containing the oligo database.
        :type file_database: str
        :param region_ids: List of region IDs to filter the database. Defaults to None.
        :type region_ids: Union[str, List[str]], optional
        :param database_overwrite: If True, overwrite the existing database. Defaults to False.
        :type database_overwrite: bool, optional

        :raises ValueError: If the database file has an incorrect format or does not exist.
        """
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
            if region_id not in database_tmp2:
                database_tmp2[region_id] = {}
            database_tmp2[region_id][oligo_id] = entry

        if not database_overwrite and self.database:
            database_tmp2 = merge_databases(
                database1=self.database,
                database2=database_tmp2,
                dir_cache_files=self._dir_cache_files,
                lru_db_max_in_memory=self.lru_db_max_in_memory,
            )

        # Filter for region ids
        region_ids = check_if_list(region_ids)
        if region_ids:
            database_tmp2 = filter_dabase_for_region(database=database_tmp2, region_ids=region_ids)
            check_if_region_in_database(
                database=database_tmp2,
                region_ids=region_ids,
                write_regions_with_insufficient_oligos=self.write_regions_with_insufficient_oligos,
                file_removed_regions=self.file_removed_regions,
            )

        self.database = database_tmp2

    ############################################
    # Save Functions
    ############################################

    def save_database(
        self,
        dir_database: str = "db_oligo",
        region_ids: list[str] = None,
    ):
        """Save the oligo database to files.

        This function saves the oligo database to the specified directory. Each region's oligo data is saved as a
        separate pickle file in the directory. Hence, one file per region is created.

        :param dir_database: Directory name where the database files will be saved, defaults to "db_oligo".
        :type dir_database: str, optional
        :param region_ids: A list of region IDs to save. If None, all regions will be saved.
        :type region_ids: list[str], optional
        :return: The directory where the database files are saved.
        :rtype: str
        """
        if region_ids:
            region_ids = check_if_list(region_ids)
        else:
            region_ids = self.database.keys()

        dir_database = os.path.join(self.dir_output, dir_database)
        Path(dir_database).mkdir(parents=True, exist_ok=True)

        for region_id in region_ids:
            database_region = self.database[region_id]
            oligoset_region = self.oligosets[region_id]
            file_output = os.path.join(dir_database, region_id)
            with open(file_output, "wb") as file:
                pickle.dump({"database_region": database_region, "oligoset_region": oligoset_region}, file)

        return dir_database

    def write_database_to_fasta(
        self,
        filename: str = "db_oligo",
        region_ids: list[str] = None,
        sequence_type: _TYPES_SEQ = "oligo",
        save_description: bool = False,
    ):
        """Write oligo sequences from the database to a FASTA file.

        This function writes the sequences from the oligo database to a FASTA file. Each sequence is represented as a
        SeqRecord in the FASTA file.

        :param filename: The base filename for the output FASTA file (without extension), defaults to "db_oligo".
        :type filename: str, optional
        :param region_ids: A list of region IDs to include in the saved database. If None, all regions are included.
        :type region_ids: list[str], optional
        :param sequence_type: The type of sequence to write (e.g., "oligo" or "target").
        :type sequence_type: str
        :param save_description: Flag to include sequence descriptions in the FASTA file.
        :type save_description: bool, optional
        :return: Path to the generated FASTA file.
        :rtype: str
        """
        # Check if sequence type is correct
        options = get_args(_TYPES_SEQ)
        assert (
            sequence_type in options
        ), f"Sequence type not supported! '{sequence_type}' is not in {options}."

        # check if region ids are list
        if region_ids:
            region_ids = check_if_list(region_ids)
        else:
            region_ids = self.database.keys()

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
        filename: str = "db_oligo",
        region_ids: list[str] = None,
    ):
        """Write the oligo database to TSV files.

        This function saves the oligo database to a TSV (tab-separated values) file containing the
        oligo database entries. The files are saved in the specified output directory with the provided
        filenames. The order of the columns in the TSV file is:

        +-----------+----------+----------+------------+-------+-----+--------+--------+------------------+
        | region_id | oligo_id | sequence | chromosome | start | end | strand | length | additional feat. |
        +-----------+----------+----------+------------+-------+-----+--------+--------+------------------+

        Additional feat. includes additional information from the header of the fasta file and additional information
        computed by the filtering classes.

        :param filename: The base filename for the output files (without extensions), defaults to "db_oligo".
        :type filename: str, optional
        :param region_ids: A list of region IDs to include in the saved database. If None, all regions are included.
        :type region_ids: list[str], optional
        :return: Paths to the saved TSV file.
        :rtype: tuple
        """
        # TODO: this should be split into write database (without brackets formatting) and save database (with brackets formatting that can be reloaded)
        if region_ids:
            region_ids = check_if_list(region_ids)
        else:
            region_ids = self.database.keys()

        file_database = os.path.join(self.dir_output, filename + ".tsv")

        first_entry = True
        for region_id in region_ids:
            database_region = self.database[region_id]
            file_tsv_content = []
            for oligo_id, oligo_attributes in database_region.items():
                entry = {"region_id": region_id, "oligo_id": oligo_id}
                entry.update(oligo_attributes)
                file_tsv_content.append(entry)
            file_tsv_content = pd.DataFrame(data=file_tsv_content)
            file_tsv_content.to_csv(file_database, sep="\t", index=False, mode="a", header=first_entry)
            first_entry = False

        return file_database

    def write_oligosets_to_yaml(
        self,
        attributes: list[str],
        top_n_sets: int,
        ascending: bool,
        filename: str = "oligos",
        region_ids: list[str] = None,
    ):
        """Write the top N oligosets to a YAML file with specified attributes.

        This function writes the specified attributes of the top N oligosets from the database to a YAML file. The
        oligosets are sorted based on their scores in ascending or descending order. If region IDs are specified,
        only those regions are included in the output.

        :param attributes: List of attributes to include for each oligo in the YAML file.
        :type attributes: list[str]
        :param top_n_sets: Number of top oligosets to include in the output.
        :type top_n_sets: int
        :param ascending: Whether to sort the oligosets in ascending order.
        :type ascending: bool
        :param filename: Name of the output YAML file, defaults to "oligos".
        :type filename: str, optional
        :param region_ids: List of region IDs to include in the output. If None, all regions are included.
        :type region_ids: list[str], optional
        """
        file_yaml = os.path.join(os.path.dirname(self.dir_output), filename)

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
                            yaml_dict_oligo_entry[attribute] = (
                                str(self.database[region_id][oligo_id][attribute])
                                .replace("'", "")
                                .replace("[[", "[")
                                .replace("]]", "]")
                            )

                    oligo_id_yaml = f"Oligo {oligo_idx + 1}"
                    yaml_dict[region_id][oligoset_id][oligo_id_yaml] = yaml_dict_oligo_entry

        with open(file_yaml, "w") as handle:
            yaml.dump(yaml_dict, handle, default_flow_style=False, sort_keys=False)

        return file_yaml

    def write_oligosets_to_table(self, foldername_out: str = "sets_of_oligos"):
        """Write oligo sets to individual TSV files.

        This function writes the oligo sets to individual TSV files, with each file representing the oligo sets
        for a specific region.

        :param foldername_out: The name of the folder to store the oligo set files, defaults to "sets_of_oligos".
        :type foldername_out: str, optional
        :return: Path to the folder containing the generated oligo set files.
        :rtype: str
        """
        dir_oligosets = os.path.join(os.path.dirname(self.dir_output), foldername_out)
        Path(dir_oligosets).mkdir(parents=True, exist_ok=True)

        for region_id in self.oligosets.keys():
            file_oligosets = os.path.join(dir_oligosets, f"{region_id}_oligosets.tsv")
            self.oligosets[region_id].to_csv(file_oligosets, sep="\t", index=False)

        return dir_oligosets

    def remove_regions_with_insufficient_oligos(self, pipeline_step: str):
        """Remove regions with insufficient oligos from the oligo database.

        This function identifies regions in the oligo database that have fewer oligos than the specified
        minimum threshold (`min_oligos_per_region`). It removes those regions from the database and updates
        the associated oligo sets. If the option to write removed regions to a file is enabled, it logs the
        removed regions along with the pipeline step.

        :param pipeline_step: Step in the pipeline that led to the removal of regions.
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

    def get_oligoid_list(self):
        """Retrieve a list of all oligo IDs in the database.

        This function iterates over the database regions and collects all oligo IDs into a single list.

        :return: A list containing all oligo IDs in the database.
        :rtype: list[str]
        """
        oligo_ids = [
            oligo_id for database_region in self.database.values() for oligo_id in database_region.keys()
        ]

        return oligo_ids

    def get_sequence_list(self, sequence_type: _TYPES_SEQ = "oligo"):
        """Retrieve a list of sequences of the specified type (e.g., 'oligo' or 'target') from the oligo database.

        :param sequence_type: Type of sequences to retrieve (default is 'oligo').
        :type sequence_type: _TYPES_SEQ
        :return: List of sequences.
        :rtype: List[str]
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
        self, sequence_type: _TYPES_SEQ = "oligo", sequence_to_upper: bool = False
    ):
        """Generate a mapping between oligonucleotide IDs and their corresponding sequences, with an option to convert sequences to uppercase.

        Validates the sequence type against predefined options. If `sequence_to_upper` is True, converts all sequences to uppercase before mapping,
        ensuring case-insensitive comparisons.

        :param sequence_type: The type of sequence to use for the mapping, defaulting to "oligo".
        :type sequence_type: _TYPES_SEQ
        :param sequence_to_upper: Flag indicating whether to convert sequences to uppercase.
        :type sequence_to_upper: bool
        :return: A dictionary mapping oligonucleotide IDs to corresponding sequences.
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
        self, sequence_type: _TYPES_SEQ = "oligo", sequence_to_upper: bool = False
    ):
        """Generate a mapping between sequences and their corresponding oligonucleotide IDs, with an option to convert sequences to uppercase.

        Validates the sequence type against predefined options. If `sequence_to_upper` is True, converts all sequences to uppercase before mapping,
        ensuring case-insensitive comparisons. This function supports scenarios where multiple oligos share the same sequence, grouping their IDs in a list.

        :param sequence_type: The type of sequence to use for the mapping, defaulting to "oligo".
        :type sequence_type: _TYPES_SEQ
        :param sequence_to_upper: Flag indicating whether to convert sequences to uppercase.
        :type sequence_to_upper: bool
        :return: A dictionary mapping sequences to lists of corresponding oligonucleotide IDs.
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

    def get_oligo_attribute(self, attribute: str, region_ids: Union[str, List[str]] = None):
        """Retrieves a specified attribute for all oligos in the database and returns it as a pandas DataFrame.
        This method assumes the presence of an attribute across all oligo records in the database. If the
        attribute is not found, a KeyError is raised.

        :param attribute: The name of the attribute to retrieve for each oligo.
        :type attribute: str
        :return: A pandas DataFrame with two columns: 'oligo_id' and the specified 'attribute', where each row
                corresponds to an oligo and its attribute value.
        :rtype: pd.DataFrame
        :raises KeyError: If the specified attribute has not been computed and added to the database.
        """
        if region_ids is None:
            region_ids = self.database.keys()
        else:
            region_ids = check_if_list(region_ids)

        oligo_ids = []
        attributes = []

        for region_id in region_ids:
            database_region = self.database[region_id]
            if not check_if_key_exists(database_region, attribute):
                warnings.warn(
                    f"The {attribute} attribute has not been computed for {region_id}! Setting to None!"
                )
            for oligo_id, oligo_attributes in database_region.items():
                if attribute in oligo_attributes:
                    oligo_ids.append(oligo_id)
                    attributes.append(oligo_attributes[attribute])
                else:
                    oligo_ids.append(oligo_id)
                    attributes.append(None)

        return pd.DataFrame({"oligo_id": oligo_ids, attribute: attributes})

    ############################################
    # Manipulation Functions
    ############################################

    def update_oligo_attributes(self, new_oligo_attribute: dict):
        """Update attributes of oligonucleotides in the database with new values.

        Iterates through each oligonucleotide in the database and updates its attributes based on the provided new attribute values.
        The update is done in place, modifying the existing attributes of the oligonucleotides.

        Input format of new attribute:
        new_attribute = {"region_x::1": {"attribute_y": 110}}

        :param new_oligo_attribute: A dictionary mapping oligonucleotide IDs to their new attributes.
        :type new_oligo_attribute: dict
        """
        for region_id, database_region in self.database.items():
            for oligo_id, oligo_attributes in database_region.items():
                oligo_attributes.update(new_oligo_attribute[oligo_id])

    def filter_oligo_attribute(
        self, name_attribute: str, thr_attribute: float, keep_if_smaller_threshold: bool
    ):
        """Filter oligos in the database based on a specific attribute and threshold.

        This function iterates through the oligo database and removes oligos that do not meet the specified threshold
        for the given attribute.

        :param name_attribute: The name of the attribute to filter oligos by.
        :type name_attribute: str
        :param thr_attribute: The threshold value for the attribute.
        :type thr_attribute: float
        :param keep_if_smaller_threshold: If True, keep oligos with attribute values smaller than the threshold.
                                        If False, keep oligos with attribute values larger than the threshold.
        :type keep_if_smaller_threshold: bool
        """
        oligos_to_delete = []
        for region_id, database_region in self.database.items():
            for oligo_id, oligo_attributes in database_region.items():
                if (keep_if_smaller_threshold) and (oligo_attributes[name_attribute] > thr_attribute):
                    oligos_to_delete.append((region_id, oligo_id))
                elif (not keep_if_smaller_threshold) and (oligo_attributes[name_attribute] < thr_attribute):
                    oligos_to_delete.append((region_id, oligo_id))

        for region_id, oligo_id in oligos_to_delete:
            del self.database[region_id][oligo_id]
