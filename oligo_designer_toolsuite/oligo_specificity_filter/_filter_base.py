############################################
# imports
############################################

import os
import re
from abc import ABC, abstractmethod
from pathlib import Path
from typing import List, Union, get_args

import pandas as pd
from Bio import Seq
from joblib import Parallel, delayed

from oligo_designer_toolsuite._constants import (
    _TYPES_SEQ,
    SEPARATOR_FASTA_HEADER_FIELDS,
    SEPARATOR_OLIGO_ID,
)
from oligo_designer_toolsuite.database import OligoDatabase, ReferenceDatabase

############################################
# Oligo Specificity Filter Classes
############################################


class SpecificityFilterBase(ABC):
    """An abstract base class for implementing specificity filters for oligonucleotide sequences.
    Specificity filters are used to ensure that oligos do not have unintended matches within a given
    reference database. This base class provides a common structure for such filters, including an output
    directory setup and an abstract method for applying the filter to an oligo database.

    :param filter_name: Subdirectory path for the output, i.e. <dir_output>/<filter_name>.
    :type filter_name: str
    :param dir_output: The directory where intermediate files will be saved.
    :type dir_output: str
    """

    def __init__(self, filter_name: str, dir_output: str):
        """Constructor for the SpecificityFilterBase class."""
        # folder where we write the intermediate files
        self.filter_name = filter_name
        self.dir_output = os.path.abspath(os.path.join(dir_output, self.filter_name))
        Path(self.dir_output).mkdir(parents=True, exist_ok=True)

    @abstractmethod
    def apply(
        self,
        sequence_type: _TYPES_SEQ,
        oligo_database: OligoDatabase,
        reference_database: ReferenceDatabase = None,
        n_jobs: int = 1,
    ):
        """Abstract method to apply the specificity filter to an oligo database.

        :param sequence_type: The type of sequences being filtered, must be one of the predefined sequence types.
        :type sequence_type: _TYPES_SEQ
        :param oligo_database: The oligo database to which the filter will be applied.
        :type oligo_database: OligoDatabase
        :param reference_database: The reference database to compare against for specificity.
            For non-alignment based specificity filter reference_database is not used, i.e. set to None.
        :type reference_database: ReferenceDatabase, optional
        :param n_jobs: The number of parallel jobs to run.
        :type n_jobs: int
        """

    def _filter_hits_from_database(self, database_region: dict, oligos_with_hits: list[str]):
        """Removes oligos identified with hits in the reference database from a given region of the oligo database.

        :param database_region: A region of the oligo database to filter.
        :type database_region: dict
        :param oligos_with_hits: A list of oligo IDs that have matches in the reference database and should be removed.
        :type oligos_with_hits: list[str]
        :return: The filtered region of the oligo database.
        :rtype: dict
        """
        oligo_ids = list(database_region.keys())
        for oligo_id in oligo_ids:
            if oligo_id in oligos_with_hits:
                del database_region[oligo_id]
        return database_region


class AlignmentSpecificityFilter(SpecificityFilterBase):
    """A class that implements specificity filtering for oligonucleotides through alignments against a reference database.
    This filter creates an index of the reference database and then aligns oligonucleotide sequences to identify and exclude
    sequences with significant hits in the reference, ensuring specificity of the oligonucleotide sequences.

    :param filter_name: Subdirectory path for the output, i.e. <dir_output>/<filter_name>.
    :type filter_name: str
    :param dir_output: Directory for saving intermediate files.
    :type dir_output: str
    """

    def __init__(
        self,
        filter_name: str,
        dir_output: str,
    ):
        """Constructor for the AlignmentSpecificityFilter class."""
        # folder where we write the intermediate files
        self.filter_name = filter_name
        self.dir_output = os.path.abspath(os.path.join(dir_output, self.filter_name))
        Path(self.dir_output).mkdir(parents=True, exist_ok=True)

    def apply(
        self,
        sequence_type: _TYPES_SEQ,
        oligo_database: OligoDatabase,
        reference_database: ReferenceDatabase,
        n_jobs: int = 1,
    ):
        """Applies the alignment-based specificity filter to an oligonucleotide database.

        :param sequence_type: The type of sequences being filtered, must be one of the predefined sequence types.
        :type sequence_type: _TYPES_SEQ
        :param oligo_database: The oligo database to which the filter will be applied.
        :type oligo_database: OligoDatabase
        :param reference_database: The reference database to compare against for specificity.
        :type reference_database: ReferenceDatabase
        :param n_jobs: The number of parallel jobs to run.
        :type n_jobs: int
        :return: The filtered oligo database with sequences having significant hits removed.
        :rtype: OligoDatabase

        """
        region_ids = list(oligo_database.database.keys())
        table_hits = self._get_table_hits(
            sequence_type=sequence_type,
            oligo_database=oligo_database,
            reference_database=reference_database,
            region_ids=region_ids,
            n_jobs=n_jobs,
        )

        for region_id, table_hits_region in zip(region_ids, table_hits):
            oligos_with_hits_region = table_hits_region["query"].unique()
            database_region_filtered = self._filter_hits_from_database(
                database_region=oligo_database.database[region_id],
                oligos_with_hits=oligos_with_hits_region,
            )
            oligo_database.database[region_id] = database_region_filtered

        return oligo_database

    def get_oligo_pair_hits(
        self,
        sequence_type: _TYPES_SEQ,
        oligo_database: OligoDatabase,
        reference_database: ReferenceDatabase,
        n_jobs: int,
    ):
        """Identifies oligonucleotide pairs with significant hits in the reference database.

        :param sequence_type: The type of sequences being filtered, must be one of the predefined sequence types.
        :type sequence_type: _TYPES_SEQ
        :param oligo_database: The oligo database to which the filter will be applied.
        :type oligo_database: OligoDatabase
        :param reference_database: The reference database to compare against for specificity.
        :type reference_database: ReferenceDatabase
        :param n_jobs: The number of parallel jobs to run.
        :type n_jobs: int
        :return: List of oligo pairs with hits in the reference database.
        :rtype: list[tuple]
        """
        region_ids = list(oligo_database.database.keys())
        table_hits = self._get_table_hits(
            sequence_type=sequence_type,
            oligo_database=oligo_database,
            reference_database=reference_database,
            region_ids=region_ids,
            n_jobs=n_jobs,
        )

        table_hits = pd.concat(table_hits, ignore_index=True)
        oligo_pair_hits = list(zip(table_hits["query"].values, table_hits["reference"].values))

        return oligo_pair_hits

    def _get_table_hits(
        self,
        sequence_type: _TYPES_SEQ,
        oligo_database: OligoDatabase,
        reference_database: ReferenceDatabase,
        region_ids: List[str],
        n_jobs: int,
    ) -> List[pd.DataFrame]:
        """Applies the alignment-based specificity filter to an oligonucleotide database and return a DataFrame containing all hits.

        :param sequence_type: The type of sequences being filtered, must be one of the predefined sequence types.
        :type sequence_type: _TYPES_SEQ
        :param oligo_database: The oligo database to which the filter will be applied.
        :type oligo_database: OligoDatabase
        :param reference_database: The reference database to compare against for specificity.
        :type reference_database: ReferenceDatabase
        :param region_ids: Regions for which to generate the table hits.
        :type region_ids: List[str]
        :param n_jobs: The number of parallel jobs to run.
        :type n_jobs: int
        :return: List of tables containing all hits. Each table contain the hits of one single region, and the list is ordered by the list region_ids.
        :rtype: List[pd.DataFrame]
        """
        options = get_args(_TYPES_SEQ)
        assert (
            sequence_type in options
        ), f"Sequence type not supported! '{sequence_type}' is not in {options}."

        # Create index file for search
        file_reference = reference_database.write_database_to_fasta(
            filename=f"db_reference_{self.filter_name}"
        )
        file_index = self._create_index(file_reference=file_reference, n_jobs=n_jobs)

        # Run search for each region in parallel
        table_hits = Parallel(n_jobs=n_jobs)(
            delayed(self._run_filter)(
                sequence_type=sequence_type,
                region_id=region_id,
                oligo_database=oligo_database,
                file_index=file_index,
                consider_hits_from_input_region=True,
            )
            for region_id in region_ids
        )
        os.remove(file_reference)
        self._remove_index(file_index)
        return table_hits

    @abstractmethod
    def _create_index(self, file_reference: str, n_jobs: int):
        """Abstract method to create an index of the reference database for alignment.

        :param file_reference: Path to the reference database fasta file.
        :type file_reference: str
        :param n_jobs: The number of parallel jobs to run.
        :type n_jobs: int
        """

    def _run_filter(
        self,
        sequence_type: _TYPES_SEQ,
        region_id: str,
        oligo_database: OligoDatabase,
        file_index: str,
        consider_hits_from_input_region: bool,
    ):
        """Executes the filtering process for a specific region of the oligonucleotide database based on search results.

        :param sequence_type: The type of sequences being filtered, must be one of the predefined sequence types.
        :type sequence_type: _TYPES_SEQ
        :param region_id: The identifier for the region within the database to filter.
        :type region_id: str
        :param oligo_database: The oligonucleotide database to apply the filter on.
        :type oligo_database: OligoDatabase
        :param file_index: Path to the index file used for the reference database.
        :type file_index: str
        :param consider_hits_from_input_region: Flag to indicate whether hits from the input region should be considered.
        :type consider_hits_from_input_region: bool
        :return: A tuple containing a table of hits and a list of oligos with those hits.
        :rtype: (pd.DataFrame, list)
        """
        search_results = self._run_search(
            sequence_type=sequence_type,
            oligo_database=oligo_database,
            file_index=file_index,
            region_ids=region_id,
        )
        table_hits = self._find_hits(
            oligo_database=oligo_database,
            search_results=search_results,
            consider_hits_from_input_region=consider_hits_from_input_region,
        )

        return table_hits

    @abstractmethod
    def _run_search(
        self,
        sequence_type: _TYPES_SEQ,
        oligo_database: OligoDatabase,
        file_index: str,
        region_ids: Union[str, List[str]] = None,
    ):
        """Abstract method to execute a search against a reference database using a specific indexing strategy.

        :param sequence_type: The type of sequences being filtered, must be one of the predefined sequence types.
        :type sequence_type: _TYPES_SEQ
        :param oligo_database: The oligonucleotide database to search against.
        :type oligo_database: OligoDatabase
        :param region_ids: Identifiers for the regions within the database to be searched.
        :type region_ids: str
        :param file_index: Path to the index file of the reference database.
        :type file_index: str
        """

    def _read_search_output(self, file_search_results: str, names_search_output: list, usecols: list = None):
        """Reads the output from a search operation into a pandas DataFrame.

        :param file_search_results: Path to the file containing search results.
        :type file_search_results: str
        :param names_search_output: Column names for the search result data frame.
        :type names_search_output: list[str]
        :param usecols: Specific columns to use from the search results. If None, all columns are used.
        :type usecols: list, optional
        :return: A pandas DataFrame containing the search results.
        :rtype: pd.DataFrame
        """
        search_results = pd.read_csv(
            filepath_or_buffer=file_search_results,
            header=None,
            sep="\t",
            low_memory=False,
            engine="c",
            usecols=usecols,
            names=names_search_output,
        )

        search_results["query_region_id"] = search_results["query"].str.split(SEPARATOR_OLIGO_ID).str[0]
        search_results["reference_region_id"] = (
            search_results["reference"].str.split(SEPARATOR_FASTA_HEADER_FIELDS).str[0]
        )

        return search_results

    @abstractmethod
    def _find_hits(
        self,
        oligo_database: OligoDatabase,
        search_results: pd.DataFrame,
        consider_hits_from_input_region: bool,
    ):
        """Abstract method to identify hits from search results within the oligonucleotide database.

        :param oligo_database: The oligonucleotide database.
        :type oligo_database: OligoDatabase
        :param search_results: DataFrame containing search results.
        :type search_results: pd.DataFrame
        :param consider_hits_from_input_region: Flag to indicate whether hits from the input region should be considered.
        :type consider_hits_from_input_region: bool
        """

    def _remove_index(self, file_index: str):
        """Remove all the temporary index files generated by the alignment method.

        :param file_index: Path to the index files (the extension is not specified).
        :type file_index: str
        """
        file_index_basename = os.path.basename(file_index)
        regex = re.compile(file_index_basename + "\..*")
        for root, _, files in os.walk(self.dir_output):
            for file in files:
                if regex.match(file):
                    os.remove(os.path.join(root, file))

    def _get_queries(
        self,
        sequence_type: _TYPES_SEQ,
        table_hits: pd.DataFrame,
        oligo_database: OligoDatabase,
        region_id: str,
    ) -> List[Seq.Seq]:
        """Abstract method to retrieve the queries sequences from the search results.

        :param sequence_type: The type of sequences being filtered, must be one of the predefined sequence types.
        :type sequence_type: _TYPES_SEQ
        :param table_hits: Dataframe containing the search results.
        :type table_hits: pd.DataFrame
        :param oligo_database: The oligo database to which the filter will be applied.
        :type oligo_database: OligoDatabase
        :param region_id: The identifier for the region within the database to filter.
        :type region_id: str
        :return: Queries sequences.
        :rtype: list
        """
        queries = [
            oligo_database.database[region_id][query_id][sequence_type] for query_id in table_hits["query"]
        ]
        return queries

    @abstractmethod
    def _get_references(self, table_hits: pd.DataFrame, file_reference: str, region_id: str):
        """Abstract method to retrieve the reference sequences from the search results.

        :param table_hits: Dataframe containing the search results.
        :type table_hits: pd.DataFrame
        :param file_reference: Path to the fasta file used as reference for the search.
        :type file_reference: str
        :param region_id: The identifier for the region within the database to filter.
        :type region_id: str
        """

    @abstractmethod
    def _add_alignment_gaps(
        self, table_hits: pd.DataFrame, queries: List[Seq.Seq], references: List[Seq.Seq]
    ):
        """Abstract method to add gaps to the references and queries sequences.

        :param table_hits: Dataframe containing the search results.
        :type table_hits: pd.DataFrame
        :param queries: List of the queries sequences.
        :type queries: List[Seq.Seq]
        :param references: List of the references sequences.
        :type references: List[Seq.Seq]
        """
