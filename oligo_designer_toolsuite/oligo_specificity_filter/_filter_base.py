############################################
# imports
############################################

import os
import re
from abc import ABC, abstractmethod
from pathlib import Path
from typing import List, Union, Tuple

import pandas as pd
from Bio import Seq
from joblib import Parallel, delayed
from joblib_progress import joblib_progress

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
    """
    A base class for implementing specificity filters that operate on OligoDatabase.

    The `SpecificityFilterBase` class provides the structure for creating filters
    that assess the specificity of oligonucleotides. These filters can be customized and extended
    to apply various criteria to an OligoDatabase, helping to refine and select optimal oligos for
    specific applications.

    :param filter_name: Name of the filter for identification purposes.
    :type filter_name: str
    :param dir_output: Directory to store output files and temporary data.
    :type dir_output: str
    """

    def __init__(self, filter_name: str, dir_output: str) -> None:
        """Constructor for the SpecificityFilterBase class."""
        # folder where we write the intermediate files
        self.filter_name = filter_name
        self.dir_output = os.path.abspath(os.path.join(dir_output, self.filter_name))
        Path(self.dir_output).mkdir(parents=True, exist_ok=True)

    @abstractmethod
    def apply(
        self,
        oligo_database: OligoDatabase,
        reference_database: ReferenceDatabase,
        sequence_type: _TYPES_SEQ,
        n_jobs: int = 1,
    ) -> OligoDatabase:
        """
        Apply the specificity filter to the given OligoDatabase.

        This abstract method defines the interface for applying the specificity filter, using the oligonucleotide sequences
        in the OligoDatabase as well as a ReferenceDatabase and a specified sequence type. The implementation should return
        a filtered version of the OligoDatabase.

        :param oligo_database: The OligoDatabase containing the oligonucleotides and their associated attributes.
        :type oligo_database: OligoDatabase
        :param reference_database: The ReferenceDatabase used for alignment.
        :type reference_database: ReferenceDatabase
        :param sequence_type: The type of sequence to be used for the filter calculations.
        :type sequence_type: _TYPES_SEQ["oligo", "target"]
        :param n_jobs: The number of parallel jobs to use for processing.
        :type n_jobs: int
        :return: The filtered OligoDatabase.
        :rtype: OligoDatabase
        """

    def _filter_hits_from_database(
        self, oligo_database: OligoDatabase, region_id: str, oligos_with_hits: list[str]
    ) -> None:
        """
        Remove oligonucleotides with hits from the database for a specific region.

        This method iterates over oligonucleotides in a specified region and removes any that have been identified
        as having hits (i.e., matches) based on the filtering criteria.

        :param oligo_database: The OligoDatabase containing the oligonucleotides and their associated attributes.
        :type oligo_database: OligoDatabase
        :param region_id: Region ID to process.
        :type region_id: str
        :param oligos_with_hits: A list of oligonucleotide IDs that have been identified as hits and should be removed.
        :type oligos_with_hits: list[str]
        """
        oligo_ids = list(oligo_database.database[region_id].keys())
        for oligo_id in oligo_ids:
            if oligo_id in oligos_with_hits:
                del oligo_database.database[region_id][oligo_id]


class AlignmentSpecificityFilter(SpecificityFilterBase):
    """
    A base class for implementing filters that utilize sequence alignment methods to evaluate oligonucleotide specificity.

    The `AlignmentSpecificityFilter` class provides a framework for developing filters that assess the potential off-target effects
    of oligonucleotides by aligning them against reference sequences.

    :param filter_name: Name of the filter for identification purposes.
    :type filter_name: str
    :param dir_output: Directory to store output files and temporary data.
    :type dir_output: str
    """

    def __init__(
        self,
        filter_name: str,
        dir_output: str,
    ) -> None:
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
    ) -> OligoDatabase:
        """
        Applies the alignment-based specificity filter to an OligoDatabase, filtering out sequences with off-target hits.

        :param sequence_type: The type of sequence to be used for the filter calculations.
        :type sequence_type: _TYPES_SEQ["oligo", "target"]
        :param oligo_database: The OligoDatabase containing the oligonucleotides and their associated attributes.
        :type oligo_database: OligoDatabase
        :param reference_database: The ReferenceDatabase used for alignment.
        :type reference_database: ReferenceDatabase
        :param n_jobs: The number of parallel jobs to use for processing.
        :type n_jobs: int
        :return: The filtered OligoDatabase.
        :rtype: OligoDatabase
        """
        # when applying filters we don't want to consider hits within the same region
        consider_hits_from_input_region = False

        # Create index file for search
        file_reference = reference_database.write_database_to_fasta(
            filename=f"db_reference_{self.filter_name}"
        )
        file_index = self._create_index(file_reference=file_reference, n_jobs=n_jobs)

        # run search in parallel for each region
        region_ids = list(oligo_database.database.keys())
        name = " ".join(string.capitalize() for string in self.filter_name.split("_"))
        with joblib_progress(description=f"Specificity Filter: {name}", total=len(region_ids)):
            Parallel(n_jobs=n_jobs, prefer="threads", require="sharedmem")(
                delayed(self._apply_region)(
                    sequence_type=sequence_type,
                    oligo_database=oligo_database,
                    file_index=file_index,
                    region_id=region_id,
                    consider_hits_from_input_region=consider_hits_from_input_region,
                )
                for region_id in region_ids
            )

        os.remove(file_reference)
        self._remove_index(file_index)

        return oligo_database

    def _apply_region(
        self,
        sequence_type: _TYPES_SEQ,
        oligo_database: OligoDatabase,
        file_index: str,
        region_id: str,
        consider_hits_from_input_region: bool,
    ) -> None:
        """
        Applies the filter to a specific region in the OligoDatabase.

        :param sequence_type: The type of sequence to be used for the filter calculations.
        :type sequence_type: _TYPES_SEQ["oligo", "target"]
        :param oligo_database: The OligoDatabase containing the oligonucleotides and their associated attributes.
        :type oligo_database: OligoDatabase
        :param file_index: The index file used for alignment filtering.
        :type file_index: str
        :param region_id: Region ID to process.
        :type region_id: str
        :param consider_hits_from_input_region: Whether to consider hits from the same region, defaults to False.
        :type consider_hits_from_input_region: bool
        """
        table_hits_region = self._run_filter(
            sequence_type=sequence_type,
            region_id=region_id,
            oligo_database=oligo_database,
            file_index=file_index,
            consider_hits_from_input_region=consider_hits_from_input_region,
        )

        oligos_with_hits_region = table_hits_region["query"].unique()
        self._filter_hits_from_database(
            oligo_database=oligo_database,
            region_id=region_id,
            oligos_with_hits=oligos_with_hits_region,
        )

    def get_oligo_pair_hits(
        self,
        sequence_type: _TYPES_SEQ,
        oligo_database: OligoDatabase,
        reference_database: ReferenceDatabase,
        n_jobs: int,
    ) -> list:
        """
        Retrieves pairs of oligonucleotides that have significant alignment hits within the same region.

        :param sequence_type: The type of sequence to be used for the filter calculations.
        :type sequence_type: _TYPES_SEQ["oligo", "target"]
        :param oligo_database: The OligoDatabase containing the oligonucleotides and their associated attributes.
        :type oligo_database: OligoDatabase
        :param reference_database: The ReferenceDatabase used for alignment.
        :type reference_database: ReferenceDatabase
        :param n_jobs: The number of parallel jobs to use for processing.
        :type n_jobs: int
        :return: A list of tuples representing oligo pairs that have significant hits.
        :rtype: list
        """
        # when getting oligo pair hits we want to consider hits within the same region
        consider_hits_from_input_region = True

        # Create index file for search
        file_reference = reference_database.write_database_to_fasta(
            filename=f"db_reference_{self.filter_name}"
        )
        file_index = self._create_index(file_reference=file_reference, n_jobs=n_jobs)

        region_ids = list(oligo_database.database.keys())
        name = " ".join(string.capitalize() for string in self.filter_name.split("_"))
        with joblib_progress(description=f"Specificity Filter: {name}", total=len(region_ids)):
            table_hits = Parallel(n_jobs=n_jobs, prefer="threads", require="sharedmem")(
                delayed(self._run_filter)(
                    sequence_type=sequence_type,
                    region_id=region_id,
                    oligo_database=oligo_database,
                    file_index=file_index,
                    consider_hits_from_input_region=consider_hits_from_input_region,
                )
                for region_id in region_ids
            )

        table_hits = pd.concat(table_hits, ignore_index=True)
        oligo_pair_hits = list(zip(table_hits["query"].values, table_hits["reference"].values))

        os.remove(file_reference)
        self._remove_index(file_index)

        return oligo_pair_hits

    @abstractmethod
    def _create_index(self, file_reference: str, n_jobs: int) -> str:
        """
        Abstract method to create an index for the reference file used in the alignment process.

        :param file_reference: Path to the reference file that needs to be indexed.
        :type file_reference: str
        :param n_jobs: Number of parallel jobs to use during the indexing process.
        :type n_jobs: int
        :return: The name of the created index file.
        :rtype: str
        """

    def _run_filter(
        self,
        sequence_type: _TYPES_SEQ,
        region_id: str,
        oligo_database: OligoDatabase,
        file_index: str,
        consider_hits_from_input_region: bool,
    ) -> pd.DataFrame:
        """
        Executes the filtering process for a specific region by running the search and identifying significant hits.

        :param sequence_type: The type of sequence to be used for the filter calculations.
        :type sequence_type: _TYPES_SEQ["oligo", "target"]
        :param region_id: Region ID to process.
        :type region_id: str
        :param oligo_database: The OligoDatabase containing the oligonucleotides and their associated attributes.
        :type oligo_database: OligoDatabase
        :param file_index: The index file used for searching.
        :type file_index: str
        :param consider_hits_from_input_region: Flag to consider hits within the same input region.
        :type consider_hits_from_input_region: bool
        :return: A DataFrame containing the significant hits found in the region.
        :rtype: pd.DataFrame
        """
        search_results = self._run_search(
            oligo_database=oligo_database,
            file_index=file_index,
            sequence_type=sequence_type,
            region_ids=region_id,
        )
        table_hits = self._find_hits(
            oligo_database=oligo_database,
            search_results=search_results,
            consider_hits_from_input_region=consider_hits_from_input_region,
            region_ids=region_id,
        )

        return table_hits

    @abstractmethod
    def _run_search(
        self,
        oligo_database: OligoDatabase,
        file_index: str,
        sequence_type: _TYPES_SEQ,
        region_ids: Union[str, List[str]] = None,
    ) -> pd.DataFrame:
        """
        Abstract method to run a search against a ReferenceDatabase using a specified index file.

        :param oligo_database: The OligoDatabase containing the oligonucleotides and their associated attributes.
        :type oligo_database: OligoDatabase
        :param file_index: The path to the index file used for searching.
        :type file_index: str
        :param sequence_type: The type of sequence to be used for the filter calculations.
        :type sequence_type: _TYPES_SEQ["oligo", "target"]
        :param region_ids: List of region IDs to process. If None, all regions in the OligoDatabase are processed, defaults to None.
        :type region_ids: Union[str, List[str]], optional
        :return: A DataFrame containing the search results.
        :rtype: pd.DataFrame
        """

    def _read_search_output(
        self, file_search_results: str, names_search_output: list, usecols: list = None
    ) -> pd.DataFrame:
        """
        Reads and processes the output of a search, converting it into a structured DataFrame.

        :param file_search_results: Path to the file containing the raw search results.
        :type file_search_results: str
        :param names_search_output: List of column names for the search output.
        :type names_search_output: list
        :param usecols: Specific columns to read from the search output file, defaults to None.
        :type usecols: list, optional
        :return: A DataFrame containing the processed search results.
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
        region_ids: Union[str, List[str]],
    ) -> pd.DataFrame:
        """
        Abstract method to identify significant hits from search results, potentially excluding hits within the same region.

        :param oligo_database: The OligoDatabase containing the oligonucleotides and their associated attributes.
        :type oligo_database: OligoDatabase
        :param search_results: DataFrame containing the raw search results.
        :type search_results: pd.DataFrame
        :param consider_hits_from_input_region: Flag to consider or ignore hits from the same region.
        :type consider_hits_from_input_region: bool
        :param region_ids: List of region IDs to process. If None, all regions in the OligoDatabase are processed, defaults to None.
        :type region_ids: Union[str, List[str]], optional
        :return: A DataFrame with the identified significant hits.
        :rtype: pd.DataFrame
        """

    def _remove_index(self, file_index: str) -> None:
        """
        Removes the index files generated during the alignment process.

        :param file_index: The base name of the index file to be removed.
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
        oligo_database: OligoDatabase,
        table_hits: pd.DataFrame,
        sequence_type: _TYPES_SEQ,
        region_id: str,
    ) -> List[str]:
        """
        Retrieves the query sequences from the OligoDatabase based on the hit information.

        :param oligo_database: The OligoDatabase containing the oligonucleotides and their associated attributes.
        :type oligo_database: OligoDatabase
        :param table_hits: DataFrame containing the hits with query IDs.
        :type table_hits: pd.DataFrame
        :param sequence_type: The type of sequence to be used for the filter calculations.
        :type sequence_type: _TYPES_SEQ["oligo", "target"]
        :param region_id: Region ID to process.
        :type region_id: str
        :return: A list of query sequences corresponding to the hits.
        :rtype: List[str]
        """
        queries = [
            oligo_database.database[region_id][query_id][sequence_type] for query_id in table_hits["query"]
        ]
        return queries

    @abstractmethod
    def _get_references(self, table_hits: pd.DataFrame, file_reference: str, region_id: str) -> list:
        """
        Abstract method to retrieve the reference sequences from the ReferenceDatabase based on hit information.

        :param table_hits: DataFrame containing the hits with reference IDs and positions.
        :type table_hits: pd.DataFrame
        :param file_reference: The reference file containing the sequences.
        :type file_reference: str
        :param region_id: Region ID to process.
        :type region_id: str
        :return: A list of reference sequences corresponding to the hits.
        :rtype: list
        """

    @abstractmethod
    def _add_alignment_gaps(
        self,
        table_hits: pd.DataFrame,
        queries: list,
        references: list,
    ) -> Tuple[list, list]:
        """
        Abstract method to align query and reference sequences by adding gaps.

        :param table_hits: DataFrame containing the hits with alignment information.
        :type table_hits: pd.DataFrame
        :param queries: List of query sequences.
        :type queries: list
        :param references: List of reference sequences.
        :type references: list
        :return: A tuple containing lists of gapped query and reference sequences.
        :rtype: Tuple[list, list]
        """
