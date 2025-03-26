############################################
# imports
############################################

import os
import re
from abc import ABC, abstractmethod
from pathlib import Path
from typing import List, Tuple, get_args

import pandas as pd
from joblib import Parallel, delayed
from joblib_progress import joblib_progress

from oligo_designer_toolsuite._constants import (
    _TYPES_SEQ,
    SEPARATOR_FASTA_HEADER_FIELDS,
    SEPARATOR_OLIGO_ID,
)
from oligo_designer_toolsuite.database import OligoDatabase, ReferenceDatabase
from oligo_designer_toolsuite.utils import check_if_list

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

        self.sequence_type = None

    @abstractmethod
    def apply(
        self,
        oligo_database: OligoDatabase,
        sequence_type: _TYPES_SEQ,
        n_jobs: int = 1,
    ) -> OligoDatabase:
        """
        Apply the specificity filter to the given OligoDatabase.

        This abstract method defines the interface for applying the specificity filter to the oligonucleotide sequences
        in the OligoDatabase specified by the sequence type. The implementation should return a filtered version of the
        OligoDatabase.

        :param oligo_database: The OligoDatabase containing the oligonucleotides and their associated attributes.
        :type oligo_database: OligoDatabase
        :param sequence_type: The type of sequence to be used for filter calculations.
        :type sequence_type: _TYPES_SEQ["oligo", "target"]
        :param n_jobs: The number of parallel jobs to use for processing.
        :type n_jobs: int
        :return: The filtered OligoDatabase.
        :rtype: OligoDatabase
        """

    def _filter_hits_from_database(
        self,
        oligo_database: OligoDatabase,
        region_ids: str,
        oligos_with_hits: dict,
    ) -> None:
        """
        Remove oligonucleotides with hits from the database for a specific region.

        This method iterates over oligonucleotides in a specified region and removes any that have been identified
        as having hits (i.e., matches) based on the filtering criteria.

        :param oligo_database: The OligoDatabase containing the oligonucleotides and their associated attributes.
        :type oligo_database: OligoDatabase
        :param region_ids: List of region IDs to retrieve. If None, all regions in the database are retrieved, defaults to None.
        :type region_ids: Union[str, List[str]], optional
        :param oligos_with_hits: A dictionary of regio_ids associated with oligo_ids that have been identified as hits and should be removed.
        :type oligos_with_hits: dict
        """
        region_ids = check_if_list(region_ids) if region_ids else oligo_database.database.keys()

        for region_id in region_ids:
            oligo_ids = list(oligo_database.database[region_id].keys())
            for oligo_id in oligo_ids:
                if oligo_id in oligos_with_hits[region_id]:
                    del oligo_database.database[region_id][oligo_id]

    def _flag_hits_in_database(
        self,
        oligo_database: OligoDatabase,
        region_ids: str,
        oligos_with_hits: dict,
        oligos_with_hits_attributes: dict,
    ) -> None:
        """
        Flags oligos in the database based on whether they have hits, and assigns hit attributes.

        This method iterates over oligonucleotides in a specified region and flags each oligo based on
        whether it appears in the list of hits. If an oligo has a hit, the corresponding reference is assigned.
        Otherwise, the attribute is set to None.

        :param oligo_database: The OligoDatabase containing the oligonucleotides and their associated attributes.
        :type oligo_database: OligoDatabase
        :param region_ids: List of region IDs to retrieve. If None, all regions in the database are retrieved, defaults to None.
        :type region_ids: Union[str, List[str]], optional
        :param oligos_with_hits: A dictionary of regio_ids associated with oligo_ids that have been identified as hits and should be removed.
        :type oligos_with_hits: dict
        :param oligos_with_hits_attributes: Dictionary mapping oligo IDs to attributes related to their hits.
        :type oligos_with_hits_attributes: dict
        """
        region_ids = check_if_list(region_ids) if region_ids else oligo_database.database.keys()

        for region_id in region_ids:
            oligo_ids = list(oligo_database.database[region_id].keys())
            for oligo_id in oligo_ids:
                if oligo_id in oligos_with_hits[region_id]:
                    oligo_database.database[region_id][oligo_id][self.filter_name] = (
                        oligos_with_hits_attributes[oligo_id]
                    )
                else:
                    oligo_database.database[region_id][oligo_id][self.filter_name] = None


class SpecificityFilterReference(SpecificityFilterBase):
    """
    A base class for implementing specificity filters using a reference database.

    The `SpecificityFilterReference` class provides a framework for developing filters that
    assess the potential off-target effects of oligonucleotides wrt reference sequences.

    :param remove_hits: If True, oligos overlapping variants are removed. If False, they are flagged.
    :type remove_hits: bool
    :param filter_name: Name of the filter for identification purposes.
    :type filter_name: str
    :param dir_output: Directory to store output files and temporary data.
    :type dir_output: str
    """

    def __init__(
        self,
        remove_hits: bool,
        filter_name: str,
        dir_output: str,
    ) -> None:
        """Constructor for the SpecificityFilterReference class."""
        # folder where we write the intermediate files
        self.filter_name = filter_name
        self.dir_output = os.path.abspath(os.path.join(dir_output, self.filter_name))
        Path(self.dir_output).mkdir(parents=True, exist_ok=True)

        self.remove_hits = remove_hits
        self.reference_database = None
        self.sequence_type = None

    def set_reference_database(self, reference_database: ReferenceDatabase) -> None:
        """
        Set the ReferenceDatabase for reference based specificity filters.

        :param reference_database: The ReferenceDatabase used for alignment.
        :type reference_database: ReferenceDatabase
        """
        self.reference_database = reference_database

    @abstractmethod
    def create_reference(self, n_jobs: int) -> str:
        """
        Abstract method to write a reference database to file and create an index in case of alignment based methods.

        :param n_jobs: Number of parallel jobs to use during the indexing process.
        :type n_jobs: int
        :return: The name of the created reference file.
        :rtype: str
        """

    def remove_reference(self, file_reference: str) -> None:
        """
        Removes the reference files created for the filter.

        :param file_reference: The base name of the reference files to be removed.
        :type file_reference: str
        """
        file_reference_basename = os.path.basename(file_reference)
        regex = re.compile(file_reference_basename + "\..*")
        for root, _, files in os.walk(self.dir_output):
            for file in files:
                if regex.match(file):
                    os.remove(os.path.join(root, file))
        os.remove(file_reference)

    def _read_search_output(
        self,
        file_search_results: str,
        names_search_output: list,
        usecols: list = None,
        parse_query: bool = True,
        parse_reference: bool = True,
    ) -> pd.DataFrame:
        """
        Reads and processes the output of a search, converting it into a structured DataFrame.

        :param file_search_results: Path to the file containing the raw search results.
        :type file_search_results: str
        :param names_search_output: List of column names for the search output.
        :type names_search_output: list
        :param usecols: List of column indices to read from the file. If None, all columns are read.
        :type usecols: list, optional
        :param parse_query: Whether to parse and extract the region ID from the query sequence identifier.
        :type parse_query: bool, optional
        :param parse_reference: Whether to parse and extract the region ID from the reference sequence identifier.
        :type parse_reference: bool, optional
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

        if parse_query:
            search_results["query_region_id"] = search_results["query"].str.split(SEPARATOR_OLIGO_ID).str[0]
        if parse_reference:
            search_results["reference_region_id"] = (
                search_results["reference"].str.split(SEPARATOR_FASTA_HEADER_FIELDS).str[0]
            )

        return search_results


class SpecificityFilterAlignment(SpecificityFilterReference):
    """
    A base class for implementing filters that utilize sequence alignment methods to evaluate oligonucleotide specificity.

    The `AlignmentSpecificityFilter` class provides a framework for developing filters that assess the potential
    off-target effects of oligonucleotides by aligning them against reference sequences.

    :param remove_hits: If True, oligos overlapping variants are removed. If False, they are flagged.
    :type remove_hits: bool
    :param filter_name: Name of the filter for identification purposes.
    :type filter_name: str
    :param dir_output: Directory to store output files and temporary data.
    :type dir_output: str
    """

    def __init__(
        self,
        remove_hits: bool,
        filter_name: str,
        dir_output: str,
    ) -> None:
        """Constructor for the SpecificityFilterAlignment class."""

        # folder where we write the intermediate files
        self.filter_name = filter_name
        self.dir_output = os.path.abspath(os.path.join(dir_output, self.filter_name))
        Path(self.dir_output).mkdir(parents=True, exist_ok=True)

        self.remove_hits = remove_hits
        self.reference_database = None
        self.sequence_type = None

    def apply(
        self,
        oligo_database: OligoDatabase,
        sequence_type: _TYPES_SEQ,
        n_jobs: int = 1,
    ) -> OligoDatabase:
        """
        Applies the alignment-based specificity filter to an OligoDatabase, filtering out sequences with off-target hits.

        This function creates a reference index, runs the alignment-based specificity filter in parallel
        across all regions in the oligo database, and removes or flags oligos with off-target hits.

        :param oligo_database: The OligoDatabase containing the oligonucleotides and their associated attributes.
        :type oligo_database: OligoDatabase
        :param sequence_type: The type of sequence to be used for filter calculations.
        :type sequence_type: _TYPES_SEQ["oligo", "target"]
        :param n_jobs: The number of parallel jobs to use for processing.
        :type n_jobs: int
        :return: The filtered OligoDatabase.
        :rtype: OligoDatabase
        """
        self.sequence_type = sequence_type

        # when applying filters we don't want to consider hits within the same region
        consider_hits_from_input_region = False

        file_reference = self.create_reference(n_jobs=n_jobs)

        # run search in parallel for each region
        region_ids = list(oligo_database.database.keys())
        name = " ".join(string.capitalize() for string in self.filter_name.split("_"))
        with joblib_progress(description=f"Specificity Filter: {name}", total=len(region_ids)):
            Parallel(n_jobs=n_jobs, prefer="threads", require="sharedmem")(
                delayed(self._run_filter)(
                    region_id=region_id,
                    oligo_database=oligo_database,
                    file_reference=file_reference,
                    consider_hits_from_input_region=consider_hits_from_input_region,
                    mode=int(self.remove_hits),
                )
                for region_id in region_ids
            )

        self.remove_reference(file_reference)

        return oligo_database

    def get_oligo_pair_hits(
        self,
        oligo_database: OligoDatabase,
        sequence_type: _TYPES_SEQ,
        n_jobs: int,
    ) -> list:
        """
        Retrieves pairs of oligonucleotides that have significant alignment hits within the same region.

        This function creates a reference index, runs the alignment-based specificity filter in parallel
        across all regions in the oligo database, and returns a list of oligo pairs with hits.

        :param oligo_database: The OligoDatabase containing the oligonucleotides and their associated attributes.
        :type oligo_database: OligoDatabase
        :param sequence_type: The type of sequence to be used for filter calculations.
        :type sequence_type: _TYPES_SEQ["oligo", "target"]
        :param n_jobs: The number of parallel jobs to use for processing.
        :type n_jobs: int
        :return: A list of tuples representing oligo pairs that have significant hits.
        :rtype: list
        """
        self.sequence_type = sequence_type

        # when getting oligo pair hits we want to consider hits within the same region
        consider_hits_from_input_region = True

        file_reference = self.create_reference(n_jobs=n_jobs)

        region_ids = list(oligo_database.database.keys())
        name = " ".join(string.capitalize() for string in self.filter_name.split("_"))
        with joblib_progress(description=f"Specificity Filter: {name}", total=len(region_ids)):
            table_hits = Parallel(n_jobs=n_jobs, prefer="threads", require="sharedmem")(
                delayed(self._run_filter)(
                    region_id=region_id,
                    oligo_database=oligo_database,
                    file_reference=file_reference,
                    consider_hits_from_input_region=consider_hits_from_input_region,
                    mode=2,
                )
                for region_id in region_ids
            )

        table_hits = pd.concat(table_hits, ignore_index=True)
        oligo_pair_hits = list(zip(table_hits["query"].values, table_hits["reference"].values))

        self.remove_reference(file_reference)

        return oligo_pair_hits

    def _run_filter(
        self,
        region_id: str,
        oligo_database: OligoDatabase,
        file_reference: str,
        consider_hits_from_input_region: bool,
        mode: int,
    ) -> pd.DataFrame:
        """
        Executes the filtering process for a specific region by running the search and identifying significant hits.

        :param region_id: Region ID to process.
        :type region_id: str
        :param oligo_database: The OligoDatabase containing the oligonucleotides and their associated attributes.
        :type oligo_database: OligoDatabase
        :param file_reference: Path to the reference file used for alignment filtering.
        :type file_reference: str
        :param consider_hits_from_input_region: Flag to consider hits within the same input region.
        :type consider_hits_from_input_region: bool
        :param mode: Operation mode — 0: flag hits, 1: remove hits, 2: return hits table.
        :type mode: int
        :return: DataFrame of hit results for the region if mode is 2, otherwise None.
        :rtype: pd.DataFrame
        """
        search_results = self._run_search(
            oligo_database=oligo_database,
            file_reference=file_reference,
            region_id=region_id,
        )
        table_hits = self._find_hits(
            oligo_database=oligo_database,
            search_results=search_results,
            consider_hits_from_input_region=consider_hits_from_input_region,
            region_id=region_id,
        )

        if mode == 0:
            oligos_with_hits_region = {region_id: table_hits["query"].unique()}
            oligos_with_hits_region_attributes = (
                table_hits.groupby("query")["reference"].apply(list).to_dict()
            )
            self._flag_hits_in_database(
                oligo_database=oligo_database,
                region_ids=region_id,
                oligos_with_hits=oligos_with_hits_region,
                oligos_with_hits_attributes=oligos_with_hits_region_attributes,
            )
            return None
        elif mode == 1:
            oligos_with_hits_region = {region_id: table_hits["query"].unique()}
            self._filter_hits_from_database(
                oligo_database=oligo_database,
                region_ids=region_id,
                oligos_with_hits=oligos_with_hits_region,
            )
            return None
        elif mode == 2:
            return table_hits
        else:
            raise ValueError(
                f"Mode {mode} not available. Choose mode=0 for removing hits from the database, mode=1 for flagging the hits in the database or mode=2 for returning the hits table."
            )

    @abstractmethod
    def _run_search(
        self,
        oligo_database: OligoDatabase,
        file_reference: str,
        region_id: str,
    ) -> pd.DataFrame:
        """
        Abstract method to run a search against a ReferenceDatabase using a specified indexed reference file.

        :param oligo_database: The OligoDatabase containing the oligonucleotides and their associated attributes.
        :type oligo_database: OligoDatabase
        :param file_reference: Path to the reference file used for alignment filtering.
        :type file_reference: str
        :param region_id: Region ID to process.
        :type region_id: str
        :return: A DataFrame containing the search results.
        :rtype: pd.DataFrame
        """

    @abstractmethod
    def _find_hits(
        self,
        oligo_database: OligoDatabase,
        search_results: pd.DataFrame,
        consider_hits_from_input_region: bool,
        region_id: str,
    ) -> pd.DataFrame:
        """
        Abstract method to identify significant hits from search results, potentially excluding hits within the same region.

        :param oligo_database: The OligoDatabase containing the oligonucleotides and their associated attributes.
        :type oligo_database: OligoDatabase
        :param search_results: DataFrame containing the raw search results.
        :type search_results: pd.DataFrame
        :param consider_hits_from_input_region: Flag to consider or ignore hits from the same region.
        :type consider_hits_from_input_region: bool
        :param region_id: Region ID to process.
        :type region_id: str
        :return: A DataFrame with the identified significant hits.
        :rtype: pd.DataFrame
        """

    def _get_queries(
        self,
        oligo_database: OligoDatabase,
        table_hits: pd.DataFrame,
        region_id: str,
    ) -> List[str]:
        """
        Retrieves the query sequences from the OligoDatabase based on the hit information.

        :param oligo_database: The OligoDatabase containing the oligonucleotides and their associated attributes.
        :type oligo_database: OligoDatabase
        :param table_hits: DataFrame containing the hits with query IDs.
        :type table_hits: pd.DataFrame
        :param region_id: Region ID to process.
        :type region_id: str
        :return: A list of query sequences corresponding to the hits.
        :rtype: List[str]
        """
        queries = [
            oligo_database.database[region_id][query_id][self.sequence_type]
            for query_id in table_hits["query"]
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
