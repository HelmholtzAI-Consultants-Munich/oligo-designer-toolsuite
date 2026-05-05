############################################
# imports
############################################

import os
from abc import ABC, abstractmethod
from pathlib import Path

import pandas as pd
from joblib import Parallel, delayed
from joblib_progress import joblib_progress

from oligo_designer_toolsuite._constants import SEPARATOR_FASTA_HEADER_FIELDS, SEPARATOR_OLIGO_ID
from oligo_designer_toolsuite._exceptions import ConfigurationError
from oligo_designer_toolsuite.database import OligoDatabase, ReferenceDatabase
from oligo_designer_toolsuite.utils import cast_to_list, remove_index_files

############################################
# Oligo Specificity Filter Classes
############################################


class BaseSpecificityFilter(ABC):
    """
    A base class for implementing specificity filters that operate on OligoDatabase.

    The `BaseSpecificityFilter` class provides the structure for creating filters
    that assess the specificity of oligonucleotides. These filters can be customized and extended
    to apply various criteria to an OligoDatabase, helping to refine and select optimal oligos for
    specific applications.

    :param filter_name: Name of the filter for identification purposes.
    :type filter_name: str
    :param dir_output: Directory path where output files will be saved.
    :type dir_output: str
    """

    def __init__(self, filter_name: str, dir_output: str) -> None:
        """Constructor for the BaseSpecificityFilter class."""
        # folder where we write the intermediate files
        self.filter_name = filter_name
        self.dir_output = os.path.abspath(os.path.join(dir_output, self.filter_name))
        Path(self.dir_output).mkdir(parents=True, exist_ok=True)

        self.sequence_type: str | None = None

    @abstractmethod
    def apply(
        self,
        oligo_database: OligoDatabase,
        sequence_type: str | None,
        n_jobs: int = 1,
    ) -> OligoDatabase:
        """
        Apply the specificity filter to the given OligoDatabase.

        This abstract method defines the interface for applying the specificity filter to the oligonucleotide sequences
        in the OligoDatabase specified by the sequence type. The implementation should return a filtered version of the
        OligoDatabase.

        :param oligo_database: The OligoDatabase instance containing oligonucleotide sequences and their associated properties. This database stores oligo data organized by genomic regions and can be used for filtering, property calculations, set generation, and output operations.
        :type oligo_database: OligoDatabase
        :param sequence_type: Type of sequence being processed.
        :type sequence_type: str | None
        :param n_jobs: Number of parallel jobs to use for processing.
        :type n_jobs: int
        :return: The filtered OligoDatabase.
        :rtype: OligoDatabase
        """

    def _filter_hits_from_database(
        self,
        oligo_database: OligoDatabase,
        region_ids: str | list[str],
        oligos_with_hits: dict,
    ) -> None:
        """
        Remove oligonucleotides with hits from the database for a specific region.

        This method iterates over oligonucleotides in a specified region and removes any that have been identified
        as having hits (i.e., matches) based on the filtering criteria.

        :param oligo_database: The OligoDatabase instance containing oligonucleotide sequences and their associated properties. This database stores oligo data organized by genomic regions and can be used for filtering, property calculations, set generation, and output operations.
        :type oligo_database: OligoDatabase
        :param region_ids: Region identifier(s) to process. Can be a single region ID (str) or a list of region IDs (List[str]). If None, all regions in the database are processed, defaults to None.
        :type region_ids: Union[str, List[str]], optional
        :param oligos_with_hits: A dictionary of regio_ids associated with oligo_ids that have been identified as hits and should be removed.
        :type oligos_with_hits: dict
        """
        region_ids = cast_to_list(region_ids) if region_ids else oligo_database.database.keys()

        for region_id in region_ids:
            oligo_ids = list(oligo_database.database[region_id].keys())
            for oligo_id in oligo_ids:
                if oligo_id in oligos_with_hits[region_id]:
                    del oligo_database.database[region_id][oligo_id]

    def _flag_hits_in_database(
        self,
        oligo_database: OligoDatabase,
        region_ids: str | list[str],
        oligos_with_hits: dict,
        oligos_with_hits_properties: dict,
    ) -> None:
        """
        Flags oligos in the database based on whether they have hits, and assigns hit properties.

        This method iterates over oligonucleotides in a specified region and flags each oligo based on
        whether it appears in the list of hits. If an oligo has a hit, the corresponding reference is assigned.
        Otherwise, the property is set to None.

        :param oligo_database: The OligoDatabase instance containing oligonucleotide sequences and their associated properties. This database stores oligo data organized by genomic regions and can be used for filtering, property calculations, set generation, and output operations.
        :type oligo_database: OligoDatabase
        :param region_ids: Region identifier(s) to process. Can be a single region ID (str) or a list of region IDs (List[str]). If None, all regions in the database are processed, defaults to None.
        :type region_ids: Union[str, List[str]], optional
        :param oligos_with_hits: A dictionary of regio_ids associated with oligo_ids that have been identified as hits and should be removed.
        :type oligos_with_hits: dict
        :param oligos_with_hits_properties: Dictionary mapping oligo IDs to properties related to their hits.
        :type oligos_with_hits_properties: dict
        """
        region_ids = cast_to_list(region_ids) if region_ids else oligo_database.database.keys()

        for region_id in region_ids:
            oligo_ids = list(oligo_database.database[region_id].keys())
            for oligo_id in oligo_ids:
                if oligo_id in oligos_with_hits[region_id]:
                    oligo_database.database[region_id][oligo_id][self.filter_name] = (
                        oligos_with_hits_properties[oligo_id]
                    )
                else:
                    oligo_database.database[region_id][oligo_id][self.filter_name] = None


class ReferenceSpecificityFilter(BaseSpecificityFilter):
    """
    A base class for implementing specificity filters using a reference database.

    The `ReferenceSpecificityFilter` class provides a framework for developing filters that
    assess the potential off-target effects of oligonucleotides wrt reference sequences.

    :param remove_hits: If True, oligos overlapping variants are removed. If False, they are flagged.
    :type remove_hits: bool
    :param filter_name: Name of the filter for identification purposes.
    :type filter_name: str
    :param dir_output: Directory path where output files will be saved.
    :type dir_output: str
    """

    def __init__(
        self,
        remove_hits: bool,
        filter_name: str,
        dir_output: str,
    ) -> None:
        """Constructor for the ReferenceSpecificityFilter class."""
        # folder where we write the intermediate files
        self.filter_name = filter_name
        self.dir_output = os.path.abspath(os.path.join(dir_output, self.filter_name))
        Path(self.dir_output).mkdir(parents=True, exist_ok=True)

        self.remove_hits = remove_hits
        self.reference_database: ReferenceDatabase | None = None
        self.sequence_type: str | None = None

    def set_reference_database(self, reference_database: ReferenceDatabase | None) -> None:
        """
        Set the ReferenceDatabase for reference based specificity filters.

        :param reference_database: The ReferenceDatabase instance containing reference sequences for alignment or comparison.
        :type reference_database: ReferenceDatabase | None
        """
        self.reference_database = reference_database

    @abstractmethod
    def _create_reference(self, n_jobs: int) -> str:
        """
        Abstract method to write a reference database to file and create an index in case of alignment based methods.

        :param n_jobs: Number of parallel jobs to use for processing.
        :type n_jobs: int
        :return: The name of the created reference file.
        :rtype: str
        """

    def _remove_reference(self, file_reference: str) -> None:
        """
        Removes the reference files created for the filter.

        This method removes the reference file and all associated index files.
        For FASTA files, it removes .fai index files.
        For VCF files, it removes .csi and .tbi index files.
        For BLAST databases, it removes .nhr, .nin, .nsq index files.
        For Bowtie/Bowtie2 indexes, it removes .ebwt and .bt2 index files.

        :param file_reference: The base name of the reference files to be removed.
        :type file_reference: str
        """
        remove_index_files(file_reference=file_reference, dir_output=self.dir_output)

        if os.path.exists(file_reference):
            os.remove(file_reference)

    def _read_search_output(
        self,
        file_search_results: str,
        names_search_output: list,
        usecols: list | None = None,
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


class AlignmentSpecificityFilter(ReferenceSpecificityFilter):
    """
    A base class for implementing filters that utilize sequence alignment methods to evaluate oligonucleotide specificity.

    The `AlignmentSpecificityFilter` class provides a framework for developing filters that assess the potential
    off-target effects of oligonucleotides by aligning them against reference sequences.

    :param remove_hits: If True, oligos overlapping variants are removed. If False, they are flagged.
    :type remove_hits: bool
    :param filter_name: Name of the filter for identification purposes.
    :type filter_name: str
    :param dir_output: Directory path where output files will be saved.
    :type dir_output: str
    """

    def __init__(
        self,
        remove_hits: bool,
        filter_name: str,
        dir_output: str,
    ) -> None:
        """Constructor for the AlignmentSpecificityFilter class."""

        # folder where we write the intermediate files
        self.filter_name = filter_name
        self.dir_output = os.path.abspath(os.path.join(dir_output, self.filter_name))
        Path(self.dir_output).mkdir(parents=True, exist_ok=True)

        self.remove_hits = remove_hits
        self.reference_database: ReferenceDatabase | None = None
        self.sequence_type: str | None = None

    def apply(
        self,
        oligo_database: OligoDatabase,
        sequence_type: str | None,
        n_jobs: int = 1,
    ) -> OligoDatabase:
        """
        Applies the alignment-based specificity filter to an OligoDatabase, filtering out sequences with off-target hits.

        This function creates a reference index, runs the alignment-based specificity filter in parallel
        across all regions in the oligo database, and removes or flags oligos with off-target hits.

        :param oligo_database: The OligoDatabase instance containing oligonucleotide sequences and their associated properties. This database stores oligo data organized by genomic regions and can be used for filtering, property calculations, set generation, and output operations.
        :type oligo_database: OligoDatabase
        :param sequence_type: Type of sequence being processed.
        :type sequence_type: str | None
        :param n_jobs: Number of parallel jobs to use for processing.
        :type n_jobs: int
        :return: The filtered OligoDatabase.
        :rtype: OligoDatabase
        """
        self.sequence_type = sequence_type

        # when applying filters we don't want to consider hits within the same region
        consider_hits_from_input_region = False

        file_reference = self._create_reference(n_jobs=n_jobs)

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

        self._remove_reference(file_reference)

        return oligo_database

    def get_oligo_pair_hits(
        self,
        oligo_database: OligoDatabase,
        sequence_type: str,
        n_jobs: int,
    ) -> list[tuple[str, str]]:
        """
        Retrieves pairs of oligonucleotides that have significant alignment hits within the same region.

        This function creates a reference index, runs the alignment-based specificity filter in parallel
        across all regions in the oligo database, and returns a list of oligo pairs with hits.

        :param oligo_database: The OligoDatabase instance containing oligonucleotide sequences and their associated properties. This database stores oligo data organized by genomic regions and can be used for filtering, property calculations, set generation, and output operations.
        :type oligo_database: OligoDatabase
        :param sequence_type: Type of sequence being processed.
        :type sequence_type: str
        :param n_jobs: Number of parallel jobs to use for processing.
        :type n_jobs: int
        :return: A list of tuples representing oligo pairs that have significant hits.
        :rtype: list[tuple[str, str]]
        """
        self.sequence_type = sequence_type

        # when getting oligo pair hits we want to consider hits within the same region
        consider_hits_from_input_region = True

        file_reference = self._create_reference(n_jobs=n_jobs)

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

        self._remove_reference(file_reference)

        return oligo_pair_hits

    def _run_filter(
        self,
        region_id: str,
        oligo_database: OligoDatabase,
        file_reference: str,
        consider_hits_from_input_region: bool,
        mode: int,
    ) -> pd.DataFrame | None:
        """
        Executes the filtering process for a specific region by running the search and identifying significant hits.

        :param region_id: Region ID to process.
        :type region_id: str
        :param oligo_database: The OligoDatabase instance containing oligonucleotide sequences and their associated properties. This database stores oligo data organized by genomic regions and can be used for filtering, property calculations, set generation, and output operations.
        :type oligo_database: OligoDatabase
        :param file_reference: Path to the reference file used for alignment filtering.
        :type file_reference: str
        :param consider_hits_from_input_region: Flag to consider hits within the same input region.
        :type consider_hits_from_input_region: bool
        :param mode: Operation mode — 0: flag hits, 1: remove hits, 2: return hits table.
        :type mode: int
        :return: DataFrame of hit results for the region if mode is 2, otherwise None.
        :rtype: pd.DataFrame | None
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
            oligos_with_hits_region_properties = (
                table_hits.groupby("query")["reference"].apply(list).to_dict()
            )
            self._flag_hits_in_database(
                oligo_database=oligo_database,
                region_ids=region_id,
                oligos_with_hits=oligos_with_hits_region,
                oligos_with_hits_properties=oligos_with_hits_region_properties,
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
            raise ConfigurationError(
                f"Mode '{mode}' is not available. Choose mode=0 for removing hits from the database, "
                f"mode=1 for flagging the hits in the database, or mode=2 for returning the hits table."
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

        :param oligo_database: The OligoDatabase instance containing oligonucleotide sequences and their associated properties. This database stores oligo data organized by genomic regions and can be used for filtering, property calculations, set generation, and output operations.
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

        :param oligo_database: The OligoDatabase instance containing oligonucleotide sequences and their associated properties. This database stores oligo data organized by genomic regions and can be used for filtering, property calculations, set generation, and output operations.
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
    ) -> list[str]:
        """
        Retrieves the query sequences from the OligoDatabase based on the hit information.

        :param oligo_database: The OligoDatabase instance containing oligonucleotide sequences and their associated properties. This database stores oligo data organized by genomic regions and can be used for filtering, property calculations, set generation, and output operations.
        :type oligo_database: OligoDatabase
        :param table_hits: DataFrame containing the hits with query IDs.
        :type table_hits: pd.DataFrame
        :param region_id: Region ID to process.
        :type region_id: str
        :return: A list of query sequences corresponding to the hits.
        :rtype: list[str]
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
    ) -> tuple[list, list]:
        """
        Abstract method to align query and reference sequences by adding gaps.

        :param table_hits: DataFrame containing the hits with alignment information.
        :type table_hits: pd.DataFrame
        :param queries: List of query sequences.
        :type queries: list
        :param references: List of reference sequences.
        :type references: list
        :return: A tuple containing lists of gapped query and reference sequences.
        :rtype: tuple[list, list]
        """
