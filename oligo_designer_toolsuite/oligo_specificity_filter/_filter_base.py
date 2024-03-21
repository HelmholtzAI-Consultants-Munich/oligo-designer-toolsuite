############################################
# imports
############################################

import os
import re
from abc import ABC, abstractmethod
from pathlib import Path
from typing import List, Union, get_args
import pandas as pd
import numpy as np

from pathlib import Path
from joblib import Parallel, delayed
from abc import ABC, abstractmethod
from typing import get_args, List
from Bio import SeqIO, Seq
from odt_ai_filters.api import generate_api

import pandas as pd
from joblib import Parallel, delayed

from .._constants import _TYPES_SEQ, SEPARATOR_FASTA_HEADER_FIELDS, SEPARATOR_OLIGO_ID
from ..database import OligoDatabase, ReferenceDatabase

#TODO: add imports for specificity filters here, e.g. from . import specificity_filter_1, specificity_filter_2, ...


############################################
# Oligo Specificity Filter Classes
############################################


class SpecificityFilterBase(ABC):
    """An abstract base class for implementing specificity filters for oligonucleotide sequences.
    Specificity filters are used to ensure that oligos do not have unintended matches within a given
    reference database. This base class provides a common structure for such filters, including an output
    directory setup and an abstract method for applying the filter to an oligo database.

    :param dir_output: The directory where intermediate files will be saved. Defaults to "output".
    :type dir_output: str
    """

    def __init__(self, dir_output: str = "output"):
        """Constructor for the SpecificityFilterBase class."""

        self.dir_output = dir_output
        Path(self.dir_output).mkdir(parents=True, exist_ok=True)

    @abstractmethod
    def apply(
        self,
        sequence_type: _TYPES_SEQ,
        database: OligoDatabase,
        n_jobs: int,
        reference_database: ReferenceDatabase = None,
    ):
        """Abstract method to apply the specificity filter to an oligo database.

        :param sequence_type: The type of sequences being filtered, must be one of the predefined sequence types.
        :type sequence_type: _TYPES_SEQ
        :param database: The oligo database to which the filter will be applied.
        :type database: OligoDatabase
        :param n_jobs: The number of parallel jobs to run.
        :type n_jobs: int
        :param reference_database: The reference database to compare against for specificity.
            For non-alignment based specificity filter reference_database is not used, i.e. set to None.
        :type reference_database: ReferenceDatabase, optional
        """

    def _filter_hits_from_database(
        self, database_region: dict, oligos_with_hits: list[str]
    ):
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

    :param dir_output: Directory for saving intermediate files generated during the filtering process.
    :type dir_output: str
    :param ai_filter: The machine learning model used to filter the oligos {None, 'hybridization_probability'}, defaults to None
    :type ai_filter: str, optional
    :param ai_filter_thershold: The threshold below which the oligos are filtered, defaults to None
    :type ai_filter_thershold: float, optional
    :param ai_filter_path: The path to the machine learning model used to filter the oligos, if None the pretrained model provided will be used, defaults to None
    :type ai_filter_path: str, optional
    """

    def __init__(
            self, 
            dir_output: str = "output",
            ai_filter: str = None,
            ai_filter_threshold: float = None,
            ai_filter_path: str = None,
    ):
        """Construnctor for the AlignmentSpecificityFilter class."""
        # folder where we write the intermediate files
        self.dir_output = dir_output
        Path(self.dir_output).mkdir(parents=True, exist_ok=True)

        if ai_filter is not None and ai_filter_threshold is  None:
            raise ValueError(
                "ai_filter_threshold must be specified if ai_filter is specified."
            )
        # create ai API
        self.ai_filter_threshold  = ai_filter_threshold
        if ai_filter is None:
            self.ai_filter_api = None
        else:
            self.ai_filter_api = generate_api(
                ai_filter=ai_filter, ai_filter_path=ai_filter_path
            )
        self.filtered = {}

    def apply(
        self,
        sequence_type: _TYPES_SEQ,
        oligo_database: OligoDatabase,
        n_jobs: int,
        reference_database: ReferenceDatabase,
    ):
        """Applies the alignment-based specificity filter to an oligonucleotide database. TODO add documentation about ai filters

        :param sequence_type: The type of sequences being filtered, must be one of the predefined sequence types.
        :type sequence_type: _TYPES_SEQ
        :param database: The oligo database to which the filter will be applied.
        :type database: OligoDatabase
        :param n_jobs: The number of parallel jobs to run.
        :type n_jobs: int
        :param reference_database: The reference database to compare against for specificity.
        :type reference_database: ReferenceDatabase
        :return: The filtered oligo database with sequences having significant hits removed.
        :rtype: OligoDatabase

        """
        options = get_args(_TYPES_SEQ)
        assert (
            sequence_type in options
        ), f"Sequence type not supported! '{sequence_type}' is not in {options}."

        # Create index file for search
        file_reference = reference_database.write_database_to_fasta(
            filename="reference_db"
        )
        file_index = self._create_index(
            file_reference=file_reference, n_jobs=n_jobs
        )

        # Run search for each region in parallel
        region_ids = list(oligo_database.database.keys())
        table_hits = Parallel(n_jobs=n_jobs)(
            delayed(self._run_filter)(
                sequence_type=sequence_type,
                region_id=region_id,
                oligo_database=oligo_database,
                file_index=file_index,
                file_reference=file_reference,
                consider_hits_from_input_region=True,
            )
            for region_id in region_ids
        )
        # delete the temporary files 
        os.remove(file_reference)
        file_index_basename = os.path.basename(file_index)
        regex = re.compile(file_index_basename + "\..*")
        for root, _,files in os.walk(self.dir_output):
            for file in files:
                if regex.match(file):
                    os.remove(os.path.join(root, file))

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
        n_jobs: int,
        reference_database: ReferenceDatabase,
    ):
        """Identifies oligonucleotide pairs with significant hits in the reference database.

        :param sequence_type: The type of sequences being filtered, must be one of the predefined sequence types.
        :type sequence_type: _TYPES_SEQ
        :param database: The oligo database to which the filter will be applied.
        :type database: OligoDatabase
        :param n_jobs: The number of parallel jobs to run.
        :type n_jobs: int
        :param reference_database: The reference database to compare against for specificity.
        :type reference_database: ReferenceDatabase
        :return: List of oligo pairs with hits in the reference database.
        :rtype: list[tuple]
        """
        options = get_args(_TYPES_SEQ)
        assert (
            sequence_type in options
        ), f"Sequence type not supported! '{sequence_type}' is not in {options}."

        # Create index file for search
        file_reference = reference_database.write_database_to_fasta(
            filename="reference_db"
        )
        file_index = self._create_index(
            file_reference=file_reference, n_jobs=n_jobs
        )

        # Run search for each region in parallel
        region_ids = list(oligo_database.database.keys())
        table_hits = Parallel(n_jobs=n_jobs)(
            delayed(self._run_filter)(
                sequence_type=sequence_type,
                region_id=region_id,
                oligo_database=oligo_database,
                file_index=file_index,
                file_reference=file_reference,
                consider_hits_from_input_region=True,
            )
            for region_id in region_ids
        )
        os.remove(file_reference)
        file_index_basename = os.path.basename(file_index)
        regex = re.compile(file_index_basename + "\..*")
        for root, _,files in os.walk(self.dir_output):
            for file in files:
                if regex.match(file):
                    os.remove(os.path.join(root, file))

        table_hits = pd.concat(table_hits, ignore_index=True)
        oligo_pair_hits = list(zip(table_hits["query"].values, table_hits["reference"].values))

        return oligo_pair_hits

    @abstractmethod
    def _create_index(self, file_reference: str, n_jobs: int):
        """Abstract method to create an index of the reference database for alignment.

        :param file_reference: Path to the reference database fasta file.
        :type file_reference: str
        :param n_jobs: The number of parallel jobs to run.
        :type n_jobs: int
        """

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
        :param database: The oligonucleotide database to search against.
        :type database: OligoDatabase
        :param region_ids: Identifiers for the regions within the database to be searched.
        :type region_ids: str
        :param file_index: Path to the index file of the reference database.
        :type file_index: str
        """

    def _read_search_output(
        self, file_search_results: str, names_search_output: list, usecols: list = None
    ):
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

        search_results["query_region_id"] = (
            search_results["query"].str.split(SEPARATOR_OLIGO_ID).str[0]
        )
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

    def _run_filter(
        self,
        sequence_type: _TYPES_SEQ,
        region_id: str,
        oligo_database: OligoDatabase,
        file_index: str,
        file_reference: str,
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
        :param file_reference: Path to the reference database fasta file.
        :type file_reference: str
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
        if self.ai_filter_api is not None:
            table_hits = self._ai_filter_hits(
                sequence_type=sequence_type,
                table_hits=table_hits,
                file_reference=file_reference,
                oligo_database=oligo_database,
                region_id=region_id,
            )

        return table_hits

    def _ai_filter_hits(
        self,
        sequence_type: _TYPES_SEQ,
        table_hits: pd.DataFrame,
        file_reference: str,
        oligo_database: OligoDatabase,
        region_id: str,
    ) -> pd.DataFrame:
        """Filters the hits from a search operation using Machine Learning models. The Hits that recieve a score form the machine learning model lower that the given threshold are filtered out.

        :param sequence_type: The type of sequences being filtered, must be one of the predefined sequence types.
        :type sequence_type: _TYPES_SEQ
        :param table_hits: Dataframe containing the strue hits of the search results.
        :type table_hits: pd.DataFrame
        :param file_reference: Path to the fasta file used as reference for the search.
        :type file_reference: str
        :param oligo_database: The oligo database to which the filter will be applied.
        :type oligo_database: OligoDatabase
        :param region_id: The identifier for the region within the database to filter.
        :type region_id: str
        :return: Dataframe containing the filtered true hits
        :rtype: pd.DataFrame
        """

        required_fields = [
            
        ]
        # check if there are any oligos to filter
        if len(table_hits) == 0:
            return table_hits
        
        # generate the references and queries sequences
        references = self._get_references(table_hits, file_reference, region_id)
        queries = self._get_queries(sequence_type, table_hits, oligo_database, region_id)
        assert len(queries) == len(
            references
        ), "The reference sequences haven't been correctly retrieved."

        # align the references and queries by adding gaps
        gapped_queries, gapped_references = self._add_alignement_gaps(
            table_hits=table_hits, queries=queries, references=references
        )

        # predict the scores for each hit
        predictions = self.ai_filter_api.predict(
            queries=queries,
            gapped_queries=gapped_queries,
            references=references,
            gapped_references=gapped_references,
        )

        # filter the database, keep only the oligos above the threshold
        table_hits.reset_index(drop=True, inplace=True)
        ids_vanilla = len(table_hits["query"].unique())
        below_threshold = np.where(predictions < self.ai_filter_threshold)[0]
        table_hits.drop(index=below_threshold, inplace=True)
        ids_ai_filter = len(table_hits["query"].unique())
        self.filtered[region_id] = [
            ids_vanilla,
            ids_ai_filter,
            ids_vanilla - ids_ai_filter,
        ]

        return table_hits
    
    def _get_queries(self, sequence_type, table_hits: pd.DataFrame, oligo_database: OligoDatabase, region_id: str)->List[Seq.Seq]:
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
            oligo_database.database[region_id][query_id][sequence_type]
            for query_id in table_hits["query"]
        ]
        return queries
    
    @abstractmethod
    def _get_references(self, table_hits: pd.DataFrame, file_reference: str , region_id: str):
        """Abstract method to retrieve the reference sequences from the search results.

        :param table_hits: Dataframe containing the search results.
        :type table_hits: pd.DataFrame
        :param file_reference: Path to the fasta file used as reference for the search.
        :type file_reference: str
        :param region_id: The identifier for the region within the database to filter.
        :type region_id: str
        """

    @abstractmethod
    def _add_alignement_gaps(self, table_hits: pd.DataFrame, queries: List[Seq.Seq], references: List[Seq.Seq]):
        """Abstract method to add gaps to the references and queries sequences.

        :param table_hits: Dataframe containing the search results.
        :type table_hits: pd.DataFrame
        :param queries: List of the queries sequences.
        :type queries: List[Seq.Seq]
        :param references: List of the references sequences.
        :type references: List[Seq.Seq]
        """

    
        
