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
        """ """

    def _filter_hits_from_database(
        self, oligo_database: OligoDatabase, region_id: str, oligos_with_hits: list[str]
    ) -> None:
        oligo_ids = list(oligo_database.database[region_id].keys())
        for oligo_id in oligo_ids:
            if oligo_id in oligos_with_hits:
                del oligo_database.database[region_id][oligo_id]


class AlignmentSpecificityFilter(SpecificityFilterBase):
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
        region_id: List[str],
        consider_hits_from_input_region: bool,
    ) -> None:
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
        """ """

    def _run_filter(
        self,
        sequence_type: _TYPES_SEQ,
        region_id: str,
        oligo_database: OligoDatabase,
        file_index: str,
        consider_hits_from_input_region: bool,
    ) -> pd.DataFrame:
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
        """ """

    def _read_search_output(
        self, file_search_results: str, names_search_output: list, usecols: list = None
    ) -> pd.DataFrame:
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
        """ """

    def _remove_index(self, file_index: str) -> None:
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
        queries = [
            oligo_database.database[region_id][query_id][sequence_type] for query_id in table_hits["query"]
        ]
        return queries

    @abstractmethod
    def _get_references(self, table_hits: pd.DataFrame, file_reference: str, region_id: str) -> list:
        """ """

    @abstractmethod
    def _add_alignment_gaps(
        self,
        table_hits: pd.DataFrame,
        queries: list,
        references: list,
    ) -> Tuple[list, list]:
        """ """
