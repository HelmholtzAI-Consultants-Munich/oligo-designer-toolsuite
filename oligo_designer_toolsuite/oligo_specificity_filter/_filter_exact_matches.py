############################################
# imports
############################################

from typing import List, get_args

import iteration_utilities
import pandas as pd
from joblib import Parallel, delayed
from joblib_progress import joblib_progress

from oligo_designer_toolsuite._constants import _TYPES_SEQ
from oligo_designer_toolsuite.database import OligoDatabase, ReferenceDatabase
from oligo_designer_toolsuite.oligo_specificity_filter import SpecificityFilterBase

from ._policies import FilterPolicyBase, RemoveAllPolicy

############################################
# Oligo Exact Match Filter Classes
############################################


class ExactMatchFilter(SpecificityFilterBase):
    """
    A filter class for identifying and handling exact sequence matches within an OligoDatabase.

    The `ExactMatchFilter` class is designed to apply a specific policy to oligonucleotides that have exact sequence matches in the OligoDatabase.
    This filter can be used to remove or handle exact matches based on the provided policy.

    :param policy: The policy to apply to exact matches. If not provided, a default policy (`RemoveAllPolicy`) will be used.
    :type policy: FilterPolicyBase
    :param filter_name: Name of the filter for identification purposes.
    :type filter_name: str
    """

    def __init__(
        self,
        policy: FilterPolicyBase = None,
        filter_name: str = "exact_match_filter",
    ) -> None:
        """Constructor for the ExactMatches class."""
        if not policy:
            policy = RemoveAllPolicy()
        self.policy = policy
        self.filter_name = filter_name

    def apply(
        self,
        sequence_type: _TYPES_SEQ,
        oligo_database: OligoDatabase,
        reference_database: ReferenceDatabase = None,  # not utilized in this filter
        n_jobs: int = 1,
    ) -> OligoDatabase:
        """
        Applies the exact match filter to an OligoDatabase.

        This function identifies sequences within the OligoDatabase that have exact matches according to the specified `sequence_type`.
        It then uses the provided policy to determine how these exact matches should be handled.
        The filter operates in parallel across regions in the OligoDatabase to improve performance.

        :param sequence_type: The type of sequence to be used for the filter calculations.
        :type sequence_type: _TYPES_SEQ["oligo", "target"]
        :param oligo_database: The OligoDatabase containing the oligonucleotides and their associated attributes.
        :type oligo_database: OligoDatabase
        :param reference_database: The ReferenceDatabase used for alignment (not utilized in this filter).
        :type reference_database: ReferenceDatabase
        :param n_jobs: The number of parallel jobs to use for processing.
        :type n_jobs: int
        :return: The filtered OligoDatabase with exact matching sequences removed.
        :rtype: OligoDatabase
        """
        # extract all the sequences
        sequences = oligo_database.get_sequence_list(sequence_type=sequence_type)
        search_results = self._get_duplicated_sequences(sequences)

        # get mapping from sequence to oligo_id
        # use same sequence as sequence_type
        sequence_oligoids_mapping = oligo_database.get_sequence_oligoid_mapping(
            sequence_type=sequence_type, sequence_to_upper=True
        )

        region_ids = list(oligo_database.database.keys())
        name = " ".join(string.capitalize() for string in self.filter_name.split("_"))
        with joblib_progress(description=f"Specificity Filter: {name}", total=len(region_ids)):
            table_hits = Parallel(n_jobs=n_jobs, prefer="threads", require="sharedmem")(
                delayed(self._run_filter)(
                    oligo_database=oligo_database,
                    search_results=search_results,
                    sequence_oligoids_mapping=sequence_oligoids_mapping,
                    consider_hits_from_input_region=False,
                    sequence_type=sequence_type,
                    region_id=region_id,
                )
                for region_id in region_ids
            )

        table_hits = pd.concat(table_hits, ignore_index=True)
        oligo_pair_hits = list(zip(table_hits["query"].values, table_hits["reference"].values))
        oligos_with_hits = self.policy.apply(oligo_pair_hits=oligo_pair_hits, oligo_database=oligo_database)

        for region_id in region_ids:
            self._filter_hits_from_database(
                oligo_database=oligo_database,
                region_id=region_id,
                oligos_with_hits=oligos_with_hits[region_id],
            )
        return oligo_database

    def get_oligo_pair_hits(
        self,
        oligo_database: OligoDatabase,
        reference_database: ReferenceDatabase,  # not utilized in this filter
        sequence_type: _TYPES_SEQ,
        n_jobs: int = 1,
    ) -> list:
        """
        Identifies oligonucleotide pairs that have exact sequence matches within the OligoDatabase.

        This function compares oligo sequences within the OligoDatabase to identify pairs of oligos that have exact matches,
        including their reverse complements. It generates a list of these matching oligo pairs, which can be used for further filtering or analysis.

        :param oligo_database: The OligoDatabase containing the oligonucleotides and their associated attributes.
        :type oligo_database: OligoDatabase
        :param reference_database: The ReferenceDatabase used for alignment (not utilized in this filter).
        :type reference_database: ReferenceDatabase
        :param sequence_type: The type of sequence to be used for the filter calculations.
        :type sequence_type: _TYPES_SEQ["oligo", "target"]
        :param n_jobs: The number of parallel jobs to use for processing.
        :type n_jobs: int
        :return: A list of tuples representing oligo pairs that have exact matches.
        :rtype: list
        """
        options = get_args(_TYPES_SEQ)
        assert (
            sequence_type in options
        ), f"Sequence type not supported! '{sequence_type}' is not in {options}."
        sequence_type_reverse_complement = options[0] if options[0] != sequence_type else options[1]

        # extract all the sequences
        sequences = oligo_database.get_sequence_list(sequence_type=sequence_type)
        sequences_rc = oligo_database.get_sequence_list(sequence_type=sequence_type_reverse_complement)
        search_results = self._get_duplicated_sequences(sequences + sequences_rc)

        # get mapping from sequence to oligo_id
        sequence_oligoids_mapping = oligo_database.get_sequence_oligoid_mapping(
            sequence_type=sequence_type_reverse_complement, sequence_to_upper=True
        )

        region_ids = list(oligo_database.database.keys())
        name = " ".join(string.capitalize() for string in self.filter_name.split("_"))
        with joblib_progress(description=f"Specificity Filter: {name}", total=len(region_ids)):
            table_hits = Parallel(n_jobs=n_jobs, prefer="threads", require="sharedmem")(
                delayed(self._run_filter)(
                    sequence_type=sequence_type,
                    region_id=region_id,
                    oligo_database=oligo_database,
                    search_results=search_results,
                    sequence_oligoids_mapping=sequence_oligoids_mapping,
                    consider_hits_from_input_region=True,
                )
                for region_id in region_ids
            )

        # Process results
        table_hits = pd.concat(table_hits, ignore_index=True)
        oligo_pair_hits = list(zip(table_hits["query"].values, table_hits["reference"].values))

        return oligo_pair_hits

    def _get_duplicated_sequences(self, sequences: list) -> list:
        """
        Identifies duplicated sequences within a list of sequences.

        This function takes a list of sequences, converts them to uppercase, and identifies sequences that appear more than once.
        It returns a list of these duplicated sequences, which can be used for further analysis or filtering.

        :param sequences: A list of sequences to be checked for duplicates.
        :type sequences: list
        :return: A list of duplicated sequences.
        :rtype: list
        """
        # convert to upper sequence
        sequences = [sequence.upper() for sequence in sequences]

        # find the duplicates within the database
        duplicated_sequences = list(
            iteration_utilities.unique_everseen(iteration_utilities.duplicates(sequences))
        )

        return duplicated_sequences

    def _run_filter(
        self,
        oligo_database: OligoDatabase,
        search_results: List,
        sequence_oligoids_mapping: dict,
        consider_hits_from_input_region: bool,
        sequence_type: _TYPES_SEQ,
        region_id: str,
    ) -> pd.DataFrame:
        """
        Runs a filter on the OligoDatabase to identify and record any matching sequences found in the provided search results.

        This function checks each sequence in a specified region of the oligo OligoDatabase against a list of search results.
        It then maps these sequences to their corresponding oligo IDs and records any hits that are either within or outside the input region, depending on the filter settings.

        :param oligo_database: The OligoDatabase containing the oligonucleotides and their associated attributes.
        :type oligo_database: OligoDatabase
        :param search_results: A list of sequences that were identified as duplicates.
        :type search_results: List
        :param sequence_oligoids_mapping: A dictionary mapping sequences to their respective oligo IDs.
        :type sequence_oligoids_mapping: dict
        :param consider_hits_from_input_region: Flag indicating whether to consider hits from the same region as the input sequences.
        :type consider_hits_from_input_region: bool
        :param sequence_type: The type of sequence to be used for the filter calculations.
        :type sequence_type: _TYPES_SEQ["oligo", "target"]
        :param region_id: Region ID to process.
        :type region_id: str
        :return: A DataFrame containing the oligo ID pairs for sequences that matched the search results.
        :rtype: pd.DataFrame
        """
        database_region = oligo_database.database[region_id]
        hit_dict = {}
        for oligo_id in database_region.keys():
            oligo_seq = database_region[oligo_id][sequence_type].upper()
            if oligo_seq in search_results:
                # find all reverse complements with the same sequence
                if oligo_seq in sequence_oligoids_mapping:
                    if consider_hits_from_input_region:
                        hits = [oligo_id for oligo_id in sequence_oligoids_mapping[oligo_seq]]
                    else:
                        hits = [
                            oligo_id
                            for oligo_id in sequence_oligoids_mapping[oligo_seq]
                            if not oligo_id.startswith(region_id)
                        ]
                else:
                    hits = []
                # don't store when all hits come from the same region and are filtered out
                if hits:
                    hit_dict[oligo_id] = hits

        table_hits = pd.DataFrame(
            [(key, value) for key, value in hit_dict.items()],
            columns=["query", "reference"],
        )
        table_hits = table_hits.explode("reference", ignore_index=True)

        return table_hits
