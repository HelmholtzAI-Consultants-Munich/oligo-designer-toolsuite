############################################
# imports
############################################

import iteration_utilities
import pandas as pd

from typing import List, get_args
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
    """A filter that identifies and removes oligonucleotides with exact match sequences within the database,
    to prevent potential off-target effects. It leverages parallel processing to efficiently analyze large oligo databases.

    :param policy: The filter policy to apply for minimizing cross-hybridization.
        If no policy is provided (i.e. policy = None) the RemoveAllPolicy() is applied. Defaults to None.
    """

    def __init__(self, policy: FilterPolicyBase = None):
        """Constructor for the ExactMatches class."""
        if policy:
            self.policy = policy
        else:
            self.policy = RemoveAllPolicy()

    def apply(
        self,
        sequence_type: _TYPES_SEQ,
        oligo_database: OligoDatabase,
        reference_database: ReferenceDatabase = None,
        n_jobs: int = 1,
    ):
        """Applies the exact match filter to an oligonucleotide database.

        :param sequence_type: The type of sequences being filtered, must be one of the predefined sequence types.
        :type sequence_type: _TYPES_SEQ
        :param database: The oligo database to which the filter will be applied.
        :type database: OligoDatabase
        :param reference_database: The reference database to compare against for specificity.
            For non-alignment based specificity filter reference_database is not used, i.e. set to None, defaults to None.
        :type reference_database: ReferenceDatabase, optional
        :param n_jobs: The number of parallel jobs to run.
        :type n_jobs: int
        :return: The filtered oligo database.
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
        with joblib_progress(description="Exact Matches", total=len(region_ids)):
            table_hits = Parallel(n_jobs=n_jobs, prefer="threads")(
                delayed(self._run_filter)(
                    sequence_type=sequence_type,
                    region_id=region_id,
                    oligo_database=oligo_database,
                    search_results=search_results,
                    sequence_oligoids_mapping=sequence_oligoids_mapping,
                    consider_hits_from_input_region=False,
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
        sequence_type: _TYPES_SEQ,
        oligo_database: OligoDatabase,
        reference_database: ReferenceDatabase = None,  # not used in this filter but needed for API
        n_jobs: int = 1,
    ):
        """Retrieves pairs of oligonucleotides with exact matches within the oligo database.
        Here we match the sequenecs to their reverse complements as a basis for the cross-hybridization filter.

        :param sequence_type: The type of sequences being filtered, must be one of the predefined sequence types.
        :type sequence_type: _TYPES_SEQ
        :param database: The oligo database to which the filter will be applied.
        :type database: OligoDatabase
        :param reference_database: The reference database to compare against for specificity.
            Not used in this filter.
        :type reference_database: ReferenceDatabase
        :param n_jobs: The number of parallel jobs to run.
        :type n_jobs: int
        :return: List of oligo pairs with hits in the reference database.
        :rtype: list[tuple]
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
        with joblib_progress(description="Exact Matches", total=len(region_ids)):
            table_hits = Parallel(n_jobs=n_jobs, prefer="threads")(
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

    def _get_duplicated_sequences(self, sequences: list):
        """Identifies duplicated sequences within an oligonucleotide database for a specified sequence type.

        :param sequence_type: The type of sequences to analyze within the database.
        :type sequence_type: _TYPES_SEQ
        :param oligo_database: The oligonucleotide database being analyzed.
        :type oligo_database: OligoDatabase
        :return: A list of duplicated sequences found in the database.
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
        sequence_type: _TYPES_SEQ,
        region_id: str,
        oligo_database: OligoDatabase,
        search_results: List,
        sequence_oligoids_mapping: dict,
        consider_hits_from_input_region: bool,
    ):
        """Filters oligonucleotides based on sequence duplications, identifying hits within and across regions.

        :param sequence_type: The type of sequences to filter.
        :type sequence_type: _TYPES_SEQ
        :param region_id: The region within the oligo database to apply the filter.
        :type region_id: str
        :param oligo_database: The oligo database containing sequences to be filtered.
        :type oligo_database: OligoDatabase
        :param search_results: A list of sequences identified as duplicates.
        :type search_results: List
        :param sequence_oligoids_mapping: Mapping from sequences to oligo IDs.
        :type sequence_oligoids_mapping: dict
        :param consider_hits_from_input_region: Flag to indicate whether hits from the input region should be considered.
        :type consider_hits_from_input_region: bool
        :return: A tuple containing a table of hits and a list of oligos with those hits.
        :rtype: (pd.DataFrame, list)
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
