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
<<<<<<< HEAD
    """A filter that identifies and removes oligonucleotides with exact match sequences within the database,
    to prevent potential off-target effects. It leverages parallel processing to efficiently analyze large oligo databases.

    :param policy: The filter policy to apply for minimizing cross-hybridization.
        If no policy is provided (i.e. policy = None) the RemoveAllPolicy() is applied. Defaults to None.
    """

=======
>>>>>>> origin/pipelines
    def __init__(
        self,
        policy: FilterPolicyBase = None,
        filter_name: str = "exact_match_filter",
<<<<<<< HEAD
    ):
=======
    ) -> None:
>>>>>>> origin/pipelines
        """Constructor for the ExactMatches class."""
        if not policy:
            policy = RemoveAllPolicy()
        self.policy = policy
        self.filter_name = filter_name

        self.filter_name = filter_name

    def apply(
        self,
        oligo_database: OligoDatabase,
        reference_database: ReferenceDatabase,  # not used in this filter but needed for API
        sequence_type: _TYPES_SEQ,
        n_jobs: int = 1,
    ) -> OligoDatabase:
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
        reference_database: ReferenceDatabase,  # not used in this filter but needed for API
        sequence_type: _TYPES_SEQ,
        n_jobs: int = 1,
    ) -> list:
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
