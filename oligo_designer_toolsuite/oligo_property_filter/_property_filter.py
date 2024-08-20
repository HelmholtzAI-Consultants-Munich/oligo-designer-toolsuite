############################################
# imports
############################################

from typing import get_args
from Bio.SeqUtils import Seq

from joblib import Parallel, delayed
from joblib_progress import joblib_progress

from oligo_designer_toolsuite._constants import _TYPES_SEQ
from oligo_designer_toolsuite.database import OligoDatabase
from oligo_designer_toolsuite.oligo_property_filter import PropertyFilterBase

############################################
# Property Filter Class
############################################


class PropertyFilter:
    """
    A class for applying multiple property filters to sequences in an oligo database.

    The `PropertyFilter` class allows you to apply a list of sequence filters (subclasses of `PropertyFilterBase`) to an oligo database.
    The filters are applied in parallel across all regions of the database, and sequences that do not meet all filter criteria are removed.

    :param filters: A list of property filters to apply to sequences.
    :type filters: list[PropertyFilterBase]
    """

    def __init__(
        self,
        filters: list[PropertyFilterBase],
    ) -> None:
        """Constructor for the PropertyFilter class."""
        self.filters = filters

    def apply(
        self, oligo_database: OligoDatabase, sequence_type: _TYPES_SEQ, n_jobs: int = 1
    ) -> OligoDatabase:
        """
        Apply the property filters to all sequences in the oligo database and filter
        sequences in the oligo database based on the specified property filters.
        Sequences that do not meet the criteria of all filters are removed.

        :param oligo_database: The Oligo Database containing the oligonucleotides and their associated attributes.
        :type oligo_database: OligoDatabase
        :param sequence_type: The type of sequence to be used for filter calculations.
        :type sequence_type: _TYPES_SEQ
        :param n_jobs: The number of jobs to run in parallel, default is 1.
        :type n_jobs: int
        :return: The filtered oligo database.
        :rtype: OligoDatabase
        """
        options = get_args(_TYPES_SEQ)
        assert (
            sequence_type in options
        ), f"Sequence type not supported! '{sequence_type}' is not in {options}."

        region_ids = list(oligo_database.database.keys())
        with joblib_progress(description="Property Filter", total=len(region_ids)):
            Parallel(n_jobs=n_jobs, prefer="threads", require="sharedmem")(
                delayed(self._filter_region)(oligo_database, region_id, sequence_type)
                for region_id in region_ids
            )

        oligo_database.remove_regions_with_insufficient_oligos(pipeline_step="Property Filters")
        return oligo_database

    def _filter_region(
        self, oligo_database: OligoDatabase, region_id: str, sequence_type: _TYPES_SEQ
    ) -> None:

        oligo_ids = list(oligo_database.database[region_id].keys())
        for oligo_id in oligo_ids:
            fulfills_all_filter = self._filter_sequence(
                oligo_database.database[region_id][oligo_id][sequence_type]
            )
            if not fulfills_all_filter:
                del oligo_database.database[region_id][oligo_id]

    def _filter_sequence(self, sequence: Seq) -> bool:

        for filter in self.filters:
            fulfills_filter = filter.apply(sequence)
            if not fulfills_filter:  # stop at the first false we obtain
                return False
        return True
