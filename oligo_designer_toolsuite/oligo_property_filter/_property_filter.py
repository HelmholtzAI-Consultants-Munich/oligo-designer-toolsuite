############################################
# imports
############################################

from typing import get_args

from effidict import LRUDict
from joblib import Parallel, delayed
from joblib_progress import joblib_progress

from oligo_designer_toolsuite._constants import _TYPES_SEQ
from oligo_designer_toolsuite.database import OligoDatabase
from oligo_designer_toolsuite.oligo_property_filter import PropertyFilterBase

############################################
# Property Filter Class
############################################


class PropertyFilter:
    """A class representing a collection of filters to be applied on oligo databases.

    This class manages the application of multiple PropertyFilterBase instances on an oligo database,
    allowing for complex filtering strategies. Oligos which don't fulfill the filtering criteria are
    removed from the database.

    :param filters: A list of filter objects derived from PropertyFilterBase.
    :type filters: list[PropertyFilterBase]
    """

    def __init__(
        self,
        filters: list[PropertyFilterBase],
    ):
        """Constructor for the PropertyFilter class."""
        self.filters = filters

    def apply(self, sequence_type: _TYPES_SEQ, oligo_database: OligoDatabase, n_jobs: int = 1):
        """Applies all filters to the oligo database in parallel, modifying it in place.

        :param sequence_type: The type of sequences being filtered, must be one of the predefined sequence types.
        :type sequence_type: _TYPES_SEQ
        :param oligo_database: The oligo database to which the filters are applied.
        :type oligo_database: OligoDatabase
        :param n_jobs: The number of parallel jobs to run. Default is 1.
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
                delayed(self._filter_region)(sequence_type, region_id, oligo_database)
                for region_id in region_ids
            )

        oligo_database.remove_regions_with_insufficient_oligos(pipeline_step="Property Filters")
        return oligo_database

    def _filter_region(self, sequence_type: _TYPES_SEQ, region_id: str, oligo_database: OligoDatabase):
        """Applies filters to a specific region of the database.

        :param sequence_type: The type of sequences being filtered, must be one of the predefined sequence types.
        :type sequence_type: _TYPES_SEQ
        :param database_region: A region from the oligo database.
        :type database_region: dict
        :return: The filtered region.
        :rtype: dict
        """
        oligo_ids = list(oligo_database.database[region_id].keys())
        for oligo_id in oligo_ids:
            fulfills_all_filter = self._filter_sequence(
                oligo_database.database[region_id][oligo_id][sequence_type]
            )
            if not fulfills_all_filter:
                del oligo_database.database[region_id][oligo_id]

    def _filter_sequence(self, sequence):
        """Applies filters to a single oligo sequence and returns the filtering outcome.

        :param sequence: The oligo sequence to be filtered.
        :type sequence: Bio.SeqUtils.Seq
        :return: Filtering result.
        :rtype: bool
        """
        for filter in self.filters:
            fulfills_filter = filter.apply(sequence)
            if not fulfills_filter:  # stop at the first false we obtain
                return False
        return True
