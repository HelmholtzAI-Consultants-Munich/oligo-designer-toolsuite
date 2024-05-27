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

        database = oligo_database.database
        region_ids = list(database.keys())
        with joblib_progress(description="Property Filter", total = len(region_ids)):
            database_regions = Parallel(n_jobs=n_jobs)(
                delayed(self._filter_region)(sequence_type, database[region]) for region in region_ids
            )

        database = LRUDict(
            max_in_memory=oligo_database.lru_db_max_in_memory,
            storage_path=oligo_database._dir_cache_files,
        )
        for database_region, region_id in zip(database_regions, region_ids):
            database[region_id] = database_region

        oligo_database.database = database
        oligo_database.remove_regions_with_insufficient_oligos(pipeline_step="Property Filters")
        return oligo_database

    def _filter_region(self, sequence_type: _TYPES_SEQ, database_region: dict):
        """Applies filters to a specific region of the database.

        :param sequence_type: The type of sequences being filtered, must be one of the predefined sequence types.
        :type sequence_type: _TYPES_SEQ
        :param database_region: A region from the oligo database.
        :type database_region: dict
        :return: The filtered region.
        :rtype: dict
        """
        oligo_ids = list(database_region.keys())
        for oligo_id in oligo_ids:
            fulfills_all_filter = self._filter_sequence(database_region[oligo_id][sequence_type])
            if not fulfills_all_filter:
                del database_region[oligo_id]
        return database_region

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
