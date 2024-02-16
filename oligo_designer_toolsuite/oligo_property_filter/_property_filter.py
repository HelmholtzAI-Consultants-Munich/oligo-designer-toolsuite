############################################
# imports
############################################

from joblib import Parallel, delayed

from ..database import OligoDatabase
from . import PropertyFilterBase

############################################
# Property Filter Class
############################################


class PropertyFilter:
    """
    A class representing a collection of filters to be applied on oligo databases.

    This class manages the application of multiple PropertyFilterBase instances on an oligo database,
    allowing for complex filtering strategies. Oligos which don't fulfill the filtering criteria are
    removed from the database.

    :param filters: A list of filter objects derived from PropertyFilterBase.
    :type filters: list[PropertyFilterBase]
    """

    def __init__(
        self,
        filters: list[PropertyFilterBase],
    ) -> None:
        """Constructor for the PropertyFilter class."""
        self.filters = filters

    def apply(self, oligo_database: OligoDatabase, n_jobs: int = 1):
        """
        Applies all filters to the oligo database in parallel, modifying it in place.

        :param oligo_database: The oligo database to which the filters are applied.
        :type oligo_database: OligoDatabase
        :param n_jobs: The number of parallel jobs to run. Default is 1.
        :type n_jobs: int
        :return: The filtered oligo database.
        :rtype: OligoDatabase
        """
        database = oligo_database.database
        region_ids = list(database.keys())
        database_regions = Parallel(n_jobs=n_jobs)(
            delayed(self._filter_region)(database[region]) for region in region_ids
        )
        database = {}
        for database_region, region_id in zip(database_regions, region_ids):
            database[region_id] = database_region
        oligo_database.database = database
        oligo_database.remove_regions_with_insufficient_oligos(
            pipeline_step="Property Filters"
        )
        return oligo_database

    def _filter_region(self, database_region):
        """
        Applies filters to a specific region of the database.

        :param database_region: A region from the oligo database.
        :type database_region: dict
        :return: The filtered region.
        :rtype: dict
        """
        oligo_ids = list(database_region.keys())
        for oligo_id in oligo_ids:
            fulfills, additional_features = self._filter_sequence(
                database_region[oligo_id]["sequence"]
            )
            if fulfills:
                database_region[oligo_id].update(additional_features)
            else:
                del database_region[oligo_id]
        return database_region

    def _filter_sequence(self, sequence):
        """
        Applies filters to a single oligo sequence and returns the filtering outcome.

        :param sequence: The oligo sequence to be filtered.
        :type sequence: Bio.SeqUtils.Seq
        :return: Tuple of filtering result and additional features.
        :rtype: (bool, dict)
        """
        fulfills = True
        additional_features = {}
        for filter in self.filters:
            out, feat = filter.apply(sequence)
            if not out:  # stop at the first false we obtain
                return False, {}
            additional_features.update(feat)
        return fulfills, additional_features
