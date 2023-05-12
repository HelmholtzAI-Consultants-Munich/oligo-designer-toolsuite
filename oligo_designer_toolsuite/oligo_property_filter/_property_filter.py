############################################
# imports
############################################

from joblib import Parallel, delayed

from . import PropertyFilterBase
from ..database import OligoDatabase

############################################
# Property Filter Class
############################################


class PropertyFilter:
    """
    Applies sequentially all the specified pre-filters based on the sequences features. All the oligos not fulfilling all the constraints are deleted form the database.

    :param filters: list of filters classes already initialized
    :type filters: list of classes
    """

    def __init__(
        self,
        filters: list[PropertyFilterBase],
    ) -> None:
        """
        Constructor.
        """
        self.filters = filters

    def apply(self, oligo_database: OligoDatabase, n_jobs: int = None):
        """Filters the database of oligos based on the given filters.

        :param oligo_database: database class containig the oligos and their features
        :type oligo_database: OligoDatabase
        :param n_jobs: nr of cores used, if None the value set in database class is used, defaults to None
        :type n_jobs: int
        :return: oligo database class cointainig only oligos that passed the filters
        :rtype: OligoDatabase
        """
        if n_jobs is None:
            n_jobs = oligo_database.n_jobs

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
            pipeline_step="Property Filter"
        )
        return oligo_database

    def _filter_region(self, database_region):
        """Apply filters to alll oligos of a specific region.

        :param database_region: Oligo database entry for one region.
        :type database_region: dict
        :return: Oligo database entry for one region only contaning oligos that passed the filter.
        :rtype: dict
        """
        oligos_id = list(database_region.keys())
        for oligo_id in oligos_id:
            fulfills, additional_features = self._filter_sequence(
                database_region[oligo_id]["sequence"]
            )
            if fulfills:
                database_region[oligo_id].update(additional_features)
            else:
                del database_region[oligo_id]
        return database_region

    def _filter_sequence(self, sequence):
        """Applies the user-defined filters and returns the result and the additional computed features.

        :param sequence: sequence to check
        :type sequence: Bio.Seq
        :return: if the filters are fulfilled and the additional features computed
        :rtype: bool, dict
        """
        fulfills = True
        additional_features = {}
        for filter in self.filters:
            out, feat = filter.apply(sequence)
            if not out:  # stop at the first false we obtain
                return False, {}
            additional_features.update(feat)
        return fulfills, additional_features
