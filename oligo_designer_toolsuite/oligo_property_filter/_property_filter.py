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

    def __init__(
        self,
        filters: list[PropertyFilterBase],
    ) -> None:
        """Constructor for the PropertyFilter class."""
        self.filters = filters

    def apply(
        self, oligo_database: OligoDatabase, sequence_type: _TYPES_SEQ, n_jobs: int = 1
    ) -> OligoDatabase:

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

<<<<<<< HEAD
        :param sequence_type: The type of sequences being filtered, must be one of the predefined sequence types.
        :type sequence_type: _TYPES_SEQ
        :param database_region: A region from the oligo database.
        :type database_region: dict
        :return: None
        :rtype: None
        """
=======
>>>>>>> origin/pipelines
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
