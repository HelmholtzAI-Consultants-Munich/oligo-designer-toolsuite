############################################
# imports
############################################

from typing import get_args

from oligo_designer_toolsuite._constants import _TYPES_SEQ
from oligo_designer_toolsuite.database import OligoDatabase, ReferenceDatabase
from oligo_designer_toolsuite.oligo_specificity_filter import SpecificityFilterBase

############################################
# Specificity Filter Classe
############################################


class SpecificityFilter:
    def __init__(
        self,
        filters: list[SpecificityFilterBase],
    ) -> None:
        """Constructor for the SpecificityFilter class."""
        self.filters = filters

    def apply(
        self,
        oligo_database: OligoDatabase,
        reference_database: ReferenceDatabase = None,
        sequence_type: _TYPES_SEQ = "oligo",
        n_jobs: int = 1,
    ) -> OligoDatabase:
        options = get_args(_TYPES_SEQ)
        assert (
            sequence_type in options
        ), f"Sequence type not supported! '{sequence_type}' is not in {options}."

        for specificity_filter in self.filters:
            oligo_database = specificity_filter.apply(
                oligo_database=oligo_database,
                reference_database=reference_database,
                sequence_type=sequence_type,
                n_jobs=n_jobs,
            )
            oligo_database.remove_regions_with_insufficient_oligos("Specificity Filters")

        return oligo_database
