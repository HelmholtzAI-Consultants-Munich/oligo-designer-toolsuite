############################################
# imports
############################################

from oligo_designer_toolsuite.database import OligoDatabase
from oligo_designer_toolsuite.oligo_specificity_filter import SpecificityFilterBase

############################################
# Specificity Filter Classe
############################################


class SpecificityFilter:
    """
    A class for managing and applying a sequence of specificity filters to an OligoDatabase.

    The `SpecificityFilter` class aggregates multiple filters, each of which applies specific criteria to determine
    the suitability of oligonucleotides based on their sequence specificity. This allows for the sequential application
    of various filtering methods to refine the OligoDatabase according to specificity requirements.

    :param filters: A list of filters, each inheriting from the `SpecificityFilterBase`, that will be applied to the OligoDatabase.
    :type filters: list[SpecificityFilterBase]
    """

    def __init__(
        self,
        filters: list[SpecificityFilterBase],
    ) -> None:
        """Constructor for the SpecificityFilter class."""
        self.filters = filters

    def apply(
        self,
        oligo_database: OligoDatabase,
        n_jobs: int = 1,
    ) -> OligoDatabase:
        """
        Applies a sequence of specificity filters to an OligoDatabase.

        The `apply` method processes the provided OligoDatabase through each specificity filter in sequence.
        It evaluates the database against reference sequences, if provided, and ensures that only oligonucleotides
        meeting all specificity criteria are retained.

        :param oligo_database: The OligoDatabase containing the oligonucleotides and their associated attributes.
        :type oligo_database: OligoDatabase
        :param n_jobs: The number of parallel jobs to use for processing.
        :type n_jobs: int
        :return: The filtered OligoDatabase.
        :rtype: OligoDatabase
        """
        for specificity_filter in self.filters:
            oligo_database = specificity_filter.apply(
                oligo_database=oligo_database,
                n_jobs=n_jobs,
            )
            oligo_database.remove_regions_with_insufficient_oligos("Specificity Filters")

        return oligo_database
