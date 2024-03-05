############################################
# imports
############################################

from typing import get_args
from .._constants import _TYPES_SEQ

from . import SpecificityFilterBase
from ..database import OligoDatabase, ReferenceDatabase

############################################
# Specificity Filter Classe
############################################


class SpecificityFilter:
    """A class to apply a series of specificity filters to an oligonucleotide database to ensure that oligos
    do not bind to off-targets of a given reference databse or cross-hybridize with other oligos in the oligo databse.

    :param filters: A list of filter instances derived from SpecificityFilterBase that define the specificity criteria.
    :type filters: list[SpecificityFilterBase]
    """

    def __init__(
        self,
        filters: list[SpecificityFilterBase],
    ):
        """Constructor for the SpecificityFilter class."""
        self.filters = filters

    def apply(
        self,
        sequence_type: _TYPES_SEQ,
        oligo_database: OligoDatabase,
        n_jobs: int = None,
        reference_database: ReferenceDatabase = None,
    ):
        """Applies all provided specificity filters to the oligo database against a reference database.

        :param sequence_type: The type of sequences being filtered, must be one of the predefined sequence types.
        :type sequence_type: _TYPES_SEQ
        :param oligo_database: The database of oligonucleotides to be filtered.
        :type oligo_database: OligoDatabase
        :param n_jobs: The number of parallel jobs to run. Defaults to the number of jobs defined in oligo_database.
        :type n_jobs: int, optional
        :param reference_database: The reference database to compare against for specificity.
            For non-alignment based specificity filter reference_database is not used, i.e. set to None.
        :type reference_database: ReferenceDatabase, optional
        :return: The filtered oligo database.
        :rtype: OligoDatabase
        """
        options = get_args(_TYPES_SEQ)
        assert (
            sequence_type in options
        ), f"Sequence type not supported! '{sequence_type}' is not in {options}."

        if n_jobs is None:
            n_jobs = oligo_database.n_jobs

        for filter in self.filters:
            oligo_database = filter.apply(sequence_type, oligo_database, reference_database, n_jobs)

        oligo_database.remove_regions_with_insufficient_oligos("Specificity Filters")
        return oligo_database
