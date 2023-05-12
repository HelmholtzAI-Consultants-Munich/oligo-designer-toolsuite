############################################
# imports
############################################

from ..database import OligoDatabase, ReferenceDatabase
from . import SpecificityFilterBase

############################################
# Specificity Filter Classe
############################################

class SpecificityFilter:
    """
    This class applies the specificity filters given in input sequentially to a ``CustomDB`` class.

    :param filters: List of all the filter classes that we want to apply to the database
    :type filters: list of ``SpecificityFilterBase`` class
    """

    def __init__(
        self,
        filters: list[SpecificityFilterBase],
    ):
        """
        Constructor.
        """

        self.filters = filters

    def apply(
        self,
        oligo_database: OligoDatabase,
        reference_database: ReferenceDatabase,
        n_jobs: int = None,
    ):
        """
        Applies all the class filter sequentially to the ``oligos_DB`` stored in the data_base class given in input. Each filter parallelized and ``n_jobs`` represents the maximum number of cores available.

        :param oligo_database: Oligos Database class to filter.
        :type oligo_database: CustomOligoDB class
        :param oligo_database: Database class containing the reference region for the alignement methods.
        :type oligo_database: CustomReferenceDB class
        :param n_jobs: maximum number of available jobs, if ``None`` the number of jobs of the ``data_base`` class is used, defaults to None
        :type n_jobs: int, optional

        :return: filtered data_base class
        :rtype: CustomDB class
        """

        if n_jobs is None:
            n_jobs = oligo_database.n_jobs

        database = oligo_database.database
        for filter in self.filters:
            database = filter.apply(database, reference_database.file_fasta, n_jobs)

        oligo_database.database = database
        oligo_database.remove_regions_with_insufficient_oligos("Specificity filter")
        return oligo_database
