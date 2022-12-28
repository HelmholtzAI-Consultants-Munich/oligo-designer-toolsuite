from ..IO import CustomOligoDB, CustomReferenceDB
from . import SpecificityFilterBase


class SpecificityFilter:
    """
    This class applies the specificity filters given in input sequentially to a ``CustomDB`` class.

    :param filters: List of all the filter classes that we want to apply to the database
    :type filters: list of ``SpecificityFilterBase`` class
    :param write_genes_with_insufficient_probes: if True genes with insufficient probes are written in a file, defaults to True
    :type write_genes_with_insufficient_probes: bool, optional
    """

    def __init__(
        self,
        filters: list[SpecificityFilterBase],
        write_genes_with_insufficient_probes: bool = True,
    ):
        """
        Constructor.
        """

        self.filters = filters
        self.write_genes_with_insufficient_probes = write_genes_with_insufficient_probes

    def apply(
        self,
        oligo_database: CustomOligoDB,
        reference_database: CustomReferenceDB,
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

        oligos_DB = oligo_database.oligos_DB
        for filter in self.filters:
            oligos_DB = filter.apply(
                oligos_DB, reference_database.file_reference_DB, n_jobs
            )

        oligo_database.oligos_DB = oligos_DB
        oligo_database.remove_genes_with_insufficient_probes(
            "Specificity filter", self.write_genes_with_insufficient_probes
        )
        return oligo_database
