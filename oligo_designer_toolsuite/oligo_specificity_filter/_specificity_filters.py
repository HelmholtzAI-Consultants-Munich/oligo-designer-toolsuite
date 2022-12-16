from ..IO import CustomDB
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

    def apply(self, database: CustomDB, n_jobs: int = None):
        """
        Applies all the class filter sequentially to the ``oligos_DB`` stored in the data_base class given in input. Each filter parallelized and ``n_jobs`` represents the maximum number of cores available.

        :param data_base: Database class to filter.
        :type data_base: CustomDB class
        :param n_jobs: maximum number of available jobs, if ``None`` the number of jobs of the ``data_base`` class is used, defaults to None
        :type n_jobs: int, optional

        :return: filtered data_base class
        :rtype: CustomDB class
        """

        if n_jobs is None:
            n_jobs = database.n_jobs

        oligos_DB = database.oligos_DB
        for filter in self.filters:
            oligos_DB = filter.apply(oligos_DB, database.file_reference_DB, n_jobs)

        database.oligos_DB = oligos_DB
        database.remove_genes_with_insufficient_probes(
            "Specificity filter", self.write_genes_with_insufficient_probes
        )
        return database
