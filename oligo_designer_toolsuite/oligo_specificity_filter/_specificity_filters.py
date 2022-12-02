class SpecificityFilter:
    """This class applies the specificity filters given in input sequentially to a ``CustomDB`` class.

    :param filters: List of all the filter classes that we want to apply to the database
    :type filters: list of ``SpecificityFilterBase`` class
    """

    def __init__(self, filters):
        self.filters = filters

    def apply(self, data_base, n_jobs=None, write_empty_genes=True):
        """Applies all the class filter sequentially to the ``oligos_DB`` stored in the data_base class given in input. Each filter parallelized and ``n_jobs`` represents the maximum number of cores available.

        :param data_base: Database class to filter.
        :type data_base: CustomDB class
        :param n_jobs: maximum number of available jobs, if ``None`` the number of jobs of the ``data_base`` class is used, defaults to None
        :type n_jobs: int, optional
        :param write_empty_genes: if True genes with insufficient probes are written in a file, defaults to True
        :type write_empty_genes: bool, optional
        :return: filtered data_base class
        :rtype: CustomDB class
        """

        if n_jobs is None:
            n_jobs = data_base.n_jobs

        oligos_DB = data_base.oligos_DB
        for filter in self.filters:
            oligos_DB = filter.apply(oligos_DB, data_base.file_reference_DB, n_jobs)

        data_base.oligos_DB = oligos_DB
        data_base.remove_genes_with_insufficient_probes(
            "Specificity filter", write_empty_genes
        )
        return data_base
