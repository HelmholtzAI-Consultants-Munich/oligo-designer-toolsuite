class SpecificityFilter:
    """Thus class applies the specificity filters given in input sequentally to a ``CustomDB`` class.

    :param filters: List of all the filter classes that we want to apply to the database
    :type filters: list of ``SpecificityFilterBase`` class
    """

    def __init__(self, filters):
        self.filters = filters

    def apply(self, DB, n_jobs=None, write_empty_genes=True):
        """Applies all the class filter sequentailly to the ``oligos_DB`` sotored in the DB class given in input. Each filter parallelized and ``n_jobs`` represents the maximum number of cores available.

        :param DB: Database class to filter.
        :type DB: CustomDB class
        :param n_jobs: maximum number of available jobs, if ``None`` hte number of jobs of the ``DB`` class is used, defaults to None
        :type n_jobs: int, optional
        :param write_empty_genes: if True genes with insufficient probes are written in a file, defaults to True
        :type write_empty_genes: bool, optional
        :return: filtered DB class
        :rtype: CustomDB class
        """

        if n_jobs is None:
            n_jobs = DB.n_jobs

        oligos_DB = DB.oligos_DB
        for filter in self.filters:
            oligos_DB = filter.apply(oligos_DB, DB.file_reference_DB, n_jobs)

        DB.oligos_DB = oligos_DB
        DB.remove_genes_with_insufficient_probes(
            "Specificity filter", write_empty_genes
        )
        return DB
