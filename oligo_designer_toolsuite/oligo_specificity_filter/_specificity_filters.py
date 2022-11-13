class SpecificityFilter:
    def __init__(self, filters):
        self.filters = filters

    def apply(self, DB, n_jobs=None, write_empty_genes=True):

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
