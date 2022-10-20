import os

from joblib import Parallel, delayed


class PreFilter:
    """
    Applies sequentially all the specified pre-filters based on the sequences features. All the probes not fulfilling all the constraints are deleted form the database.

    :param filters: list of filters classes already initialized
    :type filters: list of classes
    """

    def __init__(self, filters) -> None:
        """
        Constructor.
        """

        self.filters = filters

    def apply(self, DB, n_jobs=None):
        """Filters the database of probes based on the given filters

        :param DB: database class containig the probes and their features
        :type DB: CustomDB class
        :param n_jobs: nr of cores used, if None the value set in DB class is used, defaults to None
        :type n_jobs: int
        :return: database classs cointainig filtered oligos
        :rtype: CustomDB class
        """
        # TODO make it parallel and take into account that some genes might disappear
        if n_jobs is None:
            n_jobs = DB.n_jobs
        file_removed_genes = os.path.join(
            DB.dir_output, "genes_with_insufficient_probes.txt"
        )
        oligos_DB = DB.oligos_DB
        gene_ids = list(oligos_DB.keys())
        filtered_probes = Parallel(n_jobs=n_jobs)(
            delayed(self._filter_gene)(oligos_DB[gene]) for gene in gene_ids
        )
        oligos_DB = {}
        for probes_gene, gene in zip(filtered_probes, gene_ids):
            if probes_gene == {}:
                # all teh probes were filtered out
                with open(DB.file_removed_genes, "a") as handle:
                    handle.write(f"{gene}\tPre filter\n")
            else:
                oligos_DB[gene] = probes_gene

        DB.oligos_DB = oligos_DB
        return DB

    def _filter_gene(self, probes_gene):
        probes_id = list(probes_gene.keys())
        for probe_id in probes_id:
            fulfills, additional_features = self._filter_sequence(
                probes_gene[probe_id]["probe_sequence"]
            )
            if fulfills:
                probes_gene[probe_id].update(additional_features)
            else:
                del probes_gene[probe_id]
        return probes_gene

    def _filter_sequence(self, sequence):
        """Applies the used-defined filters and returns the result and the additional computed features

        :param sequence: sequence to check
        :type sequence: Bio.Seq
        :return: if the filters are fulfilled and the additional features computed
        :rtype: bool, dict
        """

        fulfills = True
        additional_features = {}
        for filter in self.filters:
            out, feat = filter.apply(sequence)
            if not out:  # stop at the first false we obtain
                return False, {}
            additional_features.update(feat)
        return fulfills, additional_features
