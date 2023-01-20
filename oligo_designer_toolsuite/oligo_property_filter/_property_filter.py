from joblib import Parallel, delayed

from ..database import CustomOligoDB
from . import PropertyFilterBase


class PropertyFilter:
    """
    Applies sequentially all the specified pre-filters based on the sequences features. All the oligos not fulfilling all the constraints are deleted form the database.

    :param filters: list of filters classes already initialized
    :type filters: list of classes
    :param write_genes_with_insufficient_oligos: if True genes with insufficient oligos are written in a file, defaults to True
    :type write_genes_with_insufficient_oligos: bool, optional
    """

    def __init__(
        self,
        filters: list[PropertyFilterBase],
        write_genes_with_insufficient_oligos: bool = True,
    ) -> None:
        """
        Constructor.
        """

        self.filters = filters
        self.write_genes_with_insufficient_oligos = write_genes_with_insufficient_oligos

    def apply(self, oligo_database: CustomOligoDB, n_jobs: int = None):
        """Filters the database of oligos based on the given filters

        :param database: database class containig the oligos and their features
        :type database: CustomDB class
        :param n_jobs: nr of cores used, if None the value set in database class is used, defaults to None
        :type n_jobs: int
        :return: database classs cointainig filtered oligos
        :rtype: CustomDB class
        """
        # TODO make it parallel and take into account that some genes might disappear
        if n_jobs is None:
            n_jobs = oligo_database.n_jobs

        oligos_DB = oligo_database.oligos_DB
        gene_ids = list(oligos_DB.keys())
        filtered_oligos = Parallel(n_jobs=n_jobs)(
            delayed(self._filter_gene)(oligos_DB[gene]) for gene in gene_ids
        )
        oligos_DB = {}
        for oligos_gene, gene in zip(filtered_oligos, gene_ids):
            oligos_DB[gene] = oligos_gene

        oligo_database.remove_genes_with_insufficient_oligos(
            pipeline_step="property filter",
            write=self.write_genes_with_insufficient_oligos,
        )

        oligo_database.oligos_DB = oligos_DB
        return oligo_database

    def _filter_gene(self, oligos_gene):
        oligos_id = list(oligos_gene.keys())
        for oligo_id in oligos_id:
            fulfills, additional_features = self._filter_sequence(
                oligos_gene[oligo_id]["sequence"]
            )
            if fulfills:
                oligos_gene[oligo_id].update(additional_features)
            else:
                del oligos_gene[oligo_id]
        return oligos_gene

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
