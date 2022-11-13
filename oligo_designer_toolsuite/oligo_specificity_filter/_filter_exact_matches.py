import iteration_utilities
from joblib import Parallel, delayed

from . import SpecificityFilterBase


class ExactMatches(SpecificityFilterBase):
    def __init__(self, dir_specificity):
        """This class filters probes based on exact matches."""

        super().__init__(dir_specificity)

    def apply(self, oligo_DB, file_reference_DB, n_jobs):
        """Parallelize filtering of exact matches from input dictionary"""

        duplicated_sequences = self._get_duplicated_sequences(oligo_DB)

        # run filter with joblib
        genes = list(oligo_DB.keys())
        filtered_oligo_DBs = Parallel(n_jobs=n_jobs)(
            delayed(self._filter_exactmatch_gene)(oligo_DB[gene], duplicated_sequences)
            for gene in genes
        )
        # reconstruct the database
        for gene, filtered_oligo_DB in zip(genes, filtered_oligo_DBs):
            oligo_DB[gene] = filtered_oligo_DB

        return oligo_DB

    def _get_duplicated_sequences(self, oligo_DB):
        """Get a list of probe sequences that have an exact match within the pool of all
        possible probe sequences for the list of input genes.

        :return: List of probe sequences with exact matches in the pool of probes.
        :rtype: list
        """
        # extract all the sequences
        sequences = []
        for gene in oligo_DB.keys():
            for probe_id in oligo_DB[gene].keys():
                sequences.append(
                    oligo_DB[gene][probe_id]["probe_sequence"].upper()
                )  # sequences might be also written in lower letters
        # find the duplicates within the database
        duplicated_sequences = list(
            iteration_utilities.unique_everseen(
                iteration_utilities.duplicates(sequences)
            )
        )

        return duplicated_sequences

    def _filter_exactmatch_gene(self, gene_oligo_DB, duplicated_sequences):
        """Remove sequences with exact matches within the pool of all possible probe sequences for the list of input genes.

        :param batch_id: Batch ID.
        :type batch_id: int
        """
        propbes_ids = list(gene_oligo_DB.keys())
        for probe_id in propbes_ids:
            if (
                gene_oligo_DB[probe_id]["probe_sequence"].upper()
                in duplicated_sequences
            ):
                del gene_oligo_DB[probe_id]

        return gene_oligo_DB
