import iteration_utilities
from joblib import Parallel, delayed

from . import SpecificityFilterBase


class ExactMatches(SpecificityFilterBase):
    """This class filters probes based duplicates found in the ``oligos_DB``. That is, probes with the same sequences but belonging to different genes are filtered out.

    :param dir_specificity: directory where alignement temporary files can be written
    :type dir_specificity: str
    """

    def __init__(self, dir_specificity: str):
        """Counstructor."""

        super().__init__(dir_specificity)

    def apply(self, oligo_DB: dict, file_reference_DB: str, n_jobs: int):
        """Apply the filter in parallel on the given ``oligo_DB``. Each jobs filters a single gene, and at the same time are generated at most ``n_job`` jobs.
        The filtered database is returned.

        :param oligo_DB: database containing the probes and their features
        :type oligo_DB: dict
        :param file_reference_DB: path to the file that will be used as reference for the alignement tools
        :type file_reference_DB: str
        :param n_jobs: number of simultaneous parallel computations
        :type n_jobs: int
        :return: probe info of user-specified genes
        :rtype : dict
        """

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
        """Get a list of probe sequences that have an exact match within the oligos_DB.

        :param oligo_DB: database with all the probes and their features
        :type oligo_DB: dict
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
        """Remove sequences with exact matches.

        :param gene_oligo_DB: database with all the probes from one gene
        :type gene_oligo_DB: dict
        :param duplicated_sequences: list of the sequences which have duplicates
        :type duplicated_sequences: list
        """
        propbes_ids = list(gene_oligo_DB.keys())
        for probe_id in propbes_ids:
            if (
                gene_oligo_DB[probe_id]["probe_sequence"].upper()
                in duplicated_sequences
            ):
                del gene_oligo_DB[probe_id]

        return gene_oligo_DB
