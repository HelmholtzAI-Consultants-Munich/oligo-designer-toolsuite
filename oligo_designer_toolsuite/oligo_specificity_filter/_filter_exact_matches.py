############################################
# imports
############################################

import iteration_utilities
from joblib import Parallel, delayed

from . import SpecificityFilterBase

############################################
# Oligo Exact Match Filter Classes
############################################


# TODO provide get_oligo_pair_hits function to integrate into crosshybridization filter
class ExactMatches(SpecificityFilterBase):
    """This class filters oligos based duplicates found in the ``database``. That is, oligos with the same sequences but belonging to different regions are filtered out.

    :param dir_specificity: directory where alignement temporary files can be written
    :type dir_specificity: str
    """

    def __init__(self):
        """Constructor for the ExactMatches class."""

    # TODO: rewrite function such that we save oligos as "reference" file into fasta file amd then compare all oligos to themselves but remove hits from same region
    def apply(self, database: dict, file_reference: str, n_jobs: int):
        """Apply the filter in parallel on the given ``database``. Each jobs filters a single region, and at the same time are generated at most ``n_job`` jobs.
        The filtered database is returned.

        :param database: database containing the oligos and their features
        :type database: dict
        :param file_reference: path to the file that will be used as reference for the alignement tools
        :type file_reference: str
        :param n_jobs: number of simultaneous parallel computations
        :type n_jobs: int
        :return: oligo info of user-specified regions
        :rtype: dict
        """

        duplicated_sequences = self._get_duplicated_sequences(database)

        # run filter with joblib
        regions = list(database.keys())
        filtered_oligo_DBs = Parallel(n_jobs=n_jobs)(
            delayed(self._filter_exactmatch_gene)(database[region], duplicated_sequences)
            for region in regions
        )
        # reconstruct the database
        for region, filtered_oligo_DB in zip(regions, filtered_oligo_DBs):
            database[region] = filtered_oligo_DB

        return database

    def _get_duplicated_sequences(self, database):
        """Get a list of oligo sequences that have an exact match within the oligos_DB.

        :param database: database with all the oligos and their features
        :type database: dict
        :return: List of oligo sequences with exact matches in the pool of oligos.
        :rtype: list
        """
        # extract all the sequences
        sequences = []
        for region in database.keys():
            for oligo_id in database[region].keys():
                sequences.append(
                    database[region][oligo_id]["sequence"].upper()
                )  # sequences might be also written in lower letters
        # find the duplicates within the database
        duplicated_sequences = list(
            iteration_utilities.unique_everseen(iteration_utilities.duplicates(sequences))
        )

        return duplicated_sequences

    def _filter_exactmatch_gene(self, database_region, duplicated_sequences):
        """Remove sequences with exact matches.

        :param database_region: database with all the oligos from one region
        :type database_region: dict
        :param duplicated_sequences: list of the sequences which have duplicates
        :type duplicated_sequences: list
        """
        propbes_ids = list(database_region.keys())
        for oligo_id in propbes_ids:
            if database_region[oligo_id]["sequence"].upper() in duplicated_sequences:
                del database_region[oligo_id]

        return database_region
