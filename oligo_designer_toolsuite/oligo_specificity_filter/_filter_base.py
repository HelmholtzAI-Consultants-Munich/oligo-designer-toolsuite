############################################
# imports
############################################

import os
from abc import ABC, abstractmethod
from pathlib import Path

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


############################################
# Oligo Specificity Filter Classes
############################################


# TODO: specify sequence type for database input (i.e. target or oligo)
class SpecificityFilterBase(ABC):
    "This is the base class for all specificity filter classes"

    def __init__(self, dir_specificity: str):
        """Construnctor"""
        # folder where we write the intermediate files
        self.dir_specificity = dir_specificity
        Path(self.dir_specificity).mkdir(parents=True, exist_ok=True)

    @abstractmethod
    def apply(self, database: dict, file_reference: str, n_jobs: int):
        """Apply filter to list of all possible oligos in oligo_info dictionary and filter out the oligos which don't fulfill the requirements.
        Temporary files can be written in the ``self.dir_specificiy`` folder, but they must be removed.

        :param database: database containing the oligos and their features
        :type database: dict
        :param file_reference: path to the file that is used as an reference for the alignment
        :type file_reference: str
        :param n_jobs: number of simultaneous parallel computations
        :type n_jobs: int
        :return: oligo info of user-specified genes
        :rtype: dict
        """

    def _filter_hits_from_database(self, database_region: dict, oligo_hits: list[str]):
        """Filer out form the database the sequences with a match.

        :param database_region: dictionary with all the oligos belonging to the current gene
        :type database_region: dict
        :param matching_oligos: list of the oligos with a match
        :type matching_oligos: list
        :return: database_region without the matching oligos
        :rtype: dict
        """
        oligo_ids = list(database_region.keys())
        for oligo_id in oligo_ids:
            if oligo_id in oligo_hits:
                del database_region[oligo_id]
        return database_region


class AlignmentSpecificityFilter(SpecificityFilterBase):
    """
    This is the base class for all specificity filter classes which are based on an alignment tool
    :param dir_specificity: directory where alignment temporary files can be written
    :type dir_specificity: str
    """

    def __init__(self, dir_specificity: str):
        """Construnctor"""
        # folder where we write the intermediate files
        self.dir_specificity = dir_specificity
        Path(self.dir_specificity).mkdir(parents=True, exist_ok=True)

    @abstractmethod
    def get_oligo_pair_hits(self, database: dict, file_reference: str, n_jobs: int):
        """
        Retrieve matching oligo pairs between a reference FASTA and a database. It returns a list of pairs, where each pair
        contains the name of the oligo from the database and its corresponding match from the reference.

        :param database: database containing the oligos.
        :type database: dict
        :param reference_fasta: path to the file that is used as an reference for the alignment
        :type reference_fasta: str

        :return: A list of matching oligo pairs.
        :rtype: list of tuple
        """

    @abstractmethod
    def _create_index(self, file_reference: str, n_jobs: int):
        """
        Abstract method for creating an index based on the provided file reference.
        The index serves as a database specific to the alignment tool, used for facilitating the alignment processes.

        :param file_reference: path to the file that will be used as reference for the alignment.
        :type file_reference: str
        :param n_jobs: number of jobs for parallel processing during index creation.
        :type n_jobs: int
        :returns: name of the created or initialized database.
        :rtype: str
        :raises NotImplementedError: If not overridden in a subclass.
        """

    @abstractmethod
    def _run_search(self, database: dict, region: str, file_index: str, **kwargs):
        """
        Abstract method for running a search of database region and the against the provided index.

        :param database: database containing the oligos
        :type database: dict
        :param region: id of the region processed
        :type region: str
        :param index_name: path of the database or index to search against
        :type index_name: str
        :param filter_same_region_matches: flag to indicate whether to filter matches from the same gene
        :type filter_same_region_matches: bool
        :returns: A tuple containing: an array of oligos with matches, and a dataframe containing alignment match data
        :rtype: (numpy.ndarray, pandas.DataFrame)
        """

    @abstractmethod
    def _find_hits(self, search_results, filter_hits_from_input_region):
        """ """

    def _run_filter(
        self, database: dict, region: str, file_index: str, filter_hits_from_input_region=True, **kwargs
    ):
        """
        Filters a database region based on alignment matches.

        :param database: database containing the oligos.
        :type database: dict
        :param region: id of the region processed
        :type region: str
        :param index_name: path of the database or index to search against
        :type index_name: str
        :param filter_same_region_matches: Whether to filter matches within the same region (default is True).
        :type filter_same_region_matches: bool
        :param kwargs: Additional keyword arguments to customize the search.

        :return: The filtered database region containing matching oligos.
        :rtype: dict
        """
        search_results = self._run_search(database, region, file_index, **kwargs)
        oligo_hits = self._find_hits(search_results, filter_hits_from_input_region)
        database_region_filtered = self._filter_hits_from_database(database[region], oligo_hits)

        return database_region_filtered

    # TODO: Both these functions are temporary, should be solved with database.write_to_fasta
    def _create_fasta_file(self, database, directory, region):
        file_fasta_region = os.path.join(directory, f"oligos_{region}.fna")
        output = []
        for oligo_id in database[region].keys():
            output.append(SeqRecord(database[region][oligo_id]["sequence"], oligo_id, "", ""))
        with open(file_fasta_region, "w") as handle:
            SeqIO.write(output, handle, "fasta")
        return file_fasta_region

    def _create_fasta_multiple_regions(self, database, directory, regions):
        """Temporary function"""
        file_fasta = os.path.join(directory, "oligos_database.fna")
        output = []
        for region in regions:
            for oligo_id in database[region].keys():
                output.append(SeqRecord(database[region][oligo_id]["sequence"], oligo_id, "", ""))
        with open(file_fasta, "w") as handle:
            SeqIO.write(output, handle, "fasta")
        return file_fasta
