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
class SpecificityFilterBase(ABC):
    @abstractmethod
    def apply(self, database: dict, file_reference: str, n_jobs: int):
        """Apply filter to list of all possible oligos in oligo_info dictionary and filter out the oligos which don't fulfill the requirements.
        Temporary files can be written in the ``self.dir_specificiy`` folder, but they must be removed.

        :param database: database containing the oligos and their features
        :type database: dict
        :param index_name: path to the file that is used as an index for the alignment
        :type index_name: str
        :param n_jobs: number of simultaneous parallel computations
        :type n_jobs: int
        :return: oligo info of user-specified genes
        :rtype: dict
        """

    def _filter_matching_oligos(
        self, database_region: dict, matching_oligos: list[str]
    ):
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
            if oligo_id in matching_oligos:
                del database_region[oligo_id]
        return database_region


class AlignmentSpecificityFilter(SpecificityFilterBase):
    """This is the base class for all specificity filter classes

    :param dir_specificity: directory where alignement temporary files can be written
    :type dir_specificity: str
    """

    def __init__(self, dir_specificity: str):
        """Construnctor"""
        # folder where we write the intermediate files
        self.dir_specificity = dir_specificity
        Path(self.dir_specificity).mkdir(parents=True, exist_ok=True)

    @abstractmethod
    def _create_index(self, file_reference: str, n_jobs: int):
        """
        Abstract method for creating an index based on the provided file reference.

        :param file_reference: path to the file that will be used as reference for the alignement.
        :type file_reference: str
        :param n_jobs: number of jobs for parallel processing during index creation.
        :type n_jobs: int
        :returns: name of the created or initialized database.
        :rtype: str
        :raises NotImplementedError: If not overridden in a subclass.
        """

    @abstractmethod
    def _run_search(self, database, region, index_name, filter_same_region_matches):
        """
        Abstract method for running a search of database region and the against the provided index, and optionally filters matches from the same region.

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

    def _run_filter(
        self, database, region, database_name, filter_same_region_matches=True, **kwargs
    ):
        matching_oligos, _ = self._run_search(
            database, region, database_name, filter_same_region_matches, **kwargs
        )
        filtered_database_region = self._filter_matching_oligos(
            database[region], matching_oligos
        )
        return filtered_database_region

    def get_matching_oligo_pairs(self, database: dict, reference_fasta: str, **kwargs):
        """_summary_

        :param database: _description_
        :type database: dict
        :param reference_fasta: _description_
        :type reference_fasta: str
        :return: _description_
        :rtype: _type_
        """
        database_name = self._create_index(reference_fasta, n_jobs=1)
        matches = self._run_search(
            database,
            region=None,
            index_name=database_name,
            filter_same_region_matches=False,
            **kwargs,
        )
        matches = matches[1]
        return list(zip(matches["query"].values, matches["reference"].values))

    # TODO: Both these functions are temporary, should be solved with database.write_to_fasta
    def _create_fasta_file(self, database, directory, region):
        """
        Temporary function
        Creates a fasta file with all the oligos of a specific gene. The fasta files are then used by the alignement tools.

        :param database: database with all the oligos contained
        :type database: dict
        :param dir: path to the directory where the files will be stored
        :type dir: str
        :param region: id of the region we are processing
        :type region: str
        :return: path of the fasta file generated
        :rtype: str
        """
        file_fasta_region = os.path.join(directory, f"oligos_{region}.fna")
        output = []
        for oligo_id in database[region].keys():
            output.append(
                SeqRecord(database[region][oligo_id]["sequence"], oligo_id, "", "")
            )
        with open(file_fasta_region, "w") as handle:
            SeqIO.write(output, handle, "fasta")
        return file_fasta_region

    def _create_fasta_multiple_regions(self, database, directory, regions):
        """Temporary function"""
        file_fasta = os.path.join(directory, "oligos_database.fna")
        output = []
        for region in regions:
            for oligo_id in database[region].keys():
                output.append(
                    SeqRecord(database[region][oligo_id]["sequence"], oligo_id, "", "")
                )
        with open(file_fasta, "w") as handle:
            SeqIO.write(output, handle, "fasta")
        return file_fasta
