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
    def apply(self, oligo_database: dict, file_reference: str, n_jobs: int):
        """Apply filter to list of all possible oligos in oligo_info dictionary and filter out the oligos which don't fulfill the requirements.
        Temporary files can be written in the ``self.dir_specificiy`` folder, but they must be removed.

        :param oligo_database: database containing the oligos and their features
        :type oligo_database: dict
        :param file_reference: path to the file that will be used as reference for the alignement tools
        :type file_reference: str
        :param n_jobs: number of simultaneous parallel computations
        :type n_jobs: int
        :return: oligo info of user-specified genes
        :rtype: dict
        """

    def _create_fasta_file(self, database_region, dir, gene):
        """Creates a fasta file with all the oligos of a specific gene. The fasta files are then used by the alignement tools.

        :param database_region: database with all the oligos contained in one gene
        :type database_region: dict
        :param dir: path to the directory where the files will be stored
        :type dir: str
        :param gene: id of the gene we are processing
        :type gene: str
        :return: path of the fasta file generated
        :rtype: str
        """
        file_fasta_gene = os.path.join(dir, f"oligos_{gene}.fna")
        output = []
        for oligo_id in database_region.keys():
            output.append(
                SeqRecord(database_region[oligo_id]["sequence"], oligo_id, "", "")
            )
        with open(file_fasta_gene, "w") as handle:
            SeqIO.write(output, handle, "fasta")
        return file_fasta_gene

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
