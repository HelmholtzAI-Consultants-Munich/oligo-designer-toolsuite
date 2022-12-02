import os
from abc import ABC, abstractmethod

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


class SpecificityFilterBase(ABC):
    """This is the base class for all specificity filter classes

    :param dir_specificity: directory where alignement temporary files can be written
    :type dir_specificity: str
    """

    def __init__(self, dir_specificity):
        """Construnctor"""

        self.dir_specificity = (
            dir_specificity  # folder where we write the intermediate files
        )
        os.makedirs(self.dir_specificity, exist_ok=True)

    @abstractmethod
    def apply(self, oligo_DB, file_reference_DB, n_jobs):
        """Apply filter to list of all possible probes in probe_info dictionary and filter out the probes which don't fulfill the requirements.
        Temporary files can be written in the ``self.dir_specificiy``folder, but they must be removed.

        :param oligo_DB: database containing the probes and their features
        :type oligo_DB: dict
        :param file_reference_DB: path to the file that will be used as reference for the alignement tools
        :type file_reference_DB: str
        :param n_jobs: number of simultaneous parallel computations
        :type n_jobs: int
        :return: probe info of user-specified genes
        :rtype : dict
        """

    def _create_fasta_file(self, gene_DB, dir, gene):
        """Creates a fasta file with all the probes of a specific gene. The fasta files are then used by the alignement tools.

        :param gene_DB: database with all the probes contained in one gene
        :type gene_DB: dict
        :param dir: path to the directory where the files will be stored
        :type dir: str
        :param gene: id of the gene we are processing
        :type gene: str
        :return: path of the fasta file generated
        :rtype: str
        """
        file_fasta_gene = os.path.join(dir, f"probes_{gene}.fna")
        output = []
        for probe_id in gene_DB.keys():
            output.append(
                SeqRecord(gene_DB[probe_id]["probe_sequence"], probe_id, "", "")
            )
        with open(file_fasta_gene, "w") as handle:
            SeqIO.write(output, handle, "fasta")
        return file_fasta_gene

    def _filter_matching_probes(self, gene_DB, matching_probes):
        """Filer out form the database the sequences with a match.

        :param gene_DB: dictionary with all the probes belonging to the current gene
        :type gene_DB: dict
        :param matching_probes: list of the probes with a match
        :type matching_probes: list
        :return: gene_DB without the matching probes
        :rtype: dict
        """
        probe_ids = list(gene_DB.keys())
        for probe_id in probe_ids:
            if probe_id in matching_probes:
                del gene_DB[probe_id]
        return gene_DB
