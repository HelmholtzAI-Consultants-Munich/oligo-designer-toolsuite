import os
from abc import ABC, abstractmethod

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


class SpecificityFilterBase(ABC):
    """This is the base class for all probe filter classes

    :param n_jobs: numbers of jobs to run in parallel when filtering the probes
    :type n_jobs: int
    :param dir_output: output directory to write file containing filtered probes to
    :type dir_output: str
    :param dir_annotations: directory storing file annotations of probes, defaults to None
    :type dir_annotations: str, optional
    :param number_subbatches: number of subbatches to split batched data into to control how many genes are in one batch, defaults to None
    :type number_subbatches: int, optional
    :param max_genes_in_batch: max number of genes allowed in a batch, defaults to 300
    :type max_genes_in_batch: int, optional
    """

    def __init__(self, dir_specificity):
        """Construnctor"""

        self.dir_specificity = (
            dir_specificity  # folder where we write the intermediate files
        )
        os.makedirs(self.dir_specificity, exist_ok=True)

    @abstractmethod
    def apply(self, oligo_DB, file_reference_DB, n_jobs):
        """Apply filter to list of all possible probes in probe_info dictionary given user-specified genes and save output in dir_output;
        Temporary files can be written in the `self.dir_specificiy` folder, but they must be removed.

        :param probe_info: probe info of user-specified genes
        :type probe_info: dict
        """

    def _create_fasta_file(self, gene_DB, dir, gene):

        file_fasta_gene = os.path.join(dir, f"probes_{gene}.fna")
        output = []
        for probe_id in gene_DB.keys():
            output.append(
                SeqRecord(gene_DB[probe_id]["probe_sequence"], probe_id, "", "")
            )
        with open(file_fasta_gene, "w") as handle:
            SeqIO.write(output, handle, "fasta")
        return file_fasta_gene

    def __del__(self):
        if os.path.exists(self.dir_specificity):
            pass  # shutil.rmtree(self.dir_specificity)
