import logging
import os
from abc import ABC, abstractmethod
from pathlib import Path


class ProbeFilterBase(ABC):
    """This is the base class for all filter classes

    :param no_batches: How many sub-processes to use
    :type no_batches: int
    :param ligation_region: If ligation region present
    :type ligation_region: bool"""

    def __init__(
        self,
        number_batches,
        ligation_region,
        dir_output,
        file_transcriptome_fasta,
        genes,
        min_probes_per_gene,
        dir_annotations=None,
        number_subbatches=1,
        max_genes_in_batch=300,
    ):
        self.number_batches = number_batches
        self.ligation_region = ligation_region
        self.dir_output = dir_output
        self.number_subbatches = number_subbatches  # if number of genes in batch > max_genes_in_batch then split batch into multiple subbatches
        self.max_genes_in_batch = max_genes_in_batch  # if more than 300 genes in one batch split into subbatches to reduce required memory for loading blast results
        self.min_probes_per_gene = min_probes_per_gene
        self.file_transcriptome_fasta = file_transcriptome_fasta
        self.removed_genes = genes
        # set logger
        self.logging = logging.getLogger("filter_probes")

        # set directory
        if dir_annotations == None:
            self.dir_annotations = os.path.join(self.dir_output, "annotations")
            Path(self.dir_annotations).mkdir(parents=True, exist_ok=True)
        else:
            self.dir_annotations = dir_annotations

    @abstractmethod
    def apply(self):
        """Apply filter to file_transcriptome_fasta and save output in dir_output"""
