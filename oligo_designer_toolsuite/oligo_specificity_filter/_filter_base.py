import logging
import os
from abc import abstractmethod
from pathlib import Path

import joblib


class ProbeFilterBase:

    """This is the base class for all filter classes

    :param number_batches: How many batches of files to split data into for filtering
    :type number_batches: int
    :param ligation_region: If ligation region is present
    :type ligation_region: bool"""

    def __init__(
        self,
        number_batches,
        ligation_region,
        dir_output,
        file_probe_info,
        genes,
        n_jobs=None,
        dir_annotations=None,
        number_subbatches=None,
        max_genes_in_batch=300,
    ):
        """This is the base class for all filter classes

        :param number_batches: How many batches of files to split data into for filtering
        :type number_batches: int
        :param ligation_region: If ligation region is present
        :type ligation_region: bool
        :param dir_output: output directory to write file containing filtered probes to
        :type dir_output: str
        :param file_probe_info: gtf file containing probe info
        :type file_probe_info: str
        :param genes: list of genes given by user
        :type genes: list
        :param n_jobs: number of processes to run in parallel when filtering probes, defaults to number_batches
        :type n_jobs: int, optional
        :param dir_annotations: directory storing file annotations of probes, defaults to None
        :type dir_annotations: str, optional
        :param number_subbatches: number of subbatches to split batched data into to control how many genes are in one batch, defaults to None
        :type dir_annotations: int, optional
        :param max_genes_in_batch: max number of genes allowed in a batch, defaults to 300
        :type max_genes_in_batch: int, optional
        """
        self.number_batches = number_batches
        self.ligation_region = ligation_region
        self.dir_output = dir_output
        self.max_genes_in_batch = max_genes_in_batch  # if more than 300 genes in one batch split into subbatches to reduce required memory for loading blast results
        self.file_probe_info = file_probe_info
        self.genes = genes
        self.removed_genes = genes

        if number_subbatches == None:
            self.number_subbatches = (
                len(genes) // self.max_genes_in_batch
            ) + 1  # if number of genes in batch > max_genes_in_batch then split batch into multiple subbatches
        else:
            self.number_subbatches = number_subbatches

        if n_jobs == None:
            self.n_jobs = joblib.cpu_count()
        else:
            self.n_jobs = n_jobs

        # set logger
        log_file = os.path.join(self.dir_output, "logs/filter_probes.log")
        logging.basicConfig(filename=log_file, filemode="w", level=logging.INFO)
        self.logging = logging.getLogger("filter_probes")

        # set directory
        if dir_annotations == None:
            self.dir_annotations = os.path.join(self.dir_output, "annotations")
            Path(self.dir_annotations).mkdir(parents=True, exist_ok=True)
        else:
            self.dir_annotations = dir_annotations

    @abstractmethod
    def apply(self):
        """Apply filter to list of all possible probes given user-specified genes and save output in dir_output"""
