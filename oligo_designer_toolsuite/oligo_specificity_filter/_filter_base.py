import logging
import os
import sys
from abc import ABC, abstractmethod
from pathlib import Path

import joblib
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from ..IO._data_parser import write_oligos_DB_tsv


class ProbeFilterBase(ABC):
    def __init__(
        self,
        n_jobs,
        dir_output,
        dir_annotations=None,
        number_subbatches=None,
        max_genes_in_batch=300,
    ):
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
        self.dir_output = dir_output
        self.max_genes_in_batch = max_genes_in_batch  # if more than 300 genes in one batch split into subbatches to reduce required memory for loading blast results
        self.duplicated_sequences = None
        self.number_subbatches = number_subbatches

        if n_jobs == None:
            self.n_jobs = joblib.cpu_count()
        else:
            self.n_jobs = n_jobs

        # set logger
        log_file = os.path.join(self.dir_output, "logs/filter_probes.log")
        file_handler = logging.FileHandler(filename=log_file)
        stdout_handler = logging.StreamHandler(stream=sys.stdout)
        handlers = [file_handler, stdout_handler]

        logging.basicConfig(level=logging.INFO, handlers=handlers)
        self.logging = logging.getLogger("filter_probes")

        # set directory
        if dir_annotations == None:
            self.dir_annotations = os.path.join(self.dir_output, "annotations")
            Path(self.dir_annotations).mkdir(parents=True, exist_ok=True)
        else:
            self.dir_annotations = dir_annotations

    @abstractmethod
    def apply(self, probe_info):
        """Apply filter to list of all possible probes in probe_info dictionary given user-specified genes and save output in dir_output

        :param probe_info: probe info of user-specified genes
        :type probe_info: dict
        """

    def create_batches(self, probe_info):
        """Create batches of files that subdivide the probe info such that all probes of a gene are in one batch. These batches are then used as input to the filters. Number of batches is equal to n_jobs.

        :param probe_info: probe info of user-specified genes
        :type probe_info: dict
        """

        write_oligos_DB_tsv(probe_info, self.dir_annotations + "/tmp.tsv")
        probeinfo_tsv = pd.read_csv(
            self.dir_annotations + "/tmp.tsv", sep="\t", header=0
        )

        os.remove(self.dir_annotations + "/tmp.tsv")

        self.genes = list(probe_info.keys())

        self.removed_genes = self.genes

        if self.number_subbatches == None:
            self.number_subbatches = (
                len(self.genes) // self.max_genes_in_batch
            ) + 1  # if number of genes in batch > max_genes_in_batch then split batch into multiple subbatches

        batch_size = int(len(self.genes) / self.n_jobs) + (
            len(self.genes) % self.n_jobs > 0
        )

        for batch_id in range(self.n_jobs):

            file_probe_info_batch = os.path.join(
                self.dir_annotations, "probes_info_batch{}.txt".format(batch_id)
            )

            # batch probe info
            genes_batch = self.genes[
                (batch_size * batch_id) : (
                    min(batch_size * (batch_id + 1), len(self.genes) + 1)
                )
            ]

            probeinfo_batch = probeinfo_tsv.loc[
                probeinfo_tsv["gene_id"].isin(genes_batch)
            ].copy()
            probeinfo_batch.reset_index(inplace=True, drop=True)
            probeinfo_batch.to_csv(
                file_probe_info_batch, sep="\t", header=True, index=False
            )

            # sequences = probeinfo[probeinfo.columns[1]]
            # with open(file_probe_sequence_batch, "w") as handle_out:
            #     for seq in sequences:
            #         handle_out.write(seq + "\n")

            for subbatch_id in range(self.number_subbatches):

                file_probe_sequence_subbatch = os.path.join(
                    self.dir_annotations,
                    "probes_sequence_batch{}_{}.fna".format(batch_id, subbatch_id),
                )

                genes_subbatch = self.genes[
                    (subbatch_id * self.max_genes_in_batch) : (
                        (subbatch_id + 1) * self.max_genes_in_batch
                    )
                ]
                probes_info_subbatch = probeinfo_batch.loc[
                    probeinfo_batch["gene_id"].isin(genes_subbatch)
                ].copy()
                probes_info_subbatch.reset_index(inplace=True, drop=True)

                output = []
                for row in probes_info_subbatch.index:
                    header = probes_info_subbatch.iloc[
                        row, probes_info_subbatch.columns.get_loc("probe_id")
                    ]

                    sequence = Seq(
                        probes_info_subbatch.iloc[
                            row,
                            probes_info_subbatch.columns.get_loc("probe_sequence"),
                        ]
                    )
                    output.append(SeqRecord(sequence, header, "", ""))

                with open(file_probe_sequence_subbatch, "w") as handle:
                    SeqIO.write(output, handle, "fasta")
