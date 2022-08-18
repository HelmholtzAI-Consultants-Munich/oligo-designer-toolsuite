import logging
import os
import sys
from abc import ABC, abstractmethod
from pathlib import Path

import joblib


class ProbeFilterBase(ABC):

    """This is the base class for all filter classes

    :param number_batches: How many batches of files to split data into for filtering
    :type number_batches: int
    :param ligation_region: If ligation region is present
    :type ligation_region: bool"""

    def __init__(
        self,
        n_jobs,
        dir_output,
        genes,
        dir_annotations=None,
        number_subbatches=None,
        max_genes_in_batch=300,
        write_to_file=True,
    ):
        """This is the base class for all filter classes

        :param n_jobs: numbers of jobs to run in parallel when filtering the probes
        :type n_jobs: int
        :param dir_output: output directory to write file containing filtered probes to
        :type dir_output: str
        :param genes: list of genes given by user
        :type genes: list
        :param dir_annotations: directory storing file annotations of probes, defaults to None
        :type dir_annotations: str, optional
        :param number_subbatches: number of subbatches to split batched data into to control how many genes are in one batch, defaults to None
        :type dir_annotations: int, optional
        :param max_genes_in_batch: max number of genes allowed in a batch, defaults to 300
        :type max_genes_in_batch: int, optional
        """
        self.dir_output = dir_output
        self.max_genes_in_batch = max_genes_in_batch  # if more than 300 genes in one batch split into subbatches to reduce required memory for loading blast results
        self.genes = genes
        self.removed_genes = genes
        self.duplicated_sequences = None
        self.write_to_file = write_to_file

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
        """
        :param probe_info: probe info of user-specified genes
        :type probe_info: dict
        Apply filter to list of all possible probes in probe_info dictionary given user-specified genes and save output in dir_output"""

    '''def create_batches(self, probe_info):

        df = pd.DataFrame.from_dict(probe_info)
        probeinfo_tsv = pd.read_csv(df, sep="\t", header=0)

        batch_size = int(len(self.genes) / self.n_jobs) + (
            len(self.genes) % self.n_jobs > 0
        )

        for batch_id in range(self.n_jobs):

            file_probe_info_batch = os.path.join(
                self.dir_annotations, "probes_info_batch{}.txt".format(batch_id)
            )

            file_probe_sequence_batch = os.path.join(
                self.dir_annotations, "probes_sequence_batch{}.txt".format(batch_id)
            )

            # batch probe info
            genes_batch = self.genes[
                (batch_size * batch_id) : (
                    min(batch_size * (batch_id + 1), len(self.genes) + 1)
                )
            ]

            probeinfo = probeinfo_tsv.loc[
                probeinfo_tsv["gene_id"].isin(genes_batch)
            ].copy()
            probeinfo.reset_index(inplace=True, drop=True)
            probeinfo.to_csv(file_probe_info_batch, sep="\t", header=True, index=False)

            sequences = probeinfo[probeinfo.columns[1]]
            with open(file_probe_sequence_batch, "w") as handle_out:
                for seq in sequences:
                    handle_out.write(seq + "\n")

    def _get_duplicated_sequences(self):
        """Get a list of probe sequences that have a exact match within the pool of all
        possible probe sequences for the list of input genes.
        :return: List of probe sequences with exact matches in the pool of probes.
        :rtype: list
        """
        sequences = []

        for batch_id in range(self.n_jobs):
            file_probe_sequence_batch = os.path.join(
                self.dir_annotations, "probes_sequence_batch{}.txt".format(batch_id)
            )

            with open(file_probe_sequence_batch, "r") as handle:
                sequences_batch = [line.rstrip() for line in handle]
            sequences.extend(sequences_batch)
            os.remove(file_probe_sequence_batch)

        duplicated_sequences = list(
            iteration_utilities.unique_everseen(
                iteration_utilities.duplicates(sequences)
            )
        )

        return duplicated_sequences

    def _filter_exactmatch_batch(self, batch_id):
        """Remove sequences with exact matches within the pool of all possible probe sequences for the list of input genes.
        :param batch_id: Batch ID.
        :type batch_id: int
        """
        file_probe_info_batch = os.path.join(
            self.dir_annotations, "probes_info_batch{}.txt".format(batch_id)
        )
        file_probe_fasta_batch = os.path.join(
            self.dir_annotations, "probes_sequence_batch{}.fna".format(batch_id)
        )

        probes_info = pd.read_csv(file_probe_info_batch, sep="\t")

        os.remove(file_probe_info_batch)

        probes_info_filtered = probes_info[
            ~probes_info["probe_sequence"].isin(self.duplicated_sequences)
        ]
        probes_info_filtered.reset_index(inplace=True, drop=True)
        probe_ids = [
            f"{g_id}_pid{i}" for i, g_id in enumerate(probes_info_filtered["gene_id"])
        ]
        probes_info_filtered.insert(0, "probe_id", probe_ids)

        if self.write_probes==True:
            self._write_probes(
                probes_info_filtered, file_probe_info_batch, file_probe_fasta_batch
            )

    def _write_probes(
        self, probes_info_filtered, file_probe_info_batch, file_probe_fasta_batch
    ):
        """Save filtered probe information in tsv file. Save probe sequences as fasta file.
        :param probes_info_filtered: Dataframe with probe information, where exact matches have been filtered out.
        :type probes_info_filtered: pandas.DataFrame
        :param file_probe_info_batch: Path to tsv file with probe info.
        :type file_probe_info_batch: string
        :param file_probe_fasta_batch: Path to fasta file with probe sequences.
        :type file_probe_fasta_batch: string
        """
        # save info table
        probes_info_filtered.to_csv(
            file_probe_info_batch, sep="\t", header=True, index=False
        )

        # save sequence of probes in fasta format
        genes = probes_info_filtered.gene_id.unique()
        for subbatch_id in range(self.number_subbatches):
            file_probe_fasta_subbatch = file_probe_fasta_batch.replace(
                ".fna", "_{}.fna".format(subbatch_id)
            )

            genes_subbatch = genes[
                (subbatch_id * self.max_genes_in_batch) : (
                    (subbatch_id + 1) * self.max_genes_in_batch
                )
            ]
            probes_info_filtered_subbatch = probes_info_filtered.loc[
                probes_info_filtered["gene_id"].isin(genes_subbatch)
            ].copy()
            probes_info_filtered_subbatch.reset_index(inplace=True, drop=True)

            output = []
            for row in probes_info_filtered_subbatch.index:
                header = probes_info_filtered_subbatch.iloc[
                    row, probes_info_filtered_subbatch.columns.get_loc("probe_id")
                ]

                sequence = Seq(
                    probes_info_filtered_subbatch.iloc[
                        row,
                        probes_info_filtered_subbatch.columns.get_loc("probe_sequence"),
                    ]
                )
                output.append(SeqRecord(sequence, header, "", ""))

            with open(file_probe_fasta_subbatch, "w") as handle:
                SeqIO.write(output, handle, "fasta")

    def filter_probes_exactmatch(self, probe_info):
        self.logging.info("Creating batches")
        self.create_batches(probe_info)

        self.duplicated_sequences = self._get_duplicated_sequences()

        # run filter with joblib

        start_time = time.perf_counter()
        with parallel_backend("loky"):
            Parallel()(
                delayed(self._filter_exactmatch_batch)(batch_id)
                for batch_id in range(self.n_jobs)
            )

        finish_time = time.perf_counter()
        self.logging.info(f"Exact matches filtered in {finish_time-start_time} seconds")'''
