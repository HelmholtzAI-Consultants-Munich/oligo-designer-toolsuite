import os
import time

import iteration_utilities
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from joblib import Parallel, delayed, parallel_backend

from ._filter_base import ProbeFilterBase


class ProbeFilterExact(ProbeFilterBase):
    """This class filters probes based on exact matches. This process can be executed in a parallel fashion where no_batches = number threads.

    :param"""

    def __init__(
        self,
        number_batches,
        ligation_region,
        dir_output,
        file_probe_info,
        genes,
        dir_annotations,
    ):
        super().__init__(
            number_batches,
            ligation_region,
            dir_output,
            file_probe_info,
            genes,
            dir_annotations,
        )

        self.duplicated_sequences = None

    def _get_duplicated_sequences(self, number_batches):
        """Get a list of probe sequences that have a exact match within the pool of all
        possible probe sequences for the list of input genes.
        :return: List of probe sequences with exact matches in the pool of probes.
        :rtype: list
        """
        sequences = []

        for batch_id in range(number_batches):
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

    def _filter_probes_exactmatch(self, batch_id):
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

        self._write_probes(
            probes_info_filtered, file_probe_info_batch, file_probe_fasta_batch
        )

    def _write_probes(
        self, probes_info_filtered, file_probe_info_batch, file_probe_fasta_batch
    ):
        """Save filtered probe information in tsv file. Save probe sequences as fasta file.
        :param probes_info_filtered: Dataframe with probe information, filtered based on sequence properties.
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
        self.number_subbatches = (len(genes) // self.max_genes_in_batch) + 1
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

    def apply(self, number_batches):
        self.create_batches(number_batches)

        self.duplicated_sequences = self._get_duplicated_sequences(number_batches)

        # run filter with joblib

        start_time = time.perf_counter()
        with parallel_backend("loky"):
            Parallel(n_jobs=2)(
                delayed(self._filter_probes_exactmatch)(batch_id)
                for batch_id in range(self.number_batches)
            )
            # Parallel(n_jobs=3, prefer="threads")(delayed(self._filter_probes_exactmatch)(batch_id) for batch_id in range(self.number_batches))
        finish_time = time.perf_counter()
        print(f"Program finished in {finish_time-start_time} seconds")

        """proc = multiprocessing.Process(
            target=self._filter_probes_exactmatch, args=(batch_id,)
        )
        jobs.append(proc)
        proc.start()

    for job in jobs:
        job.join()"""
