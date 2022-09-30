import os
import time

import iteration_utilities
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from joblib import Parallel, delayed, parallel_backend

from oligo_designer_toolsuite.IO._data_parser import read_oligos_DB_tsv

from ._filter_base import ProbeFilterBase


class ProbeFilterExact(ProbeFilterBase):
    def __init__(self, n_jobs, dir_output, dir_annotations):
        """This class filters probes based on exact matches."""

        super().__init__(n_jobs, dir_output, dir_annotations)

    def _get_duplicated_sequences(self):
        """Get a list of probe sequences that have an exact match within the pool of all
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

        probes_info = pd.read_csv(file_probe_info_batch, sep="\t")

        os.remove(file_probe_info_batch)

        probes_info_filtered = probes_info[
            ~probes_info["probe_sequence"].isin(self.duplicated_sequences)
        ]
        probes_info_filtered.reset_index(inplace=True, drop=True)

        for subbatch_id in range(self.number_subbatches):
            file_probe_fasta_batch = os.path.join(
                self.dir_annotations,
                "probes_sequence_batch{}_{}.fna".format(batch_id, subbatch_id),
            )

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

    def apply(self, probe_info):
        """Parallelize filtering of exact matches from input dictionary"""

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

        # Get probe dictionary
        file_probe_info_batch = os.path.join(
            self.dir_annotations, "probes_info_batch0.txt"
        )
        probe_info_combined = pd.read_csv(file_probe_info_batch, sep="\t")

        for batch_id in range(1, self.n_jobs):
            file_probe_info_batch = os.path.join(
                self.dir_annotations, "probes_info_batch{}.txt".format(batch_id)
            )
            probe_info_batch = pd.read_csv(file_probe_info_batch, sep="\t")
            probe_info_combined.append(probe_info_batch)

        # Check that column order is preserved
        assert probe_info_combined.columns == list(probe_info[0][0].keys())
        file_probe_info_dict = read_oligos_DB_tsv(probe_info_combined)

        finish_time = time.perf_counter()
        self.logging.info(f"Exact matches filtered in {finish_time-start_time} seconds")

        return file_probe_info_dict
