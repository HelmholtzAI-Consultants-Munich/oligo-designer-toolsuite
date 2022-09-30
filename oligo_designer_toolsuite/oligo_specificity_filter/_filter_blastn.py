import os
import re
import shutil
import time
from pathlib import Path

import pandas as pd
from Bio.Blast.Applications import NcbiblastnCommandline, NcbimakeblastdbCommandline
from joblib import Parallel, delayed, parallel_backend

from ._filter_base import ProbeFilterBase
from ._utils import _load_probes_info, _write_output, _write_removed_genes


class ProbeFilterBlastn(ProbeFilterBase):
    def __init__(
        self,
        n_jobs,
        file_transcriptome_fasta,
        dir_output,
        dir_annotations,
        word_size,
        percent_identity,
        coverage,
        probe_length_min,
        probe_length_max,
        ligation_region=0,
    ):
        """This class filters probes based on the blast alignment tool

        :param file_transcriptome_fasta: path to fasta file containing all probes
        :type file_transcriptome_fasta: str
        :param word_size: word size for the blastn seed (exact match to target)
        :type word_size: int
        :param percent_identity: maximum similarity between probes and target sequences, ranging from 0 to 100% (no missmatch)
        :type percent_identity: int
        :param coverage: minimum coverage between probes and target sequence, ranging from 0 to 100% (full coverage)
        :type coverage: int
        :param probe_length_min: minimum length of the probes created
        :type probe_length_min: int
        :param probe_length_max: maximum length of the probes created
        :type probe_length_max: int
        :param ligation_region: coverage between probes and target sequence should not span region around ligation site (e.g. ligation_region = 5 would correspond to -4 to +5 nt around ligation site), if ligation_region = 0, omit this requirement
        :type ligation_region: int
        """
        super().__init__(n_jobs, dir_output, dir_annotations)

        self.word_size = word_size
        self.percent_identity = percent_identity
        self.coverage = coverage
        self.file_transcriptome_fasta = file_transcriptome_fasta
        self.probe_length_min = probe_length_min
        self.probe_length_max = probe_length_max
        self.ligation_region = ligation_region

        self.dir_blast = os.path.join(self.dir_output, "blast")
        Path(self.dir_blast).mkdir(parents=True, exist_ok=True)

        self.dir_probes = os.path.join(self.dir_output, "probes_blast")
        Path(self.dir_probes).mkdir(parents=True, exist_ok=True)

        self.file_removed_genes = os.path.join(
            self.dir_probes, "genes_with_insufficient_probes.txt"
        )

    def run_blast_search(self):
        """Run BlastN alignment tool to find regions of local similarity between sequences, where sequences are probes and transcripts.
        BlastN identifies the transcript regions where probes match with a certain coverage and similarity.
        """

        def _run_blast(batch_id):
            """Run BlastN alignment search for all probes of one batch
            :param batch_id: Batch ID.
            :type batch_id: int
            """
            for subbatch_id in range(self.number_subbatches):
                file_probe_fasta_batch = os.path.join(
                    self.dir_annotations,
                    "probes_sequence_batch{}_{}.fna".format(batch_id, subbatch_id),
                )
                file_blast_batch = os.path.join(
                    self.dir_blast, "blast_batch{}_{}.txt".format(batch_id, subbatch_id)
                )

                cmd = NcbiblastnCommandline(
                    query=file_probe_fasta_batch,
                    db=self.file_transcriptome_fasta,
                    outfmt="10 qseqid sseqid length qstart qend qlen",
                    out=file_blast_batch,
                    strand="plus",
                    word_size=self.word_size,
                    perc_identity=self.percent_identity,
                    num_threads=self.n_jobs,
                )
                out, err = cmd()

        # create blast database
        self.logging.info("Creating blast database")
        start_time = time.perf_counter()
        cmd = NcbimakeblastdbCommandline(
            input_file=self.file_transcriptome_fasta, dbtype="nucl"
        )
        out, err = cmd()

        finish_time = time.perf_counter()

        self.logging.info(f"Blast database created in {finish_time-start_time} seconds")

        # run blast in parallel
        start_time = time.perf_counter()
        with parallel_backend("loky"):
            Parallel()(delayed(_run_blast)(batch_id) for batch_id in range(self.n_jobs))

        finish_time = time.perf_counter()
        self.logging.info(
            f"Blast alignmnet search finished in {finish_time-start_time} seconds"
        )

    def filter_probes_by_blast_results(self):
        """Process the output of the BlastN alignment search"""

        def _process_blast_results(batch_id):

            probes_info = _load_probes_info(self, batch_id)

            num_probes_wo_match = 0
            for subbatch_id in range(self.number_subbatches):
                blast_results = _read_blast_output(batch_id, subbatch_id)
                num_probes_wo_match += _filter_probes_blast(probes_info, blast_results)

        def _read_blast_output(batch_id, subbatch_id):
            """Load the output of the BlastN alignment search into a DataFrame and process the results.
            :return: DataFrame with processed blast alignment search results.
            :rtype: pandas.DataFrame
            """
            file_blast_batch = os.path.join(
                self.dir_blast, "blast_batch{}_{}.txt".format(batch_id, subbatch_id)
            )
            blast_results = pd.read_csv(
                file_blast_batch,
                header=None,
                sep=",",
                low_memory=False,
                names=[
                    "query",
                    "target",
                    "alignment_length",
                    "query_start",
                    "query_end",
                    "query_length",
                ],
                engine="c",
                dtype={
                    "query": str,
                    "target": str,
                    "alignment_length": int,
                    "query_start": int,
                    "query_end": int,
                    "query_length": int,
                },
            )

            blast_results["query_gene_id"] = (
                blast_results["query"].str.split("_").str[0]
            )
            blast_results["target_gene_id"] = (
                blast_results["target"].str.split("::").str[0]
            )

            return blast_results

        def _filter_probes_blast(probes_info, blast_results):
            """Use the results of the BlastN alignement search to remove probes with high similarity,
            probe coverage and ligation site coverage based on user defined thresholds.
            :param probes_info: Dataframe with probe information, filtered based on sequence properties.
            :type probes_info: pandas.DataFrame
            :param blast_results: DataFrame with processed blast alignment search results.
            :type blast_results: pandas.DataFrame
            """
            blast_results_matches = blast_results[
                blast_results["query_gene_id"] != blast_results["target_gene_id"]
            ]
            blast_results_matches = blast_results_matches.merge(
                probes_info[["probe_id", "ligation_site"]],
                left_on="query",
                right_on="probe_id",
                how="inner",
            )
            blast_results_matches.drop(columns=["probe_id"], inplace=True)
            blast_results_matches_filtered = []

            for probe_length in range(self.probe_length_min, self.probe_length_max + 1):

                min_alignment_length = probe_length * self.coverage / 100
                blast_results_matches_probe_length = blast_results_matches[
                    blast_results_matches.query_length == probe_length
                ]

                if self.ligation_region == 0:
                    blast_results_matches_probe_length = (
                        blast_results_matches_probe_length[
                            blast_results_matches_probe_length.alignment_length
                            > min_alignment_length
                        ]
                    )
                else:
                    blast_results_matches_probe_length = (
                        blast_results_matches_probe_length[
                            (
                                blast_results_matches_probe_length.alignment_length
                                > min_alignment_length
                            )
                            | (
                                (
                                    blast_results_matches_probe_length.query_start
                                    < blast_results_matches_probe_length.ligation_site
                                    - (self.ligation_region - 1)
                                )
                                & (
                                    blast_results_matches_probe_length.query_end
                                    > (
                                        blast_results_matches_probe_length.ligation_site
                                        + self.ligation_region
                                    )
                                )
                            )
                        ]
                    )

                blast_results_matches_filtered.append(
                    blast_results_matches_probe_length
                )

            blast_results_matches_filtered = pd.concat(blast_results_matches_filtered)

            probes_with_match = blast_results_matches_filtered["query"].unique()
            probes_wo_match = blast_results[
                ~blast_results["query"].isin(probes_with_match)
            ]

            for gene_id in blast_results["query_gene_id"].unique():
                probes_wo_match_gene = probes_wo_match[
                    probes_wo_match.query_gene_id == gene_id
                ]
                probes_wo_match_gene = probes_wo_match_gene["query"].unique()

                if len(probes_wo_match_gene) > 0:  # gene has to have at least one probe
                    _write_output(self, probes_info, gene_id, probes_wo_match_gene)

            return len(probes_wo_match["query"].unique())

        # Parallelize filtering of probes
        start_time = time.perf_counter()
        with parallel_backend("loky"):
            Parallel()(
                delayed(_process_blast_results)(batch_id)
                for batch_id in range(self.n_jobs)
            )

        finish_time = time.perf_counter()

        self.logging.info(
            f"Blast results processed in {finish_time-start_time} seconds"
        )

        _write_removed_genes(self)

        # remove intermediate files
        shutil.rmtree(self.dir_blast)
        for file in os.listdir(self.dir_annotations):
            if re.search("probes_*", file):
                os.remove(os.path.join(self.dir_annotations, file))

        for file in os.listdir(self.dir_annotations):
            if re.search(
                "reference_DB_unknown_unknown_Custom_release_unknown_genome_False_gene_transcript_True.n*",
                file,
            ):
                os.remove(os.path.join(self.dir_annotations, file))

    def apply(self, probe_info):
        """Apply blastN filter to all batches in parallel"""

        # # Filter out exact matches
        # self.filter_probes_exactmatch(probe_info)

        self.logging.info("Creating batches")
        self.create_batches(probe_info)

        self.run_blast_search()

        self.filter_probes_by_blast_results()
