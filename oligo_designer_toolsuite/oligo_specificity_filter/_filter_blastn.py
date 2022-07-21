import os
import re
import shutil
import time
from pathlib import Path

import pandas as pd
from Bio.Blast.Applications import NcbiblastnCommandline, NcbimakeblastdbCommandline
from joblib import Parallel, delayed, parallel_backend

from ._filter_base import ProbeFilterBase


class ProbeFilterBlastn(ProbeFilterBase):
    def __init__(
        self,
        number_batches,
        file_transcriptome_fasta,
        ligation_region,
        dir_output,
        file_probe_info,
        genes,
        n_jobs,
        word_size,
        percent_identity,
        coverage,
        probe_length_min,
        probe_length_max,
    ):
        """This class filters probes based on the blast alignment tool

        :param file_transcriptome_fasta: path to fasta file containing all probes
        :type file_transcriptome_fasta: str
        :param word_size: #word size for the blastn seed (exact match to target)
        :type word_size: int
        :param percent_identity: maximum similarity between probes and target sequences, ranging from 0 to 100% (no missmatch)
        :type percent_identity: int
        :param coverage: minimum coverage between probes and target sequence, ranging from 0 to 100% (full coverage)
        :type coverage: int
        :param probe_length_min: minimum length of the probes created
        :type probe_length_min: int
        :param probe_length_max: maximum length of the probes created
        :type probe_length_max: int

        """
        super().__init__(
            number_batches, ligation_region, dir_output, file_probe_info, genes, n_jobs
        )
        self.number_batches = number_batches
        self.word_size = word_size
        self.percent_identity = percent_identity
        self.coverage = coverage
        self.file_transcriptome_fasta = file_transcriptome_fasta
        self.probe_length_min = probe_length_min
        self.probe_length_max = probe_length_max

        self.dir_blast = os.path.join(self.dir_output, "blast")
        Path(self.dir_blast).mkdir(parents=True, exist_ok=True)

        self.dir_probes = os.path.join(self.dir_output, "probes")
        Path(self.dir_probes).mkdir(parents=True, exist_ok=True)

        self.file_removed_genes = os.path.join(
            self.dir_output, "genes_with_insufficient_probes.txt"
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
                    num_threads=1,
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
            """Parallel(n_jobs=self.n_jobs)(
                delayed(_run_blast)(batch_id)
                for batch_id in range(self.number_batches)
            )"""
            Parallel(n_jobs=self.n_jobs, prefer="threads")(
                delayed(_run_blast)(batch_id) for batch_id in range(self.number_batches)
            )
        finish_time = time.perf_counter()
        self.logging.info(
            f"Blast alignmnet search finished in {finish_time-start_time} seconds"
        )

    def filter_probes_by_blast_results(self):
        """Process the output of the BlastN alignment search
        :param batch_id: Batch ID.
        :type batch_id: int
        """

        def _process_blast_results(batch_id):

            probes_info = _load_probes_info(batch_id)

            num_probes_wo_match = 0
            for subbatch_id in range(self.number_subbatches):
                blast_results = _read_blast_output(batch_id, subbatch_id)
                num_probes_wo_match += _filter_probes_blast(probes_info, blast_results)

        def _load_probes_info(batch_id):
            """Load filtered probe information from tsv file
            :param batch_id: Batch ID.
            :type batch_id: int
            :return: Dataframe with probe information, filtered based on sequence properties.
            :rtype: pandas.DataFrame
            """
            file_probe_info_batch = os.path.join(
                self.dir_annotations, "probes_info_batch{}.txt".format(batch_id)
            )
            probes_info = pd.read_csv(
                file_probe_info_batch,
                sep="\t",
                dtype={
                    "probe_id": str,
                    "gene_id": str,
                    "probe_sequence": str,
                    "transcript_id": str,
                    "exon_id": str,
                    "chromosome": str,
                    "start": str,
                    "end": str,
                    "strand": str,
                    "GC_content": float,
                    "melting_temperature": float,
                    "melt_temp_arm1": float,
                    "melt_temp_arm2": float,
                    "dif_melt_temp_arms": float,
                    "ligation_site": int,
                },
            )

            return probes_info

        def _read_blast_output(batch_id, subbatch_id):
            """Load the output of the BlastN alignment search into a DataFrame and process the results.
            :param batch_id: Batch ID.
            :type batch_id: int
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
                blast_results["query"].str.split("_pid").str[0]
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
                    _write_output(probes_info, gene_id, probes_wo_match_gene)

            return len(probes_wo_match["query"].unique())

        def _write_output(probes_info, gene_id, probes_wo_match):
            """Write results of probe design pipeline to file and create one file with suitable probes per gene.
            :param probes_info: Dataframe with probe information, filtered based on sequence properties.
            :type probes_info: pandas.DataFrame
            :param gene_id: Gene ID of processed gene.
            :type gene_id: string
            :param probes_wo_match: List of suitable probes that don't have matches in the transcriptome.
            :type probes_wo_match: list
            """
            file_output = os.path.join(self.dir_probes, "probes_{}.txt".format(gene_id))
            valid_probes = probes_info[probes_info["probe_id"].isin(probes_wo_match)]
            valid_probes.to_csv(file_output, sep="\t", index=False)

        def _write_removed_genes():
            """Write list of genes for which not enough probes could be designed for."""

            # create file where removed genes are saved
            _, _, probe_files = next(os.walk(self.dir_probes))
            for probe_file in probe_files:
                gene_id = probe_file[len("probes_") : -len(".txt")]
                if gene_id in self.removed_genes:
                    self.removed_genes.remove(gene_id)

            with open(self.file_removed_genes, "w") as output:
                for gene_id in self.removed_genes:
                    output.write(f"{gene_id}\t0\n")

        start_time = time.perf_counter()
        with parallel_backend("loky"):
            """Parallel(n_jobs=self.n_jobs)(
                delayed(_process_blast_results)(batch_id)
                for batch_id in range(self.number_batches)
            )"""
            Parallel(n_jobs=self.n_jobs, prefer="threads")(
                delayed(_process_blast_results)(batch_id)
                for batch_id in range(self.number_batches)
            )
        finish_time = time.perf_counter()

        self.logging.info(
            f"Blast results processed in {finish_time-start_time} seconds"
        )

        _write_removed_genes()

        # remove intermediate files
        shutil.rmtree(self.dir_blast)
        for file in os.listdir(self.dir_annotations):
            if re.search("probes_*", file):
                os.remove(os.path.join(self.dir_annotations, file))

    def apply(self):

        self.run_blast_search()

        self.filter_probes_by_blast_results()
