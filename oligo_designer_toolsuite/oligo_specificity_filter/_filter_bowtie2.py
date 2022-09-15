import os
import re
import time
from pathlib import Path
from subprocess import Popen

import pandas as pd
from joblib import Parallel, delayed, parallel_backend

from ._filter_base import ProbeFilterBase
from ._utils import _load_probes_info, _write_output, _write_removed_genes

# TODO: Add optional user defined arguments for bowtie command line
# Throw error when wrong min_score format is used


class ProbeFilterBowtie2(ProbeFilterBase):
    def __init__(
        self,
        n_jobs,
        dir_output,
        dir_annotations,
        file_transcriptome_fasta,
        min_score=None,
    ):

        super().__init__(n_jobs, dir_output, dir_annotations)

        self.n_jobs = n_jobs
        self.file_transcriptome_fasta = file_transcriptome_fasta
        self.min_score = min_score

        self.dir_bowtie2 = os.path.join(self.dir_output, "bowtie2")
        Path(self.dir_bowtie2).mkdir(parents=True, exist_ok=True)

        self.dir_probes = os.path.join(self.dir_output, "probes_bowtie2")
        Path(self.dir_probes).mkdir(parents=True, exist_ok=True)

        self.file_removed_genes = os.path.join(
            self.dir_probes, "genes_with_insufficient_probes.txt"
        )

    def run_bowtie2(self):
        """Run Bowtie alignment tool to find regions of local similarity between sequences, where sequences are probes and transcripts.
        Bowtie identifies all allignments between the probes and transcripts and returns the number of mismatches and mismatch position for each alignment.

        :return: DataFrame with processed bowtie alignment search results.
        :rtype: pandas.DataFrame
        """

        def _run_bowtie2_batch(batch_id):
            """Run Bowtie alignment tool for all probes of one batch
            :param batch_id: Batch ID.
            :type batch_id: int
            """

            for subbatch_id in range(self.number_subbatches):
                file_probe_fasta_batch = os.path.join(
                    self.dir_annotations,
                    "probes_sequence_batch{}_{}.fna".format(batch_id, subbatch_id),
                )
                file_bowtie2_batch = os.path.join(
                    self.dir_bowtie2,
                    "bowtie2_batch{}_{}.sam".format(batch_id, subbatch_id),
                )

                if self.min_score is not None:
                    command = (
                        "bowtie2 -f -a --score-min "
                        + self.min_score
                        + " --no-hd -x transcriptome_bowtie2"
                        + " -U "
                        + file_probe_fasta_batch
                        + " -S "
                        + file_bowtie2_batch
                    )
                else:
                    command = (
                        "bowtie2 -f -a "
                        + " --no-hd -x transcriptome_bowtie2"
                        + " -U "
                        + file_probe_fasta_batch
                        + " -S "
                        + file_bowtie2_batch
                    )

                process = Popen(command, shell=True, cwd=self.dir_bowtie2).wait()

        # Check if bowtie2 index exists
        index_exists = False
        for file in os.listdir(self.dir_bowtie2):
            if re.search("^transcriptome_bowtie2*", file):
                index_exists = True
                break

        # Create bowtie index if none exists
        if not index_exists:
            command1 = (
                "bowtie2-build -f "
                + self.file_transcriptome_fasta
                + " transcriptome_bowtie2"
            )
            process = Popen(command1, shell=True, cwd=self.dir_bowtie2).wait()

        start_time = time.perf_counter()

        with parallel_backend("loky"):
            Parallel()(
                delayed(_run_bowtie2_batch)(batch_id) for batch_id in range(self.n_jobs)
            )

        finish_time = time.perf_counter()
        self.logging.info(
            f"Bowtie2 alignmnet search finished in {finish_time-start_time} seconds"
        )

    def filter_probes_by_bowtie2_results(self):
        """Use bowtie 2 results to filter probes based on min_score parameter

        :return: Text file for each gene containing the filtered probes for that gene. Additionally, a text file containing a list of genes with insufficient number of probes is returned.
        :rtype: Txt files"""

        def _process_bowtie2_results(batch_id):

            probes_info = _load_probes_info(self, batch_id)

            num_probes_wo_match = 0
            for subbatch_id in range(self.number_subbatches):
                bowtie2_results = _read_bowtie2_output(batch_id, subbatch_id)
                num_probes_wo_match += _filter_probes_bowtie2(
                    probes_info, bowtie2_results
                )

        def _read_bowtie2_output(batch_id, subbatch_id):
            """Load the output of the bowtie 2 alignment search into a DataFrame and process the results.

            :return: DataFrame with processed bowtie 2 alignment search results.
            :rtype: pandas.DataFrame
            """
            file_bowtie2_batch = os.path.join(
                self.dir_bowtie2, "bowtie2_batch{}_{}.sam".format(batch_id, subbatch_id)
            )
            bowtie2_results = pd.read_csv(
                file_bowtie2_batch,
                header=None,
                sep="\t",
                low_memory=False,
                names=[
                    "query",
                    "flags",
                    "reference",
                    "ref_start",
                    "mapping_quality",
                    "CIGAR_alignment",
                    "ref_seq_mate",
                    "offset",
                    "fragment_length",
                    "sequence",
                    "read_qualities",
                    "alignment_score",
                    "_1",
                    "_2",
                    "_3",
                    "_4",
                    "_5",
                    "_6",
                    "_7",
                    "_8",
                    "_9",
                    "_10",
                ],
                engine="c",
                dtype={
                    "query": str,
                    "flags": str,
                    "reference": str,
                    "ref_start": str,
                    "mapping_quality": str,
                    "CIGAR_alignment": str,
                    "ref_seq_mate": str,
                    "offset": str,
                    "fragment_length": str,
                    "sequence": str,
                    "read_qualities": str,
                    "alignment_score": str,
                    "_1": str,
                    "_2": str,
                    "_3": str,
                    "_4": str,
                    "_5": str,
                    "_6": str,
                    "_7": str,
                    "_8": str,
                    "_9": str,
                    "_10": str,
                },
            )

            bowtie2_results["query_gene_id"] = (
                bowtie2_results["query"].str.split("_").str[0]
            )
            bowtie2_results["reference_gene_id"] = (
                bowtie2_results["reference"].str.split("::").str[0]
            )

            return bowtie2_results

        def _filter_probes_bowtie2(probes_info, bowtie2_results):
            """Use the results of the Bowtie 2 alignement search to remove probes with high similarity (i.e. low number of mismatches represented by a high alignment score) based on user-defined min_score.
            :param probes_info: Dataframe with probe information, filtered based on sequence properties.
            :type probes_info: pandas.DataFrame
            :param blast_results: DataFrame with processed bowtie alignment search results.
            :type blast_results: pandas.DataFrame
            """

            bowtie2_results_matches = bowtie2_results[
                bowtie2_results["query_gene_id"] != bowtie2_results["reference_gene_id"]
            ]

            probes_with_match = bowtie2_results_matches["query"].unique()
            probes_wo_match = probes_info[
                ~probes_info["probe_id"].isin(probes_with_match)
            ]

            for gene_id in probes_info["gene_id"].unique():
                probes_wo_match_gene = probes_wo_match[
                    probes_wo_match.gene_id == gene_id
                ]
                probes_wo_match_gene = probes_wo_match_gene["probe_id"].unique()

                if len(probes_wo_match_gene) > 0:  # gene has to have at least one probe
                    _write_output(self, probes_info, gene_id, probes_wo_match_gene)

            return len(probes_wo_match["probe_id"].unique())

        # Parallelize filtering of probes
        start_time = time.perf_counter()
        with parallel_backend("loky"):
            Parallel()(
                delayed(_process_bowtie2_results)(batch_id)
                for batch_id in range(self.n_jobs)
            )

        finish_time = time.perf_counter()

        self.logging.info(
            f"Bowtie2 results processed in {finish_time-start_time} seconds"
        )

        _write_removed_genes(self)

        # remove intermediate files
        for file in os.listdir(self.dir_bowtie2):
            if re.search("bowtie2_*", file):
                os.remove(os.path.join(self.dir_bowtie2, file))

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
        """Apply bowtie 2 filter to all batches in parallel"""

        # Filter out exact matches
        self.filter_probes_exactmatch(probe_info)

        self.run_bowtie2()
        self.filter_probes_by_bowtie2_results()
