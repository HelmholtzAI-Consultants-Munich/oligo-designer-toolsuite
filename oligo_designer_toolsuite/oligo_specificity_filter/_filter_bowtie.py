import os
import re
import time
from pathlib import Path
from subprocess import Popen

import pandas as pd
from joblib import Parallel, delayed, parallel_backend

from ._filter_base import ProbeFilterBase


class ProbeFilterBowtie(ProbeFilterBase):
    def __init__(
        self,
        n_jobs,
        dir_output,
        file_probe_info,
        genes,
        probe_length_min,
        probe_length_max,
        file_transcriptome_fasta,
        min_mismatches=4,
        mismatch_region=None,
    ):
        """This class filters probes based on the Bowtie short read alignment tool.
        User can customize the filtering by specifying the min_mismatches per probe and mismatch_region, the region that should only be considered for counting mismatches.
        That is, all probes with number mismatches less than min_mismatches inside the mismatch_region are filtered out.

        :param file_transcriptome_fasta: path to fasta file containing all probes
        :type file_transcriptome_fasta: str
        :param min_mismatches: Threshhold value on the number of mismatches required for each probe. Probes where the number of mismatches are greater than or equal to this threshhold are considered valid. Possible values range from 0 to 4.
        :type min_mismatches: int
        :param mismatch_region: The region of the probe where the mismatches are considered. Probes that have less than min_mismatches in the first L bases (where L is a number 5 or greater) are filtered out
        :type mismatch_region: int


        """
        super().__init__(n_jobs, dir_output, file_probe_info, genes)

        self.n_jobs = n_jobs
        self.file_transcriptome_fasta = file_transcriptome_fasta
        self.probe_length_min = probe_length_min
        self.probe_length_max = probe_length_max
        self.genes = genes

        if min_mismatches > 4:
            raise ValueError(
                "Choice of min_mismatches out of range for bowtie allignment tool. Please choose a value no greater than 4"
            )
        else:
            self.min_mismatches = min_mismatches

        if mismatch_region < 5:
            raise ValueError(
                "Choice of mismatch_region out of range for bowtie allignment tool. Please choose a value no less than 5"
            )
        else:
            self.mismatch_region = mismatch_region

        self.dir_bowtie = os.path.join(self.dir_output, "bowtie")
        Path(self.dir_bowtie).mkdir(parents=True, exist_ok=True)

        self.dir_probes = os.path.join(self.dir_output, "probes_bowtie")
        Path(self.dir_probes).mkdir(parents=True, exist_ok=True)

        self.file_removed_genes = os.path.join(
            self.dir_probes, "genes_with_insufficient_probes.txt"
        )

    def run_bowtie(self):
        """Run Bowtie alignment tool to find regions of local similarity between sequences, where sequences are probes and transcripts.
        Bowtie identifies all allignments between the probes and transcripts and returns the number of mismatches and mismatch position for each alignment.

        :return: DataFrame with processed bowtie alignment search results.
        :rtype: pandas.DataFrame
        """

        def _run_bowtie_batch(batch_id):
            """Run Bowtie alignment tool for all probes of one batch
            :param batch_id: Batch ID.
            :type batch_id: int
            """

            for subbatch_id in range(self.number_subbatches):
                file_probe_fasta_batch = os.path.join(
                    self.dir_annotations,
                    "probes_sequence_batch{}_{}.fna".format(batch_id, subbatch_id),
                )
                file_bowtie_batch = os.path.join(
                    self.dir_bowtie,
                    "bowtie_batch{}_{}.txt".format(batch_id, subbatch_id),
                )

                if self.mismatch_region is not None:
                    command = (
                        "bowtie transcriptome -f -a -n "
                        + str(self.min_mismatches - 1)
                        + "-l"
                        + str(self.mismatch_region)
                        + " "
                        + file_probe_fasta_batch
                        + " "
                        + file_bowtie_batch
                    )
                else:
                    command = (
                        "bowtie transcriptome -f -a -v "
                        + str(self.min_mismatches - 1)
                        + " "
                        + file_probe_fasta_batch
                        + " "
                        + file_bowtie_batch
                    )

                process = Popen(command, shell=True, cwd=self.dir_bowtie).wait()

        # Check if bowtie index exists
        for file in os.listdir(self.dir_bowtie):
            if re.search("^transcriptome*", file):
                index_exists = True
                break

        # Create bowtie index if none exists
        if not index_exists:
            command1 = (
                "bowtie-build " + self.file_transcriptome_fasta + " transcriptome"
            )
            process = Popen(command1, shell=True, cwd=self.dir_bowtie).wait()

        start_time = time.perf_counter()

        with parallel_backend("loky"):
            Parallel()(
                delayed(_run_bowtie_batch)(batch_id) for batch_id in range(self.n_jobs)
            )

        finish_time = time.perf_counter()
        self.logging.info(
            f"Bowtie alignmnet search finished in {finish_time-start_time} seconds"
        )

    def filter_probes_by_bowtie_results(self):
        """Use bowtie results to filter probes based on number of mismatches and mismatch position

        :return: Text file for each gene containing the filtered probes for that gene. Additionally, a text file containing a list of genes with insufficient number of probes is returned.
        :rtype: Text files"""

        def _process_bowtie_results(batch_id):

            probes_info = _load_probes_info(batch_id)

            num_probes_wo_match = 0
            for subbatch_id in range(self.number_subbatches):
                bowtie_results = _read_bowtie_output(batch_id, subbatch_id)
                num_probes_wo_match += _filter_probes_bowtie(
                    probes_info, bowtie_results
                )

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

        def _read_bowtie_output(batch_id, subbatch_id):
            """Load the output of the bowtie alignment search into a DataFrame and process the results.
            :param batch_id: Batch ID.
            :type batch_id: int
            :return: DataFrame with processed bowtie alignment search results.
            :rtype: pandas.DataFrame
            """
            file_bowtie_batch = os.path.join(
                self.dir_bowtie, "bowtie_batch{}_{}.txt".format(batch_id, subbatch_id)
            )
            bowtie_results = pd.read_csv(
                file_bowtie_batch,
                header=None,
                sep="\t",
                low_memory=False,
                names=[
                    "query",
                    "strand",
                    "reference",
                    "ref_start",
                    "query_sequence",
                    "read_quality",
                    "num_instances",
                    "mismatch_positions",
                ],
                engine="c",
                dtype={
                    "query": str,
                    "strand": str,
                    "reference": str,
                    "ref_start": int,
                    "query_sequence": str,
                    "read_quality": str,
                    "num_instances": int,
                    "mismatch_positions": str,
                },
            )

            bowtie_results["query_gene_id"] = (
                bowtie_results["query"].str.split("_pid").str[0]
            )
            bowtie_results["reference_gene_id"] = (
                bowtie_results["reference"].str.split("::").str[0]
            )

            return bowtie_results

        def _filter_probes_bowtie(probes_info, bowtie_results):
            """Use the results of the Bowtie alignement search to remove probes with high similarity (i.e. low number of mismatches) based on user defined thresholds.
            :param probes_info: Dataframe with probe information, filtered based on sequence properties.
            :type probes_info: pandas.DataFrame
            :param blast_results: DataFrame with processed bowtie alignment search results.
            :type blast_results: pandas.DataFrame
            """

            bowtie_results_matches = bowtie_results[
                bowtie_results["query_gene_id"] != bowtie_results["reference_gene_id"]
            ]

            probes_with_match = bowtie_results_matches["query"].unique()
            probes_wo_match = probes_info[
                ~probes_info["probe_id"].isin(probes_with_match)
            ]

            for gene_id in probes_info["gene_id"].unique():
                probes_wo_match_gene = probes_wo_match[
                    probes_wo_match.gene_id == gene_id
                ]
                probes_wo_match_gene = probes_wo_match_gene["probe_id"].unique()

                if len(probes_wo_match_gene) > 0:  # gene has to have at least one probe
                    _write_output(probes_info, gene_id, probes_wo_match_gene)

            return len(probes_wo_match["probe_id"].unique())

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

        # Parallelize filtering of probes
        start_time = time.perf_counter()
        with parallel_backend("loky"):
            Parallel()(
                delayed(_process_bowtie_results)(batch_id)
                for batch_id in range(self.n_jobs)
            )

        finish_time = time.perf_counter()

        self.logging.info(
            f"Bowtie results processed in {finish_time-start_time} seconds"
        )

        _write_removed_genes()

        # remove intermediate files
        for file in os.listdir(self.dir_bowtie):
            if re.search("bowtie_*", file):
                os.remove(os.path.join(self.dir_bowtie, file))

        for file in os.listdir(self.dir_annotations):
            if re.search("probes_*", file):
                os.remove(os.path.join(self.dir_annotations, file))

        for file in os.listdir(self.dir_annotations):
            if re.search(
                "reference_DB_unknown_unknown_Custom_release_unknown_genome_False_gene_transcript_True.n*",
                file,
            ):
                os.remove(os.path.join(self.dir_annotations, file))

    def apply(self):
        self.run_bowtie()
        self.filter_probes_by_bowtie_results()
