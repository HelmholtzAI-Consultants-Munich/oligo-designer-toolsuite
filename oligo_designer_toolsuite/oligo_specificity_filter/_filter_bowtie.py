import os
import re
import time
from pathlib import Path
from subprocess import Popen

import pandas as pd
from joblib import Parallel, delayed, parallel_backend

from ._filter_base import ProbeFilterBase
from ._utils import (
    _load_probes_info,
    _write_output,
    _write_removed_genes,
    mismatch_in_ligation,
)

# TODO: Add ligation region


class ProbeFilterBowtie(ProbeFilterBase):
    def __init__(
        self,
        n_jobs,
        dir_output,
        dir_annotations,
        file_transcriptome_fasta,
        min_mismatches=4,
        mismatch_region=None,
        ligation_region=0,
    ):
        """This class filters probes based on the Bowtie short read alignment tool.
        The user can customize the filtering by specifying the min_mismatches per probe and mismatch_region, the region that should only be considered for counting mismatches.
        That is, all probes with number mismatches less than min_mismatches inside the mismatch_region are filtered out.

        :param file_transcriptome_fasta: path to fasta file containing all probes
        :type file_transcriptome_fasta: str
        :param min_mismatches: Threshhold value on the number of mismatches required for each probe. Probes where the number of mismatches are greater than or equal to this threshhold are considered valid. Possible values range from 0 to 4.
        :type min_mismatches: int
        :param mismatch_region: The region of the probe where the mismatches are considered. Probes that have less than min_mismatches in the first L bases (where L is a number 5 or greater) are filtered out
        :type mismatch_region: int
        :param ligation_region: coverage between probes and target sequence should not span region around ligation site (e.g. ligation_region = 5 would correspond to -4 to +5 nt around ligation site), if ligation_region = 0, omit this requirement
        :type ligation_region: int


        """
        super().__init__(n_jobs, dir_output, dir_annotations)

        self.n_jobs = n_jobs
        self.file_transcriptome_fasta = file_transcriptome_fasta
        self.ligation_region = ligation_region

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

        index_exists = False
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
        :rtype: Txt files"""

        def _process_bowtie_results(batch_id):

            probes_info = _load_probes_info(self, batch_id)

            num_probes_wo_match = 0
            for subbatch_id in range(self.number_subbatches):
                bowtie_results = _read_bowtie_output(batch_id, subbatch_id)
                num_probes_wo_match += _filter_probes_bowtie(
                    probes_info, bowtie_results
                )

        def _read_bowtie_output(batch_id, subbatch_id):
            """Load the output of the bowtie alignment search into a DataFrame and process the results.

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
                bowtie_results["query"].str.split("_").str[0]
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

            if self.ligation_region > 0:

                # Get probes where mismatch positions are given
                bowtie_results_matches_mismatch_positions = bowtie_results_matches[
                    bowtie_results_matches["mismatch_positions"].notna()
                ]
                # Extract positions
                positions = (
                    bowtie_results_matches_mismatch_positions["mismatch_positions"]
                    .str.split(",")
                    .apply(lambda x: [elem.split(":")[0] for elem in x])
                )
                # add positions to bowtie results
                bowtie_results_positions = bowtie_results_matches.merge(
                    positions, left_index=True, right_index=True, how="left"
                )

                bowtie_results_add_ligation = bowtie_results_positions.merge(
                    probes_info[["probe_id", "ligation_site"]],
                    left_on="query",
                    right_on="probe_id",
                    how="inner",
                )

                # Calculate ligation region and search if mismatch is in this region
                bowtie_results_add_ligation[
                    "ligation_region_start"
                ] = bowtie_results_add_ligation.ligation_site - (
                    self.ligation_region - 1
                )
                bowtie_results_add_ligation["ligation_region_end"] = (
                    bowtie_results_add_ligation.ligation_site + self.ligation_region
                )

                bowtie_results_add_ligation.rename(
                    columns={"mismatch_positions_y": "positions"}, inplace=True
                )
                bowtie_results_add_ligation["mismatch_in_ligation"] = [
                    False for i in range(len(bowtie_results_add_ligation))
                ]

                # Calculate number of mismatches in ligation region
                bowtie_results_add_ligation_notna = bowtie_results_add_ligation[
                    bowtie_results_add_ligation["positions"].notna()
                ]
                bowtie_results_add_ligation_notna = (
                    bowtie_results_add_ligation_notna.apply(
                        mismatch_in_ligation, axis=1
                    )
                )

                # filter out probes where there is mismatch in ligation region
                bowtie_results_matches = bowtie_results_add_ligation_notna[
                    bowtie_results_add_ligation_notna["mismatch_in_ligation"] == False
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
                    _write_output(self, probes_info, gene_id, probes_wo_match_gene)

            return len(probes_wo_match["probe_id"].unique())

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

        _write_removed_genes(self)

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

    def apply(self, probe_info):
        """Apply bowtie filter to all batches in parallel"""

        # Filter out exact matches
        self.filter_probes_exactmatch(probe_info)

        self.run_bowtie()
        self.filter_probes_by_bowtie_results()
