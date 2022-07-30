import os
import time
from pathlib import Path
from subprocess import Popen

from joblib import Parallel, delayed, parallel_backend

from ._filter_base import ProbeFilterBase


class ProbeFilterBowtie(ProbeFilterBase):
    def __init__(
        self,
        n_jobs,
        file_transcriptome_fasta,
        dir_output,
        file_probe_info,
        genes,
        probe_length_min,
        probe_length_max,
        min_mismatches=3,
        mismatch_region=None,
    ):
        """This class filters probes based on the Bowtie short read alignment tool.
        User can customize the filtering by specifying the min_mismatches per probe and mismatch_region, the region that should only be considered for counting mismatches.
        That is, all probes with number mismatches less than min_mismatches inside mismatch_region are filtered out.

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
        :param min_mismatches: minimum mismatches required for each probe to be considered valid
        :type min_mismatches: int
        :param mismatch_region: The region of the probe where the mismatches are considered
        :type mismatch_region: int


        """
        super().__init__(n_jobs, dir_output, file_probe_info, genes)
        self.n_jobs = n_jobs
        self.file_transcriptome_fasta = file_transcriptome_fasta
        self.probe_length_min = probe_length_min
        self.probe_length_max = probe_length_max
        self.min_mismatches = min_mismatches
        self.mismatch_region = mismatch_region

        self.dir_bowtie = os.path.join(self.dir_output, "bowtie")
        Path(self.dir_bowtie).mkdir(parents=True, exist_ok=True)

        self.dir_probes = os.path.join(self.dir_output, "probes")
        Path(self.dir_probes).mkdir(parents=True, exist_ok=True)

        self.file_removed_genes = os.path.join(
            self.dir_output, "genes_with_insufficient_probes.txt"
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

            for subbatch_id in range(self.n_jobs):
                file_probe_fasta_batch = os.path.join(
                    self.dir_annotations,
                    "probes_sequence_batch{}_{}.fna".format(batch_id, subbatch_id),
                )
                file_bowtie_batch = os.path.join(
                    self.dir_bowtie,
                    "bowtie_batch{}_{}.txt".format(batch_id, subbatch_id),
                )

                command = (
                    "./bowtie" + self.file_transcriptome_fasta + file_probe_fasta_batch
                )

                with open(file_bowtie_batch, "w") as file:
                    process = Popen(command, shell=True, stdout=file)

        self.logging.info("Creating blast database")

        start_time = time.perf_counter()

        with parallel_backend("loky"):
            Parallel()(
                delayed(_run_bowtie_batch)(batch_id) for batch_id in range(self.n_jobs)
            )

        finish_time = time.perf_counter()
        self.logging.info(
            f"Bowtie alignmnet search finished in {finish_time-start_time} seconds"
        )

    def filter_probes_by_bowtie_results():
        """Use bowtie results to filter probes based on number of mismatches and mismatch position

        :return: Text file for each gene containing the filtered probes for that gene. Additionally, a text file containing a list of genes with insufficient number of probes is returned.
        :rtype: Text files"""

    def apply(self):

        self.run_bowtie()
        # self.filter_probes_by_bowtie_results()
