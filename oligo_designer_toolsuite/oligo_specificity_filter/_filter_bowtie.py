import os
import re
from pathlib import Path
from subprocess import Popen

import pandas as pd
from joblib import Parallel, delayed

from . import SpecificityFilterBase


class Bowtie(SpecificityFilterBase):
    def __init__(self, dir_specificity, min_mismatches=4, mismatch_region=None):
        """This class filters probes based on the Bowtie short read alignment tool.
        The user can customize the filtering by specifying the min_mismatches per probe and mismatch_region, the region that should be considered for counting mismatches.
        That is, all probes with number mismatches less than min_mismatches inside the mismatch_region are filtered out.

        Use conda install -c bioconda bowtie to install Bowtie package

        :param file_transcriptome_fasta: path to fasta file containing all probes
        :type file_transcriptome_fasta: str
        :param min_mismatches: Threshhold value on the number of mismatches required for each probe. Probes where the number of mismatches are greater than or equal to this threshhold are considered valid. Possible values range from 0 to 4.
        :type min_mismatches: int
        :param mismatch_region: The region of the probe where the mismatches are considered. Probes that have less than min_mismatches in the first L bases (where L is 5 or greater) are filtered out
        :type mismatch_region: int


        """
        super().__init__(dir_specificity)

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

        self.dir_bowtie = os.path.join(self.dir_specificity, "bowtie")
        Path(self.dir_bowtie).mkdir(parents=True, exist_ok=True)

        self.dir_fasta = os.path.join(self.dir_specificity, "fasta")
        Path(self.dir_fasta).mkdir(parents=True, exist_ok=True)

    def apply(self, oligo_DB, file_reference_DB, n_jobs):
        """Apply bowtie filter to all batches in parallel"""

        # Some bowtie initializations, change the names
        index_exists = False
        # Check if bowtie index exists
        for file in os.listdir(self.dir_bowtie):
            if re.search("^reference*", file):
                index_exists = True
                break

        # Create bowtie index if none exists
        if not index_exists:
            command1 = "bowtie-build " + file_reference_DB + " reference"
            process = Popen(command1, shell=True, cwd=self.dir_bowtie).wait()

        genes = list(oligo_DB.keys())
        filtered_oligo_DBs = Parallel(n_jobs=n_jobs)(
            delayed(self._run_bowtie)(oligo_DB[gene], gene) for gene in genes
        )

        # reconstruct the oligos_DB
        for gene, filtered_oligo_DB in zip(genes, filtered_oligo_DBs):
            oligo_DB[gene] = filtered_oligo_DB
        # remove the files
        for file in os.listdir(self.dir_bowtie):
            os.remove(os.path.join(self.dir_bowtie, file))
        for file in os.listdir(self.dir_fasta):
            os.remove(os.path.join(self.dir_fasta, file))
        return oligo_DB

    def _run_bowtie(self, gene_DB, gene):
        """Run Bowtie alignment tool to find regions of local similarity between sequences, where sequences are probes and transcripts.
        Bowtie identifies all allignments between the probes and transcripts and returns the number of mismatches and mismatch position for each alignment.

        :return: DataFrame with processed bowtie alignment search results.
        :rtype: pandas.DataFrame
        """

        file_probe_fasta_gene = self._create_fasta_file(gene_DB, self.dir_fasta, gene)
        file_bowtie_gene = os.path.join(
            self.dir_bowtie,
            f"bowtie_{gene}.txt",
        )
        if self.mismatch_region is not None:
            command = (
                "bowtie reference -f -a -n "
                + str(self.min_mismatches - 1)
                + "-l"
                + str(self.mismatch_region)
                + " "
                + file_probe_fasta_gene
                + " "
                + file_bowtie_gene
            )
        else:
            command = (
                "bowtie reference -f -a -v "
                + str(self.min_mismatches - 1)
                + " "
                + file_probe_fasta_gene
                + " "
                + file_bowtie_gene
            )

        process = Popen(command, shell=True, cwd=self.dir_bowtie).wait()

        # read the results of the bowtie search
        bowtie_results = self._read_bowtie_output(file_bowtie_gene)
        # filter the DB based on the bowtie results
        matching_probes = self._find_matching_probes(gene_DB, bowtie_results)
        filtered_gene_DB = self._filter_matching_probes(gene_DB, matching_probes)
        return filtered_gene_DB

    def _read_bowtie_output(self, file_bowtie_gene):
        """Load the output of the bowtie alignment search into a DataFrame and process the results.

        :return: DataFrame with processed bowtie alignment search results.
        :rtype: pandas.DataFrame
        """
        bowtie_results = pd.read_csv(
            file_bowtie_gene,
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
        # return the real matches, that is the ones not belonging to the same gene of the query probe
        bowtie_results["query_gene_id"] = bowtie_results["query"].str.split("_").str[0]
        bowtie_results["reference_gene_id"] = (
            bowtie_results["reference"].str.split("::").str[0]
        )

        return bowtie_results

    def _find_matching_probes(self, gene_DB, bowtie_results):
        """Use the results of the Bowtie alignment search to identify probes with high similarity (i.e. low number of mismatches) based on user defined thresholds.
        :param probes_info: Dataframe with probe information, filtered based on sequence properties.
        :type probes_info: pandas.DataFrame
        :param blast_results: DataFrame with processed bowtie alignment search results.
        :type blast_results: pandas.DataFrame
        """
        bowtie_matches = bowtie_results[
            bowtie_results["query_gene_id"] != bowtie_results["reference_gene_id"]
        ]
        probes_with_match = bowtie_matches["query"].unique()
        return probes_with_match

    def _filter_matching_probes(self, gene_DB, matching_probes):
        probe_ids = list(gene_DB.keys())
        for probe_id in probe_ids:
            if probe_id in matching_probes:
                del gene_DB[probe_id]
        return gene_DB
