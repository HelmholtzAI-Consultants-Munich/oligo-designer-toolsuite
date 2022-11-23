import os
import re
from pathlib import Path
from subprocess import Popen

import pandas as pd
from joblib import Parallel, delayed

from . import SpecificityFilterBase


class BowtieSeedRegion(SpecificityFilterBase):
    """This class filters probes based on the Bowtie short read alignment tool on a specific region of the probe. The region take in consideration is created according to the ``seed_region_creation`` class.
    The user can customize the filtering by specifying the num_mismatches, and all probes with number mismatches lower or equal to min_mismatches inside the mismatch_region are filtered out.

    Use conda install -c bioconda bowtie to install Bowtie package

    :param file_transcriptome_fasta: path to fasta file containing all probes
    :type file_transcriptome_fasta: str
    :param seed_region_creation: The region of the probe where the mismatches are considered. Probes that have less than min_mismatches in the first L bases (where L is 5 or greater) are filtered out
    :type seed_region_creation: SeedRegionCreationBase class
    :param num_mismatches: Threshhold value on the number of mismatches required for each probe. Probes where the number of mismatches are greater than or equal to this threshhold are considered valid. Possible values range from 0 to 4.
    :type num_mismatches: int
    """

    def __init__(
        self, dir_specificity: str, seed_region_creation, num_mismatches: int = 0
    ):
        """Constructor."""
        super().__init__(dir_specificity)

        if num_mismatches > 3:
            raise ValueError(
                "Choice of num_mismatches out of range for bowtie allignment tool. Please choose a value no greater than 3"
            )
        else:
            self.num_mismatches = num_mismatches

        self.seed_region_creation = seed_region_creation

        self.dir_seed_region = os.path.join(
            self.dir_specificity, "bowtie"
        )  # ecplot some possible already computed indices
        Path(self.dir_seed_region).mkdir(parents=True, exist_ok=True)

        self.dir_fasta = os.path.join(self.dir_specificity, "fasta")
        Path(self.dir_fasta).mkdir(parents=True, exist_ok=True)

    def apply(self, oligo_DB, file_reference_DB, n_jobs):
        """Apply bowtie filter to all batches in parallel"""
        # generater the seed region coordinates
        oligo_DB = self.seed_region_creation.apply(oligo_DB)

        # Some bowtie initializations, change the names
        index_exists = False
        index_name = os.path.basename(file_reference_DB)
        # Check if bowtie index exists
        for file in os.listdir(self.dir_seed_region):
            if re.search(f"^{index_name}.*", file):
                index_exists = True
                break
        # Create bowtie index if none exists
        if not index_exists:
            command1 = "bowtie-build " + file_reference_DB + " " + index_name
            process = Popen(command1, shell=True, cwd=self.dir_seed_region).wait()

        oligo_DB_seed = self._extract_seed_regions(oligo_DB)
        genes = list(oligo_DB_seed.keys())
        filtered_oligo_DBs = Parallel(n_jobs=n_jobs)(
            delayed(self._run_bowtie_seed_region)(
                oligo_DB_seed[gene], oligo_DB[gene], gene, index_name
            )
            for gene in genes
        )

        # reconstruct the oligos_DB
        for gene, filtered_oligo_DB in zip(genes, filtered_oligo_DBs):
            oligo_DB[gene] = filtered_oligo_DB

        return oligo_DB

    def _run_bowtie_seed_region(self, gene_DB_seed, gene_DB, gene, index_name):
        """Run Bowtie alignment tool to find regions of local similarity between sequences, where sequences are probes and transcripts.
        Bowtie identifies all allignments between the probes and transcripts and returns the number of mismatches and mismatch position for each alignment.

        :return: DataFrame with processed bowtie alignment search results.
        :rtype: pandas.DataFrame
        """

        file_probe_fasta_gene = self._create_fasta_file(
            gene_DB_seed, self.dir_fasta, gene
        )
        file_bowtie_gene = os.path.join(
            self.dir_seed_region,
            f"bowtie_{gene}.txt",
        )
        command = (
            "bowtie -x "
            + index_name
            + " -f -a -v "
            + str(self.num_mismatches)
            + " "
            + file_probe_fasta_gene
            + " "
            + file_bowtie_gene
        )

        process = Popen(command, shell=True, cwd=self.dir_seed_region).wait()

        # read the results of the bowtie search
        bowtie_results = self._read_bowtie_output(file_bowtie_gene)
        # filter the DB based on the bowtie results
        matching_probes = self._find_matching_probes(bowtie_results)
        filtered_gene_DB = self._filter_matching_probes(gene_DB, matching_probes)
        # remove the temporary files
        os.remove(os.path.join(self.dir_seed_region, file_bowtie_gene))
        os.remove(os.path.join(self.dir_fasta, file_probe_fasta_gene))
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

    def _find_matching_probes(self, bowtie_results):
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

    def _extract_seed_regions(self, oligo_DB):
        """geneate a new oligos DB containing only the seed regions of the probes."""
        oligo_DB_seed = {}
        for gene in oligo_DB.keys():
            oligo_DB_seed[gene] = {}
            for probe_id in oligo_DB[gene].keys():
                oligo_DB_seed[gene][probe_id] = {}
                start, end = (
                    oligo_DB[gene][probe_id]["seed_region_start"],
                    oligo_DB[gene][probe_id]["seed_region_end"],
                )
                seed_region_seq = oligo_DB[gene][probe_id]["probe_sequence"][
                    start : end + 1
                ]  # end must be included
                oligo_DB_seed[gene][probe_id]["probe_sequence"] = seed_region_seq
        return oligo_DB_seed
