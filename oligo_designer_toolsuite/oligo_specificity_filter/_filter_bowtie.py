############################################
# imports
############################################

import os
import re
import subprocess
import pandas as pd

from pathlib import Path
from joblib import Parallel, delayed

from . import SpecificityFilterBase

############################################
# Oligo Bowtie Filter Classes
############################################

class Bowtie(SpecificityFilterBase):
    """This class filters oligos based on the Bowtie short read alignment tool.
    The user can customize the filtering by specifying the num_mismatches per oligo and mismatch_region, the region that should be considered for counting mismatches.
    That is, all oligos with number mismatches higher than num_mismatches inside the mismatch_region are filtered out.

    Use ``conda install -c bioconda bowtie`` to install Bowtie package

    :param dir_specificity: directory where alignement temporary files can be written
    :type dir_specificity: str
    :param num_mismatches: Threshold value on the number of mismatches required for each oligo. ligos where the number of mismatches greater than this threshhold are considered valid. Possible values range from 0 to 3.
    :type num_mismatches: int
    :param mismatch_region: The region of the oligo where the mismatches are considered. Oligos that have less than or equal to num_mismatches in the first L bases (where L is 5 or greater) are filtered out. If ``None`` then the whole sequence is considered, defaults to None
    :type mismatch_region: int
    """

    def __init__(
        self,
        dir_specificity: str,
        num_mismatches: int = 3,
        mismatch_region: int = None,
        strand: str = None,
    ):
        super().__init__(dir_specificity)

        if num_mismatches > 3:
            raise ValueError(
                "Choice of min_mismatches out of range for bowtie allignment tool. Please choose a value no greater than 3"
            )
        else:
            self.num_mismatches = num_mismatches

        if mismatch_region is not None and mismatch_region < 5:
            raise ValueError(
                "Choice of mismatch_region out of range for bowtie allignment tool. Please choose a value no less than 5"
            )
        else:
            self.mismatch_region = mismatch_region

        self.strand = strand

        self.dir_bowtie = os.path.join(self.dir_specificity, "bowtie")
        Path(self.dir_bowtie).mkdir(parents=True, exist_ok=True)

        self.dir_fasta = os.path.join(self.dir_specificity, "fasta")
        Path(self.dir_fasta).mkdir(parents=True, exist_ok=True)

    def apply(self, database: dict, file_reference: str, n_jobs: int):
        """Apply the bowtie filter in parallel on the given ``database``. Each jobs filters a single region, and  at the same time are generated at most ``n_job`` jobs.
        The filtered database is returned.

        :param database: database containing the oligos and their features
        :type database: dict
        :param file_reference: path to the file that will be used as reference for the alignement
        :type file_reference: str
        :param n_jobs: number of simultaneous parallel computations
        :type n_jobs: int
        :return: oligo info of user-specified regions
        :rtype: dict
        """
        # Some bowtie initializations, change the names
        index_exists = False
        index_name = os.path.basename(file_reference)
        # Check if bowtie index exists
        for file in os.listdir(self.dir_bowtie):
            if re.search(f"^{index_name}.*", file):
                index_exists = True
                break

        # Create bowtie index if none exists
        if not index_exists:
            command1 = (
                "bowtie-build --quiet --threads "
                + str(n_jobs)
                + " -f "
                + file_reference
                + " "
                + index_name
            )
            process = subprocess.Popen(command1, shell=True, cwd=self.dir_bowtie).wait()

        regions = list(database.keys())
        filtered_database_regions = Parallel(n_jobs=n_jobs)(
            delayed(self._run_bowtie)(database[region], region, index_name)
            for region in regions
        )

        # reconstruct the oligos_DB
        for region, filtered_database_region in zip(regions, filtered_database_regions):
            database[region] = filtered_database_region

        return database

    def _run_bowtie(self, databse_region, region, index_name):
        """Run Bowtie alignment tool to find regions of local similarity between sequences, where sequences are oligos and transcripts.
        Bowtie identifies all alignments between the oligos and transcripts and returns the number of mismatches and mismatch position for each alignment.

        :param databse_region: database containing the oligos form one region
        :type databse_region: dict
        :param region: id of the region processed
        :type region: str
        """

        file_oligo_fasta_gene = self._create_fasta_file(
            databse_region, self.dir_fasta, region
        )
        file_bowtie_gene = os.path.join(
            self.dir_bowtie,
            f"bowtie_{region}.txt",
        )
        if self.mismatch_region is not None:
            command = (
                "bowtie --quiet -x "
                + index_name
                + " -f -a -n "
                + str(self.num_mismatches)
                + " -l "
                + str(self.mismatch_region)
                + " "
                + file_oligo_fasta_gene
                + " "
                + file_bowtie_gene
            )
        else:
            command = (
                "bowtie --quiet -x "
                + index_name
                + " -f -a -v "
                + str(self.num_mismatches)
                + " "
                + file_oligo_fasta_gene
                + " "
                + file_bowtie_gene
            )

        process = subprocess.Popen(
            command,
            shell=True,
            cwd=self.dir_bowtie,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        ).wait()

        # read the results of the bowtie search
        bowtie_results = self._read_bowtie_output(file_bowtie_gene)
        # filter the DB based on the bowtie results
        matching_oligos = self._find_matching_oligos(bowtie_results)
        filtered_database_region = self._filter_matching_oligos(
            databse_region, matching_oligos
        )
        # remove the temporary files
        os.remove(os.path.join(self.dir_bowtie, file_bowtie_gene))
        os.remove(os.path.join(self.dir_fasta, file_oligo_fasta_gene))
        return filtered_database_region

    def _read_bowtie_output(self, file_bowtie_gene):
        """Load the output of the bowtie alignment search into a DataFrame and process the results."""
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
        # return the real matches, that is the ones not belonging to the same region of the query oligo
        bowtie_results["query_gene_id"] = bowtie_results["query"].str.split("_").str[0]
        bowtie_results["reference_gene_id"] = (
            bowtie_results["reference"].str.split("::").str[0]
        )

        return bowtie_results

    def _find_matching_oligos(self, bowtie_results):
        """Use the results of the Bowtie alignment search to identify oligos with high similarity (i.e. low number of mismatches) based on user defined thresholds.

        :param bowtie_results: DataFrame with processed bowtie alignment search results.
        :type bowtie_results: pandas.DataFrame
        """

        bowtie_matches = bowtie_results[
            bowtie_results["query_gene_id"] != bowtie_results["reference_gene_id"]
        ]

        if self.strand == "plus":
            bowtie_matches = bowtie_matches[bowtie_matches["strand"] == "+"]
        elif self.strand == "minus":
            bowtie_matches = bowtie_matches[bowtie_matches["strand"] == "-"]

        oligos_with_match = bowtie_matches["query"].unique()
        return oligos_with_match
