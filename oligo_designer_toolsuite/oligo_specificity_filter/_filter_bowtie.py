############################################
# imports
############################################

import os
import re
import subprocess
from pathlib import Path

import pandas as pd
from joblib import Parallel, delayed

from oligo_designer_toolsuite.constants import (
    REGION_SEPARATOR_ANNOTATION,
    REGION_SEPARATOR_OLIGO,
)

from . import AlignmentSpecificityFilter

############################################
# Oligo Bowtie Filter Classes
############################################


class Bowtie(AlignmentSpecificityFilter):
    """This class filters oligos based on the Bowtie short read alignment tool.
    The user can customize the filtering by specifying the num_mismatches per oligo and mismatch_region, the region that should be considered for counting mismatches.
    That is, all oligos with number mismatches higher than num_mismatches inside the mismatch_region are filtered out.

    Use ``conda install -c bioconda bowtie`` to install Bowtie package

    :param dir_specificity: directory where alignment temporary files can be written
    :type dir_specificity: str
    :param num_mismatches: Threshold value on the number of mismatches required for each oligo. oligos where the number of mismatches greater than this threshold are considered valid. Possible values range from 0 to 3.
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
        """
        Apply the Bowtie filter in parallel on the given ``database``. Each job filters a single region, and up to ``n_jobs`` are run simultaneously.
        The filtered database is returned.

        :param database: Database containing the oligos and their features.
        :type database: dict
        :param file_reference: Path to the file that will be used as a reference for the alignment.
        :type file_reference: str
        :param n_jobs: Number of simultaneous parallel computations.
        :type n_jobs: int
        :return: Oligo info of user-specified regions.
        :rtype: dict
        """
        index_name = self._create_index(file_reference, n_jobs=n_jobs)

        regions = list(database.keys())
        filtered_database_regions = Parallel(n_jobs=n_jobs)(
            delayed(self._run_filter)(database, region, index_name)
            for region in regions
        )

        # reconstruct the oligos_DB
        for region, filtered_database_region in zip(regions, filtered_database_regions):
            database[region] = filtered_database_region

        return database

    def _create_index(self, file_reference: str, n_jobs: int):
        """
        Create a Bowtie database index.

        :param file_reference: Path to the reference file used for creating the index.
        :type file_reference: str
        :param n_jobs: Number of simultaneous parallel computations (currently unused)
        :type n_jobs: int
        :return: The name of the created database index.
        :rtype: str
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

        return index_name

    def _run_search(self, database, region, index_name, filter_same_region_matches):
        """Run Bowtie alignment tool to find regions of local similarity between sequences, where sequences are oligos and transcripts.
        Bowtie identifies all alignments between the oligos and transcripts and returns the number of mismatches and mismatch position for each alignment.


        :param database: Database containing the oligos.
        :type database: dict
        :param region: ID of the region processed.
        :type region: str
        :param index_name: Name of the Blastn database index.
        :type index_name: str
        :param filter_same_region_matches: Whether to filter out results within the same.
        :type filter_same_region_matches: bool
        :return: A tuple containing: an array of oligos with matches, and a dataframe containing alignment match data
        :rtype: (numpy.ndarray, pandas.DataFrame)
        """

        # TODO: This part has to change
        if region is not None:
            file_oligo_fasta_gene = self._create_fasta_file(
                database, self.dir_fasta, region
            )
        else:
            file_oligo_fasta_gene = self._create_fasta_multiple_regions(
                database, self.dir_fasta, regions=database.keys()
            )
            region = "multiple_regions"

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

        # remove the temporary files
        os.remove(os.path.join(self.dir_bowtie, file_bowtie_gene))
        os.remove(os.path.join(self.dir_fasta, file_oligo_fasta_gene))

        # filter the DB based on the bowtie results
        return self._find_matching_oligos(
            bowtie_results, filter_same_gene_matches=filter_same_region_matches
        )

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
        bowtie_results["query_gene_id"] = (
            bowtie_results["query"].str.split(REGION_SEPARATOR_OLIGO).str[0]
        )
        bowtie_results["reference_gene_id"] = (
            bowtie_results["reference"].str.split(REGION_SEPARATOR_ANNOTATION).str[0]
        )

        return bowtie_results

    def _find_matching_oligos(self, bowtie_results, filter_same_gene_matches=True):
        """Use the results of the Bowtie alignment search to identify oligos with high similarity (i.e. low number of mismatches) based on user defined thresholds.

        :param bowtie_results: DataFrame with processed blast alignment search results.
        :type bowtie_results: pandas.DataFrame:
        param filter_same_region_matches: Whether to filter out results within the same region (default is True)
        :type filter_same_region_matches: bool
        :return: A tuple containing: an array of oligos with matches, and a dataframe containing alignment match data
        :rtype: (numpy.ndarray, pandas.DataFrame)
        """

        if filter_same_gene_matches:
            bowtie_matches = bowtie_results[
                bowtie_results["query_gene_id"] != bowtie_results["reference_gene_id"]
            ]
        else:
            bowtie_matches = bowtie_results

        if self.strand == "plus":
            bowtie_matches = bowtie_matches[bowtie_matches["strand"] == "+"]
        elif self.strand == "minus":
            bowtie_matches = bowtie_matches[bowtie_matches["strand"] == "-"]

        oligos_with_match = bowtie_matches["query"].unique()

        return oligos_with_match, bowtie_matches

    def get_matching_oligo_pairs(
        self, database: dict, reference_fasta: str, n_jobs: int
    ):
        """
        Retrieve matching oligo pairs between a reference FASTA and a database. It returns a list of pairs, where each pair
        contains the name of the oligo from the database and its corresponding match from the reference.

        :param database: database containing the oligos.
        :type database: dict
        :param reference_fasta: path to the file that is used as an reference for the alignment
        :type reference_fasta: str

        :return: A list of matching oligo pairs.
        :rtype: list of tuple
        """
        database_name = self._create_index(reference_fasta, n_jobs=n_jobs)
        matches = self._run_search(
            database,
            region=None,
            index_name=database_name,
            filter_same_region_matches=False,
        )
        matches = matches[1]
        return list(zip(matches["query"].values, matches["reference"].values))
