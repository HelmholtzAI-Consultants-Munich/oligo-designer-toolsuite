############################################
# imports
############################################

import os
import re
from pathlib import Path

import pandas as pd
from Bio.Blast.Applications import NcbiblastnCommandline, NcbimakeblastdbCommandline
from joblib import Parallel, delayed

from oligo_designer_toolsuite.constants import (
    REGION_SEPARATOR_ANNOTATION,
    REGION_SEPARATOR_OLIGO,
)

from . import AlignmentSpecificityFilter

############################################
# Oligo Blast Filter Classes
############################################


class Blastn(AlignmentSpecificityFilter):
    """This class filters oligos based on the blast alignment tool. All the oligos which have a match with a percentage identity higher
    than the one given in input are filtered out.

    :param dir_specificity: directory where alignment temporary files can be written
    :type dir_specificity: str
    :param word_size: word size for the blastn seed (exact match to target)
    :type word_size: int
    :param percent_identity: maximum similarity between oligos and target sequences, ranging from 0 to 100% (no missmatch)
    :type percent_identity: int
    :param strand: strand of the query sequence to search
    :type strand: str
    :param coverage: minimum coverage between oligos and target sequence, ranging from 0 to 100% (full coverage)
    :type coverage: float
    :param min_alignment_length: minimum alignment length between oligos and target sequence
    :type min_alignment_length: int
    """

    def __init__(
        self,
        dir_specificity: str,
        coverage: float = None,
        min_alignment_length: int = None,
        percent_identity: float = None,
        blast_search_parameters: dict = {}
    ):
        """Constructor."""
        super().__init__(dir_specificity)
        if coverage is None and min_alignment_length is None:
            raise Exception(
                "Please provide either coverage or a minimum alignment length"
            )
        self.coverage = coverage
        self.min_alignment_length = min_alignment_length
        self.perc_identity = percent_identity
        self.blast_search_parameters = blast_search_parameters

        self.dir_blast = os.path.join(self.dir_specificity, "blast")
        Path(self.dir_blast).mkdir(parents=True, exist_ok=True)

    def apply(self, database: dict, file_reference: str, n_jobs: int):
        """
        Apply the blastn filter in parallel on the given ``database``. Each job filters a
        single region, and at the same time, at most ``n_jobs`` jobs are generated.
        The filtered database is returned.

        :param database: Database containing the oligos and their features.
        :type database: dict
        :param file_reference: Path to the file that will be used as reference for the alignment.
        :type file_reference: str
        :param n_jobs: Number of simultaneous parallel computations.
        :type n_jobs: int
        :return: Oligo info of user-specified regions.
        :rtype: dict
        """

        # run the balst search
        database_name = self._create_index(file_reference, n_jobs=n_jobs)
        kwargs = {
            "outfmt": "6 qseqid sseqid length qstart qend qlen",
            "num_threads": 1,  # ????
        }
        regions = list(database.keys())
        filtered_database_regions = Parallel(n_jobs=n_jobs)(
            delayed(self._run_filter)(database, region, database_name, **kwargs)
            for region in regions
        )

        # reconstruct the oligos db and return it
        for region, filtered_database_region in zip(regions, filtered_database_regions):
            database[region] = filtered_database_region

        return database

    def _compute_percent_identity_from_alignment_length(alignment_length):
        if alignment_length > 50: percent_identity = 50
        elif alignment_length 


    def _create_index(self, file_reference: str, n_jobs: int):
        """
        Create a Blastn database index.

        :param file_reference: Path to the reference file used for creating the index.
        :type file_reference: str
        :param n_jobs: Number of simultaneous parallel computations (currently unused)
        :type n_jobs: int
        :return: The name of the created database index.
        :rtype: str
        """
        # create blast database
        database_exists = False
        database_name = os.path.basename(file_reference)
        # Check if blast database exists
        for file in os.listdir(self.dir_blast):
            if re.search(f"^{database_name}.*", file):
                database_exists = True
                break

        if not database_exists:
            cmd = NcbimakeblastdbCommandline(
                input_file=file_reference,
                dbtype="nucl",
                out=os.path.join(self.dir_blast, database_name),
            )
            out, err = cmd()
        return database_name

    def _run_search(
        self, database, region, index_name, filter_same_region_matches, **kwargs
    ):
        """
        Run the BlastN alignment tool to identify regions of local similarity between sequences.
        These sequences include oligos from the database and those in the index (e.g., transcript, genome, etc.).
        BlastN determines the transcript regions where oligos have matching sequences.

        :param database: Database containing the oligos.
        :type database: dict
        :param region: ID of the region being processed.
        :type region: str
        :param index_name: Name of the BlastN database index.
        :type index_name: str
        :param filter_same_region_matches: Whether to filter out matches within the same region.
        :type filter_same_region_matches: bool
        :param kwargs: Additional keyword arguments for the BlastN command.
        :type kwargs: dict
        :return: A tuple containing: an array of oligos with matches, and a dataframe containing alignment match data.
        :rtype: (numpy.ndarray, pandas.DataFrame)
        """

        # TODO: This part has to change
        if region is not None:
            file_oligo_fasta_region = self._create_fasta_file(
                database, self.dir_fasta, region
            )
        else:
            file_oligo_fasta_region = self._create_fasta_multiple_regions(
                database, self.dir_fasta, regions=database.keys()
            )
            region = "multiple_regions"

        file_blast_region = os.path.join(self.dir_blast, f"blast_{region}.txt")
        cmd = NcbiblastnCommandline(
            query=file_oligo_fasta_region,
            out=file_blast_region,
            db=os.path.join(self.dir_blast, index_name),
            word_size=self.word_size,
            perc_identity=self.perc_identity,
            strand=self.strand,
            **kwargs,
        )
        out, err = cmd()

        # read the reuslts of the blast seatch
        blast_results = self._read_blast_output(file_blast_region)

        # remove temporary files
        os.remove(file_blast_region)
        os.remove(file_oligo_fasta_region)
        # filter the DB based on the blast results
        return self._find_matching_oligos(
            blast_results, filter_same_region_matches=filter_same_region_matches
        )

    def _read_blast_output(self, file_blast_region):
        """
        Load the output of the BlastN alignment search into a DataFrame and process the results.

        :param file_blast_region: Path to the BlastN output file.
        :type file_blast_region: str
        :return: Processed BlastN results as a DataFrame.
        :rtype: pandas.DataFrame
        """

        blast_results = pd.read_csv(
            file_blast_region,
            header=None,
            sep="\t",
            low_memory=False,
            names=[
                "query",
                "reference",
                "alignment_length",
                "query_start",
                "query_end",
                "query_length",
            ],
            engine="c",
            dtype={
                "query": str,
                "reference": str,
                "alignment_length": int,
                "query_start": int,
                "query_end": int,
                "query_length": int,
            },
        )
        # return the real matches, that is the ones not belonging to the same region of the query oligo
        blast_results["query_region_id"] = (
            blast_results["query"].str.split(REGION_SEPARATOR_OLIGO).str[0]
        )
        blast_results["reference_region_id"] = (
            blast_results["reference"].str.split(REGION_SEPARATOR_ANNOTATION).str[0]
        )
        return blast_results

    def _find_matching_oligos(self, blast_results, filter_same_region_matches=True):
        """Use the results of the BlastN alignement search to remove oligos.

        :param blast_results: DataFrame with processed blast alignment search results.
        :type blast_results: pandas.DataFrame:
        param filter_same_region_matches: Whether to filter out results within the same region (default is True)
        :type filter_same_region_matches: bool
        :return: A tuple containing: an array of oligos with matches, and a dataframe containing alignment match data
        :rtype: (numpy.ndarray, pandas.DataFrame)
        """

        if filter_same_region_matches:
            blast_matches = blast_results[
                blast_results["query_region_id"] != blast_results["reference_region_id"]
            ]
        else:
            blast_matches = blast_results
        if self.coverage is not None:
            values = blast_matches["query_length"] * self.coverage / 100
            blast_matches.insert(
                len(blast_matches.columns), "min_alignment_length", values
            )
        else:
            blast_matches["min_alignment_length"] = self.min_alignment_length

        blast_matches_filtered = blast_matches.loc[
            blast_matches.alignment_length > blast_matches.min_alignment_length
        ]

        oligos_with_match = blast_matches_filtered["query"].unique()

        return oligos_with_match, blast_matches_filtered

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
            outfmt="6 qseqid sseqid length qstart qend qlen",
            num_threads=1,
        )
        matches = matches[1]
        return list(zip(matches["query"].values, matches["reference"].values))
