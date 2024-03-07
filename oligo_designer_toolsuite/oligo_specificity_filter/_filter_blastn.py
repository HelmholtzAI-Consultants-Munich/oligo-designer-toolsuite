############################################
# imports
############################################

import os
import warnings
from abc import abstractmethod
from pathlib import Path
from typing import List, Union

import pandas as pd
from Bio.Blast.Applications import NcbiblastnCommandline, NcbimakeblastdbCommandline

from .._constants import _TYPES_SEQ
from ..database import OligoDatabase
from ..utils._checkers import check_if_list
from . import AlignmentSpecificityFilter

############################################
# Oligo Blast Filter Classes
############################################


class BlastNFilter(AlignmentSpecificityFilter):
    """A class that implements specificity filtering using BLASTN to align oligonucleotides against a reference database.
    This class allows for customization of BLAST search and hit parameters, and manages output directories for BLAST results.

    Blast (2.12 or higher)  can be installed via NCBI webpage (https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
    or via Bioconda (http://bioconda.github.io/recipes/blast/README.html) installation of Bowtie2 with:
    ``conda config --add channels conda-forge``
    ``conda config --add channels bioconda``
    ``conda update --all``
    ``conda install "blast>=2.12"``

    Useful BlastN search parameter:
    - word_size: Length of initial exact match. Default: 11
    - strand: Query strand(s) to search against database/subject. Choice of both, minus, or plus. Default: both
    - perc_identity: Percent identity cutoff. Default: 0
    All available BlastN search parameters are listed on the NCBI webpage (https://www.ncbi.nlm.nih.gov/books/NBK279684/).

    :param blast_search_parameters: Custom parameters for the BLAST search command.
    :type blast_search_parameters: dict
    :param blast_hit_parameters: Criteria to consider a BLAST hit significant for filtering.
    :type blast_hit_parameters: dict
    :param dir_output: Base directory for saving output files and BLAST databases. Defaults to "output".
    :type dir_output: str
    :param names_search_output: Column names for parsing BLAST search output.
    :type names_search_output: list
    """

    def __init__(
        self,
        blast_search_parameters: dict = {},
        blast_hit_parameters: dict = {},
        dir_output: str = "output",
        names_search_output: list = [
            "query",
            "reference",
            "alignment_length",
            "query_start",
            "query_end",
            "query_length",
        ],
    ):
        """Constructor for the BlastNFilter class."""
        super().__init__(dir_output)
        self.blast_search_parameters = blast_search_parameters
        self.blast_hit_parameters = blast_hit_parameters
        self.names_search_output = names_search_output

        # Define default output format for blast search filter. The fields are:
        # query, reference, alignment_length, query_start, query_end, query_length
        if "outfmt" not in self.blast_search_parameters.keys():
            self.blast_search_parameters[
                "outfmt"
            ] = "6 qseqid sseqid length qstart qend qlen"

        self.dir_blast = os.path.join(dir_output, "blast")
        Path(self.dir_blast).mkdir(parents=True, exist_ok=True)

    def _create_index(self, file_reference: str, n_jobs: int):
        """Creates a BLAST index for the reference database if it does not exist.

        :param file_reference: Path to the reference database file.
        :type file_reference: str
        :param n_jobs: The number of parallel jobs to run (currently not utilized in index creation but reserved for future use).
        :type n_jobs: int
        :return: The basename of the reference file used as the index name.
        :rtype: str
        """
        ## Create blast index
        filename_reference_index = os.path.basename(file_reference)

        # Check if blast database exists -> check for any of the blast index files, e.g. ".nhr" file
        if not os.path.exists(
            os.path.join(self.dir_blast, filename_reference_index + ".nhr")
        ):
            cmd = NcbimakeblastdbCommandline(
                input_file=file_reference,
                dbtype="nucl",
                out=os.path.join(self.dir_blast, filename_reference_index),
            )
            out, err = cmd()
        return filename_reference_index

    def _run_search(
        self,
        sequence_type: _TYPES_SEQ,
        oligo_database: OligoDatabase,
        filename_reference_index: str,
        region_ids: Union[str, List[str]] = None,
    ):
        """Executes a BLASTN search for an oligo database against a reference database index.

        :param sequence_type: The type of sequences being filtered, must be one of the predefined sequence types.
        :type sequence_type: _TYPES_SEQ
        :param oligo_database: The database of oligonucleotides to search.
        :type oligo_database: OligoDatabase
        :param filename_reference_index: The filename of the reference database index for BLASTN search.
        :type filename_reference_index: str
        :param region_ids: Specific region IDs within the oligo database to search. If None, searches all regions.
        :type region_ids: Union[str, List[str]], optional
        :return: A DataFrame containing the BLASTN search results.
        :rtype: pd.DataFrame
        """
        region_ids = check_if_list(region_ids)
        if region_ids:
            region_name = "_".join(region_ids)
        else:
            region_name = "all_regions"

        file_oligo_database = oligo_database.write_database_to_fasta(
            filename=f"oligo_database_blast_{region_name}",
            region_ids=region_ids,
            sequence_type=sequence_type,
        )
        file_blast_results = os.path.join(
            self.dir_blast, f"blast_results_{region_name}.txt"
        )

        cmd = NcbiblastnCommandline(
            query=file_oligo_database,
            out=file_blast_results,
            db=os.path.join(self.dir_blast, filename_reference_index),
            **self.blast_search_parameters,
        )
        out, err = cmd()

        # read the reuslts of the blast seatch
        blast_results = self._read_search_output(
            file_search_results=file_blast_results,
            names_search_output=self.names_search_output,
        )

        # remove temporary files
        os.remove(file_oligo_database)
        os.remove(file_blast_results)

        # return loaded results
        return blast_results

    def _find_hits(
        self,
        oligo_database: OligoDatabase,  # not utilized in this filter
        search_results: pd.DataFrame,
        consider_hits_from_input_region: bool,
    ):
        """Identifies significant hits from BLASTN results based on alignment length or coverage criteria.

        :param oligo_database: The oligonucleotide database.
        :type oligo_database: OligoDatabase
        :param search_results: The DataFrame containing results from a BLASTN search.
        :type search_results: pd.DataFrame
        :param consider_hits_from_input_region: Flag to indicate whether hits from the input region should be considered.
        :type consider_hits_from_input_region: bool
        :return: A tuple containing a DataFrame of significant BLAST hits and a list of oligos with these hits.
        :rtype: (pd.DataFrame, list)
        """
        if "min_alignment_length" in self.blast_hit_parameters.keys():
            if "coverage" in self.blast_hit_parameters.keys():
                warnings.warn(
                    "Both, 'min_alignment_length' and 'coverage' parameters were provided. Using 'min_alignment_length' parameter."
                )
            min_alignment_length = self.blast_hit_parameters["min_alignment_length"]
        elif "coverage" in self.blast_hit_parameters.keys():
            min_alignment_length = (
                search_results["query_length"]
                * self.blast_hit_parameters["coverage"]
                / 100
            )
        else:
            raise KeyError(
                "Please provide either 'coverage' or a 'min_alignment_length' in blast_hit_parameters!"
            )

        search_results["min_alignment_length"] = min_alignment_length

        if consider_hits_from_input_region:
            # remove all hits where query and reference come from the same region
            search_results = search_results[
                search_results["query_region_id"]
                != search_results["reference_region_id"]
            ]

        blast_table_hits = search_results.loc[
            search_results.alignment_length > search_results.min_alignment_length
        ]
        oligos_with_hits = blast_table_hits["query"].unique()

        return blast_table_hits, oligos_with_hits


############################################
# Oligo Blast Filter with Seedregion Classes
############################################


class BlastNSeedregionFilterBase(BlastNFilter):
    """Base class for BLASTN-based filters focusing on seed regions within oligonucleotides. It extends the BlastNFilter class,
    adding functionality to incorporate seed region information into the specificity filtering process. This class is designed
    to be subclassed with a concrete implementation of the method to add seed region information to BLAST results.

    :param blast_search_parameters: Custom parameters for the BLAST search command.
    :type blast_search_parameters: dict
    :param blast_hit_parameters: Criteria to consider a BLAST hit significant for filtering.
    :type blast_hit_parameters: dict
    :param dir_output: Directory for saving output files and BLAST databases.
    :type dir_output: str
    """

    def __init__(
        self,
        blast_search_parameters: dict = {},
        blast_hit_parameters: dict = {},
        dir_output: str = "output",
    ):
        """Constructor for the BlastNSeedregionFilterBase class."""
        super().__init__(blast_search_parameters, blast_hit_parameters, dir_output)

    @abstractmethod
    def _add_seed_region_information(
        self, oligo_database: OligoDatabase, blast_results: pd.DataFrame
    ):
        """Abstract method to add seed region information to BLAST results. This method must be implemented in subclasses to define
        how seed region data is incorporated into the filtering logic based on BLAST results.

        :param oligo_database: The database of oligonucleotides being analyzed.
        :type oligo_database: OligoDatabase
        :param blast_results: DataFrame containing results from a BLASTN search.
        :type blast_results: pd.DataFrame
        """

    def _find_hits(
        self,
        oligo_database: OligoDatabase,
        search_results: pd.DataFrame,
        consider_hits_from_input_region: bool,
    ):
        """Identifies oligonucleotides with significant hits based on BLAST results, taking into account seed region information. Adjusts hit
        criteria based on alignment length or coverage, and applies logic to exclude or include hits from the same input region.

        :param oligo_database: The oligonucleotide database.
        :type oligo_database: OligoDatabase
        :param search_results: DataFrame containing BLAST search results, enhanced with seed region information.
        :type search_results: pd.DataFram
        :param consider_hits_from_input_region: Flag to indicate whether hits from the input region should be considered.
        :type consider_hits_from_input_region: bool
        :return: A tuple of a DataFrame of significant BLAST hits and a list of oligo IDs with these hits.
        :rtype: (pd.DataFrame, list)
        """
        if "min_alignment_length" in self.blast_hit_parameters.keys():
            if "coverage" in self.blast_hit_parameters.keys():
                warnings.warn(
                    "Both, 'min_alignment_length' and 'coverage' parameters were provided. Using 'min_alignment_length' parameter."
                )
            min_alignment_length = self.blast_hit_parameters["min_alignment_length"]
        elif "coverage" in self.blast_hit_parameters.keys():
            min_alignment_length = (
                search_results["query_length"]
                * self.blast_hit_parameters["coverage"]
                / 100
            )
        else:
            raise KeyError(
                "Please provide either 'coverage' or a 'min_alignment_length' in blast_hit_parameters!"
            )

        search_results["min_alignment_length"] = min_alignment_length

        search_results = self._add_seed_region_information(
            oligo_database=oligo_database, search_results=search_results
        )

        if consider_hits_from_input_region:
            # remove all hits where query and reference come from the same region
            search_results = search_results[
                search_results["query_region_id"]
                != search_results["reference_region_id"]
            ]

        blast_table_hits = search_results.loc[
            (search_results.alignment_length > search_results.min_alignment_length)
            | (
                (search_results.query_start < search_results.seedregion_start)
                & (search_results.query_end > search_results.seedregion_end)
            )
        ]

        oligos_with_hits = blast_table_hits["query"].unique()

        return blast_table_hits, oligos_with_hits


class BlastNSeedregionFilter(BlastNSeedregionFilterBase):
    """Extends BlastNSeedregionFilterBase to implement specificity filtering with a focus on seed regions within oligonucleotides.
    This class allows for the integration of seed region considerations into the specificity filtering process against a reference database.

    :param seedregion_start: The start of the seed region, either as an absolute position (int) or as a percentage of the sequence length (float).
    :type seedregion_start: Union[int, float]
    :param seedregion_end: The end of the seed region, with the same type as seedregion_start.
    :type seedregion_end: Union[int, float]
    :param blast_search_parameters: Custom parameters for the BLAST search command.
    :type blast_search_parameters: dict
    :param blast_hit_parameters: Criteria to consider a BLAST hit significant for filtering.
    :type blast_hit_parameters: dict
    :param dir_output: Directory for saving output files and BLAST databases.
    :type dir_output: str
    """

    def __init__(
        self,
        seedregion_start: Union[int, float],
        seedregion_end: Union[int, float],
        blast_search_parameters: dict = {},
        blast_hit_parameters: dict = {},
        dir_output: str = "output",
    ):
        """Constructor for the BlastNSeedregionFilter class."""
        super().__init__(blast_search_parameters, blast_hit_parameters, dir_output)

        self.seedregion_start = seedregion_start
        self.seedregion_end = seedregion_end

    def _add_seed_region_information(
        self, oligo_database: OligoDatabase, search_results: pd.DataFrame
    ):
        """Adds seed region information to BLAST results for further filtering.

        :param oligo_database: The oligonucleotide database being analyzed.
        :type oligo_database: OligoDatabase
        :param search_results: DataFrame containing results from a BLASTN search.
        :type search_results: pd.DataFrame
        :return: Updated BLAST results with seed region information.
        :rtype: pd.DataFrame
        """
        oligo_database.calculate_seedregion(
            start=self.seedregion_start, end=self.seedregion_end
        )

        seedregion = pd.merge(
            left=oligo_database.get_oligo_attribute("seedregion_start"),
            right=oligo_database.get_oligo_attribute("seedregion_end"),
            on="oligo_id",
        )
        search_results = pd.merge(
            left=search_results,
            right=seedregion,
            left_on="query",
            right_on="oligo_id",
            how="left",
        )
        return search_results


class BlastNSeedregionLigationsiteFilter(BlastNSeedregionFilterBase):
    """Extends BlastNSeedregionFilterBase to implement specificity filtering focused on the ligation site's seed regions within oligonucleotides.
    This class considers the seed region around the ligation site for BLASTN-based specificity filtering against a reference database.

    :param seedregion_size: The size of the seed region around the ligation site to be considered.
    :type seedregion_size: int
    :param blast_search_parameters: Custom parameters for the BLAST search command.
    :type blast_search_parameters: dict
    :param blast_hit_parameters: Criteria to consider a BLAST hit significant for filtering.
    :type blast_hit_parameters: dict
    :param dir_output: Directory for saving output files and BLAST databases.
    :type dir_output: str
    """

    def __init__(
        self,
        seedregion_size: int,
        blast_search_parameters: dict = {},
        blast_hit_parameters: dict = {},
        dir_output: str = "output",
    ):
        """Constructor for the BlastNSeedregionLigationsiteFilter class."""
        super().__init__(blast_search_parameters, blast_hit_parameters, dir_output)
        self.seedregion_size = seedregion_size

    def _add_seed_region_information(
        self, oligo_database: OligoDatabase, search_results: pd.DataFrame
    ):
        """Adds seed region information around the ligation site to BLAST results for further filtering.

        :param oligo_database: The oligonucleotide database being analyzed.
        :type oligo_database: OligoDatabase
        :param search_results: DataFrame containing results from a BLASTN search.
        :type search_results: pd.DataFrame
        :return: Updated BLAST results with seed region information around the ligation site.
        :rtype: pd.DataFrame
        """
        oligo_database.calculate_seedregion_ligationsite(
            seedregion_size=self.seedregion_size
        )

        seedregion = pd.merge(
            left=oligo_database.get_oligo_attribute("seedregion_start"),
            right=oligo_database.get_oligo_attribute("seedregion_end"),
            on="oligo_id",
        )
        search_results = pd.merge(
            left=search_results,
            right=seedregion,
            left_on="query",
            right_on="oligo_id",
            how="left",
        )
        return search_results
