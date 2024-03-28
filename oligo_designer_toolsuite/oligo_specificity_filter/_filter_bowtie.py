############################################
# imports
############################################

import os
import subprocess
from pathlib import Path
from typing import List, Union

import pandas as pd

from .._constants import _TYPES_SEQ
from ..database import OligoDatabase
from ..utils._checkers import check_if_list
from . import AlignmentSpecificityFilter

############################################
# Oligo Bowtie Filter Classes
############################################


class BowtieFilter(AlignmentSpecificityFilter):
    """A class that implements specificity filtering using the Bowtie alignment tool to align oligonucleotides against a reference database.
    This class manages Bowtie search parameters and parses the output for specificity analysis.

    Bowtie (1.3 or higher) can be installed via Bowtie webpage (https://bowtie-bio.sourceforge.net/manual.shtml#obtaining-bowtie)
    or via Bioconda (http://bioconda.github.io/recipes/bowtie/README.html) installation of Bowtie with:
    ``conda config --add channels conda-forge``
    ``conda config --add channels bioconda``
    ``conda update --all``
    ``conda install "bowtie>=1.3.1"``

    Useful bowtie search parameter:
    --norc: only report hits on the plus strand
    --nofw: only report hits on the minus strand
    -v <int>: Report alignments with at most <int> mismatches.
    -n <int>: Maximum number of mismatches permitted in the “seed”, i.e. the first L base pairs of the read (where L is set with -l/--seedlen). This may be 0, 1, 2 or 3 and the default is 2.
    -l <int>: The “seed length”; i.e., the number of bases on the high-quality end of the read to which the -n ceiling applies. The lowest permitted setting is 5 and the default is 28. bowtie is faster for larger values of -l.
    All available Bowtie search parameters are listed on the Bowtie webpage (https://bowtie-bio.sourceforge.net/manual.shtml#obtaining-bowtie).

    :param bowtie_search_parameters: Custom parameters for the Bowtie search command.
    :type bowtie_search_parameters: dict
    :param dir_output: Base directory for saving output files and Bowtie databases. Defaults to "output".
    :type dir_output: str
    :param names_search_output: Column names for parsing Bowtie search output.
    :type names_search_output: list
    """

    def __init__(
        self,
        bowtie_search_parameters: dict = {},
        # bowtie_hit_parameters: dict = {},
        dir_output: str = "output",
        names_search_output: list = [
            "query",
            "strand",
            "reference",
            "reference_start",
            "query_sequence",
            "read_quality",
            "num_instances",
            "mismatch_positions",
        ],
    ):
        """Constructor for the BowtieFilter class."""
        super().__init__(dir_output)

        self.bowtie_search_parameters = bowtie_search_parameters
        # self.bowtie_hit_parameters = bowtie_hit_parameters
        self.names_search_output = names_search_output

        self.dir_bowtie = os.path.join(self.dir_output, "bowtie")
        Path(self.dir_bowtie).mkdir(parents=True, exist_ok=True)

    def _create_index(self, file_reference: str, n_jobs: int):
        """Creates a Bowtie index for the reference database. The index facilitates
        fast alignment searches against the reference.

        :param file_reference: Path to the reference database file.
        :type file_reference: str
        :param n_jobs: The number of parallel jobs to run.
        :type n_jobs: int
        :return: The basename of the reference file used as the index name.
        :rtype: str
        """
        ## Create bowtie index
        file_reference = os.path.abspath(file_reference)
        filename_reference_index = os.path.basename(file_reference)

        cmd = (
            "bowtie-build --quiet --offrate 4"
            + " --threads "
            + str(n_jobs)
            + " -f "
            + file_reference
            + " "
            + filename_reference_index
        )
        process = subprocess.Popen(cmd, shell=True, cwd=self.dir_bowtie).wait()

        return filename_reference_index

    def _run_search(
        self,
        sequence_type: _TYPES_SEQ,
        oligo_database: OligoDatabase,
        file_index: str,
        region_ids: Union[str, List[str]] = None,
    ):
        """Performs a Bowtie search of oligonucleotide sequences against a reference database using a previously
        created Bowtie index.

        :param sequence_type: The type of sequences being filtered, must be one of the predefined sequence types.
        :type sequence_type: _TYPES_SEQ
        :param oligo_database: The database of oligonucleotides to search.
        :type oligo_database: OligoDatabase
        :param file_index: The filename of the reference database index for Bowtie search.
        :type file_index: str
        :param region_ids: Specific region IDs within the oligo database to search. If None, searches all regions.
        :type region_ids: Union[str, List[str]], optional
        :return: A DataFrame containing the Bowtie search results.
        :rtype: pd.DataFrame
        """
        region_ids = check_if_list(region_ids)
        if region_ids:
            region_name = "_".join(region_ids)
        else:
            region_name = "all_regions"

        file_oligo_database = oligo_database.write_database_to_fasta(
            filename=f"oligo_database_bowtie_{region_name}",
            region_ids=region_ids,
            sequence_type=sequence_type,
        )
        file_bowtie_results = os.path.join(self.dir_bowtie, f"bowtie_results_{region_name}.txt")

        cmd_parameters = ""
        for parameter, value in self.bowtie_search_parameters.items():
            cmd_parameters += f" {parameter} {value}"

        cmd = (
            "bowtie --quiet"
            + " -x "
            + file_index
            + " -f"  # fast file is input
            + " -a"  # report all alignments -> TODO: does this make sense or set e.g. -k 100
            + cmd_parameters
            + " "
            + file_oligo_database
            + " "
            + file_bowtie_results
        )
        process = subprocess.Popen(cmd, shell=True, cwd=self.dir_bowtie).wait()

        # read the reuslts of the bowtie search
        bowtie_results = self._read_search_output(
            file_search_results=file_bowtie_results,
            names_search_output=self.names_search_output,
        )

        # remove temporary files
        os.remove(file_oligo_database)
        os.remove(file_bowtie_results)

        # return loaded results
        return bowtie_results

    def _find_hits(
        self,
        oligo_database: OligoDatabase,  # not used in this filter
        search_results: pd.DataFrame,
        consider_hits_from_input_region: bool,
    ):
        """Identifies oligonucleotides with significant alignment hits from Bowtie search results, optionally
        excluding hits from the same input region.

        :param oligo_database: The oligonucleotide database.
        :type oligo_database: OligoDatabase
        :param search_results: The DataFrame containing results from a Bowtie search.
        :type search_results: pd.DataFrame
        :param consider_hits_from_input_region: Flag to indicate whether hits from the input region should be considered.
        :type consider_hits_from_input_region: bool
        :return: A tuple containing a DataFrame of significant Bowtie hits and a list of oligos with these hits.
        :rtype: (pd.DataFrame, list)
        """
        if consider_hits_from_input_region:
            # remove all hits where query and reference come from the same region
            search_results = search_results[
                search_results["query_region_id"] != search_results["reference_region_id"]
            ]

        oligos_with_hits = search_results["query"].unique()

        return search_results, oligos_with_hits


############################################
# Oligo Bowtie2 Filter Classes
############################################


class Bowtie2Filter(AlignmentSpecificityFilter):
    """A class that implements specificity filtering using the Bowtie2 alignment tool to align oligonucleotides against a reference database.
    This class manages Bowtie2 search parameters and parses the output for specificity analysis.

    Bowtie2 (2.5 or higher) can be installed via Bowtie2 webpage (https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#obtaining-bowtie-2)
    or via Bioconda (http://bioconda.github.io/recipes/bowtie2/README.html) installation of Bowtie2 with:
    ``conda config --add channels conda-forge``
    ``conda config --add channels bioconda``
    ``conda update --all``
    ``conda install "bowtie2>=2.5"``

    Useful Bowtie2 search parameter:
    --norc: only report hits on the plus strand
    --nofw: only report hits on the minus strand
    -N <int>: Sets the number of mismatches to allowed in a seed alignment during multiseed alignment. Can be set to 0 or 1. Setting this higher makes alignment slower (often much slower) but increases sensitivity. Default: 0.
    -L <int>: Sets the length of the seed substrings to align during multiseed alignment. Smaller values make alignment slower but more sensitive. Default: the --sensitive preset is used by default, which sets -L to 22 and 20 in --end-to-end mode and in --local mode.
    All available Bowtie2 search parameters are listed on the Bowtie2 webpage (https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#obtaining-bowtie-2).

    :param bowtie_search_parameters: Custom parameters for the Bowtie2 search command.
    :type bowtie_search_parameters: dict
    :param dir_output: Base directory for saving output files and Bowtie2 databases. Defaults to "output".
    :type dir_output: str
    :param names_search_output: Column names for parsing Bowtie2 search output.
    :type names_search_output: list
    """

    def __init__(
        self,
        bowtie_search_parameters: dict = {},
        # bowtie_hit_parameters: dict = {},
        dir_output: str = "output",
        names_search_output: list = [
            "query",
            "flags",
            "reference",
            "reference_start",
            "mapping_quality",
            "CIGAR_alignment",
            "mate_sequence_name",
            "mate_sequence_offset",
            "mate_sequence_fragment_length",
            "sequence",
            "read_qualities",
        ],
    ):
        """Constructor for the Bowtie2Filter class."""
        super().__init__(dir_output)

        self.bowtie_search_parameters = bowtie_search_parameters
        # self.bowtie_hit_parameters = bowtie_hit_parameters
        self.names_search_output = names_search_output

        self.dir_bowtie = os.path.join(self.dir_output, "bowtie2")
        Path(self.dir_bowtie).mkdir(parents=True, exist_ok=True)

    def _create_index(self, file_reference: str, n_jobs: int):
        """Creates a Bowtie2 index for the reference database. The index facilitates
        fast alignment searches against the reference.

        :param file_reference: Path to the reference database file.
        :type file_reference: str
        :param n_jobs: The number of parallel jobs to run.
        :type n_jobs: int
        :return: The basename of the reference file used as the index name.
        :rtype: str
        """
        ## Create bowtie index
        file_reference = os.path.abspath(file_reference)
        filename_reference_index = os.path.basename(file_reference)

        # Check if bowtie database exists -> check for any of the bowtie index files, e.g. ".1.bt2" file
        cmd = (
            "bowtie2-build --quiet --offrate 4"
            + " --threads "
            + str(n_jobs)
            + " -f "
            + file_reference
            + " "
            + filename_reference_index
        )
        process = subprocess.Popen(cmd, shell=True, cwd=self.dir_bowtie).wait()

        return filename_reference_index

    def _run_search(
        self,
        sequence_type: _TYPES_SEQ,
        oligo_database: OligoDatabase,
        file_index: str,
        region_ids: Union[str, List[str]] = None,
    ):
        """Performs a Bowtie2 search of oligonucleotide sequences against a reference database using a previously
        created Bowtie2 index.

        :param sequence_type: The type of sequences being filtered, must be one of the predefined sequence types.
        :type sequence_type: _TYPES_SEQ
        :param oligo_database: The database of oligonucleotides to search.
        :type oligo_database: OligoDatabase
        :param file_index: The filename of the reference database index for Bowtie2 search.
        :type file_index: str
        :param region_ids: Specific region IDs within the oligo database to search. If None, searches all regions.
        :type region_ids: Union[str, List[str]], optional
        :return: A DataFrame containing the Bowtie2 search results.
        :rtype: pd.DataFrame
        """
        region_ids = check_if_list(obj=region_ids)
        if region_ids:
            region_name = "_".join(region_ids)
        else:
            region_name = "all_regions"

        file_oligo_database = oligo_database.write_database_to_fasta(
            filename=f"oligo_database_bowtie2_{region_name}",
            region_ids=region_ids,
            sequence_type=sequence_type,
        )
        file_bowtie_results = os.path.join(self.dir_bowtie, f"bowtie2_results_{region_name}.txt")

        cmd_parameters = ""
        for parameter, value in self.bowtie_search_parameters.items():
            cmd_parameters += f" {parameter} {value}"

        cmd = (
            "bowtie2 --quiet"
            + " --no-hd --no-unal"
            + " -x "
            + file_index
            + " -f"  # fast file is input
            + " -a"  # report all alignments -> TODO: does this make sense or set e.g. -k 100
            + cmd_parameters
            + " -U "
            + file_oligo_database
            + " -S "
            + file_bowtie_results
        )
        process = subprocess.Popen(cmd, shell=True, cwd=self.dir_bowtie).wait()

        # read the reuslts of the bowtie seatch
        bowtie_results = self._read_search_output(
            file_search_results=file_bowtie_results,
            names_search_output=self.names_search_output,
            usecols=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
        )

        # remove temporary files
        os.remove(file_oligo_database)
        os.remove(file_bowtie_results)

        # return loaded results
        return bowtie_results

    def _find_hits(
        self,
        oligo_database: OligoDatabase,  # not used in this filter
        search_results: pd.DataFrame,
        consider_hits_from_input_region: bool,
    ):
        """Identifies oligonucleotides with significant alignment hits from Bowtie2 search results, optionally
        excluding hits from the same input region.

        :param oligo_database: The oligonucleotide database.
        :type oligo_database: OligoDatabase
        :param search_results: The DataFrame containing results from a Bowtie2 search.
        :type search_results: pd.DataFrame
        :param consider_hits_from_input_region: Flag to indicate whether hits from the input region should be considered.
        :type consider_hits_from_input_region: bool
        :return: A tuple containing a DataFrame of significant Bowtie2 hits and a list of oligos with these hits.
        :rtype: (pd.DataFrame, list)
        """
        if consider_hits_from_input_region:
            # remove all hits where query and reference come from the same region
            search_results = search_results[
                search_results["query_region_id"] != search_results["reference_region_id"]
            ]

        oligos_with_hits = search_results["query"].unique()

        return search_results, oligos_with_hits
