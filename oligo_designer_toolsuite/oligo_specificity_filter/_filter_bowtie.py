############################################
# imports
############################################

import os
import subprocess
from typing import List, Union

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq

from oligo_designer_toolsuite._constants import _TYPES_SEQ
from oligo_designer_toolsuite.database import OligoDatabase
from oligo_designer_toolsuite.oligo_specificity_filter import AlignmentSpecificityFilter

from ..utils._checkers import check_if_list
from ..utils._sequence_processor import get_sequence_from_annotation

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

    The hits returned by Bowtie can be further filtered using machine learning models. For more information regarding which filters are available
    refer to https://github.com/HelmholtzAI-Consultants-Munich/oligo-designer-toolsuite-AI-filters.

    :param search_parameters: Custom parameters for the Bowtie search command.
    :type search_parameters: dict
    :param hit_parameters: Criteria to consider a Bowtie hit significant for filtering.
    :type hit_parameters: dict
    :param names_search_output: Column names for parsing Bowtie search output.
    :type names_search_output: list
    :param filter_name: Subdirectory path for the output, i.e. <dir_output>/<filter_name>, defaults to "bowtie_filter".
    :type filter_name: str, optional
    :param dir_output: Directory for saving intermediate files, defaults to "output"
    :type dir_output: str, optional
    """

    def __init__(
        self,
        search_parameters: dict = {},
        hit_parameters: dict = {},
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
        filter_name: str = "bowtie_filter",
        dir_output: str = "output",
    ):
        """Constructor for the BowtieFilter class."""
        super().__init__(filter_name, dir_output)

        self.search_parameters = search_parameters
        self.hit_parameters = hit_parameters  # currently not used
        self.names_search_output = names_search_output

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
        process = subprocess.Popen(cmd, shell=True, cwd=self.dir_output).wait()

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
        file_bowtie_results = os.path.join(self.dir_output, f"bowtie_results_{region_name}.txt")

        cmd_parameters = ""
        for parameter, value in self.search_parameters.items():
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
        process = subprocess.Popen(cmd, shell=True, cwd=self.dir_output).wait()

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
        :rtype: pd.DataFrame
        """
        if not consider_hits_from_input_region:
            # remove all hits where query and reference come from the same region
            search_results = search_results[
                search_results["query_region_id"] != search_results["reference_region_id"]
            ]

        return search_results

    def _get_references(self, table_hits: pd.DataFrame, file_reference: str, region_id: str):
        """
        Retrieve the references sequences from the search results.

        :param table_hits: DataFrame with the oligos that match the blast search.
        :type table_hits: pd.DataFrame
        :param file_reference: Path to the fasta file used as reference for the search.
        :type file_reference: str
        :param region_id: The identifier for the region within the database to filter.
        :type region_id: str
        :return: Reference sequences
        :rtype: list
        """
        required_fields = [
            "query",
            "strand",
            "reference",
            "reference_start",
            "query_sequence",
        ]
        if not all(field in self.names_search_output for field in required_fields):
            raise ValueError(
                f"Some of the required fields {required_fields} are missing in the search results."
            )
        table_hits["reference_end"] = table_hits.apply(
            lambda x: x["reference_start"] + len(x["query_sequence"]), axis=1
        )
        bed = pd.DataFrame(
            {
                "chr": table_hits["reference"],
                "start": table_hits["reference_start"],
                "end": table_hits["reference_end"],
                "name": table_hits["query"],
                "score": 0,
                "strand": table_hits["strand"],
            }
        )
        file_bed = os.path.join(self.dir_output, f"references_{region_id}.bed")
        bed.to_csv(file_bed, sep="\t", index=False, header=False)

        references_fasta_file = os.path.join(self.dir_output, f"references_{region_id}.fasta")

        get_sequence_from_annotation(
            file_bed, file_reference, references_fasta_file, strand=True, nameOnly=True
        )
        references = [off_reference.seq for off_reference in SeqIO.parse(references_fasta_file, "fasta")]
        os.remove(references_fasta_file)
        os.remove(file_bed)
        return references

    def _add_alignement_gaps(self, table_hits: pd.DataFrame, queries: list, references: list):
        """Adjust the sequences of the oligos and the gaps found by the alignement search.
        The gapped references and queries are are defined such that nucleotides in the same index are binding.

        :param table_hits: DataFrame with the oligos that have a match with the query oligo.
        :type table_hits: pandas.DataFrame
        :param database_region: Dictionary with the information of the oligo sequences.
        :type database_region: dict
        """

        # bowtie does not support gaps
        return queries, references


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

    The hits returned by Bowtie2 can be further filtered using machine learning models. For more information regarding which filters are available
    refer to https://github.com/HelmholtzAI-Consultants-Munich/oligo-designer-toolsuite-AI-filters

    :param search_parameters: Custom parameters for the Bowtie2 search command.
    :type search_parameters: dict
    :param hit_parameters: Criteria to consider a Bowtie hit significant for filtering.
    :type hit_parameters: dict
    :param names_search_output: Column names for parsing Bowtie2 search output.
    :type names_search_output: list
    :param filter_name: Subdirectory path for the output, i.e. <dir_output>/<filter_name>, defaults to "bowtie2_filter".
    :type filter_name: str, optional
    :param dir_output: Directory for saving intermediate files, defaults to "output"
    :type dir_output: str, optional
    """

    def __init__(
        self,
        search_parameters: dict = {},
        hit_parameters: dict = {},
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
        filter_name: str = "bowtie2_filter",
        dir_output: str = "output",
    ):
        """Constructor for the Bowtie2Filter class."""
        super().__init__(filter_name, dir_output)

        self.search_parameters = search_parameters
        self.hit_parameters = hit_parameters  # currently not used
        self.names_search_output = names_search_output

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
        process = subprocess.Popen(cmd, shell=True, cwd=self.dir_output).wait()

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
        file_bowtie_results = os.path.join(self.dir_output, f"bowtie2_results_{region_name}.txt")

        cmd_parameters = ""
        for parameter, value in self.search_parameters.items():
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
        process = subprocess.Popen(cmd, shell=True, cwd=self.dir_output).wait()

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
        :return: A tuple containing a DataFrame of significant Bowtie2 hits.
        :rtype: pd.DataFrame
        """
        if not consider_hits_from_input_region:
            # remove all hits where query and reference come from the same region
            search_results = search_results[
                search_results["query_region_id"] != search_results["reference_region_id"]
            ]

        return search_results

    def _get_references(self, search_results: pd.DataFrame, file_reference: str, region_id: str):
        raise NotImplementedError("AI filters not supported for Bowtie2.")

    def _add_alignement_gaps(self, search_results: pd.DataFrame, queries: List[Seq], references: List[Seq]):
        raise NotImplementedError("AI filters not supported for Bowtie2.")
