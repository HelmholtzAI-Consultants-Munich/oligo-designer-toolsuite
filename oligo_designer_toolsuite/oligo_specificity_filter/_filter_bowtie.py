############################################
# Imports
############################################

import os
import subprocess
from typing import List, Union, Tuple

import pandas as pd
from Bio import SeqIO

from oligo_designer_toolsuite._constants import _TYPES_SEQ
from oligo_designer_toolsuite.database import OligoDatabase
from oligo_designer_toolsuite.oligo_specificity_filter import AlignmentSpecificityFilter

from ..utils._checkers_and_helpers import check_if_list
from ..utils._sequence_processor import get_sequence_from_annotation

############################################
# Oligo Bowtie Filter Classes
############################################


class BowtieFilter(AlignmentSpecificityFilter):
    """
    A filter class that uses Bowtie for sequence alignment and filtering.

    The `BowtieFilter` class is designed to align sequences against a reference database using Bowtie.
    This class manages Bowtie search parameters and interprets the alignment results according to user-defined criteria.
    It can be used to filter sequences based on their specificity and alignment properties.

    Bowtie (1.3 or higher) can be installed via the Bowtie webpage (https://bowtie-bio.sourceforge.net/manual.shtml#obtaining-bowtie)
    or via Bioconda (http://bioconda.github.io/recipes/bowtie/README.html) with:
    ``conda config --add channels conda-forge``
    ``conda config --add channels bioconda``
    ``conda update --all``
    ``conda install "bowtie>=1.3.1"``

    Useful Bowtie search parameters:
    --norc: only report hits on the plus strand
    --nofw: only report hits on the minus strand
    -v <int>: Report alignments with at most <int> mismatches.
    -n <int>: Maximum number of mismatches permitted in the “seed”, i.e., the first L base pairs of the read (where L is set with -l/--seedlen). This may be 0, 1, 2, or 3 and the default is 2.
    -l <int>: The “seed length”; i.e., the number of bases on the high-quality end of the read to which the -n ceiling applies. The lowest permitted setting is 5 and the default is 28. Bowtie is faster for larger values of -l.
    All available Bowtie search parameters are listed on the Bowtie webpage.

    The hits returned by Bowtie can be further filtered using machine learning models. For more information regarding which filters are available
    refer to https://github.com/HelmholtzAI-Consultants-Munich/oligo-designer-toolsuite-AI-filters.

    :param search_parameters: Parameters to configure the Bowtie search.
    :type search_parameters: dict
    :param hit_parameters: Criteria for interpreting Bowtie hits (not utilized in this filter).
    :type hit_parameters: dict
    :param names_search_output: List of names for the Bowtie search output fields.
    :type names_search_output: list
    :param filter_name: Name of the filter for identification purposes.
    :type filter_name: str
    :param dir_output: Directory to store output files and temporary data.
    :type dir_output: str
    """

    def __init__(
        self,
        search_parameters: dict = {},
        hit_parameters: dict = {},  # not utilized in this filter
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
    ) -> None:
        """Constructor for the BowtieFilter class."""
        super().__init__(filter_name, dir_output)

        self.search_parameters = search_parameters
        self.hit_parameters = hit_parameters  # not utilized in this filter
        self.names_search_output = names_search_output

    def _create_index(self, file_reference: str, n_jobs: int) -> str:
        """
        Creates a Bowtie index for a given reference file.

        This function generates a Bowtie index from a reference file, which is necessary for performing sequence alignment.
        The index creation is parallelized across multiple threads to optimize performance.

        :param file_reference: The path to the reference file for which the Bowtie index will be created.
        :type file_reference: str
        :param n_jobs: The number of parallel jobs to use for creating the index.
        :type n_jobs: int
        :return: The name of the created Bowtie index file.
        :rtype: str
        """
        ## Create bowtie index
        file_reference = os.path.abspath(file_reference)
        filename_reference_index = os.path.basename(file_reference)

        cmd = (
            "bowtie-build --offrate 4"
            + " --threads "
            + str(n_jobs)
            + " -f "
            + file_reference
            + " "
            + filename_reference_index
        )
        process = subprocess.Popen(cmd, shell=True, cwd=self.dir_output, stdout=subprocess.DEVNULL).wait()
        return filename_reference_index

    def _run_search(
        self,
        oligo_database: OligoDatabase,
        file_index: str,
        sequence_type: _TYPES_SEQ,
        region_ids: Union[str, List[str]] = None,
    ) -> pd.DataFrame:
        """
        Runs a Bowtie search against a reference index using sequences from the oligo database.

        This function performs a Bowtie search to align oligonucleotide sequences from the specified
        oligo database to a reference genome index. The results are processed and returned as a DataFrame.

        :param oligo_database: The Oligo Database containing the oligonucleotides and their associated attributes.
        :type oligo_database: OligoDatabase
        :param file_index: The name of the Bowtie database index file to search against.
        :type file_index: str
        :param sequence_type: The type of sequence to be used for the filter calculations.
        :type sequence_type: _TYPES_SEQ["oligo", "target"]
        :param region_ids: List of region IDs to process. If None, all regions in the database are processed, default is None.
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
            save_description=False,
            region_ids=region_ids,
            sequence_type=sequence_type,
        )
        file_bowtie_results = os.path.join(self.dir_output, f"bowtie_results_{region_name}.txt")

        cmd_parameters = ""
        for parameter, value in self.search_parameters.items():
            cmd_parameters += f" {parameter} {value}"

        cmd = (
            "bowtie"
            + " -x "
            + file_index
            + " -f"  # fasta file is input
            + " -a"  # report all alignments -> TODO: does this make sense or set e.g. -k 100
            + cmd_parameters
            + " "
            + file_oligo_database
            + " "
            + file_bowtie_results
        )
        process = subprocess.Popen(cmd, shell=True, cwd=self.dir_output, stdout=subprocess.DEVNULL).wait()

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
        region_ids: Union[str, List[str]],  # not used in this filter
    ) -> pd.DataFrame:
        """
        Filters Bowtie search results based on whether the query and reference sequences come from different regions.

        This method processes the Bowtie alignment results to remove hits where the query and reference sequences
        originate from the same region, if specified. This is useful for excluding self-hits or hits within the same region,
        depending on the configuration.

        :param oligo_database: The Oligo Database containing the oligonucleotides and their associated attributes (not utilized in this filter).
        :type oligo_database: OligoDatabase
        :param search_results: DataFrame containing the results of the Bowtie search.
        :type search_results: pd.DataFrame
        :param consider_hits_from_input_region: Whether to include hits from the same region as the query.
        :type consider_hits_from_input_region: bool
        :param region_ids: List of region IDs to process (not utilized in this filter).
        :type region_ids: Union[str, List[str]]
        :return: A DataFrame containing the filtered Bowtie search hits.
        :rtype: pd.DataFrame
        """
        if not consider_hits_from_input_region:
            # remove all hits where query and reference come from the same region
            search_results = search_results[
                search_results["query_region_id"] != search_results["reference_region_id"]
            ]

        return search_results

    def _get_references(self, table_hits: pd.DataFrame, file_reference: str, region_id: str) -> list:
        """
        Extracts reference sequences from the Bowtie search results and returns them as a list.

        This function processes Bowtie alignment results to determine the start and end positions of the reference sequences aligned to the queries.
        It then uses these coordinates to extract the corresponding sequences from the reference genome file.

        :param table_hits: DataFrame containing Bowtie search hits with alignment information.
        :type table_hits: pd.DataFrame
        :param file_reference: Path to the reference genome file.
        :type file_reference: str
        :param region_id: Region ID to process.
        :type region_id: str
        :return: A list of reference sequences.
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

    def _add_alignment_gaps(
        self, table_hits: pd.DataFrame, queries: list, references: list
    ) -> Tuple[list, list]:
        """
        Handles the addition of alignment gaps for queries and references.

        This function is designed to process alignment gaps in sequences to ensure that
        the query and reference sequences are properly aligned with gaps inserted where necessary.
        However, since Bowtie does not support gaps in alignments, this implementation simply
        returns the original query and reference sequences without modification.

        :param table_hits: DataFrame containing information about the alignment hits, including gap positions (not utilized in this filter).
        :type table_hits: pd.DataFrame
        :param queries: List of query sequences to be aligned.
        :type queries: list
        :param references: List of reference sequences to be aligned.
        :type references: list
        :return: Unmodified lists of query and reference sequences.
        :rtype: Tuple[list, list]
        """
        # bowtie does not support gaps
        return queries, references


############################################
# Oligo Bowtie2 Filter Classes
############################################


class Bowtie2Filter(AlignmentSpecificityFilter):
    """
    A filter class that utilizes Bowtie2 for alignment-based specificity filtering.

    The `Bowtie2Filter` class is designed to align sequences against a ReferenceDatabase using Bowtie2.
    This class manages Bowtie2 search parameters and interprets the alignment results according to user-defined criteria.
    It can be used to filter sequences based on their specificity and alignment properties.

    Bowtie2 (2.5 or higher) can be installed via the Bowtie2 webpage (https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#obtaining-bowtie-2)
    or via Bioconda (http://bioconda.github.io/recipes/bowtie2/README.html) with:
    ``conda config --add channels conda-forge``
    ``conda config --add channels bioconda``
    ``conda update --all``
    ``conda install "bowtie2>=2.5"``

    Useful Bowtie2 search parameters:
    --norc: only report hits on the plus strand
    --nofw: only report hits on the minus strand
    -N <int>: Sets the number of mismatches allowed in a seed alignment during multiseed alignment. Can be set to 0 or 1. Setting this higher makes alignment slower (often much slower) but increases sensitivity. Default: 0.
    -L <int>: Sets the length of the seed substrings to align during multiseed alignment. Smaller values make alignment slower but more sensitive. Default: the --sensitive preset is used by default, which sets -L to 22 and 20 in --end-to-end mode and in --local mode.
    All available Bowtie2 search parameters are listed on the Bowtie2 webpage.

    The hits returned by Bowtie2 can be further filtered using machine learning models. For more information regarding which filters are available
    refer to https://github.com/HelmholtzAI-Consultants-Munich/oligo-designer-toolsuite-AI-filters.

    :param search_parameters: Parameters to configure the Bowtie2 search.
    :type search_parameters: dict
    :param hit_parameters: Criteria for interpreting Bowtie2 hits (not utilized in this filter).
    :type hit_parameters: dict
    :param names_search_output: List of names for the Bowtie2 search output fields.
    :type names_search_output: list
    :param filter_name: Name of the filter for identification purposes.
    :type filter_name: str
    :param dir_output: Directory to store output files and temporary data.
    :type dir_output: str
    """

    def __init__(
        self,
        search_parameters: dict = {},
        hit_parameters: dict = {},  # not utilized in this filter
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
    ) -> None:
        """Constructor for the Bowtie2Filter class."""
        super().__init__(filter_name, dir_output)

        self.search_parameters = search_parameters
        self.hit_parameters = hit_parameters  # not utilized in this filter
        self.names_search_output = names_search_output

    def _create_index(self, file_reference: str, n_jobs: int) -> str:
        """
        Creates an index for the Bowtie2 alignment tool using a specified reference file.

        This method generates a Bowtie2 index from a reference file, which is necessary for performing sequence alignment.
        The index creation is parallelized across multiple threads to optimize performance.

        :param file_reference: The path to the reference file for which the Bowtie2 index will be created.
        :type file_reference: str
        :param n_jobs: The number of parallel jobs to use for creating the index.
        :type n_jobs: int
        :return: The name of the created Bowtie2 index file.
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
        oligo_database: OligoDatabase,
        file_index: str,
        sequence_type: _TYPES_SEQ,
        region_ids: Union[str, List[str]] = None,
    ) -> pd.DataFrame:
        """
        Runs a Bowtie2 search against a reference index using sequences from the OligoDatabase.

        This function performs a Bowtie2 search to align oligonucleotide sequences from the specified
        OligoDatabase to a reference genome index. The results are processed and returned as a DataFrame.

        :param oligo_database: The OligoDatabase containing the oligonucleotides and their associated attributes.
        :type oligo_database: OligoDatabase
        :param file_index: The name of the Bowtie2 database index file to search against.
        :type file_index: str
        :param sequence_type: The type of sequence to be used for the filter calculations.
        :type sequence_type: _TYPES_SEQ["oligo", "target"]
        :param region_ids: List of region IDs to process. If None, all regions in the OligoDatabase are processed, default is None.
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
            save_description=False,
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
        region_ids: Union[str, List[str]],  # not used in this filter
    ) -> pd.DataFrame:
        """
        Filters Bowtie2 search results based on whether the query and reference sequences come from different regions.

        This method processes the Bowtie2 alignment results to remove hits where the query and reference sequences
        originate from the same region, if specified. This is useful for excluding self-hits or hits within the same region,
        depending on the configuration.

        :param oligo_database: The OligoDatabase containing the oligonucleotides and their associated attributes (not utilized in this filter).
        :type oligo_database: OligoDatabase
        :param search_results: DataFrame containing the results of the Bowtie2 search.
        :type search_results: pd.DataFrame
        :param consider_hits_from_input_region: Whether to include hits from the same region as the query.
        :type consider_hits_from_input_region: bool
        :param region_ids: List of region IDs to process (not utilized in this filter).
        :type region_ids: Union[str, List[str]]
        :return: A DataFrame containing the filtered Bowtie2 search hits.
        :rtype: pd.DataFrame
        """
        if not consider_hits_from_input_region:
            # remove all hits where query and reference come from the same region
            search_results = search_results[
                search_results["query_region_id"] != search_results["reference_region_id"]
            ]

        return search_results

    def _get_references(self, table_hits: pd.DataFrame, file_reference: str, region_id: str) -> list:
        """
        Raises an error indicating that AI filters are not supported for Bowtie2.

        This method is intended to retrieve reference sequences based on search results.
        However, in the context of Bowtie2, this method is not implemented and raises a `NotImplementedError`.

        :param table_hits: DataFrame containing Bowtie2 search hits with alignment information.
        :type table_hits: pd.DataFrame
        :param file_reference: Path to the reference genome file.
        :type file_reference: str
        :param region_id: Region ID to process.
        :type region_id: str
        :raises NotImplementedError: Always, because AI filters are not supported for Bowtie2.
        """
        raise NotImplementedError("AI filters not supported for Bowtie2.")

    def _add_alignment_gaps(
        self, search_results: pd.DataFrame, queries: list, references: list
    ) -> Tuple[list, list]:
        """
        Raises an error indicating that AI filters are not supported for Bowtie2.

        This method is intended to add alignment gaps to sequences based on search results.
        However, in the context of Bowtie2, this method is not implemented and raises a `NotImplementedError`.

        :param table_hits: DataFrame containing information about the alignment hits, including gap positions.
        :type table_hits: pd.DataFrame
        :param queries: List of query sequences to be aligned.
        :type queries: list
        :param references: List of reference sequences to be aligned.
        :type references: list
        :raises NotImplementedError: Always, because AI filters are not supported for Bowtie2.
        """

        raise NotImplementedError("AI filters not supported for Bowtie2.")
