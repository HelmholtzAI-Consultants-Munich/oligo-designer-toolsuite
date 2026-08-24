############################################
# imports
############################################

import os
import subprocess
from abc import abstractmethod

import pandas as pd

from oligo_designer_toolsuite._exceptions import ConfigurationError, DatabaseError
from oligo_designer_toolsuite.database import OligoDatabase
from oligo_designer_toolsuite.oligo_property_calculator import (
    BaseProperty,
    PropertyCalculator,
    SeedregionProperty,
    SeedregionSiteProperty,
)
from oligo_designer_toolsuite.oligo_specificity_filter import AlignmentSpecificityFilter
from oligo_designer_toolsuite.utils import logger
from oligo_designer_toolsuite.utils._checkers_and_helpers import safe_append_filename

############################################
# Oligo Blast Filter Classes
############################################


class BlastNFilter(AlignmentSpecificityFilter):
    """
    A class for filtering oligonucleotide sequences using BLASTN alignments.

    The `BlastNFilter` class is designed to align sequences against a ReferenceDatabase using BLASTN.
    This class manages BLAST search parameters and interprets the alignment results according to user-defined criteria.
    It can be used to filter sequences based on their specificity and alignment properties.

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

    :param remove_hits: If True, oligos overlapping variants are removed. If False, they are flagged.
    :type remove_hits: bool
    :param search_parameters: Parameters to configure the BLAST search.
    :type search_parameters: dict
    :param hit_parameters: Criteria for interpreting BLAST hits.
    :type hit_parameters: dict
    :param names_search_output: List of names for the BLAST search output fields.
    :type names_search_output: list
    :param filter_name: Name of the filter for identification purposes.
    :type filter_name: str
    :param dir_output: Directory path where output files will be saved.
    :type dir_output: str
    """

    def __init__(
        self,
        remove_hits: bool = True,
        search_parameters: dict = {},
        hit_parameters: dict = {},
        names_search_output: list = [
            "query",
            "reference",
            "alignment_length",
            "query_start",
            "query_end",
            "query_length",
        ],
        filter_name: str = "blast_filter",
        dir_output: str = "output",
    ) -> None:
        """Constructor for the BlastNFilter class."""
        super().__init__(remove_hits, filter_name, dir_output)

        self.search_parameters = search_parameters
        self.hit_parameters = hit_parameters
        self.names_search_output = names_search_output

        # Define default output format for blast search filter. The fields are:
        # query, reference, alignment_length, query_start, query_end, query_length
        if "-outfmt" not in self.search_parameters.keys():
            self.search_parameters["-outfmt"] = "6 qseqid sseqid length qstart qend qlen"

    def _create_reference(
        self,
        n_jobs: int,  # not utilized in this filter
    ) -> str:
        """
        Creates a BLAST index for a given reference file.

        This method generates an index for a nucleotide BLAST database from the specified reference file.
        The index is stored in the specified output directory.

        :param n_jobs: Number of parallel jobs to use for processing. Note: This parameter is not utilized in this filter.
        :type n_jobs: int
        :return: The name of the created BLAST reference file.
        :rtype: str
        """
        if self.reference_database is None:
            raise DatabaseError("reference_database must be set before calling _create_reference")

        # write refrence database to fasta
        file_reference = self.reference_database.write_database_to_file(
            filename=f"db_reference_{self.filter_name}",
            dir_output=self.dir_output,
        )

        ## Create blast index
        args = ["makeblastdb", "-dbtype", "nucl", "-out", file_reference, "-in", file_reference]

        subprocess.run(args, cwd=self.dir_output, check=True, stdout=subprocess.DEVNULL)

        return file_reference

    def _run_search(
        self,
        oligo_database: OligoDatabase,
        file_reference: str,
        region_id: str,
    ) -> pd.DataFrame:
        """
        Runs a BLAST search for the specified oligo sequences against a ReferenceDatabase.

        This function writes the sequences from the OligoDatabase to a FASTA file and
        performs a BLAST search against the provided ReferenceDatabase index.
        The results are read into a DataFrame for further analysis.

        :param oligo_database: The OligoDatabase instance containing oligonucleotide sequences and their associated properties. This database stores oligo data organized by genomic regions and can be used for filtering, property calculations, set generation, and output operations.
        :type oligo_database: OligoDatabase
        :param file_reference: Path to the reference file used for alignment filtering.
        :type file_reference: str
        :param region_id: Region ID to process.
        :type region_id: str
        :return: A DataFrame containing the BLAST search results.
        :rtype: pd.DataFrame
        """
        if self.sequence_type is None:
            raise ConfigurationError("sequence_type must be set before calling _run_search")

        file_oligo_database = oligo_database.write_database_to_fasta(
            sequence_type=self.sequence_type,
            filename=f"oligo_database_blast_{region_id}",
            dir_output=self.dir_output,
            save_description=False,
            region_ids=region_id,
        )
        file_blast_results = safe_append_filename(self.dir_output, f"blast_results_{region_id}.txt")

        args = [
            "blastn",
            "-query",
            file_oligo_database,
            "-out",
            file_blast_results,
            "-db",
            file_reference,
        ]

        for parameter, value in self.search_parameters.items():
            args.append(parameter)
            if str(value) != "":
                args.append(str(value))

        subprocess.run(args, cwd=self.dir_output, check=True, stdout=subprocess.DEVNULL)

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
        oligo_database: OligoDatabase,  # not used in this filter
        search_results: pd.DataFrame,
        consider_hits_from_input_region: bool,
        region_id: str,  # not used in this filter
    ) -> pd.DataFrame:
        """
        Identifies significant hits from BLAST search results based on alignment length or coverage.

        This function processes the BLAST search results to identify significant hits,
        considering either a minimum alignment length or a coverage percentage of the query sequence.
        Additionally, it can exclude hits where the query and reference sequences originate from the same region.

        :param oligo_database: The OligoDatabase instance containing oligonucleotide sequences and their associated properties. This database stores oligo data organized by genomic regions and can be used for filtering, property calculations, set generation, and output operations. Note: This parameter is not utilized in this filter.
        :type oligo_database: OligoDatabase
        :param search_results: DataFrame containing the results of the BLAST search.
        :type search_results: pd.DataFrame
        :param consider_hits_from_input_region: Whether to include hits from the same region as the query.
        :type consider_hits_from_input_region: bool
        :param region_id: Region ID to process (not utilized in this filter).
        :type region_id: str
        :return: A DataFrame containing the filtered BLAST search hits.
        :rtype: pd.DataFrame
        """
        if "min_alignment_length" in self.hit_parameters.keys():
            if "coverage" in self.hit_parameters.keys():
                logger.warning(
                    "Both, 'min_alignment_length' and 'coverage' parameters were provided. Using 'min_alignment_length' parameter."
                )
            min_alignment_length = self.hit_parameters["min_alignment_length"]
        elif "coverage" in self.hit_parameters.keys():
            min_alignment_length = search_results["query_length"] * self.hit_parameters["coverage"] / 100
        else:
            raise ConfigurationError(
                "Please provide either 'coverage' or a 'min_alignment_length' in hit_parameters!"
            )

        search_results["min_alignment_length"] = min_alignment_length

        if not consider_hits_from_input_region:
            # remove all hits where query and reference come from the same region
            search_results = search_results[
                search_results["query_region_id"] != search_results["reference_region_id"]
            ]

        blast_table_hits = search_results.loc[
            search_results.alignment_length > search_results.min_alignment_length
        ]

        return blast_table_hits


############################################
# Oligo Blast Filter with Seedregion Classes
############################################


class BlastNSeedregionFilterBase(BlastNFilter):
    """
    A base class for BLASTN-based seed region filters in oligonucleotide design.

    The `BlastNSeedregionFilterBase` class extends the `BlastNFilter` class,
    providing a foundation for implementing specific filters that utilize BLASTN
    to assess the specificity of seed regions in oligonucleotide sequences.

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

    :param remove_hits: If True, oligos overlapping variants are removed. If False, they are flagged.
    :type remove_hits: bool
    :param search_parameters: Parameters to configure the BLAST search.
    :type search_parameters: dict
    :param hit_parameters: Criteria for interpreting BLAST hits.
    :type hit_parameters: dict
    :param names_search_output: List of names for the BLAST search output fields.
    :type names_search_output: list
    :param filter_name: Name of the filter for identification purposes.
    :type filter_name: str
    :param dir_output: Directory path where output files will be saved.
    :type dir_output: str
    """

    def __init__(
        self,
        remove_hits: bool = True,
        search_parameters: dict | None = None,
        hit_parameters: dict | None = None,
        names_search_output: list | None = None,
        filter_name: str = "blast_filter",
        dir_output: str = "output",
    ) -> None:
        """Constructor for the BlastNSeedregionFilterBase class."""
        if not search_parameters:
            search_parameters = {}
        if not hit_parameters:
            hit_parameters = {}
        if not names_search_output:
            names_search_output = [
                "query",
                "reference",
                "alignment_length",
                "query_start",
                "query_end",
                "query_length",
            ]
        super().__init__(
            remove_hits,
            search_parameters,
            hit_parameters,
            names_search_output,
            filter_name,
            dir_output,
        )

    @abstractmethod
    def _add_seed_region_information(
        self, oligo_database: OligoDatabase, search_results: pd.DataFrame, region_id: str
    ) -> pd.DataFrame:
        """
        An abstract method to add seed region information to BLAST search results.

        The `_add_seed_region_information` method is intended to be implemented in subclasses.
        It processes the BLAST search results and integrates specific seed region data into the provided OligoDatabase.

        :param oligo_database: The OligoDatabase instance containing oligonucleotide sequences and their associated properties. This database stores oligo data organized by genomic regions and can be used for filtering, property calculations, set generation, and output operations.
        :type oligo_database: OligoDatabase
        :param search_results: The DataFrame containing the results of the BLAST search.
        :type search_results: pd.DataFrame
        :param region_id: Region ID to process.
        :type region_id: str
        :return: The BLAST search results with added seed region information.
        :rtype: pd.DataFrame
        """

    def _find_hits(
        self,
        oligo_database: OligoDatabase,
        search_results: pd.DataFrame,
        consider_hits_from_input_region: bool,
        region_id: str,
    ) -> pd.DataFrame:
        """
        Finds and filters hits from BLAST search results based on alignment length, coverage, and seed region information.

        This method processes the search results to identify significant hits by applying user-defined thresholds for alignment length or coverage.
        It also incorporates seed region information into the search results and optionally filters out hits that originate from the same region as the input sequence.

        :param oligo_database: The OligoDatabase instance containing oligonucleotide sequences and their associated properties. This database stores oligo data organized by genomic regions and can be used for filtering, property calculations, set generation, and output operations.
        :type oligo_database: OligoDatabase
        :param search_results: DataFrame containing BLAST search results.
        :type search_results: pd.DataFrame
        :param consider_hits_from_input_region: Whether to include hits from the same region as the input sequence.
        :type consider_hits_from_input_region: bool
        :param region_id: Region ID to process.
        :type region_id: str
        :return: Filtered BLAST search results containing significant hits.
        :rtype: pd.DataFrame
        """
        if "min_alignment_length" in self.hit_parameters.keys():
            if "coverage" in self.hit_parameters.keys():
                logger.warning(
                    "Both, 'min_alignment_length' and 'coverage' parameters were provided. Using 'min_alignment_length' parameter."
                )
            min_alignment_length = self.hit_parameters["min_alignment_length"]
        elif "coverage" in self.hit_parameters.keys():
            min_alignment_length = search_results["query_length"] * self.hit_parameters["coverage"] / 100
        else:
            raise ConfigurationError(
                "Please provide either 'coverage' or a 'min_alignment_length' in hit_parameters!"
            )

        search_results["min_alignment_length"] = min_alignment_length

        search_results = self._add_seed_region_information(
            oligo_database=oligo_database, search_results=search_results, region_id=region_id
        )

        # if seedregion not given
        search_results["seedregion_start"] = search_results["seedregion_start"].fillna(0)
        search_results["seedregion_end"] = search_results["seedregion_end"].fillna(
            search_results["query_length"]
        )

        if not consider_hits_from_input_region:
            # remove all hits where query and reference come from the same region
            search_results = search_results[
                search_results["query_region_id"] != search_results["reference_region_id"]
            ]

        blast_table_hits = search_results.loc[
            (search_results.alignment_length > search_results.min_alignment_length)
            | (
                (search_results.query_start < search_results.seedregion_start)
                & (search_results.query_end > search_results.seedregion_end)
            )
        ]

        return blast_table_hits


class BlastNSeedregionFilter(BlastNSeedregionFilterBase):
    """
    A filter class for BLAST searches that considers seed regions in oligonucleotide sequences.

    The `BlastNSeedregionFilter` class extends the functionality of BLAST searches by incorporating seed region information.
    It allows for precise filtering based on the alignment within specific regions of the oligonucleotide sequences,
    ensuring that only relevant alignments are retained.

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

    :param seedregion_start: The start position of the seed region within the sequence.
    :type seedregion_start: Union[int, float]
    :param seedregion_end: The end position of the seed region within the sequence.
    :type seedregion_end: Union[int, float]
    :param remove_hits: If True, oligos overlapping variants are removed. If False, they are flagged.
    :type remove_hits: bool
    :param search_parameters: Parameters to configure the BLAST search.
    :type search_parameters: dict
    :param hit_parameters: Criteria for interpreting BLAST hits.
    :type hit_parameters: dict
    :param names_search_output: List of names for the BLAST search output fields.
    :type names_search_output: list
    :param filter_name: Name of the filter for identification purposes.
    :type filter_name: str
    :param dir_output: Directory path where output files will be saved.
    :type dir_output: str
    """

    def __init__(
        self,
        seedregion_start: int | float,
        seedregion_end: int | float,
        remove_hits: bool = True,
        search_parameters: dict | None = None,
        hit_parameters: dict | None = None,
        names_search_output: list | None = None,
        filter_name: str = "blast_filter",
        dir_output: str = "output",
    ) -> None:
        """Constructor for the BlastNSeedregionFilter class."""
        if not search_parameters:
            search_parameters = {}
        if not hit_parameters:
            hit_parameters = {}
        if not names_search_output:
            names_search_output = [
                "query",
                "reference",
                "alignment_length",
                "query_start",
                "query_end",
                "query_length",
            ]
        super().__init__(
            remove_hits,
            search_parameters,
            hit_parameters,
            names_search_output,
            filter_name,
            dir_output,
        )

        self.seedregion_start = seedregion_start
        self.seedregion_end = seedregion_end

    def _add_seed_region_information(
        self, oligo_database: OligoDatabase, search_results: pd.DataFrame, region_id: str
    ) -> pd.DataFrame:
        """
        Adds seed region information to the BLAST search results.

        This method enhances the BLAST search results by incorporating seed region start and end positions from the OligoDatabase.
        It calculates the seed region positions for each oligonucleotide and merges this information with the BLAST search results to refine the output.

        :param oligo_database: The OligoDatabase instance containing oligonucleotide sequences and their associated properties. This database stores oligo data organized by genomic regions and can be used for filtering, property calculations, set generation, and output operations.
        :type oligo_database: OligoDatabase
        :param search_results: The DataFrame containing the results of the BLAST search.
        :type search_results: pd.DataFrame
        :param region_id: Region ID to process.
        :type region_id: str
        :return: The BLAST search results with added seed region information.
        :rtype: pd.DataFrame
        """
        if self.sequence_type is None:
            raise ConfigurationError("sequence_type must be set before calling _add_seed_region_information")

        properties: list[BaseProperty] = [
            SeedregionProperty(start=self.seedregion_start, end=self.seedregion_end)
        ]
        calculator = PropertyCalculator(properties=properties)
        oligo_database = calculator.apply(
            oligo_database=oligo_database, sequence_type=self.sequence_type, n_jobs=1
        )

        seedregion = oligo_database.get_oligo_property_table(
            properties=["seedregion_start", "seedregion_end"],
            flatten=True,
            region_ids=region_id,
        )
        search_results = pd.merge(
            left=search_results,
            right=seedregion,
            left_on="query",
            right_on="oligo_id",
            how="left",
        )
        return search_results


class BlastNSeedregionSiteFilter(BlastNSeedregionFilterBase):
    """
    A filter class for BLASTN alignment considering seed region and seed region site proximity.

    The `BlastNSeedregionSiteFilter` class extends the `BlastNSeedregionFilterBase` to focus on evaluating oligonucleotide
    sequences based on their alignment within a specified seed region size, particularly around the seed region site.
    This is useful for ensuring that sequences align in a biologically meaningful region.

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

    :param seedregion_size: The size of the seed region around the seed region site to consider.
    :type seedregion_size: int
    :param seedregion_site_name: The property name of the seed region site stored in the OligoDatabase.
    :type seedregion_site_name: str
    :param remove_hits: If True, oligos overlapping variants are removed. If False, they are flagged.
    :type remove_hits: bool
    :param search_parameters: Parameters to configure the BLAST search.
    :type search_parameters: dict
    :param hit_parameters: Criteria for interpreting BLAST hits.
    :type hit_parameters: dict
    :param names_search_output: List of names for the BLAST search output fields.
    :type names_search_output: list
    :param filter_name: Name of the filter for identification purposes.
    :type filter_name: str
    :param dir_output: Directory path where output files will be saved.
    :type dir_output: str
    """

    def __init__(
        self,
        seedregion_size: int,
        seedregion_site_name: str,
        remove_hits: bool = True,
        search_parameters: dict = {},
        hit_parameters: dict = {},
        names_search_output: list = [
            "query",
            "reference",
            "alignment_length",
            "query_start",
            "query_end",
            "query_length",
        ],
        filter_name: str = "blast_filter",
        dir_output: str = "output",
    ) -> None:
        """Constructor for the BlastNSeedregionSiteFilter class."""
        super().__init__(
            remove_hits,
            search_parameters,
            hit_parameters,
            names_search_output,
            filter_name,
            dir_output,
        )
        self.seedregion_size = seedregion_size
        self.seedregion_site_name = seedregion_site_name

    def _add_seed_region_information(
        self, oligo_database: OligoDatabase, search_results: pd.DataFrame, region_id: str
    ) -> pd.DataFrame:
        """
        Adds seed region information to the BLASTN search results based on the seed region site.

        This method calculates the seed region around the seed region site for each oligonucleotide in the OligoDatabase and
        merges this information with the BLASTN search results. The seed region is defined by the `seedregion_size` parameter.

        :param oligo_database: The OligoDatabase instance containing oligonucleotide sequences and their associated properties. This database stores oligo data organized by genomic regions and can be used for filtering, property calculations, set generation, and output operations.
        :type oligo_database: OligoDatabase
        :param search_results: The DataFrame containing the results of the BLASTN search.
        :type search_results: pd.DataFrame
        :param region_id: Region ID to process.
        :type region_id: str
        :return: The BLAST search results with added seed region information.
        :rtype: pd.DataFrame
        """
        if self.sequence_type is None:
            raise ConfigurationError("sequence_type must be set before calling _add_seed_region_information")

        properties: list[BaseProperty] = [
            SeedregionSiteProperty(
                seedregion_size=self.seedregion_size,
                seedregion_site_name=self.seedregion_site_name,
            )
        ]
        calculator = PropertyCalculator(properties=properties)
        oligo_database = calculator.apply(
            oligo_database=oligo_database, sequence_type=self.sequence_type, n_jobs=1
        )

        seedregion = oligo_database.get_oligo_property_table(
            properties=["seedregion_start", "seedregion_end"],
            flatten=True,
            region_ids=region_id,
        )
        search_results = pd.merge(
            left=search_results,
            right=seedregion,
            left_on="query",
            right_on="oligo_id",
            how="left",
        )
        return search_results
