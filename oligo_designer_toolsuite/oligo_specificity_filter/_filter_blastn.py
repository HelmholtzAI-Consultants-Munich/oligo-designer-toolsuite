############################################
# imports
############################################

import os
import subprocess
import warnings
from abc import abstractmethod
from typing import Tuple, Union

import numpy as np
import pandas as pd

from abc import abstractmethod
from Bio import SeqIO

from oligo_designer_toolsuite._constants import _TYPES_SEQ
from oligo_designer_toolsuite.database import OligoAttributes, OligoDatabase
from oligo_designer_toolsuite.oligo_specificity_filter import SpecificityFilterAlignment

from ..utils._sequence_processor import get_sequence_from_annotation

############################################
# Oligo Blast Filter Classes
############################################


class BlastNFilter(SpecificityFilterAlignment):
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

    The hits returned by BLASTN can be further filtered using machine learning models. For more information regarding which filters are available
    refer to https://github.com/HelmholtzAI-Consultants-Munich/oligo-designer-toolsuite-AI-filters.

    :param sequence_type: The type of sequence to be used for the filter calculations.
    :type sequence_type: _TYPES_SEQ["oligo", "target"]
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
    :param dir_output: Directory to store output files and temporary data.
    :type dir_output: str
    """

    def __init__(
        self,
        sequence_type: _TYPES_SEQ,
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
        super().__init__(sequence_type, remove_hits, filter_name, dir_output)

        self.search_parameters = search_parameters
        self.hit_parameters = hit_parameters
        self.names_search_output = names_search_output

        # Define default output format for blast search filter. The fields are:
        # query, reference, alignment_length, query_start, query_end, query_length
        if "outfmt" not in self.search_parameters.keys():
            self.search_parameters["outfmt"] = "6 qseqid sseqid length qstart qend qlen"

    def create_reference(
        self,
        n_jobs: int,  # not utilized in this filter
    ) -> str:
        """
        Creates a BLAST index for a given reference file.

        This method generates an index for a nucleotide BLAST database from the specified reference file.
        The index is stored in the specified output directory.

        :param n_jobs: The number of parallel jobs to use for creating the index (not utilized in this filter).
        :type n_jobs: int
        :return: The name of the created BLAST reference file.
        :rtype: str
        """
        # write refrence database to fasta
        file_reference = self.reference_database.write_database_to_file(
            filename=f"db_reference_{self.filter_name}",
            dir_output=self.dir_output,
        )
        ## Create blast index
        cmd = "makeblastdb -dbtype nucl" + " -out " + file_reference + " -in " + file_reference
        process = subprocess.Popen(cmd, shell=True, cwd=self.dir_output, stdout=subprocess.DEVNULL).wait()

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

        :param oligo_database: The OligoDatabase containing the oligonucleotides and their associated attributes.
        :type oligo_database: OligoDatabase
        :param file_reference: Path to the reference file used for alignment filtering.
        :type file_reference: str
        :param region_id: Region ID to process.
        :type region_id: str
        :return: A DataFrame containing the BLAST search results.
        :rtype: pd.DataFrame
        """
        file_oligo_database = oligo_database.write_database_to_fasta(
            sequence_type=self.sequence_type,
            filename=f"oligo_database_blast_{region_id}",
            dir_output=self.dir_output,
            save_description=False,
            region_ids=region_id,
        )
        file_blast_results = os.path.join(self.dir_output, f"blast_results_{region_id}.txt")

        cmd_parameters = ""
        for parameter, value in self.search_parameters.items():
            # add quotes if list of strings seperated by whitespace
            value = f'"{value}"' if " " in str(value) else value
            cmd_parameters += f" -{parameter} {value}"

        cmd = (
            "blastn"
            + " -query "
            + file_oligo_database
            + " -out "
            + file_blast_results
            + " -db "
            + os.path.join(self.dir_output, file_reference)
            + " "
            + cmd_parameters
        )
        process = subprocess.Popen(cmd, shell=True, cwd=self.dir_output, stdout=subprocess.DEVNULL).wait()

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

        :param oligo_database: The OligoDatabase containing the oligonucleotides and their associated attributes (not utilized in this filter).
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
                warnings.warn(
                    "Both, 'min_alignment_length' and 'coverage' parameters were provided. Using 'min_alignment_length' parameter."
                )
            min_alignment_length = self.hit_parameters["min_alignment_length"]
        elif "coverage" in self.hit_parameters.keys():
            min_alignment_length = search_results["query_length"] * self.hit_parameters["coverage"] / 100
        else:
            raise KeyError("Please provide either 'coverage' or a 'min_alignment_length' in hit_parameters!")

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

    def _get_references(self, table_hits: pd.DataFrame, file_reference: str, region_id: str) -> list:
        """
        Generates reference sequences from BLAST hits and returns them in a padded format.

        This function takes a DataFrame of BLAST hits and processes it to generate a list of reference sequences.
        The function creates a BED file based on the hits, retrieves the corresponding sequences from the reference genome,
        and then pads the sequences as necessary to ensure they align correctly with the query sequences.

        :param table_hits: DataFrame containing BLAST search hits.
        :type table_hits: pd.DataFrame
        :param file_reference: Path to the reference genome file.
        :type file_reference: str
        :param region_id: Region ID to process.
        :type region_id: str
        :return: A list of padded reference sequences.
        :rtype: list
        """
        required_fields = [
            "query",
            "reference",
            "alignment_length",
            "query_start",
            "query_end",
            "query_length",
            "query_sequence",
            "reference_start",
            "reference_end",
            "reference_sequence",
            "reference_strand",
        ]
        if not all(field in self.names_search_output for field in required_fields):
            raise ValueError(
                f"Some of the required fields {required_fields} are missing in the search results."
            )
        # set the positions to a 0-based index
        table_hits = self._0_index_coordinates(table_hits)
        # Calculate adjusted "start" and "end" for BED format based on strand
        table_hits = self._extend_reference_start_end_coordinates(table_hits)

        # create bed file
        bed = pd.DataFrame(
            {
                "chr": table_hits["reference"],
                "start": table_hits["start"],
                "end": table_hits["end"],
                "name": table_hits["query"],
                "score": 0,
                "strand": table_hits["reference_strand"].map({"plus": "+", "minus": "-"}),
            }
        )

        # adjust for possible overflows (e.g. new coordinates are not included in the gene boundaries)
        # additionally we store how muchpadding we have to do to have two seqeunces of the same length
        bed = self._remove_overflows(bed, file_reference)
        file_bed = os.path.join(self.dir_output, f"references_{region_id}.bed")
        bed.to_csv(
            file_bed,
            sep="\t",
            index=False,
            header=False,
            columns=["chr", "start", "end", "name", "score", "strand"],
        )
        # generate the fasta file
        references_fasta_file = os.path.join(self.dir_output, f"references_{region_id}.fasta")
        get_sequence_from_annotation(
            file_bed, file_reference, references_fasta_file, strand=True, nameOnly=True
        )

        references = [off_reference.seq for off_reference in SeqIO.parse(references_fasta_file, "fasta")]
        references_padded = self._pad_overflows(bed, references)

        os.remove(references_fasta_file)
        os.remove(file_bed)
        return references_padded

    def _0_index_coordinates(self, table_hits: pd.DataFrame) -> pd.DataFrame:
        """
        Adjusts BLAST hit coordinates to a 0-based index system.

        This method modifies the coordinates in the BLAST hits DataFrame to follow a 0-based indexing system,
        which is commonly used in bioinformatics formats like BED files.
        The adjustments are made based on the strand orientation of the sequences,
        ensuring that the coordinates are accurate for both plus and minus strands.

        :param table_hits: DataFrame containing the BLAST hits with coordinates to be adjusted.
        :type table_hits: pd.DataFrame
        :return: DataFrame with adjusted 0-based coordinates.
        :rtype: pd.DataFrame
        """
        table_hits["query_start_corr"] = table_hits["query_start"] - 1
        table_hits["reference_start_corr"] = np.where(
            table_hits["reference_strand"] == "plus",
            table_hits["reference_start"] - 1,
            table_hits["reference_start"],
        )
        table_hits["reference_end_corr"] = np.where(
            table_hits["reference_strand"] == "minus",
            table_hits["reference_end"] - 1,
            table_hits["reference_end"],
        )
        return table_hits

    def _extend_reference_start_end_coordinates(self, table_hits: pd.DataFrame) -> pd.DataFrame:
        """
        Extends the reference start and end coordinates for BLAST hits.

        This function calculates and extends the start and end coordinates of reference sequences based on the alignment's strand orientation.
        The adjustments ensure the coordinates fully encompass the aligned query sequence in the reference genome.

        :param table_hits: DataFrame containing BLAST hits with corrected 0-based reference coordinates.
        :type table_hits: pd.DataFrame
        :return: DataFrame with extended start and end coordinates for reference sequences.
        :rtype: pd.DataFrame
        """
        table_hits["start"] = np.where(
            table_hits["reference_strand"] == "plus",
            table_hits["reference_start_corr"] - table_hits["query_start_corr"],
            table_hits["reference_end_corr"] - (table_hits["query_length"] - table_hits["query_end"]),
        )
        table_hits["end"] = np.where(
            table_hits["reference_strand"] == "plus",
            table_hits["reference_end_corr"] + (table_hits["query_length"] - table_hits["query_end"]),
            table_hits["reference_start_corr"] + table_hits["query_start_corr"],
        )
        return table_hits

    def _remove_overflows(self, bed: pd.DataFrame, file_reference: str) -> pd.DataFrame:
        """
        Removes coordinate overflows in BED format entries relative to reference sequences.

        This function adjusts the start and end coordinates in a BED file to ensure they fit within the boundaries of the reference sequence.
        It calculates and corrects overflows that occur when coordinates extend beyond the reference sequence's length, ensuring the adjusted coordinates are valid.

        :param bed: DataFrame containing BED format data with start and end coordinates.
        :type bed: pd.DataFrame
        :param file_reference: Path to the reference sequence file in FASTA format.
        :type file_reference: str
        :return: DataFrame with corrected start and end coordinates and calculated overflows.
        :rtype: pd.DataFrame
        """
        bed["overflow_start"] = bed["start"].apply(lambda x: -x if x < 0 else 0)
        bed["start"] = bed["start"].apply(lambda x: x if x >= 0 else 0)

        records = SeqIO.parse(file_reference, "fasta")
        regions_length = {record.id: len(record.seq) for record in records}
        bed["len_region"] = bed["chr"].map(regions_length)

        bed["overflow_end"] = bed[["end", "len_region"]].apply(
            lambda x: x["end"] - x["len_region"] if x["end"] > x["len_region"] else 0,
            axis=1,
        )
        bed["end"] = bed[["end", "len_region"]].apply(
            lambda x: x["end"] if x["end"] <= x["len_region"] else x["len_region"],
            axis=1,
        )
        return bed

    def _pad_overflows(self, bed: pd.DataFrame, references: list) -> list:
        """
        Pads sequence references to account for coordinate overflows.

        This function adjusts sequence references by adding padding to account for any overflows that
        occur when sequence coordinates exceed the boundaries of the reference sequence.
        The padding is added at the beginning or end of the sequence depending on the strand direction.

        :param bed: DataFrame containing BED format data, including overflow calculations.
        :type bed: pd.DataFrame
        :param references: List of sequence references extracted from the reference file.
        :type references: list
        :return: List of padded sequence references with overflow accounted for.
        :rtype: list
        """
        references_padded = []
        for reference, overflow_start, overflow_end, strand in zip(
            references, bed["overflow_start"], bed["overflow_end"], bed["strand"]
        ):
            padding_start = "-" * overflow_start
            padding_end = "-" * overflow_end
            if strand == "+":
                reference = padding_start + reference + padding_end
            elif strand == "-":
                reference = padding_end + reference + padding_start

            references_padded.append(reference)
        return references_padded

    def _add_alignment_gaps(
        self, table_hits: pd.DataFrame, queries: list, references: list
    ) -> Tuple[list, list]:
        """
        Aligns query and reference sequences by introducing gaps based on alignment hits.

        This function adjusts query and reference sequences by inserting gaps at specific positions,
        aligning them according to the provided alignment hits.
        The gaps are added to match the alignment patterns derived from the sequence comparisons.

        :param table_hits: DataFrame containing information about the alignment hits, including gap positions.
        :type table_hits: pd.DataFrame
        :param queries: List of query sequences to be aligned.
        :type queries: list
        :param references: List of reference sequences to be aligned.
        :type references: list
        :return: A tuple containing the gapped query and reference sequences.
        :rtype: Tuple[list, list]
        """

        def add_gaps(seq, gaps):
            for gap in gaps:
                seq = seq[:gap] + "-" + seq[gap:]
            return seq

        table_hits["query_gaps"] = (
            table_hits["query_sequence"].apply(lambda x: np.where(np.array(list(x)) == "-")[0])
            + table_hits["query_start"]
            - 1  # blastn has 1-based indices
        )
        table_hits["reference_gaps"] = (
            table_hits["reference_sequence"].apply(lambda x: np.where(np.array(list(x)) == "-")[0])
            + table_hits["query_start"]
            - 1
        )
        gapped_queries = [add_gaps(query, gaps) for query, gaps in zip(queries, table_hits["query_gaps"])]
        gapped_references = [
            add_gaps(reference, gaps) for reference, gaps in zip(references, table_hits["reference_gaps"])
        ]
        return gapped_queries, gapped_references


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

    :param sequence_type: The type of sequence to be used for the filter calculations.
    :type sequence_type: _TYPES_SEQ["oligo", "target"]
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
    :param dir_output: Directory to store output files and temporary data.
    :type dir_output: str
    """

    def __init__(
        self,
        sequence_type: _TYPES_SEQ,
        remove_hits: bool = True,
        search_parameters: dict = None,
        hit_parameters: dict = None,
        names_search_output: list = None,
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
            sequence_type,
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

        :param oligo_database: The OligoDatabase containing the oligonucleotides and their associated attributes.
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

        :param oligo_database: The OligoDatabase containing the oligonucleotides and their associated attributes.
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
                warnings.warn(
                    "Both, 'min_alignment_length' and 'coverage' parameters were provided. Using 'min_alignment_length' parameter."
                )
            min_alignment_length = self.hit_parameters["min_alignment_length"]
        elif "coverage" in self.hit_parameters.keys():
            min_alignment_length = search_results["query_length"] * self.hit_parameters["coverage"] / 100
        else:
            raise KeyError("Please provide either 'coverage' or a 'min_alignment_length' in hit_parameters!")

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
    :param sequence_type: The type of sequence to be used for the filter calculations.
    :type sequence_type: _TYPES_SEQ["oligo", "target"]
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
    :param dir_output: Directory to store output files and temporary data.
    :type dir_output: str
    """

    def __init__(
        self,
        seedregion_start: Union[int, float],
        seedregion_end: Union[int, float],
        sequence_type: _TYPES_SEQ,
        remove_hits: bool = True,
        search_parameters: dict = None,
        hit_parameters: dict = None,
        names_search_output: list = None,
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
            sequence_type,
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

        :param oligo_database: The OligoDatabase containing the oligonucleotides and their associated attributes.
        :type oligo_database: OligoDatabase
        :param search_results: The DataFrame containing the results of the BLAST search.
        :type search_results: pd.DataFrame
        :param region_id: Region ID to process.
        :type region_id: str
        :return: The BLAST search results with added seed region information.
        :rtype: pd.DataFrame
        """
        oligo_attributes_calculator = OligoAttributes()
        oligo_database = oligo_attributes_calculator.calculate_seedregion(
            oligo_database=oligo_database,
            region_ids=region_id,
            start=self.seedregion_start,
            end=self.seedregion_end,
        )

        seedregion = oligo_database.get_oligo_attribute_table(
            ["seedregion_start", "seedregion_end"], region_id
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
    """
    A filter class for BLASTN alignment considering seed region and ligation site proximity.

    The `BlastNSeedregionLigationsiteFilter` class extends the `BlastNSeedregionFilterBase` to focus on evaluating oligonucleotide
    sequences based on their alignment within a specified seed region size, particularly around the ligation site.
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

    :param seedregion_size: The size of the seed region around the ligation site to consider.
    :type seedregion_size: int
    :param sequence_type: The type of sequence to be used for the filter calculations.
    :type sequence_type: _TYPES_SEQ["oligo", "target"]
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
    :param dir_output: Directory to store output files and temporary data.
    :type dir_output: str
    """

    def __init__(
        self,
        seedregion_size: int,
        sequence_type: _TYPES_SEQ,
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
        """Constructor for the BlastNSeedregionLigationsiteFilter class."""
        super().__init__(
            sequence_type,
            remove_hits,
            search_parameters,
            hit_parameters,
            names_search_output,
            filter_name,
            dir_output,
        )
        self.seedregion_size = seedregion_size

    def _add_seed_region_information(
        self, oligo_database: OligoDatabase, search_results: pd.DataFrame, region_id: str
    ) -> pd.DataFrame:
        """
        Adds seed region information to the BLASTN search results based on the ligation site.

        This method calculates the seed region around the ligation site for each oligonucleotide in the OligoDatabase and
        merges this information with the BLASTN search results. The seed region is defined by the `seedregion_size` parameter.

        :param oligo_database: The OligoDatabase containing the oligonucleotides and their associated attributes.
        :type oligo_database: OligoDatabase
        :param search_results: The DataFrame containing the results of the BLASTN search.
        :type search_results: pd.DataFrame
        :param region_id: Region ID to process.
        :type region_id: str
        :return: The BLAST search results with added seed region information.
        :rtype: pd.DataFrame
        """
        oligo_attributes_calculator = OligoAttributes()
        oligo_database = oligo_attributes_calculator.calculate_seedregion_ligationsite(
            oligo_database=oligo_database, region_ids=region_id, seedregion_size=self.seedregion_size
        )

        seedregion = oligo_database.get_oligo_attribute_table(
            ["seedregion_start", "seedregion_end"], region_id
        )
        search_results = pd.merge(
            left=search_results,
            right=seedregion,
            left_on="query",
            right_on="oligo_id",
            how="left",
        )
        return search_results
