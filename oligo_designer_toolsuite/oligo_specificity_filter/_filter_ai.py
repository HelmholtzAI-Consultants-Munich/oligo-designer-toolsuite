############################################
# imports
############################################

import os

import pandas as pd
from joblib import Parallel, delayed
from joblib_progress import joblib_progress
from oligo_designer_toolsuite_ai_filters.api import APIHybridizationProbability

from oligo_designer_toolsuite._constants import _TYPES_SEQ
from oligo_designer_toolsuite.database import OligoDatabase, ReferenceDatabase
from oligo_designer_toolsuite.oligo_specificity_filter import (
    AlignmentSpecificityFilter,
    BlastNFilter,
    BlastNSeedregionFilter,
    BlastNSeedregionLigationsiteFilter,
    SpecificityFilterBase,
)

############################################
# Oligo AI Filter Classes
############################################


class HybridizationProbabilityFilter(SpecificityFilterBase):
    """
    A filter for assessing hybridization probabilities of oligonucleotides using AI models.

    The `HybridizationProbabilityFilter` class evaluates the likelihood of hybridization between oligonucleotides and potential off-target sequences based on a threshold.
    It utilizes an AI model to predict hybridization probabilities and filters out oligonucleotides that exceed the specified threshold.

    To define the hybridization probability we consider an experiment where the oligo sequence, the on-target site, and the off-target site are introduced with the
    same concentration. Next, the hybridization probability is the ratio of the final concentration of the DNA complex composed of the oligo and the off-target site,
    and all the DNA complexes that contain the oligo (oligo + off-target and oligo + on-target):

    $$C_{oligo + off-t} /C_{oligo + off-t} + C_{oligo + on-t}$$.

    In simpler terms, the hybridization probability is the frequency with which the oligo binds to the off-target site, compared to the number of successful bindings of the oligo.

    :param alignment_method: The alignment method used to identify potential off-targets with reference sequences.
    :type alignment_method: AlignmentSpecificityFilter
    :param threshold: The probability threshold above which sequences will be filtered out.
    :type threshold: float
    :param ai_filter_path: Path to the AI model used for hybridization probability predictions. If None, the default path will be used.
    :type ai_filter_path: str, optional
    :param filter_name: Name of the filter for identification purposes.
    :type filter_name: str
    :param dir_output: Directory to store output files and temporary data.
    :type dir_output: str
    """

    def __init__(
        self,
        alignment_method: AlignmentSpecificityFilter,
        threshold: float,
        ai_filter_path: str = None,
        filter_name: str = "hybridization_probability_filter",
        dir_output: str = "output",
    ) -> None:
        """Constructor for the HybridizationProbabilityFilter class."""
        super().__init__(filter_name, dir_output)

        self.alignment_method = alignment_method
        self.overwrite_output_format()

        # instantiate ai model
        self.threshold = threshold
        self.model = APIHybridizationProbability(ai_filter_path=ai_filter_path)

    def apply(
        self,
        oligo_database: OligoDatabase,
        reference_database: ReferenceDatabase,
        sequence_type: _TYPES_SEQ,
        n_jobs: int = 1,
    ) -> OligoDatabase:
        """
        Applies the hybridization probability filter to the OligoDatabase.

        This function aligns the oligonucleotides in the given OligoDatabase against a reference database and then
        uses an AI model that predicts hybridization probabilities. The filter removes oligonucleotides that are
        likely to hybridize with off-target sequences, based on a predefined threshold.

        :param oligo_database: The OligoDatabase containing the oligonucleotides and their associated attributes.
        :type oligo_database: OligoDatabase
        :param reference_database: The ReferenceDatabase used for alignment.
        :type reference_database: ReferenceDatabase
        :param sequence_type: The type of sequence to be used for the filter calculations.
        :type sequence_type: _TYPES_SEQ["oligo", "target"]
        :param n_jobs: The number of parallel jobs to use for processing.
        :type n_jobs: int
        :return: The filtered OligoDatabase with off-targets removed.
        :rtype: OligoDatabase
        """
        # When applying the filter we don't want to consider hits within the same region
        consider_hits_from_input_region = False

        # Create reference and index file for search
        # Defined in advance to avoid writing the file multiple times
        file_reference = reference_database.write_database_to_fasta(
            filename=f"db_reference_{self.filter_name}"
        )
        file_index = self.alignment_method._create_index(file_reference=file_reference, n_jobs=n_jobs)

        # run search in parallel for each region
        region_ids = list(oligo_database.database.keys())
        name = " ".join(string.capitalize() for string in self.filter_name.split("_"))
        with joblib_progress(description=f"Specificity Filter: {name}", total=len(region_ids)):
            Parallel(n_jobs=n_jobs, prefer="threads", require="sharedmem")(
                delayed(self._apply_region)(
                    sequence_type=sequence_type,
                    oligo_database=oligo_database,
                    file_reference=file_reference,
                    file_index=file_index,
                    region_id=region_id,
                    consider_hits_from_input_region=consider_hits_from_input_region,
                )
                for region_id in region_ids
            )

        os.remove(file_reference)
        os.remove(file_reference + ".fai")
        self.alignment_method._remove_index(file_index)

        return oligo_database

    def _apply_region(
        self,
        sequence_type: _TYPES_SEQ,
        oligo_database: OligoDatabase,
        file_reference: str,
        file_index: str,
        region_id: str,
        consider_hits_from_input_region: bool,
    ) -> None:
        """
        Applies the hybridization probability filter to a specific region in the OligoDatabase.

        This function processes a specific region of the OligoDatabase by filtering sequences based on their
        likelihood to hybridize with off-target sequences in a reference database. The filtering process uses
        an alignment method and AI model to evaluate the sequences, removing those that meet the criteria for potential cross-hybridization.

        :param sequence_type: The type of sequence to be used for filter calculations.
        :type sequence_type: _TYPES_SEQ["oligo", "target"]
        :param oligo_database: The OligoDatabase containing the oligonucleotides and their associated attributes.
        :type oligo_database: OligoDatabase
        :param file_reference: The path to the reference FASTA file used for the alignment.
        :type file_reference: str
        :param file_index: The index file for the ReferenceDatabase.
        :type file_index: str
        :param region_id: Region ID to process.
        :type region_id: str
        :param consider_hits_from_input_region: Whether to consider hits from the input region in the filtering process.
        :type consider_hits_from_input_region: bool
        """
        table_hits_region = self.alignment_method._run_filter(
            sequence_type=sequence_type,
            region_id=region_id,
            oligo_database=oligo_database,
            file_index=file_index,
            consider_hits_from_input_region=consider_hits_from_input_region,
        )

        table_hits_region = self._filter_table_hits(
            sequence_type=sequence_type,
            table_hits=table_hits_region,
            oligo_database=oligo_database,
            file_reference=file_reference,
            region_id=region_id,
        )

        oligos_with_hits_region = table_hits_region["query"].unique()
        self._filter_hits_from_database(
            oligo_database=oligo_database,
            region_id=region_id,
            oligos_with_hits=oligos_with_hits_region,
        )

    def _filter_table_hits(
        self,
        sequence_type: _TYPES_SEQ,
        table_hits: pd.DataFrame,
        oligo_database: OligoDatabase,
        file_reference: str,
        region_id: str,
    ) -> pd.DataFrame:
        """
        Filters the hits from the alignment table based on hybridization probability.

        This function processes the alignment results by filtering out hits with a hybridization probability
        below a specified threshold. It generates aligned sequences (queries and references) with gaps,
        predicts their hybridization probability, and retains only those hits that exceed the threshold.

        :param sequence_type: The type of sequence to be used for filter calculations.
        :type sequence_type: _TYPES_SEQ["oligo", "target"]
        :param table_hits: A DataFrame containing the alignment results.
        :type table_hits: pd.DataFrame
        :param oligo_database: The OligoDatabase containing the oligonucleotides and their associated attributes.
        :type oligo_database: OligoDatabase
        :param file_reference: The path to the reference FASTA file used for the alignment.
        :type file_reference: str
        :param region_id: Region ID to process.
        :type region_id: str
        :return: A filtered DataFrame containing only the hits with hybridization probabilities above the threshold.
        :rtype: pd.DataFrame
        """

        # check if there are any oligos to filter
        if len(table_hits) == 0:
            return table_hits
        # generate the references and queries sequences
        references = self.alignment_method._get_references(table_hits, file_reference, region_id)
        queries = self.alignment_method._get_queries(
            oligo_database=oligo_database,
            table_hits=table_hits,
            sequence_type=sequence_type,
            region_id=region_id,
        )
        # align the references and queries by adding gaps
        gapped_queries, gapped_references = self.alignment_method._add_alignment_gaps(
            table_hits=table_hits, queries=queries, references=references
        )

        # predict the scores for each hit
        predictions = self.model.predict(
            queries=queries,
            gapped_queries=gapped_queries,
            references=references,
            gapped_references=gapped_references,
        )

        # filter the database, keep only the oligos above the threshold
        table_hits = table_hits[predictions >= self.threshold]
        return table_hits

    def overwrite_output_format(self) -> None:
        """
        Overwrites the output format for the alignment method if it is a BlastN method.

        This function modifies the output format of the alignment method to include specific fields required for subsequent processing.
        It checks if the alignment method is an instance of any of the BlastN-related classes, and if so, updates
        the `names_search_output` and `search_parameters["outfmt"]` to include detailed alignment information.

        :param alignment_method: The alignment method object whose output format may need to be overwritten.
        :type alignment_method: AlignmentSpecificityFilter
        """
        # if the alignment method is a  Blastn method overwrite the
        if type(self.alignment_method) in [
            BlastNFilter,
            BlastNSeedregionFilter,
            BlastNSeedregionLigationsiteFilter,
        ]:
            self.alignment_method.names_search_output = [
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
            self.alignment_method.search_parameters[
                "outfmt"
            ] = "6 qseqid sseqid length qstart qend qlen qseq sstart send sseq sstrand"
