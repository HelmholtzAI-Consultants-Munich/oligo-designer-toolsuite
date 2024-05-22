############################################
# imports
############################################

import os
import pandas as pd

from joblib import Parallel, delayed
from oligo_designer_toolsuite_ai_filters.api import APIHybridizationProbability

from oligo_designer_toolsuite._constants import _TYPES_SEQ
from oligo_designer_toolsuite.database import OligoDatabase, ReferenceDatabase

from oligo_designer_toolsuite.oligo_specificity_filter import (
    SpecificityFilterBase,
    AlignmentSpecificityFilter,
    BlastNFilter,
    BlastNSeedregionFilter,
    BlastNSeedregionLigationsiteFilter,
)

############################################
# Oligo AI Filter Classes
############################################


class HybridizationProbabilityFilter(SpecificityFilterBase):
    """
    This class filters oligos based on the result of a user-specified aligment method and the hibridization probability of the oligos with the off-target sites.
    In particiular, the method considers as hits only the off-target sites with an hybridization probability wiht the oligos higher than a user-defined threshold.

    To define the hybridization probability we cosider an experiment where the oligo sequence, the on-target site and the off-target site are intruduced with the
    same concentration. Next, the hybridization probability is the ration of the final concentration of DNA complex composed by the the oligo and the off-target site,
    and all the DNA complexes that contain the oligo (oligo + off-target and oligo + on-target):

    $$C_{oligo + off-t} /C_{oligo + off-t}  + C_{oligo + on-t}$$.

    In simpler terms, the hybridization probability the frequency with which the oligo binds to the off-target site, compared to the number of succesfull bindinbgs of the oligo.

    :param alignment_method: The alignment specificity filter to use.
    :type alignment_method: AlignmentSpecificityFilter
    :param threshold: The threshold below which the oligos are filtered.
    :type threshold: float, optional
    :param ai_filter_path: The path to the machine learning model used to filter the oligos, if None the pretrained model provided will be used, defaults to None
    :type ai_filter_path: str, optional
    :param filter_name: Subdirectory path for the output, i.e. <dir_output>/<filter_name>, defaults to "ai_filter".
    :type filter_name: str, optional
    :param dir_output: Directory for saving intermediate files, defaults to "output"
    :type dir_output: str, optional
    """

    def __init__(
        self,
        alignment_method: AlignmentSpecificityFilter,
        threshold: float,
        ai_filter_path: str = None,
        filter_name: str = "ai_filter",
        dir_output: str = "output",
    ) -> None:
        """Constructor for the HybridizationProbabilityFilter class."""
        super().__init__(filter_name, dir_output)

        self.alignment_method = alignment_method
        self.overwrite_output_format()

        # instatiate ai model
        self.threshold = threshold
        self.model = APIHybridizationProbability(ai_filter_path=ai_filter_path)

    def apply(
        self,
        sequence_type: _TYPES_SEQ,
        oligo_database: OligoDatabase,
        reference_database: ReferenceDatabase,
        n_jobs: int = 1,
    ):
        """
        Applies the alignment-based specificity filter and the  machine leanrning based filter to an oligonucleotide database.

        :param sequence_type: The type of sequences being filtered, must be one of the predefined sequence types.
        :type sequence_type: _TYPES_SEQ
        :param database: The oligo database to which the filter will be applied.
        :type database: OligoDatabase
        :param reference_database: The reference database to compare against for specificity.
        :type reference_database: ReferenceDatabase
        :param n_jobs: The number of parallel jobs to run.
        :type n_jobs: int
        :return: The filtered oligo database with sequences having significant hits removed.
        :rtype: OligoDatabase
        """
        region_ids = list(oligo_database.database.keys())
        table_hits = self.alignment_method._get_table_hits(
            sequence_type=sequence_type,
            oligo_database=oligo_database,
            reference_database=reference_database,
            region_ids=region_ids,
            n_jobs=n_jobs,
        )
        # filter the table hits
        file_reference = reference_database.write_database_to_fasta(
            filename=f"db_reference_{self.filter_name}"
        )  # defined in advance to avoid writing the file multiple times
        table_hits = Parallel(n_jobs=n_jobs)(
            delayed(self._filter_table_hits)(
                sequence_type=sequence_type,
                table_hits=table_hits_region,
                oligo_database=oligo_database,
                file_reference=file_reference,
                region_id=region_id,
            )
            for table_hits_region, region_id in zip(table_hits, region_ids)
        )

        # filter the oligo database
        for region_id, table_hits_region in zip(region_ids, table_hits):
            oligos_with_hits_region = table_hits_region["query"].unique()
            database_region_filtered = self._filter_hits_from_database(
                database_region=oligo_database.database[region_id],
                oligos_with_hits=oligos_with_hits_region,
            )
            oligo_database.database[region_id] = database_region_filtered

        os.remove(file_reference)
        os.remove(file_reference + ".fai")
        return oligo_database

    def _filter_table_hits(
        self,
        sequence_type: _TYPES_SEQ,
        table_hits: pd.DataFrame,
        oligo_database: OligoDatabase,
        file_reference: str,
        region_id: str,
    ) -> pd.DataFrame:
        """Filters the hits from a search operation using Machine Learning models. The Hits that recieve a score form the machine learning model lower that the given threshold are filtered out.

        :param sequence_type: The type of sequences being filtered, must be one of the predefined sequence types.
        :type sequence_type: _TYPES_SEQ
        :param table_hits: Dataframe containing the strue hits of the search results.
        :type table_hits: pd.DataFrame
        :param file_reference: Path to the fasta file used as reference for the search.
        :type file_reference: str
        :param oligo_database: The oligo database to which the filter will be applied.
        :type oligo_database: OligoDatabase
        :param region_id: The identifier for the region within the database to filter.
        :type region_id: str
        :return: Dataframe containing the filtered true hits
        :rtype: pd.DataFrame
        """

        # check if there are any oligos to filter
        if len(table_hits) == 0:
            return table_hits

        # generate the references and queries sequences
        references = self.alignment_method._get_references(table_hits, file_reference, region_id)
        queries = self.alignment_method._get_queries(sequence_type, table_hits, oligo_database, region_id)

        # align the references and queries by adding gaps
        gapped_queries, gapped_references = self.alignment_method._add_alignement_gaps(
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

    def overwrite_output_format(self):
        """
        The output format of the alignment method is overwritten to be compatible with the hybridization probability filter.
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
            self.alignment_method.search_parameters["outfmt"] = (
                "6 qseqid sseqid length qstart qend qlen qseq sstart send sseq sstrand"
            )
