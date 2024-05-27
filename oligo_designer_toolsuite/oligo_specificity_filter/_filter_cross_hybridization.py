############################################
# Imports
############################################

import os

from oligo_designer_toolsuite._constants import _TYPES_SEQ
from oligo_designer_toolsuite.database import OligoDatabase, ReferenceDatabase
from oligo_designer_toolsuite.oligo_specificity_filter import (
    AlignmentSpecificityFilter,
    SpecificityFilterBase,
)

from ._policies import FilterPolicyBase

############################################
# Oligo Crosshybridization Filter Class
############################################


class CrossHybridizationFilter(SpecificityFilterBase):
    """A class that implements a filter to minimize cross-hybridization within an oligonucleotide database by applying a specified policy.
    This class integrates a cross-hybridization policy and a specificity filter to identify and mitigate potential
    cross-hybridization events between oligonucleotides.

    :param policy: The filter policy to apply for minimizing cross-hybridization.
    :type policy: FilterPolicyBase
    :param alignment_method: The alignment specificity filter used to identify potential cross-hybridization events.
    :type alignment_method: AlignmentSpecificityFilter
    :param database_name_reference: Subdirectory path for the reference database, i.e., <dir_output>/<database_name_reference>, defaults to "db_reference".
    :type database_name_reference: str, optional
    :param filter_name: Subdirectory path for the output, i.e., <dir_output>/<filter_name>, defaults to "crosshybridization_filter".
    :type filter_name: str, optional
    :param dir_output: Directory for saving intermediate files, defaults to "output".
    :type dir_output: str, optional
    """

    def __init__(
        self,
        policy: FilterPolicyBase,
        alignment_method: AlignmentSpecificityFilter,
        database_name_reference: str = "db_reference",
        filter_name: str = "crosshybridization_filter",
        dir_output: str = "output",
    ):
        """Constructor for the CrossHybridizationFilter class."""
        super().__init__(filter_name, dir_output)

        self.database_name_reference = database_name_reference
        self.dir_output_reference = dir_output

        self.policy = policy
        self.alignment_method = alignment_method

    def apply(
        self,
        sequence_type: _TYPES_SEQ,
        oligo_database: OligoDatabase,
        reference_database: ReferenceDatabase = None,
        n_jobs: int = 1,
    ):
        """Applies the cross-hybridization filter to an oligonucleotide database.

        :param sequence_type: The type of sequences being filtered, must be one of the predefined sequence types.
        :type sequence_type: _TYPES_SEQ
        :param oligo_database: The database of oligonucleotides to be filtered.
        :type oligo_database: OligoDatabase
        :param reference_database: The reference database to compare against for specificity.
        :type reference_database: ReferenceDatabase, optional
        :param n_jobs: The number of parallel jobs to run.
        :type n_jobs: int
        :return: The oligo database with cross-hybridization minimized according to the policy.
        :rtype: OligoDatabase
        """
        regions = list(oligo_database.database.keys())

        reference_database = self._create_reference_database(
            sequence_type=sequence_type, oligo_database=oligo_database
        )
        oligo_pair_hits = self.alignment_method.get_oligo_pair_hits(
            sequence_type=sequence_type,
            oligo_database=oligo_database,
            reference_database=reference_database,
            n_jobs=n_jobs,
        )
        oligos_with_hits = self.policy.apply(oligo_pair_hits=oligo_pair_hits, oligo_database=oligo_database)

        for region in regions:
            database_region_filtered = self._filter_hits_from_database(
                database_region=oligo_database.database[region],
                oligos_with_hits=oligos_with_hits[region],
            )
            oligo_database.database[region] = database_region_filtered
        return oligo_database

    def _create_reference_database(self, sequence_type: _TYPES_SEQ, oligo_database: OligoDatabase):
        """Creates a reference database from the given oligo database. This involves writing the oligo database to a FASTA file
        and loading it into a new ReferenceDatabase instance.

        :param sequence_type: The type of sequences being filtered, must be one of the predefined sequence types.
        :type sequence_type: _TYPES_SEQ
        :param oligo_database: The oligo database from which the reference database is created.
        :type oligo_database: OligoDatabase
        :return: A ReferenceDatabase instance populated with sequences from the oligo database.
        :rtype: ReferenceDatabase
        """
        file_reference = oligo_database.write_database_to_fasta(
            filename=f"db_reference_{self.filter_name}",
            region_ids=None,
            sequence_type=sequence_type,
        )
        reference_database = ReferenceDatabase(
            database_name=self.database_name_reference, dir_output=self.dir_output_reference
        )
        reference_database.load_sequences_from_fasta(files_fasta=file_reference, database_overwrite=True)

        os.remove(file_reference)

        return reference_database
