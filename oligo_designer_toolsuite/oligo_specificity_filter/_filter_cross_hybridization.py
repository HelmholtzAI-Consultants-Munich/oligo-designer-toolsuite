############################################
# imports
############################################

import os
from abc import ABC, abstractmethod
from pathlib import Path

import networkx as nx

from .._constants import _TYPES_SEQ, SEPARATOR_OLIGO_ID
from ..database import OligoDatabase, ReferenceDatabase
from ._filter_base import AlignmentSpecificityFilter, SpecificityFilterBase
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
    :param dir_output: Directory for saving output files related to cross-hybridization filtering.
    :type dir_output: str
    """

    def __init__(
        self,
        policy: FilterPolicyBase,
        alignment_method: AlignmentSpecificityFilter,
        dir_output: str = "output",
    ):
        """Constructor for the CrossHybridizationFilter class."""
        self.policy = policy
        self.alignment_method = alignment_method

        self.dir_output = os.path.join(dir_output, "crosshybridization")
        Path(self.dir_output).mkdir(parents=True, exist_ok=True)

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
        :param database: The database of oligonucleotides to be filtered.
        :type database: OligoDatabase
        :param reference_database: The reference database to compare against for specificity.
        :type reference_database: ReferenceDatabase
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
            filename=f"oligo_database_crosshybridization_with_{sequence_type}",
            region_ids=None,
            sequence_type=sequence_type,
        )
        reference_database = ReferenceDatabase(dir_output=self.dir_output)
        reference_database.load_sequences_from_fasta(files_fasta=file_reference, database_overwrite=True)

        os.remove(file_reference)

        return reference_database
