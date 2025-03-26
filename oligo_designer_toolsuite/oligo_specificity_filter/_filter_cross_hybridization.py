############################################
# Imports
############################################

import os
from typing import List


from oligo_designer_toolsuite._constants import _TYPES_SEQ
from oligo_designer_toolsuite.database import OligoDatabase, ReferenceDatabase
from oligo_designer_toolsuite.oligo_specificity_filter import (
    SpecificityFilterBase,
    SpecificityFilterAlignment,
)

from ._policies import FilterPolicyBase
from oligo_designer_toolsuite.utils import FastaParser


############################################
# Oligo Crosshybridization Filter Class
############################################


class CrossHybridizationFilter(SpecificityFilterBase):
    """
    A filter class designed to identify and filter out oligonucleotides that are likely to cross-hybridize with other sequences in the OligoDatabase.

    The `CrossHybridizationFilter` class uses a specified alignment method to compare oligonucleotides against a ReferenceDatabase,
    i.e. in the case of cross-hybridization, the ReferenceDatabase consists of all sequences in the OligoDatabase.
    The filter removes sequences that show potential cross-hybridization, based on the provided filtering policy.

    :param policy: The policy that defines how to handle oligonucleotides that meet or violate the cross-hybridization criteria.
    :type policy: FilterPolicyBase
    :param alignment_method: The alignment method used to identify potential cross-hybridization with reference sequences (i.e. OligoDatabase sequences).
    :type alignment_method: SpecificityFilterAlignment
    :param filter_name: Name of the filter for identification purposes.
    :type filter_name: str
    :param dir_output: Directory to store output files and temporary data.
    :type dir_output: str
    """

    def __init__(
        self,
        policy: FilterPolicyBase,
        alignment_method: SpecificityFilterAlignment,
        sequence_type_reference: _TYPES_SEQ = None,
        filter_name: str = "cross_hybridization_filter",
        dir_output: str = "output",
    ) -> None:
        """Constructor for the CrossHybridizationFilter class."""
        super().__init__(filter_name, dir_output)

        self.policy = policy
        self.alignment_method = alignment_method
        self.sequence_type_reference = sequence_type_reference

    def set_reference_database(self, oligo_database: OligoDatabase) -> None:
        """
        Creates a ReferenceDatabase for cross-hybridization filtering based on the sequences in the provided OligoDatabase.

        This function generates a FASTA file from the OligoDatabase, which is then loaded into a new ReferenceDatabase object.
        The ReferenceDatabase is used to identify potential cross-hybridization with the oligo sequences.

        :param oligo_database: The OligoDatabase containing the oligonucleotides and their associated attributes.
        :type oligo_database: OligoDatabase
        :param sequence_type: The type of sequence to be used for filter calculations.
        :type sequence_type: _TYPES_SEQ["oligo", "target"]
        :return: The generated ReferenceDatabase.
        :rtype: ReferenceDatabase
        """
        file_reference = oligo_database.write_database_to_fasta(
            filename=f"db_oligos_{self.filter_name}_{self.sequence_type_reference}",
            dir_output=self.dir_output,
            save_description=True,
            region_ids=None,
            sequence_type=self.sequence_type_reference,
        )

        reference_database = ReferenceDatabase(
            database_name=f"db_reference_{self.filter_name}_{self.sequence_type_reference}",
            dir_output=self.dir_output,
        )
        reference_database.load_database_from_file(
            files=file_reference, file_type="fasta", database_overwrite=True
        )
        os.remove(file_reference)

        return reference_database

    def apply(
        self,
        oligo_database: OligoDatabase,
        sequence_type: _TYPES_SEQ,
        n_jobs: int = 1,
    ) -> OligoDatabase:
        """
        Applies the cross-hybridization filter to the OligoDatabase, removing sequences that may cross-hybridize with other sequences in the OligoDatabase.

        This function compares oligonucleotides in the OligoDatabase with reference sequences consisting of all sequences in the OligoDatabase using the specified alignment method.
        Based on the results, it filters out oligonucleotides that meet the cross-hybridization criteria defined by the policy.

        :param oligo_database: The OligoDatabase containing the oligonucleotides and their associated attributes.
        :type oligo_database: OligoDatabase
        :param sequence_type: The type of sequence to be used for filter calculations.
        :type sequence_type: _TYPES_SEQ["oligo", "target"]
        :param n_jobs: The number of parallel jobs to use for processing.
        :type n_jobs: int
        :return: The filtered OligoDatabase with potential cross-hybridizing sequences removed.
        :rtype: OligoDatabase
        """
        if not self.sequence_type_reference:
            self.sequence_type_reference = sequence_type

        print(self.dir_output)

        region_ids = list(oligo_database.database.keys())

        reference_database = self.set_reference_database(oligo_database=oligo_database)
        self.alignment_method.set_reference_database(reference_database=reference_database)

        print(self.alignment_method.reference_database.database_file)

        oligo_pair_hits = self.alignment_method.get_oligo_pair_hits(
            oligo_database=oligo_database,
            sequence_type=sequence_type,
            n_jobs=n_jobs,
        )
        oligos_with_hits = self.policy.apply(oligo_pair_hits=oligo_pair_hits, oligo_database=oligo_database)

        self._filter_hits_from_database(
            oligo_database=oligo_database,
            region_ids=region_ids,
            oligos_with_hits=oligos_with_hits,
        )
        # need to reset reference database, oterwise it causes problems when using the filter twice
        self.alignment_method.set_reference_database(reference_database=None)

        return oligo_database
