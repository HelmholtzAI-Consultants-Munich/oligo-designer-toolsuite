############################################
# Imports
############################################

import os
from typing import get_args


from oligo_designer_toolsuite._constants import _TYPES_SEQ
from oligo_designer_toolsuite.database import OligoDatabase, ReferenceDatabase
from oligo_designer_toolsuite.oligo_specificity_filter import (
    SpecificityFilterBase,
    SpecificityFilterAlignment,
)

from ._policies import FilterPolicyBase

############################################
# Oligo Crosshybridization Filter Class
############################################


class CrossHybridizationFilter(SpecificityFilterBase):
    """
    A filter class designed to identify and filter out oligonucleotides that are likely to cross-hybridize with other sequences in the OligoDatabase.

    The `CrossHybridizationFilter` class uses a specified alignment method to compare oligonucleotides against a ReferenceDatabase,
    i.e. in the case of cross-hybridization, the ReferenceDatabase consists of all sequences in the OligoDatabase.
    The filter removes sequences that show potential cross-hybridization, based on the provided filtering policy.

    :param sequence_type: The type of sequence to be used for the filter calculations.
    :type sequence_type: _TYPES_SEQ["oligo", "target"]
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
        sequence_type: _TYPES_SEQ,
        policy: FilterPolicyBase,
        alignment_method: SpecificityFilterAlignment,
        filter_name: str = "cross_hybridization_filter",
        dir_output: str = "output",
    ) -> None:
        """Constructor for the CrossHybridizationFilter class."""
        super().__init__(filter_name, dir_output)

        options = get_args(_TYPES_SEQ)
        assert (
            sequence_type in options
        ), f"Sequence type not supported! '{sequence_type}' is not in {options}."

        self.sequence_type = sequence_type
        self.policy = policy
        self.alignment_method = alignment_method
        self.alignment_method.sequence_type = sequence_type  # make sure sequence_types match
        self.reference_database = None

    def set_reference_database(self, oligo_database: OligoDatabase) -> None:
        """
        Creates a ReferenceDatabase for cross-hybridization filtering based on the sequences in the provided OligoDatabase.

        This function generates a FASTA file from the OligoDatabase, which is then loaded into a new ReferenceDatabase object.
        The ReferenceDatabase is used to identify potential cross-hybridization with the oligo sequences.

        :param oligo_database: The OligoDatabase containing the oligonucleotides and their associated attributes.
        :type oligo_database: OligoDatabase
        """
        file_reference = oligo_database.write_database_to_fasta(
            filename=f"db_oligos_{self.filter_name}",
            dir_output=self.dir_output,
            save_description=False,
            region_ids=None,
            sequence_type=self.sequence_type,
        )
        self.reference_database = ReferenceDatabase(
            database_name=f"db_reference_{self.filter_name}", dir_output=self.dir_output
        )
        self.reference_database.load_database_from_file(
            files=file_reference, file_type="fasta", database_overwrite=True
        )
        self.alignment_method.set_reference_database(self.reference_database)
        os.remove(file_reference)

    def apply(
        self,
        oligo_database: OligoDatabase,
        n_jobs: int = 1,
    ) -> OligoDatabase:
        """
        Applies the cross-hybridization filter to the OligoDatabase, removing sequences that may cross-hybridize with other sequences in the OligoDatabase.

        This function compares oligonucleotides in the OligoDatabase with reference sequences consisting of all sequences in the OligoDatabase using the specified alignment method.
        Based on the results, it filters out oligonucleotides that meet the cross-hybridization criteria defined by the policy.

        :param oligo_database: The OligoDatabase containing the oligonucleotides and their associated attributes.
        :type oligo_database: OligoDatabase
        :param n_jobs: The number of parallel jobs to use for processing.
        :type n_jobs: int
        :return: The filtered OligoDatabase with potential cross-hybridizing sequences removed.
        :rtype: OligoDatabase
        """
        region_ids = list(oligo_database.database.keys())

        self.set_reference_database(oligo_database=oligo_database)
        oligo_pair_hits = self.alignment_method.get_oligo_pair_hits(
            oligo_database=oligo_database,
            n_jobs=n_jobs,
        )
        oligos_with_hits = self.policy.apply(oligo_pair_hits=oligo_pair_hits, oligo_database=oligo_database)

        self._filter_hits_from_database(
            oligo_database=oligo_database,
            region_ids=region_ids,
            oligos_with_hits=oligos_with_hits,
        )

        return oligo_database
