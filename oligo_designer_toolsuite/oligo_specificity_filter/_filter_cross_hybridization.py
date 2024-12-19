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
    def __init__(
        self,
        policy: FilterPolicyBase,
        alignment_method: AlignmentSpecificityFilter,
        database_name_reference: str = "db_reference",
        filter_name: str = "crosshybridization_filter",
        dir_output: str = "output",
    ) -> None:
        """Constructor for the CrossHybridizationFilter class."""
        super().__init__(filter_name, dir_output)

        self.database_name_reference = database_name_reference
        self.dir_output_reference = dir_output

        self.policy = policy
        self.alignment_method = alignment_method

    def apply(
        self,
        oligo_database: OligoDatabase,
        reference_database: ReferenceDatabase,
        sequence_type: _TYPES_SEQ,
        n_jobs: int = 1,
    ) -> OligoDatabase:
        region_ids = list(oligo_database.database.keys())

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

        for region_id in region_ids:
            self._filter_hits_from_database(
                oligo_database=oligo_database,
                region_id=region_id,
                oligos_with_hits=oligos_with_hits[region_id],
            )
        return oligo_database

    def _create_reference_database(
        self, sequence_type: _TYPES_SEQ, oligo_database: OligoDatabase
    ) -> ReferenceDatabase:
        file_reference = oligo_database.write_database_to_fasta(
            filename=f"db_reference_{self.filter_name}",
            save_description=False,
            region_ids=None,
            sequence_type=sequence_type,
        )
        reference_database = ReferenceDatabase(
            database_name=self.database_name_reference, dir_output=self.dir_output_reference
        )
        reference_database.load_database_from_fasta(files_fasta=file_reference, database_overwrite=True)

        os.remove(file_reference)

        return reference_database
