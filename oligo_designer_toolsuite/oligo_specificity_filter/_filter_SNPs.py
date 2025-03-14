############################################
# Imports
############################################

import os
import pandas as pd

from oligo_designer_toolsuite._constants import _TYPES_SEQ, SEPARATOR_OLIGO_ID
from oligo_designer_toolsuite.database import OligoDatabase, ReferenceDatabase
from oligo_designer_toolsuite.oligo_specificity_filter import (
    SpecificityFilterBase,
)

from oligo_designer_toolsuite.utils import get_intersection

############################################
# Oligo Crosshybridization Filter Class
############################################


class SNPFilter(SpecificityFilterBase):

    def __init__(
        self,
        filter_name: str = "SNP_filter",
        dir_output: str = "output",
    ) -> None:
        """Constructor for the SNPFilter class."""
        super().__init__(filter_name, dir_output)

    def apply(
        self,
        oligo_database: OligoDatabase,
        reference_database: ReferenceDatabase,
        sequence_type: _TYPES_SEQ,
        n_jobs: int = 1,
    ) -> OligoDatabase:
        """
        Applies the cross-hybridization filter to the OligoDatabase, removing sequences that may cross-hybridize with other sequences in the OligoDatabase.

        This function compares oligonucleotides in the OligoDatabase with reference sequences consisting of all sequences in the OligoDatabase using the specified alignment method.
        Based on the results, it filters out oligonucleotides that meet the cross-hybridization criteria defined by the policy.

        :param sequence_type: The type of sequence to be used for the filter calculations.
        :type sequence_type: _TYPES_SEQ["oligo", "target"]
        :param oligo_database: The OligoDatabase containing the oligonucleotides and their associated attributes.
        :type oligo_database: OligoDatabase
        :param reference_database: The ReferenceDatabase used for alignment (not utilized in this filter).
        :type reference_database: ReferenceDatabase
        :param n_jobs: The number of parallel jobs to use for processing.
        :type n_jobs: int
        :return: The filtered OligoDatabase with potential cross-hybridizing sequences removed.
        :rtype: OligoDatabase
        """
        region_ids = list(oligo_database.database.keys())

        file_vcf = reference_database.write_database_to_vcf(filename=f"db_reference_{self.filter_name}")

        for region_id in region_ids:

            file_bed = oligo_database.write_database_to_bed(
                filename=f"{region_id}.bed", dir_output=self.dir_output, region_ids=region_id
            )
            file_bed_out = os.path.join(self.dir_output, f"SNP_results_{region_id}.txt")

            get_intersection(file_bed, file_vcf, file_bed_out)

            search_results = pd.read_csv(
                filepath_or_buffer=file_bed_out,
                header=None,
                sep="\t",
                low_memory=False,
                engine="c",
                usecols=[3, 8],
                names=["oligo_id", "snp_id"],
            )

            search_results["region_id"] = search_results["oligo_id"].str.split(SEPARATOR_OLIGO_ID).str[0]

            search_results = search_results.groupby("oligo_id", as_index=False).agg(
                {"snp_id": list, "region_id": "first"}
            )

        ## -> add new function which flags instead of filters hit
        for region_id in region_ids:
            self._filter_hits_from_database(
                oligo_database=oligo_database,
                region_id=region_id,
                oligos_with_hits=oligos_with_hits[region_id],
            )
        return oligo_database
