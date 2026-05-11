############################################
# Imports
############################################

import os

from joblib import Parallel, delayed
from joblib_progress import joblib_progress

from oligo_designer_toolsuite._exceptions import ConfigurationError, DatabaseError
from oligo_designer_toolsuite.database import OligoDatabase
from oligo_designer_toolsuite.oligo_specificity_filter import ReferenceSpecificityFilter
from oligo_designer_toolsuite.utils import get_intersection

############################################
# Oligo Variants Filter Class
############################################


class VariantsFilter(ReferenceSpecificityFilter):
    """
    A specificity filter that removes or flags oligos overlapping known variants (e.g., SNPs)
    by intersecting oligo positions with a reference database using BEDTools.

    This class extends the `ReferenceSpecificityFilter` and implements a BED-based
    filtering pipeline. Depending on the `remove_hits` flag, oligos overlapping
    entries in the variant reference are either removed from the database or flagged
    with the associated hit information.

    :param remove_hits: If True, oligos overlapping variants are removed. If False, they are flagged.
    :type remove_hits: bool
    :param filter_name: Name used to label the filter and create subdirectories for outputs.
    :type filter_name: str
    :param dir_output: Directory path where output files will be saved.
    :type dir_output: str
    :param usecols_search_output: Column indices to extract from the BEDTools intersection output.
    :type usecols_search_output: list
    :param names_search_output: Column names corresponding to the selected columns from the BEDTools output.
    :type names_search_output: list
    """

    def __init__(
        self,
        remove_hits: bool = True,
        filter_name: str = "variants_filter",
        dir_output: str = "output",
        usecols_search_output: list[int] = [3, 8],
        names_search_output: list[str] = [
            "query",
            "reference",
        ],
    ) -> None:
        """Constructor for the SNPFilter class."""
        super().__init__(remove_hits, filter_name, dir_output)

        self.usecols_search_output = usecols_search_output
        self.names_search_output = names_search_output

    def _create_reference(
        self,
        n_jobs: int,  # not utilized in this filter
    ) -> str:
        """
        Create a reference file from the variant reference database.

        :param n_jobs: Number of parallel jobs to use for processing. Note: This parameter is not utilized in this filter.
        :type n_jobs: int
        :return: Path to the written reference file.
        :rtype: str
        """
        if self.reference_database is None:
            raise DatabaseError("reference_database must be set before calling _create_reference")

        file_reference = self.reference_database.write_database_to_file(
            filename=f"db_reference_{self.filter_name}",
            dir_output=self.dir_output,
        )
        return file_reference

    def apply(
        self,
        oligo_database: OligoDatabase,
        sequence_type: str | None = None,
        n_jobs: int = 1,
    ) -> OligoDatabase:
        """
        Apply the variant filter to the given oligo database.

        Intersects oligo positions with the variant reference file using BEDTools,
        and either removes or flags oligos depending on the filter mode.

        :param oligo_database: The OligoDatabase instance containing oligonucleotide sequences and their associated properties. This database stores oligo data organized by genomic regions and can be used for filtering, property calculations, set generation, and output operations.
        :type oligo_database: OligoDatabase
        :param sequence_type: Type of sequence being processed.  Note: This parameter is not utilized in this filter.
        :type sequence_type: str | None
        :param n_jobs: Number of parallel jobs to use for processing.
        :type n_jobs: int
        :return: The filtered OligoDatabase.
        :rtype: OligoDatabase
        """
        file_reference = self._create_reference(n_jobs=n_jobs)

        # run search in parallel for each region
        region_ids = list(oligo_database.database.keys())
        name = " ".join(string.capitalize() for string in self.filter_name.split("_"))
        with joblib_progress(description=f"Specificity Filter: {name}", total=len(region_ids)):
            Parallel(n_jobs=n_jobs, prefer="threads", require="sharedmem")(
                delayed(self._run_filter)(
                    region_id=region_id,
                    oligo_database=oligo_database,
                    file_reference=file_reference,
                    mode=int(self.remove_hits),
                )
                for region_id in region_ids
            )

        self._remove_reference(file_reference)

        return oligo_database

    def _run_filter(
        self,
        region_id: str,
        oligo_database: OligoDatabase,
        file_reference: str,
        mode: int,
    ) -> None:
        """
        Execute the filtering logic for a single region.

        Performs a BEDTools intersection of the oligos from the region with the variant reference,
        and applies the appropriate filtering or annotation logic based on the mode.

        :param region_id: Region ID to process.
        :type region_id: str
        :param oligo_database: The OligoDatabase instance containing oligonucleotide sequences and their associated properties. This database stores oligo data organized by genomic regions and can be used for filtering, property calculations, set generation, and output operations.
        :type oligo_database: OligoDatabase
        :param file_reference: Path to the reference file used for BED intersection.
        :type file_reference: str
        :param mode: Operation mode — 0: flag hits, 1: remove hits.
        :type mode: int
        """
        file_oligo_database = oligo_database.write_database_to_bed(
            filename=f"{region_id}.bed", dir_output=self.dir_output, region_ids=region_id
        )
        file_bed_results = os.path.join(self.dir_output, f"bed_results_{region_id}.txt")

        get_intersection(file_A=file_oligo_database, file_B=file_reference, file_bed_out=file_bed_results)

        # read the reuslts of the bed seatch
        table_hits = self._read_search_output(
            file_search_results=file_bed_results,
            names_search_output=self.names_search_output,
            usecols=self.usecols_search_output,
            parse_query=True,
            parse_reference=False,
        )

        if mode == 0:
            oligos_with_hits_region = {region_id: table_hits["query"].unique()}
            oligos_with_hits_region_properties = (
                table_hits.groupby("query")["reference"].apply(list).to_dict()
            )
            self._flag_hits_in_database(
                oligo_database=oligo_database,
                region_ids=region_id,
                oligos_with_hits=oligos_with_hits_region,
                oligos_with_hits_properties=oligos_with_hits_region_properties,
            )
        elif mode == 1:
            oligos_with_hits_region = {region_id: table_hits["query"].unique()}
            self._filter_hits_from_database(
                oligo_database=oligo_database,
                region_ids=region_id,
                oligos_with_hits=oligos_with_hits_region,
            )
        else:
            raise ConfigurationError(
                f"Mode '{mode}' is not available. Choose mode=0 for flagging the hits in the database, "
                f"or mode=1 for removing hits from the database."
            )

        # remove temporary files
        os.remove(file_oligo_database)
        os.remove(file_bed_results)
