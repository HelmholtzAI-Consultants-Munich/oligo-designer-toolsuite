############################################
# imports
############################################

from abc import abstractmethod

from joblib import Parallel, delayed
from joblib_progress import joblib_progress

from oligo_designer_toolsuite.database import OligoDatabase

############################################
# Oligoset Selection Classes
############################################


class BaseOligoSelection:
    """
    Base class for oligo set generators. Defines the common structure for applying set generation
    across all regions of an oligo database in parallel. Subclasses must implement
    :meth:`_get_oligo_sets_for_region` to define how sets are generated for a single region.

    The :meth:`apply` method iterates over all regions in the database, runs the region-level
    generator in parallel via joblib, and returns the updated database with
    ``oligo_database.oligosets[region_id]`` populated for each region.
    """

    def apply(
        self,
        oligo_database: OligoDatabase,
        sequence_type: str,
        n_sets: int = 1,
        n_jobs: int = 1,
    ) -> OligoDatabase:
        """
        Applies the oligo set generation process to an oligo database. For each region in the
        database, the subclass implementation of :meth:`_get_oligo_sets_for_region` is called to
        generate oligo sets. The process is parallelized across regions using joblib.

        Oligosets are stored in ``oligo_database.oligosets``, which is a dictionary with region
        names as keys and oligoset DataFrames as values. The structure of each DataFrame is:

        +-------------+----------+----------+----------+-------+----------+-------------+-------+
        | oligoset_id | oligo_0  | oligo_1  | oligo_2  |  ...  | oligo_n  | set_score_1 |  ...  |
        +-------------+----------+----------+----------+-------+----------+-------------+-------+

        :param oligo_database: The OligoDatabase instance containing oligonucleotide sequences and
            their associated properties, organized by genomic regions.
        :type oligo_database: OligoDatabase
        :param sequence_type: Type of sequence being processed (e.g. ``"oligo"``).
        :type sequence_type: str
        :param n_sets: The number of oligo sets to generate per region.
        :type n_sets: int, optional
        :param n_jobs: Number of parallel jobs to use for processing regions.
        :type n_jobs: int, optional
        :return: The updated oligo database with generated oligo sets for each region.
        :rtype: OligoDatabase
        """
        region_ids = list(oligo_database.database.keys())
        with joblib_progress(description="Find Oligosets", total=len(region_ids)):
            Parallel(n_jobs=n_jobs, prefer="threads", require="sharedmem")(
                delayed(self._get_oligo_sets_for_region)(
                    oligo_database=oligo_database,
                    sequence_type=sequence_type,
                    region_id=region_id,
                    n_sets=n_sets,
                )
                for region_id in region_ids
            )
        return oligo_database

    @abstractmethod
    def _get_oligo_sets_for_region(
        self,
        oligo_database: OligoDatabase,
        sequence_type: str,
        region_id: str,
        n_sets: int,
    ) -> None:
        """
        Generates oligo sets for a single region and stores the result in
        ``oligo_database.oligosets[region_id]``. Must be implemented by subclasses.

        :param oligo_database: The OligoDatabase instance; will be updated in place with
            ``oligo_database.oligosets[region_id]`` set to a DataFrame of oligo sets.
        :type oligo_database: OligoDatabase
        :param sequence_type: Type of sequence being processed.
        :type sequence_type: str
        :param region_id: Region ID to process.
        :type region_id: str
        :param n_sets: The number of oligo sets to generate for this region.
        :type n_sets: int
        """
