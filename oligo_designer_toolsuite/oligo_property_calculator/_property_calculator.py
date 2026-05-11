############################################
# imports
############################################

from typing import Any

from joblib import Parallel, delayed
from joblib_progress import joblib_progress

from oligo_designer_toolsuite.database import OligoDatabase
from oligo_designer_toolsuite.oligo_property_calculator import BaseProperty
from oligo_designer_toolsuite.utils import check_if_key_in_database

############################################
# Property Calculator Class
############################################


class PropertyCalculator:
    """
    A class for applying multiple property calculators to oligonucleotides in an OligoDatabase.

    The `PropertyCalculator` class allows you to apply a list of property calculators (subclasses of `BaseProperty`) to an OligoDatabase.
    The properties are calculated in parallel across all regions of the database, and the calculated values are stored as properties in the database.

    :param properties: A list of property calculators to apply to oligonucleotides.
    :type properties: list[BaseProperty]
    """

    def __init__(self, properties: list[BaseProperty]) -> None:
        """Constructor for the PropertyCalculator class."""
        self.properties = properties

    def apply(self, oligo_database: OligoDatabase, sequence_type: str, n_jobs: int = 1) -> OligoDatabase:
        """
        Apply the property calculators to all oligonucleotides in the OligoDatabase and update
        the database with the calculated property values.

        :param oligo_database: The OligoDatabase instance containing oligonucleotide sequences and their associated properties. This database stores oligo data organized by genomic regions and can be used for filtering, property calculations, set generation, and output operations.
        :type oligo_database: OligoDatabase
        :param sequence_type: Type of sequence being processed.
        :type sequence_type: str
        :param n_jobs: Number of parallel jobs to use for processing. Defaults to 1.
        :type n_jobs: int
        :return: The updated OligoDatabase with the calculated properties.
        :rtype: OligoDatabase
        """
        assert check_if_key_in_database(
            oligo_database.database, sequence_type
        ), f"Sequence type '{sequence_type}' not found in database."

        region_ids = list(oligo_database.database.keys())
        with joblib_progress(description="Property Calculator", total=len(region_ids)):
            Parallel(n_jobs=n_jobs, prefer="threads", require="sharedmem")(
                delayed(self._calculate_region)(oligo_database, region_id, sequence_type)
                for region_id in region_ids
            )

        return oligo_database

    def _calculate_region(self, oligo_database: OligoDatabase, region_id: str, sequence_type: str) -> None:
        """
        Calculate properties for all oligonucleotides in a specific region of the OligoDatabase.

        This method iterates through the oligonucleotides in a given region of the database,
        applying all property calculators to each oligo and updating the database with the calculated values.

        :param oligo_database: The OligoDatabase instance containing oligonucleotide sequences and their associated properties. This database stores oligo data organized by genomic regions and can be used for filtering, property calculations, set generation, and output operations.
        :type oligo_database: OligoDatabase
        :param region_id: Region ID to process.
        :type region_id: str
        :param sequence_type: Type of sequence being processed.
        :type sequence_type: str
        """
        new_oligo_property: dict[str, dict[str, Any]] = {}

        for oligo_id in oligo_database.database[region_id].keys():
            # Calculate all properties for this oligo
            for property_calc in self.properties:
                property_result = property_calc.apply(
                    oligo_database=oligo_database,
                    region_id=region_id,
                    oligo_id=oligo_id,
                    sequence_type=sequence_type,
                )
                # Merge results into the property dictionary
                if oligo_id not in new_oligo_property:
                    new_oligo_property[oligo_id] = {}
                new_oligo_property[oligo_id].update(property_result)

        # Update only this region (avoids O(regions × total_oligos) full-database scan per worker)
        oligo_database.update_oligo_properties(new_oligo_property, region_ids=region_id)
