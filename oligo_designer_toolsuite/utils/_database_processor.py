############################################
# imports
############################################

import hashlib
import warnings
from typing import Any

from effidict import EffiDict, LRUReplacement, PickleBackend

from ._checkers_and_helpers import cast_to_list

############################################
# Collection of utility functions
############################################


def merge_databases(
    database1: EffiDict,
    database2: EffiDict,
    database_sequence_types: list[str],
    dir_cache_files: str,
    max_entries_in_memory: int,
) -> EffiDict:
    """
    Merges two oligo databases by combining their content based on regions and oligo IDs,
    ensuring that oligos with the same oligo id are merged, and avoiding duplicates.

    :param database1: The first database to be merged.
    :type database1: EffiDict
    :param database2: The second database to be merged.
    :type database2: EffiDict
    :param database_sequence_types: List of sequence type keys in the database.
    :type database_sequence_types: list[str]
    :param dir_cache_files: Directory to store cache files used for merging.
    :type dir_cache_files: str
    :param max_entries_in_memory: Maximum number of entries to keep in memory for the LRU (Least Recently Used) cache.
    :type max_entries_in_memory: int
    :return: The merged database.
    :rtype: EffiDict
    """
    backend = PickleBackend(storage_path=dir_cache_files)
    strategy = LRUReplacement(disk_backend=backend, max_in_memory=max_entries_in_memory)
    database_merged = EffiDict(disk_backend=backend, replacement_strategy=strategy)

    regions_union = database1.keys() | database2.keys()
    for region in regions_union:
        database_merged[region] = {}
        if region in database1:
            for oligo_id, oligo_properties in database1[region].items():
                database_merged[region][oligo_id] = oligo_properties
        if region in database2:
            for oligo_id, oligo_properties in database2[region].items():
                if oligo_id in database_merged[region]:
                    database_merged[region][oligo_id] = collapse_properties_for_duplicated_coordinates(
                        database_merged[region][oligo_id],
                        oligo_properties,
                        database_sequence_types,
                    )
                else:
                    database_merged[region][oligo_id] = oligo_properties

    return database_merged


def collapse_properties_for_duplicated_coordinates(
    oligo_properties1: dict[str, Any], oligo_properties2: dict[str, Any], database_sequence_types: list[str]
) -> dict[str, Any]:
    """
    Merges two dictionaries of oligo properties, combining values for non-sequence keys and issuing warnings if sequences for the same oligo ID.

    :param oligo_properties1: The first dictionary of oligo properties.
    :type oligo_properties1: dict
    :param oligo_properties2: The second dictionary of oligo properties.
    :type oligo_properties2: dict
    :param database_sequence_types: List of sequence type keys in the database.
    :type database_sequence_types: list[str]
    :return: A merged dictionary with combined oligo properties.
    :rtype: dict
    """
    oligo_properties = {}

    if oligo_properties1 == oligo_properties2:
        return oligo_properties1

    for d in (oligo_properties1, oligo_properties2):
        for key, values in d.items():
            if key not in oligo_properties:
                oligo_properties[key] = values
            else:
                if key in database_sequence_types and oligo_properties[key] != values:
                    warnings.warn(
                        f"Values for key {key} are different in the two oligo_properties dictionaries."
                    )
                elif key not in database_sequence_types:
                    oligo_properties[key].extend(values)
                    oligo_properties[key] = list(dict.fromkeys(oligo_properties[key]))

    return oligo_properties


def check_if_region_in_database(
    database: dict[str, Any],
    region_ids: list[str],
    write_regions_with_insufficient_oligos: bool,
    file_removed_regions: str,
) -> None:
    """
    Checks if specified regions are present in the database and logs missing regions.

    :param database: The database.
    :type database: dict
    :param region_ids: Region identifier(s) to process. Can be a single region ID (str) or a list of region IDs (List[str]). If None, all regions in the database are processed.
    :type region_ids: list
    :param write_regions_with_insufficient_oligos: Whether to write regions with insufficient oligos to a file.
    :type write_regions_with_insufficient_oligos: bool
    :param file_removed_regions: The file to which missing regions should be logged.
    :type file_removed_regions: str
    """
    keys = list(database.keys())
    for region_id in region_ids:
        if region_id not in keys:
            warnings.warn(f"Region {region_id} not available in reference file.")
            if write_regions_with_insufficient_oligos:
                with open(file_removed_regions, "a") as hanlde:
                    hanlde.write(f"{region_id}\t{'Not in Annotation'}\n")


def check_if_key_in_database(database: EffiDict, key: str, region_ids: str | list[str] | None = None) -> bool:
    """
    Checks if a key exists in the database.

    :param database: The database to check.
    :type database: EffiDict
    :param key: The key to check.
    :type key: str
    :param region_ids: Optional region ID(s) to check in specific region(s). Can be a single region ID (str) or a list of region IDs (list[str]). If provided, checks that the key exists in ALL specified regions. If None, checks if it exists in at least one region.
    :type region_ids: str | list[str] | None
    :return: True if the key exists in the database (in all specified regions if region_ids is provided, or in at least one region if region_ids is None), False otherwise.
    :rtype: bool
    """

    # Helper for recursive search
    def recursive_contains(d: Any, target: str) -> bool:
        """
        Recursively checks if `target` exists as a key in a nested structure.
        Top-level is EffiDict, deeper levels are plain dicts.
        """
        try:
            if isinstance(d, dict):  # this will match EffiDict at the top AND nested dicts
                if target in d:
                    return True
                return any(recursive_contains(v, target) for v in d.values())
            return False
        except Exception:
            return False

    # --- Case: region restriction ---
    if region_ids is not None:
        region_ids = cast_to_list(region_ids)
        regions_found = False  # check if at least one region is found in the database
        for region_id in region_ids:
            if region_id not in database:
                continue
            else:
                regions_found = True
            region_data = database[region_id]
            if not recursive_contains(region_data, key):
                return False
        return regions_found

    # --- Case: no region restriction → any region may match ---
    return any(recursive_contains(region_data, key) for region_data in database.values())


def format_oligo_properties(
    oligo_properties: dict[str, Any], database_sequence_types: list[str]
) -> dict[str, Any]:
    """
    Ensures that the values in an oligo properties dictionary are formatted as lists of lists.

    :param oligo_properties: The dictionary of oligo properties to format.
    :type oligo_properties: dict
    :param database_sequence_types: List of sequence type keys in the database.
    :type database_sequence_types: list[str]
    :return: The formatted dictionary with lists of lists for non-sequence keys.
    :rtype: dict
    """
    for key, value in oligo_properties.items():
        if key not in database_sequence_types:
            oligo_properties[key] = cast_to_list(value)
    return oligo_properties


def flatten_property_list(property: list[Any]) -> list[Any]:
    """
    Flattens a nested list of properties into a single list, or returns the item if only one element remains.

    :param property: The list or nested list of properties to flatten.
    :type property: list[Any]
    :return: A flattened list or a single item if only one element exists.
    :rtype: list[Any]
    """
    flattened_property_list = [
        item
        for sublist in (property if isinstance(property, list) else [property])
        for item in (sublist if isinstance(sublist, list) else [sublist])
    ]
    if len(flattened_property_list) == 1:
        return cast_to_list(flattened_property_list[0])
    return cast_to_list(flattened_property_list)


def sequence_to_id(seq: str) -> str:
    seq = seq.upper()
    return hashlib.sha256(seq.encode()).hexdigest()[:16]
