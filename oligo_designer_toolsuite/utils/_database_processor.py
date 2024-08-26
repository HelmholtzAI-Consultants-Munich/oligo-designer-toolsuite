############################################
# imports
############################################

import warnings

from typing import Union, get_args
from effidict import LRUPickleDict

from oligo_designer_toolsuite._constants import SEPARATOR_OLIGO_ID, _TYPES_SEQ
from ._checkers_and_helpers import check_if_list_of_lists
from .._constants import _TYPES_SEQ

############################################
# Collection of utility functions
############################################


def merge_databases(
    database1: dict, database2: dict, dir_cache_files: str, lru_db_max_in_memory: int
) -> dict:
    """
    Merges two oligo databases by combining their content based on sequence keys,
    ensuring that sequences with the same oligo are merged, and avoiding duplicates.

    :param database1: The first database to be merged.
    :type database1: dict
    :param database2: The second database to be merged.
    :type database2: dict
    :param dir_cache_files: Directory to store cache files used for merging.
    :type dir_cache_files: str
    :param lru_db_max_in_memory: Maximum number of entries to keep in memory for the LRU (Least Recently Used) cache.
    :type lru_db_max_in_memory: int
    :return: The merged database.
    :rtype: dict
    """

    def _get_sequence_as_key(database: dict, regions: list) -> dict:
        """
        Converts oligo sequences to dictionary keys, grouping oligo attributes by sequence for each specified region.

        :param database: The database containing sequences and their attributes.
        :type database: dict
        :param regions: List of regions within the database to process.
        :type regions: list
        :return: A dictionary with sequences as keys and oligo attributes as values.
        :rtype: dict
        """
        database_modified = LRUPickleDict(
            max_in_memory=lru_db_max_in_memory,
            storage_path=dir_cache_files,
        )
        for region in regions:
            database_modified[region] = {}
            database_region = database[region]
            for oligo_id, oligo_attributes in database_region.items():
                oligo_sequence = oligo_attributes["oligo"]
                oligo_attributes.pop("oligo")
                database_modified[region][oligo_sequence] = oligo_attributes
        return database_modified

    def _add_database_content(database_merged_tmp: dict, database_in_tmp: dict) -> dict:
        """
        Merges oligo attributes from two databases, ensuring sequences with the same oligo are combined and attributes are updated.

        :param database_merged_tmp: The dictionary to which content is added.
        :type database_merged_tmp: dict
        :param database_in_tmp: The dictionary containing new content to merge.
        :type database_in_tmp: dict
        :return: The updated dictionary with merged oligo attributes.
        :rtype: dict
        """
        for region, database_region in database_in_tmp.items():
            for oligo_sequence, oligo_attributes in database_region.items():
                if oligo_sequence in database_merged_tmp[region]:
                    oligo_attributes_merged = collapse_attributes_for_duplicated_sequences(
                        database_merged_tmp[region][oligo_sequence], oligo_attributes
                    )
                    database_merged_tmp[region][oligo_sequence] = oligo_attributes_merged
                else:
                    database_merged_tmp[region][oligo_sequence] = oligo_attributes
        return database_merged_tmp

    # keys that are in both dicts
    regions_intersection = list(set(database1) & set(database2))

    database_merged = LRUPickleDict(
        max_in_memory=lru_db_max_in_memory,
        storage_path=dir_cache_files,
    )
    for region in regions_intersection:
        database_merged[region] = {}

    # only loop over entries that have keys in both dicts
    db1_sequences_as_keys = _get_sequence_as_key(database1, regions_intersection)
    db2_sequences_as_keys = _get_sequence_as_key(database2, regions_intersection)

    database_merged = _add_database_content(database_merged, db1_sequences_as_keys)
    database_merged = _add_database_content(database_merged, db2_sequences_as_keys)

    database_concat = LRUPickleDict(
        max_in_memory=lru_db_max_in_memory,
        storage_path=dir_cache_files,
    )
    for region in regions_intersection:
        database_concat[region] = {}

    for region, database_merged_region in database_merged.items():
        i = 1
        for oligo_sequence, oligo_attributes in database_merged_region.items():
            oligo_id = f"{region}{SEPARATOR_OLIGO_ID}{i}"
            oligo_seq_info = {"oligo": oligo_sequence} | oligo_attributes
            database_concat[region][oligo_id] = oligo_seq_info
            i += 1

    # add entries with keys in only one dict
    for region in set(database1) - set(database2):
        database_concat[region] = database1[region]
    for region in set(database2) - set(database1):
        database_concat[region] = database2[region]

    return database_concat


def collapse_attributes_for_duplicated_sequences(oligo_attributes1: dict, oligo_attributes2: dict) -> dict:
    """
    Merges two dictionaries of oligo attributes, combining values for non-sequence keys and issuing warnings if sequences for the same oligo ID.

    :param oligo_attributes1: The first dictionary of oligo attributes.
    :type oligo_attributes1: dict
    :param oligo_attributes2: The second dictionary of oligo attributes.
    :type oligo_attributes2: dict
    :return: A merged dictionary with combined oligo attributes.
    :rtype: dict
    """
    oligo_attributes = {}

    if oligo_attributes1 == oligo_attributes2:
        return oligo_attributes1

    for d in (oligo_attributes1, oligo_attributes2):
        for key, values in d.items():
            if key not in oligo_attributes:
                oligo_attributes[key] = values
            else:
                if key in get_args(_TYPES_SEQ) and oligo_attributes[key] != values:
                    warnings.warn(
                        f"Values for key {key} are different in the two oligo_attributes dictionaries."
                    )
                elif key not in get_args(_TYPES_SEQ):
                    oligo_attributes[key].extend(values)

    return oligo_attributes


def check_if_region_in_database(
    database: dict, region_ids: list, write_regions_with_insufficient_oligos: bool, file_removed_regions: str
) -> None:
    """
    Checks if specified regions are present in the database and logs missing regions.

    :param database: The database.
    :type database: dict
    :param region_ids: A list of region IDs to check.
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


def format_oligo_attributes(oligo_attributes: dict) -> dict:
    """
    Ensures that the values in an oligo attributes dictionary are formatted as lists of lists.

    :param oligo_attributes: The dictionary of oligo attributes to format.
    :type oligo_attributes: dict
    :return: The formatted dictionary with lists of lists for non-sequence keys.
    :rtype: dict
    """
    for key, value in oligo_attributes.items():
        if key not in get_args(_TYPES_SEQ):
            oligo_attributes[key] = check_if_list_of_lists(value)
    return oligo_attributes


def flatten_attribute_list(attribute: list) -> Union[list, str, int, float, bool]:
    """
    Flattens a nested list of attributes into a single list, or returns the item if only one element remains.

    :param attribute: The list or nested list of attributes to flatten.
    :type attribute: list
    :return: A flattened list or a single item if only one element exists.
    :rtype: Union[list, str, int, float, bool]
    """
    flattened_attribute_list = [
        item
        for sublist in (attribute if isinstance(attribute, list) else [attribute])
        for item in (sublist if isinstance(sublist, list) else [sublist])
    ]
    if len(flattened_attribute_list) == 1:
        return flattened_attribute_list[0]
    return flattened_attribute_list
