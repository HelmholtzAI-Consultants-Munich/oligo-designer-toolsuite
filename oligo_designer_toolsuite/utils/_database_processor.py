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

    def _get_sequence_as_key(database: dict, regions: list) -> dict:
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


def format_oligo_attributes(oligo_attributes: dict) -> dict:
    for key, value in oligo_attributes.items():
        if key not in get_args(_TYPES_SEQ):
            oligo_attributes[key] = check_if_list_of_lists(value)
    return oligo_attributes


def check_if_region_in_database(
    database: dict, region_ids: list, write_regions_with_insufficient_oligos: bool, file_removed_regions: str
) -> None:

    keys = list(database.keys())
    for region_id in region_ids:
        if region_id not in keys:
            warnings.warn(f"Region {region_id} not available in reference file.")
            if write_regions_with_insufficient_oligos:
                with open(file_removed_regions, "a") as hanlde:
                    hanlde.write(f"{region_id}\t{'Not in Annotation'}\n")


def flatten_attribute_list(attribute: list) -> Union[list, str, int, float, bool]:
    flattened_attribute_list = [
        item
        for sublist in (attribute if isinstance(attribute, list) else [attribute])
        for item in (sublist if isinstance(sublist, list) else [sublist])
    ]
    if len(flattened_attribute_list) == 1:
        return flattened_attribute_list[0]
    return flattened_attribute_list
