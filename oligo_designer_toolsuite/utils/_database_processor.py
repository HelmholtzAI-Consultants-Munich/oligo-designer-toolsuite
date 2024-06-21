############################################
# imports
############################################

import warnings
from collections import defaultdict

from effidict import LRUPickleDict

from oligo_designer_toolsuite._constants import SEPARATOR_OLIGO_ID

from ._checkers import check_if_list_of_lists

############################################
# Collection of utility functions
############################################


def merge_databases(database1, database2, dir_cache_files, lru_db_max_in_memory):
    """Merge two databases, combining their content while handling potential overlapping oligo sequences.

    :param database1: The first database.
    :type database1: dict
    :param database2: The second database.
    :type database2: dict
    :param dir_cache_files: Directory path for the cache files.
    :type dir_cache_files: str
    :param lru_db_max_in_memory: Maximum number of dictionary entries stored in RAM, defaults to 100.
    :type lru_db_max_in_memory: int, optional
    :return: The merged database.
    :rtype: dict
    """

    def _get_sequence_as_key(database, regions):
        """Modify the structure of the input database by using oligo sequences as keys.

        :param database: The input database with regions and associated oligo information.
        :type database: dict
        :return: A modified database with oligo sequences as keys.
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

    def _add_database_content(database_tmp, database_in):
        """Add the content of a second database to an existing database, merging information for overlapping oligo sequences.

        :param database_tmp: The existing database to which the content will be added.
        :type database_tmp: dict
        :param database_in: The database containing additional content.
        :type database_in: dict
        :return: The modified database with the added content.
        :rtype: dict
        """
        for region, database_region in database_in.items():
            for oligo_sequence, oligo_attributes in database_region.items():
                if oligo_sequence in database_tmp[region]:
                    oligo_attributes_merged = collapse_attributes_for_duplicated_sequences(
                        database_tmp[region][oligo_sequence], oligo_attributes
                    )
                    database_tmp[region][oligo_sequence] = oligo_attributes_merged
                else:
                    database_tmp[region][oligo_sequence] = oligo_attributes
        return database_tmp

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


def collapse_attributes_for_duplicated_sequences(oligo_attributes1, oligo_attributes2):
    """Collapse information for duplicated sequences by combining information from two dictionaries.

    :param oligo_attributes1: The first dictionary of oligo attributes.
    :type oligo_attributes1: dict
    :param oligo_attributes2: The second dictionary of oligo attributes.
    :type oligo_attributes2: dict
    :return: A dictionary containing combined attributes from the input dictionaries.
    :rtype: dict
    """

    oligo_attributes = defaultdict(list)

    for d in (oligo_attributes1, oligo_attributes2):
        for key, values in d.items():
            if key not in oligo_attributes:
                oligo_attributes[key] = values
            else:
                if key in ["oligo", "target"] and oligo_attributes[key] != values:
                    warnings.warn(
                        f"Values for key {key} are different in the two oligo_attributes dictionaries."
                    )
                else:
                    oligo_attributes[key].extend(values)

    oligo_attributes = dict(oligo_attributes)

    return oligo_attributes


def format_oligo_attributes(oligo_attributes: dict):
    """Format the entries of an oligo_attributes dictionary to be list of lists.

    :param oligo_attributes: Dictionary containing the features of the oligo
    :type oligoligo_attributeso_info: dict
    :return: Formatted oligo_attributes dictionary.
    :rtype: dict
    """
    for key, value in oligo_attributes.items():
        if not check_if_list_of_lists(value):
            oligo_attributes[key] = [value]
    return oligo_attributes


def filter_dabase_for_region(database, region_ids):
    """Filter the provided database to include only specified region IDs.

    This internal method filters the given database to retain only the entries corresponding to the provided list
    of region IDs. If a region ID is not in the specified list, it is removed from the database.

    :param database: The database to filter.
    :type database: dict
    :param region_ids: The list of region IDs to retain in the filtered database.
    :type region_ids: list
    :return: The filtered database.
    :rtype: dict
    """
    regions = list(database.keys())
    for region in regions:
        if region not in region_ids:
            database.pop(region)
    return database


def check_if_region_in_database(
    database, region_ids, write_regions_with_insufficient_oligos, file_removed_regions
):
    """Check if specified regions exist in the provided database.

    This internal method checks whether all regions provided in the region_ids list exist in the given database.
    If a region is not found, a warning is issued, and if enabled, the information is recorded in the log file.

    :param database: The database to check for region existence.
    :type database: dict
    :param region_ids: The list of region IDs to check.
    :type region_ids: list
    :param write_regions_with_insufficient_oligos: Flag to enable writing regions with insufficient oligos to a file.
    :type write_regions_with_insufficient_oligos: bool
    :param file_removed_regions: Path to file writing regions with insufficient oligos.
    :type file_removed_regions: str
    """
    keys = list(database.keys())
    for region_id in region_ids:
        if region_id not in keys:
            warnings.warn(f"Region {region_id} not available in reference file.")
            if write_regions_with_insufficient_oligos:
                with open(file_removed_regions, "a") as hanlde:
                    hanlde.write(f"{region_id}\t{'Not in Annotation'}\n")
