############################################
# imports
############################################

import warnings
from collections import defaultdict
from itertools import chain

from effidict import LRUDict

from .._constants import SEPARATOR_OLIGO_ID

############################################
# Collection of utility functions
############################################


def merge_databases(database1, database2, max_in_memory=10):
    """Merge two databases, combining their content while handling potential overlapping oligo sequences.

    :param database1: The first database.
    :type database1: dict
    :param database2: The second database.
    :type database2: dict
    :return: The merged database.
    :rtype: dict
    """

    def _get_sequence_as_key(database):
        """Modify the structure of the input database by using oligo sequences as keys.

        :param database: The input database with regions and associated oligo information.
        :type database: dict
        :return: A modified database with oligo sequences as keys.
        :rtype: dict
        """
        database_modified = LRUDict(
            max_in_memory=max_in_memory,
            storage_path=database.storage_path,
        )
        for region in database.keys():
            database_modified[region] = {}

        for region, values in database.items():
            for oligo_id, oligo_info in values.items():
                oligo_sequence = oligo_info["oligo"]
                oligo_info.pop("oligo")
                database_modified[region][oligo_sequence] = oligo_info
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
        for region, values in database_in.items():
            for oligo_sequence, oligo_info in values.items():
                if oligo_sequence in database_tmp[region]:
                    oligo_info_merged = collapse_info_for_duplicated_sequences(
                        database_tmp[region][oligo_sequence], oligo_info
                    )
                    database_tmp[region][oligo_sequence] = oligo_info_merged
                else:
                    database_tmp[region][oligo_sequence] = oligo_info
        return database_tmp

    database_tmp = LRUDict(
        max_in_memory=max_in_memory,
        storage_path=database1.storage_path,
    )
    for region in chain(database1.keys(), database2.keys()):
        database_tmp[region] = {}

    db1_sequences_as_keys = _get_sequence_as_key(database1)
    db2_sequences_as_keys = _get_sequence_as_key(database2)

    database_tmp = _add_database_content(database_tmp, db1_sequences_as_keys)
    database_tmp = _add_database_content(database_tmp, db2_sequences_as_keys)

    database_merged = LRUDict(
        max_in_memory=max_in_memory,
        storage_path=database1.storage_path,
    )
    for region in database_tmp.keys():
        database_merged[region] = {}

    for region, value in database_tmp.items():
        i = 1
        for oligo_sequence, oligo_info in value.items():
            oligo_id = f"{region}{SEPARATOR_OLIGO_ID}{i}"
            oligo_seq_info = {"oligo": oligo_sequence} | oligo_info
            database_merged[region][oligo_id] = oligo_seq_info
            i += 1

    return database_merged


def collapse_info_for_duplicated_sequences(oligo_info1, oligo_info2):
    """Collapse information for duplicated sequences by combining information from two dictionaries.

    :param oligo_info1: The first dictionary of information.
    :type oligo_info1: dict
    :param oligo_info2: The second dictionary of information.
    :type oligo_info2: dict
    :return: A dictionary containing combined information from the input dictionaries.
    :rtype: dict
    """

    def _is_list_of_lists(item):
        """Check if the given item is a list of lists.

        :param item: The item to check.
        :type item: Any
        :return: True if the item is a list containing only lists, False otherwise.
        :rtype: bool
        """
        return isinstance(item, list) and all(isinstance(subitem, list) for subitem in item)

    oligo_info = defaultdict(list)

    for d in (oligo_info1, oligo_info2):
        for key, values in d.items():
            if key not in oligo_info:
                if _is_list_of_lists(values):
                    oligo_info[key] = values
                else:
                    oligo_info[key] = [values]
            else:
                if _is_list_of_lists(values):
                    oligo_info[key].extend(values)
                else:
                    oligo_info[key].extend([values])

    oligo_info = dict(oligo_info)

    return oligo_info


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


def make_entries_list_of_list(oligo_info):
    """Tranform the entries of a disctionary into a list of lists.

    :param oligo_info: The dictionary to convert.
    :type oligo_info: dict
    :return: Dictionary tranformed.
    :rtype: list
    """
    return {key: [value] for key, value in oligo_info.items()}