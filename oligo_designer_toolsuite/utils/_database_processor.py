############################################
# imports
############################################

from collections import defaultdict
from itertools import chain

from .._constants import SEPARATOR_OLIGO_ID

############################################
# Collection of utility functions
############################################


def merge_databases(database1, database2):
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
        database_modified = {region: {} for region in database.keys()}
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

    database_tmp = {region: {} for region in chain(database1.keys(), database2.keys())}
    database_tmp = _add_database_content(database_tmp, _get_sequence_as_key(database1))
    database_tmp = _add_database_content(database_tmp, _get_sequence_as_key(database2))

    database_merged = {region: {} for region in database_tmp.keys()}
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
                oligo_info[key].extend([values])

    oligo_info = dict(oligo_info)

    return oligo_info