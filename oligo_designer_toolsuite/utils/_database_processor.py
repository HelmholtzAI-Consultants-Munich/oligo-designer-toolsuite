############################################
# imports
############################################

from itertools import chain
from collections import defaultdict

############################################
# Collection of utility functions
############################################


def merge_databases(database1, database2):
    def _get_sequence_as_key(database):
        database_modified = {region: {} for region in database.keys()}
        for region, values in database.items():
            for oligo_id, oligo_info in values.items():
                oligo_sequence = oligo_info["oligo"]
                oligo_info.pop("oligo")
                database_modified[region][oligo_sequence] = oligo_info
        return database_modified

    def _add_database_content(database_tmp, database_in):
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
            oligo_id = f"{region}::{i}"
            oligo_seq_info = {"oligo": oligo_sequence} | oligo_info
            database_merged[region][oligo_id] = oligo_seq_info
            i += 1

    return database_merged


def collapse_info_for_duplicated_sequences(oligo_info1, oligo_info2):
    def _is_list_of_lists(item):
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
