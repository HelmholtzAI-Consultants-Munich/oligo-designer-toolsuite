############################################
# imports
############################################

import csv

############################################
# Collection of utility functions
############################################


def check_if_list(obj):
    if obj:
        obj = [obj] if not isinstance(obj, list) else obj
    return obj


def check_tsv_format(file):
    """File check. Is the file a tsv file?

    :param file: Path to file.
    :type file: str
    :return: Returns True if file has tsv format.
    :rtype: bool
    """
    with open(file, "r") as tsv:
        read_tsv = csv.reader(tsv, delimiter="\t")
        return any(read_tsv)
