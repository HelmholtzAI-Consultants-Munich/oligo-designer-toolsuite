############################################
# imports
############################################

import csv

############################################
# Collection of utility functions
############################################


def check_if_list(obj):
    """
    Check if the input is a list, otherwise ensure that the input is converted to a list.

    :param obj: The input object.
    :type obj: Any
    :return: A list containing the input object or the input object itself if it's already a list.
    :rtype: list
    """
    if obj:
        obj = [obj] if not isinstance(obj, list) else obj
    return obj


def check_tsv_format(file):
    """Check if a given TSV file has any content.

    :param file: The path to the TSV file.
    :type file: str
    :return: True if the TSV file has content, False otherwise.
    :rtype: bool
    """
    with open(file, "r") as tsv:
        read_tsv = csv.reader(tsv, delimiter="\t")
        return any(read_tsv)
