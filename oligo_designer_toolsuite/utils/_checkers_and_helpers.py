############################################
# imports
############################################

import csv
import time
import uuid
import warnings

############################################
# Collection of utility functions
############################################


def check_if_dna_sequence(seq: str, valid_characters=["A", "C", "T", "G"]):
    """Verifies if a given sequence consists solely of valid DNA nucleotide characters (A, C, T, G).
    This function is case-insensitive and validates the entire sequence against a set of allowed characters.

    :param seq: The DNA sequence to validate.
    :type seq: str
    :param valid_characters: A list of characters representing valid nucleotides. Defaults to ["A", "C", "T", "G"].
    :type valid_characters: list
    :return: True if the sequence is valid DNA, False otherwise.
    :rtype: bool
    """
    if any(len(char) > 1 for char in valid_characters):
        raise ValueError("Valid characters must be single characters.")

    valid_characters_upper = [char.upper() for char in valid_characters]
    if not all(char.upper() in ["A", "C", "T", "G", "U"] for char in valid_characters_upper):
        warnings.warn("Valid characters should be A, C, T, G, or U.")

    if seq == "":
        return False
    return all(char.upper() in valid_characters_upper for char in seq)


def check_if_key_exists(nested_dict: dict, key: str):
    """Recursively check if a given key exists in a nested dictionary. The function searches through all levels of the nested dictionary to find the key.

    :param nested_dict: The nested dictionary to search for the key.
    :type nested_dict: dict
    :param key: The key to search for within the nested dictionary.
    :type key: str
    :return: True if the key is found anywhere within the nested dictionary, False otherwise.
    :rtype: bool
    """
    try:
        if key in nested_dict.keys():
            return True
        else:
            for value in nested_dict.values():
                if check_if_key_exists(value, key):
                    return True
    except:
        return False


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


def check_if_list_of_lists(item):
    """Check if the given item is a list of lists.

    :param item: The item to check.
    :type item: Any
    :return: True if the item is a list containing only lists, False otherwise.
    :rtype: bool
    """
    return isinstance(item, list) and all(isinstance(subitem, list) for subitem in item)


def check_tsv_format(file: str):
    """Check if a given TSV file has any content.

    :param file: The path to the TSV file.
    :type file: str
    :return: True if the TSV file has content, False otherwise.
    :rtype: bool
    """
    with open(file, "r") as tsv:
        read_tsv = csv.reader(tsv, delimiter="\t")
        return any(read_tsv)


def generate_unique_filename(base_name: str, extension: str = "") -> str:
    """Generate a unique filename using a base name, current timestamp, and a UUID.

    :param base_name: The base name for the file.
    :type base_name: str
    :param extension: The file extension, e.g., ".txt". (default is an empty string)
    :type extension: str
    :return: A unique filename.
    :rtype: str
    """
    timestamp = time.strftime("%Y%m%d-%H%M%S")
    unique_id = uuid.uuid4().hex
    filename = f"{base_name}_{timestamp}_{unique_id}.{extension}"
    return filename
