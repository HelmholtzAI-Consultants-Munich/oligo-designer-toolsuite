############################################
# imports
############################################

import os
import csv
import yaml
import time
import uuid
import warnings

############################################
# Collection of utility functions
############################################


class CustomYamlDumper(yaml.SafeDumper):
    """
    A custom YAML dumper class to modify the default indentation and formatting behavior in PyYAML.

    This class extends the `yaml.SafeDumper` to provide custom behaviors for YAML dumping,
    including increased indentation levels and specific formatting for lists and dictionaries.

    :param flow: Indicates if the current structure is in flow style.
    :type flow: bool

    :param data: The list or dictionary data to be represented in YAML format.
    :type data: list or dict
    """

    def increase_indent(self, flow=False):
        return super(CustomYamlDumper, self).increase_indent(flow, False)

    def represent_list(self, data):
        return self.represent_sequence("tag:yaml.org,2002:seq", data, flow_style=True)

    def represent_dict(self, data):
        return self.represent_mapping("tag:yaml.org,2002:map", data, flow_style=False)

    # Disables the use of aliases by returning True for all data
    def ignore_aliases(self, data):
        return True


def check_if_dna_sequence(seq: str, valid_characters: list = ["A", "C", "T", "G"]) -> bool:
    """
    Checks if a given sequence is a valid DNA sequence containing only specified characters.

    :param seq: The DNA sequence to check.
    :type seq: str
    :param valid_characters: A list of valid characters that the sequence should contain, defaults to ["A", "C", "T", "G"].
    :type valid_characters: list
    :return: `True` if the sequence is valid, `False` otherwise.
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


def check_if_key_exists(nested_dict: dict, key: str) -> bool:
    """
    Checks if a given key exists within a nested dictionary.

    :param nested_dict: The nested dictionary to search.
    :type nested_dict: dict
    :param key: The key to search for within the nested dictionary.
    :type key: str
    :return: `True` if the key exists, `False` otherwise.
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


def check_if_list(obj: any) -> list:
    """
    Ensures that the given object is returned as a list. Wraps non-list objects in a list.

    :param obj: The object to check and possibly convert.
    :type obj: any
    :return: The object wrapped in a list if it wasn't already a list.
    :rtype: list
    """
    if obj:
        obj = [obj] if not isinstance(obj, list) else obj
    return obj


def check_if_list_of_lists(obj: any) -> list:
    """
    Ensures that the given object is returned as a list of lists.

    :param obj: The object to check and possibly convert.
    :type obj: any
    :return: The object as a list of lists.
    :rtype: list
    """
    if isinstance(obj, list):
        # Check if it's a list of lists
        if all(isinstance(subitem, list) for subitem in obj):
            # Already a list of lists
            return obj
        else:
            # Convert the single list into a list of lists
            return [obj]
    else:
        # Wrap the non-list obj in a list of lists
        return [[obj]]


def check_tsv_format(file: str) -> list:
    """
    Checks if a given file is in valid TSV (Tab-Separated Values) format.

    :param file: The path to the TSV file to check.
    :type file: str
    :return: `True` if the file is in valid TSV format, `False` otherwise.
    :rtype: list
    """
    with open(file, "r") as tsv:
        read_tsv = csv.reader(tsv, delimiter="\t")
        return any(read_tsv)


def generate_unique_filename(dir_output: str, base_name: str, extension: str = "") -> str:
    """
    Generates a unique filename based on the current timestamp and a random UUID.

    :param dir_output: The directory where the file will be saved.
    :type dir_output: str
    :param base_name: The base name for the file.
    :type base_name: str
    :param extension: The file extension to use, defaults to an empty string.
    :type extension: str
    :return: The unique filename with the specified directory, base name, and extension.
    :rtype: str
    """
    timestamp = time.strftime("%Y%m%d-%H%M%S")
    unique_id = uuid.uuid4().hex
    filename = f"{base_name}_{timestamp}_{unique_id}.{extension}"
    filename = os.path.join(dir_output, filename)
    return filename
