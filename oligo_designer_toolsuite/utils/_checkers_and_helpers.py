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
    def increase_indent(self, flow=False):
        return super(CustomYamlDumper, self).increase_indent(flow, False)

    def represent_list(self, data):
        return self.represent_sequence("tag:yaml.org,2002:seq", data, flow_style=True)

    def represent_dict(self, data):
        return self.represent_mapping("tag:yaml.org,2002:map", data, flow_style=False)


def check_if_dna_sequence(seq: str, valid_characters: list = ["A", "C", "T", "G"]) -> bool:
    if any(len(char) > 1 for char in valid_characters):
        raise ValueError("Valid characters must be single characters.")

    valid_characters_upper = [char.upper() for char in valid_characters]
    if not all(char.upper() in ["A", "C", "T", "G", "U"] for char in valid_characters_upper):
        warnings.warn("Valid characters should be A, C, T, G, or U.")

    if seq == "":
        return False
    return all(char.upper() in valid_characters_upper for char in seq)


def check_if_key_exists(nested_dict: dict, key: str) -> bool:
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
    if obj:
        obj = [obj] if not isinstance(obj, list) else obj
    return obj


def check_if_list_of_lists(item: any) -> list:
    if isinstance(item, list):
        # Check if it's a list of lists
        if all(isinstance(subitem, list) for subitem in item):
            # Already a list of lists
            return item
        else:
            # Convert the single list into a list of lists
            return [item]
    else:
        # Wrap the non-list item in a list of lists
        return [[item]]


def check_tsv_format(file: str) -> list:
    with open(file, "r") as tsv:
        read_tsv = csv.reader(tsv, delimiter="\t")
        return any(read_tsv)


def generate_unique_filename(dir_output: str, base_name: str, extension: str = "") -> str:
    timestamp = time.strftime("%Y%m%d-%H%M%S")
    unique_id = uuid.uuid4().hex
    filename = f"{base_name}_{timestamp}_{unique_id}.{extension}"
    filename = os.path.join(dir_output, filename)
    return filename
