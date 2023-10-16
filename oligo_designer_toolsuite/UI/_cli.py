import os
import yaml
from argparse import ArgumentParser


def initialize_parameters(parser: ArgumentParser, config_dir):
    """
    Initialize parser arguments based on configurations from provided YAML files.

    :param parser: The ArgumentParser object to add arguments to.
    :type  parser: ArgumentParser
    :param config_dir: Directory containing the YAML configuration files ('basic.yaml' and 'advanced.yaml').
    :type  config_dir: str
    """
    basic_config, basic_help_texts = read_yaml(os.path.join(config_dir, "basic.yaml"))
    advanced_config, advanced_help_texts = read_yaml(
        os.path.join(config_dir, "advanced.yaml")
    )

    flat_basic_config = flatten_dict(basic_config)
    flat_advanced_config = flatten_dict(advanced_config)

    combined_config = {**flat_basic_config, **flat_advanced_config}
    combined_help_texts = {**basic_help_texts, **advanced_help_texts}

    for key, value in combined_config.items():
        arg_name = key.replace("_", "-")
        parser.add_argument(
            f"--{arg_name}",
            help=combined_help_texts.get(key, ""),
            type=type(value) if value is not None else str,
            default=value,
            metavar="",
        )


def extract_arguments_to_dict(parser: ArgumentParser, args):
    """
    Extracts the parsed arguments from the parser and converts it to a nested dictionary.

    :param parser: The ArgumentParser object containing the arguments.
    :type  parser: ArgumentParser
    :param args: The parsed arguments returned by `parser.parse_args()`.
    :type  args: Namespace
    :return: A nested dictionary containing the parsed arguments.
    :rtype: dict
    """
    # Convert the Namespace object to a dictionary
    arg_dict = vars(args)
    # Unflatten the dictionary to get the nested structure
    nested_dict = unflatten_dict(arg_dict)
    return nested_dict


# Helper functions


def read_yaml(file_path):
    """
    Helper function that parses a YAML file and extract help texts from comments.

    :param file_path: Path to the YAML file to be read.
    :type  file_path: str
    :return: A tuple containing two dictionaries: The parsed YAML content and the extracted help texts.
    :rtype: tuple
    """
    help_texts = {}
    with open(file_path, "r") as file:
        for line in file:
            if "#" in line:
                key = line.split(":")[0].strip()
                help_text = line.split("#")[1].strip()
                help_texts[key] = help_text

        # Go back to the start of the file and parse it as YAML
        file.seek(0)
        config = yaml.safe_load(file)
    return config, help_texts


def flatten_dict(d, parent_key="", sep="__"):
    """
    Helper function to recursively flatten a dictionary.

    :param d: Dictionary to be flattened.
    :type  d: dict
    :param parent_key: Key to start with in the flattened version. Default is an empty string.
    :type  parent_key: str, optional
    :param sep: Separator to use between keys in the flattened dictionary. Default is '__'.
    :type  sep: str, optional
    :return: The flattened dictionary.
    :rtype: dict
    """
    items = {}
    for k, v in d.items():
        new_key = parent_key + sep + k if parent_key else k
        if isinstance(v, dict):
            items.update(flatten_dict(v, new_key, sep=sep))
        else:
            items[new_key] = v
    return items


def unflatten_dict(d, sep="__"):
    """
    Helper function to recursively unflatten a dictionary.

    :param d: Dictionary to be unflattened.
    :type  d: dict
    :param sep: Separator used between keys in the flattened dictionary. Default is '--'.
    :type  sep: str, optional
    :return: The unflattened dictionary.
    :rtype: dict
    """

    def set_nested_item(dct, keys, value):
        """Set item in nested dictionary."""
        for key in keys[:-1]:
            dct = dct.setdefault(key, {})
        dct[keys[-1]] = value

    result = {}
    for key, value in d.items():
        keys = key.split(sep)
        set_nested_item(result, keys, value)
    return result
