"""
    This module provides functions related to the user interface. Currently supported UIs: config files and cli arguments
"""

from ._config import generate_config
from ._cli import initialize_parameters, extract_arguments_to_dict

__all__ = ["generate_config", "initialize_parameters", "extract_arguments_to_dict"]
