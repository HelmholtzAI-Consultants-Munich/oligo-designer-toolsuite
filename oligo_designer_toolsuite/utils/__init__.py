"""
This module provides utilities for processing databases, parsing sequences and checking different file or object formats.

Classes:
- FastaParser: Parses FASTA files and provides methods for extracting sequence information.
- GffParser: Parses GFF/GTF files and provides methods for extracting annotation information.
- CustomYamlDumper: Custom YAML dumper for serializing Python objects to YAML format.

Functions for processing databases:
- check_if_region_in_database: Checks if specified regions exist in the provided OligoDatabase.
- collapse_attributes_for_duplicated_sequences: Collapses information for duplicated sequences in a OligoDatabase
- format_oligo_attributes: Formats the entries of an oligo_attributes dictionary into a list of lists.
- merge_databases: Merges two OligoDatabase by combining their content based on sequence keys.
- flatten_attribute_list: Flattens a nested list of attributes into a single list.

Functions for processing sequences:
- get_sequence_from_annotation: Retrieves sequences based on genomic coordinates from a BED file.
- get_complement_regions: Computes complement regions for a given set of genomic regions.

Functions for checking formats:
- check_if_dna_sequence: Verifies that a sequence consists only of valid DNA nucleotides (A, C, T, G), case-insensitively.
- check_if_key_exists: Recursively searches a nested dictionary to find if a specific key exists.
- check_if_list: Ensures the input is a list, converting it to one if necessary.
- check_if_list_of_lists: Ensures the input is a list of lists, converting it to one if necessary.
- check_tsv_format: Checks a TSV file for content, verifying it is not empty.
- generate_unique_filename: Generates a unique filename based on the current timestamp.
"""

from ._checkers_and_helpers import (
    CustomYamlDumper,
    check_if_dna_sequence,
    check_if_key_exists,
    check_if_list,
    check_if_list_of_lists,
    check_tsv_format,
    generate_unique_filename,
)
from ._database_processor import (
    check_if_region_in_database,
    collapse_attributes_for_duplicated_sequences,
    format_oligo_attributes,
    merge_databases,
    flatten_attribute_list,
)
from ._sequence_parser import FastaParser, GffParser
from ._sequence_processor import get_complement_regions, get_sequence_from_annotation

__all__ = [
    "FastaParser",
    "GffParser",
    "CustomYamlDumper",
    "check_if_dna_sequence",
    "check_if_key_exists",
    "check_if_list",
    "check_if_list_of_lists",
    "check_tsv_format",
    "check_if_region_in_database",
    "generate_unique_filename",
    "collapse_attributes_for_duplicated_sequences",
    "format_oligo_attributes",
    "get_complement_regions",
    "get_sequence_from_annotation",
    "merge_databases",
    "flatten_attribute_list",
]
