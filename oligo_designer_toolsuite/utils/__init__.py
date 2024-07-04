"""
This module provides utilities for processing genomic databases, parsing sequences,
and performing various sequence-related tasks.

Classes:
- FastaParser: Parses FASTA files and provides methods for extracting sequence information.
- GffParser: Parses GFF/GTF files and provides methods for extracting annotation information.

Functions for processing databases:
- merge_databases: Merges two genomic databases.
- collapse_info_for_duplicated_sequences: Collapses information for duplicated sequences
  in a genomic database.
- filter_dabase_for_region: Filters a genomic database for sequences overlapping a given region.
- format_oligo_attributes: Format the entries of an oligo_attributes dictionary to be list of lists.
- check_if_region_in_database: Checks if specified regions exist in the provided database.
- flatten_attribute_list: Flatten a nested list of attributes into a single list.

Functions for processing sequences:
- get_sequence_from_annotation: Retrieves sequences based on genomic coordinates from a BED file.
- get_complement_regions: Computes complement regions for a given set of genomic regions.

Functions for checking formats:
- check_if_dna_sequence: Verifies a sequence consists only of valid DNA nucleotides (A, C, T, G), case-insensitively.
- check_if_key_exists: Recursively searches a nested dictionary to find if a specific key exists.
- check_if_list: Ensures the input is a list, converting it to one if necessary.
- check_if_list_of_lists: Check if the given item is a list of lists.
- check_tsv_format: Checks a TSV file for content, verifying it's not empty.

"""

from ._checkers import (
    check_if_dna_sequence,
    check_if_key_exists,
    check_if_list,
    check_if_list_of_lists,
    check_tsv_format,
)
from ._database_processor import (
    check_if_region_in_database,
    collapse_attributes_for_duplicated_sequences,
    filter_dabase_for_region,
    format_oligo_attributes,
    merge_databases,
    flatten_attribute_list,
)
from ._sequence_parser import FastaParser, GffParser
from ._sequence_processor import get_complement_regions, get_sequence_from_annotation

__all__ = [
    "FastaParser",
    "GffParser",
    "check_if_dna_sequence",
    "check_if_key_exists",
    "check_if_list",
    "check_if_list_of_lists",
    "check_tsv_format",
    "check_if_region_in_database",
    "collapse_attributes_for_duplicated_sequences",
    "filter_dabase_for_region",
    "format_oligo_attributes",
    "get_complement_regions",
    "get_sequence_from_annotation",
    "merge_databases",
    "flatten_attribute_list",
]
