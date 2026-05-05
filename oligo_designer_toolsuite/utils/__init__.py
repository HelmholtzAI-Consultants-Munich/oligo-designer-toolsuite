"""
This module provides utilities for processing databases, parsing sequences and checking different file or object formats.
"""

from ._checkers_and_helpers import (
    CustomYamlDumper,
    cast_to_int,
    cast_to_list,
    cast_to_list_of_lists,
    cast_to_string,
    check_if_dna_sequence,
    check_if_key_exists,
    check_tsv_format,
    generate_unique_filename,
)
from ._database_processor import (
    check_if_key_in_database,
    check_if_region_in_database,
    collapse_properties_for_duplicated_sequences,
    flatten_property_list,
    format_oligo_properties,
    merge_databases,
)
from ._sequence_parser import FastaParser, GffParser, VCFParser
from ._sequence_processor import (
    append_nucleotide_to_sequences,
    count_kmer_abundance,
    get_complement_regions,
    get_intersection,
    get_sequence_from_annotation,
    remove_index_files,
)

__all__ = [
    "FastaParser",
    "GffParser",
    "VCFParser",
    "CustomYamlDumper",
    "check_if_dna_sequence",
    "check_if_key_exists",
    "cast_to_list",
    "cast_to_list_of_lists",
    "cast_to_int",
    "cast_to_string",
    "check_tsv_format",
    "check_if_region_in_database",
    "generate_unique_filename",
    "collapse_properties_for_duplicated_sequences",
    "format_oligo_properties",
    "check_if_key_in_database",
    "merge_databases",
    "flatten_property_list",
    "get_complement_regions",
    "get_sequence_from_annotation",
    "get_intersection",
    "append_nucleotide_to_sequences",
    "remove_index_files",
    "count_kmer_abundance",
]
