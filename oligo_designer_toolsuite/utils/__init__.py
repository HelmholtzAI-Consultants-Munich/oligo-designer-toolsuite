"""
This module provides all additional functionalities that are needed for oligo design.
"""

from ._database_processor import merge_databases, collapse_info_for_duplicated_sequences

from ._sequence_design import (
    generate_random_sequence,
    generate_binary_sequences,
    generate_codebook,
    get_barcode,
    SCRINSHOT_or_ISS_backbone_sequence,
    convert_complementary_seq_to_arms,
    create_seqfish_plus_barcodes,
)

from ._sequence_parser import FastaParser, GffParser

from ._sequence_processor import get_sequence_from_annotation, get_complement_regions

from ._utils import check_if_list, check_tsv_format


__all__ = [
    "merge_databases",
    "collapse_info_for_duplicated_sequences",
    "generate_random_sequence",
    "generate_binary_sequences",
    "generate_codebook",
    "get_barcode",
    "SCRINSHOT_or_ISS_backbone_sequence",
    "convert_complementary_seq_to_arms",
    "create_seqfish_plus_barcodes",
    "FastaParser",
    "GffParser",
    "get_sequence_from_annotation",
    "get_complement_regions",
    "check_if_list",
    "check_tsv_format",
]
