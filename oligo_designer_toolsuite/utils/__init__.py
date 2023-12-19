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

Functions for processing sequences:
- get_sequence_from_annotation: Retrieves sequences based on genomic coordinates from a BED file.
- get_complement_regions: Computes complement regions for a given set of genomic regions.

"""

from ._sequence_parser import FastaParser, GffParser

from ._sequence_processor import get_sequence_from_annotation, get_complement_regions

from ._database_processor import merge_databases, collapse_info_for_duplicated_sequences

from ._utils import check_if_list, check_tsv_format

from ._sequence_design import (
    generate_random_sequence,
    generate_binary_sequences,
    generate_codebook,
    get_barcode,
    SCRINSHOT_or_ISS_backbone_sequence,
    convert_complementary_seq_to_arms,
    create_seqfish_plus_barcodes,
)


__all__ = [
    "FastaParser",
    "GffParser",
    "get_sequence_from_annotation",
    "get_complement_regions",
    "merge_databases",
    "collapse_info_for_duplicated_sequences",
    "check_if_list",
    "check_tsv_format",
    "generate_random_sequence",
    "generate_binary_sequences",
    "generate_codebook",
    "get_barcode",
    "SCRINSHOT_or_ISS_backbone_sequence",
    "convert_complementary_seq_to_arms",
    "create_seqfish_plus_barcodes",
]
