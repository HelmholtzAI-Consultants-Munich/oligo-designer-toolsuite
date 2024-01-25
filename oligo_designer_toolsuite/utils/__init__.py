"""
This module provides all additional functionalities that are needed for oligo design, e.g. a class to download annotations from NCBI and Ensembl servers.
"""

from ._data_parser import (
    check_fasta_format,
    check_gff_format,
    check_tsv_format,
    get_sequence_from_annotation,
    merge_fasta,
    parse_fasta_header,
)
from ._ftp_loader import BaseFtpLoader, FtpLoaderEnsembl, FtpLoaderNCBI
from ._gff_parser import GffParser
from ._sequence_design import (
    SCRINSHOT_or_ISS_backbone_sequence,
    convert_complementary_seq_to_arms,
    create_seqfish_plus_barcodes,
    generate_binary_sequences,
    generate_codebook,
    generate_random_sequence,
    get_barcode,
)

__all__ = [
    "BaseFtpLoader",
    "FtpLoaderNCBI",
    "FtpLoaderEnsembl",
    "GffParser",
    "check_gff_format",
    "check_fasta_format",
    "check_tsv_format",
    "get_sequence_from_annotation",
    "merge_fasta",
    "parse_fasta_header",
    "generate_random_sequence",
    "generate_binary_sequences",
    "generate_codebook",
    "get_barcode",
    "SCRINSHOT_or_ISS_backbone_sequence",
    "convert_complementary_seq_to_arms",
    "create_seqfish_plus_barcodes",
]
