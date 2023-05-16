"""
Additional functionalities to dowload form the NCBI and Ensembl servers and to deal with sequence and annotation files
"""

from ._data_parser import (
    check_gff_format,
    check_fasta_format,
    check_tsv_format,
    get_sequence_from_annotation,
    merge_fasta,
    parse_fasta_header,
)
from ._ftp_loader import BaseFtpLoader, FtpLoaderEnsembl, FtpLoaderNCBI
from ._gff_parser import GffParser
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
