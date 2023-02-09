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
]
