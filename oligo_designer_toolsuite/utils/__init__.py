"""
Additional functionalities to dowload form the NCBI and Ensembl servers and to deal with reading and writing the files
"""

from ._data_parser import (
    read_gtf,
    read_oligos_DB_gtf,
    read_oligos_DB_tsv,
    write_oligos_DB_gtf,
    write_oligos_DB_tsv,
)
from ._ftp_loader import BaseFtpLoader, FtpLoaderEnsembl, FTPLoaderNCBI

__all__ = [
    "read_oligos_DB_gtf",
    "read_oligos_DB_tsv",
    "write_oligos_DB_gtf",
    "write_oligos_DB_tsv",
    "BaseFtpLoader",
    "FTPLoaderNCBI",
    "FtpLoaderEnsembl",
    "read_gtf",
]
