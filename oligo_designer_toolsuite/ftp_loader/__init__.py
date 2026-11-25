"""
This module provides classes for loading genomic data via FTP.
"""

from ._ftp_loader_base import GenomicEntity
from ._ftp_loader_ensembl import GenomicDatabaseDownloaderEnsembl
from ._ftp_loader_ncbi import GenomicDatabaseDownloaderNcbi, NcbiGenomicEntity

__all__ = [
    "GenomicDatabaseDownloaderEnsembl",
    "GenomicDatabaseDownloaderNcbi",
    "GenomicEntity",
    "NcbiGenomicEntity",
]

classes = __all__
