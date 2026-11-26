"""
This module provides classes for loading genomic data via the FTP structure and HTTPS.
"""

from ._genomic_database_downloader_base import GenomicEntity
from ._genomic_database_downloader_ensembl import GenomicDatabaseDownloaderEnsembl
from ._genomic_database_downloader_ncbi import GenomicDatabaseDownloaderNcbi, NcbiGenomicEntity

__all__ = [
    "GenomicDatabaseDownloaderEnsembl",
    "GenomicDatabaseDownloaderNcbi",
    "GenomicEntity",
    "NcbiGenomicEntity",
]

classes = __all__
