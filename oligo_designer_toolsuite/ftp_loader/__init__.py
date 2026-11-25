"""
This module provides classes for loading genomic data via FTP.
"""

from ._ftp_loader_base import GenomicEntity
from ._ftp_loader_ensembl import FtpLoaderEnsembl
from ._ftp_loader_ncbi import FtpLoaderNCBI, NcbiGenomicEntity

__all__ = ["FtpLoaderEnsembl", "FtpLoaderNCBI", "GenomicEntity", "NcbiGenomicEntity"]

classes = __all__
