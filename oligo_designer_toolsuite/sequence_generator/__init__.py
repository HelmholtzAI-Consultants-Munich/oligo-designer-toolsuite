"""
This module contains the sequence generation classes, which are used to generate the sequences from specific genomic regions or the oligo sequences.
The sequence generator classes implement read, create and write functionalities. 
The genomic regions can be automatically downloaded from NCBI or Ensemble or can be provided by the user. For automatic download an FTP download class is available.
"""

from ._oligo_sequence_generator import OligoSequenceGenerator
from ._genomic_region_generator import (
    CustomGenomicRegionGenerator,
    NcbiGenomicRegionGenerator,
    EnsemblGenomicRegionGenerator,
)
from ._ftp_loader import BaseFtpLoader, FtpLoaderEnsembl, FtpLoaderNCBI

__all__ = [
    "OligoSequenceGenerator",
    "CustomGenomicRegionGenerator",
    "NcbiGenomicRegionGenerator",
    "EnsemblGenomicRegionGenerator",
    "BaseFtpLoader",
    "FtpLoaderEnsembl",
    "FtpLoaderNCBI",
]

classes = __all__
