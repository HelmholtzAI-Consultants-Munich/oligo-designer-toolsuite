"""
This module provides classes for generating oligo sequences, custom genomic regions,
and loading genomic data via FTP.

Classes:
- OligoSequenceGenerator: Generates oligo sequences. Handles the generation of
  random sequences based on specified probabilities. Generates sequences with
  sliding window from input sequences.

- CustomGenomicRegionGenerator: Base class for custom genomic region generation.
  Provides a framework for creating genomic region generators with specific data
  sources and configurations.

- NcbiGenomicRegionGenerator: Generates genomic regions using NCBI as the data source.
  Downloads annotation files and sequences from NCBI FTP servers.

- EnsemblGenomicRegionGenerator: Generates genomic regions using Ensembl as the data source.
  Downloads annotation files and sequences from Ensembl FTP servers.

- BaseFtpLoader: Base class for FTP loading. Handles common functionality for
  downloading files from FTP servers.

- FtpLoaderEnsembl: FTP loader for Ensembl data. Specialized class for downloading
  Ensembl annotation and sequence files.

- FtpLoaderNCBI: FTP loader for NCBI data. Specialized class for downloading NCBI
  annotation and sequence files.
"""

from ._ftp_loader import BaseFtpLoader, FtpLoaderEnsembl, FtpLoaderNCBI
from ._genomic_region_generator import (
    CustomGenomicRegionGenerator,
    EnsemblGenomicRegionGenerator,
    NcbiGenomicRegionGenerator,
)
from ._oligo_sequence_generator import OligoSequenceGenerator

__all__ = [
    "BaseFtpLoader",
    "FtpLoaderEnsembl",
    "FtpLoaderNCBI",
    "CustomGenomicRegionGenerator",
    "NcbiGenomicRegionGenerator",
    "EnsemblGenomicRegionGenerator",
    "OligoSequenceGenerator",
]

classes = __all__
