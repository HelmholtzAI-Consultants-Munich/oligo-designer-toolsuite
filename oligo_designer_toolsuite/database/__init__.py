"""
This module contains the databases of the oligo sequences and the reference file. In particular, it contains the classes to
generate the oligos database and the reference database and allow to read ad write them to disk.
"""

from ._oligos_database import OligoDatabase
from ._reference_database import ReferenceDatabase
from ._genomic_region_generator import (
    CustomGenomicRegionGenerator,
    NcbiGenomicRegionGenerator,
    EnsemblGenomicRegionGenerator,
)

__all__ = [
    "OligoDatabase",
    "ReferenceDatabase",
    "CustomGenomicRegionGenerator",
    "NcbiGenomicRegionGenerator",
    "EnsemblGenomicRegionGenerator",
]

classes = __all__
