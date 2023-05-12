"""
This module contains the OligoDatabase and ReferenceDatabase classes, which are used to store the oligo sequences and sequence from a reference. 
The classes implement read, create and write functionalities. This module also contains a class that can be used to generate sequences from specific genomic regions. 
The genomic regions can be automatically downloaded from NCBI or Ensemble or can be provided by the user.
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
