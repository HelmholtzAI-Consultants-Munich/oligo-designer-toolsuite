"""
This module provides key classes for managing and processing oligonucleotide-related databases.

Classes:
- OligoDatabase: The OligoDatabase class facilitates the management of oligonucleotide databases.
  It includes functionality for loading, saving, and manipulating oligonucleotide data.
  Users can interact with this class to handle various tasks related to oligonucleotide design and analysis.

- ReferenceDatabase: The ReferenceDatabase class is designed to handle reference databases.
  It supports the loading and retrieval of reference sequences and associated metadata.
  Users can utilize this class to incorporate reference data into their analyses and designs.

- OligoAttributes: The OligoAttributes class handles the computation of various attributes 
  associated with each oligonucleotide, such as GC content or isoform consensus.
"""

from ._oligo_database import OligoDatabase
from ._reference_database import ReferenceDatabase
from ._oligo_database_attributes import OligoAttributes

__all__ = [
    "OligoDatabase",
    "ReferenceDatabase",
    "OligoAttributes",
]

classes = __all__
