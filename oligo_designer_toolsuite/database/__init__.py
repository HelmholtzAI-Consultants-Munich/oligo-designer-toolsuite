"""
This module provides two key classes for managing oligonucleotide-related databases.

Classes:
- OligoDatabase: The OligoDatabase class facilitates the management of oligonucleotide databases.
  It includes functionality for loading, saving, and manipulating oligonucleotide data.
  Users can interact with this class to handle various tasks related to oligonucleotide design and analysis.

- ReferenceDatabase: The ReferenceDatabase class is designed to handle reference databases.
  It supports the loading and retrieval of reference sequences and associated metadata.
  Users can utilize this class to incorporate reference data into their analyses and designs.
"""

from ._oligos_database import OligoDatabase
from ._reference_database import ReferenceDatabase

__all__ = [
    "OligoDatabase",
    "ReferenceDatabase",
]

classes = __all__
