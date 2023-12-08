"""
This module contains the OligoDatabase and ReferenceDatabase classes, which are used to store the oligo sequences and sequence from a reference. 
The classes implement load, save, get and write functionalities. 
"""

from ._oligos_database import OligoDatabase
from ._reference_database import ReferenceDatabase

__all__ = [
    "OligoDatabase",
    "ReferenceDatabase",
]

classes = __all__
