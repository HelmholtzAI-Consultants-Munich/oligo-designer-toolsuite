"""
This module provides sequence assembly methods for different oligonucleotide sequences, e.g. assembly of the complete padlock probe sequence.
"""

from ._MERFISH_sequence import MerfishSequence
from ._padlock_sequence import PadlockSequence

# from ._merfish_sequence import MerfishSequence

__all__ = ["PadlockSequence"]
