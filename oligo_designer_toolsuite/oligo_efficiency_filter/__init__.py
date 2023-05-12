"""
This module provides different scoring approaches to calculate the efficiency of oligos.
"""


from ._oligo_scoring import OligoScoringBase, PadlockOligoScoring
from ._set_scoring import (
    SetScoringBase,
    AverageSetScoring,
    MaxSetScoring,
    PadlockSetScoring,
)

__all__ = [
    "OligoScoringBase",
    "PadlockOligoScoring",
    "SetScoringBase",
    "PadlockSetScoring",
    "AverageSetScoring",
    "MaxSetScoring",
]
