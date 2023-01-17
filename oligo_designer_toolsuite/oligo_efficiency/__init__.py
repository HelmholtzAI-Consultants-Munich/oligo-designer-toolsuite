"""
This module assignes a score to the oligos and the oligosets for their (on-target) efficiency.
"""


from ._oligo_scoring import OligoScoringBase, PadlockOligoScoring
from ._set_scoring import (
    AverageSetScoring,
    MaxSetScoring,
    PadlockSetScoring,
    SetScoringBase,
)

__all__ = [
    "OligoScoringBase",
    "PadlockOligoScoring",
    "SetScoringBase",
    "PadlockSetScoring",
    "AverageSetScoring",
    "MaxSetScoring",
]
