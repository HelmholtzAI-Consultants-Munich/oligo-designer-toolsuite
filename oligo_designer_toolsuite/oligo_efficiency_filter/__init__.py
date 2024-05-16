"""
This module assignes a score to the oligos and the oligosets for their (on-target) efficiency.
"""

from ._oligo_scoring import OligoScoringBase, TmGCOligoScoring, SeqFISHOligoScoring
from ._set_scoring import (
    AverageSetScoring,
    MaxSetScoring,
    PadlockSetScoring,
    SetScoringBase,
)

__all__ = [
    "OligoScoringBase",
    "PadlockOligoScoring",
    "SeqFISHOligoScoring",
    "SetScoringBase",
    "AverageSetScoring",
    "MaxSetScoring",
    "TmGCOligoScoring",
]
