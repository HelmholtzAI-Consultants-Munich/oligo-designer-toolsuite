"""
This module provides different evaluation strategies for oligonucleotides and their sets based on various scoring criteria.
"""

from ._oligo_scoring import (
    OligoScoringBase,
    GCOligoScoring,
    WeightedGCUtrScoring,
    WeightedTmGCOligoScoring,
    WeightedIsoformTmGCOligoScoring,
)
from ._set_scoring import AverageSetScoring, LowestSetScoring, SetScoringBase

__all__ = [
    "OligoScoringBase",
    "GCOligoScoring",
    "WeightedGCUtrScoring",
    "WeightedTmGCOligoScoring",
    "WeightedIsoformTmGCOligoScoring",
    "SetScoringBase",
    "LowestSetScoring",
    "AverageSetScoring",
]
