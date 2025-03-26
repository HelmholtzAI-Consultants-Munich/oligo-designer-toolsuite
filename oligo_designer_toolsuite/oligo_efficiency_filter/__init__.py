"""
This module provides different evaluation strategies for oligonucleotides and their sets based on various scoring criteria.
"""

from ._oligo_scoring import (
    GCOligoScoring,
    OligoScoringBase,
    WeightedGCUtrScoring,
    WeightedIsoformTmScoring,
    WeightedIsoformTmGCOligoScoring,
    WeightedIsoformTmGCOligoScoringTargetedExons,
    WeightedTmGCOligoScoring,
)
from ._set_scoring import AverageSetScoring, LowestSetScoring, SetScoringBase

__all__ = [
    "OligoScoringBase",
    "GCOligoScoring",
    "WeightedGCUtrScoring",
    "WeightedIsoformTmScoring",
    "WeightedTmGCOligoScoring",
    "WeightedIsoformTmGCOligoScoring",
    "WeightedIsoformTmGCOligoScoringTargetedExons",
    "SetScoringBase",
    "LowestSetScoring",
    "AverageSetScoring",
]
