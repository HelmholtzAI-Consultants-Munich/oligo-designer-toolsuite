"""
This module provides different evaluation strategies for oligonucleotides and their sets based on various scoring criteria.
"""

from ._oligo_scoring import OligoScoring
from ._scorer_base import BaseScorer
from ._scorer_region_property import IsoformConsensusScorer, OverlapTargetedExonsScorer, OverlapUTRScorer
from ._scorer_sequence_property import (
    DeviationFromOptimalGCContentScorer,
    DeviationFromOptimalTmScorer,
    NormalizedDeviationFromOptimalGCContentScorer,
    NormalizedDeviationFromOptimalTmScorer,
)
from ._scorer_set_property import UniformDistanceScorer
from ._set_scoring import AverageSetScoring, LowestSetScoring, SetScoringBase

__all__ = [
    "BaseScorer",
    "OverlapTargetedExonsScorer",
    "OverlapUTRScorer",
    "IsoformConsensusScorer",
    "DeviationFromOptimalGCContentScorer",
    "DeviationFromOptimalTmScorer",
    "NormalizedDeviationFromOptimalGCContentScorer",
    "NormalizedDeviationFromOptimalTmScorer",
    "OligoScoring",
    "SetScoringBase",
    "LowestSetScoring",
    "AverageSetScoring",
    "UniformDistanceScorer",
]
