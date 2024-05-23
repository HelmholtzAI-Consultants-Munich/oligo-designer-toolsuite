"""
This module provides different evaluation strategies for oligonucleotides and their sets based on their on-target efficiency.

Classes:
- OligoScoringBase: This abstract base class provides a template for implementing scoring mechanisms that evaluate the efficiency of oligonucleotides based on specific criteria.
- WeightedTmGCOligoScoring: A class that scores oligos by calculating a weighted score based on their melting temperature (Tm) and GC content, allowing for optimal oligo design based on thermal stability and nucleotide composition.
- GCOligoScoring: This class scores oligonucleotides based purely on their GC content wrt. certain GC content ranges.
- SetScoringBase: An abstract base class designed for scoring sets of oligonucleotides, facilitating the selection of optimal oligo sets based on collective properties.
- LowestSetScoring: A scoring strategy that selects sets of oligonucleotides by scoring the set by it's lowest oligo score, or in case of ties by the sum all oligo scores within the set.
- AverageSetScoring: A scoring strategy that selects sets of oligonucleotides by scoring the set by the average of all oligo scores within the set, or in case of ties by the lowest oligo score.
"""

from ._oligo_scoring import OligoScoringBase, WeightedTmGCOligoScoring, GCOligoScoring
from ._set_scoring import (
    SetScoringBase,
    LowestSetScoring,
    AverageSetScoring,
)

__all__ = [
    "OligoScoringBase",
    "WeightedTmGCOligoScoring",
    "GCOligoScoring",
    "SetScoringBase",
    "LowestSetScoring",
    "AverageSetScoring",
]
