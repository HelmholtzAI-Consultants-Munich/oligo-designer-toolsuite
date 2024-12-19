"""
This module provides different evaluation strategies for oligonucleotides and their sets based on their on-target efficiency.

Classes:
- OligoScoringBase: This abstract base class provides a template for implementing scoring mechanisms that evaluate the efficiency of oligonucleotides based on specific criteria.
- GCOligoScoring: A class for scoring oligos considering the distance of the oligos GC content to an optimal user-defined GC content.
- WeightedGCUtrScoring: A class for scoring oligos considering the distance of the oligos GC content to an optimal user-defined GC content and whether the sequence originates from untranslated regions (UTRs).
- WeightedTmGCOligoScoring: A class for scoring oligos considering melting temperature (Tm), and GC content, with customizable weights for each factor.
- WeightedIsoformTmGCOligoScoring: A class for scoring oligos considering isoform consensus, melting temperature (Tm), and GC content, with customizable weights for each factor.
- SetScoringBase: An abstract base class designed for scoring sets of oligonucleotides, facilitating the selection of optimal oligo sets based on collective properties.
- LowestSetScoring: A scoring strategy that selects sets of oligonucleotides by scoring the set by it's lowest oligo score, or in case of ties by the sum all oligo scores within the set.
- AverageSetScoring: A scoring strategy that selects sets of oligonucleotides by scoring the set by the average of all oligo scores within the set, or in case of ties by the lowest oligo score.
"""

from ._oligo_scoring import (
    OligoScoringBase,
    GCOligoScoring,
    WeightedGCUtrScoring,
    WeightedTmGCOligoScoring,
    WeightedIsoformTmGCOligoScoring,
)
from ._set_scoring import (
    SetScoringBase,
    LowestSetScoring,
    AverageSetScoring,
)

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
