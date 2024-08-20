"""
This module provides different evaluation strategies for oligonucleotides and their sets based on various scoring criteria.

Classes:
- OligoScoringBase: An abstract base class providing a template for oligonucleotide scoring mechanisms.
- GCOligoScoring: Scores oligos based on their GC content relative to an optimal value.
- WeightedGCUtrScoring: Scores oligos based on GC content relative to an optimal value and whether the sequence originates from UTRs.
- WeightedTmGCOligoScoring: Scores oligos considering Tm and GC content with customizable weights.
- WeightedIsoformTmGCOligoScoring: Scores oligos considering isoform consensus, Tm, and GC content with customizable weights.
- SetScoringBase: An abstract base class for scoring sets of oligonucleotides.
- LowestSetScoring: Scores oligo sets based on the lowest individual oligo score, with ties broken by the sum of scores.
- AverageSetScoring: Scores oligo sets based on the average score of all oligos, with ties broken by the lowest individual score.
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
