"""
This module provides different approaches for a strategic selection of optimal oligo sets for genomic regions, focusing on efficiency and scoring.
"""

from ._generate_oligosets import (
    HomogeneousPropertyOligoSetGenerator,
    OligosetGeneratorIndependentSet,
)
from ._selection_methods import (
    GraphBasedSelectionPolicy,
    GreedySelectionPolicy,
    OligoSelectionPolicy,
)

__all__ = [
    "OligosetGeneratorIndependentSet",
    "HomogeneousPropertyOligoSetGenerator",
    "OligoSelectionPolicy",
    "GreedySelectionPolicy",
    "GraphBasedSelectionPolicy",
]
