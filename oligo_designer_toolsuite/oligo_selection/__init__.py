"""
This module provides different approaches for a strategic selection of optimal oligo sets for genomic regions, focusing on efficiency and scoring.
"""

from ._generate_oligosets import OligosetGeneratorIndependentSet
from ._heuristic_selection_methods import heuristic_selection_independent_set

__all__ = [
    "OligosetGeneratorIndependentSet",
    "heuristic_selection_independent_set",
]
