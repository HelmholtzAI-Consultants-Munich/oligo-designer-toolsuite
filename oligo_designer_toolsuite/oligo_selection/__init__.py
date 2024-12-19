"""
This module provides different approaches for a strategic selection of optimal oligo sets for genomic regions, focusing on efficiency and scoring.

Classes:
- OligosetGeneratorIndependentSet: A class that generates sets of oligonucleotides for genomic regions, aiming to identify the most efficient combinations based on specific scoring strategies.

Functions:
- heuristic_selection_independent_set: A function that applies heuristic methods to select optimal oligos, enhancing the selection process by focusing on non-overlapping and high-efficiency criteria.
"""

from ._generate_oligosets import OligosetGeneratorIndependentSet
from ._heuristic_selection_methods import heuristic_selection_independent_set

__all__ = [
    "OligosetGeneratorIndependentSet",
    "heuristic_selection_independent_set",
]
