"""
This module provides different approaches for a strategic selection of optimal oligo sets for genomic regions, focusing on efficiency and scoring.

Classes:
- OligosetGeneratorIndependentSet: A class that generates sets of oligonucleotides for genomic regions, aiming to identify the most efficient combinations based on specific scoring strategies.
- HomogeneousPropertyOligoSetGenerator: A class that generates sets of oligonucleotides for genomic regions, aiming to identify the most homogeneous combinations based on specific properties.
- OligoSelectionPolicy: An abstract class that defines the interface for oligo selection policies.
- GreedySelectionPolicy: A class that defines a greedy selection policy for oligo selection.
- GraphBasedSelectionPolicy: A class that defines a graph-based selection policy for oligo selection.
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
