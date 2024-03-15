"""
This module provides different approaches on how to select optimal sets of oligos for each region, e.g. sets of high scoring oligos based on oligo efficiencies.
"""

from ._generate_oligosets import OligosetGenerator
from ._heuristic_methods import padlock_heuristic_selection

__all__ = [
    "OligosetGenerator",
    "padlock_heuristic_selection",
]
