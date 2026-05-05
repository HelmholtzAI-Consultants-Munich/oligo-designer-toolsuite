"""
This module provides different approaches for a strategic selection of optimal oligo sets for genomic regions, focusing on efficiency and scoring.
"""

from ._oligo_selection_base import BaseOligoSelection
from ._oligo_selection_homogeneous_properties import HomogeneousPropertyOligoSelection
from ._oligo_selection_independent_sets import IndependentSetsOligoSelection

__all__ = [
    "BaseOligoSelection",
    "IndependentSetsOligoSelection",
    "HomogeneousPropertyOligoSelection",
]
