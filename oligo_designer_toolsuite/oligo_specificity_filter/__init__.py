"""
This module provides a comprehensive set of filters designed to assess and mitigate the off-target and cross-hybridization potential of oligonucleotide sequences, ensuring high specificity and low off-target effects.
"""

from ._filter_base import SpecificityFilterBase, AlignmentSpecificityFilter

from ._filter_exact_matches import ExactMatchFilter
from ._filter_blastn import (
    BlastNFilter,
    BlastNSeedregionFilter,
    BlastNSeedregionLigationsiteFilter,
)
from ._filter_bowtie import BowtieFilter, Bowtie2Filter
from ._filter_cross_hybridization import (
    CrossHybridizationFilter,
)
from ._filter_ai import HybridizationProbabilityFilter

from ._policies import FilterPolicyBase, RemoveAllPolicy, RemoveByDegreePolicy, RemoveByLargerRegionPolicy

from ._specificity_filter import SpecificityFilter


__all__ = [
    "SpecificityFilterBase",
    "AlignmentSpecificityFilter",
    "ExactMatchFilter",
    "BlastNFilter",
    "BlastNSeedregionFilter",
    "BlastNSeedregionLigationsiteFilter",
    "BowtieFilter",
    "Bowtie2Filter",
    "CrossHybridizationFilter",
    "HybridizationProbabilityFilter",
    "FilterPolicyBase",
    "RemoveAllPolicy",
    "RemoveByDegreePolicy",
    "RemoveByLargerRegionPolicy",
    "SpecificityFilter",
]
