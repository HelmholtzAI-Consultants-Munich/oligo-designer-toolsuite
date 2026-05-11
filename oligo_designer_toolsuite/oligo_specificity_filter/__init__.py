"""
This module provides a comprehensive set of filters designed to assess and mitigate the off-target and cross-hybridization potential of oligonucleotide sequences, ensuring high specificity and low off-target effects.
"""

from ._filter_base import AlignmentSpecificityFilter, BaseSpecificityFilter, ReferenceSpecificityFilter
from ._filter_blastn import BlastNFilter, BlastNSeedregionFilter, BlastNSeedregionSiteFilter
from ._filter_bowtie import Bowtie2Filter, BowtieFilter
from ._filter_cross_hybridization import CrossHybridizationFilter
from ._filter_exact_matches import ExactMatchFilter
from ._filter_hybridization_probability import HybridizationProbabilityFilter
from ._filter_variants import VariantsFilter
from ._policies import (
    BaseFilterPolicy,
    RemoveAllFilterPolicy,
    RemoveByDegreeFilterPolicy,
    RemoveByLargerRegionFilterPolicy,
)
from ._specificity_filter import SpecificityFilter

__all__ = [
    "BaseSpecificityFilter",
    "ReferenceSpecificityFilter",
    "AlignmentSpecificityFilter",
    "ExactMatchFilter",
    "BlastNFilter",
    "BlastNSeedregionFilter",
    "BlastNSeedregionSiteFilter",
    "BowtieFilter",
    "Bowtie2Filter",
    "CrossHybridizationFilter",
    "HybridizationProbabilityFilter",
    "VariantsFilter",
    "BaseFilterPolicy",
    "RemoveAllFilterPolicy",
    "RemoveByDegreeFilterPolicy",
    "RemoveByLargerRegionFilterPolicy",
    "SpecificityFilter",
]
