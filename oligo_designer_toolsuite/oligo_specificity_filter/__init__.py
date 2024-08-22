"""
This module provides a comprehensive set of filters designed to assess and mitigate the off-target and cross-hybridization potential of oligonucleotide sequences, ensuring high specificity and low off-target effects.

Classes:
- SpecificityFilterBase: Abstract base class for creating oligo specificity filters.
- AlignmentSpecificityFilter: Extends SpecificityFilterBase, utilizing sequence alignment tools like BLAST and Bowtie for specificity checks.
- ExactMatchFilter: Filters out oligos with exact matches to others within a database.
- BlastNFilter: Uses BLASTN for filtering sequences with potential off-target regions.
- BlastNSeedregionFilter: Focuses on seed regions during BLASTN-based specificity filtering.
- BlastNSeedregionLigationsiteFilter: Focuses on seed regions around ligation sites during BLASTN-based specificity filtering.
- BowtieFilter: Uses Bowtie for filtering sequences with potential off-target regions.
- Bowtie2Filter: Uses Bowtie2 for filtering sequences with potential off-target regions.
- CrossHybridizationFilter: Reduces cross-hybridization within oligo databases.
- HybridizationProbabilityFilter: Enhances specificity filters with AI-driven hybridization probability models.
- FilterPolicyBase: Base class for implementing filtering policies to manage cross-hybridization.
- RemoveAllPolicy: Removes all oligos that exhibit potential hybridization with off-target sequences.
- RemoveByDegreePolicy: Removes oligos based on their degree of connectivity in hybridization analysis.
- RemoveByLargerRegionPolicy: Eliminates oligos from regions with a higher oligo count to manage hybridization effectively.
- SpecificityFilter: Integrates multiple specificity filters and policies, offering a comprehensive tool for optimizing oligo sequence design.
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
