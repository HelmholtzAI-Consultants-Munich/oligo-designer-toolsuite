"""
This module provides different oligo filters that assess the off-target and crosshybridization potential of oligo sequences and filter out sequences that have potential off-target regions.

Classes:
- SpecificityFilterBase: A base class for creating filters that assess the specificity of oligo sequences against target regions.
- AlignmentSpecificityFilter: Extends SpecificityFilterBase to use sequence alignment tools for specificity checking.
- ExactMatchFilter: Filters oligo sequences by detecting exact matches among other oligo sequences.
- BlastNFilter: Utilizes BLASTN for identifying and filtering sequences with potential off-target effects.
- BlastNSeedregionFilter: A BlastNFilter that focuses on the seed regions of oligo sequences for specificity filtering.
- BlastNSeedregionLigationsiteFilter: Targets the ligation site's seed regions for more precise BLASTN-based specificity filtering.
- BowtieFilter: Employs Bowtie for alignment-based specificity filtering of oligo sequences.
- Bowtie2Filter: Leverages Bowtie2, an advanced version of Bowtie, for enhanced specificity filtering.
- CrossHybridizationFilter: Applies policies to reduce cross-hybridization within oligo databases.
- RemoveByDegreePolicy: A policy for removing oligos based on their connectivity degree in cross-hybridization analysis.
- RemoveByLargerRegionPolicy: Eliminates oligos from regions with higher number of assigned oligos to manage cross-hybridization.
- SpecificityFilter: Integrates various specificity filters and policies to optimize oligo sequence design.
"""

from ._filter_base import AlignmentSpecificityFilter, SpecificityFilterBase
from ._filter_blastn import (
    BlastNFilter,
    BlastNSeedregionFilter,
    BlastNSeedregionLigationsiteFilter,
)
from ._filter_bowtie import BowtieFilter, Bowtie2Filter
from ._filter_cross_hybridization import (
    CrossHybridizationFilter,
    RemoveByDegreePolicy,
    RemoveByLargerRegionPolicy,
)
from ._filter_exact_matches import ExactMatchFilter
from ._specificity_filter import SpecificityFilter
from ._ai_filter import HybridizationProbabilityFilter


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
    "RemoveByDegreePolicy",
    "RemoveByLargerRegionPolicy",
    "SpecificityFilter",
    "HybridizationProbabilityFilter",
]
