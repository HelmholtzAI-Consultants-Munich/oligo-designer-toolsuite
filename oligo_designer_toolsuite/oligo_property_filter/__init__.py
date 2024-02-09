"""
This module provides different oligo filters that assess sequence properties of oligonucelotides critical for successful applications.

Classes:
- PropertyFilterBase: An abstract base class for creating filters that assess various properties of oligonucleotides.
- PropertyFilter: Manages and applies a collection of property-based filters to a set of oligonucleotides, facilitating batch processing.
- SoftMaskedSequenceFilter: Identifies sequences with lowercase characters, indicating regions of lower confidence or variability.
- HardMaskedSequenceFilter: Filters sequences containing specific characters (e.g., "N") that denote uncertain or variable nucleotides.
- ProhibitedSequenceFilter: Excludes sequences containing specific motifs or sequences deemed unsuitable for experimental purposes.
- HomopolymericRunsFilter: Detects and filters out sequences with long runs of a single nucleotide, which can be problematic in certain assays.
- GCContentFilter: Ensures sequences fall within a specified range of GC content, critical for stable hybridization and amplification.
- GCClampFilter: Filters for sequences with a specific number of G or C nucleotides at one end, affecting primer binding and stability.
- MeltingTemperatureNNFilter: Assesses sequences for their melting temperatures using nearest-neighbor models, ensuring they are within an optimal range.
- SecondaryStructureFilter: Evaluates sequences for potential secondary structures that could interfere with assay performance.
- ThreePrimeSequenceFilter: Filters sequences based on the presence or absence of specified sequences at the 3' end.
- FivePrimeSequenceFilter: Similar to the ThreePrimeSequenceFilter, but focuses on the 5' end of sequences.
- PadlockArmsFilter: Specific to padlock probe design, it evaluates sequences for suitability based on arm length and melting temperature criteria, ensuring optimal probe performance.
"""

from ._filter_base import PropertyFilterBase
from ._property_filter import PropertyFilter
from ._filter_experiment_unspecific import (
    SoftMaskedSequenceFilter,
    HardMaskedSequenceFilter,
    ProhibitedSequenceFilter,
    HomopolymericRunsFilter,
    GCContentFilter,
    GCClampFilter,
    MeltingTemperatureNNFilter,
    SecondaryStructureFilter,
    ThreePrimeSequenceFilter,
    FivePrimeSequenceFilter,
)
from ._filter_experiment_specific import PadlockArmsFilter

__all__ = [
    "PropertyFilter",
    "PropertyFilterBase",
    "SoftMaskedSequenceFilter",
    "HardMaskedSequenceFilter",
    "ProhibitedSequenceFilter",
    "HomopolymericRunsFilter",
    "GCContentFilter",
    "GCClampFilter",
    "MeltingTemperatureNNFilter",
    "SecondaryStructureFilter",
    "ThreePrimeSequenceFilter",
    "FivePrimeSequenceFilter",
    "PadlockArmsFilter",
]
