"""
This module provides a comprehensive set of filters designed to evaluate the sequence properties of oligonucleotides, ensuring their suitability for various applications.

Classes:
- PropertyFilterBase: An abstract base class for creating filters that assess various properties of oligonucleotides.
- PropertyFilter: Manages and applies a collection of property-based filters to a set of oligonucleotides, facilitating batch processing.
- SoftMaskedSequenceFilter: Identifies sequences with lowercase characters, indicating regions of lower confidence or variability.
- HardMaskedSequenceFilter: Filters sequences containing specific characters (e.g., "N") that denote uncertain or variable nucleotides.
- ProhibitedSequenceFilter: Excludes sequences containing specific motifs or sequences deemed unsuitable for experimental purposes.
- HomopolymericRunsFilter: Detects and filters out sequences with long runs of a single nucleotide, which can be problematic in certain assays.
- ThreePrimeSequenceFilter: Filters sequences based on the presence or absence of specified sequences at the 3' end.
- FivePrimeSequenceFilter: Similar to the ThreePrimeSequenceFilter, but focuses on the 5' end of sequences.
- GCContentFilter: Ensures sequences fall within a specified range of GC content, critical for stable hybridization and amplification.
- GCClampFilter: Filters for sequences with a specific number of G or C nucleotides at one end, affecting primer binding and stability.
- MeltingTemperatureNNFilter: Assesses sequences for their melting temperatures using nearest-neighbor models, ensuring they are within an optimal range.
- SelfComplementFilter: Filters sequences based on their potential to form self-complementary structures, which can interfere with hybridization.
- ComplementFilter: Filters sequences based on their potential to form complementary structures with other sequences, which can interfere with hybridization.
- SecondaryStructureFilter: Evaluates sequences for potential secondary structures that could interfere with assay performance.
- PadlockArmsFilter: Specific to padlock probe design, it evaluates sequence suitability for padlock arm design based on arm length and melting temperature criteria.
- DetectionOligoFilter: Specific to padlock probe design, it evaluates sequence suitability for detection oligo design based on detection oligo length and thymine content criteria.
"""

from ._filter_base import PropertyFilterBase
from ._filter_experiment_specific import DetectionOligoFilter, PadlockArmsFilter
from ._filter_experiment_unspecific import (
    ComplementFilter,
    FivePrimeSequenceFilter,
    GCClampFilter,
    GCContentFilter,
    HardMaskedSequenceFilter,
    HomopolymericRunsFilter,
    MeltingTemperatureNNFilter,
    ProhibitedSequenceFilter,
    SecondaryStructureFilter,
    SelfComplementFilter,
    SoftMaskedSequenceFilter,
    ThreePrimeSequenceFilter,
)
from ._property_filter import PropertyFilter

__all__ = [
    "PropertyFilter",
    "PropertyFilterBase",
    "SoftMaskedSequenceFilter",
    "HardMaskedSequenceFilter",
    "ProhibitedSequenceFilter",
    "HomopolymericRunsFilter",
    "ThreePrimeSequenceFilter",
    "FivePrimeSequenceFilter",
    "GCContentFilter",
    "GCClampFilter",
    "MeltingTemperatureNNFilter",
    "SelfComplementFilter",
    "ComplementFilter",
    "SecondaryStructureFilter",
    "PadlockArmsFilter",
    "DetectionOligoFilter",
]
