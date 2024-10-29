"""
This module provides a comprehensive set of filters designed to evaluate the sequence properties of oligonucleotides, ensuring their suitability for various applications.
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
