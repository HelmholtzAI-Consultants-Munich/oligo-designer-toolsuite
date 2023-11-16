"""
This module provides different sequence-property based filters for oligos like GC content or melting temperature filters.
"""

from ._filter_base import (
    ConsecutiveRepeats,
    GCClamp,
    GCContent,
    MaskedSequences,
    MeltingTemperatureNN,
    PropertyFilterBase,
    SecondaryStructure,
)
from ._filter_padlock_oligos import PadlockArms
from ._property_filter import PropertyFilter

__all__ = [
    "PropertyFilter",
    "PropertyFilterBase",
    "MaskedSequences",
    "MeltingTemperatureNN",
    "GCContent",
    "PadlockArms",
    "GCClamp",
    "ConsecutiveRepeats",
    "SecondaryStructure",
]
