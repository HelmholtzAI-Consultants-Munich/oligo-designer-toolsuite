"""
This module provides property calculation functionality for oligonucleotides.
"""

from ._property_base import BaseProperty
from ._property_calculator import PropertyCalculator

# Export calculation functions for use by filters and scorers
from ._property_functions import (
    calc_detect_oligo,
    calc_dg_secondary_structure,
    calc_gc_content,
    calc_isoform_consensus,
    calc_length_complement,
    calc_length_selfcomplement,
    calc_num_targeted_transcripts,
    calc_oligo_length,
    calc_padlock_arms,
    calc_seedregion,
    calc_split_sequence,
    calc_tm_nn,
    calculate_reverse_complement_sequence,
    calculate_seedregion_site,
    calculate_shortened_sequence,
)
from ._property_region import IsoformConsensusProperty, NumTargetedTranscriptsProperty
from ._property_sequence import (
    DetectOligoProperty,
    DGSecondaryStructureProperty,
    GCContentProperty,
    LengthComplementProperty,
    LengthProperty,
    LengthSelfComplementProperty,
    PadlockArmsProperty,
    ReverseComplementSequenceProperty,
    SeedregionProperty,
    SeedregionSiteProperty,
    ShortenedSequenceProperty,
    SplitSequenceProperty,
    TmNNProperty,
)

__all__ = [
    "BaseProperty",
    "PropertyCalculator",
    "LengthProperty",
    "GCContentProperty",
    "TmNNProperty",
    "DGSecondaryStructureProperty",
    "LengthSelfComplementProperty",
    "LengthComplementProperty",
    "ShortenedSequenceProperty",
    "ReverseComplementSequenceProperty",
    "SplitSequenceProperty",
    "SeedregionProperty",
    "SeedregionSiteProperty",
    "PadlockArmsProperty",
    "DetectOligoProperty",
    "NumTargetedTranscriptsProperty",
    "IsoformConsensusProperty",
    # Calculation functions
    "calc_oligo_length",
    "calc_gc_content",
    "calc_tm_nn",
    "calc_dg_secondary_structure",
    "calc_length_complement",
    "calc_length_selfcomplement",
    "calculate_shortened_sequence",
    "calculate_reverse_complement_sequence",
    "calc_split_sequence",
    "calc_seedregion",
    "calculate_seedregion_site",
    "calc_padlock_arms",
    "calc_detect_oligo",
    "calc_num_targeted_transcripts",
    "calc_isoform_consensus",
]
