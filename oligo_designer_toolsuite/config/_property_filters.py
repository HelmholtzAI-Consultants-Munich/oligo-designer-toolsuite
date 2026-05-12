from typing import Annotated, Literal, Union

from pydantic import BaseModel, ConfigDict, Field, NonNegativeInt, PositiveInt

from oligo_designer_toolsuite.config._general_models import HomopolymerThresholds
from oligo_designer_toolsuite.config._types import (
    DNAT,
    FractionT,
    GCContentMaxT,
    GCContentMinT,
    SecondaryStructuresThresholdDeltaGT,
    TmMaxT,
    TmMinT,
    TSecondaryStructureT,
)


class FilterBaseConfig(BaseModel):
    model_config = ConfigDict(extra="forbid")


class IsoformConsensusFilterDisabled(FilterBaseConfig):
    enabled: Literal[False]


class IsoformConsensusFilterEnabled(FilterBaseConfig):
    enabled: Literal[True]
    isoform_consensus: Annotated[
        float,
        Field(
            description="Threshold for isoform consensus filtering (typically between 0.0 and 1.0) in %. Probes with isoform consensus values below this threshold will be filtered out. This ensures that selected probes target sequences that are conserved across multiple transcript isoforms.",
            ge=0,
            le=100,
        ),
    ]


IsoformConsensusFilterConfig = Annotated[
    Union[IsoformConsensusFilterEnabled, IsoformConsensusFilterDisabled], Field(discriminator="enabled")
]


class HardMaskedFilterConfig(FilterBaseConfig):
    enabled: bool


class SoftMaskedFilterConfig(FilterBaseConfig):
    enabled: bool


class HomopolymericRunsFilterDisabled(FilterBaseConfig):
    enabled: Literal[False]


class HomopolymericRunsFilterEnabled(FilterBaseConfig):
    enabled: Literal[True]
    homopolymeric_base_n: Annotated[
        HomopolymerThresholds,
        Field(description="minimum number of nucleotides to consider it a homopolymeric run per base"),
    ]


HomopolymericRunsFilterConfig = Annotated[
    Union[HomopolymericRunsFilterEnabled, HomopolymericRunsFilterDisabled], Field(discriminator="enabled")
]


class GCContentFilterDisabled(FilterBaseConfig):
    enabled: Literal[False]


class GCContentFilterEnabled(FilterBaseConfig):
    enabled: Literal[True]
    GC_content_min: GCContentMinT
    GC_content_max: GCContentMaxT


GCContentFilterConfig = Annotated[
    Union[GCContentFilterEnabled, GCContentFilterDisabled], Field(discriminator="enabled")
]


class ProhibitedSequencesFilterDisabled(FilterBaseConfig):
    enabled: Literal[False]


class ProhibitedSequencesFilterEnabled(FilterBaseConfig):
    enabled: Literal[True]
    prohibited_sequences: Annotated[
        list[DNAT],
        Field(
            description="The sequences to prohibit in the oligos. If an oligo contains any of these sequences, it will be filtered out."
        ),
    ]
    kmer_abundance_threshold: Annotated[
        dict[PositiveInt, FractionT],
        Field(
            description="The maximum abundance of a k-mer allowed. If a k-mer has an abundance greater than this threshold, all oligos containing this k-mer will be filtered out. This is a dictionary with the k-mer length as the key and the threshold as the value."
        ),
    ]


ProhibitedSequencesFilterConfig = Annotated[
    Union[ProhibitedSequencesFilterEnabled, ProhibitedSequencesFilterDisabled], Field(discriminator="enabled")
]


class SelfComplementarityFilterDisabled(FilterBaseConfig):
    enabled: Literal[False]


class SelfComplementarityFilterEnabled(FilterBaseConfig):
    enabled: Literal[True]
    max_len_selfcomplement: Annotated[
        NonNegativeInt,
        Field(
            description="Maximum allowable length of self-complementary sequences. Probes with longer self-complementary regions can form hairpins and reduce hybridization efficiency."
        ),
    ]


SelfComplementarityFilterConfig = Annotated[
    Union[SelfComplementarityFilterEnabled, SelfComplementarityFilterDisabled], Field(discriminator="enabled")
]


class MeltingTemperatureFilterDisabled(FilterBaseConfig):
    enabled: Literal[False]


class MeltingTemperatureFilterEnabled(FilterBaseConfig):
    enabled: Literal[True]
    Tm_min: TmMinT
    Tm_max: TmMaxT


MeltingTemperatureFilterConfig = Annotated[
    Union[MeltingTemperatureFilterEnabled, MeltingTemperatureFilterDisabled], Field(discriminator="enabled")
]


class SecondaryStructureFilterDisabled(FilterBaseConfig):
    enabled: Literal[False]


class SecondaryStructureFilterEnabled(FilterBaseConfig):
    enabled: Literal[True]
    T: TSecondaryStructureT
    thr_DG: SecondaryStructuresThresholdDeltaGT


SecondaryStructureFilterConfig = Annotated[
    Union[SecondaryStructureFilterEnabled, SecondaryStructureFilterDisabled], Field(discriminator="enabled")
]
