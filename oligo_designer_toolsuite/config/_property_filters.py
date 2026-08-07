from typing import Annotated, Literal

from pydantic import BaseModel, ConfigDict, Field, NonNegativeInt, PositiveInt, model_validator
from typing_extensions import Self

from oligo_designer_toolsuite.config._general_models import HomopolymericRunThreshold
from oligo_designer_toolsuite.config._types import (
    DRNAT,
    FractionT,
    GCContentMaxT,
    GCContentMinT,
    SecondaryStructuresThresholdDeltaGT,
    TmMaxT,
    TmMinT,
    TSecondaryStructureT,
)


class FilterBaseConfigEnabled(BaseModel):
    model_config = ConfigDict(extra="forbid")
    enabled: Literal[True]


class FilterBaseConfigDisabled(BaseModel):
    model_config = ConfigDict(extra="ignore")
    enabled: Literal[False]


class IsoformConsensusFilterDisabled(FilterBaseConfigDisabled):
    pass


class IsoformConsensusFilterEnabled(FilterBaseConfigEnabled):
    isoform_consensus: Annotated[
        float,
        Field(
            description="Threshold for isoform consensus filtering in %. Probes with isoform consensus values below this threshold will be filtered out. This ensures that selected probes target sequences that are conserved across multiple transcript isoforms.",
            ge=0,
            le=100,
        ),
    ]


IsoformConsensusFilterConfig = Annotated[
    IsoformConsensusFilterEnabled | IsoformConsensusFilterDisabled, Field(discriminator="enabled")
]


class TargetedExonsFilterDisabled(FilterBaseConfigDisabled):
    pass


class TargetedExonsFilterEnabled(FilterBaseConfigEnabled):
    targeted_exons: list[str] = Field(
        description="List of exon identifiers that should be preferentially targeted by probes, e.g. ['1', '2', '3']. Only probes in these exons are kept.",
    )


TargetedExonsFilterConfig = Annotated[
    TargetedExonsFilterEnabled | TargetedExonsFilterDisabled, Field(discriminator="enabled")
]


class HardMaskedFilterConfig(BaseModel):
    model_config = ConfigDict(extra="forbid")
    enabled: bool


class SoftMaskedFilterConfig(BaseModel):
    model_config = ConfigDict(extra="forbid")
    enabled: bool


class HomopolymericRunsFilterDisabled(FilterBaseConfigDisabled):
    pass


class HomopolymericRunsFilterEnabled(FilterBaseConfigEnabled):
    homopolymeric_base_n: Annotated[
        HomopolymericRunThreshold,
        Field(description="minimum number of nucleotides to consider it a homopolymeric run per base"),
    ]


HomopolymericRunsFilterConfig = Annotated[
    HomopolymericRunsFilterEnabled | HomopolymericRunsFilterDisabled, Field(discriminator="enabled")
]


class GCContentFilterDisabled(FilterBaseConfigDisabled):
    pass


class GCContentFilterEnabled(FilterBaseConfigEnabled):
    GC_content_min: GCContentMinT
    GC_content_max: GCContentMaxT

    @model_validator(mode="after")
    def _check_min_max(self) -> Self:
        if self.GC_content_min > self.GC_content_max:
            raise ValueError(
                f"'GC_content_min' ({self.GC_content_min}) must be <= 'GC_content_max' ({self.GC_content_max})"
            )
        return self


GCContentFilterConfig = Annotated[
    GCContentFilterEnabled | GCContentFilterDisabled, Field(discriminator="enabled")
]


class ProhibitedSequencesFilterDisabled(FilterBaseConfigDisabled):
    pass


class ProhibitedSequencesFilterEnabled(FilterBaseConfigEnabled):
    prohibited_sequences: Annotated[
        list[DRNAT],
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
    ProhibitedSequencesFilterEnabled | ProhibitedSequencesFilterDisabled, Field(discriminator="enabled")
]


class SelfComplementarityFilterDisabled(FilterBaseConfigDisabled):
    pass


class SelfComplementarityFilterEnabled(FilterBaseConfigEnabled):
    max_len_selfcomplement: Annotated[
        NonNegativeInt,
        Field(
            description="Maximum allowable length of self-complementary sequences. Probes with longer self-complementary regions can form hairpins."
        ),
    ]


SelfComplementarityFilterConfig = Annotated[
    SelfComplementarityFilterEnabled | SelfComplementarityFilterDisabled, Field(discriminator="enabled")
]


class TmFilterDisabled(FilterBaseConfigDisabled):
    pass


class TmFilterEnabled(FilterBaseConfigEnabled):
    Tm_min: TmMinT
    Tm_max: TmMaxT

    @model_validator(mode="after")
    def _check_min_max(self) -> Self:
        if self.Tm_min > self.Tm_max:
            raise ValueError(f"'Tm_min' ({self.Tm_min}) must be <= 'Tm_max' ({self.Tm_max})")
        return self


TmFilterConfig = Annotated[TmFilterEnabled | TmFilterDisabled, Field(discriminator="enabled")]


class GCClampFilterDisabled(FilterBaseConfigDisabled):
    pass


class GCClampFilterEnabled(FilterBaseConfigEnabled):
    number_GC_GCclamp: Annotated[
        PositiveInt,
        Field(description="Minimum number of G or C nucleotides required at the 3' end."),
    ]
    number_three_prime_base_GCclamp: Annotated[
        PositiveInt,
        Field(description="Number of bases from the 3' end to consider for the GC clamp."),
    ]


GCClampFilterConfig = Annotated[GCClampFilterEnabled | GCClampFilterDisabled, Field(discriminator="enabled")]


class ComplementReversePrimerFilterDisabled(FilterBaseConfigDisabled):
    pass


class ComplementReversePrimerFilterEnabled(FilterBaseConfigEnabled):
    max_len_complement: Annotated[
        NonNegativeInt,
        Field(
            description="Maximum allowable length of complementarity between the forward primer and the (known) reverse primer, to prevent primer-dimer formation."
        ),
    ]


ComplementReversePrimerFilterConfig = Annotated[
    ComplementReversePrimerFilterEnabled | ComplementReversePrimerFilterDisabled,
    Field(discriminator="enabled"),
]


class SecondaryStructureFilterDisabled(FilterBaseConfigDisabled):
    pass


class SecondaryStructureFilterEnabled(FilterBaseConfigEnabled):
    T: TSecondaryStructureT
    thr_DG: SecondaryStructuresThresholdDeltaGT


SecondaryStructureFilterConfig = Annotated[
    SecondaryStructureFilterEnabled | SecondaryStructureFilterDisabled, Field(discriminator="enabled")
]
