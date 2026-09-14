from typing import Annotated, Literal

from pydantic import BaseModel, ConfigDict, Field, PositiveInt, model_validator
from typing_extensions import Self

from oligo_designer_toolsuite.config._general_models import (
    OLIGO_GENERATION_DESC,
    PROBE_SET_SELECTION_DESC,
    PROPERTY_FILTERS_DESC,
    REQUIRED_PARAMETERS_DESC,
    SCHEMA_VERSION_DESC,
    SHARED_PARAMETERS_DESC,
    SPECIFICITY_FILTERS_DESC,
    General,
    RequiredParameters,
    TmChemCorrectionParameters,
    TmParameters,
    TmSaltCorrectionParameters,
)
from oligo_designer_toolsuite.config._oligo_scoring import (
    GCContentScoreNormalized,
    IndependentSetSelection,
    IsoformConsensusScore,
    TargetedExonsScore,
    TmScoreNormalized,
    UniformDistanceScore,
)
from oligo_designer_toolsuite.config._property_filters import (
    HARDMASKED_DESC,
    SOFTMASKED_DESC,
    GCContentFilterConfig,
    HardMaskedFilterConfig,
    HomopolymericRunsFilterConfig,
    IsoformConsensusFilterConfig,
    ProhibitedSequencesFilterConfig,
    SecondaryStructureFilterConfig,
    SelfComplementarityFilterConfig,
    SoftMaskedFilterConfig,
    TargetedExonsFilterConfig,
    TmFilterConfig,
)
from oligo_designer_toolsuite.config._specificity_filters import (
    SPECIFICITY_TARGET_DESC,
    CrossHybridizationBlastnFilterConfig,
    ReadLengthBiasFilterConfig,
    SpecificityBlastnFilterConfig,
    VariantFilterDisabled,
    VariantFilterEnabled,
)
from oligo_designer_toolsuite.config._types import (
    LengthMaxT,
    LengthMinT,
)

############################################
# Oligoseq-specific overrides
############################################


class OligoSeqVariantFilterDisabled(VariantFilterDisabled):
    pass


class OligoSeqVariantFilterEnabled(VariantFilterEnabled):
    action: Annotated[
        Literal["flag", "filter"],
        Field(description="Action for variant-overlapping oligos: 'flag' (mark only) or 'filter' (exclude)."),
    ]


OligoSeqVariantFilterConfig = Annotated[
    OligoSeqVariantFilterEnabled | OligoSeqVariantFilterDisabled,
    Field(discriminator="enabled"),
]


############################################
# Target probe
############################################


class TargetProbeOligoGeneration(BaseModel):
    model_config = ConfigDict(extra="forbid")

    probe_length_min: LengthMinT = Field(json_schema_extra={"x-quick-setting": True})
    probe_length_max: LengthMaxT = Field(json_schema_extra={"x-quick-setting": True})
    probe_split_region: PositiveInt = Field(
        description="Minimum number of bases covering the exon junction, i.e. the oligo should contain at least x bases upstream/downstream of the junction.",
    )

    @model_validator(mode="after")
    def _check_min_max(self) -> Self:
        if self.probe_length_min > self.probe_length_max:
            raise ValueError(
                f"'probe_length_min' ({self.probe_length_min}) must be <= 'probe_length_max' ({self.probe_length_max})"
            )
        return self


class TargetProbePropertyFilter(BaseModel):
    model_config = ConfigDict(extra="forbid")

    isoform_consensus_filter: IsoformConsensusFilterConfig
    targeted_exons_filter: TargetedExonsFilterConfig
    hard_masked_sequences_filter: HardMaskedFilterConfig = Field(description=HARDMASKED_DESC)
    soft_masked_sequences_filter: SoftMaskedFilterConfig = Field(description=SOFTMASKED_DESC)
    homopolymeric_runs_filter: HomopolymericRunsFilterConfig
    GC_content_filter: GCContentFilterConfig
    prohibited_sequences_filter: ProhibitedSequencesFilterConfig
    self_complementarity_filter: SelfComplementarityFilterConfig
    Tm_filter: TmFilterConfig
    secondary_structure_filter: SecondaryStructureFilterConfig


# can be used to be adapted in the frontend
class TargetProbeSpecificityFilterBase(BaseModel):
    model_config = ConfigDict(extra="forbid")

    read_length_bias_filter: ReadLengthBiasFilterConfig
    specificity_blastn_filter: SpecificityBlastnFilterConfig = Field(
        description=SPECIFICITY_TARGET_DESC,
    )
    cross_hybridization_blastn_filter: CrossHybridizationBlastnFilterConfig


class TargetProbeSpecificityFilter(TargetProbeSpecificityFilterBase):
    variant_filter: OligoSeqVariantFilterConfig


class TargetProbeProbeSetSelection(BaseModel):
    model_config = ConfigDict(extra="forbid")

    independent_set_selection: IndependentSetSelection
    uniform_distance_score: UniformDistanceScore
    isoform_consensus_score: IsoformConsensusScore
    targeted_exons_score: TargetedExonsScore
    GC_content_score: GCContentScoreNormalized
    Tm_score: TmScoreNormalized


class TargetProbeShared(BaseModel):
    model_config = ConfigDict(extra="forbid")

    Tm_parameters: TmParameters
    Tm_chem_correction_parameters: TmChemCorrectionParameters
    Tm_salt_correction_parameters: TmSaltCorrectionParameters


class TargetProbes(BaseModel):
    model_config = ConfigDict(extra="forbid")

    oligo_generation: TargetProbeOligoGeneration = Field(description=OLIGO_GENERATION_DESC)
    property_filters: TargetProbePropertyFilter = Field(description=PROPERTY_FILTERS_DESC)
    specificity_filters: TargetProbeSpecificityFilter = Field(description=SPECIFICITY_FILTERS_DESC)
    probe_set_selection: TargetProbeProbeSetSelection = Field(description=PROBE_SET_SELECTION_DESC)
    shared_parameters: TargetProbeShared = Field(description=SHARED_PARAMETERS_DESC)


############################################
# Top level
############################################


# The front end builds its form from this, so `general` stays out of it.
class OligoSeqProbeDesignerConfigBase(BaseModel):
    model_config = ConfigDict(extra="forbid")
    schema_version: Literal[2] = Field(description=SCHEMA_VERSION_DESC)
    target_probes: TargetProbes


class OligoSeqProbeDesignerConfig(OligoSeqProbeDesignerConfigBase):
    general: General

    required_parameters: RequiredParameters = Field(description=REQUIRED_PARAMETERS_DESC)
