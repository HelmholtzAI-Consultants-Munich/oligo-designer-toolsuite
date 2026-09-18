from typing import Annotated, Literal

from pydantic import (
    BaseModel,
    ConfigDict,
    Field,
    NonNegativeFloat,
    NonNegativeInt,
    PositiveInt,
    model_validator,
)
from typing_extensions import Self

from oligo_designer_toolsuite.config._general_models import (
    OLIGO_GENERATION_DESC,
    PROBE_SET_SELECTION_DESC,
    PROPERTY_FILTERS_DESC,
    REQUIRED_PARAMETERS_DESC,
    SCHEMA_VERSION_DESC,
    SHARED_PARAMETERS_DESC,
    SPECIFICITY_FILTERS_DESC,
    BlastnHitParameters,
    BlastnSearchParameters,
    General,
    RequiredParameters,
    TmChemCorrectionParameters,
    TmNNParameters,
    TmSaltCorrectionParameters,
)
from oligo_designer_toolsuite.config._oligo_scoring import (
    GCContentScoreNormalized,
    IndependentSetSelection,
    IsoformConsensusScore,
    TmScoreNormalized,
)
from oligo_designer_toolsuite.config._property_filters import (
    HARDMASKED_DESC,
    SOFTMASKED_DESC,
    GCContentFilterConfig,
    HardMaskedFilterConfig,
    HomopolymericRunsFilterConfig,
    IsoformConsensusFilterConfig,
    SoftMaskedFilterConfig,
    TmFilterConfig,
)
from oligo_designer_toolsuite.config._specificity_filters import (
    SPECIFICITY_HIT_PARAMS_DESC,
    SPECIFICITY_SEARCH_PARAMS_DESC,
    SPECIFICITY_TARGET_DESC,
    CrossHybridizationBlastnFilterConfig,
    SpecificityBlastnFilterDisabled,
    SpecificityBlastnFilterEnabled,
)
from oligo_designer_toolsuite.config._types import (
    LengthMaxT,
    LengthMinT,
    TmMaxT,
    TmMinT,
    TmOptT,
)

PADLOCK_ARMS_PROPERTIES_DESC = "Parameters that determine properties of the padlock arms."
DETECTION_OLIGO_GENERATION_DESC = "Parameters that determine length and melting temperature of the probes."

############################################
# SCRINSHOT-specific overrides
############################################


class ScrinshotSpecificityBlastnFilterEnabled(SpecificityBlastnFilterEnabled):
    search_parameters: BlastnSearchParameters = Field(
        description=SPECIFICITY_SEARCH_PARAMS_DESC,
        json_schema_extra={"x-collapsed": True},
    )
    hit_parameters: BlastnHitParameters = Field(
        description=SPECIFICITY_HIT_PARAMS_DESC,
    )
    ligation_region_size: NonNegativeInt = Field(
        description=(
            "Size of the seed region around the ligation site for BLASTN seed-region filtering. "
            "If > 0, seed-based filtering is applied around the ligation site, removing all probes"
            "where BLASTN hits cover the junction region regardless of the "
            "coverage threshold. If 0, full-length specificity filtering is performed instead."
        ),
    )


ScrinshotSpecificityBlastnFilterConfig = Annotated[
    ScrinshotSpecificityBlastnFilterEnabled | SpecificityBlastnFilterDisabled,
    Field(discriminator="enabled", description=SPECIFICITY_TARGET_DESC),
]


############################################
# Target probe
############################################


class TargetProbeOligoGeneration(BaseModel):
    model_config = ConfigDict(extra="forbid")

    probe_length_min: LengthMinT = Field(json_schema_extra={"x-quick-setting": True})
    probe_length_max: LengthMaxT = Field(json_schema_extra={"x-quick-setting": True})

    @model_validator(mode="after")
    def _check_min_max(self) -> Self:
        if self.probe_length_min > self.probe_length_max:
            raise ValueError(
                f"'probe_length_min' ({self.probe_length_min}) must be <= 'probe_length_max' ({self.probe_length_max})"
            )
        return self


class PadlockArmsProperties(BaseModel):
    model_config = ConfigDict(extra="forbid")

    length_min: PositiveInt = Field(
        description="Minimum length (bases) for each padlock arm.",
    )
    Tm_dif_max: NonNegativeFloat = Field(
        description="Maximum Tm difference (°C) between the two padlock arms.",
    )
    Tm_min: TmMinT
    Tm_max: TmMaxT

    @model_validator(mode="after")
    def _check_min_max(self) -> Self:
        if self.Tm_min > self.Tm_max:
            raise ValueError(
                f"'padlock_arm_Tm_min' ({self.Tm_min}) must be <= 'padlock_arm_Tm_max' ({self.Tm_max})"
            )
        return self


class TargetProbePropertyFilter(BaseModel):
    model_config = ConfigDict(extra="forbid")

    isoform_consensus_filter: IsoformConsensusFilterConfig
    hard_masked_sequences_filter: HardMaskedFilterConfig = Field(description=HARDMASKED_DESC)
    soft_masked_sequences_filter: SoftMaskedFilterConfig = Field(description=SOFTMASKED_DESC)
    homopolymeric_runs_filter: HomopolymericRunsFilterConfig
    GC_content_filter: GCContentFilterConfig
    Tm_filter: TmFilterConfig


class TargetProbeSpecificityFilter(BaseModel):
    model_config = ConfigDict(extra="forbid")

    specificity_blastn_filter: ScrinshotSpecificityBlastnFilterConfig
    cross_hybridization_blastn_filter: CrossHybridizationBlastnFilterConfig


class TargetProbeProbeSetSelection(BaseModel):
    model_config = ConfigDict(extra="forbid")

    independent_set_selection: IndependentSetSelection
    isoform_consensus_score: IsoformConsensusScore
    GC_content_score: GCContentScoreNormalized
    Tm_score: TmScoreNormalized


class TargetProbeTmParameters(BaseModel):
    model_config = ConfigDict(extra="forbid")

    Tm_NN_parameters: TmNNParameters
    Tm_chem_correction_parameters: TmChemCorrectionParameters
    Tm_salt_correction_parameters: TmSaltCorrectionParameters


class TargetProbes(BaseModel):
    model_config = ConfigDict(extra="forbid")

    oligo_generation: TargetProbeOligoGeneration = Field(description=OLIGO_GENERATION_DESC)
    padlock_arms_properties: PadlockArmsProperties = Field(description=PADLOCK_ARMS_PROPERTIES_DESC)
    property_filters: TargetProbePropertyFilter = Field(description=PROPERTY_FILTERS_DESC)
    specificity_filters: TargetProbeSpecificityFilter = Field(description=SPECIFICITY_FILTERS_DESC)
    probe_set_selection: TargetProbeProbeSetSelection = Field(description=PROBE_SET_SELECTION_DESC)
    Tm_parameters: TargetProbeTmParameters = Field(description=SHARED_PARAMETERS_DESC)


############################################
# Detection oligo
############################################


class DetectionOligoOligoGeneration(BaseModel):
    model_config = ConfigDict(extra="forbid")

    min_thymines: PositiveInt = Field(
        description="Minimal number of thymines (T) in the detection oligo (required for UNG cleavage after U-substitution).",
        json_schema_extra={"x-quick-setting": True},
    )
    oligo_length_min: LengthMinT = Field(json_schema_extra={"x-quick-setting": True})
    oligo_length_max: LengthMaxT = Field(json_schema_extra={"x-quick-setting": True})
    U_distance: PositiveInt = Field(
        description="Preferred minimum distance (bases) between consecutive uracils.",
        json_schema_extra={"x-quick-setting": True},
    )
    Tm_opt: TmOptT

    @model_validator(mode="after")
    def _check_min_max(self) -> Self:
        if self.oligo_length_min > self.oligo_length_max:
            raise ValueError(
                f"'oligo_length_min' ({self.oligo_length_min}) must be <= 'oligo_length_max' ({self.oligo_length_max})"
            )
        return self


class DetectionOligoTmParameters(BaseModel):
    model_config = ConfigDict(extra="forbid")

    Tm_NN_parameters: TmNNParameters
    Tm_chem_correction_parameters: TmChemCorrectionParameters
    Tm_salt_correction_parameters: TmSaltCorrectionParameters


class DetectionOligo(BaseModel):
    model_config = ConfigDict(extra="forbid")

    oligo_generation: DetectionOligoOligoGeneration = Field(description=DETECTION_OLIGO_GENERATION_DESC)
    Tm_parameters: DetectionOligoTmParameters = Field(description=SHARED_PARAMETERS_DESC)


############################################
# Top level
############################################


# The front end builds its form from this, so `general` stays out of it.
class ScrinshotProbeDesignerConfigBase(BaseModel):
    model_config = ConfigDict(extra="forbid")
    schema_version: Literal[2] = Field(description=SCHEMA_VERSION_DESC)
    target_probes: TargetProbes
    detection_oligo: DetectionOligo


class ScrinshotProbeDesignerConfig(ScrinshotProbeDesignerConfigBase):
    general: General

    required_parameters: RequiredParameters = Field(description=REQUIRED_PARAMETERS_DESC)
