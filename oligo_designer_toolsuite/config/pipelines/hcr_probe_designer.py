from typing import Annotated, Literal

from pydantic import (
    BaseModel,
    ConfigDict,
    Field,
    NonNegativeInt,
    PositiveInt,
)

from oligo_designer_toolsuite.config._general_models import (
    CODEBOOK_DESC,
    INITIATOR_TABLE_DESC,
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
    IndependentSetSelection,
    IsoformConsensusScore,
)
from oligo_designer_toolsuite.config._property_filters import (
    HARDMASKED_DESC,
    SOFTMASKED_DESC,
    GCContentFilterConfig,
    HardMaskedFilterConfig,
    HomopolymericRunsFilterConfig,
    IsoformConsensusFilterConfig,
    SecondaryStructureFilterConfig,
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
from oligo_designer_toolsuite.config._types import DRNAT

############################################
# HCR-specific overrides
############################################


class HcrSpecificityBlastnFilterEnabled(SpecificityBlastnFilterEnabled):
    junction_region_size: NonNegativeInt = Field(
        description=(
            "Size of the seed region around the junction site for BLASTN seed-region filtering. "
            "If > 0, seed-based filtering is applied around the junction site, removing all probes "
            "where BLASTN hits cover the junction region regardless of the coverage threshold. "
            "If 0, full-length specificity filtering is performed instead."
        ),
    )
    search_parameters: BlastnSearchParameters = Field(
        description=SPECIFICITY_SEARCH_PARAMS_DESC,
        json_schema_extra={"x-collapsed": True},
    )
    hit_parameters: BlastnHitParameters = Field(
        description=SPECIFICITY_HIT_PARAMS_DESC,
    )


HcrSpecificityBlastnFilterConfig = Annotated[
    HcrSpecificityBlastnFilterEnabled | SpecificityBlastnFilterDisabled,
    Field(discriminator="enabled", description=SPECIFICITY_TARGET_DESC),
]


############################################
# Target probe
############################################


class TargetProbeOligoGeneration(BaseModel):
    model_config = ConfigDict(extra="forbid")

    L_probe_sequence_length: PositiveInt = Field(
        description="Length (bases) of the L arm of the probe; L + gap + R equals the total probe length.",
        json_schema_extra={"x-quick-setting": True},
    )
    gap_sequence_length: NonNegativeInt = Field(
        description="Length (bases) of the spacer between the L and R arms (covers the ligation site).",
        json_schema_extra={"x-quick-setting": True},
    )
    R_probe_sequence_length: PositiveInt = Field(
        description="Length (bases) of the R arm of the probe; L + gap + R equals the total probe length.",
        json_schema_extra={"x-quick-setting": True},
    )


class TargetProbePropertyFilter(BaseModel):
    model_config = ConfigDict(extra="forbid")

    isoform_consensus_filter: IsoformConsensusFilterConfig
    hard_masked_sequences_filter: HardMaskedFilterConfig = Field(description=HARDMASKED_DESC)
    soft_masked_sequences_filter: SoftMaskedFilterConfig = Field(description=SOFTMASKED_DESC)
    homopolymeric_runs_filter: HomopolymericRunsFilterConfig
    GC_content_filter: GCContentFilterConfig
    Tm_filter: TmFilterConfig
    secondary_structure_filter: SecondaryStructureFilterConfig


class TargetProbeSpecificityFilter(BaseModel):
    model_config = ConfigDict(extra="forbid")

    specificity_blastn_filter: HcrSpecificityBlastnFilterConfig
    cross_hybridization_blastn_filter: CrossHybridizationBlastnFilterConfig


class TargetProbeProbeSetSelection(BaseModel):
    model_config = ConfigDict(extra="forbid")

    independent_set_selection: IndependentSetSelection
    isoform_consensus_score: IsoformConsensusScore


class TargetProbeTmParameters(BaseModel):
    model_config = ConfigDict(extra="forbid")

    Tm_NN_parameters: TmNNParameters
    Tm_chem_correction_parameters: TmChemCorrectionParameters
    Tm_salt_correction_parameters: TmSaltCorrectionParameters


class TargetProbes(BaseModel):
    model_config = ConfigDict(extra="forbid")

    oligo_generation: TargetProbeOligoGeneration = Field(description=OLIGO_GENERATION_DESC)
    property_filters: TargetProbePropertyFilter = Field(description=PROPERTY_FILTERS_DESC)
    specificity_filters: TargetProbeSpecificityFilter = Field(description=SPECIFICITY_FILTERS_DESC)
    probe_set_selection: TargetProbeProbeSetSelection = Field(description=PROBE_SET_SELECTION_DESC)
    Tm_parameters: TargetProbeTmParameters = Field(description=SHARED_PARAMETERS_DESC)


############################################
# Initiator probes
############################################


class HcrCodebookLoad(BaseModel):
    model_config = ConfigDict(extra="forbid")

    source: Literal["load"]
    file: str = Field(
        description="Path to the codebook file: columns = 'bits', rows = 'gene_name'; entries are 0/1 bit-encodings for each gene."
    )


class HcrCodebookGenerate(BaseModel):
    model_config = ConfigDict(extra="ignore")

    source: Literal["generate"]


HcrCodebook = Annotated[HcrCodebookLoad | HcrCodebookGenerate, Field(discriminator="source")]


class HcrInitiatorTableLoad(BaseModel):
    model_config = ConfigDict(extra="forbid")

    source: Literal["load"]
    file: str = Field(
        description="Path to the bit-indexed initiator table (csv/tsv) with columns 'bit', 'initiator_L_sequence', and 'initiator_R_sequence'."
    )


class HcrInitiatorTableGenerate(BaseModel):
    model_config = ConfigDict(extra="ignore")

    source: Literal["generate"]


HcrInitiatorTable = Annotated[
    HcrInitiatorTableLoad | HcrInitiatorTableGenerate, Field(discriminator="source")
]


class InitiatorProbes(BaseModel):
    model_config = ConfigDict(extra="forbid")

    codebook: HcrCodebook = Field(description=CODEBOOK_DESC)
    initiator_table: HcrInitiatorTable = Field(description=INITIATOR_TABLE_DESC)


############################################
# Hybridization probes
############################################


class HybridizationProbes(BaseModel):
    model_config = ConfigDict(extra="forbid")

    linker_sequence: DRNAT = Field(
        description="Linker sequence between the initiator and the target-binding L/R arm.",
    )


############################################
# Top level
############################################


# The front end builds its form from this, so `general` stays out of it.
class HcrProbeDesignerConfigBase(BaseModel):
    model_config = ConfigDict(extra="forbid")
    schema_version: Literal[2] = Field(description=SCHEMA_VERSION_DESC)
    target_probes: TargetProbes
    initiator_probes: InitiatorProbes
    hybridization_probes: HybridizationProbes


class HcrProbeDesignerConfig(HcrProbeDesignerConfigBase):
    general: General

    required_parameters: RequiredParameters = Field(description=REQUIRED_PARAMETERS_DESC)
