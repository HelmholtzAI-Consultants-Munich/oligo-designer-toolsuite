from typing import Annotated, Literal

from pydantic import (
    BaseModel,
    ConfigDict,
    Field,
    NonNegativeInt,
    PositiveInt,
)

from oligo_designer_toolsuite.config._general_models import (
    GLOBAL_PARAMETERS_DESC,
    OLIGO_GENERATION_DESC,
    PROBE_SET_SELECTION_DESC,
    PROPERTY_FILTERS_DESC,
    REQUIRED_PARAMETERS_DESC,
    SPECIFICITY_FILTERS_DESC,
    BlastnHitParameters,
    BlastnSearchParameters,
    General,
    RequiredParameters,
    TmChemCorrectionParameters,
    TmParameters,
    TmSaltCorrectionParameters,
)
from oligo_designer_toolsuite.config._oligo_scoring import (
    IndependentSetSelection,
    IsoformConsensusScore,
    TmScore,
)
from oligo_designer_toolsuite.config._property_filters import (
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
    CrossHybridizationBlastnFilterCoverageConfig,
    SpecificityBlastnFilterDisabled,
    SpecificityBlastnFilterEnabled,
)
from oligo_designer_toolsuite.config._types import DRNAT

############################################
# CycleHCR-specific overrides
############################################


class CycleHcrSpecificityBlastnFilterEnabled(SpecificityBlastnFilterEnabled):
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


CycleHcrSpecificityBlastnFilterConfig = Annotated[
    CycleHcrSpecificityBlastnFilterEnabled | SpecificityBlastnFilterDisabled,
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
        description="Length (bases) of the spacer between the L and R arms (covers the junction site).",
        json_schema_extra={"x-quick-setting": True},
    )
    R_probe_sequence_length: PositiveInt = Field(
        description="Length (bases) of the R arm of the probe; L + gap + R equals the total probe length.",
        json_schema_extra={"x-quick-setting": True},
    )


class TargetProbePropertyFilter(BaseModel):
    model_config = ConfigDict(extra="forbid")

    isoform_consensus_filter: IsoformConsensusFilterConfig
    hard_masked_sequences_filter: HardMaskedFilterConfig
    soft_masked_sequences_filter: SoftMaskedFilterConfig
    homopolymeric_runs_filter: HomopolymericRunsFilterConfig
    GC_content_filter: GCContentFilterConfig
    Tm_filter: TmFilterConfig
    secondary_structure_filter: SecondaryStructureFilterConfig


class TargetProbeSpecificityFilter(BaseModel):
    model_config = ConfigDict(extra="forbid")

    specificity_blastn_filter: CycleHcrSpecificityBlastnFilterConfig
    cross_hybridization_blastn_filter: CrossHybridizationBlastnFilterCoverageConfig


class TargetProbeProbeSetSelection(BaseModel):
    model_config = ConfigDict(extra="forbid")

    independent_set_selection: IndependentSetSelection
    isoform_consensus_score: IsoformConsensusScore
    Tm_score: TmScore


class TargetProbeGlobal(BaseModel):
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
    global_parameters: TargetProbeGlobal = Field(description=GLOBAL_PARAMETERS_DESC)


############################################
# Readout probes
############################################

# The codebook supports both source="load" and source="generate" for CycleHCR. `extra="ignore"`
# lets the shipped YAML carry both `file` and `min_hamming_distance` under `codebook`.
# The readout probe table only supports source="load" (generate_readout_probe_table raises
# FeatureNotImplementedError), so `source` is modelled as Literal["load"] to reject "generate".


class CycleHcrCodebookLoad(BaseModel):
    model_config = ConfigDict(extra="ignore")

    source: Literal["load"]
    file: str = Field(
        description="Path to the codebook file (csv/tsv): columns = 'bits', rows = 'gene_name'; entries are 0/1 bit-encodings for each gene. An example file can be found at: data/functional_oligos/cycle_hcr_codebook.tsv; codebook file."
    )


class CycleHcrCodebookGenerate(BaseModel):
    model_config = ConfigDict(extra="ignore")

    source: Literal["generate"] = "generate"
    min_hamming_distance: Literal[0, 2, 4] = Field(
        description=(
            "Required minimum Hamming distance between codewords. Codewords have weight 2, so only "
            "0, 2, or 4 are achievable. 4 enables single-bit error detection but limits capacity to "
            "n_readout_probes_LR * n_channels regions."
        ),
    )


CycleHcrCodebook = Annotated[
    CycleHcrCodebookLoad | CycleHcrCodebookGenerate,
    Field(discriminator="source"),
]


class CycleHcrReadoutProbeTable(BaseModel):
    model_config = ConfigDict(extra="forbid")

    source: Literal["load"]
    file: str = Field(
        description="Path to the readout probe table (csv/tsv) with columns 'channel', 'readout_probe_id', 'readout_probe_sequence', and 'L/R'. Bit handling depends on codebook.source: when codebook.source = 'load' the file MUST also contain a 'bit' column whose values match the codebook columns (the user is responsible for that mapping); when codebook.source = 'generate' any 'bit' column is ignored and bits are reassigned deterministically by sorting on (readout_probe_id, channel, L/R)."
    )


class ReadoutProbes(BaseModel):
    model_config = ConfigDict(extra="forbid")

    codebook: CycleHcrCodebook
    readout_probe_table: CycleHcrReadoutProbeTable


############################################
# Primers
############################################

# Primer generation is not implemented (generate_forward_primer/generate_reverse_primer raise
# FeatureNotImplementedError), so `source` is modelled as Literal["load"] to reject "generate".


class CycleHcrPrimer(BaseModel):
    model_config = ConfigDict(extra="forbid")

    source: Literal["load"]
    sequence: DRNAT


class Primers(BaseModel):
    model_config = ConfigDict(extra="forbid")

    forward_primer: CycleHcrPrimer = Field(
        description="Forward PCR primer placed at the 5' end of the DNA template probe. The default is the T7 promoter sequence.",
    )
    reverse_primer: CycleHcrPrimer = Field(
        description="Reverse PCR primer placed at the 3' end of the DNA template probe.",
    )


############################################
# Hybridization probes
############################################


class HybridizationProbes(BaseModel):
    model_config = ConfigDict(extra="forbid")

    linker_sequence: DRNAT = Field(
        description="Linker sequence between the target-binding L/R arm and the readout-probe barcode (used by both hybridization and DNA template assembly).",
    )


############################################
# Top level
############################################


# The front end builds its form from this, so `general` stays out of it.
class CycleHcrProbeDesignerConfigBase(BaseModel):
    model_config = ConfigDict(extra="forbid")
    schema_version: Literal[2]
    target_probes: TargetProbes
    readout_probes: ReadoutProbes
    primers: Primers
    hybridization_probes: HybridizationProbes


class CycleHcrProbeDesignerConfig(CycleHcrProbeDesignerConfigBase):
    general: General

    required_parameters: RequiredParameters = Field(description=REQUIRED_PARAMETERS_DESC)
