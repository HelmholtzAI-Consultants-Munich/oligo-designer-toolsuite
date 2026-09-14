from typing import Annotated, Literal

from pydantic import (
    BaseModel,
    ConfigDict,
    Field,
    PositiveInt,
    model_validator,
)
from typing_extensions import Self

from oligo_designer_toolsuite.config._general_models import (
    FORWARD_PRIMER_DESC,
    OLIGO_GENERATION_DESC,
    PROBE_SET_SELECTION_DESC,
    PROPERTY_FILTERS_DESC,
    REQUIRED_PARAMETERS_DESC,
    REVERSE_PRIMER_DESC,
    SCHEMA_VERSION_DESC,
    SHARED_PARAMETERS_DESC,
    SPECIFICITY_FILTERS_DESC,
    BaseProbabilities,
    General,
    RequiredParameters,
    TmChemCorrectionParameters,
    TmParameters,
    TmSaltCorrectionParameters,
)
from oligo_designer_toolsuite.config._oligo_scoring import (
    GCContentScore,
    IndependentSetSelection,
    UTRScore,
)
from oligo_designer_toolsuite.config._property_filters import (
    HARDMASKED_DESC,
    SOFTMASKED_DESC,
    ComplementReversePrimerFilterConfig,
    GCClampFilterConfig,
    GCContentFilterConfig,
    HardMaskedFilterConfig,
    HomopolymericRunsFilterConfig,
    IsoformConsensusFilterConfig,
    SecondaryStructureFilterConfig,
    SelfComplementarityFilterConfig,
    SoftMaskedFilterConfig,
    TmFilterConfig,
)
from oligo_designer_toolsuite.config._specificity_filters import (
    SPECIFICITY_PRIMER_DESC,
    SPECIFICITY_READOUT_DESC,
    SPECIFICITY_TARGET_DESC,
    CrossHybridizationBlastnFilterConfig,
    HybridizationProbesBlastnFilterConfig,
    SpecificityBlastnFilterConfig,
)
from oligo_designer_toolsuite.config._types import (
    DRNAT,
    LengthMaxT,
    LengthMinT,
)

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


class TargetProbePropertyFilter(BaseModel):
    model_config = ConfigDict(extra="forbid")

    isoform_consensus_filter: IsoformConsensusFilterConfig
    hard_masked_sequences_filter: HardMaskedFilterConfig = Field(description=HARDMASKED_DESC)
    soft_masked_sequences_filter: SoftMaskedFilterConfig = Field(description=SOFTMASKED_DESC)
    homopolymeric_runs_filter: HomopolymericRunsFilterConfig
    GC_content_filter: GCContentFilterConfig
    secondary_structure_filter: SecondaryStructureFilterConfig


class TargetProbeSpecificityFilter(BaseModel):
    model_config = ConfigDict(extra="forbid")

    specificity_blastn_filter: SpecificityBlastnFilterConfig = Field(
        description=SPECIFICITY_TARGET_DESC,
    )
    cross_hybridization_blastn_filter: CrossHybridizationBlastnFilterConfig


class TargetProbeProbeSetSelection(BaseModel):
    model_config = ConfigDict(extra="forbid")

    independent_set_selection: IndependentSetSelection
    GC_content_score: GCContentScore
    UTR_score: UTRScore


class TargetProbes(BaseModel):
    model_config = ConfigDict(extra="forbid")

    oligo_generation: TargetProbeOligoGeneration = Field(description=OLIGO_GENERATION_DESC)
    property_filters: TargetProbePropertyFilter = Field(description=PROPERTY_FILTERS_DESC)
    specificity_filters: TargetProbeSpecificityFilter = Field(description=SPECIFICITY_FILTERS_DESC)
    probe_set_selection: TargetProbeProbeSetSelection = Field(description=PROBE_SET_SELECTION_DESC)


############################################
# Readout probes
############################################

# The codebook and readout-probe table support both source="load" and
# source="generate", modelled as discriminated unions on `source`. Both use
# extra="ignore" so the shipped YAML validates in either branch.


class SeqfishPlusCodebookBase(BaseModel):
    model_config = ConfigDict(extra="ignore")

    n_barcode_rounds: PositiveInt = Field(
        description="Number of barcode rounds. Equals active bits per gene and readout overhangs per encoding probe (first half 5' of target, second half 3').",
        json_schema_extra={"x-quick-setting": True},
    )
    n_pseudocolors: PositiveInt = Field(
        description="Number of pseudocolors per round.",
        json_schema_extra={"x-quick-setting": True},
    )
    channels_ids: list[str] = Field(
        description="Fluorescence channels used in the experiment.",
        json_schema_extra={"x-quick-setting": True},
    )


class SeqfishPlusCodebookLoad(SeqfishPlusCodebookBase):
    source: Literal["load"]
    file: str = Field(
        description="Only used when source = load. Path to the codebook file (csv/tsv): columns = 'bits', rows = 'gene_name'; entries are 0/1 bit-encodings for each gene."
    )


class SeqfishPlusCodebookGenerate(SeqfishPlusCodebookBase):
    source: Literal["generate"]


SeqfishPlusCodebook = Annotated[
    # use this order so that react-jsonschema-forms in ODT-Cloud selects the first
    # model as a default when no defaults are specified (as is currently the case)
    SeqfishPlusCodebookGenerate | SeqfishPlusCodebookLoad,
    Field(discriminator="source"),
]


class SeqfishPlusReadoutProbeTableLoad(BaseModel):
    model_config = ConfigDict(extra="ignore")

    source: Literal["load"]
    file: str = Field(
        description="Only used when source = load. Path to the bit-indexed readout probe table (csv/tsv) with columns 'barcode_round', 'pseudocolor', 'channel', and 'readout_probe_sequence'."
    )


class SeqfishPlusReadoutProbeOligoGeneration(BaseModel):
    model_config = ConfigDict(extra="forbid")

    probe_length: PositiveInt = Field(
        description="Length (bases) of each random readout probe.",
    )
    base_probabilities: BaseProbabilities = Field(
        description="Per-base probability used to generate random readout-probe sequences.",
    )
    initial_num_sequences: PositiveInt = Field(
        description="Number of random sequences to generate before filtering.",
    )


class SeqfishPlusReadoutProbePropertyFilter(BaseModel):
    model_config = ConfigDict(extra="forbid")

    homopolymeric_runs_filter: HomopolymericRunsFilterConfig
    GC_content_filter: GCContentFilterConfig


class SeqfishPlusReadoutProbeSpecificityFilter(BaseModel):
    model_config = ConfigDict(extra="forbid")

    specificity_blastn_filter: SpecificityBlastnFilterConfig = Field(
        description=SPECIFICITY_READOUT_DESC,
    )
    cross_hybridization_blastn_filter: CrossHybridizationBlastnFilterConfig


class SeqfishPlusReadoutProbeTableGenerate(BaseModel):
    model_config = ConfigDict(extra="ignore")

    source: Literal["generate"]
    oligo_generation: SeqfishPlusReadoutProbeOligoGeneration = Field(json_schema_extra={"x-collapsed": True})
    property_filters: SeqfishPlusReadoutProbePropertyFilter = Field(
        description=PROPERTY_FILTERS_DESC,
        json_schema_extra={"x-collapsed": True},
    )
    specificity_filters: SeqfishPlusReadoutProbeSpecificityFilter = Field(
        description=SPECIFICITY_FILTERS_DESC,
        json_schema_extra={"x-collapsed": True},
    )


SeqfishPlusReadoutProbeTable = Annotated[
    # use this order so that react-jsonschema-forms in ODT-Cloud selects the first
    # model as a default when no defaults are specified (as is currently the case)
    SeqfishPlusReadoutProbeTableGenerate | SeqfishPlusReadoutProbeTableLoad,
    Field(discriminator="source"),
]


class ReadoutProbes(BaseModel):
    model_config = ConfigDict(extra="forbid")

    codebook: SeqfishPlusCodebook
    readout_probe_table: SeqfishPlusReadoutProbeTable


############################################
# Primers
############################################

# The forward primer supports source="load" (provide a sequence) and
# source="generate" (design one), modelled as a discriminated union on `source`
# with extra="ignore". The reverse primer only supports source="load"
# (generate_reverse_primer raises FeatureNotImplementedError).


class SeqfishPlusForwardPrimerLoad(BaseModel):
    model_config = ConfigDict(extra="ignore")

    source: Literal["load"]
    sequence: DRNAT = Field(description="Only used when source = load.")


class SeqfishPlusForwardPrimerOligoGeneration(BaseModel):
    model_config = ConfigDict(extra="forbid")

    probe_length: PositiveInt = Field(
        description="Length (bases) of each generated primer.",
    )
    base_probabilities: BaseProbabilities = Field(
        description="Per-base probability used to generate random primer sequences.",
    )
    initial_num_sequences: PositiveInt = Field(
        description="Number of random sequences to generate before filtering.",
    )


class SeqfishPlusForwardPrimerPropertyFilter(BaseModel):
    model_config = ConfigDict(extra="forbid")

    homopolymeric_runs_filter: HomopolymericRunsFilterConfig
    GC_content_filter: GCContentFilterConfig
    GC_clamp_filter: GCClampFilterConfig
    self_complementarity_filter: SelfComplementarityFilterConfig
    complement_reverse_primer_filter: ComplementReversePrimerFilterConfig
    Tm_filter: TmFilterConfig
    secondary_structure_filter: SecondaryStructureFilterConfig


class SeqfishPlusForwardPrimerSpecificityFilter(BaseModel):
    model_config = ConfigDict(extra="forbid")

    specificity_blastn_filter: SpecificityBlastnFilterConfig = Field(
        description=SPECIFICITY_PRIMER_DESC,
    )
    hybridization_probes_blastn_filter: HybridizationProbesBlastnFilterConfig


class SeqfishPlusForwardPrimerShared(BaseModel):
    model_config = ConfigDict(extra="forbid")

    Tm_parameters: TmParameters
    Tm_chem_correction_parameters: TmChemCorrectionParameters
    Tm_salt_correction_parameters: TmSaltCorrectionParameters


class SeqfishPlusForwardPrimerGenerate(BaseModel):
    model_config = ConfigDict(extra="ignore")

    source: Literal["generate"]
    oligo_generation: SeqfishPlusForwardPrimerOligoGeneration = Field(json_schema_extra={"x-collapsed": True})
    property_filters: SeqfishPlusForwardPrimerPropertyFilter = Field(
        description=PROPERTY_FILTERS_DESC,
        json_schema_extra={"x-collapsed": True},
    )
    specificity_filters: SeqfishPlusForwardPrimerSpecificityFilter = Field(
        description=SPECIFICITY_FILTERS_DESC,
        json_schema_extra={"x-collapsed": True},
    )
    shared_parameters: SeqfishPlusForwardPrimerShared = Field(
        description=SHARED_PARAMETERS_DESC,
        json_schema_extra={"x-collapsed": True},
    )


SeqfishPlusForwardPrimer = Annotated[
    # use this order so that react-jsonschema-forms in ODT-Cloud selects the first
    # model as a default when no defaults are specified (as is currently the case)
    SeqfishPlusForwardPrimerGenerate | SeqfishPlusForwardPrimerLoad,
    Field(discriminator="source"),
]


class SeqfishPlusReversePrimer(BaseModel):
    model_config = ConfigDict(extra="forbid")

    source: Literal["load"]
    sequence: DRNAT = Field(
        description="Reverse complement of 20 nt T7 promoter sequence. Only used when source = load.",
    )


class Primers(BaseModel):
    model_config = ConfigDict(extra="forbid")

    forward_primer: SeqfishPlusForwardPrimer = Field(description=FORWARD_PRIMER_DESC)
    reverse_primer: SeqfishPlusReversePrimer = Field(description=REVERSE_PRIMER_DESC)


############################################
# Top level
############################################


# The front end builds its form from this, so `general` stays out of it.
class SeqfishPlusProbeDesignerConfigBase(BaseModel):
    model_config = ConfigDict(extra="forbid")
    schema_version: Literal[2] = Field(description=SCHEMA_VERSION_DESC)
    target_probes: TargetProbes
    readout_probes: ReadoutProbes
    primers: Primers


class SeqfishPlusProbeDesignerConfig(SeqfishPlusProbeDesignerConfigBase):
    general: General

    required_parameters: RequiredParameters = Field(description=REQUIRED_PARAMETERS_DESC)
