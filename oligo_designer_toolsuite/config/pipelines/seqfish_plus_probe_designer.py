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
    BaseProbabilities,
    BlastnHitParameters,
    BlastnSearchParameters,
    General,
    HomopolymericRunThreshold,
    TmChemCorrectionParameters,
    TmChemCorrectionParametersDisabled,
    TmParameters,
    TmSaltCorrectionParameters,
    TmSaltCorrectionParametersDisabled,
)
from oligo_designer_toolsuite.config._oligo_scoring import (
    PROBE_SET_SELECTION_DESC,
    GCContentScore,
    IndependentSetSelection,
    UTRScore,
)
from oligo_designer_toolsuite.config._property_filters import (
    ComplementReversePrimerFilterConfig,
    ComplementReversePrimerFilterEnabled,
    GCClampFilterConfig,
    GCClampFilterEnabled,
    GCContentFilterConfig,
    GCContentFilterEnabled,
    HardMaskedFilterConfig,
    HomopolymericRunsFilterConfig,
    HomopolymericRunsFilterEnabled,
    IsoformConsensusFilterConfig,
    IsoformConsensusFilterEnabled,
    SecondaryStructureFilterConfig,
    SecondaryStructureFilterEnabled,
    SelfComplementarityFilterConfig,
    SelfComplementarityFilterEnabled,
    SoftMaskedFilterConfig,
    TmFilterConfig,
    TmFilterEnabled,
)
from oligo_designer_toolsuite.config._specificity_filters import (
    SPECIFICITY_PRIMER_DESC,
    SPECIFICITY_READOUT_DESC,
    SPECIFICITY_TARGET_DESC,
    CrossHybridizationBlastnFilterConfig,
    CrossHybridizationBlastnFilterEnabled,
    HybridizationProbesBlastnFilterConfig,
    HybridizationProbesBlastnFilterEnabled,
    SpecificityBlastnFilterDisabled,
    SpecificityBlastnFilterEnabled,
)
from oligo_designer_toolsuite.config._types import (
    DRNAT,
    FilesFastaDatabaseT,
    LengthMaxT,
    LengthMinT,
    RegionListT,
)

############################################
# seqFISH+-specific overrides
############################################

# The specificity BLASTN filter carries a required reference-database path
# (files_fasta_reference_database, no default). These subclasses only set the
# seqFISH+ default search/hit parameters; the path stays required, so the whole
# filter must be supplied by the user.


class SeqfishPlusSpecificityBlastnFilterDisabled(SpecificityBlastnFilterDisabled):
    pass


class SeqfishPlusTargetSpecificityBlastnFilterEnabled(SpecificityBlastnFilterEnabled):
    search_parameters: BlastnSearchParameters = BlastnSearchParameters(
        perc_identity=100,
        strand="minus",
        word_size=7,
        dust="no",
        soft_masking=False,
        max_target_seqs=10,
        max_hsps=1000,
    )
    hit_parameters: BlastnHitParameters = BlastnHitParameters(min_alignment_length=15)


SeqfishPlusTargetSpecificityBlastnFilterConfig = Annotated[
    SeqfishPlusTargetSpecificityBlastnFilterEnabled | SeqfishPlusSpecificityBlastnFilterDisabled,
    Field(discriminator="enabled", description=SPECIFICITY_TARGET_DESC),
]


class SeqfishPlusReadoutSpecificityBlastnFilterEnabled(SpecificityBlastnFilterEnabled):
    search_parameters: BlastnSearchParameters = BlastnSearchParameters(
        perc_identity=100,
        strand="minus",
        word_size=7,
        dust="no",
        soft_masking=False,
        max_target_seqs=10,
        max_hsps=1000,
    )
    hit_parameters: BlastnHitParameters = BlastnHitParameters(min_alignment_length=10)


SeqfishPlusReadoutSpecificityBlastnFilterConfig = Annotated[
    SeqfishPlusReadoutSpecificityBlastnFilterEnabled | SeqfishPlusSpecificityBlastnFilterDisabled,
    Field(discriminator="enabled", description=SPECIFICITY_READOUT_DESC),
]


class SeqfishPlusPrimerSpecificityBlastnFilterEnabled(SpecificityBlastnFilterEnabled):
    search_parameters: BlastnSearchParameters = BlastnSearchParameters(
        perc_identity=100,
        strand="minus",
        word_size=7,
        dust="no",
        soft_masking=False,
        max_target_seqs=10,
        max_hsps=1000,
    )
    hit_parameters: BlastnHitParameters = BlastnHitParameters(min_alignment_length=14)


SeqfishPlusPrimerSpecificityBlastnFilterConfig = Annotated[
    SeqfishPlusPrimerSpecificityBlastnFilterEnabled | SeqfishPlusSpecificityBlastnFilterDisabled,
    Field(discriminator="enabled", description=SPECIFICITY_PRIMER_DESC),
]


############################################
# Target probe
############################################


class TargetProbeOligoGeneration(BaseModel):
    model_config = ConfigDict(extra="forbid")

    file_region_ids: RegionListT
    files_fasta_probe_database: FilesFastaDatabaseT
    probe_length_min: LengthMinT = Field(default=28, json_schema_extra={"x-quick-setting": True})
    probe_length_max: LengthMaxT = Field(default=28, json_schema_extra={"x-quick-setting": True})

    @model_validator(mode="after")
    def _check_min_max(self) -> Self:
        if self.probe_length_min > self.probe_length_max:
            raise ValueError(
                f"'probe_length_min' ({self.probe_length_min}) must be <= 'probe_length_max' ({self.probe_length_max})"
            )
        return self


class TargetProbePropertyFilter(BaseModel):
    model_config = ConfigDict(extra="forbid")

    isoform_consensus_filter: IsoformConsensusFilterConfig = IsoformConsensusFilterEnabled(
        enabled=True, isoform_consensus=100
    )
    hard_masked_sequences_filter: HardMaskedFilterConfig = HardMaskedFilterConfig(enabled=True)
    soft_masked_sequences_filter: SoftMaskedFilterConfig = SoftMaskedFilterConfig(enabled=True)
    homopolymeric_runs_filter: HomopolymericRunsFilterConfig = HomopolymericRunsFilterEnabled(
        enabled=True, homopolymeric_base_n=HomopolymericRunThreshold(A=5, T=5, C=5, G=5)
    )
    GC_content_filter: GCContentFilterConfig = GCContentFilterEnabled(
        enabled=True, GC_content_min=45, GC_content_max=65
    )
    secondary_structure_filter: SecondaryStructureFilterConfig = SecondaryStructureFilterEnabled(
        enabled=True, T=76, thr_DG=0
    )


class TargetProbeSpecificityFilter(BaseModel):
    model_config = ConfigDict(extra="forbid")

    specificity_blastn_filter: SeqfishPlusTargetSpecificityBlastnFilterConfig
    cross_hybridization_blastn_filter: CrossHybridizationBlastnFilterConfig = (
        CrossHybridizationBlastnFilterEnabled(
            enabled=True,
            search_parameters=BlastnSearchParameters(
                perc_identity=80,
                strand="minus",
                word_size=7,
                dust="no",
                soft_masking=False,
                max_target_seqs=10,
            ),
            hit_parameters=BlastnHitParameters(min_alignment_length=17),
        )
    )


class TargetProbeProbeSetSelection(BaseModel):
    model_config = ConfigDict(extra="forbid")

    independent_set_selection: IndependentSetSelection = IndependentSetSelection(
        n_sets=3,
        set_size_min=24,
        set_size_opt=24,
        distance_between_probes=2,
        n_attempts_graph=50,
        n_attempts_clique_enum=50,
        diversification_fraction=0.1,
        jaccard_opt=0.5,
        jaccard_step=0.1,
    )
    GC_content_score: GCContentScore = GCContentScore(weight=1, GC_content_opt=55)
    UTR_score: UTRScore = UTRScore(weight=10)


class TargetProbes(BaseModel):
    model_config = ConfigDict(extra="forbid")

    oligo_generation: TargetProbeOligoGeneration
    property_filters: TargetProbePropertyFilter
    specificity_filters: TargetProbeSpecificityFilter
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
        default=4,
        json_schema_extra={"x-quick-setting": True},
    )
    n_pseudocolors: PositiveInt = Field(
        description="Number of pseudocolors per round.",
        default=4,
        json_schema_extra={"x-quick-setting": True},
    )
    channels_ids: list[str] = Field(
        description="Fluorescence channels used in the experiment.",
        default=["Alexa488", "Cy3b", "Alexa647"],
        json_schema_extra={"x-quick-setting": True},
    )


class SeqfishPlusCodebookLoad(SeqfishPlusCodebookBase):

    source: Literal["load"]
    file: str = Field(
        description="Only used when source = load. Path to the codebook file (csv/tsv): columns = 'bits', rows = 'gene_name'; entries are 0/1 bit-encodings for each gene."
    )


class SeqfishPlusCodebookGenerate(SeqfishPlusCodebookBase):

    source: Literal["generate"] = "generate"


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
        default=15,
    )
    base_probabilities: BaseProbabilities = Field(
        default=BaseProbabilities(A=0.25, C=0.25, G=0.25, T=0.25),
        description="Per-base probability used to generate random readout-probe sequences.",
    )
    initial_num_sequences: PositiveInt = Field(
        description="Number of random sequences to generate before filtering.",
        default=100000,
    )


class SeqfishPlusReadoutProbePropertyFilter(BaseModel):
    model_config = ConfigDict(extra="forbid")

    homopolymeric_runs_filter: HomopolymericRunsFilterConfig = HomopolymericRunsFilterEnabled(
        enabled=True, homopolymeric_base_n=HomopolymericRunThreshold(G=3)
    )
    GC_content_filter: GCContentFilterConfig = GCContentFilterEnabled(
        enabled=True, GC_content_min=40, GC_content_max=60
    )


class SeqfishPlusReadoutProbeSpecificityFilter(BaseModel):
    model_config = ConfigDict(extra="forbid")

    specificity_blastn_filter: SeqfishPlusReadoutSpecificityBlastnFilterConfig
    cross_hybridization_blastn_filter: CrossHybridizationBlastnFilterConfig = (
        CrossHybridizationBlastnFilterEnabled(
            enabled=True,
            search_parameters=BlastnSearchParameters(
                perc_identity=100,
                strand="minus",
                word_size=7,
                dust="no",
                soft_masking=False,
                max_target_seqs=10,
            ),
            hit_parameters=BlastnHitParameters(min_alignment_length=10),
        )
    )


class SeqfishPlusReadoutProbeTableGenerate(BaseModel):
    model_config = ConfigDict(extra="ignore")

    source: Literal["generate"] = "generate"
    oligo_generation: SeqfishPlusReadoutProbeOligoGeneration
    property_filters: SeqfishPlusReadoutProbePropertyFilter
    specificity_filters: SeqfishPlusReadoutProbeSpecificityFilter


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
        default=20,
    )
    base_probabilities: BaseProbabilities = Field(
        default=BaseProbabilities(A=0.25, C=0.25, G=0.25, T=0.25),
        description="Per-base probability used to generate random primer sequences.",
    )
    initial_num_sequences: PositiveInt = Field(
        description="Number of random sequences to generate before filtering.",
        default=1000000,
    )


class SeqfishPlusForwardPrimerPropertyFilter(BaseModel):
    model_config = ConfigDict(extra="forbid")

    homopolymeric_runs_filter: HomopolymericRunsFilterConfig = HomopolymericRunsFilterEnabled(
        enabled=True, homopolymeric_base_n=HomopolymericRunThreshold(A=4, T=4, C=4, G=4)
    )
    GC_content_filter: GCContentFilterConfig = GCContentFilterEnabled(
        enabled=True, GC_content_min=50, GC_content_max=65
    )
    GC_clamp_filter: GCClampFilterConfig = GCClampFilterEnabled(
        enabled=True, number_GC_GCclamp=1, number_three_prime_base_GCclamp=2
    )
    self_complementarity_filter: SelfComplementarityFilterConfig = SelfComplementarityFilterEnabled(
        enabled=True, max_len_selfcomplement=6
    )
    complement_reverse_primer_filter: ComplementReversePrimerFilterConfig = (
        ComplementReversePrimerFilterEnabled(enabled=True, max_len_complement=5)
    )
    Tm_filter: TmFilterConfig = TmFilterEnabled(enabled=True, Tm_min=60, Tm_max=75)
    secondary_structure_filter: SecondaryStructureFilterConfig = SecondaryStructureFilterEnabled(
        enabled=True, T=76, thr_DG=0
    )


class SeqfishPlusForwardPrimerSpecificityFilter(BaseModel):
    model_config = ConfigDict(extra="forbid")

    specificity_blastn_filter: SeqfishPlusPrimerSpecificityBlastnFilterConfig
    hybridization_probes_blastn_filter: HybridizationProbesBlastnFilterConfig = (
        HybridizationProbesBlastnFilterEnabled(
            enabled=True,
            search_parameters=BlastnSearchParameters(
                perc_identity=100,
                strand="minus",
                word_size=7,
                dust="no",
                soft_masking=False,
                max_target_seqs=10,
                max_hsps=1000,
            ),
            hit_parameters=BlastnHitParameters(min_alignment_length=11),
        )
    )


class SeqfishPlusForwardPrimerGlobal(BaseModel):
    model_config = ConfigDict(extra="forbid")

    Tm_parameters: TmParameters = TmParameters(
        nn_table="DNA_NN4",
        tmm_table="DNA_TMM1",
        imm_table="DNA_IMM1",
        de_table="DNA_DE1",
        dnac1=250,
        dnac2=250,
        saltcorr=5,
        Na=300,
        K=0,
        Tris=0,
        Mg=0,
        dNTPs=0,
    )
    Tm_chem_correction_parameters: TmChemCorrectionParameters = TmChemCorrectionParametersDisabled(
        enabled=False
    )
    Tm_salt_correction_parameters: TmSaltCorrectionParameters = TmSaltCorrectionParametersDisabled(
        enabled=False
    )


class SeqfishPlusForwardPrimerGenerate(BaseModel):
    model_config = ConfigDict(extra="ignore")

    source: Literal["generate"] = "generate"
    oligo_generation: SeqfishPlusForwardPrimerOligoGeneration
    property_filters: SeqfishPlusForwardPrimerPropertyFilter
    specificity_filters: SeqfishPlusForwardPrimerSpecificityFilter
    global_parameters: SeqfishPlusForwardPrimerGlobal


SeqfishPlusForwardPrimer = Annotated[
    # use this order so that react-jsonschema-forms in ODT-Cloud selects the first
    # model as a default when no defaults are specified (as is currently the case)
    SeqfishPlusForwardPrimerGenerate | SeqfishPlusForwardPrimerLoad,
    Field(discriminator="source"),
]


class SeqfishPlusReversePrimer(BaseModel):
    model_config = ConfigDict(extra="forbid")

    source: Literal["load"] = "load"
    sequence: DRNAT = Field(
        description="Reverse complement of 20 nt T7 promoter sequence. Only used when source = load.",
        default="CCCTATAGTGAGTCGTATTA",
    )


class Primers(BaseModel):
    model_config = ConfigDict(extra="forbid")

    forward_primer: SeqfishPlusForwardPrimer
    reverse_primer: SeqfishPlusReversePrimer


############################################
# Top level
############################################


class SeqfishPlusProbeDesignerConfig(BaseModel):
    model_config = ConfigDict(extra="forbid")
    schema_version: Literal[2] = 2
    general: General = General(
        n_jobs=4,
        dir_output="output_SeqfishPlusplus_probe_designer",
        write_intermediate_steps=True,
    )

    target_probes: TargetProbes
    readout_probes: ReadoutProbes
    primers: Primers
