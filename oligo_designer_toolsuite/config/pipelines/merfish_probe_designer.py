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
    IsoformConsensusScore,
    TmScore,
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
    HomogeneousPropertiesWeightsT,
    LengthMaxT,
    LengthMinT,
    RegionListT,
)

############################################
# MERFISH-specific overrides
############################################

# The specificity BLASTN filter carries a required reference-database path
# (files_fasta_reference_database, no default). These subclasses only set the
# MERFISH default search/hit parameters; the path stays required, so the whole
# filter must be supplied by the user.


class MerfishSpecificityBlastnFilterDisabled(SpecificityBlastnFilterDisabled):
    pass


class MerfishTargetSpecificityBlastnFilterEnabled(SpecificityBlastnFilterEnabled):
    search_parameters: BlastnSearchParameters = Field(
        default=BlastnSearchParameters(
            perc_identity=80,
            strand="minus",
            word_size=10,
            dust="no",
            soft_masking=False,
            max_target_seqs=10,
            max_hsps=1000,
        ),
        description="BLAST options for the specificity search.",
    )
    hit_parameters: BlastnHitParameters = Field(
        default=BlastnHitParameters(min_alignment_length=17),
        description="Hit criteria. Hits satisfying these lead to oligo rejection.",
    )


MerfishTargetSpecificityBlastnFilterConfig = Annotated[
    MerfishTargetSpecificityBlastnFilterEnabled | MerfishSpecificityBlastnFilterDisabled,
    Field(discriminator="enabled", description=SPECIFICITY_TARGET_DESC),
]


class MerfishReadoutSpecificityBlastnFilterEnabled(SpecificityBlastnFilterEnabled):
    search_parameters: BlastnSearchParameters = BlastnSearchParameters(
        perc_identity=100,
        strand="minus",
        word_size=7,
        dust="no",
        soft_masking=False,
        max_target_seqs=10,
        max_hsps=1000,
    )
    hit_parameters: BlastnHitParameters = BlastnHitParameters(min_alignment_length=11)


MerfishReadoutSpecificityBlastnFilterConfig = Annotated[
    MerfishReadoutSpecificityBlastnFilterEnabled | MerfishSpecificityBlastnFilterDisabled,
    Field(discriminator="enabled"),
]


class MerfishPrimerSpecificityBlastnFilterEnabled(SpecificityBlastnFilterEnabled):
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


MerfishPrimerSpecificityBlastnFilterConfig = Annotated[
    MerfishPrimerSpecificityBlastnFilterEnabled | MerfishSpecificityBlastnFilterDisabled,
    Field(discriminator="enabled"),
]


############################################
# Target probe
############################################


class TargetProbeOligoGeneration(BaseModel):
    model_config = ConfigDict(extra="forbid")

    file_region_ids: RegionListT
    files_fasta_probe_database: FilesFastaDatabaseT
    probe_length_min: LengthMinT = Field(default=30, json_schema_extra={"x-quick-setting": True})
    probe_length_max: LengthMaxT = Field(default=30, json_schema_extra={"x-quick-setting": True})

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
        enabled=True, isoform_consensus=50
    )
    hard_masked_sequences_filter: HardMaskedFilterConfig = HardMaskedFilterConfig(enabled=True)
    soft_masked_sequences_filter: SoftMaskedFilterConfig = SoftMaskedFilterConfig(enabled=True)
    homopolymeric_runs_filter: HomopolymericRunsFilterConfig = HomopolymericRunsFilterEnabled(
        enabled=True, homopolymeric_base_n=HomopolymericRunThreshold(A=6, T=6, C=6, G=6)
    )
    GC_content_filter: GCContentFilterConfig = GCContentFilterEnabled(
        enabled=True, GC_content_min=43, GC_content_max=63
    )
    Tm_filter: TmFilterConfig = TmFilterEnabled(enabled=True, Tm_min=66, Tm_max=76)
    secondary_structure_filter: SecondaryStructureFilterConfig = SecondaryStructureFilterEnabled(
        enabled=True, T=76, thr_DG=0
    )


class TargetProbeSpecificityFilter(BaseModel):
    model_config = ConfigDict(extra="forbid")

    specificity_blastn_filter: MerfishTargetSpecificityBlastnFilterConfig
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
        set_size_min=50,
        set_size_opt=50,
        distance_between_probes=0,
        n_attempts_graph=50,
        n_attempts_clique_enum=50,
        diversification_fraction=0.1,
        jaccard_opt=0.5,
        jaccard_step=0.1,
    )
    GC_content_score: GCContentScore = GCContentScore(weight=1, GC_content_opt=55)
    Tm_score: TmScore = TmScore(weight=1, Tm_opt=72)
    isoform_consensus_score: IsoformConsensusScore = IsoformConsensusScore(weight=2)


class TargetProbeGlobal(BaseModel):
    model_config = ConfigDict(extra="forbid")

    Tm_parameters: TmParameters = TmParameters(
        nn_table="DNA_NN4",
        tmm_table="DNA_TMM1",
        imm_table="DNA_IMM1",
        de_table="DNA_DE1",
        dnac1=5,
        dnac2=0,
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


class TargetProbes(BaseModel):
    model_config = ConfigDict(extra="forbid")

    oligo_generation: TargetProbeOligoGeneration
    property_filters: TargetProbePropertyFilter
    specificity_filters: TargetProbeSpecificityFilter
    probe_set_selection: TargetProbeProbeSetSelection = Field(description=PROBE_SET_SELECTION_DESC)
    global_parameters: TargetProbeGlobal


############################################
# Readout probes
############################################

# The codebook and readout-probe table support both source="load" and
# source="generate", modelled as discriminated unions on `source`. Both use
# extra="ignore" so the shipped YAML validates in either branch.


class MerfishCodebookLoad(BaseModel):
    model_config = ConfigDict(extra="ignore")

    source: Literal["load"]
    file: str = Field(
        description="Only used when source = load. Path to the codebook file (csv/tsv): columns = 'bits', rows = 'gene_name'; entries are 0/1 bit-encodings for each gene."
    )
    hamming_weight: PositiveInt = Field(
        description="Number of bits set to 1 in each barcode (MERFISH standard = 2). Equals readout overhangs per encoding probe (first half 5' of target, second half 3').",
        default=2,
    )


class MerfishCodebookGenerate(BaseModel):
    model_config = ConfigDict(extra="ignore")

    source: Literal["generate"]
    n_bits: PositiveInt = Field(
        description="Number of bits in each barcode.",
        default=16,
    )
    min_hamming_dist: PositiveInt = Field(
        description="Minimum Hamming distance between any two valid barcodes (MHD4 = 4 for single-bit error correction).",
        default=4,
    )
    hamming_weight: PositiveInt = Field(
        description="Number of bits set to 1 in each barcode (MERFISH standard = 2). Equals readout overhangs per encoding probe (first half 5' of target, second half 3').",
        default=2,
    )


MerfishCodebook = Annotated[MerfishCodebookLoad | MerfishCodebookGenerate, Field(discriminator="source")]


class MerfishReadoutProbeTableLoad(BaseModel):
    model_config = ConfigDict(extra="ignore")

    source: Literal["load"]
    file: str = Field(
        description="Path to the bit-indexed readout probe table (csv/tsv) with columns 'channel', 'readout_probe_id', 'readout_probe_sequence'."
    )


class MerfishReadoutProbeOligoGeneration(BaseModel):
    model_config = ConfigDict(extra="forbid")

    probe_length: PositiveInt = Field(
        description="Length (bases) of each random readout probe.",
        default=20,
    )
    base_probabilities: BaseProbabilities = Field(
        default=BaseProbabilities(A=0.25, C=0.00, G=0.50, T=0.25),
        description="Per-base probability used to generate random readout-probe sequences.",
    )
    initial_num_sequences: PositiveInt = Field(
        description="Number of random sequences to generate before filtering.",
        default=100000,
    )


class MerfishReadoutProbePropertyFilter(BaseModel):
    model_config = ConfigDict(extra="forbid")

    homopolymeric_runs_filter: HomopolymericRunsFilterConfig = Field(
        default=HomopolymericRunsFilterEnabled(
            enabled=True, homopolymeric_base_n=HomopolymericRunThreshold(G=3)
        ),
        description="Exclude probes containing long homopolymeric runs.",
    )
    GC_content_filter: GCContentFilterConfig = Field(
        default=GCContentFilterEnabled(enabled=True, GC_content_min=40, GC_content_max=50),
        description="Enforce GC content within a specified range.",
    )


class MerfishReadoutProbeSpecificityFilter(BaseModel):
    model_config = ConfigDict(extra="forbid")

    specificity_blastn_filter: MerfishReadoutSpecificityBlastnFilterConfig = Field(
        description="Remove readout probes with significant hits to the reference (BLAST-based)."
    )
    cross_hybridization_blastn_filter: CrossHybridizationBlastnFilterConfig = Field(
        default=CrossHybridizationBlastnFilterEnabled(
            enabled=True,
            search_parameters=BlastnSearchParameters(
                perc_identity=100,
                strand="minus",
                word_size=7,
                dust="no",
                soft_masking=False,
                max_target_seqs=10,
            ),
            hit_parameters=BlastnHitParameters(min_alignment_length=11),
        ),
        description="Remove readout probes that cross-hybridize to each other (BLAST-based).",
    )


class MerfishReadoutProbeSetSelection(BaseModel):
    model_config = ConfigDict(extra="forbid")

    set_size: PositiveInt = Field(
        default=16,
        description="Total number of readout probes to select (must equal n_bits so each bit gets a probe).",
    )
    n_combinations: PositiveInt = Field(
        default=100000,
        description="Number of candidate sets to evaluate; higher values may find better sets but increase computation time.",
    )
    homogeneous_properties_weights: HomogeneousPropertiesWeightsT = {"TmNN_oligo": 1, "GC_content_oligo": 1.0}


class MerfishReadoutProbesGlobal(BaseModel):
    model_config = ConfigDict(extra="forbid")

    Tm_parameters: TmParameters = TmParameters(
        nn_table="DNA_NN4",
        tmm_table="DNA_TMM1",
        imm_table="DNA_IMM1",
        de_table="DNA_DE1",
        dnac1=25,
        dnac2=25,
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


class MerfishReadoutProbeTableGenerate(BaseModel):
    model_config = ConfigDict(extra="ignore")

    source: Literal["generate"]
    channels_ids: list[str] = Field(
        description="Fluorescence channels used for round-robin channel assignment across bits.",
        default=["Alexa488", "Cy3b", "Alexa647"],
    )
    oligo_generation: MerfishReadoutProbeOligoGeneration
    property_filters: MerfishReadoutProbePropertyFilter
    specificity_filters: MerfishReadoutProbeSpecificityFilter
    probe_set_selection: MerfishReadoutProbeSetSelection
    global_parameters: MerfishReadoutProbesGlobal


MerfishReadoutProbeTable = Annotated[
    MerfishReadoutProbeTableLoad | MerfishReadoutProbeTableGenerate, Field(discriminator="source")
]


class ReadoutProbes(BaseModel):
    model_config = ConfigDict(extra="forbid")

    codebook: MerfishCodebook
    readout_probe_table: MerfishReadoutProbeTable


############################################
# Primers
############################################

# The forward primer supports source="load" (provide a sequence) and
# source="generate" (design one), modelled as a discriminated union on `source`
# with extra="ignore". The reverse primer only supports source="load"
# (generate_reverse_primer raises FeatureNotImplementedError).


class MerfishForwardPrimerLoad(BaseModel):
    model_config = ConfigDict(extra="ignore")

    source: Literal["load"]
    sequence: DRNAT = Field(
        description="Forward primer sequence placed at the 5' end of the DNA template probe."
    )


class MerfishForwardPrimerOligoGeneration(BaseModel):
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


class MerfishForwardPrimerPropertyFilter(BaseModel):
    model_config = ConfigDict(extra="forbid")

    homopolymeric_runs_filter: HomopolymericRunsFilterConfig = HomopolymericRunsFilterEnabled(
        enabled=True, homopolymeric_base_n=HomopolymericRunThreshold(A=4, T=4, C=4, G=4)
    )
    GC_content_filter: GCContentFilterConfig = Field(
        default=GCContentFilterEnabled(enabled=True, GC_content_min=50, GC_content_max=65),
        description="Enforce GC content within a specified range.",
    )
    GC_clamp_filter: GCClampFilterConfig = Field(
        default=GCClampFilterEnabled(enabled=True, number_GC_GCclamp=1, number_three_prime_base_GCclamp=2),
        description="Require G/C nucleotides at the 3' end of the primer for stable binding.",
    )
    self_complementarity_filter: SelfComplementarityFilterConfig = Field(
        default=SelfComplementarityFilterEnabled(enabled=True, max_len_selfcomplement=6),
        description="Limit self-complementarity to avoid hairpins.",
    )
    complement_reverse_primer_filter: ComplementReversePrimerFilterConfig = Field(
        default=ComplementReversePrimerFilterEnabled(enabled=True, max_len_complement=5),
        description="Limit complementarity between the forward primer and the (known) reverse primer to prevent primer-dimer formation.",
    )
    Tm_filter: TmFilterConfig = TmFilterEnabled(enabled=True, Tm_min=60, Tm_max=75)
    secondary_structure_filter: SecondaryStructureFilterConfig = SecondaryStructureFilterEnabled(
        enabled=True, T=76, thr_DG=0
    )


class MerfishForwardPrimerSpecificityFilter(BaseModel):
    model_config = ConfigDict(extra="forbid")

    specificity_blastn_filter: MerfishPrimerSpecificityBlastnFilterConfig = Field(
        description="Ensure primer specificity against a reference database (BLAST-based)."
    )
    hybridization_probes_blastn_filter: HybridizationProbesBlastnFilterConfig = Field(
        default=HybridizationProbesBlastnFilterEnabled(
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
        ),
        description="Remove primers that match the assembled hybridization probes (BLAST-based).",
    )


class MerfishForwardPrimerGlobal(BaseModel):
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


class MerfishForwardPrimerGenerate(BaseModel):
    model_config = ConfigDict(extra="ignore")

    source: Literal["generate"] = "generate"
    oligo_generation: MerfishForwardPrimerOligoGeneration = MerfishForwardPrimerOligoGeneration()
    property_filters: MerfishForwardPrimerPropertyFilter = MerfishForwardPrimerPropertyFilter()
    specificity_filters: MerfishForwardPrimerSpecificityFilter
    global_parameters: MerfishForwardPrimerGlobal = MerfishForwardPrimerGlobal()


MerfishForwardPrimer = Annotated[
    MerfishForwardPrimerLoad | MerfishForwardPrimerGenerate, Field(discriminator="source")
]


class MerfishReversePrimer(BaseModel):
    model_config = ConfigDict(extra="forbid")

    source: Literal["load"] = "load"
    sequence: DRNAT = Field(
        description="Reverse primer sequence placed at the 3' end of the DNA template probe. The default is the reverse complement of the 20 nt T7 promoter sequence.",
        default="CCCTATAGTGAGTCGTATTA",
    )


class Primers(BaseModel):
    model_config = ConfigDict(extra="forbid")

    forward_primer: MerfishForwardPrimer
    reverse_primer: MerfishReversePrimer = MerfishReversePrimer()


############################################
# Top level
############################################


class MerfishProbeDesignerConfig(BaseModel):
    model_config = ConfigDict(extra="forbid")
    schema_version: Literal[2] = 2
    general: General = General(
        n_jobs=4,
        dir_output="output_Merfishplus_probe_designer",
        write_intermediate_steps=True,
    )

    target_probes: TargetProbes
    readout_probes: ReadoutProbes
    primers: Primers
