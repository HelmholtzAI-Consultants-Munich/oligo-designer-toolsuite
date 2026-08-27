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
    PROPERTY_FILTERS_DESC,
    REQUIRED_PARAMETERS_DESC,
    SPECIFICITY_FILTERS_DESC,
    BlastnHitParameters,
    BlastnHitParametersCoverage,
    BlastnSearchParameters,
    General,
    HomopolymericRunThreshold,
    RequiredParameters,
    TmChemCorrectionParameters,
    TmChemCorrectionParametersDisabled,
    TmParameters,
    TmSaltCorrectionParameters,
    TmSaltCorrectionParametersDisabled,
)
from oligo_designer_toolsuite.config._oligo_scoring import (
    PROBE_SET_SELECTION_DESC,
    IndependentSetSelection,
    IsoformConsensusScore,
    TmScore,
)
from oligo_designer_toolsuite.config._property_filters import (
    GCContentFilterConfig,
    GCContentFilterEnabled,
    HardMaskedFilterConfig,
    HomopolymericRunsFilterConfig,
    HomopolymericRunsFilterEnabled,
    IsoformConsensusFilterConfig,
    IsoformConsensusFilterEnabled,
    SecondaryStructureFilterConfig,
    SecondaryStructureFilterEnabled,
    SoftMaskedFilterConfig,
    TmFilterConfig,
    TmFilterEnabled,
)
from oligo_designer_toolsuite.config._specificity_filters import (
    SPECIFICITY_TARGET_DESC,
    CrossHybridizationBlastnFilterConfig,
    CrossHybridizationBlastnFilterEnabled,
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
        default=13,
    )
    search_parameters: BlastnSearchParameters = Field(
        default=BlastnSearchParameters(
            perc_identity=100,
            strand="minus",
            word_size=10,
            dust="no",
            soft_masking=False,
            max_target_seqs=10,
            max_hsps=1000,
        ),
        description=SpecificityBlastnFilterEnabled.model_fields["search_parameters"].description,
        json_schema_extra={"x-collapsed": True},
    )
    hit_parameters: BlastnHitParameters = BlastnHitParametersCoverage(value=90)


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
        default=45,
        json_schema_extra={"x-quick-setting": True},
    )
    gap_sequence_length: NonNegativeInt = Field(
        description="Length (bases) of the spacer between the L and R arms (covers the junction site).",
        default=2,
        json_schema_extra={"x-quick-setting": True},
    )
    R_probe_sequence_length: PositiveInt = Field(
        description="Length (bases) of the R arm of the probe; L + gap + R equals the total probe length.",
        default=45,
        json_schema_extra={"x-quick-setting": True},
    )


class TargetProbePropertyFilter(BaseModel):
    model_config = ConfigDict(extra="forbid")

    isoform_consensus_filter: IsoformConsensusFilterConfig = IsoformConsensusFilterEnabled(
        enabled=True, isoform_consensus=0
    )
    hard_masked_sequences_filter: HardMaskedFilterConfig = HardMaskedFilterConfig(enabled=True)
    soft_masked_sequences_filter: SoftMaskedFilterConfig = SoftMaskedFilterConfig(enabled=False)
    homopolymeric_runs_filter: HomopolymericRunsFilterConfig = HomopolymericRunsFilterEnabled(
        enabled=True,
        homopolymeric_base_n=HomopolymericRunThreshold(A=6, T=6, C=6, G=6),
    )
    GC_content_filter: GCContentFilterConfig = GCContentFilterEnabled(
        enabled=True, GC_content_min=30, GC_content_max=90
    )
    Tm_filter: TmFilterConfig = TmFilterEnabled(enabled=True, Tm_min=75, Tm_max=100)
    secondary_structure_filter: SecondaryStructureFilterConfig = SecondaryStructureFilterEnabled(
        enabled=True, T=90, thr_DG=0
    )


class TargetProbeSpecificityFilter(BaseModel):
    model_config = ConfigDict(extra="forbid")

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
            hit_parameters=BlastnHitParametersCoverage(value=90),
        )
    )
    specificity_blastn_filter: CycleHcrSpecificityBlastnFilterConfig


class TargetProbeProbeSetSelection(BaseModel):
    model_config = ConfigDict(extra="forbid")

    independent_set_selection: IndependentSetSelection = IndependentSetSelection(
        n_sets=3,
        set_size_min=10,
        set_size_opt=25,
        distance_between_probes=2,
        n_attempts_graph=50,
        n_attempts_clique_enum=50,
        diversification_fraction=0.1,
        jaccard_opt=0.5,
        jaccard_step=0.1,
    )
    isoform_consensus_score: IsoformConsensusScore = IsoformConsensusScore(weight=10)
    Tm_score: TmScore = TmScore(weight=1, Tm_opt=75)


class TargetProbeGlobal(BaseModel):
    model_config = ConfigDict(extra="forbid")

    Tm_parameters: TmParameters = TmParameters(
        nn_table="DNA_NN3",
        tmm_table="DNA_TMM1",
        imm_table="DNA_IMM1",
        de_table="DNA_DE1",
        dnac1=25,
        dnac2=25,
        saltcorr=0,
        Na=50,
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
        default=4,
    )


CycleHcrCodebook = Annotated[
    CycleHcrCodebookLoad | CycleHcrCodebookGenerate,
    Field(discriminator="source"),
]


class CycleHcrReadoutProbeTable(BaseModel):
    model_config = ConfigDict(extra="forbid")

    source: Literal["load"] = "load"
    file: str = Field(
        description="Path to the readout probe table (csv/tsv) with columns 'channel', 'readout_probe_id', 'readout_probe_sequence', and 'L/R'. Bit handling depends on codebook.source: when codebook.source = 'load' the file MUST also contain a 'bit' column whose values match the codebook columns (the user is responsible for that mapping); when codebook.source = 'generate' any 'bit' column is ignored and bits are reassigned deterministically by sorting on (readout_probe_id, channel, L/R)."
    )


class ReadoutProbes(BaseModel):
    model_config = ConfigDict(extra="forbid")

    codebook: CycleHcrCodebook = CycleHcrCodebookGenerate()
    readout_probe_table: CycleHcrReadoutProbeTable


############################################
# Primers
############################################

# Primer generation is not implemented (generate_forward_primer/generate_reverse_primer raise
# FeatureNotImplementedError), so `source` is modelled as Literal["load"] to reject "generate".


class CycleHcrPrimer(BaseModel):
    model_config = ConfigDict(extra="forbid")

    source: Literal["load"] = "load"
    sequence: DRNAT


class Primers(BaseModel):
    model_config = ConfigDict(extra="forbid")

    forward_primer: CycleHcrPrimer = Field(
        default=CycleHcrPrimer(sequence="TAATACGACTCACTATAGCGTCATC"),
        description="Forward PCR primer placed at the 5' end of the DNA template probe. The default is the T7 promoter sequence.",
    )
    reverse_primer: CycleHcrPrimer = Field(
        default=CycleHcrPrimer(sequence="CGACACCGAACGTGCGACAA"),
        description="Reverse PCR primer placed at the 3' end of the DNA template probe.",
    )


############################################
# Hybridization probes
############################################


class HybridizationProbes(BaseModel):
    model_config = ConfigDict(extra="forbid")

    linker_sequence: DRNAT = Field(
        description="Linker sequence between the target-binding L/R arm and the readout-probe barcode (used by both hybridization and DNA template assembly).",
        default="TT",
    )


############################################
# Top level
############################################


# The front end builds its form from this, so `general` stays out of it.
class CycleHcrProbeDesignerConfigBase(BaseModel):
    model_config = ConfigDict(extra="forbid")
    schema_version: Literal[2] = 2
    target_probes: TargetProbes
    readout_probes: ReadoutProbes
    primers: Primers
    hybridization_probes: HybridizationProbes


class CycleHcrProbeDesignerConfig(CycleHcrProbeDesignerConfigBase):
    general: General = General(
        n_jobs=4,
        dir_output="output_cyclehcr_probe_designer",
        write_intermediate_steps=True,
    )

    required_parameters: RequiredParameters = Field(description=REQUIRED_PARAMETERS_DESC)
