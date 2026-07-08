from typing import Annotated, Literal

from pydantic import (
    BaseModel,
    ConfigDict,
    Field,
    NonNegativeInt,
    PositiveInt,
)

from oligo_designer_toolsuite.config._general_models import (
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
    IndependentSetSelection,
    IsoformConsensusScore,
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
    CrossHybridizationBlastnFilterConfig,
    CrossHybridizationBlastnFilterEnabled,
    SpecificityBlastnFilterDisabled,
    SpecificityBlastnFilterEnabled,
)
from oligo_designer_toolsuite.config._types import (
    DRNAT,
    FilesFastaDatabaseT,
    RegionListT,
)

############################################
# HCR-specific overrides
############################################


class HcrSpecificityBlastnFilterDisabled(SpecificityBlastnFilterDisabled):
    pass


class HcrSpecificityBlastnFilterEnabled(SpecificityBlastnFilterEnabled):
    junction_region_size: NonNegativeInt = Field(
        description=(
            "Size of the seed region around the junction site for BLASTN seed-region filtering. "
            "If > 0, seed-based filtering is applied around the junction site, removing all probes "
            "where BLASTN hits cover the junction region regardless of the coverage threshold. "
            "If 0, full-length specificity filtering is performed instead."
        ),
        default=0,
    )
    search_parameters: BlastnSearchParameters = BlastnSearchParameters(
        perc_identity=100,
        strand="minus",
        word_size=10,
        dust="no",
        soft_masking=False,
        max_target_seqs=10,
        max_hsps=1000,
    )
    hit_parameters: BlastnHitParameters = BlastnHitParameters(coverage=90)


HcrSpecificityBlastnFilterConfig = Annotated[
    HcrSpecificityBlastnFilterEnabled | HcrSpecificityBlastnFilterDisabled,
    Field(discriminator="enabled"),
]


############################################
# Target probe
############################################


class TargetProbeOligoGeneration(BaseModel):
    model_config = ConfigDict(extra="forbid")

    file_region_ids: RegionListT
    files_fasta_probe_database: FilesFastaDatabaseT
    L_probe_sequence_length: PositiveInt = Field(
        description="Length (bases) of the L arm of the probe; L + gap + R equals the total probe length.",
        default=25,
    )
    gap_sequence_length: NonNegativeInt = Field(
        description="Length (bases) of the spacer between the L and R arms (covers the ligation site).",
        default=2,
    )
    R_probe_sequence_length: PositiveInt = Field(
        description="Length (bases) of the R arm of the probe; L + gap + R equals the total probe length.",
        default=25,
    )


class TargetProbePropertyFilter(BaseModel):
    model_config = ConfigDict(extra="forbid")

    isoform_consensus_filter: IsoformConsensusFilterConfig = IsoformConsensusFilterEnabled(
        enabled=True, isoform_consensus=0
    )
    hard_masked_sequences_filter: HardMaskedFilterConfig = HardMaskedFilterConfig(enabled=True)
    soft_masked_sequences_filter: SoftMaskedFilterConfig = SoftMaskedFilterConfig(enabled=False)
    homopolymeric_runs_filter: HomopolymericRunsFilterConfig = HomopolymericRunsFilterEnabled(
        enabled=True, homopolymeric_base_n=HomopolymericRunThreshold(A=4, T=4, C=4, G=4)
    )
    GC_content_filter: GCContentFilterConfig = GCContentFilterEnabled(
        enabled=True, GC_content_min=40, GC_content_max=65
    )
    Tm_filter: TmFilterConfig = TmFilterEnabled(enabled=True, Tm_min=50, Tm_max=80)
    secondary_structure_filter: SecondaryStructureFilterConfig = SecondaryStructureFilterEnabled(
        enabled=True, T=37, thr_DG=-7
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
            hit_parameters=BlastnHitParameters(coverage=90),
        )
    )
    specificity_blastn_filter: HcrSpecificityBlastnFilterConfig


class TargetProbeProbeSetSelection(BaseModel):
    model_config = ConfigDict(extra="forbid")

    independent_set_selection: IndependentSetSelection = IndependentSetSelection(
        n_sets=30,
        set_size_min=10,
        set_size_opt=25,
        distance_between_probes=2,
        n_attempts_graph=50,
        n_attempts_clique_enum=50,
        diversification_fraction=0.1,
        jaccard_opt=0.5,
        jaccard_step=0.1,
    )
    isoform_consensus_score: IsoformConsensusScore = IsoformConsensusScore(weight=1)


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


class TargetProbe(BaseModel):
    model_config = ConfigDict(extra="forbid")

    oligo_generation: TargetProbeOligoGeneration
    property_filters: TargetProbePropertyFilter
    specificity_filters: TargetProbeSpecificityFilter
    probe_set_selection: TargetProbeProbeSetSelection
    global_parameters: TargetProbeGlobal


############################################
# Initiator probes
############################################

# HCR only supports source="load" at the moment: InitiatorDesigner.generate_codebook and
# generate_initiator_table both raise FeatureNotImplementedError. Modelling `source` as
# Literal["load"] rejects source="generate" at validation. When generation is implemented,
# introduce a discriminated union on `source`.


class HcrCodebook(BaseModel):
    model_config = ConfigDict(extra="forbid")

    source: Literal["load"] = "load"
    file: str = Field(
        description="Path to the codebook file: columns = 'bits', rows = 'gene_name'; entries are 0/1 bit-encodings for each gene."
    )


class HcrInitiatorTable(BaseModel):
    model_config = ConfigDict(extra="forbid")

    source: Literal["load"] = "load"
    file: str = Field(
        description="Path to the bit-indexed initiator table (csv/tsv) with columns 'bit', 'initiator_L_sequence', and 'initiator_R_sequence'."
    )


class InitiatorProbes(BaseModel):
    model_config = ConfigDict(extra="forbid")

    codebook: HcrCodebook
    initiator_table: HcrInitiatorTable


############################################
# Hybridization probes
############################################


class HybridizationProbes(BaseModel):
    model_config = ConfigDict(extra="forbid")

    linker_sequence: DRNAT = Field(
        description="Linker sequence between the initiator and the target-binding L/R arm.",
        default="AA",
    )


############################################
# Top level
############################################


class HcrProbeDesignerConfig(BaseModel):
    model_config = ConfigDict(extra="forbid")
    schema_version: Literal[2] = 2
    general: General = General(
        n_jobs=4,
        dir_output="output_hcr_probe_designer",
        write_intermediate_steps=True,
    )

    target_probe: TargetProbe
    initiator_probes: InitiatorProbes
    hybridization_probes: HybridizationProbes
