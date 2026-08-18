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
    BlastnHitParameters,
    BlastnHitParametersCoverage,
    BlastnSearchParameters,
    General,
    HomopolymericRunThreshold,
    RequiredParameters,
    TmChemCorrectionParameters,
    TmChemCorrectionParametersDetails,
    TmChemCorrectionParametersEnabled,
    TmParameters,
    TmSaltCorrectionParameters,
    TmSaltCorrectionParametersDisabled,
)
from oligo_designer_toolsuite.config._oligo_scoring import (
    GCContentScoreNormalized,
    IndependentSetSelection,
    IsoformConsensusScore,
    TmScoreNormalized,
)
from oligo_designer_toolsuite.config._property_filters import (
    GCContentFilterConfig,
    GCContentFilterEnabled,
    HardMaskedFilterConfig,
    HomopolymericRunsFilterConfig,
    HomopolymericRunsFilterEnabled,
    IsoformConsensusFilterConfig,
    IsoformConsensusFilterEnabled,
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
from oligo_designer_toolsuite.config._types import (
    LengthMaxT,
    LengthMinT,
    TmMaxT,
    TmMinT,
    TmOptT,
)

############################################
# SCRINSHOT-specific overrides
############################################


class ScrinshotSpecificityBlastnFilterEnabled(SpecificityBlastnFilterEnabled):
    search_parameters: BlastnSearchParameters = BlastnSearchParameters(
        perc_identity=80,
        strand="minus",
        word_size=10,
        dust="no",
        soft_masking=False,
        max_target_seqs=10,
        max_hsps=1000,
    )
    hit_parameters: BlastnHitParameters = BlastnHitParametersCoverage(value=50)
    ligation_region_size: NonNegativeInt = Field(
        description=(
            "Size of the seed region around the ligation site for BLASTN seed-region filtering. "
            "If > 0, seed-based filtering is applied around the ligation site, removing all probes"
            "where BLASTN hits cover the junction region regardless of the "
            "coverage threshold. If 0, full-length specificity filtering is performed instead."
        ),
        default=5,
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

    probe_length_min: LengthMinT = 40
    probe_length_max: LengthMaxT = 45

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
        default=10,
    )
    Tm_dif_max: NonNegativeFloat = Field(
        description="Maximum Tm difference (°C) between the two padlock arms.",
        default=2,
    )
    Tm_min: TmMinT = 50
    Tm_max: TmMaxT = 60

    @model_validator(mode="after")
    def _check_min_max(self) -> Self:
        if self.Tm_min > self.Tm_max:
            raise ValueError(
                f"'padlock_arm_Tm_min' ({self.Tm_min}) must be <= 'padlock_arm_Tm_max' ({self.Tm_max})"
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
        enabled=True, homopolymeric_base_n=HomopolymericRunThreshold(A=5, T=5, C=5, G=5)
    )
    GC_content_filter: GCContentFilterConfig = GCContentFilterEnabled(
        enabled=True, GC_content_min=40, GC_content_max=60
    )
    Tm_filter: TmFilterConfig = TmFilterEnabled(enabled=True, Tm_min=65, Tm_max=75)


class TargetProbeSpecificityFilter(BaseModel):
    model_config = ConfigDict(extra="forbid")

    cross_hybridization_blastn_filter: CrossHybridizationBlastnFilterConfig = (
        CrossHybridizationBlastnFilterEnabled(
            enabled=True,
            search_parameters=BlastnSearchParameters(
                perc_identity=80,
                strand="minus",
                word_size=10,
                dust="no",
                soft_masking=False,
                max_target_seqs=10,
            ),
            hit_parameters=BlastnHitParametersCoverage(value=80),
        )
    )
    specificity_blastn_filter: ScrinshotSpecificityBlastnFilterConfig


class TargetProbeProbeSetSelection(BaseModel):
    model_config = ConfigDict(extra="forbid")

    independent_set_selection: IndependentSetSelection = IndependentSetSelection(
        n_sets=3,
        set_size_min=3,
        set_size_opt=5,
        distance_between_probes=0,
        n_attempts_graph=20,
        n_attempts_clique_enum=70,
        diversification_fraction=0.1,
        jaccard_opt=0.5,
        jaccard_step=0.1,
    )
    isoform_consensus_score: IsoformConsensusScore = IsoformConsensusScore(weight=2)
    GC_content_score: GCContentScoreNormalized = GCContentScoreNormalized(
        weight=1, GC_content_min=40, GC_content_opt=50, GC_content_max=60
    )
    Tm_score: TmScoreNormalized = TmScoreNormalized(weight=1, Tm_min=65, Tm_opt=70, Tm_max=75)


class TargetProbeGlobal(BaseModel):
    model_config = ConfigDict(extra="forbid")

    Tm_parameters: TmParameters = TmParameters(
        nn_table="DNA_NN3",
        tmm_table="DNA_TMM1",
        imm_table="DNA_IMM1",
        de_table="DNA_DE1",
        dnac1=50,
        dnac2=0,
        saltcorr=7,
        Na=39,
        K=75,
        Tris=20,
        Mg=10,
        dNTPs=0,
    )
    Tm_chem_correction_parameters: TmChemCorrectionParameters = TmChemCorrectionParametersEnabled(
        enabled=True,
        parameters=TmChemCorrectionParametersDetails(
            DMSO=0, DMSOfactor=0.75, fmd=20, fmdfactor=0.65, fmdmethod=1, GC=None
        ),
    )
    Tm_salt_correction_parameters: TmSaltCorrectionParameters = TmSaltCorrectionParametersDisabled(
        enabled=False
    )


class TargetProbes(BaseModel):
    model_config = ConfigDict(extra="forbid")

    oligo_generation: TargetProbeOligoGeneration
    padlock_arms_properties: PadlockArmsProperties
    property_filters: TargetProbePropertyFilter
    specificity_filters: TargetProbeSpecificityFilter
    probe_set_selection: TargetProbeProbeSetSelection
    global_parameters: TargetProbeGlobal


############################################
# Detection oligo
############################################


class DetectionOligoOligoGeneration(BaseModel):
    model_config = ConfigDict(extra="forbid")

    min_thymines: PositiveInt = Field(
        description="Minimal number of thymines (T) in the detection oligo (required for UNG cleavage after U-substitution).",
        default=2,
    )
    oligo_length_min: LengthMinT = 15
    oligo_length_max: LengthMaxT = 40
    U_distance: PositiveInt = Field(
        description="Preferred minimum distance (bases) between consecutive uracils.",
        default=5,
    )
    Tm_opt: TmOptT = 56

    @model_validator(mode="after")
    def _check_min_max(self) -> Self:
        if self.oligo_length_min > self.oligo_length_max:
            raise ValueError(
                f"'oligo_length_min' ({self.oligo_length_min}) must be <= 'oligo_length_max' ({self.oligo_length_max})"
            )
        return self


class DetectionOligoGlobal(BaseModel):
    model_config = ConfigDict(extra="forbid")

    Tm_parameters: TmParameters = TmParameters(
        nn_table="DNA_NN3",
        tmm_table="DNA_TMM1",
        imm_table="DNA_IMM1",
        de_table="DNA_DE1",
        dnac1=50,
        dnac2=0,
        saltcorr=7,
        Na=39,
        K=0,
        Tris=0,
        Mg=0,
        dNTPs=0,
    )
    Tm_chem_correction_parameters: TmChemCorrectionParameters = TmChemCorrectionParametersEnabled(
        enabled=True,
        parameters=TmChemCorrectionParametersDetails(
            DMSO=0, DMSOfactor=0.75, fmd=30, fmdfactor=0.65, fmdmethod=1, GC=None
        ),
    )
    Tm_salt_correction_parameters: TmSaltCorrectionParameters = TmSaltCorrectionParametersDisabled(
        enabled=False
    )


class DetectionOligo(BaseModel):
    model_config = ConfigDict(extra="forbid")

    oligo_generation: DetectionOligoOligoGeneration
    global_parameters: DetectionOligoGlobal


############################################
# Top level
############################################


class ScrinshotProbeDesignerConfig(BaseModel):
    model_config = ConfigDict(extra="forbid")
    schema_version: Literal[2] = 2
    general: General = General(
        n_jobs=4,
        dir_output="output_scrinshot_probe_designer",
        write_intermediate_steps=True,
    )

    required_parameters: RequiredParameters
    target_probes: TargetProbes
    detection_oligo: DetectionOligo
