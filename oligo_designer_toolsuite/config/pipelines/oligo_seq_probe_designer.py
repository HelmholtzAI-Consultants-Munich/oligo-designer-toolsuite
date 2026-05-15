from pydantic import BaseModel, ConfigDict, Field, PositiveInt

from oligo_designer_toolsuite.config._general_models import (
    BlastnHitParameters,
    BlastnSearchParameters,
    General,
    HomopolymerThresholds,
    TmChemCorrectionParameters,
    TmChemCorrectionParametersDetails,
    TmChemCorrectionParametersEnabled,
    TmParameters,
    TmSaltCorrectionParameters,
    TmSaltCorrectionParametersDisabled,
)
from oligo_designer_toolsuite.config._oligo_scoring import (
    GCContentScore,
    IndependentSetSelection,
    IsoformConsensusScore,
    TargetedExonsScore,
    TmScore,
    UniformDistanceScore,
)
from oligo_designer_toolsuite.config._property_filters import (
    GCContentFilterConfig,
    GCContentFilterEnabled,
    HardMaskedFilterConfig,
    HomopolymericRunsFilterConfig,
    HomopolymericRunsFilterEnabled,
    IsoformConsensusFilterConfig,
    IsoformConsensusFilterEnabled,
    ProhibitedSequencesFilterConfig,
    ProhibitedSequencesFilterEnabled,
    SecondaryStructureFilterConfig,
    SecondaryStructureFilterEnabled,
    SelfComplementarityFilterConfig,
    SelfComplementarityFilterEnabled,
    SoftMaskedFilterConfig,
    TargetedExonsFilterConfig,
    TargetedExonsFilterDisabled,
    TmFilterConfig,
    TmFilterEnabled,
)
from oligo_designer_toolsuite.config._specificity_filters import (
    CrossHybridizationBlastnFilterConfig,
    CrossHybridizationBlastnFilterEnabled,
    ReadLengthBiasFilterConfig,
    ReadLengthBiasFilterEnabled,
    SpecificityBlastnFilterConfig,
    SpecificityBlastnFilterEnabled,
    VariantFilterConfig,
    VariantFilterEnabled,
)
from oligo_designer_toolsuite.config._types import (
    FilesFastaDatabaseT,
    LengthMaxT,
    LengthMinT,
)


class TargetProbeOligoGeneration(BaseModel):
    model_config = ConfigDict(extra="forbid")

    file_region_ids: str | None = Field(
        description="file with a list the genes used to generate the probe sequences, leave empty if all the genes are used",
        default="data/genes/custom_100.txt",
    )
    files_fasta_probe_database: FilesFastaDatabaseT = FilesFastaDatabaseT(
        [
            "data/genomic_regions/exon_annotation_source-NCBI_species-Homo_sapiens_annotation_release-110_genome_assemly-GRCh38.fna",
            "data/genomic_regions/exon_exon_junction_annotation_source-NCBI_species-Homo_sapiens_annotation_release-110_genome_assemly-GRCh38.fna",
        ]
    )
    probe_length_min: LengthMinT = 26
    probe_length_max: LengthMaxT = 30
    probe_split_region: PositiveInt = Field(
        description="Minimum number of bases covering the exon junction, i.e. the oligo should contain at least x bases upstream/downstream of the junction.",
        default=4,
    )


class TargetProbePropertyFilter(BaseModel):
    model_config = ConfigDict(extra="forbid")

    isoform_consensus_filter: IsoformConsensusFilterConfig = IsoformConsensusFilterEnabled(
        enabled=True, isoform_consensus=0
    )
    targeted_exons_filter: TargetedExonsFilterConfig = TargetedExonsFilterDisabled(enabled=False)
    hard_masked_sequences_filter: HardMaskedFilterConfig = HardMaskedFilterConfig(enabled=True)
    soft_masked_sequences_filter: SoftMaskedFilterConfig = SoftMaskedFilterConfig(enabled=True)
    homopolymeric_runs_filter: HomopolymericRunsFilterConfig = HomopolymericRunsFilterEnabled(
        enabled=True, homopolymeric_base_n=HomopolymerThresholds(A=6, T=6, C=6, G=6)
    )
    GC_content_filter: GCContentFilterConfig = GCContentFilterEnabled(
        enabled=True, GC_content_min=45, GC_content_max=65
    )
    prohibited_sequences_filter: ProhibitedSequencesFilterConfig = ProhibitedSequencesFilterEnabled(
        enabled=True,
        prohibited_sequences=["TCT", "CTC"],
        kmer_abundance_threshold={3: 0.0469, 4: 0.0117, 5: 0.0029, 6: 0.00073},
    )
    self_complementarity_filter: SelfComplementarityFilterConfig = SelfComplementarityFilterEnabled(
        enabled=True, max_len_selfcomplement=10
    )
    Tm_filter: TmFilterConfig = TmFilterEnabled(enabled=True, Tm_min=50, Tm_max=70)
    secondary_structure_filter: SecondaryStructureFilterConfig = SecondaryStructureFilterEnabled(
        enabled=True, T=37, thr_DG=0
    )


class TargetProbeSpecificityFilter(BaseModel):
    model_config = ConfigDict(extra="forbid")

    read_length_bias_filter: ReadLengthBiasFilterConfig = ReadLengthBiasFilterEnabled(
        enabled=True, read_length_bias=20
    )
    cross_hybridization_blastn_filter: CrossHybridizationBlastnFilterConfig = (
        CrossHybridizationBlastnFilterEnabled(
            enabled=True,
            search_parameters=BlastnSearchParameters(perc_identity=80, strand="minus", word_size=10),
            hit_parameters=BlastnHitParameters(coverage=50),
        )
    )
    specificity_blastn_filter: SpecificityBlastnFilterConfig = SpecificityBlastnFilterEnabled(
        enabled=True,
        search_parameters=BlastnSearchParameters(perc_identity=80, strand="minus", word_size=10),
        hit_parameters=BlastnHitParameters(coverage=50),
        files_fasta_reference_database=[
            "data/genomic_regions/exon_annotation_source-NCBI_species-Homo_sapiens_annotation_release-110_genome_assemly-GRCh38.fna",
            "data/genomic_regions/exon_exon_junction_annotation_source-NCBI_species-Homo_sapiens_annotation_release-110_genome_assemly-GRCh38.fna",
        ],
    )
    variant_filter: VariantFilterConfig = VariantFilterEnabled(
        enabled=True,
        files_vcf_reference_database=["data/annotations/custom_GCF_000001405.40.chr16.vcf"],
        action="flag",
    )


class TargetProbeProbeSetSelection(BaseModel):
    model_config = ConfigDict(extra="forbid")

    independent_set_selection: IndependentSetSelection = IndependentSetSelection(
        n_sets=3,
        set_size_min=3,
        set_size_opt=5,
        distance_between_target_probes=0,
        n_attempts_graph=50,
        n_attempts_clique_enum=50,
        diversification_fraction=0.1,
        jaccard_opt=0.5,
        jaccard_step=0.1,
    )
    uniform_distance_score: UniformDistanceScore = UniformDistanceScore(weight=1)
    isoform_consensus_score: IsoformConsensusScore = IsoformConsensusScore(weight=1)
    targeted_exons_score: TargetedExonsScore = TargetedExonsScore(weight=0, targeted_exons=["1", "2", "3"])
    GC_content_score: GCContentScore = GCContentScore(
        weight=1, GC_content_min=45, GC_content_opt=55, GC_content_max=65
    )
    Tm_score: TmScore = TmScore(weight=1, Tm_min=50, Tm_opt=60, Tm_max=70)


class TargetProbeGlobal(BaseModel):
    model_config = ConfigDict(extra="forbid")

    Tm_parameters: TmParameters = TmParameters(
        nn_table="DNA_NN3",
        imm_table="DNA_IMM1",
        de_table="DNA_DE1",
        dnac1=50,
        dnac2=0,
        saltcorr=7,
        Na=1000,
        K=0,
        Tris=0,
        Mg=0,
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


class TargetProbe(BaseModel):
    model_config = ConfigDict(extra="forbid")

    oligo_generation: TargetProbeOligoGeneration = TargetProbeOligoGeneration()
    property_filters: TargetProbePropertyFilter = TargetProbePropertyFilter()
    specificity_filters: TargetProbeSpecificityFilter = TargetProbeSpecificityFilter()
    probe_set_selection: TargetProbeProbeSetSelection = TargetProbeProbeSetSelection()
    global_parameters: TargetProbeGlobal = TargetProbeGlobal()


class OligoSeqProbeDesignerConfig(BaseModel):
    model_config = ConfigDict(extra="forbid")
    schema_version: PositiveInt = 2
    general: General = General(
        n_jobs=4,
        dir_output="output_oligo_seq_probe_designer",
        write_intermediate_steps=True,
    )

    target_probe: TargetProbe = TargetProbe()
