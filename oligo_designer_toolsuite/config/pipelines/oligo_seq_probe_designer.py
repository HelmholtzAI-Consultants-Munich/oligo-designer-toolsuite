from pydantic import BaseModel, ConfigDict, Field, PositiveInt

from oligo_designer_toolsuite.config._general_models import (
    BlastnHitParameters,
    BlastnSearchParameters,
    General,
    HomopolymerThresholds,
    TmChemCorrectionParameters,
    TmChemCorrectionParametersCustom,
    TmChemCorrectionParametersDetails,
    TmParameters,
    TmParametersCustom,
    TmParametersDetails,
    TmSaltCorrectionParameters,
    TmSaltCorrectionParametersDisabled,
)
from oligo_designer_toolsuite.config._property_filters import (
    GCContentFilterConfig,
    GCContentFilterEnabled,
    HardMaskedFilterConfig,
    HomopolymericRunsFilterConfig,
    HomopolymericRunsFilterEnabled,
    IsoformConsensusFilterConfig,
    IsoformConsensusFilterEnabled,
    MeltingTemperatureFilterConfig,
    MeltingTemperatureFilterEnabled,
    ProhibitedSequencesFilterConfig,
    ProhibitedSequencesFilterEnabled,
    SecondaryStructureFilterConfig,
    SecondaryStructureFilterEnabled,
    SelfComplementarityFilterConfig,
    SelfComplementarityFilterEnabled,
    SoftMaskedFilterConfig,
)
from oligo_designer_toolsuite.config._specificity_filters import (
    CrossHybridizationFilterConfig,
    CrossHybridizationFilterEnabled,
    ReadLengthBiasFilterConfig,
    ReadLengthBiasFilterEnabled,
    SpecificityBlastnFilterConfig,
    SpecificityBlastnFilterEnabled,
    VariantFilterConfig,
    VariantFilterEnabled,
)
from oligo_designer_toolsuite.config._types import (
    FilesFastaDatabaseT,
    GCContentOptT,
    LengthMaxT,
    LengthMinT,
    TmOptT,
)


class TargetProbeOligoGeneration(BaseModel):
    model_config = ConfigDict(extra="forbid")

    file_regions: str | None = Field(
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
    melting_temperature_filter: MeltingTemperatureFilterConfig = MeltingTemperatureFilterEnabled(
        enabled=True, Tm_min=50, Tm_max=70
    )
    secondary_structure_filter: SecondaryStructureFilterConfig = SecondaryStructureFilterEnabled(
        enabled=True, T=37, thr_DG=0
    )


class TargetProbeSpecificityFilter(BaseModel):
    model_config = ConfigDict(extra="forbid")

    read_length_bias_filter: ReadLengthBiasFilterConfig = ReadLengthBiasFilterEnabled(
        enabled=True, read_length_bias=20
    )
    cross_hybridization_filter: CrossHybridizationFilterConfig = CrossHybridizationFilterEnabled(
        enabled=True,
        search_parameters=BlastnSearchParameters(perc_identity=80, strand="minus", word_size=10),
        hit_parameters=BlastnHitParameters(coverage=50),
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
        enabled=True, files_vcf_reference_database=["data/annotations/custom_GCF_000001405.40.chr16.vcf"]
    )


class TargetProbeProbeSetSelection(BaseModel):
    model_config = ConfigDict(extra="forbid")

    n_sets: PositiveInt = Field(
        description="Number of oligo sets to generate per region. Multiple sets allow for redundancy and selection of the best-performing set based on scoring criteria.",
        default=3,
    )
    set_size_min: PositiveInt = Field(
        description="Minimum size (number of probes) required for each oligo set. Sets with fewer probes than this value will be rejected, and regions that cannot generate sets meeting this minimum will be removed. Regions with less oligos will be filtered out and stored in regions_with_insufficient_oligos_for_db_probes.",
        default=3,
    )
    set_size_opt: PositiveInt = Field(
        description="Optimal size (number of probes) for each oligo set. The set selection algorithm will attempt to generate sets of this size, but may produce sets with fewer probes if constraints cannot be met.",
        default=5,
    )
    distance_between_target_probes: int = Field(
        description="How much overlap should be allowed between oligos, e.g. if oligos can overlap x bases choose -x, if oligos can be next to one another choose 0, if oligos should be x bases apart choose x.",
        default=0,
    )
    uniform_distance_weight: float = Field(
        description="Weight assigned to uniform distance in the scoring function.", default=1
    )
    isoform_weight: float = Field(
        description="Weight assigned to isoform consensus in the scoring function.", default=1
    )
    targeted_exons_weight: float = Field(
        description="Weight assigned to targeted exons overlap in the scoring function.", default=1
    )
    targeted_exons: list[str] = Field(
        description="List of exon identifiers that should be preferentially targeted by probes. Probes overlapping these exons receive higher scores.",
        default=["1", "2", "3"],
    )
    GC_weight: float = Field(description="Weight assigned to GC content in the scoring function.", default=1)
    GC_content_opt: GCContentOptT = 55
    Tm_weight: float = Field(
        description="Weight assigned to melting temperature in the scoring function.", default=1
    )
    Tm_opt: TmOptT = 60
    n_attempts_graph: PositiveInt = Field(
        description="Number of randomized graph attempts. In each attempt, a fraction of nodes is randomly removed from the compatibility graph to create diversity; more attempts increase diversity at the cost of runtime.",
        default=50,
    )
    n_attempts_clique_enum: PositiveInt = Field(
        description="Maximum number of cliques enumerated per graph attempt. Limits how many cliques are explored before stopping enumeration for the current graph.",
        default=50,
    )
    diversification_fraction: float = Field(
        ge=0,
        le=1,
        description="Fraction of oligos to remove from the graph per attempt to create diversity in the set selection.",
        default=0.1,
    )
    jaccard_opt: float = Field(
        ge=0,
        le=1,
        description="Optimal maximum Jaccard overlap allowed between selected sets. Lower values enforce more diversity between sets.",
        default=0.5,
    )
    jaccard_step: float = Field(
        ge=0,
        le=1,
        description="Step size used to relax the Jaccard constraint when not enough sets are found.",
        default=0.1,
    )


class TargetProbeAdvanced(BaseModel):
    model_config = ConfigDict(extra="forbid")

    Tm_parameters: TmParameters = TmParametersCustom(
        mode="custom",
        parameters=TmParametersDetails(
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
        ),
    )
    Tm_chem_correction_parameters: TmChemCorrectionParameters = TmChemCorrectionParametersCustom(
        mode="custom",
        parameters=TmChemCorrectionParametersDetails(
            DMSO=0, DMSOfactor=0.75, fmd=20, fmdfactor=0.65, fmdmethod=1, GC=None
        ),
    )
    Tm_salt_correction_parameters: TmSaltCorrectionParameters = TmSaltCorrectionParametersDisabled(
        mode="disabled"
    )


class TargetProbe(BaseModel):
    model_config = ConfigDict(extra="forbid")

    oligo_generation: TargetProbeOligoGeneration = TargetProbeOligoGeneration()
    property_filters: TargetProbePropertyFilter = TargetProbePropertyFilter()
    specificity_filters: TargetProbeSpecificityFilter = TargetProbeSpecificityFilter()
    probe_set_selection: TargetProbeProbeSetSelection = TargetProbeProbeSetSelection()
    advanced: TargetProbeAdvanced = TargetProbeAdvanced()


class OligoSeqProbeDesignerConfig(BaseModel):
    model_config = ConfigDict(extra="forbid")
    schema_version: PositiveInt = 2
    general: General = General(
        n_jobs=4,
        dir_output="output_oligo_seq_probe_designer",
        write_intermediate_steps=True,
    )

    target_probe: TargetProbe = TargetProbe()
