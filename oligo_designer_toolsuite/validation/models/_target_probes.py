from typing import Annotated

from pydantic import BaseModel, ConfigDict, Field, NonNegativeFloat, NonNegativeInt, PositiveInt

from oligo_designer_toolsuite.validation._types import (
    DNAT,
    FilesFastaDatabaseT,
    FilesFastaReferenceDatabaseT,
    FractionT,
    GCContentMaxT,
    GCContentMinT,
    GCContentOptT,
    LengthMaxT,
    LengthMinT,
    SecondaryStructuresThresholdDeltaGT,
    TmMaxT,
    TmMinT,
    TmOptT,
    TSecondaryStructureT,
    WeightT,
)
from oligo_designer_toolsuite.validation.models._general import HomopolymerThresholds


class TargetProbeBase(BaseModel):
    model_config = ConfigDict(extra="forbid")

    file_regions: Annotated[
        str | None,
        Field(
            description="file with a list the genes used to generate the probe sequences, leave empty if all the genes are used",
        ),
    ]
    files_fasta_database: FilesFastaDatabaseT
    files_fasta_reference_database: FilesFastaReferenceDatabaseT

    isoform_consensus: Annotated[
        float,
        Field(
            description="Threshold for isoform consensus filtering (typically between 0.0 and 1.0). Probes with isoform consensus values below this threshold will be filtered out. This ensures that selected probes target sequences that are conserved across multiple transcript isoforms.",
            ge=0,
            le=100,
        ),
    ]
    GC_content_min: GCContentMinT
    GC_content_max: GCContentMaxT
    homopolymeric_base_n: Annotated[
        HomopolymerThresholds,
        Field(
            description="Specifying the maximum allowed length of homopolymeric runs for each nucleotide base. Keys should be 'A', 'T', 'G', 'C' and values are the maximum run length. For example: A: 3, T: 3, G: 3, C: 3 allows up to 3 consecutive identical bases."
        ),
    ]

    set_size_min: Annotated[
        PositiveInt,
        Field(
            description="Minimum size (number of probes) required for each oligo set. Sets with fewer probes than this value will be rejected, and regions that cannot generate sets meeting this minimum will be removed. Regions with less oligos will be filtered out and stored in regions_with_insufficient_oligos_for_db_probes."
        ),
    ]
    set_size_opt: Annotated[
        PositiveInt,
        Field(
            description="Optimal size (number of probes) for each oligo set. The set selection algorithm will attempt to generate sets of this size, but may produce sets with fewer probes if constraints cannot be met."
        ),
    ]
    distance_between_target_probes: Annotated[
        int,
        Field(
            description="How much overlap should be allowed between oligos, e.g. if oligos can overlap x bases choose -x, if oligos can be next to one another choose 0, if oligos should be x bases apart choose x."
        ),
    ]
    n_sets: Annotated[
        PositiveInt,
        Field(
            description="Number of oligo sets to generate per region. Multiple sets allow for redundancy and selection of the best-performing set based on scoring criteria."
        ),
    ]


class TargetProbeCycleHCR(TargetProbeBase):

    L_probe_sequence_length: Annotated[
        PositiveInt,
        Field(
            description="Length of the left probe sequence in nucleotides. This is the 5' portion of the target probe that binds to the 3' end of the target RNA sequence. L + spacer + R sequence should equal the total probe length, e.g. 45 + 2 + 45 = 92",
        ),
    ]
    gap_sequence_length: Annotated[
        NonNegativeInt,
        Field(
            description="Length of the gap sequence between left and right probes in nucleotides. This gap is not included in the probe sequences but represents the spacing between the two probe halves on the target transcript.",
        ),
    ]
    R_probe_sequence_length: Annotated[
        PositiveInt,
        Field(
            description="Length of the right probe sequence in nucleotides. This is the 3' portion of the target probe that binds to the 5' end of the target RNA sequence. L + spacer + R sequence should equal the total probe length, e.g. 45 + 2 + 45 = 92",
        ),
    ]
    Tm_min: TmMinT
    Tm_max: TmMaxT
    T_secondary_structure: TSecondaryStructureT
    junction_region_size: Annotated[
        NonNegativeInt,
        Field(
            description="Size of the junction region (in nucleotides) used for seed-based specificity filtering. If set to 0, full-length specificity filtering is used instead of seed-based filtering. When seed-based filtering is enabled, any probe with a BLASTN hit covering the junction region between the left and right probe halves will be removed, regardless of the alignment coverage percentage.",
        ),
    ]
    Tm_weight: Annotated[
        WeightT,
        Field(
            description="Weight assigned to melting temperature (Tm) in the probe scoring function. Higher values prioritize probes with Tm closer to the optimal value (Tm_max). This weight is used in combination with isoform_weight to calculate a composite score for each probe.",
        ),
    ]
    isoform_weight: Annotated[
        WeightT,
        Field(
            description="Weight assigned to isoform consensus in the probe scoring function. Higher values prioritize probes with higher isoform consensus values (probes that are conserved across multiple transcript isoforms). This weight is used in combination with Tm_weight to calculate a composite score for each probe.",
        ),
    ]
    linker_sequence: Annotated[
        DNAT,
        Field(
            description="DNA sequence used to link target probes and readout probes in the hybridization probe. This sequence is inserted between the target probe sequence and the readout probe sequence during assembly. Typically a short spacer sequence (e.g., 'TT').",
        ),
    ]


class TargetProbeMerfish(TargetProbeBase):

    GC_content_opt: GCContentOptT

    length_min: LengthMinT
    length_max: LengthMaxT

    Tm_min: TmMinT
    Tm_opt: TmOptT
    Tm_max: TmMaxT

    T_secondary_structure: TSecondaryStructureT

    GC_weight: Annotated[
        WeightT,
        Field(
            description="Weight assigned to GC content in the scoring function.",
        ),
    ]
    Tm_weight: Annotated[
        WeightT,
        Field(
            description="Weight assigned to melting temperature (Tm) in the probe scoring function. Higher values prioritize probes with Tm closer to the optimal value (Tm_max). This weight is used in combination with isoform_weight to calculate a composite score for each probe.",
        ),
    ]
    isoform_weight: WeightT


class TargetProbeSeqFishPlus(TargetProbeBase):

    GC_content_opt: GCContentOptT

    length_min: LengthMinT
    length_max: LengthMaxT

    T_secondary_structure: TSecondaryStructureT

    GC_weight: Annotated[
        WeightT,
        Field(
            description="Weight assigned to GC content in the scoring function.",
        ),
    ]
    UTR_weight: Annotated[
        WeightT,
        Field(
            description="Weight assigned to untranslated region (UTR) targeting in the scoring function. Probes that target UTRs are prioritized.",
        ),
    ]


class TargetProbeOligoSeq(TargetProbeBase):

    files_vcf_reference_database: Annotated[
        list[str],
        Field(
            description="List of paths to VCF files containing variant information used for filtering probes that overlap with known single nucleotide polymorphisms (SNPs) or other variants. Probes overlapping variants may have reduced specificity.",
        ),
    ]
    length_min: LengthMinT
    length_max: LengthMaxT
    split_region: Annotated[
        PositiveInt,
        Field(
            description="Size of regions (in nucleotides) to split large genomic regions into when generating probe candidates. This helps manage memory usage for very long sequences."
        ),
    ]
    Tm_min: TmMinT
    Tm_max: TmMaxT
    T_secondary_structure: TSecondaryStructureT
    secondary_structures_threshold_deltaG: SecondaryStructuresThresholdDeltaGT
    kmer_abundance_threshold: Annotated[
        dict[PositiveInt, FractionT],
        Field(
            description="The maximum abundance of a k-mer allowed. If a k-mer has an abundance greater than this threshold, all oligos containing this k-mer will be filtered out. This is a dictionary with the k-mer length as the key and the threshold as the value."
        ),
    ]
    prohibited_sequences: Annotated[
        list[DNAT],
        Field(
            description="The sequences to prohibit in the oligos. If an oligo contains any of these sequences, it will be filtered out."
        ),
    ]
    max_len_selfcomplement: Annotated[
        NonNegativeInt,
        Field(
            description="Maximum allowable length of self-complementary sequences. Probes with longer self-complementary regions can form hairpins and reduce hybridization efficiency."
        ),
    ]
    read_length_bias: Annotated[
        NonNegativeInt,
        Field(
            description="Number of nucleotides from the 5' end of probes to check for read length bias. Probes where the first N bases match exactly with other probes are removed to prevent sequencing read length biases."
        ),
    ]
    uniform_distance_weight: Annotated[
        WeightT, Field(description="Weight assigned to uniform distance in the scoring function.")
    ]
    isoform_weight: Annotated[
        WeightT, Field(description="Weight assigned to isoform consensus in the scoring function.")
    ]
    targeted_exons_weight: Annotated[
        WeightT, Field(description="Weight assigned to targeted exons overlap in the scoring function.")
    ]
    targeted_exons: Annotated[
        list[str],
        Field(
            description="List of exon identifiers that should be preferentially targeted by probes. Probes overlapping these exons receive higher scores."
        ),
    ]
    GC_weight: Annotated[WeightT, Field(description="Weight assigned to GC content in the scoring function.")]
    GC_content_opt: GCContentOptT
    Tm_weight: Annotated[
        WeightT, Field(description="Weight assigned to melting temperature in the scoring function.")
    ]
    Tm_opt: TmOptT


class TargetProbeScrinshot(TargetProbeBase):

    length_min: LengthMinT
    length_max: LengthMaxT

    GC_content_opt: GCContentOptT
    Tm_min: TmMinT
    Tm_opt: TmOptT
    Tm_max: TmMaxT

    padlock_arm_Tm_dif_max: Annotated[
        NonNegativeInt,
        Field(
            description="Maximum allowed difference in melting temperature (in degrees Celsius) between the two padlock arms. This ensures balanced binding of both arms. The difference shouldn't be higher than 5. But the range is not super important, the lower the better"
        ),
    ]
    padlock_arm_length_min: Annotated[
        PositiveInt,
        Field(
            description="Minimum length (in nucleotides) for each padlock probe arm. Each arm must meet this requirement for the probe to pass filtering."
        ),
    ]
    padlock_arm_Tm_min: Annotated[
        NonNegativeFloat,
        Field(description="Minimum acceptable melting temperature (Tm) for padlock arms in degrees Celsius."),
    ]
    padlock_arm_Tm_max: Annotated[
        NonNegativeFloat,
        Field(description="Maximum acceptable melting temperature (Tm) for padlock arms in degrees Celsius."),
    ]
    ligation_region_size: Annotated[
        PositiveInt,
        Field(
            description="Size of the ligation region (in nucleotides) around the ligation site. This parameter is used for seed-based specificity filtering around the junction region where padlock arms meet."
        ),
    ]
    isoform_weight: Annotated[
        WeightT, Field(description="Weight assigned to isoform consensus in the scoring function.")
    ]
    GC_weight: Annotated[WeightT, Field(description="Weight assigned to GC content in the scoring function.")]
    Tm_weight: Annotated[
        WeightT, Field(description="Weight assigned to melting temperature in the scoring function.")
    ]
