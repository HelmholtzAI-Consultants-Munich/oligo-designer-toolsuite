from typing import Annotated

from pydantic import Field, NonNegativeFloat, NonNegativeInt, PositiveInt

RegionListT = Annotated[
    str | None,
    Field(
        description="Path to file listing gene/region IDs to design probes for. Empty = use all genes from the FASTA."
    ),
]
FastaFileListT = Annotated[list[str], Field(min_length=1)]
FilesFastaDatabaseT = Annotated[
    FastaFileListT,
    Field(
        description="FASTA file(s) from which oligo sequences are generated. Use genomic_region_generator for custom regions.."
    ),
]

FilesFastaReferenceDatabaseT = Annotated[
    FastaFileListT,
    Field(description="FASTA file(s) used as reference for specificity."),
]

VCFReferenceDatabaseT = Annotated[
    list[str],
    Field(
        description="List of paths to VCF files containing variant information used for filtering probes that overlap with known single nucleotide polymorphisms (SNPs) or other variants. Probes overlapping variants may have reduced specificity.",
    ),
]

LengthMinT = Annotated[NonNegativeInt, Field(description="Minimum allowed oligo length (bases)")]
LengthMaxT = Annotated[NonNegativeInt, Field(description="Maximum allowed oligo length (bases).")]

GCContentMinT = Annotated[
    float,
    Field(
        description="Minimum GC content (%). Oligos below this are rejected.",
        ge=0,
        le=100,
    ),
]
GCContentMaxT = Annotated[
    float,
    Field(
        description="Maximum GC content (%). Oligos above this are rejected.",
        ge=0,
        le=100,
    ),
]
GCContentOptT = Annotated[
    float,
    Field(
        description="Target (optimal) GC content (%) for scoring.",
        ge=0,
        le=100,
    ),
]

DRNAT = Annotated[str, Field(pattern=r"^[ATUGC]+$")]

FractionT = Annotated[float, Field(ge=0, le=1)]

TmMinT = Annotated[
    NonNegativeFloat,
    Field(description="Minimum melting temperature (°C). Oligos below this are rejected."),
]
TmMaxT = Annotated[
    NonNegativeFloat,
    Field(description="Maximum melting temperature (°C). Oligos above this are rejected."),
]
TmOptT = Annotated[
    NonNegativeFloat,
    Field(description="Target (optimal) Tm (°C) for scoring."),
]

TSecondaryStructureT = Annotated[
    PositiveInt,
    Field(description="Temperature (°C) at which secondary structure free energy is evaluated."),
]

SecondaryStructuresThresholdDeltaGT = Annotated[
    float,
    Field(
        description="Free-energy threshold (kcal/mol). Oligos with ΔG ≤ thr_DG (stable structure) are rejected.",
    ),
]

HomogeneousPropertiesWeightsT = Annotated[
    dict[str, float],
    Field(
        min_length=1,
        description=(
            "Weights for property homogeneity in the score minimised across the set (weighted sum of variances)."
        ),
    ),
]
