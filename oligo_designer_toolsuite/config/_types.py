from typing import Annotated

from pydantic import Field, NonNegativeFloat, NonNegativeInt, PositiveInt

RegionListT = Annotated[
    str | None,
    Field(
        description="File with a list the genes used to generate the probe sequences, leave empty if all the genes are used. For examples, see 'data/genes'."
    ),
]
FastaFileListT = Annotated[list[str], Field(min_length=1)]
FilesFastaDatabaseT = Annotated[
    FastaFileListT,
    Field(
        description="List of paths to FASTA file(s) containing sequences from which probes will be generated. These files should contain genomic regions of interest (e.g., exons, exon-exon junctions). Hint: use the genomic_region_generator pipeline to create FASTA files of genomic regions of interest. For examples, see 'data/genomic_regions'."
    ),
]

FilesFastaReferenceDatabaseT = Annotated[
    FastaFileListT,
    Field(
        description="List of paths to FASTA file(s) containing reference sequences against which specificity will be evaluated. These typically include the entire genome or transcriptome to identify off-target binding sites. Hint: use the genomic_region_generator pipeline to create FASTA files of genomic regions of interest. For examples, see 'data/genomic_regions'."
    ),
]

VCFReferenceDatabaseT = Annotated[
    list[str],
    Field(
        description="List of paths to VCF files containing variant information used for filtering probes that overlap with known single nucleotide polymorphisms (SNPs) or other variants. Probes overlapping variants may have reduced specificity.",
    ),
]

LengthMinT = Annotated[
    NonNegativeInt, Field(description="Minimum length (in nucleotides) for oligo sequences.")
]
LengthMaxT = Annotated[
    NonNegativeInt, Field(description="Maximum length (in nucleotides) for probe sequences.")
]

GCContentMinT = Annotated[
    float,
    Field(
        description="Minimum GC content (in %) for probes. Probes with GC content below this value will be filtered out/rejected.",
        ge=0,
        le=100,
    ),
]
GCContentMaxT = Annotated[
    float,
    Field(
        description="Maximum GC content (in %) for probes. Probes with GC content above this value will be filtered out/rejected.",
        ge=0,
        le=100,
    ),
]
GCContentOptT = Annotated[
    float,
    Field(
        description="Optimal GC content (in %) for target probes. Used in scoring to prioritize probes closer to this value.",
        ge=0,
        le=100,
    ),
]

DRNAT = Annotated[str, Field(pattern=r"^[ATUGC]+$")]

FractionT = Annotated[float, Field(ge=0, le=1)]

TmMinT = Annotated[
    NonNegativeFloat,
    Field(
        description="Minimum melting temperature (Tm) in degrees Celsius for probes. Probes with calculated Tm below this value will be filtered out."
    ),
]
TmMaxT = Annotated[
    NonNegativeFloat,
    Field(
        description="Maximum melting temperature (Tm) in degrees Celsius for probes. Probes with calculated Tm above this value will be filtered out. This value is also used as the optimal Tm target in probe scoring."
    ),
]
TmOptT = Annotated[
    NonNegativeFloat,
    Field(
        description="Optimal melting temperature (Tm) for oligos in degrees Celsius. Used in scoring to prioritize probes closer to this value."
    ),
]

TSecondaryStructureT = Annotated[
    PositiveInt,
    Field(
        description="Temperature in degrees Celsius at which to evaluate secondary structure formation (free energy calculation). Secondary structures that form at this temperature can interfere with probe binding."
    ),
]

SecondaryStructuresThresholdDeltaGT = Annotated[
    float,
    Field(
        description="DeltaG threshold (in kcal/mol) for secondary structure stability. Probes with secondary structures having deltaG values more negative (more stable) than this threshold will be filtered out.",
    ),
]
