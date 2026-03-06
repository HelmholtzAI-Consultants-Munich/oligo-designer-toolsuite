from typing import Annotated

from pydantic import Field, NonNegativeFloat, NonNegativeInt, PositiveInt

FastaFileListT = Annotated[list[str], Field(min_length=1)]

FilesFastaDatabaseT = Annotated[
    FastaFileListT,
    Field(
        description="List of paths to FASTA file(s) containing sequences from which probes will be generated. These files should contain genomic regions of interest (e.g., exons, exon-exon junctions). Hint: use the genomic_region_generator pipeline to create FASTA files of genomic regions of interest."
    ),
]

FilesFastaReferenceDatabaseT = Annotated[
    FastaFileListT,
    Field(
        description="List of paths to FASTA file(s) containing reference sequences against which specificity will be evaluated. These typically include the entire genome or transcriptome to identify off-target binding sites. Hint: use the genomic_region_generator pipeline to create FASTA files of genomic regions of interest."
    ),
]

DNAT = Annotated[str, Field(pattern=r"^[ATGC]+$")]

GCContentMinT = Annotated[
    float,
    Field(
        description="Minimum GC content (as a fraction between 0.0 and 1.0) for probes. Probes with GC content below this value will be filtered out.",
        ge=0,
        le=100,
    ),
]
GCContentMaxT = Annotated[
    float,
    Field(
        description="Maximum GC content (as a fraction between 0.0 and 1.0) for probes. Probes with GC content above this value will be filtered out.",
        ge=0,
        le=100,
    ),
]

GCContentOptT = Annotated[
    float,
    Field(
        description="Optimal GC content for target probes, expressed as a fraction between 0.0 and 1.0. Used in scoring to prioritize probes closer to this value.",
        ge=0,
        le=100,
    ),
]

TmMinT = Annotated[
    NonNegativeFloat,
    Field(
        description="Minimum melting temperature (Tm) in degrees Celsius for probes. Probes with calculated Tm below this value will be filtered out."
    ),
]
TmOptT = Annotated[
    NonNegativeFloat,
    Field(
        description="Optimal melting temperature (Tm) for oligos in degrees Celsius. Used in scoring to prioritize probes closer to this value."
    ),
]
TmMaxT = Annotated[
    NonNegativeFloat,
    Field(
        description="Maximum melting temperature (Tm) in degrees Celsius for target probes. Probes with calculated Tm above this value will be filtered out. This value is also used as the optimal Tm target in probe scoring."
    ),
]

TSecondaryStructureT = Annotated[
    PositiveInt,
    Field(
        description="Temperature in degrees Celsius at which to evaluate secondary structure formation (free energy calculation). Secondary structures that form at this temperature can interfere with probe binding."
    ),
]

WeightT = Annotated[float, Field(description="Weight in the efficiency score of the respective measure")]

SecondaryStructuresThresholdDeltaGT = Annotated[
    float,
    Field(
        description="DeltaG threshold (in kcal/mol) for secondary structure stability. Probes with secondary structures having deltaG values more negative (more stable) than this threshold will be filtered out.",
    ),
]

LengthMinT = Annotated[
    NonNegativeInt, Field(description="Minimum length (in nucleotides) for probe sequences.")
]
LengthMaxT = Annotated[
    NonNegativeInt, Field(description="Maximum length (in nucleotides) for probe sequences.")
]

FractionT = Annotated[float, Field(ge=0, le=1)]
