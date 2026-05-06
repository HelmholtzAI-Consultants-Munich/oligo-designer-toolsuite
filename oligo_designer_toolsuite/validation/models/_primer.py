from typing import Annotated

from pydantic import BaseModel, ConfigDict, Field, NonNegativeInt, PositiveInt

from oligo_designer_toolsuite.validation._types import (
    DNAT,
    FilesFastaReferenceDatabaseT,
    GCContentMaxT,
    GCContentMinT,
    SecondaryStructuresThresholdDeltaGT,
    TmMaxT,
    TmMinT,
    TSecondaryStructureT,
)
from oligo_designer_toolsuite.validation.models._general import BaseProbabilities, HomopolymerThresholds


class PrimerCycleHCR(BaseModel):
    model_config = ConfigDict(extra="forbid")

    forward_primer_sequence: Annotated[
        DNAT,
        Field(
            description="DNA sequence of the forward primer. Should be a string containing valid nucleotide characters (A, T, G, C). This primer will be placed at the 5' end of the DNA template probe during assembly. Defaults to T7 promoter sequence, change if different sequence desired.",
        ),
    ]
    reverse_primer_sequence: Annotated[
        DNAT,
        Field(
            description="DNA sequence of the reverse primer. Should be a string containing valid nucleotide characters (A, T, G, C). This primer will be placed at the 3' end of the DNA template probe during assembly. Defaults to oligo sequence, change if different sequence desired.",
        ),
    ]


class PrimerFish(BaseModel):
    model_config = ConfigDict(extra="forbid")

    files_fasta_reference_database: Annotated[
        FilesFastaReferenceDatabaseT,
        Field(
            description="List of paths to FASTA files containing reference sequences used for specificity filtering. These files are used to identify off-target binding sites (e.g., whole genome or transcriptome sequences)."
        ),
    ]
    reverse_primer_sequence: Annotated[
        DNAT,
        Field(
            description="DNA sequence of the reverse primer that will be used for complementarity filtering. This prevents the forward and reverse primers from binding to each other. Defaults to reverse complement of 20 nt T7 promoter sequence, change if different sequence desired",
        ),
    ]
    length: Annotated[
        PositiveInt, Field(description="Length (in nucleotides) of each primer sequence to generate.")
    ]
    base_probabilities: Annotated[
        BaseProbabilities,
        Field(
            description="Probabilities of each base for random primer sequence generation. Parameters should be 'A', 'T', 'G', 'C' and values should sum to 1.0."
        ),
    ]
    GC_content_min: GCContentMinT
    GC_content_max: GCContentMaxT
    number_GC_GCclamp: Annotated[
        NonNegativeInt,
        Field(
            description="Minimum number of G or C nucleotides required within the specified number of bases at the 3' end (GC clamp). This improves primer binding stability.",
        ),
    ]
    number_three_prime_base_GCclamp: Annotated[
        NonNegativeInt,
        Field(description="Number of bases from the 3' end to consider for the GC clamp requirement."),
    ]
    homopolymeric_base_n: HomopolymerThresholds
    max_len_selfcomplement: Annotated[
        NonNegativeInt,
        Field(
            description="Maximum allowable length of self-complementary sequences. Primers with longer self-complementary regions can form hairpins and reduce PCR efficiency."
        ),
    ]
    max_len_complement_reverse_primer: Annotated[
        NonNegativeInt,
        Field(
            description="Maximum allowable length of complementarity to the reverse primer sequence. This prevents the forward and reverse primers from binding to each other."
        ),
    ]
    Tm_min: TmMinT
    Tm_max: TmMaxT
    T_secondary_structure: TSecondaryStructureT


class PrimerMerfish(PrimerFish):
    secondary_structures_threshold_deltaG: SecondaryStructuresThresholdDeltaGT


class PrimerSeqFishPlus(PrimerFish):
    pass
