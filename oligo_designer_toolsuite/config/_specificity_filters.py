from typing import Annotated, Literal

from pydantic import BaseModel, ConfigDict, Field, NonNegativeInt

from oligo_designer_toolsuite.config._general_models import BlastnHitParameters, BlastnSearchParameters
from oligo_designer_toolsuite.config._types import FilesFastaReferenceDatabaseT, VCFReferenceDatabaseT

SPECIFICITY_TARGET_DESC = "Ensure oligo specificity against a reference database (BLAST-based)."
SPECIFICITY_READOUT_DESC = "Remove readout probes with significant hits to the reference (BLAST-based)."
SPECIFICITY_PRIMER_DESC = "Ensure primer specificity against a reference database (BLAST-based)."


class FilterBaseConfigEnabled(BaseModel):
    model_config = ConfigDict(extra="forbid")
    enabled: Literal[True]


class FilterBaseConfigDisabled(BaseModel):
    model_config = ConfigDict(extra="ignore")
    enabled: Literal[False]


class ReadLengthBiasFilterDisabled(FilterBaseConfigDisabled):
    pass


class ReadLengthBiasFilterEnabled(FilterBaseConfigEnabled):
    read_length_bias: Annotated[
        NonNegativeInt,
        Field(
            description="Number of nucleotides from the 5' end of probes to check for read length bias. Probes where the first N bases match exactly with other probes are removed to prevent sequencing read length biases."
        ),
    ]


ReadLengthBiasFilterConfig = Annotated[
    ReadLengthBiasFilterEnabled | ReadLengthBiasFilterDisabled, Field(discriminator="enabled")
]


class CrossHybridizationBlastnFilterDisabled(FilterBaseConfigDisabled):
    pass


class CrossHybridizationBlastnFilterEnabled(FilterBaseConfigEnabled):
    """Remove oligos that may cross-hybridize to other probes (BLAST-based)."""

    search_parameters: Annotated[
        BlastnSearchParameters,
        Field(description="BLAST options for the cross-hybridization search."),
    ]
    hit_parameters: Annotated[
        BlastnHitParameters,
        Field(description="Hit criteria. Hits satisfying these lead to oligo rejection."),
    ]


CrossHybridizationBlastnFilterConfig = Annotated[
    CrossHybridizationBlastnFilterEnabled | CrossHybridizationBlastnFilterDisabled,
    Field(discriminator="enabled"),
]


class SpecificityBlastnFilterDisabled(FilterBaseConfigDisabled):
    pass


class SpecificityBlastnFilterEnabled(FilterBaseConfigEnabled):
    search_parameters: Annotated[
        BlastnSearchParameters,
        Field(description="BLAST options for the specificity search."),
    ]
    hit_parameters: Annotated[
        BlastnHitParameters,
        Field(
            description="Hit criteria (either coverage % or mininum alignment length). Hits satisfying these lead to oligo rejection."
        ),
    ]
    files_fasta_reference_database: FilesFastaReferenceDatabaseT


SpecificityBlastnFilterConfig = Annotated[
    SpecificityBlastnFilterEnabled | SpecificityBlastnFilterDisabled, Field(discriminator="enabled")
]


class HybridizationProbesBlastnFilterDisabled(FilterBaseConfigDisabled):
    pass


class HybridizationProbesBlastnFilterEnabled(FilterBaseConfigEnabled):
    search_parameters: Annotated[
        BlastnSearchParameters,
        Field(
            description="BLASTN search parameters for filtering primers that match the assembled hybridization probes."
        ),
    ]
    hit_parameters: Annotated[
        BlastnHitParameters,
        Field(
            description="Parameters for filtering BLASTN hits against the hybridization probes. Use either coverage or min_alignment_length."
        ),
    ]


HybridizationProbesBlastnFilterConfig = Annotated[
    HybridizationProbesBlastnFilterEnabled | HybridizationProbesBlastnFilterDisabled,
    Field(
        discriminator="enabled",
        description="Remove primers that match the assembled hybridization probes (BLAST-based).",
    ),
]


class VariantFilterDisabled(FilterBaseConfigDisabled):
    pass


class VariantFilterEnabled(FilterBaseConfigEnabled):
    files_vcf_reference_database: VCFReferenceDatabaseT
    action: Annotated[
        Literal["flag", "filter"],
        Field(description="Action for variant-overlapping oligos: 'flag' (mark only) or 'filter' (exclude)."),
    ]


VariantFilterConfig = Annotated[VariantFilterEnabled | VariantFilterDisabled, Field(discriminator="enabled")]
