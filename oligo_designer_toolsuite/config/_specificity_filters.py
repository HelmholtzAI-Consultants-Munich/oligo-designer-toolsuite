from typing import Annotated, Literal

from pydantic import BaseModel, ConfigDict, Field, NonNegativeInt

from oligo_designer_toolsuite.config._general_models import BlastnHitParameters, BlastnSearchParameters
from oligo_designer_toolsuite.config._types import FilesFastaReferenceDatabaseT, VCFReferenceDatabaseT


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
    search_parameters: Annotated[
        BlastnSearchParameters,
        Field(description="Parameters for BLASTN searches used in cross-hybridization filtering."),
    ]
    hit_parameters: Annotated[
        BlastnHitParameters,
        Field(
            description="Parameters for filtering BLASTN hits in cross-hybridization searches. Use either coverage or min_alignment_length."
        ),
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
        Field(
            description="BLASTN search parameters for specificity filtering. These parameters control how BLASTN searches are performed to identify off-target binding sites."
        ),
    ]
    hit_parameters: Annotated[
        BlastnHitParameters,
        Field(
            description="Parameters for filtering BLASTN hits during specificity analysis. Use either coverage or min_alignment_length."
        ),
    ]
    files_fasta_reference_database: FilesFastaReferenceDatabaseT


SpecificityBlastnFilterConfig = Annotated[
    SpecificityBlastnFilterEnabled | SpecificityBlastnFilterDisabled, Field(discriminator="enabled")
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
