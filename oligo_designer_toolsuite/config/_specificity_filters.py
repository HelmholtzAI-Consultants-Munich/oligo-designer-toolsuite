from typing import Annotated, Literal, Union

from pydantic import BaseModel, ConfigDict, Field, NonNegativeInt

from oligo_designer_toolsuite.config._general_models import BlastnHitParameters, BlastnSearchParameters
from oligo_designer_toolsuite.config._types import FilesFastaReferenceDatabaseT


class FilterBaseConfig(BaseModel):
    model_config = ConfigDict(extra="forbid")


class ReadLengthBiasFilterDisabled(FilterBaseConfig):
    enabled: Literal[False]


class ReadLengthBiasFilterEnabled(FilterBaseConfig):
    enabled: Literal[True]
    read_length_bias: Annotated[
        NonNegativeInt,
        Field(
            description="Number of nucleotides from the 5' end of probes to check for read length bias. Probes where the first N bases match exactly with other probes are removed to prevent sequencing read length biases."
        ),
    ]


ReadLengthBiasFilterConfig = Annotated[
    Union[ReadLengthBiasFilterEnabled, ReadLengthBiasFilterDisabled], Field(discriminator="enabled")
]


class CrossHybridizationFilterDisabled(FilterBaseConfig):
    enabled: Literal[False]


class CrossHybridizationFilterEnabled(FilterBaseConfig):
    enabled: Literal[True]
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


CrossHybridizationFilterConfig = Annotated[
    Union[CrossHybridizationFilterEnabled, CrossHybridizationFilterDisabled], Field(discriminator="enabled")
]


class SpecificityBlastnFilterDisabled(FilterBaseConfig):
    enabled: Literal[False]


class SpecificityBlastnFilterEnabled(FilterBaseConfig):
    enabled: Literal[True]
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
    Union[SpecificityBlastnFilterEnabled, SpecificityBlastnFilterDisabled], Field(discriminator="enabled")
]


class VariantFilterDisabled(FilterBaseConfig):
    enabled: Literal[False]


class VariantFilterEnabled(FilterBaseConfig):
    enabled: Literal[True]
    files_vcf_reference_database: Annotated[
        list[str],
        Field(
            description="List of paths to VCF files containing variant information used for filtering probes that overlap with known single nucleotide polymorphisms (SNPs) or other variants. Probes overlapping variants may have reduced specificity.",
        ),
    ]


VariantFilterConfig = Annotated[
    Union[VariantFilterEnabled, VariantFilterDisabled], Field(discriminator="enabled")
]
