from typing import Annotated, Literal

from pydantic import BaseModel, ConfigDict, Field, NonNegativeInt

from oligo_designer_toolsuite.config._general_models import (
    BlastnHitParameters,
    BlastnHitParametersCoverage,
    BlastnHitParametersMinAlignmentLength,
    BlastnSearchParameters,
)
from oligo_designer_toolsuite.config._types import VCFReferenceDatabaseT

SPECIFICITY_TARGET_DESC = "Ensure oligo specificity against a reference database (BLAST-based)."
SPECIFICITY_READOUT_DESC = "Remove readout probes with significant hits to the reference (BLAST-based)."
SPECIFICITY_PRIMER_DESC = "Ensure primer specificity against a reference database (BLAST-based)."
CROSS_HYBRIDIZATION_DESC = "Remove oligos that may cross-hybridize to other probes (BLAST-based)."
HYBRIDIZATION_PROBES_DESC = "Remove primers that match the assembled hybridization probes (BLAST-based)."


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
        Field(
            description="BLAST options for the cross-hybridization search.",
            json_schema_extra={"x-collapsed": True},
        ),
    ]
    hit_parameters: Annotated[
        BlastnHitParameters,
        Field(description="Hit criteria. Hits satisfying these lead to oligo rejection."),
    ]


CrossHybridizationBlastnFilterConfig = Annotated[
    CrossHybridizationBlastnFilterEnabled | CrossHybridizationBlastnFilterDisabled,
    Field(
        discriminator="enabled",
        description=CROSS_HYBRIDIZATION_DESC,
    ),
]


class SpecificityBlastnFilterDisabled(FilterBaseConfigDisabled):
    pass


class SpecificityBlastnFilterEnabled(FilterBaseConfigEnabled):
    search_parameters: Annotated[
        BlastnSearchParameters,
        Field(
            description="BLAST options for the specificity search.",
            json_schema_extra={"x-collapsed": True},
        ),
    ]
    hit_parameters: Annotated[
        BlastnHitParameters,
        Field(
            description="Hit criteria (either coverage % or mininum alignment length). Hits satisfying these lead to oligo rejection."
        ),
    ]


SpecificityBlastnFilterConfig = Annotated[
    SpecificityBlastnFilterEnabled | SpecificityBlastnFilterDisabled, Field(discriminator="enabled")
]


class HybridizationProbesBlastnFilterDisabled(FilterBaseConfigDisabled):
    pass


class HybridizationProbesBlastnFilterEnabled(FilterBaseConfigEnabled):
    search_parameters: Annotated[
        BlastnSearchParameters,
        Field(
            description="BLASTN search parameters for filtering primers that match the assembled hybridization probes.",
            json_schema_extra={"x-collapsed": True},
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
        description=HYBRIDIZATION_PROBES_DESC,
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


# Variants that set the hit criterion on `hit_parameters` itself, not on an enclosing config.
#
# A JSON Schema consumer picks the branch from a default on the field; one further out is not
# read, so it falls back to the first branch while the outer value still merges in. That made
# `value=15` meant as 15 bases arrive as 15 % coverage -- silently, both being plain numbers.
# Coverage variants are needed for the same reason: otherwise a pipeline is correct only while
# coverage stays declared first.
#
# Values are placeholders, overridden per call site. `description` is restated because
# redeclaring a field replaces its whole FieldInfo.


class SpecificityBlastnFilterMinAlignmentEnabled(SpecificityBlastnFilterEnabled):
    hit_parameters: BlastnHitParameters = Field(
        default=BlastnHitParametersMinAlignmentLength(value=15),
        description=SpecificityBlastnFilterEnabled.model_fields["hit_parameters"].description,
    )


SpecificityBlastnFilterMinAlignmentConfig = Annotated[
    SpecificityBlastnFilterMinAlignmentEnabled | SpecificityBlastnFilterDisabled,
    Field(discriminator="enabled"),
]


class SpecificityBlastnFilterCoverageEnabled(SpecificityBlastnFilterEnabled):
    hit_parameters: BlastnHitParameters = Field(
        default=BlastnHitParametersCoverage(value=50),
        description=SpecificityBlastnFilterEnabled.model_fields["hit_parameters"].description,
    )


SpecificityBlastnFilterCoverageConfig = Annotated[
    SpecificityBlastnFilterCoverageEnabled | SpecificityBlastnFilterDisabled,
    Field(discriminator="enabled"),
]


class CrossHybridizationBlastnFilterMinAlignmentEnabled(CrossHybridizationBlastnFilterEnabled):
    hit_parameters: BlastnHitParameters = Field(
        default=BlastnHitParametersMinAlignmentLength(value=17),
        description=CrossHybridizationBlastnFilterEnabled.model_fields["hit_parameters"].description,
    )


CrossHybridizationBlastnFilterMinAlignmentConfig = Annotated[
    CrossHybridizationBlastnFilterMinAlignmentEnabled | CrossHybridizationBlastnFilterDisabled,
    Field(discriminator="enabled", description=CROSS_HYBRIDIZATION_DESC),
]


class CrossHybridizationBlastnFilterCoverageEnabled(CrossHybridizationBlastnFilterEnabled):
    hit_parameters: BlastnHitParameters = Field(
        default=BlastnHitParametersCoverage(value=50),
        description=CrossHybridizationBlastnFilterEnabled.model_fields["hit_parameters"].description,
    )


CrossHybridizationBlastnFilterCoverageConfig = Annotated[
    CrossHybridizationBlastnFilterCoverageEnabled | CrossHybridizationBlastnFilterDisabled,
    Field(discriminator="enabled", description=CROSS_HYBRIDIZATION_DESC),
]


class HybridizationProbesBlastnFilterMinAlignmentEnabled(HybridizationProbesBlastnFilterEnabled):
    hit_parameters: BlastnHitParameters = Field(
        default=BlastnHitParametersMinAlignmentLength(value=11),
        description=HybridizationProbesBlastnFilterEnabled.model_fields["hit_parameters"].description,
    )


HybridizationProbesBlastnFilterMinAlignmentConfig = Annotated[
    HybridizationProbesBlastnFilterMinAlignmentEnabled | HybridizationProbesBlastnFilterDisabled,
    Field(discriminator="enabled", description=HYBRIDIZATION_PROBES_DESC),
]
