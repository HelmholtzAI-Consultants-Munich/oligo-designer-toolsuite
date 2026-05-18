from typing import Annotated, Literal

from pydantic import Field

from oligo_designer_toolsuite.config._general_models import BlastnHitParameters, BlastnSearchParameters
from oligo_designer_toolsuite.config._specificity_filters import (
    SpecificityBlastnFilterDisabled,
    SpecificityBlastnFilterEnabled,
    VariantFilterDisabled,
    VariantFilterEnabled,
)


class OligoSeqSpecificityBlastnFilterDisabled(SpecificityBlastnFilterDisabled):
    pass


class OligoSeqSpecificityBlastnFilterEnabled(SpecificityBlastnFilterEnabled):
    search_parameters: BlastnSearchParameters = BlastnSearchParameters(
        perc_identity=80, strand="minus", word_size=10
    )
    hit_parameters: BlastnHitParameters = BlastnHitParameters(coverage=50)


OligoSeqSpecificityBlastnFilterConfig = Annotated[
    OligoSeqSpecificityBlastnFilterEnabled | OligoSeqSpecificityBlastnFilterDisabled,
    Field(discriminator="enabled"),
]


class OligoSeqVariantFilterDisabled(VariantFilterDisabled):
    pass


class OligoSeqVariantFilterEnabled(VariantFilterEnabled):
    action: Annotated[
        Literal["flag", "filter"],
        Field(description="Action for variant-overlapping oligos: 'flag' (mark only) or 'filter' (exclude)."),
    ] = "flag"


OligoSeqVariantFilterConfig = Annotated[
    VariantFilterEnabled | VariantFilterDisabled, Field(discriminator="enabled")
]
