from typing import Annotated

from pydantic import BaseModel, ConfigDict, Field, NonNegativeInt

from oligo_designer_toolsuite.validation._types import LengthMaxT, LengthMinT, TmOptT


class DetectionProbeScrinshot(BaseModel):
    model_config = ConfigDict(extra="forbid")

    min_thymines: Annotated[
        NonNegativeInt,
        Field(
            description="Minimum number of thymine (T) nucleotides required in the detection probe sequence. These thymines will be converted to uracils (U) for UNG cleavage."
        ),
    ]
    length_min: LengthMinT
    length_max: LengthMaxT
    U_distance: Annotated[
        NonNegativeInt,
        Field(
            description="Maximum distance (in nucleotides) allowed between uracil substitutions in the detection oligo. Uracils must be spaced ≤ this distance apart."
        ),
    ]
    Tm_opt: TmOptT
