from __future__ import annotations

from pydantic import BaseModel, ConfigDict, PositiveInt

from oligo_designer_toolsuite.validation.models._developer_parameters import (
    DeveloperParametersCycleHCR,
    DeveloperParametersMerfish,
    DeveloperParametersSeqFishPlus,
)
from oligo_designer_toolsuite.validation.models._general import General
from oligo_designer_toolsuite.validation.models._primer import PrimerCycleHCR, PrimerFish, PrimerMerfish
from oligo_designer_toolsuite.validation.models._readout_probes import (
    ReadoutProbeCycleHCR,
    ReadoutProbeMerfish,
    ReadoutProbeSeqFishPlus,
)
from oligo_designer_toolsuite.validation.models._target_probes import (
    TargetProbeCycleHCR,
    TargetProbeMerfish,
    TargetProbeSeqFishPlus,
)


class CycleHCRProbeDesignerConfig(BaseModel):
    model_config = ConfigDict(extra="forbid")
    schema_version: PositiveInt
    general: General
    target_probe: TargetProbeCycleHCR
    readout_probe: ReadoutProbeCycleHCR
    primer: PrimerCycleHCR
    developer_param: DeveloperParametersCycleHCR


class MerfishProbeDesignerConfig(BaseModel):
    model_config = ConfigDict(extra="forbid")
    schema_version: PositiveInt
    general: General
    target_probe: TargetProbeMerfish
    readout_probe: ReadoutProbeMerfish
    primer: PrimerMerfish
    developer_param: DeveloperParametersMerfish


class SeqFishPlusDesignerConfig(BaseModel):
    model_config = ConfigDict(extra="forbid")
    schema_version: PositiveInt
    general: General
    target_probe: TargetProbeSeqFishPlus
    readout_probe: ReadoutProbeSeqFishPlus
    primer: PrimerFish
    developer_param: DeveloperParametersSeqFishPlus
