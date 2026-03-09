from __future__ import annotations

from typing import Annotated

from pydantic import BaseModel, ConfigDict, Field, PositiveInt

from oligo_designer_toolsuite.validation._types import DirOutputT, ExonExonJunctionBlockSizeT
from oligo_designer_toolsuite.validation.models._detection_probes import DetectionProbeScrinshot
from oligo_designer_toolsuite.validation.models._developer_parameters import (
    DeveloperParametersCycleHCR,
    DeveloperParametersMerfish,
    DeveloperParametersOligoSeq,
    DeveloperParametersScrinshot,
    DeveloperParametersSeqFishPlus,
)
from oligo_designer_toolsuite.validation.models._general import General, GenomicRegions, SourceConfigs
from oligo_designer_toolsuite.validation.models._primer import PrimerCycleHCR, PrimerFish, PrimerMerfish
from oligo_designer_toolsuite.validation.models._readout_probes import (
    ReadoutProbeCycleHCR,
    ReadoutProbeMerfish,
    ReadoutProbeSeqFishPlus,
)
from oligo_designer_toolsuite.validation.models._target_probes import (
    TargetProbeCycleHCR,
    TargetProbeMerfish,
    TargetProbeOligoSeq,
    TargetProbeScrinshot,
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


class SeqFishPlusProbeDesignerConfig(BaseModel):
    model_config = ConfigDict(extra="forbid")
    schema_version: PositiveInt
    general: General
    target_probe: TargetProbeSeqFishPlus
    readout_probe: ReadoutProbeSeqFishPlus
    primer: PrimerFish
    developer_param: DeveloperParametersSeqFishPlus


class OligoSeqProbeDesignerConfig(BaseModel):
    model_config = ConfigDict(extra="forbid")
    schema_version: PositiveInt
    general: General
    target_probe: TargetProbeOligoSeq
    developer_param: DeveloperParametersOligoSeq


class ScrinshotProbeDesignerConfig(BaseModel):
    model_config = ConfigDict(extra="forbid")
    schema_version: PositiveInt
    general: General
    target_probe: TargetProbeScrinshot
    detection_probe: DetectionProbeScrinshot
    developer_param: DeveloperParametersScrinshot


class GenomicRegionGeneratorConfig(BaseModel):
    model_config = ConfigDict(extra="forbid")
    schema_version: PositiveInt
    dir_output: DirOutputT
    source: Annotated[SourceConfigs, Field(description="Parameters for genome and gene annotation.")]
    genomic_regions: Annotated[
        GenomicRegions,
        Field(
            description="List of genomic regions that should be generated, set the genomic regions you want to generate to true."
        ),
    ]
    exon_exon_junction_block_size: ExonExonJunctionBlockSizeT
