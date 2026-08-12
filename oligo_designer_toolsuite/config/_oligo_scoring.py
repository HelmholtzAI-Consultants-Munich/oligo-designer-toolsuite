from pydantic import BaseModel, ConfigDict, Field, PositiveInt, model_validator
from typing_extensions import Self

from oligo_designer_toolsuite.config._types import (
    GCContentMaxT,
    GCContentMinT,
    GCContentOptT,
    TmMaxT,
    TmMinT,
    TmOptT,
)


class IndependentSetSelection(BaseModel):
    model_config = ConfigDict(extra="forbid")

    n_sets: PositiveInt = Field(
        description="Number of probe sets to generate per gene.",
    )
    set_size_min: PositiveInt = Field(
        description="Minimum set size. Genes that cannot form a set of at least this size are reported as insufficient.",
    )
    set_size_opt: PositiveInt = Field(
        description="Preferred (optimal) number of oligos per set.",
    )
    distance_between_probes: int = Field(
        description="Required spacing: negative = allow overlap (bases), 0 = adjacent OK, positive = minimum gap (bases).",
    )
    n_attempts_graph: PositiveInt = Field(
        description="Number of randomized graph attempts used to diversify set selection.",
    )
    n_attempts_clique_enum: PositiveInt = Field(
        description="Maximum cliques enumerated per graph attempt.",
    )
    diversification_fraction: float = Field(
        ge=0,
        le=1,
        description="Fraction of oligos removed from the graph in each attempt to increase set diversity.",
    )
    jaccard_opt: float = Field(
        ge=0,
        le=1,
        description="Target maximum Jaccard similarity between selected sets; lower = more distinct sets.",
    )
    jaccard_step: float = Field(
        ge=0,
        le=1,
        description="Step size for relaxing Jaccard threshold when too few sets are found.",
    )

    @model_validator(mode="after")
    def _check_min_max(self) -> Self:
        if self.set_size_min > self.set_size_opt:
            raise ValueError(
                f"'set_size_min' ({self.set_size_min}) must be <= 'set_size_opt' ({self.set_size_opt})"
            )
        return self


class ScoreBaseConfig(BaseModel):
    model_config = ConfigDict(extra="forbid")

    weight: float = Field(description="Weight assigned in the scoring function.")


class UniformDistanceScore(BaseModel):
    model_config = ConfigDict(extra="forbid")

    weight: float = Field(description="Weight for uniform spacing in the set-selection scoring function.")


class IsoformConsensusScore(BaseModel):
    model_config = ConfigDict(extra="forbid")

    weight: float = Field(description="Weight for isoform coverage in the set-selection scoring function.")


class TargetedExonsScore(BaseModel):
    model_config = ConfigDict(extra="forbid")

    weight: float = Field(
        description="Weight for targeted-exon coverage in the set-selection scoring function."
    )
    targeted_exons: list[str] = Field(
        description="List of exon identifiers that should be preferentially targeted by probes, e.g. ['1', '2', '3']. Probes overlapping these exons receive higher scores."
    )


class UTRScore(BaseModel):
    model_config = ConfigDict(extra="forbid")

    weight: float = Field(
        description="Weight for UTR-overlap targeting in the set-selection scoring function."
    )


class GCContentScore(BaseModel):
    model_config = ConfigDict(extra="forbid")

    weight: float = Field(description="Weight for GC content in the set-selection scoring function.")
    GC_content_opt: GCContentOptT


class GCContentScoreNormalized(GCContentScore):
    GC_content_min: GCContentMinT
    GC_content_max: GCContentMaxT

    @model_validator(mode="after")
    def _check_min_max(self) -> Self:
        if not (self.GC_content_min <= self.GC_content_opt <= self.GC_content_max):
            raise ValueError(
                f"Values must be: 'GC_content_min' ({self.GC_content_min}) <= 'GC_content_opt' ({self.GC_content_opt}) <= 'GC_content_max' ({self.GC_content_max})"
            )
        return self


class TmScore(BaseModel):
    model_config = ConfigDict(extra="forbid")

    weight: float = Field(description="Weight for Tm in the set-selection scoring function.")
    Tm_min: TmMinT
    Tm_opt: TmOptT
    Tm_max: TmMaxT

    @model_validator(mode="after")
    def _check_min_max(self) -> Self:
        if not (self.Tm_min <= self.Tm_opt <= self.Tm_max):
            raise ValueError(
                f"Values must be: 'Tm_min' ({self.Tm_min}) <= 'Tm_opt' ({self.Tm_opt}) <= 'Tm_max' ({self.Tm_max})"
            )
        return self
