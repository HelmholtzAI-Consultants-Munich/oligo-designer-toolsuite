from pydantic import BaseModel, ConfigDict, Field, PositiveInt

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
        description="Number of oligo sets to generate per region. Multiple sets allow for redundancy and selection of the best-performing set based on scoring criteria.",
    )
    set_size_min: PositiveInt = Field(
        description="Minimum size (number of probes) required for each oligo set. Sets with fewer probes than this value will be rejected, and regions that cannot generate sets meeting this minimum will be removed. Regions with less oligos will be filtered out and stored in regions_with_insufficient_oligos_for_db_probes.",
    )
    set_size_opt: PositiveInt = Field(
        description="Optimal size (number of probes) for each oligo set. The set selection algorithm will attempt to generate sets of this size, but may produce sets with fewer probes if constraints cannot be met.",
    )
    distance_between_target_probes: int = Field(
        description="How much overlap should be allowed between oligos, e.g. if oligos can overlap x bases choose -x, if oligos can be next to one another choose 0, if oligos should be x bases apart choose x.",
    )
    n_attempts_graph: PositiveInt = Field(
        description="Number of randomized graph attempts. In each attempt, a fraction of nodes is randomly removed from the compatibility graph to create diversity; more attempts increase diversity at the cost of runtime.",
    )
    n_attempts_clique_enum: PositiveInt = Field(
        description="Maximum number of cliques enumerated per graph attempt. Limits how many cliques are explored before stopping enumeration for the current graph.",
    )
    diversification_fraction: float = Field(
        ge=0,
        le=1,
        description="Fraction of oligos to remove from the graph per attempt to create diversity in the set selection.",
    )
    jaccard_opt: float = Field(
        ge=0,
        le=1,
        description="Optimal maximum Jaccard overlap allowed between selected sets. Lower values enforce more diversity between sets.",
    )
    jaccard_step: float = Field(
        ge=0,
        le=1,
        description="Step size used to relax the Jaccard constraint when not enough sets are found.",
    )


class ScoreBaseConfig(BaseModel):
    model_config = ConfigDict(extra="forbid")

    weight: float = Field(description="Weight assigned in the scoring function.")


class UniformDistanceScore(ScoreBaseConfig):
    pass


class IsoformConsensusScore(ScoreBaseConfig):
    pass


class TargetedExonsScore(ScoreBaseConfig):
    targeted_exons: list[str] = Field(
        description="List of exon identifiers that should be preferentially targeted by probes. Probes overlapping these exons receive higher scores.",
        default=["1", "2", "3"],
    )


class GCContentScore(ScoreBaseConfig):
    GC_content_min: GCContentMinT
    GC_content_opt: GCContentOptT
    GC_content_max: GCContentMaxT


class TmScore(ScoreBaseConfig):
    Tm_min: TmMinT
    Tm_opt: TmOptT
    Tm_max: TmMaxT
