from math import isclose
from typing import Annotated, Literal

from pydantic import (
    BaseModel,
    ConfigDict,
    Field,
    NonNegativeInt,
    NonPositiveInt,
    PositiveInt,
    model_validator,
)
from typing_extensions import Self


class HomopolymerThresholds(BaseModel):
    model_config = ConfigDict(extra="forbid")
    A: Annotated[PositiveInt | None, Field(default=None)]
    T: Annotated[PositiveInt | None, Field(default=None)]
    C: Annotated[PositiveInt | None, Field(default=None)]
    G: Annotated[PositiveInt | None, Field(default=None)]


class General(BaseModel):
    model_config = ConfigDict(extra="forbid")
    n_jobs: Annotated[
        PositiveInt,
        Field(
            description="number of cores used to run the pipeline and 2*n_jobs +1 of regions that should be stored in cache. If memory consumption of pipeline is too high reduce this number, if a lot of RAM is available increase this number to decrease runtime"
        ),
    ]
    dir_output: Annotated[
        str, Field(description="name of the directory where the output files will be written")
    ]
    write_intermediate_steps: Annotated[
        bool,
        Field(
            default=True,
            description="if true, writes the oligo sequences after each step of the pipeline into a csv file",
        ),
    ]


class OligoSetSelection(BaseModel):
    model_config = ConfigDict(extra="forbid")

    n_attempts_graph: Annotated[
        PositiveInt,
        Field(
            description="Number of randomized graph attempts. In each attempt, a fraction of nodes is randomly removed from the compatibility graph to create diversity; more attempts increase diversity at the cost of runtime.",
        ),
    ]
    n_attempts_clique_enum: Annotated[
        PositiveInt,
        Field(
            description="Maximum number of cliques enumerated per graph attempt. Limits how many cliques are explored before stopping enumeration for the current graph."
        ),
    ]
    diversification_fraction: Annotated[
        float,
        Field(
            ge=0,
            le=1,
            description="Fraction of oligos to remove from the graph per attempt to create diversity in the set selection.",
        ),
    ]
    jaccard_opt: Annotated[
        float,
        Field(
            ge=0,
            le=1,
            description="Optimal maximum Jaccard overlap allowed between selected sets. Lower values enforce more diversity between sets.",
        ),
    ]
    jaccard_step: Annotated[
        float,
        Field(
            ge=0,
            le=1,
            description="Step size used to relax the Jaccard constraint when not enough sets are found.",
        ),
    ]


class TmParametersDetails(BaseModel):
    model_config = ConfigDict(extra="forbid")

    nn_table: Annotated[
        Literal["DNA_NN1", "DNA_NN2", "DNA_NN3", "DNA_NN4"] | None,
        Field(default=None, description="Thermodynamic NN values."),
    ]
    tmm_table: Annotated[
        Literal["DNA_TMM1"] | None,
        Field(default=None, description="Thermodynamic values for terminal mismatches."),
    ]
    imm_table: Annotated[
        Literal["DNA_IMM1"] | None,
        Field(
            default=None,
            description="Thermodynamic values for internal mismatches, may include insosine mismatches.",
        ),
    ]
    de_table: Annotated[
        Literal["DNA_DE1"] | None, Field(default=None, description="Thermodynamic values for dangling ends.")
    ]
    dnac1: Annotated[
        NonNegativeInt | None,
        Field(default=None, description="Concentration of the higher concentrated strand [nM]."),
    ]
    dnac2: Annotated[
        NonNegativeInt | None,
        Field(default=None, description="Concentration of the lower concentrated strand [nM]."),
    ]
    saltcorr: Annotated[
        NonNegativeInt | None,
        Field(
            default=None,
            ge=0,
            le=7,
            description="Salt correction method, see Bio.SeqUtils.MeltingTemp.salt_correction.",
        ),
    ]
    Na: Annotated[NonNegativeInt | None, Field(default=None, description="Concentration of the ions [mM].")]
    K: Annotated[NonNegativeInt | None, Field(default=None, description="Concentration of the ions [mM].")]
    Tris: Annotated[NonNegativeInt | None, Field(default=None, description="Concentration of the ions [mM].")]
    Mg: Annotated[NonNegativeInt | None, Field(default=None, description="Concentration of the ions [mM].")]
    dNTPs: Annotated[
        NonNegativeInt | None, Field(default=None, description="Concentration of the ions [mM].")
    ]

    @model_validator(mode="after")
    def _at_least_one_parameter_provided(self) -> Self:
        data = self.model_dump()
        if all(value is None for value in data.values()):
            raise ValueError("At least one parameter must be provided.")
        return self


# class TmParameters(BaseModel):
#     model_config = ConfigDict(extra="forbid")

#     mode: Literal["biopython_defaults", "custom"] = Field(
#         default="biopython_defaults",
#         description="Should the defaults of the underlying BioPython function (Bio.SeqUtils.MeltingTemp.Tm_NN) be used or custom parameters?",
#     )
#     parameters: TmParametersDetails | None = Field(
#         default=None, description="Required when mode='custom'. Must be omitted otherwise."
#     )

#     @model_validator(mode="after")
#     def _check_mode_and_parameters(self) -> Self:
#         if self.mode == "custom" and self.parameters is None:
#             raise ValueError("TmParameters.parameters must be provided when mode='custom'.")
#         if self.mode == "biopython_defaults" and self.parameters is not None:
#             raise ValueError("TmParameters.parameters must be omitted when mode is 'biopython_defaults'.")
#         return self


class TmParametersBiopythonDefaults(BaseModel):
    model_config = ConfigDict(extra="forbid")

    mode: Literal["biopython_defaults"] = Field(
        default="biopython_defaults",
        description="Should the defaults of the underlying BioPython function (Bio.SeqUtils.MeltingTemp.Tm_NN) be used or custom parameters?",
    )


class TmParametersCustom(BaseModel):
    model_config = ConfigDict(extra="forbid")

    mode: Literal["custom"] = Field(
        default="custom",
        description="Should the defaults of the underlying BioPython function (Bio.SeqUtils.MeltingTemp.Tm_NN) be used or custom parameters?",
    )
    parameters: TmParametersDetails = Field(
        description="Required when mode='custom'. Must be omitted otherwise."
    )


TmParameters = Annotated[TmParametersBiopythonDefaults, TmParametersCustom, Field(discriminator="mode")]


class TmChemCorrectionParametersDetails(BaseModel):
    model_config = ConfigDict(extra="forbid")

    # defaults are from Bio.SeqUtils.MeltingTemp.chem_correction
    DMSO: Annotated[float | None, Field(default=None, ge=0, le=100, description="Percent DMSO")]
    DMSOfactor: Annotated[
        float | None, Field(default=None, description="How much Tm should decrease per percent DMSO")
    ]
    fmd: Annotated[
        float | None,
        Field(default=None, description="Formamide concentration in %(fmdmethod=1) or molar (fmdmethod=2)."),
    ]
    fmdfactor: Annotated[
        float | None, Field(default=None, description="How much Tm should decrease per percent formamide")
    ]
    fmdmethod: Annotated[
        int | None,
        Field(
            default=None,
            ge=1,
            le=2,
            description="Tm = Tm - factor(%formamide) (Default); Tm = Tm + (0.453(f(GC)) - 2.88) x [formamide]",
        ),
    ]
    GC: Annotated[float | None, Field(default=None, ge=0, le=100, description="GC content in percent.")]

    @model_validator(mode="after")
    def _check_fmd_vs_method(self) -> Self:
        # method 1: fmd is percent
        if self.fmdmethod == 1:
            if self.fmd is not None and (self.fmd < 0.0 or self.fmd > 100.0):
                raise ValueError("For fmdmethod=1, fmd must be a percentage in [0, 100].")

        # method 2: fmd is molar, and GC is required by the formula
        elif self.fmdmethod == 2:
            if self.fmd is not None and self.fmd < 0.0:
                raise ValueError("For fmdmethod=2, fmd must be a non-negative molar concentration.")
            if self.GC is None:
                raise ValueError("For fmdmethod=2, GC must be provided (0–100%) for the formula.")

        return self

    @model_validator(mode="after")
    def _at_least_one_parameter_provided(self) -> Self:
        data = self.model_dump()
        if all(value is None for value in data.values()):
            raise ValueError("At least one parameter must be provided.")
        return self


class TmChemCorrectionParameters(BaseModel):
    model_config = ConfigDict(extra="forbid")

    mode: Literal["biopython_defaults", "custom", "disabled"] = Field(
        default="disabled",
        description="Should the defaults of the underlying BioPython function (Bio.SeqUtils.MeltingTemp.chem_correction) be used, custom parameters or chem correction be disabled?",
    )
    parameters: TmChemCorrectionParametersDetails | None = Field(
        default=None, description="Required when mode='custom'. Must be omitted otherwise."
    )

    @model_validator(mode="after")
    def _check_mode_and_parameters(self) -> Self:
        if self.mode == "custom" and self.parameters is None:
            raise ValueError("TmChemCorrectionParameters.parameters must be provided when mode='custom'.")
        if self.mode in ("biopython_defaults", "disabled") and self.parameters is not None:
            raise ValueError(
                "TmChemCorrectionParameters must be omitted when mode is 'biopython_defaults' or 'disabled'."
            )
        return self


class TmSaltCorrectionParametersDetails(BaseModel):
    model_config = ConfigDict(extra="forbid")

    Na: Annotated[NonNegativeInt | None, Field(default=None, description="[mM] of ion")]
    K: Annotated[NonNegativeInt | None, Field(default=None, description="[mM] of ion")]
    Tris: Annotated[NonNegativeInt | None, Field(default=None, description="[mM] of ion")]
    Mg: Annotated[NonNegativeInt | None, Field(default=None, description="[mM] of ion")]
    dNTPs: Annotated[NonNegativeInt | None, Field(default=None, description="[mM] of ion")]
    method: Annotated[
        PositiveInt | None,
        Field(
            default=None,
            ge=1,
            le=7,
            description="Correction method to be applied. Methods 1-4 correct Tm, method 5 corrects deltaS, methods 6 and 7 correct 1/Tm.",
        ),
    ]

    @model_validator(mode="after")
    def _at_least_one_parameter_provided(self) -> Self:
        data = self.model_dump()
        if all(value is None for value in data.values()):
            raise ValueError("At least one parameter must be provided.")
        return self


class TmSaltCorrectionParameters(BaseModel):
    model_config = ConfigDict(extra="forbid")

    mode: Literal["biopython_defaults", "custom", "disabled"] = Field(
        default="disabled",
        description="Should the defaults of the underlying BioPython function (Bio.SeqUtils.MeltingTemp.salt_correction) be used, custom parameters or salt correction be disabled?",
    )
    parameters: TmSaltCorrectionParametersDetails | None = Field(
        default=None, description="Required when mode='custom'. Must be omitted otherwise."
    )

    @model_validator(mode="after")
    def _check_mode_and_parameters(self) -> Self:
        if self.mode == "custom" and self.parameters is None:
            raise ValueError("TmSaltmCorrectionParameters.parameters must be provided when mode='custom'.")
        if self.mode in ("biopython_defaults", "disabled") and self.parameters is not None:
            raise ValueError(
                "TmSaltCorrectionParameters must be omitted when mode is 'biopython_defaults' or 'disabled'."
            )
        return self


class BlastnSearchParameters(BaseModel):
    model_config = ConfigDict(extra="forbid", validate_by_name=True, validate_by_alias=True)

    # use exclude_if=lambda v: v is None on a Field level here and not exlcude_none during model dumpimg,
    # because for other config arguments, the None is actually needed
    # but for blastn we don't want to provide default parameters but have the defaults handled directly
    # by blastn (as there can be quite complicated dependencies between arguments)
    # don't allow
    # -h
    # -help
    # -version
    # -query
    query_loc: str | None = Field(
        default=None,
        alias="-query_loc",
        description="Location on the query sequence in 1-based offsets (Format: start-stop).",
    )
    strand: Literal["plus", "minus", "both"] | None = Field(
        default=None,
        alias="-strand",
        description="Query strand(s) to search against database/subject. Choice of both, minus, or plus.",
    )
    task: Literal["megablast", "dc-megablast", "blastn", "blastn-short", "rmblastn"] | None = Field(
        default=None, alias="-task", description="Supported tasks."
    )
    # don't allow
    # -db
    # -out
    evalue: float | None = Field(
        default=None,
        alias="-evalue",
        description="Expectation value (E) threshold for saving hits. Default = 10 (1000 for blastn-short)",
    )
    word_size: NonNegativeInt | None = Field(
        default=None,
        alias="-word_size",
        ge=4,
        description="Length of initial exact match.",
    )
    gapopen: NonNegativeInt | None = Field(default=None, alias="-gapopen", description="Cost to open a gap.")
    gapextend: NonNegativeInt | None = Field(
        default=None,
        alias="-gapextend",
        description="Cost to extend a gap.",
    )
    penalty: NonPositiveInt | None = Field(
        default=None,
        alias="-penalty",
        description="Penalty for a nucleotide mismatch.",
    )
    reward: NonNegativeInt | None = Field(
        default=None,
        alias="-reward",
        description="Reward for a nucleotide match.",
    )
    # don't allow use_index/index_name as another file would be needed
    # don't allow subject/subject_loc as another file would be needed
    # don't allow
    # -outfmt
    # -show_gis
    num_descriptions: NonNegativeInt | None = Field(
        default=None,
        alias="-num_descriptions",
        description="Number of database sequences to show one-line descriptions for.",
    )
    num_alignments: NonNegativeInt | None = Field(
        default=None,
        alias="-num_alignments",
        description="Number of database sequences to show alignments for.",
    )
    # don't allow
    # line_length
    # html
    sorthits: NonNegativeInt | None = Field(
        default=None,
        alias="-sorthits",
        le=4,
        description="Sorting option for hits.",
    )
    sorthsps: NonNegativeInt | None = Field(
        default=None,
        alias="-sorthsps",
        le=4,
        description="Sorting option for hps.",
    )
    dust: str | None = Field(
        default=None,
        alias="-dust",
        description="Filter query sequence with dust.",
    )
    # don't allow
    # filtering_db as another file would be needed
    # window_masker_taxid
    # window_masker_db as then another file would be needed
    soft_masking: bool | None = Field(
        default=None,
        alias="-soft_masking",
        description="Apply filtering locations as soft masks (i.e., only for finding initial matches).",
    )
    lcase_masking: Literal[""] | None = Field(
        default=None,
        alias="-lcase_masking",
        description="Use lower case filtering in query and subject sequence(s).",
    )
    # don't allow the following options as additional files would be needed
    # gilist
    # seqidlist
    # negative_gilist
    # negative_seqidlist
    # taxids (theoretically no extra file needed, but we blast against only the target genome)
    # negative_taxids
    # taxidlist
    # negative_taxidlist
    # no_taxid_expansion
    # entrez_query
    db_soft_mask: str | None = Field(
        default=None,
        alias="-db_soft_mask",
        description="Filtering algorithm ID to apply to the BLAST database as soft mask (i.e., only for finding initial matches).",
    )
    db_hard_mask: str | None = Field(
        default=None,
        alias="-db_hard_mask",
        description="Filtering algorithm ID to apply to the BLAST database as hard mask (i.e., sequence is masked for all phases of search).",
    )
    perc_identity: float | None = Field(
        default=None,
        alias="-perc_identity",
        description="Percent identity cutoff.",
        ge=0,
        le=100,
    )
    qcov_hsp_perc: float | None = Field(
        default=None,
        alias="-qcov_hsp_perc",
        description="Percent query coverage per hsp.",
        ge=0,
        le=100,
    )
    max_hsps: PositiveInt | None = Field(
        default=None,
        alias="-max_hsps",
        description="Set maximum number of HSPs per subject sequence to save for each query.",
    )
    culling_limit: NonNegativeInt | None = Field(
        default=None,
        alias="-culling_limit",
        description="If the query range of a hit is enveloped by that of at least this many higher-scoring hits, delete the hit",
    )
    best_hit_overhang: float | None = Field(
        default=None,
        alias="-best_hit_overhang",
        description="Best Hit algorithm overhang value (recommended value: 0.1).",
        gt=0,
        lt=0.5,
    )
    best_hit_score_edge: float | None = Field(
        default=None,
        alias="-best_hit_score_edge",
        description="Best Hit algorithm score edge value (recommended value: 0.1)",
        gt=0,
        lt=0.5,
    )
    subject_besthit: Literal[""] | None = Field(
        default=None,
        alias="-subject_besthit",
        description="Turn on best hit per subject sequence.",
    )
    max_target_seqs: PositiveInt | None = Field(
        default=None,
        alias="-max_target_seqs",
        description="Maximum number of aligned sequences to keep.",
    )
    template_type: Literal["coding", "coding_and_optimal", "optimal"] | None = Field(
        default=None,
        alias="-template_type",
        description="Discontiguous MegaBLAST template type. Allowed values are coding, optimal and coding_and_optimal.",
    )
    # template_length is actually int, but only 3 values, therefore implemented as literal
    template_length: Literal["16", "18", "21"] | None = Field(
        default=None,
        alias="-template_length",
        description="Discontiguous MegaBLAST template length.",
    )
    db_size: int | None = Field(
        default=None,
        alias="-db_size",
        description="Effective length of the database.",
    )
    searchsp: NonNegativeInt | None = Field(
        default=None,
        alias="-searchsp",
        description="Effective length of the search space.",
    )
    # don't allow because extra file needed
    # import_search_strategy
    # export_search_strategy
    xdrop_ungap: float | None = Field(
        default=None,
        alias="-xdrop_ungap",
        description="X-dropoff value (in bits) for ungapped extensions.",
    )
    xdrop_gap: float | None = Field(
        default=None,
        alias="-xdrop_gap",
        description="X-dropoff value (in bits) for preliminary gapped extensions.",
    )
    xdrop_gap_final: float | None = Field(
        default=None,
        alias="-xdrop_gap_final",
        description="X-dropoff value (in bits) for final gapped alignment.",
    )
    no_greedy: Literal[""] | None = Field(
        default=None,
        alias="-no_greedy",
        description="Use non-greedy dynamic programming extension.",
    )
    min_raw_gapped_score: int | None = Field(
        default=None,
        alias="-min_raw_gapped_score",
        description="Minimum raw gapped score to keep an alignment in the preliminary gapped and trace-back stages. Normally set based upon expect value.",
    )
    ungapped: Literal[""] | None = Field(
        default=None,
        alias="-ungapped",
        description="Perform ungapped alignment only?",
    )
    window_size: NonNegativeInt | None = Field(
        default=None,
        alias="-window_size",
        description="Multiple hits window size, use 0 to specify 1-hit algorithm.",
    )
    off_diagonal_range: NonNegativeInt | None = Field(
        default=None,
        alias="-off_diagonal_range",
        description="Number of off-diagonals to search for the 2nd hit, use 0 to turn off.",
    )
    # don't allow
    # parse_deflines
    # num_threads
    # mt_mode
    # remote


class BlastnHitParameters(BaseModel):
    model_config = ConfigDict(extra="forbid")

    coverage: Annotated[
        float | None,
        Field(
            default=None,
            ge=0,
            le=100,
            description="Coverage in %, alternatively, min_alignment_length can be used",
        ),
    ]
    min_alignment_length: Annotated[
        NonNegativeInt | None,
        Field(
            default=None,
            description="Number of nucleotides for alignment, alternatively, coverage can be used",
        ),
    ]

    @model_validator(mode="after")
    def _check_mutually_exclusive(self) -> Self:
        if (self.coverage is None) == (self.min_alignment_length is None):
            raise ValueError("Exactly one of 'coverage' or 'min_alignment_length' must be set.")
        return self


class BaseProbabilities(BaseModel):
    model_config = ConfigDict(extra="forbid")

    A: float = Field(default=0.25, ge=0, le=1)
    C: float = Field(default=0.25, ge=0, le=1)
    G: float = Field(default=0.25, ge=0, le=1)
    T: float = Field(default=0.25, ge=0, le=1)

    @model_validator(mode="after")
    def _check_sums_up_to_1(self) -> Self:
        sum_probabilities = self.A + self.C + self.G + self.T
        if not isclose(sum_probabilities, 1):
            raise ValueError("The probabilities for all 4 bases needs to sum up to 1.")
        return self


class OligoPropertyWeights(BaseModel):
    model_config = ConfigDict(extra="forbid")

    length_oligo: float | None = Field(default=None)
    GC_content_oligo: float | None = Field(default=None)
    TmNN_oligo: float | None = Field(default=None)
    DG_secondary_structure_oligo: float | None = Field(default=None)
    length_selfcomplement_oligo: float | None = Field(default=None)
