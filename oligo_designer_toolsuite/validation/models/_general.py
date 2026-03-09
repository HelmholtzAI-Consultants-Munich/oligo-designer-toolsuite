from math import isclose
from typing import Annotated, Literal

from pydantic import (
    AliasChoices,
    BaseModel,
    ConfigDict,
    Field,
    NonNegativeInt,
    NonPositiveInt,
    PositiveInt,
    model_validator,
)
from typing_extensions import Self

from oligo_designer_toolsuite.validation._types import FractionT


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


class TmChemCorrectionParametersDisabled(BaseModel):
    model_config = ConfigDict(extra="forbid")

    mode: Literal["disabled"] = Field(
        description="Should the defaults of the underlying BioPython function (Bio.SeqUtils.MeltingTemp.chem_correction) be used, custom parameters or chem correction be disabled?",
    )


class TmChemCorrectionParametersBiopythonDefaults(BaseModel):
    model_config = ConfigDict(extra="forbid")

    mode: Literal["biopython_defaults"] = Field(
        description="Should the defaults of the underlying BioPython function (Bio.SeqUtils.MeltingTemp.chem_correction) be used, custom parameters or chem correction be disabled?",
    )


class TmChemCorrectionParametersCustom(BaseModel):
    model_config = ConfigDict(extra="forbid")

    mode: Literal["custom"] = Field(
        description="Should the defaults of the underlying BioPython function (Bio.SeqUtils.MeltingTemp.chem_correction) be used, custom parameters or chem correction be disabled?",
    )
    parameters: TmChemCorrectionParametersDetails = Field(
        description="Required when mode='custom'. Must be omitted otherwise."
    )


TmChemCorrectionParameters = Annotated[
    TmChemCorrectionParametersDisabled
    | TmChemCorrectionParametersBiopythonDefaults
    | TmChemCorrectionParametersCustom,
    Field(discriminator="mode"),
]


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


class TmSaltCorrectionParametersDisabled(BaseModel):
    model_config = ConfigDict(extra="forbid")

    mode: Literal["disabled"] = Field(
        description="Should the defaults of the underlying BioPython function (Bio.SeqUtils.MeltingTemp.salt_correction) be used, custom parameters or salt correction be disabled?",
    )


class TmSaltCorrectionParametersBiopythonDefaults(BaseModel):
    model_config = ConfigDict(extra="forbid")

    mode: Literal["biopython_defaults"] = Field(
        description="Should the defaults of the underlying BioPython function (Bio.SeqUtils.MeltingTemp.salt_correction) be used, custom parameters or salt correction be disabled?",
    )


class TmSaltCorrectionParametersCustom(BaseModel):
    model_config = ConfigDict(extra="forbid")

    mode: Literal["custom"] = Field(
        description="Should the defaults of the underlying BioPython function (Bio.SeqUtils.MeltingTemp.salt_correction) be used, custom parameters or salt correction be disabled?",
    )
    parameters: TmSaltCorrectionParametersDetails = Field(
        description="Required when mode='custom'. Must be omitted otherwise."
    )


TmSaltCorrectionParameters = Annotated[
    TmSaltCorrectionParametersDisabled
    | TmSaltCorrectionParametersBiopythonDefaults
    | TmSaltCorrectionParametersCustom,
    Field(discriminator="mode"),
]


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


class BowtieSearchParameters(BaseModel):
    model_config = ConfigDict(extra="forbid", validate_by_name=True, validate_by_alias=True)

    # don't allow any input options

    v: int | None = Field(
        default=None,
        alias="-v",
        description="Report alignments with at most <int> mismatches. -e and -l options are ignored and quality values have no effect on what alignments are valid. -v is mutually exclusive with -n.",
    )
    seedmms: int | None = Field(
        default=None,
        ge=0,
        le=3,
        alias="-n",
        validation_alias=AliasChoices("-n", "--seedmms"),
        description="Maximum number of mismatches permitted in the 'seed', i.e. the first L base pairs of the read (where L is set with -l/--seedlen). This may be 0, 1, 2 or 3 and the default is 2. This option is mutually exclusive with the -v option.",
    )
    maqerr: NonNegativeInt | None = Field(
        default=None,
        alias="-e",
        validation_alias=AliasChoices("-e", "--maqerr"),
        description="Maximum permitted total of quality values at all mismatched read positions throughout the entire alignment, not just in the 'seed'. The default is 70. Like Maq, bowtie rounds quality values to the nearest 10 and saturates at 30; rounding can be disabled with --nomaqround.",
    )
    seedlen: NonNegativeInt | None = Field(
        default=None,
        ge=5,
        alias="-l",
        validation_alias=AliasChoices("-l", "--seedlen"),
        description="The 'seed length'; i.e., the number of bases on the high-quality end of the read to which the -n ceiling applies. The lowest permitted setting is 5 and the default is 28. bowtie is faster for larger values of -l.",
    )
    nomaqround: Literal[""] | None = Field(
        default=None,
        alias="--nomaqround",
        description="Maq accepts quality values in the Phred quality scale, but internally rounds values to the nearest 10, with a maximum of 30. By default, bowtie also rounds this way. --nomaqround prevents this rounding in bowtie.",
    )
    nofw: Literal[""] | None = Field(
        default=None,
        alias="--nofw",
        description=(
            "If --nofw is specified, bowtie will not attempt to align against the forward reference strand."
        ),
    )
    norc: Literal[""] | None = Field(
        default=None,
        alias="--norc",
        description=(
            "If --norc is specified, bowtie will not attempt to align against the reverse-complement reference strand."
        ),
    )
    maxbts: PositiveInt | None = Field(
        default=None,
        alias="--maxbts",
        description=(
            "The maximum number of backtracks permitted when aligning a read in -n 2 or -n 3 mode (default: 125 without --best, 800 with --best). A 'backtrack' is the introduction of a speculative substitution into the alignment. Without this limit, the default parameters will sometimes require that bowtie try 100s or 1,000s of backtracks to align a read, especially if the read has many low-quality bases and/or has no valid alignments, slowing bowtie down significantly. However, this limit may cause some valid alignments to be missed. Higher limits yield greater sensitivity at the expensive of longer running times. See also: -y/--tryhard."
        ),
    )
    tryhard: Literal[""] | None = Field(
        default=None,
        alias="-y",
        validation_alias=AliasChoices("-y", "--tryhard"),
        description=(
            "Try as hard as possible to find valid alignments when they exist, including paired-end alignments. This is equivalent to specifying very high values for the --maxbts and --pairtries options. This mode is generally much slower than the default settings, but can be useful for certain problems. This mode is slower when (a) the reference is very repetitive, (b) the reads are low quality, or (c) not many reads have valid alignments."
        ),
    )
    chunkmbs: PositiveInt | None = Field(
        default=None,
        alias="--chunkmbs",
        description=(
            "The number of megabytes of memory a given thread is given to store path descriptors in --best mode. Best-first search must keep track of many paths at once to ensure it is always extending the path with the lowest cumulative cost. Bowtie tries to minimize the memory impact of the descriptors, but they can still grow very large in some cases. If you receive an error message saying that chunk memory has been exhausted in --best mode, try adjusting this parameter up to dedicate more memory to the descriptors. Default: 64."
        ),
    )
    reads_per_batch: PositiveInt | None = Field(
        default=None,
        alias="--reads-per-batch",
        description=(
            "Part of bowtie's batch parsing and used to specify the number of reads that bowtie will consume from the input file at once. Default: 16"
        ),
    )
    k: PositiveInt | None = Field(
        default=None,
        alias="-k",
        description=(
            "Report up to <int> valid alignments per read or pair (default: 1). "
            "Validity of alignments is determined by the alignment policy (combined effects of -n, -v, -l, and -e). "
            "If more than one valid alignment exists and the --best and --strata options are specified, "
            "then only those alignments belonging to the best alignment 'stratum' will be reported. "
            "Bowtie is designed to be very fast for small -k but bowtie can become significantly slower as -k increases. "
            "If you would like to use Bowtie for larger values of -k, consider building an index with a denser suffix-array sample, "
            "i.e. specify a smaller -o/--offrate when invoking bowtie-build for the relevant index "
            "(see the Performance tuning section for details)."
        ),
    )
    all: Literal[""] | None = Field(
        default=None,
        alias="-a",
        validation_alias=AliasChoices("-a", "--all"),
        description=(
            "Report all valid alignments per read or pair (default: off). "
            "Validity of alignments is determined by the alignment policy (combined effects of -n, -v, -l, and -e). "
            "If more than one valid alignment exists and the --best and --strata options are specified, "
            "then only those alignments belonging to the best alignment 'stratum' will be reported. "
            "Bowtie is designed to be very fast for small -k but bowtie can become significantly slower if -a/--all is specified. "
            "If you would like to use Bowtie with -a, consider building an index with a denser suffix-array sample, "
            "i.e. specify a smaller -o/--offrate when invoking bowtie-build for the relevant index "
            "(see the Performance tuning section for details)."
        ),
    )
    m: PositiveInt | None = Field(
        default=None,
        alias="-m",
        description=(
            "Suppress all alignments for a particular read or pair if more than <int> reportable alignments exist for it. "
            "Reportable alignments are those that would be reported given the -n, -v, -l, -e, -k, -a, --best, and --strata options. "
            "Default: no limit. Bowtie is designed to be very fast for small -m but bowtie can become significantly slower "
            "for larger values of -m. If you would like to use Bowtie for larger values of -k, consider building an index "
            "with a denser suffix-array sample, i.e. specify a smaller -o/--offrate when invoking bowtie-build for the relevant index "
            "(see the Performance tuning section for details)."
        ),
    )
    M: PositiveInt | None = Field(
        default=None,
        alias="-M",
        description=(
            "Behaves like -m except that if a read has more than <int> reportable alignments, one is reported at random. "
            "In default output mode, the selected alignment's 7th column is set to <int>+1 to indicate the read has at least "
            "<int>+1 valid alignments. In -S/--sam mode, the selected alignment is given a MAPQ (mapping quality) of 0 "
            "and the XM:I field is set to <int>+1. This option requires --best; if specified without --best, "
            "--best is enabled automatically."
        ),
    )
    best: Literal[""] | None = Field(
        default=None,
        alias="--best",
        description=(
            "Make Bowtie guarantee that reported singleton alignments are 'best' in terms of stratum "
            "(i.e. number of mismatches, or mismatches in the seed in the case of -n mode) and in terms of "
            "the quality values at the mismatched position(s). Stratum always trumps quality; e.g. a 1-mismatch alignment "
            "where the mismatched position has Phred quality 40 is preferred over a 2-mismatch alignment where the mismatched "
            "positions both have Phred quality 10. When --best is not specified, Bowtie may report alignments that are "
            "sub-optimal in terms of stratum and/or quality (though an effort is made to report the best alignment). "
            "--best mode also removes all strand bias. Note that --best does not affect which alignments are considered "
            "'valid' by bowtie, only which valid alignments are reported by bowtie. When --best is specified and multiple "
            "hits are allowed (via -k or -a), the alignments for a given read are guaranteed to appear in best-to-worst "
            "order in bowtie’s output. bowtie is somewhat slower when --best is specified."
        ),
    )
    strata: Literal[""] | None = Field(
        default=None,
        alias="--strata",
        description=(
            "If many valid alignments exist and are reportable (e.g. are not disallowed via the -k option) "
            "and they fall into more than one alignment 'stratum', report only those alignments that fall "
            "into the best stratum. By default, Bowtie reports all reportable alignments regardless of "
            "whether they fall into multiple strata. When --strata is specified, --best must also be specified."
        ),
    )
    # don't allow any output parameters
    # don't allow any SAM parameters
    offrate: PositiveInt | None = Field(
        default=None,
        alias="-o",
        validation_alias=AliasChoices("-o", "--offrate"),
        description=(
            "Override the offrate of the index with <int>. If <int> is greater than the offrate used to build the index, "
            "then some row markings are discarded when the index is read into memory. This reduces the memory footprint "
            "of the aligner but requires more time to calculate text offsets. <int> must be greater than the value used "
            "to build the index."
        ),
    )
    threads: PositiveInt | None = Field(
        default=None,
        alias="-p",
        validation_alias=AliasChoices("-p", "--threads"),
        description=(
            "Launch <int> parallel search threads (default: 1). Threads will run on separate processors/cores and "
            "synchronize when parsing reads and outputting alignments. Searching for alignments is highly parallel, "
            "and speedup is fairly close to linear."
        ),
    )
    reorder: Literal[""] | None = Field(
        default=None,
        alias="--reorder",
        description=(
            "Guarantees that output SAM records are printed in an order corresponding to the order of the reads in the "
            "original input file, even when -p is set greater than 1. Specifying --reorder and setting -p greater than 1 "
            "causes Bowtie to run somewhat slower than if --reorder were not specified. Has no effect if -p is set to 1, "
            "since output order will naturally correspond to input order in that case. It is an error to specify "
            "--reorder without the -S parameter. N.B. --reorder does not affect the outputs of --al/--max/--un."
        ),
    )
    mm: Literal[""] | None = Field(
        default=None,
        alias="--mm",
        description=(
            "Use memory-mapped I/O to load the index, rather than normal C file I/O. Memory-mapping the index allows "
            "many concurrent bowtie processes on the same computer to share the same memory image of the index "
            "(i.e. you pay the memory overhead just once). This facilitates memory-efficient parallelization of bowtie "
            "in situations where using -p is not possible."
        ),
    )
    shmem: Literal[""] | None = Field(
        default=None,
        alias="--shmem",
        description=(
            "Use shared memory to load the index, rather than normal C file I/O. Using shared memory allows many "
            "concurrent bowtie processes on the same computer to share the same memory image of the index "
            "(i.e. you pay the memory overhead just once). This facilitates memory-efficient parallelization of bowtie "
            "in situations where using -p is not desirable. Unlike --mm, --shmem installs the index into shared memory "
            "permanently, or until the user deletes the shared memory chunks manually. See your operating system "
            "documentation for details on how to manually list and remove shared memory chunks (on Linux and Mac OS X, "
            "these commands are ipcs and ipcrm). You may also need to increase your OS's maximum shared-memory chunk "
            "size to accommodate larger indexes; see your OS documentation."
        ),
    )


class CrossHybridizationProbabilityFilterBlastnConfig(BaseModel):
    model_config = ConfigDict(extra="forbid")

    alignment_method: Annotated[
        Literal["blastn"],
        Field(description="Which alignment method to use; options are 'blastn' or 'bowtie'."),
    ]
    search_parameters: BlastnSearchParameters
    hit_parameters: BlastnHitParameters


class CrossHybridizationProbabilityFilterBowtieConfig(BaseModel):
    model_config = ConfigDict(extra="forbid")

    alignment_method: Annotated[
        Literal["bowtie"],
        Field(description="Which alignment method to use; options are 'blastn' or 'bowtie'."),
    ]
    search_parameters: BowtieSearchParameters


CrossHybridizationProbabilityFilterConfig = Annotated[
    CrossHybridizationProbabilityFilterBlastnConfig | CrossHybridizationProbabilityFilterBowtieConfig,
    Field(discriminator="alignment_method"),
]


class HybridizationProbabilityFilterBlastnConfig(CrossHybridizationProbabilityFilterBlastnConfig):
    threshold: Annotated[
        FractionT,
        Field(
            description="Threshold for hybridization probability filtering. Probes with hybridization probabilities above this threshold are removed, as they may bind non-specifically."
        ),
    ]


class HybridizationProbabilityFilterBowtieConfig(CrossHybridizationProbabilityFilterBowtieConfig):
    threshold: Annotated[
        FractionT,
        Field(
            description="Threshold for hybridization probability filtering. Probes with hybridization probabilities above this threshold are removed, as they may bind non-specifically."
        ),
    ]


HybridizationProbabilityFilterConfig = Annotated[
    HybridizationProbabilityFilterBlastnConfig | HybridizationProbabilityFilterBowtieConfig,
    Field(discriminator="alignment_method"),
]


class GenomicRegions(BaseModel):
    model_config = ConfigDict(extra="forbid")

    gene: Annotated[bool, Field(description="Generate gene regions.")]
    intergenic: Annotated[bool, Field(description="Generate intergenic regions.")]
    exon: Annotated[bool, Field(description="Generate exon regions.")]
    exon_exon_junction: Annotated[bool, Field(description="Generate exon–exon junction regions.")]
    utr: Annotated[bool, Field(description="Generate UTR regions.")]
    cds: Annotated[bool, Field(description="Generate coding sequence (CDS) regions.")]
    intron: Annotated[bool, Field(description="Generate intron regions.")]


class SourceParamsCustom(BaseModel):
    model_config = ConfigDict(extra="forbid")

    file_annotation: Annotated[
        str,
        Field(
            description="GTF file with gene annotation",
        ),
    ]
    file_sequence: Annotated[
        str,
        Field(
            description="FASTA file with genome sequence",
        ),
    ]
    files_source: Annotated[str | None, Field(description="original source of the genomic files")]
    species: Annotated[
        str | None,
        Field(description="species of provided annotation, set to 'None' if unknown"),
    ]
    annotation_release: Annotated[
        str | None,
        Field(description="release number of provided annotation, set to 'None' if unknown"),
    ]
    genome_assembly: Annotated[
        str | None,
        Field(description="genome assembly of provided annotation, set to 'None' if unknown"),
    ]


class SourceParamsEnsembl(BaseModel):
    model_config = ConfigDict(extra="forbid")

    species: Annotated[str, Field(description="species of provided annotation")]
    annotation_release: Annotated[str, Field(description="release number of provided annotation")]


class SourceParamsNcbi(BaseModel):
    model_config = ConfigDict(extra="forbid")

    taxon: Annotated[
        Literal[
            "archaea",
            "bacteria",
            "fungi",
            "invertebrate",
            "metagenomes",
            "plant",
            "plants",
            "protozoa",
            "unkown",
            "vertebrate_mammalian",
            "vertebrate_other",
            "viral",
        ],
        Field(description="taxon of the species"),
    ]
    species: Annotated[str, Field(description="species of provided annotation")]
    annotation_release: Annotated[str, Field(description="release number of provided annotation")]


class SourceCustom(BaseModel):
    model_config = ConfigDict(extra="forbid")

    source: Annotated[
        Literal["custom"],
        Field(
            description="Indicate which annotation should be used. Possible values are 'ensembl', 'ncbi' or 'custom'."
        ),
    ]
    parameters: Annotated[
        SourceParamsCustom,
        Field(
            description="If a custom source, the metadata of the provided genome and annotation. If ensembl or ncbi, the parameters used to retrieve the genomic data."
        ),
    ]


class SourceEnsembl(BaseModel):
    model_config = ConfigDict(extra="forbid")

    source: Annotated[
        Literal["ensembl"],
        Field(
            description="Indicate which annotation should be used. Possible values are 'ensembl', 'ncbi' or 'custom'."
        ),
    ]
    parameters: Annotated[
        SourceParamsEnsembl,
        Field(
            description="If a custom source, the metadata of the provided genome and annotation. If ensembl or ncbi, the parameters used to retrieve the genomic data."
        ),
    ]


class SourceNcbi(BaseModel):
    model_config = ConfigDict(extra="forbid")

    source: Annotated[
        Literal["ncbi"],
        Field(
            description="Indicate which annotation should be used. Possible values are 'ensembl', 'ncbi' or 'custom'."
        ),
    ]
    parameters: Annotated[
        SourceParamsNcbi,
        Field(
            description="If a custom source, the metadata of the provided genome and annotation. If ensembl or ncbi, the parameters used to retrieve the genomic data."
        ),
    ]


SourceConfigs = Annotated[SourceCustom | SourceEnsembl | SourceNcbi, Field(discriminator="source")]
