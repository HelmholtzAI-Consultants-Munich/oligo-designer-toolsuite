from typing import Annotated, Literal

from pydantic import (
    AliasChoices,
    BaseModel,
    ConfigDict,
    Field,
    NonNegativeInt,
    NonPositiveInt,
    PositiveInt,
    model_serializer,
    model_validator,
)
from typing_extensions import Self


class General(BaseModel):
    model_config = ConfigDict(extra="forbid")
    n_jobs: PositiveInt = Field(
        description="number of cores used to run the pipeline and 2*n_jobs +1 of regions that should be stored in cache. If memory consumption of pipeline is too high reduce this number, if a lot of RAM is available increase this number to decrease runtime"
    )
    dir_output: str = Field(description="name of the directory where the output files will be written")
    write_intermediate_steps: bool = Field(
        description="if true, writes the oligo sequences after each step of the pipeline into a csv file",
    )


class HomopolymerThresholds(BaseModel):
    model_config = ConfigDict(extra="forbid")
    A: PositiveInt | None = None
    T: PositiveInt | None = None
    C: PositiveInt | None = None
    G: PositiveInt | None = None


class BlastnSearchParameters(BaseModel):
    model_config = ConfigDict(extra="forbid", validate_by_name=True, validate_by_alias=True)

    # don't allow
    # -h
    # -help
    # -version
    # -query
    query_loc: str | None = Field(
        default=None,
        validation_alias=AliasChoices("-query_loc", "query_loc"),
        serialization_alias="-query_loc",
        description="Location on the query sequence in 1-based offsets (Format: start-stop).",
    )
    strand: Literal["plus", "minus", "both"] | None = Field(
        default=None,
        validation_alias=AliasChoices("-strand", "strand"),
        serialization_alias="-strand",
        description="Query strand(s) to search against database/subject. Choice of both, minus, or plus.",
    )
    task: Literal["megablast", "dc-megablast", "blastn", "blastn-short", "rmblastn"] | None = Field(
        default=None,
        validation_alias=AliasChoices("-task", "task"),
        serialization_alias="-task",
        description="Supported tasks.",
    )
    # don't allow
    # -db
    # -out
    evalue: float | None = Field(
        default=None,
        validation_alias=AliasChoices("-evalue", "evalue"),
        serialization_alias="-evalue",
        description="Expectation value (E) threshold for saving hits. Default = 10 (1000 for blastn-short)",
    )
    word_size: NonNegativeInt | None = Field(
        default=None,
        validation_alias=AliasChoices("-word_size", "word_size"),
        serialization_alias="-word_size",
        ge=4,
        description="Length of initial exact match.",
    )
    gapopen: NonNegativeInt | None = Field(
        default=None,
        validation_alias=AliasChoices("-gapopen", "gapopen"),
        serialization_alias="-gapopen",
        description="Cost to open a gap.",
    )
    gapextend: NonNegativeInt | None = Field(
        default=None,
        validation_alias=AliasChoices("-gapextend", "gapextend"),
        serialization_alias="-gapextend",
        description="Cost to extend a gap.",
    )
    penalty: NonPositiveInt | None = Field(
        default=None,
        validation_alias=AliasChoices("-penalty", "penalty"),
        serialization_alias="-penalty",
        description="Penalty for a nucleotide mismatch.",
    )
    reward: NonNegativeInt | None = Field(
        default=None,
        validation_alias=AliasChoices("-reward", "reward"),
        serialization_alias="-reward",
        description="Reward for a nucleotide match.",
    )
    # don't allow use_index/index_name as another file would be needed
    # don't allow subject/subject_loc as another file would be needed
    # don't allow
    # -outfmt
    # -show_gis
    num_descriptions: NonNegativeInt | None = Field(
        default=None,
        validation_alias=AliasChoices("-num_descriptions", "num_descriptions"),
        serialization_alias="-num_descriptions",
        description="Number of database sequences to show one-line descriptions for.",
    )
    num_alignments: NonNegativeInt | None = Field(
        default=None,
        validation_alias=AliasChoices("-num_alignments", "num_alignments"),
        serialization_alias="-num_alignments",
        description="Number of database sequences to show alignments for.",
    )
    # don't allow
    # line_length
    # html
    sorthits: NonNegativeInt | None = Field(
        default=None,
        validation_alias=AliasChoices("-sorthits", "sorthits"),
        serialization_alias="-sorthits",
        le=4,
        description="Sorting option for hits.",
    )
    sorthsps: NonNegativeInt | None = Field(
        default=None,
        validation_alias=AliasChoices("-sorthsps", "sorthsps"),
        serialization_alias="-sorthsps",
        le=4,
        description="Sorting option for hps.",
    )
    dust: str | None = Field(
        default=None,
        validation_alias=AliasChoices("-dust", "dust"),
        serialization_alias="-dust",
        description="Filter query sequence with dust.",
    )
    # don't allow
    # filtering_db as another file would be needed
    # window_masker_taxid
    # window_masker_db as then another file would be needed
    soft_masking: bool | None = Field(
        default=None,
        validation_alias=AliasChoices("-soft_masking", "soft_masking"),
        serialization_alias="-soft_masking",
        description="Apply filtering locations as soft masks (i.e., only for finding initial matches).",
    )
    lcase_masking: Literal[""] | None = Field(
        default=None,
        validation_alias=AliasChoices("-lcase_masking", "lcase_masking"),
        serialization_alias="-lcase_masking",
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
        validation_alias=AliasChoices("-db_soft_mask", "db_soft_mask"),
        serialization_alias="-db_soft_mask",
        description="Filtering algorithm ID to apply to the BLAST database as soft mask (i.e., only for finding initial matches).",
    )
    db_hard_mask: str | None = Field(
        default=None,
        validation_alias=AliasChoices("-db_hard_mask", "db_hard_mask"),
        serialization_alias="-db_hard_mask",
        description="Filtering algorithm ID to apply to the BLAST database as hard mask (i.e., sequence is masked for all phases of search).",
    )
    perc_identity: float | None = Field(
        default=None,
        validation_alias=AliasChoices("-perc_identity", "perc_identity"),
        serialization_alias="-perc_identity",
        description="Percent identity cutoff.",
        ge=0,
        le=100,
    )
    qcov_hsp_perc: float | None = Field(
        default=None,
        validation_alias=AliasChoices("-qcov_hsp_perc", "qcov_hsp_perc"),
        serialization_alias="-qcov_hsp_perc",
        description="Percent query coverage per hsp.",
        ge=0,
        le=100,
    )
    max_hsps: PositiveInt | None = Field(
        default=None,
        validation_alias=AliasChoices("-max_hsps", "max_hsps"),
        serialization_alias="-max_hsps",
        description="Set maximum number of HSPs per subject sequence to save for each query.",
    )
    culling_limit: NonNegativeInt | None = Field(
        default=None,
        validation_alias=AliasChoices("-culling_limit", "culling_limit"),
        serialization_alias="-culling_limit",
        description="If the query range of a hit is enveloped by that of at least this many higher-scoring hits, delete the hit",
    )
    best_hit_overhang: float | None = Field(
        default=None,
        validation_alias=AliasChoices("-best_hit_overhang", "best_hit_overhang"),
        serialization_alias="-best_hit_overhang",
        description="Best Hit algorithm overhang value (recommended value: 0.1).",
        gt=0,
        lt=0.5,
    )
    best_hit_score_edge: float | None = Field(
        default=None,
        validation_alias=AliasChoices("-best_hit_score_edge", "best_hit_score_edge"),
        serialization_alias="-best_hit_score_edge",
        description="Best Hit algorithm score edge value (recommended value: 0.1)",
        gt=0,
        lt=0.5,
    )
    subject_besthit: Literal[""] | None = Field(
        default=None,
        validation_alias=AliasChoices("-subject_besthit", "subject_besthit"),
        serialization_alias="-subject_besthit",
        description="Turn on best hit per subject sequence.",
    )
    max_target_seqs: PositiveInt | None = Field(
        default=None,
        validation_alias=AliasChoices("-max_target_seqs", "max_target_seqs"),
        serialization_alias="-max_target_seqs",
        description="Maximum number of aligned sequences to keep.",
    )
    template_type: Literal["coding", "coding_and_optimal", "optimal"] | None = Field(
        default=None,
        validation_alias=AliasChoices("-template_type", "template_type"),
        serialization_alias="-template_type",
        description="Discontiguous MegaBLAST template type. Allowed values are coding, optimal and coding_and_optimal.",
    )
    # template_length is actually int, but only 3 values, therefore implemented as literal
    template_length: Literal["16", "18", "21"] | None = Field(
        default=None,
        validation_alias=AliasChoices("-template_length", "template_length"),
        serialization_alias="-template_length",
        description="Discontiguous MegaBLAST template length.",
    )
    db_size: int | None = Field(
        default=None,
        validation_alias=AliasChoices("-db_size", "db_size"),
        serialization_alias="-db_size",
        description="Effective length of the database.",
    )
    searchsp: NonNegativeInt | None = Field(
        default=None,
        validation_alias=AliasChoices("-searchsp", "searchsp"),
        serialization_alias="-searchsp",
        description="Effective length of the search space.",
    )
    # don't allow because extra file needed
    # import_search_strategy
    # export_search_strategy
    xdrop_ungap: float | None = Field(
        default=None,
        validation_alias=AliasChoices("-xdrop_ungap", "xdrop_ungap"),
        serialization_alias="-xdrop_ungap",
        description="X-dropoff value (in bits) for ungapped extensions.",
    )
    xdrop_gap: float | None = Field(
        default=None,
        validation_alias=AliasChoices("-xdrop_gap", "xdrop_gap"),
        serialization_alias="-xdrop_gap",
        description="X-dropoff value (in bits) for preliminary gapped extensions.",
    )
    xdrop_gap_final: float | None = Field(
        default=None,
        validation_alias=AliasChoices("-xdrop_gap_final", "xdrop_gap_final"),
        serialization_alias="-xdrop_gap_final",
        description="X-dropoff value (in bits) for final gapped alignment.",
    )
    no_greedy: Literal[""] | None = Field(
        default=None,
        validation_alias=AliasChoices("-no_greedy", "no_greedy"),
        serialization_alias="-no_greedy",
        description="Use non-greedy dynamic programming extension.",
    )
    min_raw_gapped_score: int | None = Field(
        default=None,
        validation_alias=AliasChoices("-min_raw_gapped_score", "min_raw_gapped_score"),
        serialization_alias="-min_raw_gapped_score",
        description="Minimum raw gapped score to keep an alignment in the preliminary gapped and trace-back stages. Normally set based upon expect value.",
    )
    ungapped: Literal[""] | None = Field(
        default=None,
        validation_alias=AliasChoices("-ungapped", "ungapped"),
        serialization_alias="-ungapped",
        description="Perform ungapped alignment only?",
    )
    window_size: NonNegativeInt | None = Field(
        default=None,
        validation_alias=AliasChoices("-window_size", "window_size"),
        serialization_alias="-window_size",
        description="Multiple hits window size, use 0 to specify 1-hit algorithm.",
    )
    off_diagonal_range: NonNegativeInt | None = Field(
        default=None,
        validation_alias=AliasChoices("-off_diagonal_range", "off_diagonal_range"),
        serialization_alias="-off_diagonal_range",
        description="Number of off-diagonals to search for the 2nd hit, use 0 to turn off.",
    )
    # don't allow
    # parse_deflines
    # num_threads
    # mt_mode
    # remote

    @model_serializer
    def serialize(self) -> dict:
        return {
            self.__class__.model_fields[name].serialization_alias or name: value
            for name, value in self.__dict__.items()
            if value is not None
        }


class BlastnHitParameters(BaseModel):
    model_config = ConfigDict(extra="forbid")

    coverage: float | None = Field(
        default=None,
        ge=0,
        le=100,
        description="Coverage in %, alternatively, min_alignment_length can be used",
    )
    min_alignment_length: NonNegativeInt | None = Field(
        default=None,
        description="Number of nucleotides for alignment, alternatively, coverage can be used",
    )

    @model_validator(mode="after")
    def _check_mutually_exclusive(self) -> Self:
        if (self.coverage is None) == (self.min_alignment_length is None):
            raise ValueError("Exactly one of 'coverage' or 'min_alignment_length' must be set.")
        return self


class TmParameters(BaseModel):
    model_config = ConfigDict(extra="forbid")

    check: bool = Field(default=True, description="Checks if the sequence is valid for the given method.")
    strict: bool = Field(
        default=True,
        description="Do not allow base characters or neighbor duplex keys (e.g. 'AT/NA') that could not or not unambiguously be evaluated for the respective method.",
    )
    c_seq: str | None = Field(
        default=None,
        description="Complementary sequence. The sequence of the template/target in 3'->5' direction. c_seq is necessary for mismatch correction and dangling-ends correction. Both corrections will automatically be applied if mismatches or dangling ends are present.",
    )
    shift: int = Field(
        default=0,
        description="Shift of the primer/probe sequence on the template/target sequence. The shift parameter is necessary to align seq and c_seq if they have different lengths or if they should have dangling ends.",
    )
    selfcomp: bool = Field(
        default=False,
        description="Is the sequence self-complementary? If 'True' the primer is thought binding to itself, thus dnac2 is not considered.",
    )
    nn_table: Literal["DNA_NN1", "DNA_NN2", "DNA_NN3", "DNA_NN4"] | None = Field(
        default=None, description="Thermodynamic NN values."
    )
    tmm_table: Literal["DNA_TMM1"] | None = Field(
        default=None, description="Thermodynamic values for terminal mismatches."
    )
    imm_table: Literal["DNA_IMM1"] | None = Field(
        default=None,
        description="Thermodynamic values for internal mismatches, may include insosine mismatches.",
    )
    de_table: Literal["DNA_DE1"] | None = Field(
        default=None, description="Thermodynamic values for dangling ends."
    )
    dnac1: NonNegativeInt = Field(
        default=25, description="Concentration of the higher concentrated strand [nM]."
    )
    dnac2: NonNegativeInt = Field(
        default=25, description="Concentration of the lower concentrated strand [nM]."
    )
    saltcorr: NonNegativeInt = Field(
        default=5,
        ge=0,
        le=7,
        description="Salt correction method, see Bio.SeqUtils.MeltingTemp.salt_correction.",
    )
    Na: NonNegativeInt = Field(default=50, description="Concentration of the ions [mM].")
    K: NonNegativeInt = Field(default=0, description="Concentration of the ions [mM].")
    Tris: NonNegativeInt = Field(default=0, description="Concentration of the ions [mM].")
    Mg: NonNegativeInt = Field(default=0, description="Concentration of the ions [mM].")
    dNTPs: NonNegativeInt = Field(default=0, description="Concentration of the ions [mM].")


class TmChemCorrectionParametersDetails(BaseModel):
    model_config = ConfigDict(extra="forbid")

    # defaults are from Bio.SeqUtils.MeltingTemp.chem_correction
    DMSO: float = Field(default=0, ge=0, le=100, description="Percent DMSO")
    DMSOfactor: float = Field(default=0.75, description="How much Tm should decrease per percent DMSO")
    fmd: float = Field(
        default=0, description="Formamide concentration in %(fmdmethod=1) or molar (fmdmethod=2)."
    )
    fmdfactor: float = Field(default=0.65, description="How much Tm should decrease per percent formamide")
    fmdmethod: int = Field(
        default=1,
        ge=1,
        le=2,
        description="Tm = Tm - factor(%formamide) (Default); Tm = Tm + (0.453(f(GC)) - 2.88) x [formamide]",
    )
    GC: float | None = Field(default=None, ge=0, le=100, description="GC content in percent.")

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


class TmChemCorrectionParametersEnabled(BaseModel):
    model_config = ConfigDict(extra="forbid")

    enabled: Literal[True] = Field(
        description="Should chem correction be used in the calculation of the melting temperature?",
        default=True,
    )
    parameters: TmChemCorrectionParametersDetails = TmChemCorrectionParametersDetails()


class TmChemCorrectionParametersDisabled(BaseModel):
    model_config = ConfigDict(extra="ignore")

    enabled: Literal[False] = Field(
        description="Should chem correction be used in the calculation of the melting temperature?",
        default=False,
    )


TmChemCorrectionParameters = Annotated[
    TmChemCorrectionParametersEnabled | TmChemCorrectionParametersDisabled, Field(discriminator="enabled")
]


class TmSaltCorrectionParametersDetails(BaseModel):
    model_config = ConfigDict(extra="forbid")

    Na: NonNegativeInt = Field(default=0, description="[mM] of ion")
    K: NonNegativeInt = Field(default=0, description="[mM] of ion")
    Tris: NonNegativeInt = Field(default=0, description="[mM] of ion")
    Mg: NonNegativeInt = Field(default=0, description="[mM] of ion")
    dNTPs: NonNegativeInt = Field(default=0, description="[mM] of ion")
    method: PositiveInt = Field(
        default=1,
        ge=1,
        le=7,
        description="Correction method to be applied. Methods 1-4 correct Tm, method 5 corrects deltaS, methods 6 and 7 correct 1/Tm.",
    )


class TmSaltCorrectionParametersEnabled(BaseModel):
    model_config = ConfigDict(extra="forbid")

    enabled: Literal[True] = Field(
        description="Should salt correction be used in the calculation of the melting temperature?",
        default=True,
    )
    parameters: TmSaltCorrectionParametersDetails = TmSaltCorrectionParametersDetails()


class TmSaltCorrectionParametersDisabled(BaseModel):
    model_config = ConfigDict(extra="ignore")

    enabled: Literal[False] = Field(
        description="Should salt correction be used in the calculation of the melting temperature?",
        default=False,
    )


TmSaltCorrectionParameters = Annotated[
    TmSaltCorrectionParametersEnabled | TmSaltCorrectionParametersDisabled, Field(discriminator="enabled")
]
