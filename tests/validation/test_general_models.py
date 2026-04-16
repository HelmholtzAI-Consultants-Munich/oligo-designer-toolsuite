"""Tests for general-purpose models in _general.py."""

from typing import Any

import pytest
from pydantic import TypeAdapter, ValidationError

from oligo_designer_toolsuite.validation.models._general import (
    BaseProbabilities,
    BlastnHitParameters,
    BlastnSearchParameters,
    BowtieSearchParameters,
    CrossHybridizationProbabilityFilterBlastnConfig,
    CrossHybridizationProbabilityFilterBowtieConfig,
    CrossHybridizationProbabilityFilterConfig,
    General,
    HomopolymerThresholds,
    HybridizationProbabilityFilterBlastnConfig,
    HybridizationProbabilityFilterBowtieConfig,
    HybridizationProbabilityFilterConfig,
    OligoPropertyWeights,
    OligoSetSelection,
)

_CrossHybAdapter: TypeAdapter[CrossHybridizationProbabilityFilterConfig] = TypeAdapter(
    CrossHybridizationProbabilityFilterConfig
)
_HybAdapter: TypeAdapter[HybridizationProbabilityFilterConfig] = TypeAdapter(
    HybridizationProbabilityFilterConfig
)


# ---------------------------------------------------------------------------
# HomopolymerThresholds
# ---------------------------------------------------------------------------


class TestHomopolymerThresholds:
    def test_invalid_zero(self) -> None:
        # PositiveInt boundary: 0 is not allowed
        with pytest.raises(ValidationError):
            HomopolymerThresholds(A=0)

    def test_extra_field_forbidden(self) -> None:
        with pytest.raises(ValidationError):
            HomopolymerThresholds(A=3, X=1)  # type: ignore[call-arg]


# ---------------------------------------------------------------------------
# General
# ---------------------------------------------------------------------------


class TestGeneral:
    def test_invalid_n_jobs_zero(self) -> None:
        # PositiveInt boundary: 0 is not allowed
        with pytest.raises(ValidationError):
            General(n_jobs=0, dir_output="/out")

    def test_extra_field_forbidden(self) -> None:
        with pytest.raises(ValidationError):
            General(n_jobs=1, dir_output="/out", unknown_field=42)  # type: ignore[call-arg]


# ---------------------------------------------------------------------------
# OligoSetSelection
# ---------------------------------------------------------------------------


class TestOligoSetSelection:
    def _make(self, **overrides: Any) -> OligoSetSelection:
        fields = dict(
            n_attempts_graph=1,
            n_attempts_clique_enum=1,
            diversification_fraction=0.3,
            jaccard_opt=0.5,
            jaccard_step=0.1,
        )
        fields.update(overrides)
        return OligoSetSelection(**fields)

    def test_boundary_fractions_zero(self) -> None:
        s = self._make(diversification_fraction=0.0, jaccard_opt=0.0, jaccard_step=0.0)
        assert s.diversification_fraction == 0.0

    def test_boundary_fractions_one(self) -> None:
        s = self._make(diversification_fraction=1.0, jaccard_opt=1.0, jaccard_step=1.0)
        assert s.jaccard_opt == 1.0

    def test_invalid_n_attempts_graph_zero(self) -> None:
        with pytest.raises(ValidationError):
            self._make(n_attempts_graph=0)

    def test_invalid_diversification_fraction_above_one(self) -> None:
        with pytest.raises(ValidationError):
            self._make(diversification_fraction=1.1)

    def test_invalid_jaccard_opt_negative(self) -> None:
        with pytest.raises(ValidationError):
            self._make(jaccard_opt=-0.1)

    def test_extra_field_forbidden(self) -> None:
        with pytest.raises(ValidationError):
            self._make(extra=99)


# ---------------------------------------------------------------------------
# BaseProbabilities — custom validator: _check_sums_up_to_1
# ---------------------------------------------------------------------------


class TestBaseProbabilities:
    def test_valid_sum_up_to_one(self) -> None:
        # test that the validator that calculates the sum doesn't throw an error
        bp = BaseProbabilities(A=0.25, T=0.25, C=0.25, G=0.25)
        assert bp.A == 0.25

    def test_valid_one_base_all(self) -> None:
        # Edge case: all probability on one base still sums to 1
        bp = BaseProbabilities(A=1.0, T=0.0, C=0.0, G=0.0)
        assert bp.A == 1.0

    def test_invalid_does_not_sum_to_one(self) -> None:
        with pytest.raises(ValidationError):
            BaseProbabilities(A=0.3, T=0.3, C=0.3, G=0.3)

    def test_invalid_sum_less_than_one(self) -> None:
        with pytest.raises(ValidationError):
            BaseProbabilities(A=0.1, T=0.1, C=0.1, G=0.1)

    def test_extra_field_forbidden(self) -> None:
        with pytest.raises(ValidationError):
            BaseProbabilities(A=0.25, T=0.25, C=0.25, G=0.25, X=0.0)  # type: ignore[call-arg]


# ---------------------------------------------------------------------------
# BlastnSearchParameters
# ---------------------------------------------------------------------------


class TestBlastnSearchParameters:
    def test_extra_field_forbidden(self) -> None:
        with pytest.raises(ValidationError):
            BlastnSearchParameters(unknown=123)  # type: ignore[call-arg]

    @pytest.mark.parametrize("field", ["-query", "-db", "-out"])  # type: ignore[untyped-decorator]
    def test_forbids_unsupported_blastn_fields(self, field: str) -> None:
        with pytest.raises(ValidationError):
            BlastnSearchParameters(**{field: "x"})

    def test_accepts_alias_input(self) -> None:
        model = BlastnSearchParameters(**{"-strand": "plus"})
        assert model.strand == "plus"

    def test_accepts_field_name_input(self) -> None:
        model = BlastnSearchParameters(strand="plus")
        assert model.strand == "plus"

    def test_all_fields_none_by_default(self) -> None:
        model = BlastnSearchParameters()
        for field_name in BlastnSearchParameters.model_fields:
            assert getattr(model, field_name) is None, f"{field_name} should default to None"

    @pytest.mark.parametrize(  # type: ignore[untyped-decorator]
        ("payload", "attr", "expected"),
        [
            ({"-strand": "plus"}, "strand", "plus"),
            ({"-task": "blastn"}, "task", "blastn"),
            ({"-template_type": "optimal"}, "template_type", "optimal"),
            ({"-template_length": "18"}, "template_length", "18"),
        ],
    )
    def test_accepts_valid_literal_values(self, payload: dict[str, str], attr: str, expected: str) -> None:
        model = BlastnSearchParameters(**payload)
        assert getattr(model, attr) == expected

    @pytest.mark.parametrize(  # type: ignore[untyped-decorator]
        "payload",
        [
            {"-strand": "forward"},
            {"-task": "foo"},
            {"-template_type": "bad"},
            {"-template_length": "20"},
        ],
    )
    def test_rejects_invalid_literal_values(self, payload: dict[str, str]) -> None:
        with pytest.raises(ValidationError):
            BlastnSearchParameters(**payload)

    @pytest.mark.parametrize("field", ["gapopen", "num_alignments", "window_size"])  # type: ignore[untyped-decorator]
    def test_non_negative_int_fields_accept_zero(self, field: str) -> None:
        model = BlastnSearchParameters(**{field: 0})
        assert getattr(model, field) == 0

    @pytest.mark.parametrize("field", ["gapopen", "num_alignments", "window_size"])  # type: ignore[untyped-decorator]
    def test_non_negative_int_fields_reject_negative_values(self, field: str) -> None:
        with pytest.raises(ValidationError):
            BlastnSearchParameters(**{field: -1})

    @pytest.mark.parametrize("field", ["max_hsps", "max_target_seqs"])  # type: ignore[untyped-decorator]
    def test_positive_int_fields_accept_positive_values(self, field: str) -> None:
        model = BlastnSearchParameters(**{field: 1})
        assert getattr(model, field) == 1

    @pytest.mark.parametrize("field", ["max_hsps", "max_target_seqs"])  # type: ignore[untyped-decorator]
    def test_positive_int_fields_reject_zero(self, field: str) -> None:
        with pytest.raises(ValidationError):
            BlastnSearchParameters(**{field: 0})

    @pytest.mark.parametrize("value", [0, 4])  # type: ignore[untyped-decorator]
    def test_sorthits_accepts_bounds(self, value: int) -> None:
        model = BlastnSearchParameters(**{"-sorthits": value})
        assert model.sorthits == value

    @pytest.mark.parametrize("value", [-1, 5])  # type: ignore[untyped-decorator]
    def test_sorthits_rejects_out_of_range(self, value: int) -> None:
        with pytest.raises(ValidationError):
            BlastnSearchParameters(**{"-sorthits": value})

    @pytest.mark.parametrize(  # type: ignore[untyped-decorator]
        "field,value",
        [("perc_identity", 0), ("perc_identity", 100), ("qcov_hsp_perc", 0), ("qcov_hsp_perc", 100)],
    )
    def test_percentage_fields_accept_bounds(self, field: str, value: int) -> None:
        model = BlastnSearchParameters(**{field: value})
        assert getattr(model, field) == value

    @pytest.mark.parametrize(  # type: ignore[untyped-decorator]
        "field,value",
        [("perc_identity", -0.1), ("perc_identity", 100.1), ("qcov_hsp_perc", -1), ("qcov_hsp_perc", 101)],
    )
    def test_percentage_fields_reject_out_of_range(self, field: str, value: int | float) -> None:
        with pytest.raises(ValidationError):
            BlastnSearchParameters(**{field: value})

    @pytest.mark.parametrize("field", ["best_hit_overhang", "best_hit_score_edge"])  # type: ignore[untyped-decorator]
    def test_best_hit_fields_accept_valid_inner_value(self, field: str) -> None:
        model = BlastnSearchParameters(**{field: 0.1})
        assert getattr(model, field) == 0.1

    @pytest.mark.parametrize(  # type: ignore[untyped-decorator]
        "field,value",
        [
            ("-best_hit_overhang", 0),
            ("-best_hit_overhang", 0.5),
            ("-best_hit_score_edge", 0),
            ("-best_hit_score_edge", 0.5),
        ],
    )
    def test_best_hit_fields_reject_boundaries(self, field: str, value: int | float) -> None:
        with pytest.raises(ValidationError):
            BlastnSearchParameters(**{field: value})

    @pytest.mark.parametrize("field", ["lcase_masking", "subject_besthit", "no_greedy", "ungapped"])  # type: ignore[untyped-decorator]
    def test_empty_string_literal_flags_accept_empty_string(self, field: str) -> None:
        model = BlastnSearchParameters(**{field: ""})
        assert getattr(model, field) == ""

    @pytest.mark.parametrize("field", ["-lcase_masking", "-subject_besthit", "-no_greedy", "-ungapped"])  # type: ignore[untyped-decorator]
    def test_empty_string_literal_flags_reject_other_values(self, field: str) -> None:
        with pytest.raises(ValidationError):
            BlastnSearchParameters(**{field: "true"})


# ---------------------------------------------------------------------------
# BowtieSearchParameters
# ---------------------------------------------------------------------------


class TestBowtieSearchParameters:
    def test_extra_field_forbidden(self) -> None:
        with pytest.raises(ValidationError):
            BowtieSearchParameters(unknown=123)  # type: ignore[call-arg]

    def test_all_fields_none_by_default(self) -> None:
        model = BowtieSearchParameters()
        for field_name in BowtieSearchParameters.model_fields:
            assert getattr(model, field_name) is None, f"{field_name} should default to None"

    def test_accepts_alias_input(self) -> None:
        model = BowtieSearchParameters(**{"-n": 2})
        assert model.seedmms == 2

    def test_accepts_field_name_input(self) -> None:
        model = BowtieSearchParameters(seedmms=2)
        assert model.seedmms == 2

    @pytest.mark.parametrize("value", [0, 1, 2, 3])  # type: ignore[untyped-decorator]
    def test_seedmms_accepts_valid_range(self, value: int) -> None:
        model = BowtieSearchParameters(seedmms=value)
        assert model.seedmms == value

    def test_seedmms_rejects_above_3(self) -> None:
        with pytest.raises(ValidationError):
            BowtieSearchParameters(seedmms=4)

    def test_seedmms_rejects_negative(self) -> None:
        with pytest.raises(ValidationError):
            BowtieSearchParameters(seedmms=-1)

    def test_seedlen_rejects_below_5(self) -> None:
        with pytest.raises(ValidationError):
            BowtieSearchParameters(seedlen=4)

    def test_maqerr_rejects_negative(self) -> None:
        with pytest.raises(ValidationError):
            BowtieSearchParameters(maqerr=-1)

    @pytest.mark.parametrize(  # type: ignore[untyped-decorator]
        "alias,field",
        [
            ("--maxbts", "maxbts"),
            ("--chunkmbs", "chunkmbs"),
            ("--reads-per-batch", "reads_per_batch"),
            ("-k", "k"),
            ("-m", "m"),
            ("-M", "M"),
            ("-p", "threads"),
        ],
    )
    def test_positive_int_fields_accept_positive_values(self, alias: str, field: str) -> None:
        model = BowtieSearchParameters(**{alias: 1})
        assert getattr(model, field) == 1

    @pytest.mark.parametrize("alias", ["--maxbts", "--chunkmbs", "--reads-per-batch", "-k", "-m", "-M", "-p"])  # type: ignore[untyped-decorator]
    def test_positive_int_fields_reject_zero(self, alias: str) -> None:
        with pytest.raises(ValidationError):
            BowtieSearchParameters(**{alias: 0})

    @pytest.mark.parametrize(  # type: ignore[untyped-decorator]
        "alias,field",
        [
            ("--nomaqround", "nomaqround"),
            ("--nofw", "nofw"),
            ("--norc", "norc"),
            ("-y", "tryhard"),
            ("--best", "best"),
            ("--strata", "strata"),
            ("--reorder", "reorder"),
            ("--mm", "mm"),
            ("--shmem", "shmem"),
        ],
    )
    def test_empty_string_literal_flags_accept_empty_string(self, alias: str, field: str) -> None:
        model = BowtieSearchParameters(**{alias: ""})
        assert getattr(model, field) == ""

    @pytest.mark.parametrize(  # type: ignore[untyped-decorator]
        "alias",
        ["--nomaqround", "--nofw", "--norc", "-y", "--best", "--strata", "--reorder", "--mm", "--shmem"],
    )
    def test_empty_string_literal_flags_reject_other_values(self, alias: str) -> None:
        with pytest.raises(ValidationError):
            BowtieSearchParameters(**{alias: "true"})

    @pytest.mark.parametrize("field", ["--offbase", "-B", "--max"])  # type: ignore[untyped-decorator]
    def test_forbids_unsupported_bowtie_fields(self, field: str) -> None:
        with pytest.raises(ValidationError):
            BowtieSearchParameters(**{field: "x"})


# ---------------------------------------------------------------------------
# CrossHybridizationProbabilityFilterConfig — discriminator="alignment_method"
# ---------------------------------------------------------------------------


class TestCrossHybridizationProbabilityFilterConfig:
    def test_blastn_mode(self) -> None:
        obj = _CrossHybAdapter.validate_python(
            {
                "alignment_method": "blastn",
                "search_parameters": {},
                "hit_parameters": {"coverage": 80.0},
            }
        )
        assert isinstance(obj, CrossHybridizationProbabilityFilterBlastnConfig)

    def test_bowtie_mode(self) -> None:
        obj = _CrossHybAdapter.validate_python(
            {
                "alignment_method": "bowtie",
                "search_parameters": {},
            }
        )
        assert isinstance(obj, CrossHybridizationProbabilityFilterBowtieConfig)

    def test_invalid_mode(self) -> None:
        with pytest.raises(ValidationError):
            _CrossHybAdapter.validate_python({"alignment_method": "bwa"})

    def test_blastn_requires_hit_parameters(self) -> None:
        with pytest.raises(ValidationError):
            _CrossHybAdapter.validate_python(
                {
                    "alignment_method": "blastn",
                    "search_parameters": {},
                }
            )

    def test_bowtie_rejects_hit_parameters(self) -> None:
        with pytest.raises(ValidationError):
            _CrossHybAdapter.validate_python(
                {
                    "alignment_method": "bowtie",
                    "search_parameters": {},
                    "hit_parameters": {"coverage": 80.0},
                }
            )


# ---------------------------------------------------------------------------
# HybridizationProbabilityFilterConfig — discriminator="alignment_method"
# ---------------------------------------------------------------------------


class TestHybridizationProbabilityFilterConfig:
    def test_blastn_mode(self) -> None:
        obj = _HybAdapter.validate_python(
            {
                "alignment_method": "blastn",
                "search_parameters": {},
                "hit_parameters": {"coverage": 80.0},
                "threshold": 0.5,
            }
        )
        assert isinstance(obj, HybridizationProbabilityFilterBlastnConfig)
        assert obj.threshold == 0.5

    def test_bowtie_mode(self) -> None:
        obj = _HybAdapter.validate_python(
            {
                "alignment_method": "bowtie",
                "search_parameters": {},
                "threshold": 0.5,
            }
        )
        assert isinstance(obj, HybridizationProbabilityFilterBowtieConfig)
        assert obj.threshold == 0.5

    def test_invalid_mode(self) -> None:
        with pytest.raises(ValidationError):
            _HybAdapter.validate_python({"alignment_method": "bwa", "threshold": 0.5})

    def test_threshold_required(self) -> None:
        with pytest.raises(ValidationError):
            _HybAdapter.validate_python(
                {
                    "alignment_method": "bowtie",
                    "search_parameters": {},
                }
            )


# ---------------------------------------------------------------------------
# BlastnHitParameters — custom validator: _check_mutually_exclusive
# ---------------------------------------------------------------------------


class TestBlastnHitParameters:
    def test_valid_coverage_only(self) -> None:
        p = BlastnHitParameters(coverage=80.0)
        assert p.coverage == 80.0
        assert p.min_alignment_length is None

    def test_valid_min_alignment_length_only(self) -> None:
        p = BlastnHitParameters(min_alignment_length=20)
        assert p.min_alignment_length == 20
        assert p.coverage is None

    def test_invalid_both_set(self) -> None:
        with pytest.raises(ValidationError):
            BlastnHitParameters(coverage=80.0, min_alignment_length=20)

    def test_invalid_neither_set(self) -> None:
        with pytest.raises(ValidationError):
            BlastnHitParameters()

    def test_invalid_coverage_above_100(self) -> None:
        # boundary: coverage is [0, 100]
        with pytest.raises(ValidationError):
            BlastnHitParameters(coverage=101.0)

    def test_invalid_coverage_negative(self) -> None:
        # boundary: coverage is [0, 100]
        with pytest.raises(ValidationError):
            BlastnHitParameters(coverage=-1.0)

    def test_invalid_min_alignment_length_negative(self) -> None:
        # boundary: NonNegativeInt
        with pytest.raises(ValidationError):
            BlastnHitParameters(min_alignment_length=-1)

    def test_extra_field_forbidden(self) -> None:
        with pytest.raises(ValidationError):
            BlastnHitParameters(coverage=50.0, extra_field=1)  # type: ignore[call-arg]


# ---------------------------------------------------------------------------
# OligoPropertyWeights
# ---------------------------------------------------------------------------


class TestOligoPropertyWeights:
    def test_extra_field_forbidden(self) -> None:
        with pytest.raises(ValidationError):
            OligoPropertyWeights(unknown_weight=1.0)  # type: ignore[call-arg]
