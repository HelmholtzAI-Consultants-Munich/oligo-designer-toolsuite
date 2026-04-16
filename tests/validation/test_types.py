"""Tests for custom Annotated type constraints defined in _types.py."""

import pytest
from pydantic import BaseModel, ValidationError

from oligo_designer_toolsuite.validation._types import (
    DNAT,
    ExonExonJunctionBlockSizeT,
    FastaFileListT,
    FractionT,
    GCContentMaxT,
    GCContentMinT,
    GCContentOptT,
    LengthMaxT,
    LengthMinT,
    TmMaxT,
    TmMinT,
    TmOptT,
    TSecondaryStructureT,
)

# ---------------------------------------------------------------------------
# Wrapper models — required to exercise Annotated field constraints
# ---------------------------------------------------------------------------


class _DNA(BaseModel):
    seq: DNAT


class _GCMin(BaseModel):
    val: GCContentMinT


class _GCMax(BaseModel):
    val: GCContentMaxT


class _GCOpt(BaseModel):
    val: GCContentOptT


class _TmMin(BaseModel):
    val: TmMinT


class _TmOpt(BaseModel):
    val: TmOptT


class _TmMax(BaseModel):
    val: TmMaxT


class _TSecondary(BaseModel):
    val: TSecondaryStructureT


class _Fraction(BaseModel):
    val: FractionT


class _LengthMin(BaseModel):
    val: LengthMinT


class _LengthMax(BaseModel):
    val: LengthMaxT


class _FastaList(BaseModel):
    files: FastaFileListT


class _ExonBlock(BaseModel):
    val: ExonExonJunctionBlockSizeT


# ---------------------------------------------------------------------------
# DNAT
# ---------------------------------------------------------------------------


class TestDNAT:
    def test_valid_all_bases(self) -> None:
        assert _DNA(seq="ATGC").seq == "ATGC"

    def test_valid_single_base(self) -> None:
        for base in "ATGC":
            assert _DNA(seq=base).seq == base

    def test_invalid_lowercase(self) -> None:
        with pytest.raises(ValidationError):
            _DNA(seq="atgc")

    def test_invalid_ambiguous_base(self) -> None:
        with pytest.raises(ValidationError):
            _DNA(seq="ATGCN")

    def test_invalid_empty_string(self) -> None:
        with pytest.raises(ValidationError):
            _DNA(seq="")


# ---------------------------------------------------------------------------
# GC content types  [0, 100]
# ---------------------------------------------------------------------------


class TestGCContentTypes:
    GCModelCls = type[_GCMin] | type[_GCMax] | type[_GCOpt]

    @pytest.mark.parametrize("cls", [_GCMin, _GCMax, _GCOpt])  # type: ignore[untyped-decorator]
    def test_valid_zero(self, cls: GCModelCls) -> None:
        assert cls(val=0).val == 0

    @pytest.mark.parametrize("cls", [_GCMin, _GCMax, _GCOpt])  # type: ignore[untyped-decorator]
    def test_valid_100(self, cls: GCModelCls) -> None:
        assert cls(val=100).val == 100

    @pytest.mark.parametrize("cls", [_GCMin, _GCMax, _GCOpt])  # type: ignore[untyped-decorator]
    def test_valid_midpoint(self, cls: GCModelCls) -> None:
        assert cls(val=50.5).val == 50.5

    @pytest.mark.parametrize("cls", [_GCMin, _GCMax, _GCOpt])  # type: ignore[untyped-decorator]
    def test_invalid_below_zero(self, cls: GCModelCls) -> None:
        with pytest.raises(ValidationError):
            cls(val=-0.1)

    @pytest.mark.parametrize("cls", [_GCMin, _GCMax, _GCOpt])  # type: ignore[untyped-decorator]
    def test_invalid_above_100(self, cls: GCModelCls) -> None:
        with pytest.raises(ValidationError):
            cls(val=100.1)


# ---------------------------------------------------------------------------
# Tm types — NonNegativeFloat (no upper bound)
# ---------------------------------------------------------------------------


class TestTmTypes:
    TmModelCls = type[_TmMin] | type[_TmOpt] | type[_TmMax]

    @pytest.mark.parametrize("cls", [_TmMin, _TmOpt, _TmMax])  # type: ignore[untyped-decorator]
    def test_valid_zero(self, cls: TmModelCls) -> None:
        assert cls(val=0.0).val == 0.0

    @pytest.mark.parametrize("cls", [_TmMin, _TmOpt, _TmMax])  # type: ignore[untyped-decorator]
    def test_valid_positive(self, cls: TmModelCls) -> None:
        assert cls(val=60.0).val == 60.0

    @pytest.mark.parametrize("cls", [_TmMin, _TmOpt, _TmMax])  # type: ignore[untyped-decorator]
    def test_invalid_negative(self, cls: TmModelCls) -> None:
        with pytest.raises(ValidationError):
            cls(val=-0.1)


# ---------------------------------------------------------------------------
# TSecondaryStructureT — PositiveInt
# ---------------------------------------------------------------------------


class TestTSecondaryStructure:
    def test_valid(self) -> None:
        assert _TSecondary(val=37).val == 37

    def test_invalid_zero(self) -> None:
        with pytest.raises(ValidationError):
            _TSecondary(val=0)

    def test_invalid_negative(self) -> None:
        with pytest.raises(ValidationError):
            _TSecondary(val=-1)


# ---------------------------------------------------------------------------
# FractionT  [0, 1]
# ---------------------------------------------------------------------------


class TestFractionT:
    def test_valid_zero(self) -> None:
        assert _Fraction(val=0.0).val == 0.0

    def test_valid_one(self) -> None:
        assert _Fraction(val=1.0).val == 1.0

    def test_valid_midpoint(self) -> None:
        assert _Fraction(val=0.5).val == 0.5

    def test_invalid_below_zero(self) -> None:
        with pytest.raises(ValidationError):
            _Fraction(val=-0.01)

    def test_invalid_above_one(self) -> None:
        with pytest.raises(ValidationError):
            _Fraction(val=1.01)


# ---------------------------------------------------------------------------
# LengthMinT / LengthMaxT — NonNegativeInt
# ---------------------------------------------------------------------------


class TestLengthTypes:
    LengthModelCls = type[_LengthMin] | type[_LengthMax]

    @pytest.mark.parametrize("cls", [_LengthMin, _LengthMax])  # type: ignore[untyped-decorator]
    def test_valid_zero(self, cls: LengthModelCls) -> None:
        assert cls(val=0).val == 0

    @pytest.mark.parametrize("cls", [_LengthMin, _LengthMax])  # type: ignore[untyped-decorator]
    def test_valid_positive(self, cls: LengthModelCls) -> None:
        assert cls(val=25).val == 25

    @pytest.mark.parametrize("cls", [_LengthMin, _LengthMax])  # type: ignore[untyped-decorator]
    def test_invalid_negative(self, cls: LengthModelCls) -> None:
        with pytest.raises(ValidationError):
            cls(val=-1)


# ---------------------------------------------------------------------------
# FastaFileListT — min_length=1
# ---------------------------------------------------------------------------


class TestFastaFileListT:
    def test_valid_single_file(self) -> None:
        assert _FastaList(files=["genome.fa"]).files == ["genome.fa"]

    def test_valid_multiple_files(self) -> None:
        assert len(_FastaList(files=["a.fa", "b.fa"]).files) == 2

    def test_invalid_empty_list(self) -> None:
        with pytest.raises(ValidationError):
            _FastaList(files=[])


# ---------------------------------------------------------------------------
# ExonExonJunctionBlockSizeT — PositiveInt
# ---------------------------------------------------------------------------


class TestExonExonJunctionBlockSizeT:
    def test_valid_one(self) -> None:
        assert _ExonBlock(val=1).val == 1

    def test_valid_large(self) -> None:
        assert _ExonBlock(val=100).val == 100

    def test_invalid_zero(self) -> None:
        with pytest.raises(ValidationError):
            _ExonBlock(val=0)

    def test_invalid_negative(self) -> None:
        with pytest.raises(ValidationError):
            _ExonBlock(val=-1)
