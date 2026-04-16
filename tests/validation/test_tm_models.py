"""Tests for Tm parameter models and their discriminated unions."""

import unittest

from pydantic import TypeAdapter, ValidationError

from oligo_designer_toolsuite.validation.models._general import (
    TmChemCorrectionParameters,
    TmChemCorrectionParametersBiopythonDefaults,
    TmChemCorrectionParametersCustom,
    TmChemCorrectionParametersDetails,
    TmChemCorrectionParametersDisabled,
    TmParameters,
    TmParametersBiopythonDefaults,
    TmParametersCustom,
    TmParametersDetails,
    TmSaltCorrectionParameters,
    TmSaltCorrectionParametersBiopythonDefaults,
    TmSaltCorrectionParametersCustom,
    TmSaltCorrectionParametersDetails,
    TmSaltCorrectionParametersDisabled,
)

_TmParametersAdapter: TypeAdapter[TmParameters] = TypeAdapter(TmParameters)
_TmChemAdapter: TypeAdapter[TmChemCorrectionParameters] = TypeAdapter(TmChemCorrectionParameters)
_TmSaltAdapter: TypeAdapter[TmSaltCorrectionParameters] = TypeAdapter(TmSaltCorrectionParameters)


# ---------------------------------------------------------------------------
# TmParametersDetails
# ---------------------------------------------------------------------------


class TestTmParametersDetails(unittest.TestCase):
    def test_valid_single_field(self) -> None:
        d = TmParametersDetails(nn_table="DNA_NN3")
        assert d.nn_table == "DNA_NN3"

    def test_valid_multiple_fields(self) -> None:
        d = TmParametersDetails(check=True, strict=False, nn_table="DNA_NN1", dnac1=250, Na=None)
        assert d.dnac1 == 250

    def test_invalid_empty_input(self) -> None:
        with self.assertRaisesRegex(ValidationError, "At least one parameter must be provided"):
            TmParametersDetails()

    def test_invalid_all_none(self) -> None:
        with self.assertRaisesRegex(ValidationError, "At least one parameter must be provided"):
            TmParametersDetails(strict=None, check=None)

    def test_invalid_saltcorr_above_7(self) -> None:
        with self.assertRaises(ValidationError):
            TmParametersDetails(saltcorr=8)

    def test_invalid_saltcorr_negative(self) -> None:
        with self.assertRaises(ValidationError):
            TmParametersDetails(saltcorr=-1)

    def test_valid_saltcorr_boundaries(self) -> None:
        assert TmParametersDetails(saltcorr=0).saltcorr == 0
        assert TmParametersDetails(saltcorr=7).saltcorr == 7

    def test_invalid_nn_table(self) -> None:
        with self.assertRaises(ValidationError):
            TmParametersDetails(nn_table="INVALID_TABLE")

    def test_valid_all_nn_tables(self) -> None:
        for table in ("DNA_NN1", "DNA_NN2", "DNA_NN3", "DNA_NN4"):
            d = TmParametersDetails(nn_table=table)
            assert d.nn_table == table

    def test_invalid_tmm_table(self) -> None:
        with self.assertRaises(ValidationError):
            TmParametersDetails(tmm_table="DNA_TMM2")

    def test_valid_all_tmm_tables(self) -> None:
        for table in ("DNA_TMM1",):
            d = TmParametersDetails(tmm_table=table)
            assert d.tmm_table == table

    def test_invalid_imm_table(self) -> None:
        with self.assertRaises(ValidationError):
            TmParametersDetails(imm_table="DNA_IMM2")

    def test_valid_all_imm_tables(self) -> None:
        for table in ("DNA_IMM1",):
            d = TmParametersDetails(imm_table=table)
            assert d.imm_table == table

    def test_invalid_de_table(self) -> None:
        with self.assertRaises(ValidationError):
            TmParametersDetails(de_table="DNA_DE2")

    def test_valid_all_de_tables(self) -> None:
        for table in ("DNA_DE1",):
            d = TmParametersDetails(de_table=table)
            assert d.de_table == table

    def test_invalid_negative_ion_concentration(self) -> None:
        with self.assertRaises(ValidationError):
            TmParametersDetails(Na=-1)

    def test_extra_field_forbidden(self) -> None:
        with self.assertRaises(ValidationError):
            TmParametersDetails(check=True, unknown=99)  # type: ignore[call-arg]


# ---------------------------------------------------------------------------
# TmParameters union — discriminator="mode"
# ---------------------------------------------------------------------------


class TestTmParametersUnion(unittest.TestCase):
    def test_biopython_defaults_mode(self) -> None:
        obj = _TmParametersAdapter.validate_python({"mode": "biopython_defaults"})
        assert isinstance(obj, TmParametersBiopythonDefaults)

    def test_custom_mode(self) -> None:
        obj = _TmParametersAdapter.validate_python({"mode": "custom", "parameters": {"nn_table": "DNA_NN3"}})
        assert isinstance(obj, TmParametersCustom)
        assert obj.parameters.nn_table == "DNA_NN3"

    def test_invalid_mode(self) -> None:
        with self.assertRaises(ValidationError):
            _TmParametersAdapter.validate_python({"mode": "nonsense"})

    def test_custom_without_parameters_fails(self) -> None:
        with self.assertRaises(ValidationError):
            _TmParametersAdapter.validate_python({"mode": "custom"})


# ---------------------------------------------------------------------------
# TmChemCorrectionParametersDetails
# ---------------------------------------------------------------------------


class TestTmChemCorrectionParametersDetails(unittest.TestCase):
    def test_valid_single_field(self) -> None:
        d = TmChemCorrectionParametersDetails(DMSO=5.0)
        assert d.DMSO == 5.0

    def test_invalid_empty_input(self) -> None:
        with self.assertRaisesRegex(ValidationError, "At least one parameter must be provided"):
            TmChemCorrectionParametersDetails()

    def test_invalid_all_none(self) -> None:
        with self.assertRaisesRegex(ValidationError, "At least one parameter must be provided"):
            TmChemCorrectionParametersDetails(DMSO=None, fmd=None)

    def test_invalid_DMSO_above_100(self) -> None:
        with self.assertRaises(ValidationError):
            TmChemCorrectionParametersDetails(DMSO=101.0)

    def test_invalid_DMSO_negative(self) -> None:
        with self.assertRaises(ValidationError):
            TmChemCorrectionParametersDetails(DMSO=-1.0)

    def test_invalid_GC_above_100(self) -> None:
        with self.assertRaises(ValidationError):
            TmChemCorrectionParametersDetails(GC=101.0)

    def test_invalid_GC_negative(self) -> None:
        with self.assertRaises(ValidationError):
            TmChemCorrectionParametersDetails(GC=-1.0)

    def test_invalid_fmdmethod_zero(self) -> None:
        with self.assertRaises(ValidationError):
            TmChemCorrectionParametersDetails(fmdmethod=0)

    def test_invalid_fmdmethod_above_2(self) -> None:
        with self.assertRaises(ValidationError):
            TmChemCorrectionParametersDetails(fmdmethod=3)

    # _check_fmd_vs_method — method 1: fmd must be [0, 100]
    def test_fmdmethod1_fmd_in_range(self) -> None:
        d = TmChemCorrectionParametersDetails(fmdmethod=1, fmd=50.0)
        assert d.fmd == 50.0

    def test_fmdmethod1_fmd_above_100(self) -> None:
        with self.assertRaisesRegex(ValidationError, "fmdmethod=1"):
            TmChemCorrectionParametersDetails(fmdmethod=1, fmd=101.0)

    def test_fmdmethod1_fmd_negative(self) -> None:
        with self.assertRaisesRegex(ValidationError, "fmdmethod=1"):
            TmChemCorrectionParametersDetails(fmdmethod=1, fmd=-1.0)

    # _check_fmd_vs_method — method 2: GC required, fmd must be >= 0
    def test_fmdmethod2_requires_GC(self) -> None:
        with self.assertRaisesRegex(ValidationError, "GC must be provided"):
            TmChemCorrectionParametersDetails(fmdmethod=2)

    def test_fmdmethod2_with_GC_valid(self) -> None:
        d = TmChemCorrectionParametersDetails(fmdmethod=2, GC=50.0)
        assert d.GC == 50.0

    def test_fmdmethod2_fmd_negative(self) -> None:
        with self.assertRaisesRegex(ValidationError, "fmdmethod=2"):
            TmChemCorrectionParametersDetails(fmdmethod=2, GC=50.0, fmd=-0.1)

    def test_extra_field_forbidden(self) -> None:
        with self.assertRaises(ValidationError):
            TmChemCorrectionParametersDetails(DMSO=5.0, unknown=1)  # type: ignore[call-arg]


# ---------------------------------------------------------------------------
# TmChemCorrectionParameters union — discriminator="mode"
# ---------------------------------------------------------------------------


class TestTmChemCorrectionParametersUnion(unittest.TestCase):
    def test_disabled_mode(self) -> None:
        obj = _TmChemAdapter.validate_python({"mode": "disabled"})
        assert isinstance(obj, TmChemCorrectionParametersDisabled)

    def test_biopython_defaults_mode(self) -> None:
        obj = _TmChemAdapter.validate_python({"mode": "biopython_defaults"})
        assert isinstance(obj, TmChemCorrectionParametersBiopythonDefaults)

    def test_custom_mode(self) -> None:
        obj = _TmChemAdapter.validate_python({"mode": "custom", "parameters": {"DMSO": 5.0}})
        assert isinstance(obj, TmChemCorrectionParametersCustom)

    def test_invalid_mode(self) -> None:
        with self.assertRaises(ValidationError):
            _TmChemAdapter.validate_python({"mode": "bad"})


# ---------------------------------------------------------------------------
# TmSaltCorrectionParametersDetails
# ---------------------------------------------------------------------------


class TestTmSaltCorrectionParametersDetails(unittest.TestCase):
    def test_valid_single_field(self) -> None:
        d = TmSaltCorrectionParametersDetails(Na=50)
        assert d.Na == 50

    def test_invalid_empty_input(self) -> None:
        with self.assertRaisesRegex(ValidationError, "At least one parameter must be provided"):
            TmSaltCorrectionParametersDetails()

    def test_invalid_all_none(self) -> None:
        with self.assertRaisesRegex(ValidationError, "At least one parameter must be provided"):
            TmSaltCorrectionParametersDetails(Na=None, method=None)

    def test_invalid_method_zero(self) -> None:
        with self.assertRaises(ValidationError):
            TmSaltCorrectionParametersDetails(method=0)

    def test_invalid_method_above_7(self) -> None:
        with self.assertRaises(ValidationError):
            TmSaltCorrectionParametersDetails(method=8)

    def test_invalid_negative_ion_concentration(self) -> None:
        with self.assertRaises(ValidationError):
            TmSaltCorrectionParametersDetails(Na=-1)

    def test_extra_field_forbidden(self) -> None:
        with self.assertRaises(ValidationError):
            TmSaltCorrectionParametersDetails(Na=50, unknown=1)  # type: ignore[call-arg]


# ---------------------------------------------------------------------------
# TmSaltCorrectionParameters union — discriminator="mode"
# ---------------------------------------------------------------------------


class TestTmSaltCorrectionParametersUnion(unittest.TestCase):
    def test_disabled_mode(self) -> None:
        obj = _TmSaltAdapter.validate_python({"mode": "disabled"})
        assert isinstance(obj, TmSaltCorrectionParametersDisabled)

    def test_biopython_defaults_mode(self) -> None:
        obj = _TmSaltAdapter.validate_python({"mode": "biopython_defaults"})
        assert isinstance(obj, TmSaltCorrectionParametersBiopythonDefaults)

    def test_custom_mode(self) -> None:
        obj = _TmSaltAdapter.validate_python({"mode": "custom", "parameters": {"Na": 100}})
        assert isinstance(obj, TmSaltCorrectionParametersCustom)

    def test_invalid_mode(self) -> None:
        with self.assertRaises(ValidationError):
            _TmSaltAdapter.validate_python({"mode": "wrong"})
