"""Tests for developer parameter models in _developer_parameters.py."""

import unittest
from typing import Any

from pydantic import ValidationError

from oligo_designer_toolsuite.validation.models._developer_parameters import (
    DetectionProbeDevScrinshot,
    DeveloperParametersBase,
    PrimerDevFish,
    ReadoutProbeDevFish,
    ReadoutProbeDevMerfish,
    TargetProbeDev,
    TargetProbeDevMerfish,
    TargetProbeDevOligoSeq,
)

# ---------------------------------------------------------------------------
# Shared minimal-valid building blocks
# ---------------------------------------------------------------------------

_TM_PARAMS = {"mode": "biopython_defaults"}
_TM_CHEM = {"mode": "disabled"}
_TM_SALT = {"mode": "disabled"}
_BLASTN_SEARCH: dict[str, Any] = {}
_BLASTN_HIT = {"coverage": 80.0}

_READOUT_FISH_BASE = dict(
    initial_num_sequences=100,
    specificity_blastn_search_parameters=_BLASTN_SEARCH,
    specificity_blastn_hit_parameters=_BLASTN_HIT,
    cross_hybridization_blastn_search_parameters=_BLASTN_SEARCH,
    cross_hybridization_blastn_hit_parameters=_BLASTN_HIT,
)

_PRIMER_FISH_BASE = dict(
    initial_num_sequences=100,
    Tm_parameters=_TM_PARAMS,
    Tm_chem_correction_parameters=_TM_CHEM,
    Tm_salt_correction_parameters=_TM_SALT,
    specificity_reference_blastn_search_parameters=_BLASTN_SEARCH,
    specificity_reference_blastn_hit_parameters=_BLASTN_HIT,
    specificity_hybridization_probes_blastn_search_parameters=_BLASTN_SEARCH,
    specificity_hybridization_probes_blastn_hit_parameters=_BLASTN_HIT,
)

_OLIGO_SEQ_BASE = dict(
    Tm_parameters=_TM_PARAMS,
    Tm_chem_correction_parameters=_TM_CHEM,
    Tm_salt_correction_parameters=_TM_SALT,
    hybridization_probability={
        "alignment_method": "blastn",
        "search_parameters": {},
        "hit_parameters": {"coverage": 80.0},
        "threshold": 0.5,
    },
    cross_hybridization={
        "alignment_method": "blastn",
        "search_parameters": {},
        "hit_parameters": {"coverage": 80.0},
    },
)

_TARGET_PROBE_DEV_BASE = dict(
    specificity_blastn_search_parameters=_BLASTN_SEARCH,
    specificity_blastn_hit_parameters=_BLASTN_HIT,
    cross_hybridization_blastn_search_parameters=_BLASTN_SEARCH,
    cross_hybridization_blastn_hit_parameters=_BLASTN_HIT,
)

_DEVELOPER_PARAMETERS_BASE = dict(
    oligo_set_selection={
        "n_attempts_graph": 10,
        "n_attempts_clique_enum": 10,
        "diversification_fraction": 0.5,
        "jaccard_opt": 0.5,
        "jaccard_step": 0.1,
    }
)

_READOUT_MERFISH_BASE = {
    **_READOUT_FISH_BASE,
    "Tm_parameters": _TM_PARAMS,
    "Tm_chem_correction_parameters": _TM_CHEM,
    "Tm_salt_correction_parameters": _TM_SALT,
    "n_combinations": 10,
}

_DETECTION_PROBE_DEV_SCRINSHOT_BASE = dict(
    Tm_parameters=_TM_PARAMS,
    Tm_chem_correction_parameters=_TM_CHEM,
    Tm_salt_correction_parameters=_TM_SALT,
)

_TARGET_PROBE_DEV_MERFISH_BASE = {
    **_TARGET_PROBE_DEV_BASE,
    "secondary_structures_threshold_deltaG": -1.0,
    "Tm_parameters": _TM_PARAMS,
    "Tm_chem_correction_parameters": _TM_CHEM,
    "Tm_salt_correction_parameters": _TM_SALT,
}


# ---------------------------------------------------------------------------
# TargetProbeDev
# ---------------------------------------------------------------------------


class TestTargetProbeDev(unittest.TestCase):
    def test_extra_field_forbidden(self) -> None:
        with self.assertRaises(ValidationError):
            fields = {**_TARGET_PROBE_DEV_BASE, "unknown": 99}
            TargetProbeDev(**fields)


# ---------------------------------------------------------------------------
# TargetProbeDevMerfish
# ---------------------------------------------------------------------------


class TestTargetProbeDevMerfish(unittest.TestCase):
    def test_extra_field_forbidden(self) -> None:
        with self.assertRaises(ValidationError):
            fields = {**_TARGET_PROBE_DEV_MERFISH_BASE, "unknown": 99}
            TargetProbeDevMerfish(**fields)


# ---------------------------------------------------------------------------
# DeveloperParametersBase
# ---------------------------------------------------------------------------


class TestDeveloperParametersBase(unittest.TestCase):
    def test_extra_field_forbidden(self) -> None:
        with self.assertRaises(ValidationError):
            fields = {**_DEVELOPER_PARAMETERS_BASE, "unknown": 99}
            DeveloperParametersBase(**fields)


# ---------------------------------------------------------------------------
# ReadoutProbeDevFish
# ---------------------------------------------------------------------------


class TestReadoutProbeDevFish(unittest.TestCase):
    def _make(self, **overrides: Any) -> ReadoutProbeDevFish:
        fields = {**_READOUT_FISH_BASE, **overrides}
        return ReadoutProbeDevFish(**fields)

    def test_valid(self) -> None:
        r = self._make()
        assert r.initial_num_sequences == 100

    def test_invalid_initial_num_sequences_zero(self) -> None:
        with self.assertRaises(ValidationError):
            self._make(initial_num_sequences=0)

    def test_extra_field_forbidden(self) -> None:
        with self.assertRaises(ValidationError):
            self._make(unknown=99)


# ---------------------------------------------------------------------------
# ReadoutProbeDevMerfish
# ---------------------------------------------------------------------------


class TestReadoutProbeDevMerfish(unittest.TestCase):
    def _make(self, **overrides: Any) -> ReadoutProbeDevMerfish:
        fields = {**_READOUT_MERFISH_BASE, **overrides}
        return ReadoutProbeDevMerfish(**fields)

    def test_valid(self) -> None:
        r = self._make()
        assert r.n_combinations == 10

    def test_invalid_n_combinations_zero(self) -> None:
        with self.assertRaises(ValidationError):
            self._make(n_combinations=0)

    def test_extra_field_forbidden(self) -> None:
        with self.assertRaises(ValidationError):
            self._make(unknown=99)


# ---------------------------------------------------------------------------
# PrimerDevFish
# ---------------------------------------------------------------------------


class TestPrimerDevFish(unittest.TestCase):
    def _make(self, **overrides: Any) -> PrimerDevFish:
        fields = {**_PRIMER_FISH_BASE, **overrides}
        return PrimerDevFish(**fields)

    def test_valid(self) -> None:
        p = self._make()
        assert p.initial_num_sequences == 100

    def test_invalid_initial_num_sequences_zero(self) -> None:
        with self.assertRaises(ValidationError):
            self._make(initial_num_sequences=0)

    def test_extra_field_forbidden(self) -> None:
        with self.assertRaises(ValidationError):
            self._make(unknown=99)


# ---------------------------------------------------------------------------
# TargetProbeDevOligoSeq
# ---------------------------------------------------------------------------


class TestTargetProbeDevOligoSeq(unittest.TestCase):
    def test_extra_field_forbidden(self) -> None:
        with self.assertRaises(ValidationError):
            fields = {**_OLIGO_SEQ_BASE, "unknown": 99}
            TargetProbeDevOligoSeq(**fields)


# ---------------------------------------------------------------------------
# DetectionProbeDevScrinshot
# ---------------------------------------------------------------------------


class TestDetectionProbeDevScrinshot(unittest.TestCase):
    def test_extra_field_forbidden(self) -> None:
        with self.assertRaises(ValidationError):
            fields = {**_DETECTION_PROBE_DEV_SCRINSHOT_BASE, "unknown": 99}
            DetectionProbeDevScrinshot(**fields)
