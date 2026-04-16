"""Tests for detection probe models."""

import unittest
from typing import Any

from pydantic import ValidationError

from oligo_designer_toolsuite.validation.models._detection_probes import DetectionProbeScrinshot

_VALID_FIELDS = dict(
    min_thymines=3,
    length_min=18,
    length_max=40,
    U_distance=5,
    Tm_opt=60.0,
)


class TestDetectionProbeScrinshot(unittest.TestCase):
    def _make(self, **overrides: Any) -> DetectionProbeScrinshot:
        fields = {**_VALID_FIELDS, **overrides}
        return DetectionProbeScrinshot(**fields)

    def test_valid(self) -> None:
        p = self._make()
        assert p.min_thymines == 3

    def test_invalid_min_thymines_negative(self) -> None:
        with self.assertRaises(ValidationError):
            self._make(min_thymines=-1)

    def test_invalid_U_distance_negative(self) -> None:
        with self.assertRaises(ValidationError):
            self._make(U_distance=-1)

    def test_extra_field_forbidden(self) -> None:
        with self.assertRaises(ValidationError):
            self._make(unknown_field=99)
