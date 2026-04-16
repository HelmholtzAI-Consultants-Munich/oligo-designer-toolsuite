"""Tests for primer models in _primer.py."""

import unittest

from pydantic import ValidationError

from oligo_designer_toolsuite.validation.models._primer import (
    PrimerCycleHCR,
    PrimerFish,
    PrimerMerfish,
)

_VALID_FISH_FIELDS = dict(
    files_fasta_reference_database=["ref.fa"],
    reverse_primer_sequence="ATGCATGC",
    length=20,
    base_probabilities={"A": 0.25, "T": 0.25, "C": 0.25, "G": 0.25},
    GC_content_min=40.0,
    GC_content_max=60.0,
    number_GC_GCclamp=2,
    number_three_prime_base_GCclamp=5,
    homopolymeric_base_n={},
    max_len_selfcomplement=6,
    max_len_complement_reverse_primer=6,
    Tm_min=55.0,
    Tm_max=70.0,
    T_secondary_structure=37,
)


class TestPrimerCycleHCRExtraFields(unittest.TestCase):
    def test_extra_field_forbidden(self) -> None:
        with self.assertRaises(ValidationError):
            PrimerCycleHCR(
                forward_primer_sequence="ATGC",
                reverse_primer_sequence="ATGC",
                extra_field="x",  # type: ignore[call-arg]
            )


class TestPrimerFishExtraFields(unittest.TestCase):
    def test_extra_field_forbidden(self) -> None:
        with self.assertRaises(ValidationError):
            PrimerFish(**{**_VALID_FISH_FIELDS, "extra_field": "x"})


class TestPrimerMerfishExtraFields(unittest.TestCase):
    def test_extra_field_forbidden(self) -> None:
        with self.assertRaises(ValidationError):
            PrimerMerfish(
                **{**_VALID_FISH_FIELDS, "secondary_structures_threshold_deltaG": -1.0, "extra_field": "x"}
            )
