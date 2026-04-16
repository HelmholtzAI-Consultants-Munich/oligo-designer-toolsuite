"""Tests for readout probe models."""

import unittest
from typing import Any

from pydantic import ValidationError

from oligo_designer_toolsuite.validation.models._general import (
    BaseProbabilities,
    HomopolymerThresholds,
    OligoPropertyWeights,
)
from oligo_designer_toolsuite.validation.models._readout_probes import (
    ReadoutProbeCycleHCR,
    ReadoutProbeFish,
    ReadoutProbeMerfish,
    ReadoutProbeSeqFishPlus,
)

# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------

_BASE_FISH_FIELDS = dict(
    files_fasta_reference_database=["ref.fa"],
    length=20,
    base_probabilities=BaseProbabilities(),
    GC_content_min=40.0,
    GC_content_max=60.0,
    homopolymeric_base_n=HomopolymerThresholds(A=4, T=4, C=4, G=4),
    channels_ids=["Cy3", "Cy5"],
)

_CYCLE_HCR_FIELDS = dict(
    file_readout_probe_table="probes.csv",
    file_codebook=None,
)

_MERFISH_EXTRA = dict(
    set_size=16,
    homogeneous_properties_weights=OligoPropertyWeights(GC_content_oligo=1.0),
    n_bits=16,
    min_hamming_dist=4,
    hamming_weight=4,
)

_SEQ_FISH_PLUS_EXTRA = dict(
    n_barcode_rounds=3,
    n_pseudocolors=4,
)


# ---------------------------------------------------------------------------
# ReadoutProbeCycleHCR — must_be_csv_or_tsv validator
# ---------------------------------------------------------------------------


class TestReadoutProbeCycleHCR(unittest.TestCase):
    def _make(self, **overrides: Any) -> ReadoutProbeCycleHCR:
        fields = {**_CYCLE_HCR_FIELDS, **overrides}
        return ReadoutProbeCycleHCR(**fields)

    def test_valid_csv(self) -> None:
        r = self._make()
        assert r.file_readout_probe_table == "probes.csv"

    def test_valid_tsv(self) -> None:
        r = self._make(file_readout_probe_table="probes.tsv")
        assert r.file_readout_probe_table == "probes.tsv"

    def test_valid_with_codebook_csv(self) -> None:
        r = self._make(file_codebook="book.csv")
        assert r.file_codebook == "book.csv"

    def test_valid_with_codebook_tsv(self) -> None:
        r = self._make(file_codebook="book.tsv")
        assert r.file_codebook == "book.tsv"

    def test_invalid_xlsx_table(self) -> None:
        with self.assertRaises(ValidationError):
            self._make(file_readout_probe_table="probes.xlsx")

    def test_invalid_xlsx_codebook(self) -> None:
        with self.assertRaises(ValidationError):
            self._make(file_codebook="book.xlsx")

    def test_extra_field_forbidden(self) -> None:
        with self.assertRaises(ValidationError):
            self._make(extra="x")


# ---------------------------------------------------------------------------
# ReadoutProbeFish
# ---------------------------------------------------------------------------


class TestReadoutProbeFish(unittest.TestCase):
    def _make(self, **overrides: Any) -> ReadoutProbeFish:
        fields = {**_BASE_FISH_FIELDS, **overrides}
        return ReadoutProbeFish(**fields)

    def test_valid(self) -> None:
        r = self._make()
        assert r.length == 20

    def test_invalid_length_zero(self) -> None:
        with self.assertRaises(ValidationError):
            self._make(length=0)

    def test_extra_field_forbidden(self) -> None:
        with self.assertRaises(ValidationError):
            self._make(unexpected=99)


# ---------------------------------------------------------------------------
# ReadoutProbeMerfish
# ---------------------------------------------------------------------------


class TestReadoutProbeMerfish(unittest.TestCase):
    def _make(self, **overrides: Any) -> ReadoutProbeMerfish:
        fields = {**_BASE_FISH_FIELDS, **_MERFISH_EXTRA, **overrides}
        return ReadoutProbeMerfish(**fields)

    def test_valid(self) -> None:
        r = self._make()
        assert r.set_size == 16

    def test_invalid_set_size_zero(self) -> None:
        with self.assertRaises(ValidationError):
            self._make(set_size=0)

    def test_invalid_n_bits_zero(self) -> None:
        with self.assertRaises(ValidationError):
            self._make(n_bits=0)

    def test_invalid_min_hamming_dist_zero(self) -> None:
        with self.assertRaises(ValidationError):
            self._make(min_hamming_dist=0)

    def test_invalid_hamming_weight_zero(self) -> None:
        with self.assertRaises(ValidationError):
            self._make(hamming_weight=0)


# ---------------------------------------------------------------------------
# ReadoutProbeSeqFishPlus
# ---------------------------------------------------------------------------


class TestReadoutProbeSeqFishPlus(unittest.TestCase):
    def _make(self, **overrides: Any) -> ReadoutProbeSeqFishPlus:
        fields = {**_BASE_FISH_FIELDS, **_SEQ_FISH_PLUS_EXTRA, **overrides}
        return ReadoutProbeSeqFishPlus(**fields)

    def test_valid(self) -> None:
        r = self._make()
        assert r.n_barcode_rounds == 3
        assert r.n_pseudocolors == 4

    def test_invalid_n_barcode_rounds_zero(self) -> None:
        with self.assertRaises(ValidationError):
            self._make(n_barcode_rounds=0)

    def test_invalid_n_pseudocolors_zero(self) -> None:
        with self.assertRaises(ValidationError):
            self._make(n_pseudocolors=0)
