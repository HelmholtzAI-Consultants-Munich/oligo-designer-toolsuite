"""Tests for target probe models."""

import unittest
from typing import Any

from pydantic import ValidationError

from oligo_designer_toolsuite.validation.models._general import HomopolymerThresholds
from oligo_designer_toolsuite.validation.models._target_probes import (
    TargetProbeBase,
    TargetProbeOligoSeq,
    TargetProbeScrinshot,
    TargetProbeSplitProbeBase,
)

# ---------------------------------------------------------------------------
# Minimal valid field sets
# ---------------------------------------------------------------------------

_BASE_FIELDS = dict(
    file_regions=None,
    files_fasta_database=["db.fa"],
    files_fasta_reference_database=["ref.fa"],
    isoform_consensus=80.0,
    GC_content_min=40.0,
    GC_content_max=60.0,
    homopolymeric_base_n=HomopolymerThresholds(A=4, T=4, C=4, G=4),
    set_size_min=5,
    set_size_opt=10,
    distance_between_target_probes=2,
    n_sets=1,
)

_SPLIT_PROBE_EXTRA = dict(
    L_probe_sequence_length=25,
    gap_sequence_length=0,
    R_probe_sequence_length=25,
    Tm_min=60.0,
    Tm_max=72.0,
    T_secondary_structure=37,
    junction_region_size=0,
    isoform_weight=1.0,
    linker_sequence="ATGCATGC",
)

_OLIGOSEQ_EXTRA = dict(
    files_vcf_reference_database=[],
    length_min=18,
    length_max=22,
    split_region=100,
    Tm_min=60.0,
    Tm_max=72.0,
    T_secondary_structure=37,
    secondary_structures_threshold_deltaG=-1.0,
    kmer_abundance_threshold={},
    prohibited_sequences=[],
    max_len_selfcomplement=6,
    read_length_bias=0,
    uniform_distance_weight=1.0,
    isoform_weight=1.0,
    targeted_exons_weight=1.0,
    targeted_exons=[],
    GC_weight=1.0,
    GC_content_opt=50.0,
    Tm_weight=1.0,
    Tm_opt=66.0,
)

_SCRINSHOT_EXTRA = dict(
    length_min=18,
    length_max=40,
    GC_content_opt=50.0,
    Tm_min=60.0,
    Tm_opt=66.0,
    Tm_max=72.0,
    padlock_arm_Tm_dif_max=5,
    padlock_arm_length_min=10,
    padlock_arm_Tm_min=55.0,
    padlock_arm_Tm_max=70.0,
    ligation_region_size=5,
    isoform_weight=1.0,
    GC_weight=1.0,
    Tm_weight=1.0,
)


# ---------------------------------------------------------------------------
# TargetProbeBase — shared field validation
# ---------------------------------------------------------------------------


class TestTargetProbeBaseFields(unittest.TestCase):
    """Test TargetProbeBase constraints via TargetProbeCycleHCR."""

    def _make(self, **overrides: Any) -> TargetProbeBase:
        fields = {**_BASE_FIELDS, **overrides}
        return TargetProbeBase(**fields)

    def test_valid(self) -> None:
        p = self._make()
        assert p.set_size_min == 5

    def test_file_regions_can_be_none(self) -> None:
        p = self._make(file_regions=None)
        assert p.file_regions is None

    def test_file_regions_string(self) -> None:
        p = self._make(file_regions="genes.txt")
        assert p.file_regions == "genes.txt"

    def test_invalid_isoform_consensus_negative(self) -> None:
        with self.assertRaises(ValidationError):
            self._make(isoform_consensus=-1.0)

    def test_invalid_isoform_consensus_above_100(self) -> None:
        with self.assertRaises(ValidationError):
            self._make(isoform_consensus=101.0)

    def test_invalid_set_size_min_zero(self) -> None:
        with self.assertRaises(ValidationError):
            self._make(set_size_min=0)

    def test_invalid_set_size_opt_zero(self) -> None:
        with self.assertRaises(ValidationError):
            self._make(set_size_opt=0)

    def test_invalid_n_sets_zero(self) -> None:
        with self.assertRaises(ValidationError):
            self._make(n_sets=0)

    def test_distance_between_target_probes_accepts_negative(self) -> None:
        # Plain int: negative values are valid (negative = allowed overlap)
        p = self._make(distance_between_target_probes=-5)
        assert p.distance_between_target_probes == -5

    def test_extra_field_forbidden(self) -> None:
        with self.assertRaises(ValidationError):
            self._make(unexpected_field=99)


# ---------------------------------------------------------------------------
# TargetProbeSplitProbeBase — specific fields (don't additionally check
# TargetProbeCycleHCR as the only extra field is a weight that I haven't tested
# for any model)
# ---------------------------------------------------------------------------


class TestTargetProbeSplitProbeBase(unittest.TestCase):
    def _make(self, **overrides: Any) -> TargetProbeSplitProbeBase:
        fields = {**_BASE_FIELDS, **_SPLIT_PROBE_EXTRA, **overrides}
        return TargetProbeSplitProbeBase(**fields)

    def test_valid(self) -> None:
        p = self._make()
        assert p.linker_sequence == "ATGCATGC"

    def test_invalid_L_probe_sequence_length_zero(self) -> None:
        with self.assertRaises(ValidationError):
            self._make(L_probe_sequence_length=0)

    def test_invalid_gap_sequence_length(self) -> None:
        with self.assertRaises(ValidationError):
            self._make(gap_sequence_length=-1)

    def test_invalid_R_probe_sequence_length_zero(self) -> None:
        with self.assertRaises(ValidationError):
            self._make(R_probe_sequence_length=0)

    def test_invalid_junction_region_size(self) -> None:
        with self.assertRaises(ValidationError):
            self._make(junction_region_size=-1)


# ---------------------------------------------------------------------------
# TargetProbeMerfish — on specific fields with any types that are not tested already
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# TargetProbeScrinshot — specific fields
# ---------------------------------------------------------------------------


class TestTargetProbeScrinshot(unittest.TestCase):
    def _make(self, **overrides: Any) -> TargetProbeScrinshot:
        fields = {**_BASE_FIELDS, **_SCRINSHOT_EXTRA, **overrides}
        return TargetProbeScrinshot(**fields)

    def test_valid(self) -> None:
        p = self._make()
        assert p.padlock_arm_length_min == 10

    def test_invalid_padlock_arm_Tm_dif_max(self) -> None:
        with self.assertRaises(ValidationError):
            self._make(padlock_arm_Tm_dif_max=-1)

    def test_invalid_padlock_arm_length_min(self) -> None:
        with self.assertRaises(ValidationError):
            self._make(padlock_arm_length_min=0)

    def test_invalid_padlock_arm_Tm_min_negative(self) -> None:
        with self.assertRaises(ValidationError):
            self._make(padlock_arm_Tm_min=-1.0)

    def test_invalid_padlock_arm_Tm_max_negative(self) -> None:
        with self.assertRaises(ValidationError):
            self._make(padlock_arm_Tm_max=-1.0)

    def test_invalid_ligation_region_size_zero(self) -> None:
        with self.assertRaises(ValidationError):
            self._make(ligation_region_size=0)


# ---------------------------------------------------------------------------
# TargetProbeOligoSeq — specific fields (built-in / pydantic types only)
# ---------------------------------------------------------------------------


class TestTargetProbeOligoSeq(unittest.TestCase):
    def _make(self, **overrides: Any) -> TargetProbeOligoSeq:
        fields = {**_BASE_FIELDS, **_OLIGOSEQ_EXTRA, **overrides}
        return TargetProbeOligoSeq(**fields)

    def test_invalid_split_region_zero(self) -> None:
        with self.assertRaises(ValidationError):
            self._make(split_region=0)

    def test_invalid_kmer_abundance_threshold_key_zero(self) -> None:
        with self.assertRaises(ValidationError):
            self._make(kmer_abundance_threshold={0: 0.5})

    def test_invalid_max_len_selfcomplement_negative(self) -> None:
        with self.assertRaises(ValidationError):
            self._make(max_len_selfcomplement=-1)

    def test_invalid_read_length_bias_negative(self) -> None:
        with self.assertRaises(ValidationError):
            self._make(read_length_bias=-1)
