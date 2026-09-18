import unittest
from pathlib import Path

import yaml
from pydantic import ValidationError

from oligo_designer_toolsuite.config.pipelines.cycle_hcr_probe_designer import CycleHcrProbeDesignerConfig
from oligo_designer_toolsuite.config.pipelines.hcr_probe_designer import HcrProbeDesignerConfig
from oligo_designer_toolsuite.config.pipelines.merfish_probe_designer import MerfishProbeDesignerConfig
from oligo_designer_toolsuite.config.pipelines.oligo_seq_probe_designer import OligoSeqProbeDesignerConfig
from oligo_designer_toolsuite.config.pipelines.scrinshot_probe_designer import ScrinshotProbeDesignerConfig
from oligo_designer_toolsuite.config.pipelines.seqfish_plus_probe_designer import (
    SeqfishPlusProbeDesignerConfig,
)

_CONFIGS = Path("data/configs")


def _load(filename: str) -> dict:
    with open(_CONFIGS / filename) as fh:
        data = yaml.safe_load(fh)
        if not isinstance(data, dict):
            raise TypeError
        return data


# ---------------------------------------------------------------------------
# OligoSeqProbeDesignerConfig
# ---------------------------------------------------------------------------


class TestOligoSeqYaml(unittest.TestCase):
    def test_parses(self) -> None:
        raw = _load("oligo_seq_probe_designer.yaml")
        cfg = OligoSeqProbeDesignerConfig.model_validate(raw)
        assert cfg.schema_version == 2
        assert cfg.general.write_intermediate_steps is True

    def test_no_readout_probe_key(self) -> None:
        """OligoSeq has no readout_probe; the YAML must not contain it either."""
        raw = _load("oligo_seq_probe_designer.yaml")
        assert "readout_probe" not in raw

    def test_round_trip(self) -> None:
        raw = _load("oligo_seq_probe_designer.yaml")
        cfg = OligoSeqProbeDesignerConfig.model_validate(raw)
        cfg2 = OligoSeqProbeDesignerConfig.model_validate(cfg.model_dump())
        assert cfg == cfg2

    def test_unknown_top_level_field_forbidden(self) -> None:
        # OligoSeq has no readout_probe section; extra="forbid" must reject it
        raw = _load("oligo_seq_probe_designer.yaml")
        raw["readout_probe"] = {"file_readout_probe_table": "p.csv"}
        with self.assertRaises(ValidationError):
            OligoSeqProbeDesignerConfig.model_validate(raw)

    def test_json_schema(self) -> None:
        schema = OligoSeqProbeDesignerConfig.model_json_schema()
        assert schema["title"] == "OligoSeqProbeDesignerConfig"
        assert "required_parameters" in schema["properties"]
        assert "target_probes" in schema["properties"]


# ---------------------------------------------------------------------------
# ScrinshotProbeDesignerConfig
# ---------------------------------------------------------------------------


class TestScrinshotYaml(unittest.TestCase):
    def test_parses(self) -> None:
        raw = _load("scrinshot_probe_designer.yaml")
        cfg = ScrinshotProbeDesignerConfig.model_validate(raw)
        assert cfg.schema_version == 2
        assert cfg.general.write_intermediate_steps is True

    def test_round_trip(self) -> None:
        raw = _load("scrinshot_probe_designer.yaml")
        cfg = ScrinshotProbeDesignerConfig.model_validate(raw)
        cfg2 = ScrinshotProbeDesignerConfig.model_validate(cfg.model_dump())
        assert cfg == cfg2

    def test_unknown_top_level_field_forbidden(self) -> None:
        # Scrinshot has no readout_probe section; extra="forbid" must reject it
        raw = _load("scrinshot_probe_designer.yaml")
        raw["readout_probe"] = {"file_readout_probe_table": "p.csv"}
        with self.assertRaises(ValidationError):
            ScrinshotProbeDesignerConfig.model_validate(raw)

    def test_json_schema(self) -> None:
        schema = ScrinshotProbeDesignerConfig.model_json_schema()
        assert schema["title"] == "ScrinshotProbeDesignerConfig"
        assert "required_parameters" in schema["properties"]
        assert "target_probes" in schema["properties"]
        assert "detection_oligo" in schema["properties"]


# ---------------------------------------------------------------------------
# HcrProbeDesignerConfig
# ---------------------------------------------------------------------------


class TestHcrYaml(unittest.TestCase):
    def test_parses(self) -> None:
        raw = _load("hcr_probe_designer.yaml")
        cfg = HcrProbeDesignerConfig.model_validate(raw)
        assert cfg.schema_version == 2
        assert cfg.general.write_intermediate_steps is True

    def test_round_trip(self) -> None:
        raw = _load("hcr_probe_designer.yaml")
        cfg = HcrProbeDesignerConfig.model_validate(raw)
        cfg2 = HcrProbeDesignerConfig.model_validate(cfg.model_dump())
        assert cfg == cfg2

    def test_unknown_top_level_field_forbidden(self) -> None:
        # HCR has no readout_probe section; extra="forbid" must reject it
        raw = _load("hcr_probe_designer.yaml")
        raw["readout_probe"] = {"file_readout_probe_table": "p.csv"}
        with self.assertRaises(ValidationError):
            HcrProbeDesignerConfig.model_validate(raw)

    def test_json_schema(self) -> None:
        schema = HcrProbeDesignerConfig.model_json_schema()
        assert schema["title"] == "HcrProbeDesignerConfig"
        assert "required_parameters" in schema["properties"]
        assert "target_probes" in schema["properties"]
        assert "initiator_probes" in schema["properties"]


# ---------------------------------------------------------------------------
# CycleHcrProbeDesignerConfig
# ---------------------------------------------------------------------------


class TestCycleHcrYaml(unittest.TestCase):
    def test_parses(self) -> None:
        raw = _load("cycle_hcr_probe_designer.yaml")
        cfg = CycleHcrProbeDesignerConfig.model_validate(raw)
        assert cfg.schema_version == 2
        assert cfg.general.write_intermediate_steps is True

    def test_round_trip(self) -> None:
        raw = _load("cycle_hcr_probe_designer.yaml")
        cfg = CycleHcrProbeDesignerConfig.model_validate(raw)
        cfg2 = CycleHcrProbeDesignerConfig.model_validate(cfg.model_dump())
        assert cfg == cfg2

    def test_unknown_top_level_field_forbidden(self) -> None:
        # CycleHCR has no initiator_probes section; extra="forbid" must reject it
        raw = _load("cycle_hcr_probe_designer.yaml")
        raw["initiator_probes"] = {"codebook": {"source": "load", "file": "c.tsv"}}
        with self.assertRaises(ValidationError):
            CycleHcrProbeDesignerConfig.model_validate(raw)

    def test_codebook_generate_allowed(self) -> None:
        # CycleHCR supports codebook generation (unlike HCR); the default source is "generate"
        raw = _load("cycle_hcr_probe_designer.yaml")
        raw["readout_probes"]["codebook"]["source"] = "generate"
        cfg = CycleHcrProbeDesignerConfig.model_validate(raw)
        assert cfg.readout_probes.codebook.source == "generate"
        assert cfg.readout_probes.codebook.min_hamming_distance in (0, 2, 4)

    def test_readout_probe_table_generate_rejected(self) -> None:
        # Readout probe table generation is not implemented; only "load" is allowed
        raw = _load("cycle_hcr_probe_designer.yaml")
        raw["readout_probes"]["readout_probe_table"]["source"] = "generate"
        with self.assertRaises(ValidationError):
            CycleHcrProbeDesignerConfig.model_validate(raw)

    def test_primer_generate_rejected(self) -> None:
        # Primer generation is not implemented; only "load" is allowed
        raw = _load("cycle_hcr_probe_designer.yaml")
        raw["primers"]["forward_primer"]["source"] = "generate"
        with self.assertRaises(ValidationError):
            CycleHcrProbeDesignerConfig.model_validate(raw)

    def test_json_schema(self) -> None:
        schema = CycleHcrProbeDesignerConfig.model_json_schema()
        assert schema["title"] == "CycleHcrProbeDesignerConfig"
        assert "required_parameters" in schema["properties"]
        assert "target_probes" in schema["properties"]
        assert "readout_probes" in schema["properties"]


# ---------------------------------------------------------------------------
# SeqfishPlusProbeDesignerConfig
# ---------------------------------------------------------------------------


class TestSeqfishYaml(unittest.TestCase):
    def test_parses(self) -> None:
        raw = _load("seqfish_plus_probe_designer.yaml")
        cfg = SeqfishPlusProbeDesignerConfig.model_validate(raw)
        assert cfg.schema_version == 2
        assert cfg.general.write_intermediate_steps is True

    def test_round_trip(self) -> None:
        raw = _load("seqfish_plus_probe_designer.yaml")
        cfg = SeqfishPlusProbeDesignerConfig.model_validate(raw)
        cfg2 = SeqfishPlusProbeDesignerConfig.model_validate(cfg.model_dump())
        assert cfg == cfg2

    def test_unknown_top_level_field_forbidden(self) -> None:
        # seqFISH+ has no initiator_probes section; extra="forbid" must reject it
        raw = _load("seqfish_plus_probe_designer.yaml")
        raw["initiator_probes"] = {"linker_sequence": "AA"}
        with self.assertRaises(ValidationError):
            SeqfishPlusProbeDesignerConfig.model_validate(raw)

    def test_codebook_load_experiment_params_survive_dump(self) -> None:
        # A loaded codebook still needs n_barcode_rounds / n_pseudocolors / channels_ids
        # downstream, so they must survive model_dump() even under extra="ignore".
        raw = _load("seqfish_plus_probe_designer.yaml")
        raw["readout_probes"]["codebook"] = {"source": "load", "file": "codebook.tsv"}
        dumped = SeqfishPlusProbeDesignerConfig.model_validate(raw).model_dump()
        codebook = dumped["readout_probes"]["codebook"]
        assert codebook["n_barcode_rounds"] == 4
        assert codebook["n_pseudocolors"] == 4
        assert codebook["channels_ids"] == ["Alexa488", "Cy3b", "Alexa647"]

    def test_reverse_primer_generate_rejected(self) -> None:
        # Reverse primer generation is not implemented; only "load" is allowed
        raw = _load("seqfish_plus_probe_designer.yaml")
        raw["primers"]["reverse_primer"]["source"] = "generate"
        with self.assertRaises(ValidationError):
            SeqfishPlusProbeDesignerConfig.model_validate(raw)

    def test_base_probabilities_must_sum_to_one(self) -> None:
        raw = _load("seqfish_plus_probe_designer.yaml")
        raw["readout_probes"]["readout_probe_table"]["oligo_generation"]["base_probabilities"] = {
            "A": 0.5,
            "C": 0.5,
            "G": 0.5,
            "T": 0.5,
        }
        with self.assertRaises(ValidationError):
            SeqfishPlusProbeDesignerConfig.model_validate(raw)

    def test_json_schema(self) -> None:
        schema = SeqfishPlusProbeDesignerConfig.model_json_schema()
        assert schema["title"] == "SeqfishPlusProbeDesignerConfig"
        assert "required_parameters" in schema["properties"]
        assert "target_probes" in schema["properties"]
        assert "readout_probes" in schema["properties"]
        assert "primers" in schema["properties"]


# ---------------------------------------------------------------------------
# MerfishProbeDesignerConfig
# ---------------------------------------------------------------------------


class TestMerfishYaml(unittest.TestCase):
    def test_parses(self) -> None:
        raw = _load("merfish_probe_designer.yaml")
        cfg = MerfishProbeDesignerConfig.model_validate(raw)
        assert cfg.schema_version == 2
        assert cfg.general.write_intermediate_steps is True

    def test_round_trip(self) -> None:
        raw = _load("merfish_probe_designer.yaml")
        cfg = MerfishProbeDesignerConfig.model_validate(raw)
        cfg2 = MerfishProbeDesignerConfig.model_validate(cfg.model_dump())
        assert cfg == cfg2

    def test_unknown_top_level_field_forbidden(self) -> None:
        # MERFISH has no initiator_probes section; extra="forbid" must reject it
        raw = _load("merfish_probe_designer.yaml")
        raw["initiator_probes"] = {"linker_sequence": "AA"}
        with self.assertRaises(ValidationError):
            MerfishProbeDesignerConfig.model_validate(raw)

    def test_codebook_load_hamming_weight_survives_dump(self) -> None:
        # A loaded codebook needs hamming_weight downstream and n_bits
        # if the readoutprobe table is generated, so they must
        # survive model_dump().
        raw = _load("merfish_probe_designer.yaml")
        raw["readout_probes"]["codebook"] = {"source": "load", "file": "codebook.tsv"}
        dumped = MerfishProbeDesignerConfig.model_validate(raw).model_dump()
        codebook = dumped["readout_probes"]["codebook"]
        assert codebook["source"] == "load"
        assert codebook["hamming_weight"] == 2
        assert codebook["n_bits"] == 16

    def test_reverse_primer_generate_rejected(self) -> None:
        # Reverse primer generation is not implemented; only "load" is allowed
        raw = _load("merfish_probe_designer.yaml")
        raw["primers"]["reverse_primer"]["source"] = "generate"
        with self.assertRaises(ValidationError):
            MerfishProbeDesignerConfig.model_validate(raw)

    def test_json_schema(self) -> None:
        schema = MerfishProbeDesignerConfig.model_json_schema()
        assert schema["title"] == "MerfishProbeDesignerConfig"
        assert "required_parameters" in schema["properties"]
        assert "target_probes" in schema["properties"]
        assert "readout_probes" in schema["properties"]
        assert "primers" in schema["properties"]
