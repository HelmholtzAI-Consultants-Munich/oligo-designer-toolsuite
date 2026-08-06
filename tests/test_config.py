import unittest
from pathlib import Path

import yaml
from pydantic import ValidationError

from oligo_designer_toolsuite.config.pipelines.hcr_probe_designer import HcrProbeDesignerConfig
from oligo_designer_toolsuite.config.pipelines.oligo_seq_probe_designer import OligoSeqProbeDesignerConfig
from oligo_designer_toolsuite.config.pipelines.scrinshot_probe_designer import ScrinshotProbeDesignerConfig

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
        assert "target_probe" in schema["properties"]


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
        assert "target_probe" in schema["properties"]
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
        assert "target_probes" in schema["properties"]
        assert "initiator_probes" in schema["properties"]
