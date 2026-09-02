import unittest
from pathlib import Path

import yaml
from pydantic import AliasChoices, BaseModel, Field, ValidationError

from oligo_designer_toolsuite.config import find_missing_config_keys, format_missing_config_keys
from oligo_designer_toolsuite.config._completeness import MissingConfigKey
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


def _assert_complete(filename: str, model: type[BaseModel]) -> None:
    """Assert that a shipped config sets every parameter instead of relying on a default."""
    raw = _load(filename)
    missing = find_missing_config_keys(model.model_validate(raw), raw)
    assert missing == [], format_missing_config_keys(missing, filename)


# ---------------------------------------------------------------------------
# OligoSeqProbeDesignerConfig
# ---------------------------------------------------------------------------


class TestOligoSeqYaml(unittest.TestCase):
    def test_parses(self) -> None:
        raw = _load("oligo_seq_probe_designer.yaml")
        cfg = OligoSeqProbeDesignerConfig.model_validate(raw)
        assert cfg.schema_version == 2
        assert cfg.general.write_intermediate_steps is True

    def test_config_is_complete(self) -> None:
        _assert_complete("oligo_seq_probe_designer.yaml", OligoSeqProbeDesignerConfig)

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

    def test_config_is_complete(self) -> None:
        _assert_complete("scrinshot_probe_designer.yaml", ScrinshotProbeDesignerConfig)

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

    def test_config_is_complete(self) -> None:
        _assert_complete("hcr_probe_designer.yaml", HcrProbeDesignerConfig)

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

    def test_config_is_complete(self) -> None:
        _assert_complete("cycle_hcr_probe_designer.yaml", CycleHcrProbeDesignerConfig)

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

    def test_config_is_complete(self) -> None:
        _assert_complete("seqfish_plus_probe_designer.yaml", SeqfishPlusProbeDesignerConfig)

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

    def test_config_is_complete(self) -> None:
        _assert_complete("merfish_probe_designer.yaml", MerfishProbeDesignerConfig)

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

    def test_omitted_parameter_resolves_differently_than_an_omitted_section(self) -> None:
        """The bug this check exists for: which default wins depends on how much is left out.

        `nn_table` selects the nearest neighbour table every melting temperature is computed
        from. Leaving it out of a section that is otherwise written resolves to the default on
        the field, while leaving the whole section out resolves to the value merfish declares.
        """
        raw = _load("merfish_probe_designer.yaml")
        without_section = _load("merfish_probe_designer.yaml")
        del raw["target_probes"]["global_parameters"]["Tm_parameters"]["nn_table"]
        del without_section["target_probes"]["global_parameters"]["Tm_parameters"]

        cfg = MerfishProbeDesignerConfig.model_validate(raw)
        cfg_without_section = MerfishProbeDesignerConfig.model_validate(without_section)

        assert cfg.target_probes.global_parameters.Tm_parameters.nn_table is None
        assert cfg_without_section.target_probes.global_parameters.Tm_parameters.nn_table == "DNA_NN4"
        paths = [entry.path for entry in find_missing_config_keys(cfg, raw)]
        assert "target_probes.global_parameters.Tm_parameters.nn_table" in paths


# ---------------------------------------------------------------------------
# find_missing_config_keys
# ---------------------------------------------------------------------------


class TestConfigCompleteness(unittest.TestCase):
    def _merfish(self, raw: dict) -> list[MissingConfigKey]:
        return find_missing_config_keys(MerfishProbeDesignerConfig.model_validate(raw), raw)

    def test_missing_section_reported_once(self) -> None:
        """A whole section is one entry, not one entry per parameter below it."""
        raw = _load("merfish_probe_designer.yaml")
        del raw["general"]
        missing = self._merfish(raw)

        assert [entry.path for entry in missing] == ["general"]
        message = format_missing_config_keys(missing)
        assert "n_jobs" in message and "dir_output" in message

    def test_disabled_branch_demands_no_siblings(self) -> None:
        """A filter that is turned off only has to say so."""
        raw = _load("merfish_probe_designer.yaml")
        raw["target_probes"]["property_filters"]["Tm_filter"] = {"enabled": False}

        assert [entry for entry in self._merfish(raw) if "Tm_filter" in entry.path] == []

    def test_model_allowing_incomplete_is_exempt(self) -> None:
        """An unset BLAST option means the flag is not passed, so it needs no value."""
        raw = _load("merfish_probe_designer.yaml")
        filters = raw["target_probes"]["specificity_filters"]["specificity_blastn_filter"]
        filters["search_parameters"] = {}

        assert [entry for entry in self._merfish(raw) if "search_parameters" in entry.path] == []

    def test_branch_of_a_union_decides_what_is_demanded(self) -> None:
        """`source` selects the branch, and only that branch's parameters are required."""
        raw = _load("merfish_probe_designer.yaml")
        raw["readout_probes"]["codebook"] = {"source": "generate"}
        paths = [entry.path for entry in self._merfish(raw)]

        assert "readout_probes.codebook.n_bits" in paths
        assert "readout_probes.codebook.file" not in paths

    def test_null_counts_as_set(self) -> None:
        """Writing a parameter out as null is setting it; leaving it out is not."""
        raw = _load("merfish_probe_designer.yaml")
        assert self._merfish(raw) == []

        del raw["target_probes"]["global_parameters"]["Tm_parameters"]["c_seq"]

        assert "target_probes.global_parameters.Tm_parameters.c_seq" in [
            entry.path for entry in self._merfish(raw)
        ]

    def test_alias_is_accepted_and_suggested(self) -> None:
        """A field with an alias counts as set under either spelling, and is suggested as the alias."""

        class Aliased(BaseModel):
            evalue: int = Field(
                default=1, validation_alias=AliasChoices("-evalue", "evalue"), serialization_alias="-evalue"
            )

        assert find_missing_config_keys(Aliased.model_validate({"-evalue": 2}), {"-evalue": 2}) == []
        assert find_missing_config_keys(Aliased.model_validate({"evalue": 2}), {"evalue": 2}) == []
        missing = find_missing_config_keys(Aliased.model_validate({}), {})
        assert [(entry.path, entry.key) for entry in missing] == [("evalue", "-evalue")]

    def test_list_of_models_is_indexed(self) -> None:
        """Each element of a list is walked, and its position shows up in the path."""

        class Item(BaseModel):
            a: int = 1
            b: int = 2

        class Outer(BaseModel):
            items: list[Item]

        raw = {"items": [{"a": 1, "b": 2}, {"a": 1}]}

        assert [entry.path for entry in find_missing_config_keys(Outer.model_validate(raw), raw)] == [
            "items[1].b"
        ]

    def test_validation_alone_does_not_check(self) -> None:
        """ODT-Cloud validates the same models and must keep accepting a partial configuration."""
        raw = _load("merfish_probe_designer.yaml")
        del raw["target_probes"]["global_parameters"]["Tm_parameters"]["c_seq"]
        cfg = MerfishProbeDesignerConfig.model_validate(raw)

        # validation is happy, only the completeness check is not
        assert cfg.target_probes.global_parameters.Tm_parameters.c_seq is None
        assert self._merfish(raw) != []
