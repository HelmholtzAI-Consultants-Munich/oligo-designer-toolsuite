import unittest
from pathlib import Path

import yaml
from pydantic import ValidationError

from oligo_designer_toolsuite.config.pipelines.genomic_region_generator import (
    AnnotationNcbiAssembly,
    AnnotationNcbiSpecies,
    GenomicRegionGeneratorConfig,
)
from oligo_designer_toolsuite.config.pipelines.oligo_seq_probe_designer import OligoSeqProbeDesignerConfig

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
# GenomicRegionGeneratorConfig
# ---------------------------------------------------------------------------


class TestGenomicRegionGeneratorYaml(unittest.TestCase):
    _FILES = {
        "custom": "genomic_region_generator_custom.yaml",
        "ensembl": "genomic_region_generator_ensembl.yaml",
        "ncbi": "genomic_region_generator_ncbi.yaml",
    }

    def _cfg(self, source: str) -> GenomicRegionGeneratorConfig:
        return GenomicRegionGeneratorConfig.model_validate(_load(self._FILES[source]))

    def _ncbi_raw(self, **param_overrides: str) -> dict:
        raw = _load(self._FILES["ncbi"])
        raw["annotation"]["parameters"].update(param_overrides)
        return raw

    def test_parses_custom(self) -> None:
        cfg = self._cfg("custom")
        assert cfg.schema_version == 2
        assert cfg.annotation.source == "custom"
        assert cfg.annotation.parameters.file_annotation.endswith(".gtf")
        assert cfg.genomic_regions.gene is True

    def test_parses_ensembl(self) -> None:
        cfg = self._cfg("ensembl")
        assert cfg.schema_version == 2
        assert cfg.annotation.source == "ensembl"
        assert cfg.annotation.parameters.species == "homo_sapiens"
        assert cfg.genomic_regions.exon is True

    def test_parses_ncbi(self) -> None:
        cfg = self._cfg("ncbi")
        assert cfg.schema_version == 2
        assert cfg.annotation.source == "ncbi"
        assert cfg.annotation.mode == "species"
        assert cfg.annotation.parameters.taxon == "vertebrate_mammalian"
        assert cfg.annotation.parameters.assembly_source == "auto"

    def test_round_trip_custom(self) -> None:
        cfg = self._cfg("custom")
        cfg2 = GenomicRegionGeneratorConfig.model_validate(cfg.model_dump())
        assert cfg == cfg2

    def test_round_trip_ensembl(self) -> None:
        cfg = self._cfg("ensembl")
        cfg2 = GenomicRegionGeneratorConfig.model_validate(cfg.model_dump())
        assert cfg == cfg2

    def test_round_trip_ncbi(self) -> None:
        cfg = self._cfg("ncbi")
        cfg2 = GenomicRegionGeneratorConfig.model_validate(cfg.model_dump())
        assert cfg == cfg2

    def test_unknown_top_level_field_forbidden_custom(self) -> None:
        raw = _load(self._FILES["custom"])
        raw["bogus"] = 1
        with self.assertRaises(ValidationError):
            GenomicRegionGeneratorConfig.model_validate(raw)

    def test_unknown_top_level_field_forbidden_ensembl(self) -> None:
        raw = _load(self._FILES["ensembl"])
        raw["bogus"] = 1
        with self.assertRaises(ValidationError):
            GenomicRegionGeneratorConfig.model_validate(raw)

    def test_unknown_top_level_field_forbidden_ncbi(self) -> None:
        raw = _load(self._FILES["ncbi"])
        raw["bogus"] = 1
        with self.assertRaises(ValidationError):
            GenomicRegionGeneratorConfig.model_validate(raw)

    def test_json_schema(self) -> None:
        # The JSON schema does not vary by pipeline; assert the model shape and
        # that the annotation discriminated union exposes all three source branches.
        schema = GenomicRegionGeneratorConfig.model_json_schema()
        assert schema["title"] == "GenomicRegionGeneratorConfig"
        for field in ("annotation", "genomic_regions", "dir_output", "exon_exon_junction_block_size"):
            assert field in schema["properties"]
        sources = {
            defn["properties"]["source"]["const"]
            for defn in schema["$defs"].values()
            if "properties" in defn
            and "source" in defn["properties"]
            and "const" in defn["properties"]["source"]
        }
        assert {"custom", "ensembl", "ncbi"} <= sources

    def test_source_discriminator_selects_ncbi_species(self) -> None:
        cfg = self._cfg("ncbi")
        assert cfg.annotation.source == "ncbi"
        assert cfg.annotation.mode == "species"
        assert hasattr(cfg.annotation.parameters, "taxon")

    def test_ncbi_assembly_mode_parses(self) -> None:
        raw = _load(self._FILES["ncbi"])
        raw["annotation"]["mode"] = "assembly"
        raw["annotation"]["parameters"] = {
            "refseq_assembly_accession": "GCF_000001405.38",
            "assembly_name": "GRCh38.p12",
        }
        cfg = GenomicRegionGeneratorConfig.model_validate(raw)
        assert isinstance(cfg.annotation, AnnotationNcbiAssembly)
        assert cfg.annotation.mode == "assembly"
        assert cfg.annotation.parameters.assembly_name == "GRCh38.p12"

    def test_dir_output_default_derived(self) -> None:
        for source in self._FILES:
            raw = _load(self._FILES[source])
            raw["dir_output"] = None
            cfg = GenomicRegionGeneratorConfig.model_validate(raw)
            assert cfg.dir_output == f"output_genomic_region_generator_{source}"

    # --- test cross-field validators

    def test_assembly_source_not_available_for_taxon(self) -> None:
        # bacteria supports only {latest_assembly_versions, reference}
        raw = self._ncbi_raw(taxon="bacteria", assembly_source="annotation_releases")
        with self.assertRaises(ValidationError):
            GenomicRegionGeneratorConfig.model_validate(raw)

    def test_numeric_release_requires_annotation_releases(self) -> None:
        raw = self._ncbi_raw(
            taxon="vertebrate_mammalian",
            assembly_source="latest_assembly_versions",
            annotation_release="110",
        )
        with self.assertRaises(ValidationError):
            GenomicRegionGeneratorConfig.model_validate(raw)

    def test_numeric_release_auto_without_annotation_releases(self) -> None:
        # viral supports only {latest_assembly_versions}, so auto cannot use a numeric release
        raw = self._ncbi_raw(taxon="viral", assembly_source="auto", annotation_release="110")
        with self.assertRaises(ValidationError):
            GenomicRegionGeneratorConfig.model_validate(raw)

    def test_annotation_releases_with_numeric_release_ok(self) -> None:
        raw = self._ncbi_raw(
            taxon="vertebrate_mammalian",
            assembly_source="annotation_releases",
            annotation_release="110",
        )
        cfg = GenomicRegionGeneratorConfig.model_validate(raw)
        assert isinstance(cfg.annotation, AnnotationNcbiSpecies)
        assert cfg.annotation.parameters.annotation_release == "110"

    def test_block_size_none_rejected(self) -> None:
        # ncbi config has exon_exon_junction=True, so a null block size must be rejected.
        raw = _load(self._FILES["ncbi"])
        raw["exon_exon_junction_block_size"] = None
        with self.assertRaises(ValidationError):
            GenomicRegionGeneratorConfig.model_validate(raw)
