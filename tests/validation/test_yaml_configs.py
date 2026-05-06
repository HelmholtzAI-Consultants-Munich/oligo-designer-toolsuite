"""Validate the shipped example YAML config files against their pydantic models.

Each test loads a YAML file from data/configs/ exactly as the pipeline entry
points do (yaml.safe_load → model_validate) and asserts that validation
succeeds and that the parsed object carries expected top-level values.
"""

import unittest
from pathlib import Path

import yaml
from pydantic import ValidationError

from oligo_designer_toolsuite.validation.models.config_pipelines import (
    CycleHCRProbeDesignerConfig,
    GenomicRegionGeneratorConfig,
    HCRProbeDesignerConfig,
    MerfishProbeDesignerConfig,
    OligoSeqProbeDesignerConfig,
    ScrinshotProbeDesignerConfig,
    SeqFishPlusProbeDesignerConfig,
)

# Root of the repository, independent of where pytest is invoked from.
_REPO_ROOT = Path(__file__).parents[2]
_CONFIGS = _REPO_ROOT / "data" / "configs"


def _load(filename: str) -> dict:
    with open(_CONFIGS / filename) as fh:
        data = yaml.safe_load(fh)
        if not isinstance(data, dict):
            raise TypeError
        return data


# ---------------------------------------------------------------------------
# GenomicRegionGeneratorConfig — three source variants
# ---------------------------------------------------------------------------


class TestGenomicRegionGeneratorYaml(unittest.TestCase):
    def test_ensembl_config_parses(self) -> None:
        raw = _load("genomic_region_generator_ensembl.yaml")
        cfg = GenomicRegionGeneratorConfig.model_validate(raw)
        assert cfg.schema_version == 2
        assert cfg.source.source == "ensembl"
        assert cfg.source.parameters.species == "homo_sapiens"
        assert cfg.genomic_regions.exon is True
        assert cfg.exon_exon_junction_block_size == 50

    def test_ncbi_config_parses(self) -> None:
        raw = _load("genomic_region_generator_ncbi.yaml")
        cfg = GenomicRegionGeneratorConfig.model_validate(raw)
        assert cfg.schema_version == 2
        assert cfg.source.source == "ncbi"
        # ncbi species mode: annotation_release is a numeric string
        assert cfg.source.parameters.taxon == "vertebrate_mammalian"  # type: ignore [union-attr]
        assert cfg.genomic_regions.exon_exon_junction is True

    def test_custom_config_parses(self) -> None:
        raw = _load("genomic_region_generator_custom.yaml")
        cfg = GenomicRegionGeneratorConfig.model_validate(raw)
        assert cfg.schema_version == 2
        assert cfg.source.source == "custom"
        assert cfg.source.parameters.species == "Homo_sapiens"
        # all genomic regions enabled
        assert all(
            [
                cfg.genomic_regions.gene,
                cfg.genomic_regions.exon,
                cfg.genomic_regions.intron,
                cfg.genomic_regions.cds,
                cfg.genomic_regions.utr,
                cfg.genomic_regions.intergenic,
                cfg.genomic_regions.exon_exon_junction,
            ]
        )

    def test_ensembl_round_trip(self) -> None:
        raw = _load("genomic_region_generator_ensembl.yaml")
        cfg = GenomicRegionGeneratorConfig.model_validate(raw)
        cfg2 = GenomicRegionGeneratorConfig.model_validate(cfg.model_dump())
        assert cfg == cfg2

    def test_ncbi_round_trip(self) -> None:
        raw = _load("genomic_region_generator_ncbi.yaml")
        cfg = GenomicRegionGeneratorConfig.model_validate(raw)
        cfg2 = GenomicRegionGeneratorConfig.model_validate(cfg.model_dump())
        assert cfg == cfg2

    def test_invalid_schema_version_zero(self) -> None:
        raw = _load("genomic_region_generator_ensembl.yaml")
        raw["schema_version"] = 0
        with self.assertRaises(ValidationError):
            GenomicRegionGeneratorConfig.model_validate(raw)

    def test_extra_field_forbidden(self) -> None:
        raw = _load("genomic_region_generator_ensembl.yaml")
        raw["not_a_field"] = "x"
        with self.assertRaises(ValidationError):
            GenomicRegionGeneratorConfig.model_validate(raw)

    def test_missing_source(self) -> None:
        raw = _load("genomic_region_generator_ensembl.yaml")
        del raw["source"]
        with self.assertRaises(ValidationError):
            GenomicRegionGeneratorConfig.model_validate(raw)

    def test_invalid_source_discriminator(self) -> None:
        raw = _load("genomic_region_generator_ensembl.yaml")
        raw["source"] = {"source": "ftp", "parameters": {}}
        with self.assertRaises(ValidationError):
            GenomicRegionGeneratorConfig.model_validate(raw)


# ---------------------------------------------------------------------------
# CycleHCRProbeDesignerConfig
# ---------------------------------------------------------------------------


class TestCycleHCRYaml(unittest.TestCase):
    def test_parses(self) -> None:
        raw = _load("cycle_hcr_probe_designer.yaml")
        cfg = CycleHCRProbeDesignerConfig.model_validate(raw)
        assert cfg.schema_version == 2
        assert cfg.general.write_intermediate_steps is True
        assert cfg.target_probe.linker_sequence == "TT"
        assert cfg.primer.forward_primer_sequence == "TAATACGACTCACTATAGCGTCATC"
        assert cfg.primer.reverse_primer_sequence == "CGACACCGAACGTGCGACAA"

    def test_readout_probe_file_extension(self) -> None:
        raw = _load("cycle_hcr_probe_designer.yaml")
        cfg = CycleHCRProbeDesignerConfig.model_validate(raw)
        table = cfg.readout_probe.file_readout_probe_table
        assert table.endswith(".tsv") or table.endswith(".csv")

    def test_developer_param_tm_mode(self) -> None:
        raw = _load("cycle_hcr_probe_designer.yaml")
        cfg = CycleHCRProbeDesignerConfig.model_validate(raw)
        assert cfg.developer_param.target_probe.Tm_parameters.mode == "custom"
        assert cfg.developer_param.target_probe.Tm_chem_correction_parameters.mode == "disabled"
        assert cfg.developer_param.target_probe.Tm_salt_correction_parameters.mode == "disabled"

    def test_developer_param_blastn_hit_uses_coverage(self) -> None:
        raw = _load("cycle_hcr_probe_designer.yaml")
        cfg = CycleHCRProbeDesignerConfig.model_validate(raw)
        assert cfg.developer_param.target_probe.specificity_blastn_hit_parameters.coverage is not None
        assert cfg.developer_param.target_probe.specificity_blastn_hit_parameters.min_alignment_length is None

    def test_round_trip(self) -> None:
        raw = _load("cycle_hcr_probe_designer.yaml")
        cfg = CycleHCRProbeDesignerConfig.model_validate(raw)
        cfg2 = CycleHCRProbeDesignerConfig.model_validate(cfg.model_dump())
        assert cfg == cfg2

    def test_invalid_schema_version_zero(self) -> None:
        raw = _load("cycle_hcr_probe_designer.yaml")
        raw["schema_version"] = 0
        with self.assertRaises(ValidationError):
            CycleHCRProbeDesignerConfig.model_validate(raw)

    def test_extra_field_forbidden(self) -> None:
        raw = _load("cycle_hcr_probe_designer.yaml")
        raw["not_a_field"] = "x"
        with self.assertRaises(ValidationError):
            CycleHCRProbeDesignerConfig.model_validate(raw)


# ---------------------------------------------------------------------------
# MerfishProbeDesignerConfig
# ---------------------------------------------------------------------------


class TestMerfishYaml(unittest.TestCase):
    def test_parses(self) -> None:
        raw = _load("merfish_probe_designer.yaml")
        cfg = MerfishProbeDesignerConfig.model_validate(raw)
        assert cfg.schema_version == 2
        assert cfg.general.write_intermediate_steps is True
        assert cfg.target_probe.length_min == 30
        assert cfg.target_probe.length_max == 30

    def test_gc_content_range(self) -> None:
        raw = _load("merfish_probe_designer.yaml")
        cfg = MerfishProbeDesignerConfig.model_validate(raw)
        assert 0 <= cfg.target_probe.GC_content_min <= cfg.target_probe.GC_content_max <= 100

    def test_tm_range(self) -> None:
        raw = _load("merfish_probe_designer.yaml")
        cfg = MerfishProbeDesignerConfig.model_validate(raw)
        assert cfg.target_probe.Tm_min <= cfg.target_probe.Tm_opt <= cfg.target_probe.Tm_max

    def test_round_trip(self) -> None:
        raw = _load("merfish_probe_designer.yaml")
        cfg = MerfishProbeDesignerConfig.model_validate(raw)
        cfg2 = MerfishProbeDesignerConfig.model_validate(cfg.model_dump())
        assert cfg == cfg2

    def test_invalid_schema_version_zero(self) -> None:
        raw = _load("merfish_probe_designer.yaml")
        raw["schema_version"] = 0
        with self.assertRaises(ValidationError):
            MerfishProbeDesignerConfig.model_validate(raw)

    def test_extra_field_forbidden(self) -> None:
        raw = _load("merfish_probe_designer.yaml")
        raw["not_a_field"] = "x"
        with self.assertRaises(ValidationError):
            MerfishProbeDesignerConfig.model_validate(raw)


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

    def test_hybridization_probability_threshold(self) -> None:
        raw = _load("oligo_seq_probe_designer.yaml")
        cfg = OligoSeqProbeDesignerConfig.model_validate(raw)
        hp = cfg.developer_param.target_probe.hybridization_probability
        assert 0.0 <= hp.threshold <= 1.0

    def test_round_trip(self) -> None:
        raw = _load("oligo_seq_probe_designer.yaml")
        cfg = OligoSeqProbeDesignerConfig.model_validate(raw)
        cfg2 = OligoSeqProbeDesignerConfig.model_validate(cfg.model_dump())
        assert cfg == cfg2

    def test_invalid_schema_version_zero(self) -> None:
        raw = _load("oligo_seq_probe_designer.yaml")
        raw["schema_version"] = 0
        with self.assertRaises(ValidationError):
            OligoSeqProbeDesignerConfig.model_validate(raw)

    def test_readout_probe_field_forbidden(self) -> None:
        # OligoSeq has no readout_probe section; extra="forbid" must reject it
        raw = _load("oligo_seq_probe_designer.yaml")
        raw["readout_probe"] = {"file_readout_probe_table": "p.csv"}
        with self.assertRaises(ValidationError):
            OligoSeqProbeDesignerConfig.model_validate(raw)


# ---------------------------------------------------------------------------
# ScrinshotProbeDesignerConfig
# ---------------------------------------------------------------------------


class TestScrinshotYaml(unittest.TestCase):
    def test_parses(self) -> None:
        raw = _load("scrinshot_probe_designer.yaml")
        cfg = ScrinshotProbeDesignerConfig.model_validate(raw)
        assert cfg.schema_version == 2

    def test_detection_probe_fields(self) -> None:
        raw = _load("scrinshot_probe_designer.yaml")
        cfg = ScrinshotProbeDesignerConfig.model_validate(raw)
        dp = cfg.detection_probe
        assert dp.length_min >= 0
        assert dp.length_max >= dp.length_min
        assert dp.Tm_opt >= 0

    def test_round_trip(self) -> None:
        raw = _load("scrinshot_probe_designer.yaml")
        cfg = ScrinshotProbeDesignerConfig.model_validate(raw)
        cfg2 = ScrinshotProbeDesignerConfig.model_validate(cfg.model_dump())
        assert cfg == cfg2

    def test_invalid_schema_version_zero(self) -> None:
        raw = _load("scrinshot_probe_designer.yaml")
        raw["schema_version"] = 0
        with self.assertRaises(ValidationError):
            ScrinshotProbeDesignerConfig.model_validate(raw)

    def test_extra_field_forbidden(self) -> None:
        raw = _load("scrinshot_probe_designer.yaml")
        raw["not_a_field"] = "x"
        with self.assertRaises(ValidationError):
            ScrinshotProbeDesignerConfig.model_validate(raw)


# ---------------------------------------------------------------------------
# SeqFishPlusProbeDesignerConfig
# ---------------------------------------------------------------------------


class TestSeqFishPlusYaml(unittest.TestCase):
    def test_parses(self) -> None:
        raw = _load("seqfish_plus_probe_designer.yaml")
        cfg = SeqFishPlusProbeDesignerConfig.model_validate(raw)
        assert cfg.schema_version == 2
        assert cfg.general.write_intermediate_steps is True

    def test_round_trip(self) -> None:
        raw = _load("seqfish_plus_probe_designer.yaml")
        cfg = SeqFishPlusProbeDesignerConfig.model_validate(raw)
        cfg2 = SeqFishPlusProbeDesignerConfig.model_validate(cfg.model_dump())
        assert cfg == cfg2

    def test_invalid_schema_version_zero(self) -> None:
        raw = _load("seqfish_plus_probe_designer.yaml")
        raw["schema_version"] = 0
        with self.assertRaises(ValidationError):
            SeqFishPlusProbeDesignerConfig.model_validate(raw)

    def test_extra_field_forbidden(self) -> None:
        raw = _load("seqfish_plus_probe_designer.yaml")
        raw["not_a_field"] = "x"
        with self.assertRaises(ValidationError):
            SeqFishPlusProbeDesignerConfig.model_validate(raw)


# ---------------------------------------------------------------------------
# HCRProbeDesignerConfig
# ---------------------------------------------------------------------------


class TestHCRYaml(unittest.TestCase):
    def test_parses(self) -> None:
        raw = _load("hcr_probe_designer.yaml")
        cfg = HCRProbeDesignerConfig.model_validate(raw)
        assert cfg.schema_version == 2
        assert cfg.general.write_intermediate_steps is True
        assert cfg.target_probe.linker_sequence == "AA"

    def test_initiator_probe_file_extension(self) -> None:
        raw = _load("hcr_probe_designer.yaml")
        cfg = HCRProbeDesignerConfig.model_validate(raw)
        table = cfg.initiator_probe.file_initiator_table
        assert table.endswith(".tsv") or table.endswith(".csv")

    def test_developer_param_tm_mode(self) -> None:
        raw = _load("hcr_probe_designer.yaml")
        cfg = HCRProbeDesignerConfig.model_validate(raw)
        assert cfg.developer_param.target_probe.Tm_parameters.mode == "custom"
        assert cfg.developer_param.target_probe.Tm_chem_correction_parameters.mode == "disabled"
        assert cfg.developer_param.target_probe.Tm_salt_correction_parameters.mode == "disabled"

    def test_developer_param_blastn_hit_uses_coverage(self) -> None:
        raw = _load("hcr_probe_designer.yaml")
        cfg = HCRProbeDesignerConfig.model_validate(raw)
        assert cfg.developer_param.target_probe.specificity_blastn_hit_parameters.coverage is not None
        assert cfg.developer_param.target_probe.specificity_blastn_hit_parameters.min_alignment_length is None

    def test_round_trip(self) -> None:
        raw = _load("hcr_probe_designer.yaml")
        cfg = HCRProbeDesignerConfig.model_validate(raw)
        cfg2 = HCRProbeDesignerConfig.model_validate(cfg.model_dump())
        assert cfg == cfg2

    def test_invalid_schema_version_zero(self) -> None:
        raw = _load("hcr_probe_designer.yaml")
        raw["schema_version"] = 0
        with self.assertRaises(ValidationError):
            HCRProbeDesignerConfig.model_validate(raw)

    def test_extra_field_forbidden(self) -> None:
        raw = _load("hcr_probe_designer.yaml")
        raw["not_a_field"] = "x"
        with self.assertRaises(ValidationError):
            HCRProbeDesignerConfig.model_validate(raw)
