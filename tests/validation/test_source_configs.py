"""Tests for source configuration models and their discriminated union."""

import warnings

import pytest
from pydantic import TypeAdapter, ValidationError

from oligo_designer_toolsuite.utils._checkers_and_helpers import ConfigurationError
from oligo_designer_toolsuite.validation.models._general import (
    SourceConfigs,
    SourceCustom,
    SourceEnsembl,
    SourceNcbi,
    SourceNcbiAssembly,
    SourceNcbiSpecies,
    SourceParamsCustom,
    SourceParamsEnsembl,
    SourceParamsNcbiAssembly,
    SourceParamsNcbiSpecies,
)

_SourceConfigsAdapter: TypeAdapter[SourceConfigs] = TypeAdapter(SourceConfigs)
_SourceNcbiAdapter: TypeAdapter[SourceNcbi] = TypeAdapter(SourceNcbi)


# ---------------------------------------------------------------------------
# SourceParamsCustom
# ---------------------------------------------------------------------------


class TestSourceParamsCustom:
    def test_valid_minimal(self) -> None:
        # ignore mypy error as the default values are provided by a field validator
        p = SourceParamsCustom(file_annotation="ann.gtf", file_sequence="seq.fa")  # type: ignore[call-arg]
        assert p.file_annotation == "ann.gtf"
        assert p.files_source == "custom"  # default
        assert p.species == "unknown"  # default

    def test_valid_all_fields(self) -> None:
        p = SourceParamsCustom(
            file_annotation="ann.gtf",
            file_sequence="seq.fa",
            files_source="Ensembl",
            species="homo_sapiens",
            annotation_release="109",
            genome_assembly="GRCh38",
        )
        assert p.species == "homo_sapiens"

    def test_none_replaced_with_default_files_source(self) -> None:
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            p = SourceParamsCustom(file_annotation="a.gtf", file_sequence="s.fa", files_source=None)  # type: ignore[call-arg]
            assert p.files_source == "custom"
            assert any("files_source" in str(warning.message) for warning in w)

    def test_extra_field_forbidden(self) -> None:
        with pytest.raises(ValidationError):
            SourceParamsCustom(file_annotation="a.gtf", file_sequence="s.fa", extra_field="x")  # type: ignore[call-arg]


# ---------------------------------------------------------------------------
# SourceParamsEnsembl
# ---------------------------------------------------------------------------


class TestSourceParamsEnsembl:
    def test_valid(self) -> None:
        p = SourceParamsEnsembl(species="homo_sapiens", annotation_release="current")
        assert (
            p.annotation_release == "current"
        )  # previous default value, now needs to be provided explicitly

    def test_valid_with_release(self) -> None:
        p = SourceParamsEnsembl(species="mus_musculus", annotation_release="109")
        assert p.annotation_release == "109"

    def test_extra_field_forbidden(self) -> None:
        with pytest.raises(ValidationError):
            SourceParamsEnsembl(species="homo_sapiens", annotation_release="current", unknown="x")  # type: ignore[call-arg]


# ---------------------------------------------------------------------------
# SourceParamsNcbiAssembly
# ---------------------------------------------------------------------------


class TestSourceParamsNcbiAssembly:
    def test_valid(self) -> None:
        p = SourceParamsNcbiAssembly(
            refseq_assembly_accession="GCF_000001405.40",
            assembly_name="GRCh38.p14",
        )
        assert p.assembly_name == "GRCh38.p14"

    def test_invalid_accession_pattern(self) -> None:
        with pytest.raises(ValidationError):
            SourceParamsNcbiAssembly(
                refseq_assembly_accession="GCA_000001405.40",  # GCA not GCF
                assembly_name="GRCh38",
            )

    def test_invalid_accession_no_version(self) -> None:
        with pytest.raises(ValidationError):
            SourceParamsNcbiAssembly(
                refseq_assembly_accession="GCF_000001405",
                assembly_name="GRCh38",
            )

    def test_extra_field_forbidden(self) -> None:
        with pytest.raises(ValidationError):
            SourceParamsNcbiAssembly(
                refseq_assembly_accession="GCF_000001405.40", assembly_name="GRCh38.p14", unknown="x"  # type: ignore[call-arg]
            )


# ---------------------------------------------------------------------------
# SourceParamsNcbiSpecies — _check_supported_taxa_sources_and_annotation_release
# ---------------------------------------------------------------------------


class TestSourceParamsNcbiSpecies:
    def _make(self, **overrides: str) -> SourceParamsNcbiSpecies:
        data = {
            "taxon": "plant",
            "species": "arabidopsis_thaliana",
            "assembly_source": "auto",
            "annotation_release": "current",
        }
        data.update(overrides)
        return SourceParamsNcbiSpecies(**data)

    def test_accepts_valid_defaults(self) -> None:
        model = self._make()
        assert model.taxon == "plant"

    @pytest.mark.parametrize(  # type: ignore [untyped-decorator]
        ("taxon", "assembly_source"),
        [
            ("viral", "reference"),
            ("metagenomes", "reference"),
            ("bacteria", "annotation_releases"),
        ],
    )
    def test_rejects_unsupported_source_for_taxon(self, taxon: str, assembly_source: str) -> None:
        with pytest.raises(ConfigurationError):
            self._make(taxon=taxon, assembly_source=assembly_source)

    def test_allows_non_current_annotation_release_with_annotation_releases(self) -> None:
        model = self._make(
            taxon="plant",
            assembly_source="annotation_releases",
            annotation_release="109",
        )
        assert model.annotation_release == "109"

    def test_allows_non_current_annotation_release_with_auto_when_annotation_releases_available(self) -> None:
        model = self._make(
            taxon="plant",
            assembly_source="auto",
            annotation_release="109",
        )
        assert model.annotation_release == "109"

    @pytest.mark.parametrize(  # type: ignore[untyped-decorator]
        ("taxon", "assembly_source"),
        [
            ("bacteria", "latest_assembly_versions"),
            ("bacteria", "reference"),
            ("viral", "auto"),
        ],
    )
    def test_rejects_non_current_annotation_release_when_not_supported(
        # for simplicity, use str here instead of the Literals which will be checked by the pydantic model
        self,
        taxon: str,
        assembly_source: str,
    ) -> None:
        with pytest.raises(ConfigurationError):
            self._make(
                taxon=taxon,
                assembly_source=assembly_source,
                annotation_release="109",
            )

    def test_extra_field_forbidden(self) -> None:
        with pytest.raises(ValidationError):
            self._make(extra_field="x")


# ---------------------------------------------------------------------------
# SourceNcbi discriminated union — discriminator="mode"
# ---------------------------------------------------------------------------


class TestSourceNcbiUnion:
    def test_species_mode(self) -> None:
        obj = _SourceNcbiAdapter.validate_python(
            {
                "source": "ncbi",
                "mode": "species",
                "parameters": {
                    "taxon": "vertebrate_mammalian",
                    "species": "homo_sapiens",
                    "assembly_source": "auto",
                    "annotation_release": "current",
                },
            }
        )
        assert isinstance(obj, SourceNcbiSpecies)

    def test_assembly_mode(self) -> None:
        obj = _SourceNcbiAdapter.validate_python(
            {
                "source": "ncbi",
                "mode": "assembly",
                "parameters": {
                    "refseq_assembly_accession": "GCF_000001405.40",
                    "assembly_name": "GRCh38.p14",
                },
            }
        )
        assert isinstance(obj, SourceNcbiAssembly)

    def test_invalid_mode(self) -> None:
        with pytest.raises(ValidationError):
            _SourceNcbiAdapter.validate_python(
                {
                    "source": "ncbi",
                    "mode": "unknown",
                    "parameters": {},
                }
            )


# ---------------------------------------------------------------------------
# SourceConfigs discriminated union — discriminator="source"
# ---------------------------------------------------------------------------


class TestSourceConfigsUnion:
    def test_custom_source(self) -> None:
        obj = _SourceConfigsAdapter.validate_python(
            {
                "source": "custom",
                "parameters": {"file_annotation": "ann.gtf", "file_sequence": "seq.fa"},
            }
        )
        assert isinstance(obj, SourceCustom)

    def test_ensembl_source(self) -> None:
        obj = _SourceConfigsAdapter.validate_python(
            {"source": "ensembl", "parameters": {"species": "homo_sapiens"}}
        )
        assert isinstance(obj, SourceEnsembl)

    def test_ncbi_species_source(self) -> None:
        obj = _SourceConfigsAdapter.validate_python(
            {
                "source": "ncbi",
                "mode": "species",
                "parameters": {
                    "taxon": "vertebrate_mammalian",
                    "species": "homo_sapiens",
                    "assembly_source": "auto",
                    "annotation_release": "current",
                },
            }
        )
        assert isinstance(obj, SourceEnsembl) is False

    def test_ncbi_assembly_source(self) -> None:
        obj = _SourceConfigsAdapter.validate_python(
            {
                "source": "ncbi",
                "mode": "assembly",
                "parameters": {
                    "refseq_assembly_accession": "GCF_000001405.40",
                    "assembly_name": "GRCh38.p14",
                },
            }
        )
        assert obj.source == "ncbi"

    def test_invalid_source_discriminator(self) -> None:
        with pytest.raises(ValidationError):
            _SourceConfigsAdapter.validate_python({"source": "unknown_source"})
