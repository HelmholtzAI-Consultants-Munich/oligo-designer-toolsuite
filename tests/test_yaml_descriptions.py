from pathlib import Path

import yaml
from pydantic import BaseModel

from oligo_designer_toolsuite.config._yaml_descriptions import (
    PIPELINE_MODELS,
    add_descriptions_to_yaml,
    main,
)

CONFIGS = Path("data/configs")
PIPELINE_CONFIGS = [
    "cycle_hcr_probe_designer_gandin.yaml",
    "hcr_probe_designer.yaml",
    "merfish_probe_designer.yaml",
    "oligo_seq_probe_designer.yaml",
    "scrinshot_probe_designer.yaml",
    "seqfish_plus_probe_designer.yaml",
]


def _copy(name: str, tmp_path: Path) -> Path:
    """Copy a shipped config without its generated comments (starting with #:), so each
    test starts from a clean file."""
    target = tmp_path / name
    lines = (CONFIGS / name).read_text().splitlines(keepends=True)
    target.write_text("".join(line for line in lines if not line.lstrip().startswith("#: ")))
    return target


def _model(path: Path) -> type[BaseModel]:
    return PIPELINE_MODELS[yaml.safe_load(path.read_text())["metadata"]["protocol"]]


def _comment_above(text: str, line: str) -> list[str]:
    """Return the generated comment lines directly above the first occurrence of ``line``."""
    lines = text.splitlines()
    index = lines.index(line)
    comments: list[str] = []
    while lines[index - 1].lstrip().startswith("#: "):
        index -= 1
        comments.insert(0, lines[index].strip())
    return comments


def test_only_adds_comments_and_is_idempotent(tmp_path: Path) -> None:
    for name in PIPELINE_CONFIGS:
        path = _copy(name, tmp_path)
        original = path.read_text()

        assert add_descriptions_to_yaml(path, _model(path)) is True, name
        first = path.read_text()
        assert add_descriptions_to_yaml(path, _model(path)) is False, name
        assert path.read_text() == first, name

        assert yaml.safe_load(first) == yaml.safe_load(original), name
        kept = [line for line in first.splitlines() if not line.lstrip().startswith("#: ")]
        assert kept == original.splitlines(), name
        generated = [line for line in first.splitlines() if line.lstrip().startswith("#: ")]
        assert all(len(line) <= 110 for line in generated), name


def test_hand_written_comments_are_kept(tmp_path: Path) -> None:
    path = _copy("cycle_hcr_probe_designer_gandin.yaml", tmp_path)
    add_descriptions_to_yaml(path, _model(path))
    text = path.read_text()

    assert "      # --- set by Gandin et al. ---\n      #: " in text
    assert "### PROBE DESIGN PARAMETERS ###" in text


def test_stale_description_is_replaced(tmp_path: Path) -> None:
    path = _copy("hcr_probe_designer.yaml", tmp_path)
    path.write_text(path.read_text().replace("  n_jobs: 4", "  #: outdated text\n  n_jobs: 4"))

    add_descriptions_to_yaml(path, _model(path))
    comment = _comment_above(path.read_text(), "  n_jobs: 4")

    assert "#: outdated text" not in comment
    assert comment[0].startswith("#: number of cores")


def test_description_of_active_union_member(tmp_path: Path) -> None:
    path = _copy("cycle_hcr_probe_designer_gandin.yaml", tmp_path)
    add_descriptions_to_yaml(path, _model(path))
    text = path.read_text()

    assert _comment_above(text, "        value: 90") == [
        "#: Coverage in %, alternatively, min_alignment_length can be used"
    ]
    assert _comment_above(text, "        value: 16") == [
        "#: Number of nucleotides for alignment, alternatively, coverage can be used"
    ]


def test_key_ignored_by_active_member_borrows_description(tmp_path: Path) -> None:
    path = _copy("cycle_hcr_probe_designer_gandin.yaml", tmp_path)
    add_descriptions_to_yaml(path, _model(path))
    text = path.read_text()

    # isoform_consensus under a disabled filter, file under codebook source "generate"
    assert _comment_above(text, "      isoform_consensus: 0")[0].startswith("#: Minimum isoform consensus")
    assert _comment_above(text, "    file:")[0].startswith("#: Path to the codebook file")


def test_inline_mapping_lists_its_keys(tmp_path: Path) -> None:
    path = _copy("merfish_probe_designer.yaml", tmp_path)
    add_descriptions_to_yaml(path, _model(path))
    line = next(line for line in path.read_text().splitlines() if line.startswith("      search_parameters:"))
    comment = _comment_above(path.read_text(), line)

    assert comment[0] == "#: BLAST options for the specificity search."
    listed = [entry.split(":")[1].strip() for entry in comment[1:] if entry.startswith("#:   -")]
    assert listed == list(yaml.safe_load(line.strip())["search_parameters"])


def test_main_skips_configs_without_protocol_and_reports_changes(tmp_path: Path) -> None:
    region_config = _copy("genomic_region_generator_ncbi.yaml", tmp_path)
    original = region_config.read_text()
    _copy("hcr_probe_designer.yaml", tmp_path)

    assert main([str(tmp_path)]) == 1
    assert region_config.read_text() == original
    assert main([str(tmp_path)]) == 0


def test_main_fails_on_invalid_config_without_changing_it(tmp_path: Path) -> None:
    path = _copy("hcr_probe_designer.yaml", tmp_path)
    path.write_text(path.read_text().replace("n_jobs: 4", "n_jobs: -4"))
    original = path.read_text()
    valid = _copy("merfish_probe_designer.yaml", tmp_path)

    # the error status is kept even though the valid file processed afterwards is changed
    assert main([str(path), str(valid)]) == 2
    assert path.read_text() == original


def test_shipped_configs_are_up_to_date(tmp_path: Path) -> None:
    """Fails if a model description changed without re-running the pre-commit hook."""
    for name in PIPELINE_CONFIGS:
        path = tmp_path / name
        path.write_text((CONFIGS / name).read_text())
        assert add_descriptions_to_yaml(path, _model(path)) is False, name
