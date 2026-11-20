"""
Write the pydantic field descriptions as comments into the pipeline YAML configs.

Every key with a ``Field(description=...)`` gets its description as a ``#: `` comment on the line(s)
above it. The ``#: `` prefix marks the generated lines: they are removed and written again on each
run, everything else in the file (values, hand-written comments, blank lines) is left as it is.

Used as a pre-commit hook::

    python -m oligo_designer_toolsuite.config._yaml_descriptions data/configs
"""

############################################
# imports
############################################

import argparse
import re
import sys
import textwrap
import types
from collections.abc import Mapping
from pathlib import Path
from typing import Any, Union, get_args, get_origin

import yaml
from pydantic import BaseModel, ValidationError
from pydantic.fields import FieldInfo

from oligo_designer_toolsuite.config import (
    CycleHcrProbeDesignerConfig,
    HcrProbeDesignerConfig,
    MerfishProbeDesignerConfig,
    OligoSeqProbeDesignerConfig,
    ScrinshotProbeDesignerConfig,
    SeqfishPlusProbeDesignerConfig,
)

############################################
# constants
############################################

#: Pipeline model for each value of ``metadata.protocol``.
PIPELINE_MODELS: dict[str, type[BaseModel]] = {
    "Cycle-HCR": CycleHcrProbeDesignerConfig,
    "HCR": HcrProbeDesignerConfig,
    "MERFISH": MerfishProbeDesignerConfig,
    "Oligo-Seq": OligoSeqProbeDesignerConfig,
    "Scrinshot": ScrinshotProbeDesignerConfig,
    "SeqFish+": SeqfishPlusProbeDesignerConfig,
}

MARKER = "#: "
LINE_LENGTH = 110

_GENERATED_LINE = re.compile(r"^\s*#:( |$)")
_KEY_LINE = re.compile(r"^( *)([A-Za-z_]\w*):(.*)$")

DescriptionMap = dict[tuple[str, ...], str]

############################################
# collect descriptions from the model
############################################


def _find_field(model: type[BaseModel], key: str) -> tuple[str, FieldInfo] | None:
    """
    Find the field that a YAML key refers to, by field name or alias.

    :param model: Model that should contain the field.
    :type model: type[BaseModel]
    :param key: Key as written in the YAML file.
    :type key: str
    :return: Field name and field info, or None if the model has no such field.
    :rtype: tuple[str, FieldInfo] | None
    """
    for name, field in model.model_fields.items():
        if key in (name, field.alias, field.serialization_alias):
            return name, field
    return None


def _union_members(annotation: Any) -> list[type[BaseModel]]:
    """
    List the models of a union annotation, e.g. ``FilterEnabled | FilterDisabled``.

    :param annotation: Field annotation.
    :type annotation: Any
    :return: The models in the union, or an empty list if the annotation is not a union.
    :rtype: list[type[BaseModel]]
    """
    if get_origin(annotation) not in (Union, types.UnionType):
        return []
    return [arg for arg in get_args(annotation) if isinstance(arg, type) and issubclass(arg, BaseModel)]


def collect_descriptions(
    config: BaseModel,
    data: Mapping[str, Any],
    path: tuple[str, ...] = (),
    annotation: Any = None,
) -> DescriptionMap:
    """
    Map each key path in the YAML data to the description of its field.

    The validated config is walked alongside the raw data, so for a discriminated union the
    description of the member that is actually used is taken (e.g. coverage vs. minimum alignment
    length for ``hit_parameters.value``). A key that the active member ignores, like the threshold
    of a disabled filter, takes the description from another member of the union.

    :param config: Validated config (or nested part of it).
    :type config: BaseModel
    :param data: Raw YAML data belonging to ``config``.
    :type data: Mapping[str, Any]
    :param path: Key path of ``config`` within the YAML document.
    :type path: tuple[str, ...]
    :param annotation: Annotation of the field holding ``config``, used to find union siblings.
    :type annotation: Any
    :return: Description for each key path that has one.
    :rtype: DescriptionMap
    """
    descriptions: DescriptionMap = {}
    for key, value in data.items():
        # build the complete path so that every key is unique, e.g.
        # ("target_probes", "property_filters", "GC_content_filter")
        key_path = (*path, key)
        found = _find_field(type(config), key)

        if found is None:
            # key ignored by the active union member -> borrow the description from a sibling
            for member in _union_members(annotation):
                sibling = _find_field(member, key)
                # check FieldInfo (sibling[1]) for a description
                if sibling is not None and sibling[1].description:
                    descriptions[key_path] = sibling[1].description
                    break
            continue

        name, field = found
        if field.description:
            descriptions[key_path] = field.description
        # check if it is a nested model and if yes, find description for it
        child = getattr(config, name)
        if isinstance(child, BaseModel):
            descriptions.update(collect_descriptions(child, value, key_path, field.annotation))
    return descriptions


############################################
# insert descriptions into the YAML text
############################################


def _comment(text: str, indent: str, initial: str = "", subsequent: str = "") -> list[str]:
    """
    Render a description as generated comment lines, wrapped to the line length.

    :param text: Description text.
    :type text: str
    :param indent: Indentation of the key the comment belongs to.
    :type indent: str
    :param initial: Extra indentation of the first line after the marker.
    :type initial: str
    :param subsequent: Extra indentation of the following lines after the marker.
    :type subsequent: str
    :return: Comment lines including line breaks.
    :rtype: list[str]
    """
    lines = textwrap.wrap(
        " ".join(text.split()),
        width=LINE_LENGTH - len(indent) - len(MARKER),
        initial_indent=initial,
        subsequent_indent=subsequent,
        break_long_words=False,
        break_on_hyphens=False,
    )
    return [f"{indent}{MARKER}{line}\n" for line in lines]


def insert_descriptions(text: str, descriptions: DescriptionMap) -> str:
    """
    Write the descriptions as ``#: `` comments above their keys in a YAML document.

    Previously generated lines are removed first, so running this twice gives the same result.
    For an inline mapping (``key: {"a": 1, "b": 2}``) the descriptions of its keys are listed
    below the description of ``key``.

    :param text: YAML document.
    :type text: str
    :param descriptions: Description for each key path, see :func:`collect_descriptions`.
    :type descriptions: DescriptionMap
    :return: The YAML document with the generated comments.
    :rtype: str
    """
    output: list[str] = []
    parents: list[tuple[int, str]] = []  # (indentation, key) of the enclosing block keys

    for line in text.splitlines(keepends=True):
        if _GENERATED_LINE.match(line):
            continue
        match = _KEY_LINE.match(line)
        if match is None:
            output.append(line)
            continue

        indent, key, value = match.groups()
        # Find the full path of this key, e.g. ("target_probes", "property_filters", "Tm_filter").
        # Keys indented as much as or more than this one are not its parents, so remove them.
        # The remaining keys are its parents. Then remember this key as a possible parent of
        # the lines that follow.
        while parents and parents[-1][0] >= len(indent):
            parents.pop()
        path = (*(name for _, name in parents), key)
        parents.append((len(indent), key))

        if path in descriptions:
            output += _comment(descriptions[path], indent)
        # add comments to keys in blocks with "{}" above it
        if value.strip().startswith("{"):
            for child_path, description in descriptions.items():
                if child_path[:-1] == path:
                    output += _comment(f"{child_path[-1]}: {description}", indent, "  ", "    ")
        output.append(line)

    return "".join(output)


############################################
# files and command line
############################################


def add_descriptions_to_yaml(path: str | Path, model: type[BaseModel]) -> bool:
    """
    Add the descriptions of a pipeline model to a YAML config file in place.

    :param path: YAML config file.
    :type path: str | Path
    :param model: Pipeline model the file is validated against.
    :type model: type[BaseModel]
    :raises ValidationError: If the file does not validate against the model.
    :return: True if the file was changed.
    :rtype: bool
    """
    path = Path(path)
    text = path.read_text(encoding="utf-8")
    data = yaml.safe_load(text)
    config = model.model_validate(data)
    new_text = insert_descriptions(text, collect_descriptions(config, data))
    if new_text == text:
        return False
    path.write_text(new_text, encoding="utf-8")
    return True


def main(argv: list[str] | None = None) -> int:
    """
    Add the descriptions to the given YAML files, choosing the model by ``metadata.protocol``.

    Directories are expanded to the ``*.yaml`` files they contain. Files without
    ``metadata.protocol`` are skipped.

    :param argv: Files or directories, defaults to the command line arguments.
    :type argv: list[str] | None
    :return: 2 if a file has an unknown protocol or failed to validate, else 1 if a file was
        changed, otherwise 0.
    :rtype: int
    """
    parser = argparse.ArgumentParser(description="Add pydantic field descriptions to pipeline YAML configs.")
    parser.add_argument("paths", nargs="+", type=Path, help="YAML files or directories")
    args = parser.parse_args(argv)

    files = [file for p in args.paths for file in (sorted(p.glob("*.yaml")) if p.is_dir() else [p])]
    status = 0
    for file in files:
        data = yaml.safe_load(file.read_text(encoding="utf-8"))
        metadata = data.get("metadata") if isinstance(data, dict) else None
        protocol = metadata.get("protocol") if isinstance(metadata, dict) else None
        if protocol is None:
            continue
        if protocol not in PIPELINE_MODELS:
            print(f"{file}: unknown protocol {protocol!r}", file=sys.stderr)
            status = 2
            continue
        try:
            if add_descriptions_to_yaml(file, PIPELINE_MODELS[protocol]):
                print(f"Added descriptions to {file}")
                status = max(status, 1)  # keep 2 if an earlier file failed
        except ValidationError as error:
            print(
                f"{file}: does not validate against {PIPELINE_MODELS[protocol].__name__}\n{error}",
                file=sys.stderr,
            )
            status = 2
    return status


if __name__ == "__main__":
    sys.exit(main())
