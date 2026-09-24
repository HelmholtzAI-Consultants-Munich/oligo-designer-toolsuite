"""
Completeness check for pipeline configuration files.

A parameter that is left out of a configuration file does not have one obvious value: which
default applies depends on how much of the surrounding section was left out with it. The command
line pipelines therefore require every parameter to be written out.

Only the command line entry points call this, so validating a configuration elsewhere, as
ODT-Cloud does, behaves as it did before.
"""

############################################
# imports
############################################

from collections.abc import Mapping, Sequence
from typing import Any, NamedTuple

import yaml
from pydantic import BaseModel
from pydantic.aliases import AliasChoices
from pydantic.fields import FieldInfo

from oligo_designer_toolsuite._exceptions import ConfigurationError

############################################
# Completeness check
############################################

ALLOW_INCOMPLETE_MARKER = "x-allow-incomplete"


class MissingConfigKey(NamedTuple):
    """
    A parameter that the configuration file does not set.

    :param path: Dotted path of the parameter, e.g. ``target_probes.global_parameters``.
    :type path: str
    :param key: Key to write in the configuration file, which is the alias where a field has one.
    :type key: str
    :param value: Value the model fell back to.
    :type value: Any
    """

    path: str
    key: str
    value: Any


def _accepted_keys(name: str, field: FieldInfo) -> list[str]:
    """
    Collect every key that validation accepts for one field.

    Kept for models whose fields carry aliases. The only one today, ``BlastnSearchParameters``, is
    exempt from the check, so this currently resolves to the field name for every parameter walked.

    :param name: Name the field is declared under.
    :type name: str
    :param field: Field to collect the keys of.
    :type field: FieldInfo
    :return: The accepted keys, most specific first.
    :rtype: list[str]
    """
    keys: list[str] = []
    validation_alias = field.validation_alias
    if isinstance(validation_alias, AliasChoices):
        keys.extend(choice for choice in validation_alias.choices if isinstance(choice, str))
    elif isinstance(validation_alias, str):
        keys.append(validation_alias)
    if field.alias is not None:
        keys.append(field.alias)
    keys.append(name)
    return list(dict.fromkeys(keys))


def _preferred_key(name: str, field: FieldInfo) -> str:
    """
    Return the key to suggest when a field has to be added to a configuration file.

    :param name: Name the field is declared under.
    :type name: str
    :param field: Field to suggest a key for.
    :type field: FieldInfo
    :return: The key, which is the serialization alias where a field has one.
    :rtype: str
    """
    if isinstance(field.serialization_alias, str):
        return field.serialization_alias
    return _accepted_keys(name, field)[0]


def _allows_incomplete(model: type[BaseModel]) -> bool:
    """
    Report whether a model is allowed to leave parameters unset.

    Such a model passes its fields on to an external tool and drops the unset ones when it is
    serialized, so an omitted key means "do not pass this option" rather than "use a default".

    :param model: Model to check.
    :type model: type[BaseModel]
    :return: True if the model is exempt from the completeness check.
    :rtype: bool
    """
    json_schema_extra = model.model_config.get("json_schema_extra")
    return isinstance(json_schema_extra, dict) and bool(json_schema_extra.get(ALLOW_INCOMPLETE_MARKER))


def find_missing_config_keys(
    config: BaseModel,
    config_raw: Mapping[str, Any],
    path: str = "",
    missing: list[MissingConfigKey] | None = None,
) -> list[MissingConfigKey]:
    """
    Collect the parameters that a configuration file does not set.

    The validated configuration is walked in parallel with the mapping it was built from, which
    means the branch of a union does not have to be worked out from the annotations: it is the type
    of the attribute that validation produced. A key that is absent is reported and not descended
    into, so a missing section is one entry rather than one entry per parameter below it.

    :param config: Validated configuration, or one of the models nested in it.
    :type config: BaseModel
    :param config_raw: Mapping that ``config`` was validated from.
    :type config_raw: Mapping[str, Any]
    :param path: Dotted path of ``config`` within the whole configuration.
    :type path: str
    :param missing: Parameters found so far, used by the recursion.
    :type missing: list[MissingConfigKey] | None
    :return: The missing parameters, in the order the models declare them.
    :rtype: list[MissingConfigKey]
    """
    missing = [] if missing is None else missing
    model = type(config)
    if _allows_incomplete(model) or not isinstance(config_raw, Mapping):
        return missing

    for name, field in model.model_fields.items():
        field_path = f"{path}.{name}" if path else name
        key = next((key for key in _accepted_keys(name, field) if key in config_raw), None)
        if key is None:
            missing.append(MissingConfigKey(field_path, _preferred_key(name, field), getattr(config, name)))
            continue

        value = getattr(config, name)
        value_raw = config_raw[key]
        if isinstance(value, BaseModel):
            find_missing_config_keys(value, value_raw, field_path, missing)
        elif isinstance(value, Sequence) and not isinstance(value, (str, bytes)):
            if isinstance(value_raw, Sequence) and not isinstance(value_raw, (str, bytes)):
                for index, (item, item_raw) in enumerate(zip(value, value_raw)):
                    if isinstance(item, BaseModel):
                        find_missing_config_keys(item, item_raw, f"{field_path}[{index}]", missing)

    return missing


def format_missing_config_keys(missing: list[MissingConfigKey], source: str | None = None) -> str:
    """
    List the missing parameters and the value each of them fell back to.

    :param missing: Parameters as returned by :func:`find_missing_config_keys`.
    :type missing: list[MissingConfigKey]
    :param source: Path of the configuration file, named in the message.
    :type source: str | None
    :return: The message.
    :rtype: str
    """
    lines = [
        f"All parameters have to be set in {source or 'the configuration file'},"
        f" {len(missing)} are missing.",
        "",
        "The value each of them fell back to:",
        "",
    ]
    for entry in missing:
        value = entry.value.model_dump(mode="json") if isinstance(entry.value, BaseModel) else entry.value
        # only the first line: dumping a scalar on its own also emits YAML's "..." end marker
        rendered = yaml.safe_dump(value, default_flow_style=True).splitlines()[0]
        lines.append(f"  {entry.path}: {rendered}")
    return "\n".join(lines)


def check_config_complete(
    config: BaseModel, config_raw: Mapping[str, Any], source: str | None = None
) -> None:
    """
    Raise if a configuration file leaves any parameter unset.

    :param config: Configuration validated from ``config_raw``.
    :type config: BaseModel
    :param config_raw: Mapping loaded from the configuration file.
    :type config_raw: Mapping[str, Any]
    :param source: Path of the configuration file, named in the error message.
    :type source: str | None
    :raises ConfigurationError: If at least one parameter is not set.
    """
    missing = find_missing_config_keys(config, config_raw)
    if missing:
        raise ConfigurationError(format_missing_config_keys(missing, source))
