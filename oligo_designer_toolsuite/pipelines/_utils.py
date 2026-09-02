"""
Shared helpers for the probe-design pipelines.

Covers CLI entry, logging, and config/table checks used across the assay-specific
pipeline modules in :mod:`oligo_designer_toolsuite.pipelines`.
"""

############################################
# imports
############################################

import inspect
import sys
from argparse import ArgumentParser, RawDescriptionHelpFormatter
from typing import Any, Callable, TypeVar, cast

import pandas as pd
import yaml
from Bio.SeqUtils import MeltingTemp as mt
from pydantic import BaseModel, ValidationError

from oligo_designer_toolsuite._exceptions import ConfigurationError, FileFormatError
from oligo_designer_toolsuite.config._completeness import check_config_complete
from oligo_designer_toolsuite.database import OligoDatabase
from oligo_designer_toolsuite.utils import check_if_dna_sequence, count_kmer_abundance, logger

F = TypeVar("F", bound=Callable[..., Any])
ConfigT = TypeVar("ConfigT", bound=BaseModel)

############################################
# Utils functions
############################################


def base_parser(*, prog: str, usage: str, description: str | None = None) -> dict[str, Any]:
    """
    Read the common command-line arguments used by pipeline entry points.

    Each probe-designer command accepts a ``-c / --config`` argument that points
    to a YAML configuration file. Pass pipeline-specific ``prog``, ``usage``, and
    ``description`` so ``--help`` matches the command that was invoked.

    :param prog: Program name shown in help (for example ``"MERFISH Probe Designer"``).
    :type prog: str
    :param usage: Usage line shown in help (for example
        ``"merfish_probe_designer [options]"``).
    :type usage: str
    :param description: Longer help text, usually the calling module's ``__doc__``.
    :type description: str | None
    :return: Parsed command-line arguments as a dictionary.
    :rtype: dict[str, Any]
    """
    parser = ArgumentParser(
        prog=prog,
        usage=usage,
        description=description,
        formatter_class=RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "-c",
        "--config",
        help="Path to the config file in yaml format, str",
        default=None,
        type=str,
        metavar="",
    )
    args = parser.parse_args()
    return vars(args)


def load_config(file_config: str, model: type[ConfigT]) -> ConfigT:
    """
    Read a configuration file and check it against the configuration model of a pipeline.

    The file has to set every parameter of the model. A parameter that is left out resolves to a
    different value depending on how much of its section is left out with it, so the parameters
    that are missing are reported together with the value that was used so far.

    :param file_config: Path to the config file in yaml format.
    :type file_config: str
    :param model: Configuration model of the pipeline.
    :type model: type[ConfigT]
    :return: The validated configuration.
    :rtype: ConfigT
    :raises SystemExit: If the file does not match the model or leaves parameters unset.
    """
    with open(file_config, "r") as handle:
        config_raw = yaml.safe_load(handle)

    # Reported and exited rather than raised: both messages are written for the person running the
    # command, and a traceback would print the whole report a second time.
    try:
        config_validated = model.model_validate(config_raw)
    except ValidationError as error:
        sys.exit(f"Invalid configuration file:\n{error}")

    try:
        check_config_complete(config_validated, config_raw, source=file_config)
    except ConfigurationError as error:
        sys.exit(f"Incomplete configuration file:\n{error}")

    return config_validated


def base_log_parameters(parameters: dict[str, Any]) -> None:
    """
    Write pipeline parameters to the log.

    This is used at the start of a run to record the configuration that was
    actually used. The ``self`` entry is ignored, which makes the function safe
    to use with parameters collected from class methods.

    :param parameters: Parameter names and values to log.
    :type parameters: dict[str, Any]
    """
    for key, value in parameters.items():
        if key != "self":
            logger.info("Parameter: %s = %s", key, value)


def log_parameters_and_get_db(func: Callable[..., Any], args: tuple[Any, ...], kwargs: dict[str, Any]) -> Any:
    """
    Log the parameters passed to a pipeline step.

    The function combines positional arguments, keyword arguments, and default
    values into one complete view of the call. If the step receives an
    ``oligo_database`` argument, that database is returned so decorators can
    report how the step changed the number of regions and oligos.

    :param func: Function whose call should be logged.
    :type func: Callable[..., Any]
    :param args: Positional arguments passed to the function.
    :type args: tuple[Any, ...]
    :param kwargs: Keyword arguments passed to the function.
    :type kwargs: dict[str, Any]
    :return: The ``oligo_database`` argument if present, otherwise ``None``.
    :rtype: Any
    """
    sig = inspect.signature(func)
    bound_args = sig.bind(*args, **kwargs)
    bound_args.apply_defaults()

    logger.info("Function: %s", func.__name__)
    for name, value in bound_args.arguments.items():
        if name != "self":
            logger.info("Parameter: %s = %s", name, value)

    return bound_args.arguments.get("oligo_database")


def get_oligo_database_info(oligo_database: dict[str, dict[str, Any]]) -> tuple[int, int]:
    """
    Count how many regions and oligos are stored in a database.

    This is mainly used for logging after filtering steps. It gives a compact
    summary of how many regions are still represented and how many candidate
    oligos remain in total.

    :param oligo_database: Raw database dictionary with regions as top-level keys.
    :type oligo_database: dict[str, dict[str, Any]]
    :return: Number of regions and total number of oligos.
    :rtype: tuple[int, int]
    """
    num_genes = len(oligo_database)
    num_oligos = sum(len(oligos) for oligos in oligo_database.values())
    return num_genes, num_oligos


def get_oligo_length_min_max_from_database(oligo_database: OligoDatabase) -> tuple[int, int]:
    """
    Find the shortest and longest oligo in the database.

    Some downstream steps need to know the actual length range of the remaining
    oligos. This helper scans the current database and returns the minimum and
    maximum sequence length.

    :param oligo_database: Oligo database to inspect.
    :type oligo_database: OligoDatabase
    :return: Minimum and maximum oligo length in bases.
    :rtype: tuple[int, int]
    """
    oligo_length_min = sys.maxsize
    oligo_length_max = 0

    region_ids = oligo_database.database.keys()

    for region_id in region_ids:
        oligo_ids = oligo_database.database[region_id].keys()
        for oligo_id in oligo_ids:
            length = oligo_database.database[region_id][oligo_id]["length"]
            if length < oligo_length_min:
                oligo_length_min = length
            if length > oligo_length_max:
                oligo_length_max = length

    return oligo_length_min, oligo_length_max


def pipeline_step_basic(step_name: str) -> Callable[[F], F]:
    """
    Add standard logging around a pipeline step.

    Use this decorator for steps that return only an :class:`OligoDatabase`.
    It logs the step parameters before the function runs and reports how many
    regions and oligos are present afterwards.

    :param step_name: Name shown in the log for this pipeline step.
    :type step_name: str
    :return: Decorator for the wrapped pipeline step.
    :rtype: Callable[[F], F]
    """

    def decorator(function: F) -> F:
        """Return the instrumented replacement for ``function``."""

        def wrapper(*args: Any, **kwargs: Any) -> Any:
            """Run ``function``, logging its parameters and the resulting database size."""
            logger.info(f"Parameters {step_name}:")
            log_parameters_and_get_db(function, args, kwargs)

            oligo_database = function(*args, **kwargs)

            num_genes, num_oligos = get_oligo_database_info(oligo_database.database)
            logger.info(
                f"Step - {step_name}: database contains {num_oligos} oligos from {num_genes} regions."
            )

            return oligo_database

        return cast(F, wrapper)

    return decorator


def check_content_oligo_database(oligo_database: OligoDatabase) -> None:
    """
    Stop the pipeline if no candidate oligos are left.

    Filtering can remove all regions from the database. When that happens, the
    pipeline cannot produce useful output. This helper stops the run early and
    writes a clear message instead of letting a later step fail with a less
    helpful error.

    :param oligo_database: Oligo database to check.
    :type oligo_database: OligoDatabase
    :raises SystemExit: If the database contains no regions.
    """
    if len(oligo_database.get_regionid_list()) == 0:
        logger.error("The oligo database is empty. Exiting program...")
        print("The oligo database is empty. Exiting program...")
        sys.exit(1)


def format_sequence(database: OligoDatabase, property: str, region_id: str, oligo_id: str) -> str:
    """
    Return a sequence property as a plain string.

    This helper is used when a database property is expected to be a single DNA
    sequence. It retrieves the value and checks that it is a string before the
    sequence is used by downstream code.

    :param database: Oligo database to query.
    :type database: OligoDatabase
    :param property: Name of the sequence property to retrieve.
    :type property: str
    :param region_id: Region that contains the oligo.
    :type region_id: str
    :param oligo_id: Oligo whose property should be returned.
    :type oligo_id: str
    :return: Sequence value as a string.
    :rtype: str
    :raises ValueError: If the retrieved value is not a string.
    """
    value = database.get_oligo_property_value(
        property=property,
        region_id=region_id,
        oligo_id=oligo_id,
        flatten=True,
    )
    if not isinstance(value, str):
        raise ValueError(f"Expected string for {property}, got {type(value)}")
    return value


def preprocess_tm_parameters(tm_parameters: dict[str, Any]) -> dict[str, Any]:
    """
    Prepare melting-temperature settings from a config file.

    Config files store melting-temperature table names as readable strings, such
    as ``"DNA_NN4"``. The Biopython Tm functions expect the matching table
    objects instead. This helper replaces the configured table names with the
    objects used during calculation.

    :param tm_parameters: Melting-temperature parameter dictionary.
    :type tm_parameters: dict[str, Any]
    :return: Updated parameter dictionary.
    :rtype: dict[str, Any]
    """
    for key in ("nn_table", "tmm_table", "imm_table", "de_table"):
        if tm_parameters[key] is not None:
            tm_parameters[key] = getattr(mt, tm_parameters[key])
    return tm_parameters


def apply_required_parameters(config: dict[str, Any]) -> None:
    """
    Copy the entries of ``required_parameters`` into the sections that consume them.

    Config files collect the gene list and the FASTA files that every run needs in a
    single ``required_parameters`` section for better user experience. The helper writes
    them to the entries the pipelines expect. This only covers the cases that are the same
    across all pipelines; some pipelines might need additional copies.

    :param config: Configuration dictionary, updated in place.
    :type config: dict[str, Any]
    """
    required_parameters = config["required_parameters"]
    oligo_generation = config["target_probes"]["oligo_generation"]
    oligo_generation["file_region_ids"] = required_parameters["targets"]
    oligo_generation["files_fasta_probe_database"] = required_parameters["target_genome"]
    config["target_probes"]["specificity_filters"]["specificity_blastn_filter"][
        "files_fasta_reference_database"
    ] = required_parameters["reference_genome"]


def validate_codebook(
    codebook: pd.DataFrame,
    region_ids: list[str],
    *,
    source: str,
    expected_hamming_weight: int | None = None,
    index_name: str = "region_id",
) -> None:
    """
    Check that a codebook has the expected format.

    A codebook assigns each region to a binary barcode. Rows represent regions,
    columns represent bits, and each value must be ``0`` or ``1``. This helper
    checks that the table is complete, has no duplicated regions or bit columns,
    contains only valid binary values, and covers all requested regions.

    If an expected Hamming weight is given, every row must contain exactly that
    number of active bits. Otherwise, each row must contain at least one active
    bit.

    :param codebook: Codebook table to validate.
    :type codebook: pd.DataFrame
    :param region_ids: Region IDs that must be present in the codebook.
    :type region_ids: list[str]
    :param source: Name or path used to identify the codebook in error messages.
    :type source: str
    :param expected_hamming_weight: Required number of active bits per row, or
        ``None`` to allow variable-weight codes.
    :type expected_hamming_weight: int | None
    :param index_name: Expected name of the codebook index.
    :type index_name: str
    :raises FileFormatError: If the codebook is incomplete or incorrectly formatted.
    """
    if codebook.index.name != index_name:
        raise FileFormatError(
            f"Codebook '{source}' must use '{index_name}' as the index, got '{codebook.index.name}'."
        )

    duplicate_regions = codebook.index[codebook.index.duplicated()].unique().tolist()
    if duplicate_regions:
        raise FileFormatError(
            f"Codebook '{source}' contains duplicate {index_name}s: {sorted(duplicate_regions)}."
        )

    if len(codebook.columns) == 0:
        raise FileFormatError(f"Codebook '{source}' must contain at least one bit column.")

    non_bit_columns = [c for c in codebook.columns if not str(c).startswith("bit_")]
    if non_bit_columns:
        raise FileFormatError(
            f"Codebook '{source}' must have all columns named with the 'bit_*' pattern. "
            f"Found columns that don't match: {non_bit_columns}."
        )

    duplicate_cols = codebook.columns[codebook.columns.duplicated()].unique().tolist()
    if duplicate_cols:
        raise FileFormatError(
            f"Codebook '{source}' contains duplicate bit columns: {sorted(duplicate_cols)}."
        )

    if codebook.isna().any().any():
        rows_with_nan = codebook.index[codebook.isna().any(axis=1)].unique().tolist()
        raise FileFormatError(f"Codebook '{source}' contains NaN values in rows: {sorted(rows_with_nan)}.")

    if not codebook.isin([0, 1]).all().all():
        invalid_rows = codebook.index[~codebook.isin([0, 1]).all(axis=1)].unique().tolist()
        raise FileFormatError(
            f"Codebook '{source}' must contain only 0/1 values. "
            f"Rows with invalid values: {sorted(invalid_rows)}."
        )

    row_sums = (codebook == 1).sum(axis=1)
    if expected_hamming_weight is None:
        invalid_rows = codebook.index[row_sums < 1].unique().tolist()
        if invalid_rows:
            raise FileFormatError(
                f"Codebook '{source}' must have at least one bit set per row. "
                f"Rows with no bits set: {sorted(invalid_rows)}."
            )
    else:
        invalid_rows = codebook.index[row_sums != expected_hamming_weight].unique().tolist()
        if invalid_rows:
            raise FileFormatError(
                f"Codebook '{source}' must have exactly {expected_hamming_weight} bit(s) set per row. "
                f"Rows with wrong count: {sorted(invalid_rows)}."
            )

    missing_region_ids = set(region_ids) - set(codebook.index)
    if missing_region_ids:
        raise FileFormatError(f"Codebook '{source}' is missing {index_name}s: {sorted(missing_region_ids)}.")


def validate_bit_mapping_table(
    table: pd.DataFrame,
    codebook: pd.DataFrame | None = None,
    *,
    source: str,
    required_columns: list[str],
    sequence_columns: list[str],
) -> None:
    """
    Check a table that maps codebook bits to probe sequences.

    These tables connect each codebook bit to a sequence, such as a readout
    probe or initiator. The table must have one row per bit, contain the required
    columns, and provide valid DNA sequences in the selected sequence columns.

    When a codebook is provided, this helper also checks that every bit used in
    the codebook has a matching row in the table. Extra rows are allowed, which
    makes it possible to use a larger readout set with a smaller codebook.

    :param table: Bit mapping table to validate.
    :type table: pd.DataFrame
    :param codebook: Optional codebook used to check bit coverage.
    :type codebook: pd.DataFrame | None
    :param source: Name or path used to identify the table in error messages.
    :type source: str
    :param required_columns: Columns that must be present in the table.
    :type required_columns: list[str]
    :param sequence_columns: Columns that must contain valid DNA sequences.
    :type sequence_columns: list[str]
    :raises FileFormatError: If the table is incomplete or incorrectly formatted.
    """
    if len(table) == 0:
        raise FileFormatError(f"Table '{source}' is empty. Expected at least one bit-indexed row.")

    if table.index.name != "bit":
        raise FileFormatError(f"Table '{source}' must use 'bit' as the index, got '{table.index.name}'.")

    duplicate_bits = table.index[table.index.duplicated()].unique().tolist()
    if duplicate_bits:
        raise FileFormatError(f"Table '{source}' contains duplicate bit entries: {sorted(duplicate_bits)}.")

    missing_columns = set(required_columns) - set(table.columns)
    if missing_columns:
        raise FileFormatError(
            f"Table '{source}' is missing required columns: {sorted(missing_columns)}. "
            f"Required columns are: {required_columns}."
        )

    if codebook is not None:
        required_bits = set(codebook.columns)
        table_bits = set(table.index)
        missing_bits = required_bits - table_bits
        if missing_bits:
            raise FileFormatError(
                f"Table '{source}' is missing entries for codebook bits: {sorted(missing_bits)}."
            )

    for column in sequence_columns:
        invalid_bits = [
            bit
            for bit, value in table[column].items()
            if not isinstance(value, str) or not check_if_dna_sequence(value)
        ]
        if invalid_bits:
            raise FileFormatError(
                f"Table '{source}' column '{column}' must contain non-empty DNA sequences "
                f"(A/C/G/T only). Invalid entries at bits: {sorted(invalid_bits)}."
            )


def validate_primer_sequence(sequence: str, *, source: str) -> None:
    """
    Check that a primer is written as a valid DNA sequence.

    A valid primer must be a non-empty string containing only ``A``, ``C``,
    ``G``, and ``T``. The ``source`` label is included in error messages so it
    is clear which primer caused the problem.

    :param sequence: Primer sequence to check.
    :type sequence: str
    :param source: Primer label used in error messages.
    :type source: str
    :raises FileFormatError: If the primer is empty or contains invalid characters.
    """
    if not isinstance(sequence, str) or not check_if_dna_sequence(sequence):
        raise FileFormatError(
            f"Primer '{source}' must be a non-empty DNA sequence (A/C/G/T only), got {sequence!r}."
        )


def get_highly_abundant_kmer_sequences(
    files_fasta: str | list[str],
    kmer_abundance_threshold: dict[int, float],
) -> list[str]:
    """
    Find k-mers that occur too often in the reference sequences.

    Very common k-mers can make probes less specific because they appear in many
    places. This helper counts k-mers in one or more FASTA files and returns the
    sequences whose abundance is above the chosen threshold for their length.

    :param files_fasta: FASTA file path or list of FASTA file paths to scan.
    :type files_fasta: str | list[str]
    :param kmer_abundance_threshold: Maximum allowed abundance for each k-mer length.
    :type kmer_abundance_threshold: dict[int, float]
    :return: K-mer sequences that exceed their abundance threshold.
    :rtype: list[str]
    """
    highly_abundant_kmer_sequences: list[str] = []

    kmer_abundance = count_kmer_abundance(
        files_fasta=files_fasta,
        k=list(kmer_abundance_threshold.keys()),
    )

    for k, v in kmer_abundance.items():
        for kmer, abundance in v.items():
            if abundance > kmer_abundance_threshold[k]:
                logger.info(
                    f"K-mer {kmer} has abundance {abundance} which is greater than the threshold {kmer_abundance_threshold[k]}"
                )
                highly_abundant_kmer_sequences.append(kmer)

    return highly_abundant_kmer_sequences
