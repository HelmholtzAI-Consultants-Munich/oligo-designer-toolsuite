############################################
# imports
############################################

import inspect
import sys
from argparse import ArgumentParser, RawDescriptionHelpFormatter
from typing import Any, Callable, TypeVar, cast

import pandas as pd
from Bio.SeqUtils import MeltingTemp as mt

from oligo_designer_toolsuite._exceptions import FileFormatError
from oligo_designer_toolsuite.database import OligoDatabase
from oligo_designer_toolsuite.utils import check_if_dna_sequence, count_kmer_abundance, logger

F = TypeVar("F", bound=Callable[..., Any])

############################################
# Utils functions
############################################


def base_parser() -> dict[str, Any]:
    """
    Parse the command-line arguments shared by every pipeline entry point.

    Every ``<pipeline>_probe_designer`` CLI accepts a single ``-c / --config`` argument
    pointing at a YAML config file. This helper centralises that parsing so each pipeline's
    ``main()`` only needs a one-liner call.

    :return: Parsed arguments as a dictionary (currently a single key ``"config"``).
    :rtype: dict[str, Any]
    """
    parser = ArgumentParser(
        prog="Genomic Region Generator",
        usage="genomic_region_generation [options]",
        description=__doc__,
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


def base_log_parameters(parameters: dict[str, Any]) -> None:
    """
    Log all parameters from a dictionary, excluding 'self'.

    :param parameters: Dictionary of parameters to log.
    :type parameters: dict[str, Any]
    """
    for key, value in parameters.items():
        if key != "self":
            logger.info("Parameter: %s = %s", key, value)


def log_parameters_and_get_db(func: Callable[..., Any], args: tuple[Any, ...], kwargs: dict[str, Any]) -> Any:
    """
    Log function parameters and return the oligo_database argument if present.

    :param func: The function to inspect.
    :type func: Callable[..., Any]
    :param args: Positional arguments passed to the function.
    :type args: tuple[Any, ...]
    :param kwargs: Keyword arguments passed to the function.
    :type kwargs: dict[str, Any]
    :return: The oligo_database argument if present, otherwise None.
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
    Get information about the number of regions and oligos in a database.

    :param oligo_database: Dictionary containing region IDs as keys and oligo dictionaries as values.
    :type oligo_database: dict[str, dict[str, Any]]
    :return: Tuple containing (number of regions, total number of oligos).
    :rtype: tuple[int, int]
    """
    num_genes = len(oligo_database)
    num_oligos = sum(len(oligos) for oligos in oligo_database.values())
    return num_genes, num_oligos


def get_oligo_length_min_max_from_database(oligo_database: OligoDatabase) -> tuple[int, int]:
    """
    Get the minimum and maximum oligo lengths from the database.

    This function iterates through all oligos in the database to find the
    minimum and maximum length values.

    :param oligo_database: The OligoDatabase instance to query.
    :type oligo_database: OligoDatabase
    :return: A tuple containing (minimum_length, maximum_length).
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
    Decorator for basic pipeline steps that logs parameters and reports the resulting database size.

    Wraps a pipeline step method (typically a filter or a database-creation step) so that the
    parameters passed in are captured to the run log and the number of oligos / regions in the
    returned ``OligoDatabase`` is reported afterwards. Use when the wrapped step returns a single
    ``OligoDatabase``; for steps that also return auxiliary values use :func:`pipeline_step_advanced`.

    :param step_name: Human-readable name of the pipeline step (used in log messages).
    :type step_name: str
    :return: Decorator that instruments the wrapped function with logging.
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


def pipeline_step_advanced(step_name: str) -> Callable[[F], F]:
    """
    Decorator for pipeline steps that both mutate the database and return auxiliary values.

    Extends :func:`pipeline_step_basic` by also reporting **how many** oligos / regions were
    removed by the step (compared to the pre-call state). Use when the wrapped step returns
    a tuple ``(oligo_database, *extras)`` — the returned tuple is passed through unchanged
    after the log messages are emitted.

    :param step_name: Human-readable name of the pipeline step (used in log messages).
    :type step_name: str
    :return: Decorator that instruments the wrapped function with logging.
    :rtype: Callable[[F], F]
    """

    def decorator(function: F) -> F:
        """Return the instrumented replacement for ``function``."""

        def wrapper(*args: Any, **kwargs: Any) -> Any:
            """Run ``function``, logging parameters plus before/after database sizes."""
            logger.info(f"Parameters {step_name}:")
            oligo_database = log_parameters_and_get_db(function, args, kwargs)

            num_genes_before, num_oligos_before = get_oligo_database_info(oligo_database.database)

            oligo_database, *returned_values = function(*args, **kwargs)

            num_genes_after, num_oligos_after = get_oligo_database_info(oligo_database.database)
            logger.info(
                f"Step - {step_name}: database contains {num_oligos_after} oligos from {num_genes_after} regions, "
                f"{num_oligos_before - num_oligos_after} oligos and {num_genes_before - num_genes_after} regions removed."
            )

            return oligo_database, *returned_values

        return cast(F, wrapper)

    return decorator


def check_content_oligo_database(oligo_database: OligoDatabase) -> None:
    """
    Check if the oligo database is empty and exit if it is.

    :param oligo_database: The OligoDatabase instance to check.
    :type oligo_database: OligoDatabase
    :raises SystemExit: If the database is empty, exits with status code 1.
    """
    if len(oligo_database.get_regionid_list()) == 0:
        logger.error("The oligo database is empty. Exiting program...")
        print("The oligo database is empty. Exiting program...")
        sys.exit(1)  # Exit the program with a status code of 1


def format_sequence(database: OligoDatabase, property: str, region_id: str, oligo_id: str) -> str:
    """
    Get a sequence property as a string from the database, raising an error if not available.

    :param database: The OligoDatabase instance to query.
    :type database: OligoDatabase
    :param property: The property name to retrieve.
    :type property: str
    :param region_id: The region ID to query.
    :type region_id: str
    :param oligo_id: The oligo ID to query.
    :type oligo_id: str
    :return: The sequence as a string.
    :rtype: str
    :raises ValueError: If the property value is not a string.
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
    Preprocess melting temperature parameters by converting string table names to actual table objects.

    This function modifies the tm_parameters dictionary in place, converting the string names
    for nn_table, tmm_table, imm_table, and de_table to their corresponding table objects
    from Bio.SeqUtils.MeltingTemp.

    :param tm_parameters: Dictionary containing melting temperature parameters with string table names.
    :type tm_parameters: dict[str, Any]
    :return: The modified dictionary with table objects instead of string names.
    :rtype: dict[str, Any]
    """
    for key in ("nn_table", "tmm_table", "imm_table", "de_table"):
        if tm_parameters[key] is not None:
            tm_parameters[key] = getattr(mt, tm_parameters[key])
    return tm_parameters


def validate_codebook(
    codebook: pd.DataFrame,
    region_ids: list[str],
    *,
    source: str,
    expected_hamming_weight: int | None = None,
    index_name: str = "region_id",
) -> None:
    """
    Validate a codebook DataFrame against the expected format.

    The codebook is expected to have ``index_name`` as the index and one or more
    ``bit_*`` columns containing binary (0/1) values. All ``region_ids`` requested
    by the caller must be present in the index. If ``expected_hamming_weight`` is
    provided, each row must have exactly that many bits set; otherwise each row
    must have at least one bit set.

    :param codebook: Codebook with ``index_name`` as the index and ``bit_*`` columns.
    :type codebook: pd.DataFrame
    :param region_ids: Region IDs required to be present in the codebook index.
    :type region_ids: list[str]
    :param source: Source identifier (e.g. filename) used in error messages.
    :type source: str
    :param expected_hamming_weight: Required number of bits set per row, or ``None``
        to require at least one bit set per row.
    :type expected_hamming_weight: int | None
    :param index_name: Expected name of the codebook index column.
    :type index_name: str
    :raises FileFormatError: If any validation check fails.
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
    Validate a bit-indexed mapping table (e.g. readout probe table, initiator table), optionally
    against a codebook.

    A probe design run needs every codebook bit to have a corresponding readout / initiator
    sequence — without that, the assembled probes cannot be decoded at imaging time. This helper
    enforces that contract: every bit column present in ``codebook`` must appear as a row in
    the mapping table, and the sequence values in that row must be valid DNA. Bits present in
    the mapping table but not used by the codebook are silently accepted — a full-size readout
    set paired with a small-panel codebook is a legitimate configuration, not a mistake.

    When called with ``codebook=None`` (typically right after loading the table from a file, before
    a codebook has been obtained), the bit-coverage check is skipped and only the structural
    checks are enforced: non-empty table, ``bit`` index, no duplicate bits, required columns
    present, and sequence columns contain valid DNA. This lets a malformed file surface a clear
    :class:`FileFormatError` at load time instead of a cryptic downstream crash.

    :param table: Bit-indexed table (e.g. readout probe table, initiator table).
    :type table: pd.DataFrame
    :param codebook: Codebook whose ``bit_*`` columns drive the set of required bits. Pass ``None``
        to skip the bit-coverage check (useful for early load-time validation).
    :type codebook: pd.DataFrame | None
    :param source: Source identifier (e.g. filename) used in error messages.
    :type source: str
    :param required_columns: Columns that must be present in the table.
    :type required_columns: list[str]
    :param sequence_columns: Subset of ``required_columns`` whose values must be
        valid DNA sequences.
    :type sequence_columns: list[str]
    :raises FileFormatError: If any validation check fails.
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
    Validate that a primer sequence is a non-empty DNA string (A/C/G/T only).

    Thin wrapper around :func:`check_if_dna_sequence` that raises a
    ``FileFormatError`` with a ``source``-tagged message so all callers report
    errors consistently.

    :param sequence: The primer sequence to validate.
    :type sequence: str
    :param source: Identifier (e.g. ``"forward_primer"``) used in error messages.
    :type source: str
    :raises FileFormatError: If the value is not a non-empty DNA sequence.
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
    Get highly abundant k-mer sequences by identifying k-mers that exceed specified thresholds.

    This function counts k-mer abundances in FASTA files and identifies k-mers that exceed
    the specified abundance thresholds. These high-abundance k-mers are added to the list
    of highly abundant k-mer sequences, which can be used to filter out sequences during probe design.

    :param files_fasta: Path(s) to FASTA file(s) to analyze. Can be a single file path (str)
                        or a list of file paths (list[str]).
    :type files_fasta: str | list[str]
    :param kmer_abundance_threshold: Dictionary mapping k-mer length (int) to maximum allowed
                                     abundance threshold (float). K-mers exceeding this threshold
                                     will be added to highly abundant k-mer sequences.
    :type kmer_abundance_threshold: dict[int, float]
    :return: List of highly abundant k-mer sequences containing identified high-abundance k-mers.
    :rtype: list[str]
    """
    highly_abundant_kmer_sequences: list[str] = []

    # Count k-mer abundances
    kmer_abundance = count_kmer_abundance(
        files_fasta=files_fasta,
        k=list(kmer_abundance_threshold.keys()),
    )

    # Identify high-abundance k-mers and add them to highly abundant k-mer sequences
    for k, v in kmer_abundance.items():
        for kmer, abundance in v.items():
            if abundance > kmer_abundance_threshold[k]:
                logger.info(
                    f"K-mer {kmer} has abundance {abundance} which is greater than the threshold {kmer_abundance_threshold[k]}"
                )
                highly_abundant_kmer_sequences.append(kmer)

    return highly_abundant_kmer_sequences
