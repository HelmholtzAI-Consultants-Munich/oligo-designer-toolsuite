############################################
# imports
############################################

import inspect
import logging
import os
import sys
from argparse import ArgumentParser, RawDescriptionHelpFormatter
from datetime import datetime
from typing import Any, Callable, TypeVar, cast

from Bio.SeqUtils import MeltingTemp as mt

from oligo_designer_toolsuite.database import OligoDatabase
from oligo_designer_toolsuite.utils import count_kmer_abundance

F = TypeVar("F", bound=Callable[..., Any])

############################################
# Utils functions
############################################


def setup_logging(
    dir_output: str,
    pipeline_name: str,
    log_level: int = logging.NOTSET,
    include_console: bool = False,
    log_start_message: bool = False,
) -> None:
    """
    Set up logging configuration for a pipeline.

    This function creates a consistent logging setup across all pipelines, creating a log file
    in the output directory and optionally writing to the console.

    :param dir_output: Directory path where output files will be saved.
    :type dir_output: str
    :param pipeline_name: Name of the pipeline (used in log file name).
    :type pipeline_name: str
    :param log_level: Logging level (default: logging.NOTSET).
    :type log_level: int
    :param include_console: Whether to also log to console (default: False).
    :type include_console: bool
    :param log_start_message: Whether to log a "START PIPELINE" message (default: False).
    :type log_start_message: bool
    """
    timestamp = datetime.now()
    file_logger = os.path.join(
        dir_output,
        f"log_{pipeline_name}_{timestamp.year}-{timestamp.month}-{timestamp.day}-{timestamp.hour}-{timestamp.minute}.txt",
    )

    handlers: list[logging.Handler] = [logging.FileHandler(file_logger)]
    if include_console:
        handlers.append(logging.StreamHandler())

    logging.basicConfig(
        format="%(asctime)s [%(levelname)s] %(message)s",
        level=log_level,
        handlers=handlers,
        force=True,  # Force reconfiguration if logging was already configured
    )
    logging.captureWarnings(True)

    if log_start_message:
        logging.info("--------------START PIPELINE--------------")


def base_parser() -> dict[str, Any]:
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
            logging.info("Parameter: %s = %s", key, value)


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

    logging.info("Function: %s", func.__name__)
    for name, value in bound_args.arguments.items():
        if name != "self":
            logging.info("Parameter: %s = %s", name, value)

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
    Decorator for basic pipeline steps that logs parameters and tracks database info.

    :param step_name: Name of the pipeline step.
    :type step_name: str
    :return: Decorator function.
    :rtype: Callable[[F], F]
    """

    def decorator(function: F) -> F:
        def wrapper(*args: Any, **kwargs: Any) -> Any:
            logging.info(f"Parameters {step_name}:")
            log_parameters_and_get_db(function, args, kwargs)

            oligo_database = function(*args, **kwargs)

            num_genes, num_oligos = get_oligo_database_info(oligo_database.database)
            logging.info(
                f"Step - {step_name}: database contains {num_oligos} oligos from {num_genes} regions."
            )

            return oligo_database

        return cast(F, wrapper)

    return decorator


def pipeline_step_advanced(step_name: str) -> Callable[[F], F]:
    """
    Decorator for advanced pipeline steps that logs parameters and tracks database changes.

    :param step_name: Name of the pipeline step.
    :type step_name: str
    :return: Decorator function.
    :rtype: Callable[[F], F]
    """

    def decorator(function: F) -> F:
        def wrapper(*args: Any, **kwargs: Any) -> Any:
            logging.info(f"Parameters {step_name}:")
            oligo_database = log_parameters_and_get_db(function, args, kwargs)

            num_genes_before, num_oligos_before = get_oligo_database_info(oligo_database.database)

            oligo_database, *returned_values = function(*args, **kwargs)

            num_genes_after, num_oligos_after = get_oligo_database_info(oligo_database.database)
            logging.info(
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
        logging.error("The oligo database is empty. Exiting program...")
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
    tm_parameters["nn_table"] = getattr(mt, tm_parameters["nn_table"])
    tm_parameters["tmm_table"] = getattr(mt, tm_parameters["tmm_table"])
    tm_parameters["imm_table"] = getattr(mt, tm_parameters["imm_table"])
    tm_parameters["de_table"] = getattr(mt, tm_parameters["de_table"])
    return tm_parameters


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
                logging.info(
                    f"K-mer {kmer} has abundance {abundance} which is greater than the threshold {kmer_abundance_threshold[k]}"
                )
                highly_abundant_kmer_sequences.append(kmer)

    return highly_abundant_kmer_sequences
