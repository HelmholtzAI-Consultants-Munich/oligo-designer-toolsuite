############################################
# imports
############################################

import os
import sys
import inspect
import logging

from argparse import ArgumentParser, RawDescriptionHelpFormatter

############################################
# Utils functions
############################################


def base_parser():
    """
    Parses command-line arguments for the genomic region generator.

    This function initializes a command-line parser and defines an argument for specifying a configuration file in YAML format.
    It processes the input arguments provided when the script is run from the command line, returning them in a dictionary.

    :return: A dictionary containing the command-line arguments where the key is the argument name and the value is its specified value.
    :rtype: dict
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


def base_log_parameters(parameters):
    """Log function parameters.

    :param parameters: Dict with parameter name : parameter value pairs
    :type parameters: dict
    """
    for key, value in parameters.items():
        if key != "self":
            logging.info(f"{key} = {value}")


def log_parameters_and_get_db(func, args, kwargs):
    """Log function parameters."""
    sig = inspect.signature(func)
    bound_args = sig.bind(*args, **kwargs)
    bound_args.apply_defaults()

    logging.info("Function: %s", func.__name__)
    for name, value in bound_args.arguments.items():
        if name != "self":
            logging.info("Parameter: %s = %s", name, value)

    return bound_args.arguments.get("oligo_database")


def get_oligo_database_info(oligo_database: dict):
    """Count the number of oligos and genes in the database.

    :param oligo_database: Database with oligos.
    :type oligo_database: dict
    :return: Number of genes and oligos in the database.
    :rtype: int, int
    """
    num_genes = len(oligo_database)
    num_oligos = sum(len(oligos) for oligos in oligo_database.values())
    return num_genes, num_oligos


def get_oligo_length_min_max_from_database(oligo_database: dict):
    """Get minimum and maximum length of oligos stored in the oligo database.

    :param oligo_database: Database with oligos.
    :type oligo_database: dict
    :return: Min and max length of oligos
    :rtype: int, int
    """
    oligo_length_min = sys.maxsize
    oligo_length_max = 0

    for region in oligo_database.keys():
        for oligo in oligo_database[region].keys():
            length = oligo_database[region][oligo]["length"]
            if length < oligo_length_min:
                oligo_length_min = length
            if length > oligo_length_max:
                oligo_length_max = length

    return oligo_length_min, oligo_length_max


def pipeline_step_basic(step_name: str):
    """Decorator function to log the input parameter of a general generative step (where an oligo database is created) of any pipeline and the information of the database generated.
    This decorator requires that the first returned value of the function is the oligo database.

    :param step_name: Name identifying the step.
    :type step_name: str.
    """

    def decorator(function):
        def wrapper(*args, **kwargs):
            logging.info(f"Parameters {step_name}:")
            log_parameters_and_get_db(function, args, kwargs)

            oligo_database, *returned_values = function(*args, **kwargs)

            num_genes, num_oligos = get_oligo_database_info(oligo_database.database)
            logging.info(f"Step - {step_name}: database contains {num_oligos} oligos from {num_genes} genes.")

            return oligo_database, *returned_values

        return wrapper

    return decorator


def pipeline_step_advanced(step_name: str):
    """Decorator function to log the input parameter of a general filtering step (where an oligo database is filtered) of any pipeline and the information of the changes applied to the database.
    This decorator requires that the first returned value of the function is the oligo database.

    Note: using this function can increase runtime when database is large.

    :param step_name: Name identifying the step.
    :type step_name: str.
    """

    def decorator(function):
        def wrapper(*args, **kwargs):
            logging.info(f"Parameters {step_name}:")
            oligo_database = log_parameters_and_get_db(function, args, kwargs)

            num_genes_before, num_oligos_before = get_oligo_database_info(oligo_database.database)

            oligo_database, *returned_values = function(*args, **kwargs)

            num_genes_after, num_oligos_after = get_oligo_database_info(oligo_database.database)
            logging.info(
                f"Step - {step_name}: database contains {num_oligos_after} oligos from {num_genes_after} genes, "
                f"{num_oligos_before - num_oligos_after} oligos and {num_genes_before - num_genes_after} genes removed."
            )

            return oligo_database, *returned_values

        return wrapper

    return decorator
