############################################
# imports
############################################

import logging
import inspect
from argparse import ArgumentParser, RawDescriptionHelpFormatter

from oligo_designer_toolsuite.database import OligoDatabase

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
        prog="Genomic Region generator",
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
    args = vars(args)
    return args


def log_parameters(parameters):
    """Log function parameters.

    :param parameters: Dict with parameter name : parameter value pairs
    :type parameters: dict
    """
    for key, value in parameters.items():
        if key != "self":
            logging.info(f"{key} = {value}")


def get_oligo_database_info(oligo_database: dict):
    """Count the number of oligos and genes in the database.

    :param oligo_database: Database with oligos.
    :type oligo_database: dict
    :return: Number of genes and oligos in the database.
    :rtype: int, int
    """
    genes = oligo_database.keys()
    num_genes = len(genes)
    num_oligos = 0
    for gene in genes:
        num_oligos += len(oligo_database[gene].keys())

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

def generation_step(step_name):
    def decorator(function):
        def wrapper(*args, **kwargs):
            ##### log parameters #####
            logging.info(f"Parameters {step_name}:")
            arguments, _, _, values = inspect.getargvalues(inspect.currentframe())
            parameters = {i: values[i] for i in arguments}
            parameters.update(kwargs)
            log_parameters(parameters)

            #### call the function
            oligo_database, *returned_values = function(*args, **kwargs)

            ##### loggig database information #####
            num_genes, num_oligos = get_oligo_database_info(oligo_database.database)
            logging.info(f"Step - Generate oligos: the database contains {num_oligos} oligos from {num_genes} genes.")

            return oligo_database, *returned_values
        return wrapper
    return decorator

def filtering_step(step_name):
    def decorator(function):
        def wrapper(oligo_database: OligoDatabase, *args, **kwargs):

            ##### log parameters #####
            logging.info(f"Parameters {step_name}:")
            arguments, _, _, values = inspect.getargvalues(inspect.currentframe())
            parameters = {i: values[i] for i in arguments}
            parameters.update(kwargs)
            log_parameters(parameters)
            num_genes_before, num_oligos_before = get_oligo_database_info(oligo_database.database)

            #### call the function
            oligo_database, *returned_values = function(oligo_database, *args, **kwargs)

            ##### loggig database information #####
            num_genes_after, num_oligos_after = get_oligo_database_info(oligo_database.database)
            logging.info(
                f"Step - Filter Oligos by {step_name}: the database contains {num_oligos_after} oligos from {num_genes_after} genes, while {num_oligos_before - num_oligos_after} oligos and {num_genes_before - num_genes_after} genes have been deleted in this step."
            )
            return oligo_database, *returned_values
        return wrapper
    return decorator
        