############################################
# imports
############################################

import logging
import os
import time
from argparse import ArgumentParser, RawDescriptionHelpFormatter
from datetime import datetime
from pathlib import Path

import yaml

timestamp = datetime.now()
file_logger = f"log_padlock_probe_designer_{timestamp.year}-{timestamp.month}-{timestamp.day}-{timestamp.hour}-{timestamp.minute}.txt"
logging.getLogger("padlock_probe_designer")
logging.basicConfig(
    format="%(asctime)s [%(levelname)s] %(message)s",
    level=logging.NOTSET,
    handlers=[logging.FileHandler(file_logger), logging.StreamHandler()],
)


def generate_config_file(directory: str):
    # TODO
    config_file = os.path.join(directory, "padlock_oligo_designer.yaml")
    # write the config file

    return config_file


def initialize_parameters(parser: ArgumentParser):
    parser.add_argument(
        "-o",
        "--output",
        help="path of output folder",
        required=True,
        type=str,
        metavar="",
    )
    parser.add_argument(
        "-c",
        "--config",
        help="path to config yaml file",
        default=None,
        type=str,
        metavar="",
    )
    parser.add_argument("--n_jobs", help="number of cores used", type=int, metavar="")

    # TODO all the others

    args = parser.parse_args()
    args = vars(args)
    if args["config"] is None:
        config_file = generate_config_file()
    # read the config file
    with open(config_file, "r") as handle:
        config = yaml.safe_load(handle)

    # update the config file values with the given one in the command line
    for param, value in args:
        if value is not None and param != "config":
            # overwrite the config file
            config[param] = args[param]

    return config


def main():
    """Pipeline to designe Padlock oligos. The program supports two ways to recieve the input parameters:

    - command line input
    - configuration file (recommended)

    A standard configuration file can be generated using the following command ...

    Since the number of input parameters requested is too high to be handled only through the command line,
    the progam will use as baseline the configuration file (it will be automatically generated if it wasn't provided)
    and the parameters specified in the command line will be overwritten with the values given.

    REMARK: melting temperature parameters can be given only through the configuration file.
    """

    # get comman line arguments
    parser = ArgumentParser(
        prog="Padlock Oligo Designer",
        usage="%(prog)s [options]",
        description=__doc__,
        formatter_class=RawDescriptionHelpFormatter,
    )

    config = initialize_parameters(parser)
    dir_output = config["output"]
    Path(dir_output).mkdir(parents=True, exist_ok=True)
    logging.info(f"Results will be saved to: {config['output']}")

    logging.info("#########Start Pipeline#########")
    t_pipeline = time.time()

    logging.info("Time Pipeline: {} min".format(t_pipeline))
    logging.info("#########End Pipeline#########")


if __name__ == "__main__":
    main()
