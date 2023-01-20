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
        help="path to the output folder, str",
        required=True,
        type=str,
        metavar="",
    )
    parser.add_argument(
        "-c",
        "--config",
        help="path to the config yaml file, str",
        default=None,
        type=str,
        metavar="",
    )
    parser.add_argument(
        "-n",
        "--n_jobs",
        help="number of cores used, int",
        type=int,
        default=None,
        metavar="",
    )
    # TODO how to use all the genes? add a worning if no info is added
    parser.add_argument(
        "-g",
        "--file_genes",
        help="path to a file containing the genes used to generate the oligos, str",
        type=str,
        default=None,
        metavar="",
    )
    parser.add_argument(
        "-mo",
        "--min_oligos_per_gene",
        help="geens wit less that his oligos are removed, int",
        type=int,
        default=None,
        metavar="",
    )
    parser.add_argument(
        "-wg",
        "--write_removed_genes",
        help="write in a file the removed genes, bool",
        type=bool,
        default=None,
        metavar="",
    )
    parser.add_argument(
        "-ws",
        "--write_intermediate_steps",
        help="write the oligo sequences after each step of the pipeline, bool",
        type=bool,
        default=None,
        metavar="",
    )
    parser.add_argument(
        "-f",
        "--file_format",
        help="format of the files containing the oligos, [tsv, gtf]",
        choices=["tsv", "gtf"],
        default=None,
        metavar="",
    )
    parser.add_argument(
        "-s",
        "--source",
        help="how to obtain the genimic files: dowload them from a server [ncbi, ensembl] or provide the files directly [custom]",
        choices=["ncbi", "ensembl", "custom"],
        default=None,
        metavar="",
    )
    parser.add_argument(
        "-sp",
        "--species",
        help="[human, mouse]",
        choices=["human", "mouse"],
        default=None,
        metavar="",
    )
    parser.add_argument(
        "--annotation_release",
        help="annotation release number or 'current' for the latest version, str",
        type=str,
        default=None,
        metavar="",
    )
    parser.add_argument(
        "--annotation_file",
        help="path to the annotation file (only for custom source), str",
        type=str,
        default=None,
        metavar="",
    )
    parser.add_argument(
        "--sequence_file",
        help="path to the sequence file (only for custom source), str",
        type=str,
        default=None,
        metavar="",
    )
    parser.add_argument(
        "--genome_assembly",
        help="(only for custom source), str",
        type=str,
        default=None,
        metavar="",
    )
    parser.add_argument(
        "--files_source",
        help="original source of the files (only for custom source), [NCBI, Ensembl]",
        choices=["NCBI", "Ensembl"],
        default=None,
        metavar="",
    )
    parser.add_argument(
        "--oligo_length_min",
        help="min length of oligos, int",
        type=int,
        default=None,
        metavar="",
    )
    parser.add_argument(
        "--oligo_length_max",
        help="max length of oligos, int",
        type=int,
        default=None,
        metavar="",
    )
    # parser.add_argument("--", help="", type=, default=None, metavar="")
    # parser.add_argument("--", help="", type=, default=None, metavar="")
    # parser.add_argument("--", help="", type=, default=None, metavar="")
    # parser.add_argument("--", help="", type=, default=None, metavar="")
    # parser.add_argument("--", help="", type=, default=None, metavar="")
    # parser.add_argument("--", help="", type=, default=None, metavar="")
    # parser.add_argument("--", help="", type=, default=None, metavar="")
    # parser.add_argument("--", help="", type=, default=None, metavar="")
    # parser.add_argument("--", help="", type=, default=None, metavar="")
    # parser.add_argument("--", help="", type=, default=None, metavar="")
    # parser.add_argument("--", help="", type=, default=None, metavar="")
    # parser.add_argument("--", help="", type=, default=None, metavar="")
    # parser.add_argument("--", help="", type=, default=None, metavar="")
    # parser.add_argument("--", help="", type=, default=None, metavar="")
    # parser.add_argument("--", help="", type=, default=None, metavar="")
    # parser.add_argument("--", help="", type=, default=None, metavar="")
    # parser.add_argument("--", help="", type=, default=None, metavar="")
    # parser.add_argument("--", help="", type=, default=None, metavar="")
    # parser.add_argument("--", help="", type=, default=None, metavar="")
    # parser.add_argument("--", help="", type=, default=None, metavar="")
    # parser.add_argument("--", help="", type=, default=None, metavar="")
    # parser.add_argument("--", help="", type=, default=None, metavar="")

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
