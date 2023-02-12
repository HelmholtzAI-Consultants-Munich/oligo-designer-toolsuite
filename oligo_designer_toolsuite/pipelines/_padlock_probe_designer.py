############################################
# imports
############################################

import logging
import os
import time
import warnings
from argparse import ArgumentParser, RawDescriptionHelpFormatter
from datetime import datetime
from pathlib import Path

import yaml
from Bio.SeqUtils import MeltingTemp as mt

from oligo_designer_toolsuite.database import (
    CustomOligoDB,
    CustomReferenceDB,
    EnsemblOligoDB,
    NcbiOligoDB,
)
from oligo_designer_toolsuite.experiment_specific import PadlockSequenceDesigner
from oligo_designer_toolsuite.oligo_efficiency import (
    PadlockOligoScoring,
    PadlockSetScoring,
)
from oligo_designer_toolsuite.oligo_property_filter import (
    GCContent,
    MaskedSequences,
    MeltingTemperature,
    PadlockArms,
    PropertyFilter,
)
from oligo_designer_toolsuite.oligo_selection import (
    OligosetGenerator,
    padlock_heuristic_selection,
)
from oligo_designer_toolsuite.oligo_specificity_filter import (
    Blastn,
    BowtieSeedRegion,
    ExactMatches,
    LigationRegionCreation,
    SpecificityFilter,
)
from oligo_designer_toolsuite.pipelines._padlock_probe_designer_config import (
    generate_ncbi_config,
)


def oligo_database_info(oligo_database: dict):
    """Counts the number of oligos and genees in the database."""

    genes = oligo_database.keys()
    num_genes = len(genes)
    num_oligos = 0
    for gene in genes:
        num_oligos += len(oligo_database[gene].keys())

    return num_genes, num_oligos


def generate_config_file(directory: str):
    directory = os.path.join(directory, "config")
    Path(directory).mkdir(parents=True, exist_ok=True)
    config_file = generate_ncbi_config(directory)  # function generating the config file
    warnings.warn(
        f"No config file was given, a file is generatede automatically in '{config_file}'."
    )
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
        help="path to a file containing the genes used to generate the oligos, if empty all the genes are used, str",
        type=str,
        default=None,
        metavar="",
    )
    parser.add_argument(
        "-mog",
        "--min_oligos_per_gene",
        help="genes with less that this number of oligos are removed, int",
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
        "-af",
        "--annotation_file",
        help="path to the annotation file (only for custom source), str",
        type=str,
        default=None,
        metavar="",
    )
    parser.add_argument(
        "-sf",
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
        "-mol",
        "--oligo_length_min",
        help="minimum length of oligos, int",
        type=int,
        default=None,
        metavar="",
    )
    parser.add_argument(
        "-Mol",
        "--oligo_length_max",
        help="max length of oligos, int",
        type=int,
        default=None,
        metavar="",
    )
    parser.add_argument(
        "--block_size",
        help="size of the exon junctions in the transcript used for the alignement methods, int",
        type=int,
        default=None,
        metavar="",
    )
    parser.add_argument(
        "--GC_content_min",
        help="minimum GC content of oligos, [0, 100]",
        type=float,
        default=None,
        metavar="",
    )
    parser.add_argument(
        "--GC_content_max",
        help="maximum GC content of oligos, [0, 100]",
        type=float,
        default=None,
        metavar="",
    )
    parser.add_argument(
        "--Tm_min",
        help="minimum melting temperature of oligos, float",
        type=float,
        default=None,
        metavar="",
    )
    parser.add_argument(
        "--Tm_max",
        help="maximum melting temperature of oligos, float",
        type=float,
        default=None,
        metavar="",
    )
    parser.add_argument(
        "--min_arm_length",
        help="minimum length of each arm, int",
        type=int,
        default=None,
        metavar="",
    )
    parser.add_argument(
        "--max_arm_Tm_dif",
        help="max melting temperature difference of both arms, float",
        type=float,
        default=None,
        metavar="",
    )
    parser.add_argument(
        "--arm_Tm_min",
        help="minimum melting temperature of each arm, float",
        type=float,
        default=None,
        metavar="",
    )
    parser.add_argument(
        "--arm_Tm_max",
        help="max melting temperature of each arm, float",
        type=float,
        default=None,
        metavar="",
    )
    parser.add_argument(
        "--word_size",
        help="word size for the blastn seed, int",
        type=int,
        default=None,
        metavar="",
    )
    parser.add_argument(
        "--coverage",
        help="minimum coverage between oligos and target sequence for blastn, [0,100]",
        type=float,
        default=None,
        metavar="",
    )
    parser.add_argument(
        "--percent_identity",
        help="maximum similarity between oligos and target sequences for blastn, [0,100]",
        type=float,
        default=None,
        metavar="",
    )
    parser.add_argument(
        "--strand",
        help="strand of the query sequence to search, str",
        type=str,
        default=None,
        metavar="",
    )
    parser.add_argument(
        "--ligation_region_size",
        help="size of the seed region around the ligation site for bowtie seed region filter, int",
        type=int,
        default=None,
        metavar="",
    )  # TODO write help function
    parser.add_argument(
        "--Tm_opt",
        help="optimal melting temperature of oligos, float",
        type=float,
        default=None,
        metavar="",
    )
    parser.add_argument(
        "--GC_content_opt",
        help="optimal GC content of oligos, [0, 100]",
        type=float,
        default=None,
        metavar="",
    )
    parser.add_argument(
        "--Tm_weight",
        help="weight of the Tm of the oligo in the efficiency score, float",
        type=float,
        default=None,
        metavar="",
    )
    parser.add_argument(
        "--GC_weight",
        help="weight of the GC content of the oligo in the efficiency score, float",
        type=float,
        default=None,
        metavar="",
    )
    parser.add_argument(
        "--oligoset_size",
        help="ideal number of oligos per oligoset, int",
        type=int,
        default=None,
        metavar="",
    )
    parser.add_argument(
        "--min_oligoset_size",
        help="minimum number of oligos per oligoset, int",
        type=int,
        default=None,
        metavar="",
    )
    parser.add_argument(
        "--n_sets",
        help="maximum number of sets per gene, int",
        type=int,
        default=None,
        metavar="",
    )
    parser.add_argument(
        "--detect_oligo_length_min",
        help="minimum number of oligos per oligoset, int",
        type=int,
        default=None,
        metavar="",
    )
    parser.add_argument(
        "--detect_oligo_length_max",
        help="maximum length of detection oligo, int",
        type=int,
        default=None,
        metavar="",
    )
    parser.add_argument(
        "--detect_oligo_Tm_opt",
        help="optimal melting temperature of detection oligo, float",
        type=float,
        default=None,
        metavar="",
    )

    args = parser.parse_args()
    args = vars(args)
    if args["config"] is None:
        args["config"] = generate_config_file(args["output"])

    # read the config file
    with open(args["config"], "r") as handle:
        config = yaml.safe_load(handle)

    # update the config file values with the given one in the command line
    for param, value in args.items():
        if value is not None and param != "config":
            # overwrite the config file
            config[param] = args[param]

    return config


def padlock_probe_designer():
    """Command line tool to run a pipeline to design Padlock Probes, to run the tool use the command: ``padlock_probe_designer [options]``.

    The program supports two ways to recieve the input parameters:

    - command line input
    - configuration file (recommended)

    A standard configuration file can be generated using the following command ``padlock_probe_designer_config [options]``

    Since the number of input parameters requested is too high to be handled only through the command line,
    the progam will use as baseline the configuration file (it will be automatically generated if it wasn't provided)
    and the parameters specified in the command line will be overwritten with the values given.

    REMARK: melting temperature parameters can be given only through the configuration file.
    """
    # logging
    timestamp = datetime.now()
    file_logger = f"log_padlock_probe_designer_{timestamp.year}-{timestamp.month}-{timestamp.day}-{timestamp.hour}-{timestamp.minute}.txt"
    logging.getLogger("padlock_probe_designer")
    logging.basicConfig(
        format="%(asctime)s [%(levelname)s] %(message)s",
        level=logging.NOTSET,
        handlers=[logging.FileHandler(file_logger), logging.StreamHandler()],
    )
    logging.captureWarnings(True)

    # get comman line arguments
    parser = ArgumentParser(
        prog="Padlock Probe Designer",
        usage="padlock_probe_designer [options]",
        description=padlock_probe_designer.__doc__,
        formatter_class=RawDescriptionHelpFormatter,
    )

    config = initialize_parameters(parser)
    dir_output = config["output"]
    Path(dir_output).mkdir(parents=True, exist_ok=True)

    logging.info("#########Parameter settings#########")
    for item, value in config.items():
        logging.info("{}: {}".format(item, value))
    logging.info(f"Results will be saved to: {config['output']}")

    logging.info("######### Start Pipeline #########")
    time.time()

    # initialize the necessary classes

    ######## Oligos generation ########

    if config["source"] == "ncbi":
        # dowload the fasta files formthe NCBI server
        oligo_database = NcbiOligoDB(
            oligo_length_min=config["oligo_length_min"],
            oligo_length_max=config["oligo_length_max"],
            species=config["species"],
            annotation_release=config["annotation_release"],
            n_jobs=config["n_jobs"],
            dir_output=dir_output,
            min_oligos_per_gene=config["min_oligos_per_gene"],
        )
    elif config["source"] == "ensembl":
        # dowload the fasta files formthe NCBI server
        oligo_database = EnsemblOligoDB(
            oligo_length_min=config["oligo_length_min"],
            oligo_length_max=config["oligo_length_max"],
            species=config["species"],
            annotation_release=config["annotation_release"],
            n_jobs=config["n_jobs"],
            dir_output=dir_output,
            min_oligos_per_gene=config["min_oligos_per_gene"],
        )
    elif config["source"] == "custom":
        # use already dowloaded files
        oligo_database = CustomOligoDB(
            oligo_length_min=config["oligo_length_min"],
            oligo_length_max=config["oligo_length_max"],
            species=config["species"],
            genome_assembly=config["genome_assembly"],
            annotation_release=config["annotation_release"],
            files_source=config["files_source"],
            annotation_file=config["annotation_file"],
            sequence_file=config["sequence_file"],
            n_jobs=config["n_jobs"],
            dir_output=dir_output,
            min_oligos_per_gene=config["min_oligos_per_gene"],
        )
    if config["write_removed_genes"]:
        logging.info(
            f"The genes with less than {config['min_oligos_per_gene']} probes will be removed and their names will be stored in '{oligo_database.file_removed_genes}'."
        )

    ######## Property filters ########

    # the melting temperature params need to be preprocessed
    Tm_params = config["Tm_parameters"]["shared"].copy()
    Tm_params.update(config["Tm_parameters"]["property_filter"])
    Tm_params["nn_table"] = getattr(mt, Tm_params["nn_table"])
    Tm_params["tmm_table"] = getattr(mt, Tm_params["tmm_table"])
    Tm_params["imm_table"] = getattr(mt, Tm_params["imm_table"])
    Tm_params["de_table"] = getattr(mt, Tm_params["de_table"])

    Tm_correction_param = config["Tm_correction_parameters"]["shared"].copy()
    Tm_correction_param.update(config["Tm_correction_parameters"]["property_filter"])

    # initialize the filters clasees
    masked_sequences = MaskedSequences()
    gc_content = GCContent(
        GC_content_min=config["GC_content_min"], GC_content_max=config["GC_content_max"]
    )
    melting_temperature = MeltingTemperature(
        Tm_min=config["Tm_min"],
        Tm_max=config["Tm_max"],
        Tm_parameters=Tm_params,
        Tm_correction_parameters=Tm_correction_param,
    )
    padlock_arms = PadlockArms(
        min_arm_length=config["min_arm_length"],
        max_arm_Tm_dif=config["max_arm_Tm_dif"],
        arm_Tm_min=config["arm_Tm_min"],
        arm_Tm_max=config["arm_Tm_max"],
        Tm_parameters=Tm_params,
        Tm_correction_parameters=Tm_correction_param,
    )
    # create the list of filters
    filters = [masked_sequences, gc_content, melting_temperature, padlock_arms]
    # initialize the property filter class
    property_filter = PropertyFilter(
        filters=filters,
        write_genes_with_insufficient_oligos=config["write_removed_genes"],
    )

    ######## Specificity filters ########

    dir_specificity = os.path.join(
        dir_output, "specificity_temporary"
    )  # folder where the temporary files will be written
    # generate the reference
    reference_database = CustomReferenceDB(
        species=oligo_database.species,
        genome_assembly=oligo_database.genome_assembly,
        annotation_release=oligo_database.annotation_release,
        files_source=oligo_database.files_source,
        annotation_file=oligo_database.annotation_file,
        sequence_file=oligo_database.sequence_file,
        dir_output=dir_output,
    )
    reference_database.create_reference_DB(
        block_size=config["block_size"]
    )  # leave the standard parameters

    # intialize the filter classes
    exact_mathces = ExactMatches(dir_specificity=dir_specificity)
    seed_ligation = LigationRegionCreation(
        ligation_region_size=config["ligation_region_size"]
    )
    seed_region = BowtieSeedRegion(
        dir_specificity=dir_specificity, seed_region_creation=seed_ligation
    )
    blastn = Blastn(
        dir_specificity=dir_specificity,
        word_size=config["word_size"],
        percent_identity=config["percent_identity"],
        coverage=config["coverage"],
        strand=config["strand"],
    )
    filters = [exact_mathces, seed_region, blastn]
    # initialize the specificity filter class
    specificity_filter = SpecificityFilter(
        filters=filters,
        write_genes_with_insufficient_oligos=config["write_removed_genes"],
    )

    ######## Oligoset generation ########

    # initialize the scoring classes
    oligos_scoring = PadlockOligoScoring(
        Tm_min=config["Tm_min"],
        Tm_opt=config["Tm_opt"],
        Tm_max=config["Tm_max"],
        GC_content_min=config["GC_content_min"],
        GC_content_opt=config["GC_content_opt"],
        GC_content_max=config["GC_content_max"],
        Tm_weight=config["Tm_weight"],
        GC_weight=config["GC_weight"],
    )
    set_scoring = PadlockSetScoring()

    # initialize the oligoset generator class
    oligoset_generator = OligosetGenerator(
        oligoset_size=config["oligoset_size"],
        min_oligoset_size=config["min_oligoset_size"],
        oligos_scoring=oligos_scoring,
        set_scoring=set_scoring,
        heurustic_selection=padlock_heuristic_selection,
        write_genes_with_insufficient_oligos=config["write_removed_genes"],
    )

    ######## Padlock sequences ########

    # preprocessing of themelting temperature parameters
    Tm_params = config["Tm_parameters"]["shared"].copy()
    Tm_params.update(config["Tm_parameters"]["detection_oligo"])
    Tm_params["nn_table"] = getattr(mt, Tm_params["nn_table"])
    Tm_params["tmm_table"] = getattr(mt, Tm_params["tmm_table"])
    Tm_params["imm_table"] = getattr(mt, Tm_params["imm_table"])
    Tm_params["de_table"] = getattr(mt, Tm_params["de_table"])

    Tm_correction_param = config["Tm_correction_parameters"]["shared"].copy()
    Tm_correction_param.update(config["Tm_correction_parameters"]["detection_oligo"])
    # initilize the padlock sequence designer class
    padlock_sequence_designer = PadlockSequenceDesigner(
        detect_oligo_length_min=config["detect_oligo_length_min"],
        detect_oligo_length_max=config["detect_oligo_length_max"],
        detect_oligo_Tm_opt=config["detect_oligo_Tm_opt"],
        Tm_parameters=Tm_params,
        Tm_correction_parameters=Tm_correction_param,
    )

    ######## Pipeline ########

    # read the genes file
    if config["file_genes"] is None:
        warnings.warn(
            "No file containing the genes was provided, all the genes are ussed to generate the probes. This chioce can use a lot of resources."
        )
        genes = None
    else:
        with open(config["file_genes"]) as handle:
            lines = handle.readlines()
            genes = [line.rstrip() for line in lines]

    # generate the oligo sequences from gene transcripts
    oligo_database.create_oligos_DB(genes=genes, region="transcripts")

    num_genes_init, num_oligos_init = oligo_database_info(oligo_database.oligos_DB)
    logging.info(f"Probes generation:")
    logging.info(
        f"at the beginning the database contains {num_oligos_init} probes from {num_genes_init} genes."
    )

    # apply property filter to the database
    oligo_database = property_filter.apply(
        oligo_database=oligo_database, n_jobs=config["n_jobs"]
    )
    # write the intermediate result in a file
    if config["write_intermediate_steps"]:
        oligo_database.write_oligos_DB(
            format=config["file_format"], dir_oligos_DB="property_filter"
        )

    num_genes_prop, num_oligos_prop = oligo_database_info(oligo_database.oligos_DB)
    logging.info(f"Property filters: ")
    logging.info(
        f"after the property filters have been applied the database contains {num_oligos_prop} probes from {num_genes_prop} genes,"
    )
    logging.info(
        f"{num_oligos_init - num_oligos_prop} probes and {num_genes_init - num_genes_prop} genes have been deleted in this step."
    )

    # apply specificty filters to the database
    oligo_database = specificity_filter.apply(
        oligo_database=oligo_database,
        reference_database=reference_database,
        n_jobs=config["n_jobs"],
    )
    # write the intermediate result
    if config["write_intermediate_steps"]:
        oligo_database.write_oligos_DB(
            format=config["file_format"], dir_oligos_DB="specificity_filter"
        )
    num_genes_spec, num_oligos_spec = oligo_database_info(oligo_database.oligos_DB)
    logging.info(f"Specificity filters:")
    logging.info(
        f"after the specificity filters have been applied the database contains {num_oligos_spec} probes from {num_genes_spec} genes, "
    )
    logging.info(
        f"{num_oligos_prop - num_oligos_spec} probes and {num_genes_prop - num_genes_spec} genes have been deleted in this step."
    )

    # generate the oligoset
    oligo_database = oligoset_generator.apply(
        oligo_database=oligo_database, n_sets=config["n_sets"], n_jobs=config["n_jobs"]
    )
    # write the intermediate result
    if config["write_intermediate_steps"]:
        oligo_database.write_oligosets(dir_oligosets="oligosets")
    num_genes_set, num_oligos_set = oligo_database_info(oligo_database.oligos_DB)
    logging.info(f"Oligosets generation:")
    logging.info(
        f"after the oligosets have been generated the database contains {num_oligos_set} probes from {num_genes_set} genes, "
    )
    logging.info(
        f"{num_oligos_spec - num_oligos_set} probes and {num_genes_spec - num_genes_set} genes have been deleted in this step."
    )

    # generate the padlock sequence
    padlock_sequence_designer.design_padlocks(
        database=oligo_database, dir_padlock="padlock_sequences"
    )
    logging.info(f"Padlock sequences generation: ")
    logging.info(
        f"padlock sequences have been written in '{os.path.join(oligo_database.dir_output, 'padlock_sequences')}'."
    )

    logging.info("######### End Pipeline #########")


if __name__ == "__main__":
    padlock_probe_designer()
