############################################
# imports
############################################

import os
import yaml
import inspect
import logging
import warnings

from pathlib import Path
from datetime import datetime

from typing import List

from Bio.SeqUtils import MeltingTemp as mt

from oligo_designer_toolsuite.pipelines._utils import (
    log_parameters,
    base_parser,
    get_oligo_database_info,
)
from oligo_designer_toolsuite.sequence_generator import OligoSequenceGenerator
from oligo_designer_toolsuite.database import OligoDatabase, ReferenceDatabase, OligoAttributes
from oligo_designer_toolsuite.oligo_property_filter import (
    GCContentFilter,
    HardMaskedSequenceFilter,
    HomopolymericRunsFilter,
    MeltingTemperatureNNFilter,
    PropertyFilter,
    SecondaryStructureFilter,
    SoftMaskedSequenceFilter,
    HomodimerFilter,
)
from oligo_designer_toolsuite.oligo_specificity_filter import (
    ExactMatchFilter,
    BlastNFilter,
    BowtieFilter,
    CrossHybridizationFilter,
    HybridizationProbabilityFilter,
    RemoveByLargerRegionPolicy,
    SpecificityFilter,
)
from oligo_designer_toolsuite.oligo_selection import OligosetGenerator, padlock_heuristic_selection
from oligo_designer_toolsuite.oligo_efficiency_filter import TmGCOligoScoring, AverageSetScoring


############################################
# Oligo-seq Designer Functions
############################################


def create_oligo_database(
    regions: list,
    oligo_length_min: int,
    oligo_length_max: int,
    files_fasta_oligo_database: list[str],
    min_oligos_per_region: int,
    write_removed_genes: bool,
    write_intermediate_steps: bool,
    dir_output: str,
    n_jobs: int,
):
    ##### log parameters #####
    logging.info("Parameters Create Database:")
    args, _, _, values = inspect.getargvalues(inspect.currentframe())
    parameters = {i: values[i] for i in args}
    log_parameters(parameters)

    ##### creating the oligo sequences #####
    oligo_sequences = OligoSequenceGenerator(dir_output=dir_output)
    oligo_fasta_file = oligo_sequences.create_sequences_sliding_window(
        filename_out="oligo_seq_designer_oligos",
        files_fasta_in=files_fasta_oligo_database,
        length_interval_sequences=(oligo_length_min, oligo_length_max),
        region_ids=regions,
        n_jobs=n_jobs,
    )

    ##### creating the oligo database #####
    # oligo database
    oligo_database = OligoDatabase(
        min_oligos_per_region=min_oligos_per_region,
        write_regions_with_insufficient_oligos=True,
        dir_output=dir_output,
    )
    # load the oligo sequences
    oligo_database.load_sequences_from_fasta(
        files_fasta=[oligo_fasta_file],
        sequence_type="oligo",
        region_ids=regions,
    )

    ##### loggig database information #####
    if write_removed_genes:
        logging.info(
            f"Genes with <= {min_oligos_per_region} oligos will be removed from the oligo database and stored in '{oligo_database.file_removed_regions}'."
        )

    num_genes, num_oligos = get_oligo_database_info(oligo_database.database)
    logging.info(f"Step - Generate oligos: the database contains {num_oligos} oligos from {num_genes} genes.")

    ##### save database #####
    if write_intermediate_steps:
        file_database = oligo_database.save_database(filename="oligo_database_initial.txt")
    else:
        file_database = ""

    return oligo_database, file_database


def filter_by_property(
    oligo_database: OligoDatabase,
    GC_content_min: int,
    GC_content_max: int,
    Tm_min: int,
    Tm_max: int,
    secondary_structures_T: float,
    secondary_structures_threshold_deltaG: float,
    homopolymeric_base_n: str,
    homodimer_max_len_selfcomp: int,
    Tm_parameters: dict,
    Tm_chem_correction_parameters: dict,
    write_intermediate_steps: bool,
    n_jobs: int,
):

    ##### log parameters #####
    logging.info("Parameters Property Filters:")
    args, _, _, values = inspect.getargvalues(inspect.currentframe())
    parameters = {i: values[i] for i in args}
    log_parameters(parameters)

    num_genes_before, num_oligos_before = get_oligo_database_info(oligo_database.database)

    # define the filters
    hard_masked_sequences = HardMaskedSequenceFilter()
    soft_masked_sequences = SoftMaskedSequenceFilter()
    gc_content = GCContentFilter(GC_content_min=GC_content_min, GC_content_max=GC_content_max)
    melting_temperature = MeltingTemperatureNNFilter(
        Tm_min=Tm_min,
        Tm_max=Tm_max,
        Tm_parameters=Tm_parameters,
        Tm_chem_correction_parameters=Tm_chem_correction_parameters,
    )
    secondary_sctructure = SecondaryStructureFilter(
        T=secondary_structures_T,
        thr_DG=secondary_structures_threshold_deltaG,
    )
    homopolymeric_runs = HomopolymericRunsFilter(
        base_n=homopolymeric_base_n,
    )
    homodimer = HomodimerFilter(
        max_len_selfcomp=homodimer_max_len_selfcomp,
    )

    filters = [
        hard_masked_sequences,
        soft_masked_sequences,
        gc_content,
        melting_temperature,
        secondary_sctructure,
        homopolymeric_runs,
        homodimer,
    ]

    # initialize the preoperty filter class
    property_filter = PropertyFilter(filters=filters)

    # filter the database
    oligo_database = property_filter.apply(
        sequence_type="oligo",
        oligo_database=oligo_database,
        n_jobs=n_jobs,
    )

    # write the intermediate result in a file
    if write_intermediate_steps:
        file_database = oligo_database.save_database(filename="oligo_database_property_filter")
    else:
        file_database = ""

    ##### loggig database information #####
    num_genes_after, num_oligos_after = get_oligo_database_info(oligo_database.database)
    logging.info(
        f"Step - Filter Oligos by Sequence Property: the database contains {num_oligos_after} oligos from {num_genes_after} genes, while {num_oligos_before - num_oligos_after} oligos and {num_genes_before - num_genes_after} genes have been deleted in this step."
    )

    return oligo_database, file_database


def filter_by_specificity(
    oligo_database: OligoDatabase,
    files_fasta_reference_database: List[str],
    cross_hybridization_alignment_method: str,
    cross_hybridization_search_parameters: dict,
    cross_hybridization_hit_parameters: dict,
    hybridization_probability_alignment_method: str,
    hybridization_probability_search_parameters: dict,
    hybridization_probability_hit_parameters: dict,
    hybridization_probability_threshold: float,
    write_intermediate_steps: bool,
    dir_output: str,
    n_jobs: int = 1,
):
    def _get_alignment_method(alignment_method, search_parameters, hit_parameters, dir_output):
        if alignment_method == "blastn":
            return BlastNFilter(
                search_parameters=search_parameters,
                hit_parameters=hit_parameters,
                dir_output=dir_output,
            )
        elif alignment_method == "bowtie":
            return BowtieFilter(
                search_parameters=search_parameters,
                hit_parameters=hit_parameters,
                dir_output=dir_output,
            )
        else:
            raise ValueError(f"The alignment method {alignment_method} is not supported.")

    ##### log parameters #####
    logging.info("Parameters Specificty Filters:")
    args, _, _, values = inspect.getargvalues(inspect.currentframe())
    parameters = {i: values[i] for i in args}
    log_parameters(parameters)

    num_genes_before, num_oligos_before = get_oligo_database_info(oligo_database.database)

    ##### define reference database #####
    reference_database = ReferenceDatabase(dir_output=dir_output)
    reference_database.load_sequences_from_fasta(
        files_fasta=files_fasta_reference_database, database_overwrite=False
    )

    ##### specificity filters #####
    # removing duplicated oligos from the region with the most oligos
    exact_matches_policy = RemoveByLargerRegionPolicy()
    exact_matches = ExactMatchFilter(policy=exact_matches_policy)

    cross_hybridization_aligner = _get_alignment_method(
        alignment_method=cross_hybridization_alignment_method,
        search_parameters=cross_hybridization_search_parameters,
        hit_parameters=cross_hybridization_hit_parameters,
        dir_output=dir_output,
    )
    cross_hybridization_policy = RemoveByLargerRegionPolicy()
    cross_hybridization = CrossHybridizationFilter(
        policy=cross_hybridization_policy,
        alignment_method=cross_hybridization_aligner,
        dir_output=dir_output,
    )

    hybridization_probability_aligner = _get_alignment_method(
        alignment_method=hybridization_probability_alignment_method,
        search_parameters=hybridization_probability_search_parameters,
        hit_parameters=hybridization_probability_hit_parameters,
        dir_output=dir_output,
    )
    hybridization_probability = HybridizationProbabilityFilter(
        alignment_method=hybridization_probability_aligner,
        threshold=hybridization_probability_threshold,
        dir_output=dir_output,
    )

    # TODO check correct sequence types
    # filters_oligo_sequence = [exact_matches]
    # filters_target_sequence = []
    filters = [exact_matches, cross_hybridization, hybridization_probability]
    specificity_filter = SpecificityFilter(filters=filters)
    oligo_database = specificity_filter.apply(
        sequence_type="target",
        oligo_database=oligo_database,
        reference_database=reference_database,
        n_jobs=n_jobs,
    )

    # write the intermediate result in a file
    if write_intermediate_steps:
        file_database = oligo_database.save_database(filename="oligo_database_specificty_filter")
    else:
        file_database = ""

    ##### loggig database information #####
    num_genes_after, num_oligos_after = get_oligo_database_info(oligo_database.database)
    logging.info(
        f"Step - Filter Oligos by Sequence Specificity: the database contains {num_oligos_after} oligos from {num_genes_after} genes, while {num_oligos_before - num_oligos_after} oligos and {num_genes_before - num_genes_after} genes have been deleted in this step."
    )

    return oligo_database, file_database


def compute_oligo_attributes(
    oligo_database: OligoDatabase,
    secondary_structures_T: float,
    Tm_parameters: dict,
    Tm_chem_correction_parameters: dict,
):
    oligo_attributes = OligoAttributes()
    oligo_database = oligo_attributes.calculate_oligo_length(oligo_database=oligo_database)
    oligo_database = oligo_attributes.calculate_GC_content(
        oligo_database=oligo_database, sequence_type="oligo"
    )
    oligo_database = oligo_attributes.calculate_TmNN(
        oligo_database=oligo_database,
        sequence_type="oligo",
        Tm_parameters=Tm_parameters,
        Tm_chem_correction_parameters=Tm_chem_correction_parameters,
    )
    oligo_database = oligo_attributes.calculate_num_targeted_transcripts(oligo_database=oligo_database)
    oligo_database = oligo_attributes.calculate_isoform_consensus(oligo_database=oligo_database)
    oligo_database = oligo_attributes.calculate_length_selfcomplement(
        oligo_database=oligo_database, sequence_type="oligo"
    )
    oligo_database = oligo_attributes.calculate_DG_secondary_structure(
        oligo_database=oligo_database, sequence_type="oligo", T=secondary_structures_T
    )

    return oligo_database


def create_oligo_sets(
    oligo_database: OligoDatabase,
    Tm_min: float,
    Tm_opt: float,
    Tm_max: float,
    Tm_parameters: dict,
    Tm_chem_correction_parameters: dict,
    GC_content_min: float,
    GC_content_opt: float,
    GC_content_max: float,
    oligoset_size: int,
    min_oligoset_size: int,
    max_oligos: int,
    n_sets: int,
    write_intermediate_steps: bool,
    n_jobs: int,
):
    ##### log parameters #####
    logging.info("Parameters Oligo Selection:")
    args, _, _, values = inspect.getargvalues(inspect.currentframe())
    parameters = {i: values[i] for i in args}
    log_parameters(parameters)

    num_genes_before, num_oligos_before = get_oligo_database_info(oligo_database.database)

    oligos_scoring = TmGCOligoScoring(
        Tm_min=Tm_min,
        Tm_opt=Tm_opt,
        Tm_max=Tm_max,
        GC_content_min=GC_content_min,
        GC_content_opt=GC_content_opt,
        GC_content_max=GC_content_max,
        Tm_parameters=Tm_parameters,
        Tm_chem_correction_parameters=Tm_chem_correction_parameters,
    )
    set_scoring = AverageSetScoring()
    oligoset_generator = OligosetGenerator(
        oligoset_size=oligoset_size,
        min_oligoset_size=min_oligoset_size,
        oligos_scoring=oligos_scoring,
        set_scoring=set_scoring,
        heurustic_selection=padlock_heuristic_selection,
        max_oligos=max_oligos,
    )
    oligo_database = oligoset_generator.apply(
        oligo_database=oligo_database,
        sequence_type="oligo",
        n_sets=n_sets,
        n_jobs=n_jobs,
    )

    # write the intermediate result in a file
    if write_intermediate_steps:
        file_database = oligo_database.save_database(filename="oligo_database_oligosets")
        file_oligosets = oligo_database.write_oligosets()
    else:
        file_database = ""
        file_oligosets = ""

    ##### loggig database information #####
    num_genes_after, num_oligos_after = get_oligo_database_info(oligo_database.database)
    logging.info(
        f"Step - Filter Oligos by Sequence Efficiency: the database contains {num_oligos_after} oligos from {num_genes_after} genes, while {num_oligos_before - num_oligos_after} oligos and {num_genes_before - num_genes_after} genes have been deleted in this step."
    )

    return oligo_database, file_database, file_oligosets


############################################
# Oligo-seq Designer Pipeline
############################################


def main():

    args = base_parser()

    ##### read the config file #####
    with open(args["config"], "r") as handle:
        config = yaml.safe_load(handle)

    ##### create the output folder #####
    dir_output = os.path.abspath(config["dir_output"])
    Path(dir_output).mkdir(parents=True, exist_ok=True)

    ##### setup logger #####
    timestamp = datetime.now()
    file_logger = os.path.join(
        dir_output,
        f"log_oligo_seq_designer_{timestamp.year}-{timestamp.month}-{timestamp.day}-{timestamp.hour}-{timestamp.minute}.txt",
    )
    logging.getLogger("log_name")
    logging.basicConfig(
        format="%(asctime)s [%(levelname)s] %(message)s",
        level=logging.NOTSET,
        handlers=[logging.FileHandler(file_logger), logging.StreamHandler()],
    )
    logging.captureWarnings(True)

    ##### read the genes file #####
    if config["file_regions"] is None:
        warnings.warn(
            "No gene list file was provided! All genes from fasta file are used to generate the oligos. This chioce can use a lot of resources."
        )
        genes = None
    else:
        with open(config["file_regions"]) as handle:
            lines = handle.readlines()
            genes = [line.rstrip() for line in lines]

    ##### preprocess melting temperature params #####
    Tm_parameters = config["Tm_parameters"]
    Tm_parameters["nn_table"] = getattr(mt, Tm_parameters["nn_table"])
    Tm_parameters["tmm_table"] = getattr(mt, Tm_parameters["tmm_table"])
    Tm_parameters["imm_table"] = getattr(mt, Tm_parameters["imm_table"])
    Tm_parameters["de_table"] = getattr(mt, Tm_parameters["de_table"])

    ##### create oligo database #####
    oligo_database, file_database = create_oligo_database(
        regions=genes,
        oligo_length_min=config["oligo_length_min"],
        oligo_length_max=config["oligo_length_max"],
        files_fasta_oligo_database=config["files_fasta_oligo_database"],
        min_oligos_per_region=config["min_oligos_per_gene"],
        write_removed_genes=config["write_removed_genes"],
        write_intermediate_steps=config["write_intermediate_steps"],
        dir_output=dir_output,
        n_jobs=config["n_jobs"],
    )

    ##### filter oligos by property #####
    oligo_database, file_database = filter_by_property(
        oligo_database=oligo_database,
        GC_content_min=config["GC_content_min"],
        GC_content_max=config["GC_content_max"],
        Tm_min=config["Tm_min"],
        Tm_max=config["Tm_max"],
        secondary_structures_T=config["secondary_structures_T"],
        secondary_structures_threshold_deltaG=config["secondary_structures_threshold_deltaG"],
        homopolymeric_base_n=config["homopolymeric_base_n"],
        homodimer_max_len_selfcomp=config["homodimer_max_len_selfcomp"],
        Tm_parameters=Tm_parameters,
        Tm_chem_correction_parameters=config["Tm_chem_correction_parameters"],
        write_intermediate_steps=config["write_intermediate_steps"],
        n_jobs=config["n_jobs"],
    )

    # ##### filter oligos by specificity #####
    oligo_database, file_database = filter_by_specificity(
        oligo_database=oligo_database,
        files_fasta_reference_database=config["files_fasta_reference_database"],
        cross_hybridization_alignment_method=config["cross_hybridization_alignment_method"],
        cross_hybridization_search_parameters=config[
            f"cross_hybridization_{config['cross_hybridization_alignment_method']}_search_parameters"
        ],
        cross_hybridization_hit_parameters=config[
            f"cross_hybridization_{config['cross_hybridization_alignment_method']}_hit_parameters"
        ],
        hybridization_probability_alignment_method=config["hybridization_probability_alignment_method"],
        hybridization_probability_search_parameters=config[
            f"hybridization_probability_{config['hybridization_probability_alignment_method']}_search_parameters"
        ],
        hybridization_probability_hit_parameters=config[
            f"hybridization_probability_{config['hybridization_probability_alignment_method']}_hit_parameters"
        ],
        hybridization_probability_threshold=config["hybridization_probability_threshold"],
        write_intermediate_steps=config["write_intermediate_steps"],
        dir_output=dir_output,
        n_jobs=config["n_jobs"],
    )

    ##### compute all required attributes #####
    logging.info("Computing Oligo Attributes")
    oligo_database = compute_oligo_attributes(
        oligo_database=oligo_database,
        secondary_structures_T=config["secondary_structures_T"],
        Tm_parameters=Tm_parameters,
        Tm_chem_correction_parameters=config["Tm_chem_correction_parameters"],
    )

    ##### create oligo sets #####
    oligo_database, file_database, dir_oligosets = create_oligo_sets(
        oligo_database=oligo_database,
        Tm_min=config["Tm_min"],
        Tm_opt=config["Tm_opt"],
        Tm_max=config["Tm_max"],
        Tm_parameters=Tm_parameters,
        Tm_chem_correction_parameters=config["Tm_chem_correction_parameters"],
        GC_content_min=config["GC_content_min"],
        GC_content_opt=config["GC_content_opt"],
        GC_content_max=config["GC_content_max"],
        oligoset_size=config["oligoset_size"],
        min_oligoset_size=config["min_oligoset_size"],
        max_oligos=config["max_graph_size"],
        n_sets=config["n_sets"],
        write_intermediate_steps=config["write_intermediate_steps"],
        n_jobs=config["n_jobs"],
    )

    logging.info(f"Oligo sets were saved in {dir_oligosets}")
    logging.info("##### End of the pipeline. #####")


if __name__ == "__main__":
    main()
