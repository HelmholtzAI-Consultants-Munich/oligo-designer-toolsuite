import inspect
import logging
import os
from pathlib import Path
import warnings
from argparse import ArgumentParser, RawDescriptionHelpFormatter
from typing import List, Union, Literal
import yaml

from Bio.SeqUtils import MeltingTemp as mt

from oligo_designer_toolsuite.database import OligoDatabase, ReferenceDatabase, OligoAttributes
from oligo_designer_toolsuite.oligo_property_filter import (
    GCContentFilter,
    HardMaskedSequenceFilter,
    HomopolymericRunsFilter,
    MeltingTemperatureNNFilter,
    ProhibitedSequenceFilter,
    PropertyFilter,
    SecondaryStructureFilter,
    SoftMaskedSequenceFilter,
    HomodimerFilter,
)
from oligo_designer_toolsuite.oligo_specificity_filter import (
    ExactMatchFilter,
    CrossHybridizationFilter,
    SpecificityFilter,
    AlignmentSpecificityFilter,
    BlastNFilter,
    BowtieFilter,
    RemoveByLargerRegionPolicy,
    HybridizationProbabilityFilter,
)
from oligo_designer_toolsuite.oligo_selection import OligosetGenerator, padlock_heuristic_selection
from oligo_designer_toolsuite.oligo_efficiency_filter import TmGCOligoScoring, AverageSetScoring
from oligo_designer_toolsuite.pipelines._base_oligo_designer import BaseOligoDesigner
from oligo_designer_toolsuite.pipelines._utils import log_parameters

# _ALIGNMENT_METHODS = Literal["blastn", "blastn_seedregion", "blastn_seedregion_ligationsite", "bowtie", "bowtie2"] #this should go to constants?
# _POLICIES = Literal["larger_region", "degree"]


class OligoSeq(BaseOligoDesigner):
    """_summary_

    :param BaseOligoDesigner: _description_
    :type BaseOligoDesigner: _type_
    """

    def filter_by_property(
        self,
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
        n_jobs: int = 1,
    ):

        ##### log parameters #####
        logging.info("Parameters Property Filters:")
        args, _, _, values = inspect.getargvalues(inspect.currentframe())
        parameters = {i: values[i] for i in args}
        log_parameters(parameters)

        num_genes_before, num_oligos_before = self._get_oligo_database_info(oligo_database.database)

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
        # add homodimer TODO
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
        if self.write_intermediate_steps:
            file_database = oligo_database.save_database(filename="oligo_database_property_filter")
        else:
            file_database = ""

        ##### loggig database information #####
        num_genes_after, num_oligos_after = self._get_oligo_database_info(oligo_database.database)
        logging.info(
            f"Step - Filter Oligos by Sequence Property: the database contains {num_oligos_after} oligos from {num_genes_after} genes, while {num_oligos_before - num_oligos_after} oligos and {num_genes_before - num_genes_after} genes have been deleted in this step."
        )

        return oligo_database, file_database

    def filter_by_specificity(
        self,
        oligo_database: OligoDatabase,
        files_fasta_reference_database: List[str],
        cross_hybridization_alignment_method: str,
        cross_hybridization_search_parameters: dict,
        cross_hybridization_hit_parameters: dict,
        hybridization_probability_alignment_method: str,
        hybridization_probability_search_parameters: dict,
        hybridization_probability_hit_parameters: dict,
        hybridization_probability_threshold: float,
        n_jobs: int = 1,
    ):

        ##### log parameters #####
        logging.info("Parameters Specificty Filters:")
        args, _, _, values = inspect.getargvalues(inspect.currentframe())
        parameters = {i: values[i] for i in args}
        log_parameters(parameters)

        num_genes_before, num_oligos_before = self._get_oligo_database_info(oligo_database.database)

        # define reference database
        reference_database = ReferenceDatabase(dir_output=self.dir_output)
        reference_database.load_sequences_from_fasta(
            files_fasta=files_fasta_reference_database, database_overwrite=False
        )

        # specificity filters
        # for the future: add policy
        exact_matches = ExactMatchFilter()
        cross_hybridization_aligner = self._get_alignment_method(
            alignment_method=cross_hybridization_alignment_method,
            search_parameters=cross_hybridization_search_parameters,
            hit_parameters=cross_hybridization_hit_parameters,
        )
        cross_hybridization_policy = RemoveByLargerRegionPolicy()
        cross_hybridization = CrossHybridizationFilter(
            policy=cross_hybridization_policy,
            alignment_method=cross_hybridization_aligner,
            dir_output=self.dir_output,
        )
        hybridization_probability_aligner = self._get_alignment_method(
            alignment_method=hybridization_probability_alignment_method,
            search_parameters=hybridization_probability_search_parameters,
            hit_parameters=hybridization_probability_hit_parameters,
        )
        hybridization_probability = HybridizationProbabilityFilter(
            alignment_method=hybridization_probability_aligner,
            threshold=hybridization_probability_threshold,
            dir_output=self.dir_output,
        )

        filters = [exact_matches, cross_hybridization, hybridization_probability]  # TODO add hybrid. prob
        specificity_filter = SpecificityFilter(filters=filters)
        oligo_database = specificity_filter.apply(
            sequence_type="oligo",
            oligo_database=oligo_database,
            reference_database=reference_database,
            n_jobs=n_jobs,
        )

        # write the intermediate result in a file
        if self.write_intermediate_steps:
            file_database = oligo_database.save_database(filename="oligo_database_specificty_filter")
        else:
            file_database = ""

        ##### loggig database information #####
        num_genes_after, num_oligos_after = self._get_oligo_database_info(oligo_database.database)
        logging.info(
            f"Step - Filter Oligos by Sequence Specificity: the database contains {num_oligos_after} oligos from {num_genes_after} genes, while {num_oligos_before - num_oligos_after} oligos and {num_genes_before - num_genes_after} genes have been deleted in this step."
        )

        return oligo_database, file_database

    def create_oligo_sets(
            self,
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
            n_jobs: int = 1,
    ):
        ##### log parameters #####
        logging.info("Parameters Oligo Selection:")
        args, _, _, values = inspect.getargvalues(inspect.currentframe())
        parameters = {i: values[i] for i in args}
        log_parameters(parameters)

        num_genes_before, num_oligos_before = self._get_oligo_database_info(oligo_database.database)

        # add required fileds
        # oligo_attributes = OligoAttributes()
        # oligo_database = oligo_attributes.calculate_GC_content(oligo_database=oligo_database, sequence_type="oligo")
        # oligo_database = oligo_attributes.calculate_TmNN(
        #     oligo_database=oligo_database,
        #     sequence_type="oligo",
        #     Tm_parameters=Tm_parameters,
        #     Tm_chem_correction_parameters=Tm_chem_correction_parameters
        # )

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
        if self.write_intermediate_steps:
            file_database = oligo_database.save_database(filename="oligo_database_specificty_filter")
            file_oligosets = oligo_database.write_oligosets()
        else:
            file_database = ""
            file_oligosets = ""

        ##### loggig database information #####
        num_genes_after, num_oligos_after = self._get_oligo_database_info(oligo_database.database)
        logging.info(
            f"Step - Filter Oligos by Sequence Efficiency: the database contains {num_oligos_after} oligos from {num_genes_after} genes, while {num_oligos_before - num_oligos_after} oligos and {num_genes_before - num_genes_after} genes have been deleted in this step."
        )

        return oligo_database, file_database, file_oligosets

    def create_final_sequences():
        pass

    def _get_alignment_method(self, alignment_method, search_parameters, hit_parameters):
        if alignment_method == "blastn":
            return BlastNFilter(
                search_parameters=search_parameters,
                hit_parameters=hit_parameters,
                dir_output=self.dir_output,
            )
        elif alignment_method == "bowtie":
            return BowtieFilter(
                search_parameters=search_parameters,
                hit_parameters=hit_parameters,
                dir_output=self.dir_output,
            )
        else:
            raise ValueError(f"The alignment method {alignment_method} is not supported.")
        
    def compute_oligo_attributes(
            self,
            oligo_database: OligoDatabase,
            Tm_parameters: dict,
            Tm_chem_correction_parameters: dict,
    ):
        oligo_attributes = OligoAttributes()
        oligo_database = oligo_attributes.calculate_oligo_length(oligo_database=oligo_database)
        oligo_database = oligo_attributes.calculate_GC_content(oligo_database=oligo_database, sequence_type="oligo")
        oligo_database = oligo_attributes.calculate_TmNN(
            oligo_database=oligo_database, 
            sequence_type="oligo",
            Tm_parameters=Tm_parameters,
            Tm_chem_correction_parameters=Tm_chem_correction_parameters
        )
        oligo_database = oligo_attributes.calculate_num_targeted_transcripts(oligo_database=oligo_database)
        oligo_database = oligo_attributes.calculate_isoform_consensus(oligo_database=oligo_database)
        oligo_database = oligo_attributes.calculate_length_selfcomplement(oligo_database=oligo_database, sequence_type="oligo")
        oligo_database = oligo_attributes.calculate_DG_secondary_structure(oligo_database=oligo_database, sequence_type="oligo")





def main():
    """
    Command line tool to run the OLIGO-SEQ pipeline to design oligos for FISH experiments. To run the tool use the command: ``oligo_seq [options]``.

    The program supports two ways to recieve the input parameters:

    - command line input
    - configuration file (recommended)

    A standard configuration file can be generated using the following command ``oligo_seq_config [options]``

    Since the number of input parameters requested is too high to be handled only through the command line,
    the progam will use as baseline the configuration file (it will be automatically generated if it wasn't provided)
    and the parameters specified in the command line will be overwritten with the values given.

    REMARK: melting temperature parameters can be given only through the configuration file.
    """
    parser = ArgumentParser(
        prog="Oligo Seq Designer",
        usage="oligo_seq [options]",
        description=__doc__,
        formatter_class=RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "-c",
        "--config",
        help="path to the config yaml file, str",
        default=None,
        type=str,
        metavar="",
    )
    args = parser.parse_args()
    args = vars(args)
    # read the config file
    with open(args["config"], "r") as handle:
        config = yaml.safe_load(handle)
    dir_output = os.path.abspath(config["dir_output"])
    Path(dir_output).mkdir(parents=True, exist_ok=True)

    ##### Initialize OligoDesigner Class #####
    oligo_designer = OligoSeq(dir_output=dir_output, log_name="oligo_seq")

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
    oligo_database, file_database = oligo_designer.create_oligo_database(
        regions=genes,
        files_fasta_oligo_database=config["files_fasta_oligo_database"],
        oligo_length_min=config["oligo_length_min"],
        oligo_length_max=config["oligo_length_max"],
        min_oligos_per_region=config["min_oligos_per_gene"],
        n_jobs=config["n_jobs"],
    )

    ##### filter oligos by property #####

    oligo_database, file_database = oligo_designer.filter_by_property(
        oligo_database,
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
        n_jobs=config["n_jobs"],
    )

    # ##### filter oligos by specificity #####
    oligo_database, file_database = oligo_designer.filter_by_specificity(
        oligo_database,
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
        n_jobs=config["n_jobs"],
    )

    ##### create oligo sets #####
    oligo_database, file_database, dir_oligosets = oligo_designer.create_oligo_sets(
        oligo_database,
        Tm_min=config["Tm_min"],
        Tm_max=config["Tm_max"],
        Tm_opt=config["Tm_opt"],
        Tm_parameters=Tm_parameters,
        Tm_chem_correction_parameters=config["Tm_chem_correction_parameters"],
        GC_content_min=config["GC_content_min"],
        GC_content_max=config["GC_content_max"],
        GC_content_opt=config["GC_content_opt"],
        oligoset_size=config["oligoset_size"],
        min_oligoset_size=config["min_oligoset_size"],
        max_oligos=config["max_graph_size"],
        n_sets=config["n_sets"],
        n_jobs=config["n_jobs"],
    )

    logging.info(f"Oligo sets were saved in {dir_oligosets}")
    logging.info("##### End of the pipeline. #####")


if __name__ == "__main__":
    main()
