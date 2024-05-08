import inspect
import logging
import os
from pathlib import Path
import warnings
from argparse import ArgumentParser, RawDescriptionHelpFormatter
from typing import List, Union, Literal

from Bio.SeqUtils import MeltingTemp as mt

from oligo_designer_toolsuite.database import OligoDatabase, ReferenceDatabase
from oligo_designer_toolsuite.oligo_property_filter import (
    GCContentFilter,
    HardMaskedSequenceFilter,
    HomopolymericRunsFilter,
    MeltingTemperatureNNFilter,
    ProhibitedSequenceFilter,
    PropertyFilter,
    SecondaryStructureFilter,
    SoftMaskedSequenceFilter,
)
from oligo_designer_toolsuite.oligo_specificity_filter import (
    ExactMatchFilter,
    CrossHybridizationFilter,
    #TODO add hybridization probability
    SpecificityFilter,
    AlignmentSpecificityFilter,
    BlastNFilter,
    RemoveByDegreePolicy,
    HybridizationProbabilityFilter,
)
from oligo_designer_toolsuite.oligo_specificity_filter._filter_cross_hybridization import CrossHybridizationPolicy
from oligo_designer_toolsuite.oligo_selection import OligosetGenerator, padlock_heuristic_selection
from oligo_designer_toolsuite.oligo_efficiency_filter import PadlockOligoScoring, AverageSetScoring
from oligo_designer_toolsuite.pipelines._base_oligo_designer import BaseOligoDesigner
from oligo_designer_toolsuite.pipelines._utils import initialize_parameters

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
        GC_content_min: int = 40,
        GC_content_max: int = 60,
        Tm_min: int = 70,
        Tm_max: int = 80,
        secondary_structures_T: float = 76,
        secondary_structures_threshold_deltaG: float = 0,
        homopolymeric_base_n: str = {"A": 6},  # TODO: meaningful standard setting
        Tm_parameters: dict = {
            "check": True,
            "strict": True,
            "c_seq": None,
            "shift": 0,
            "nn_table": "DNA_NN3",
            "tmm_table": "DNA_TMM1",
            "imm_table": "DNA_IMM1",
            "de_table": "DNA_DE1",
            "dnac1": 50,
            "dnac2": 0,
            "selfcomp": False,
            "dNTPs": 0,
            "saltcorr": 7,
            "Na": 1000,
            "K": 0,
            "Tris": 0,
            "Mg": 0,
        },
        Tm_chem_correction_parameters: dict = {
            "DMSO": 0,
            "DMSOfactor": 0.75,
            "fmdfactor": 0.65,
            "fmdmethod": 1,
            "GC": None,
            "fmd": 20,
        },
        prohibited_sequence: Union[str, List[str]] = "",
        n_jobs: int = 1,
    ):

        ##### log parameters #####
        logging.info("Parameters Property Filters:")
        args, _, _, values = inspect.getargvalues(inspect.currentframe())
        parameters = {i: values[i] for i in args}
        self._log_parameters(parameters)

        num_genes_before, num_oligos_before = self._get_oligo_database_info(oligo_database.database)
        ##### preprocess melting temperature params #####
        Tm_parameters["nn_table"] = getattr(mt, Tm_parameters["nn_table"])
        Tm_parameters["tmm_table"] = getattr(mt, Tm_parameters["tmm_table"])
        Tm_parameters["imm_table"] = getattr(mt, Tm_parameters["imm_table"])
        Tm_parameters["de_table"] = getattr(mt, Tm_parameters["de_table"])

        # define the filters
        masked_sequences = HardMaskedSequenceFilter()
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
        prohibited_sequences = ProhibitedSequenceFilter(
            prohibited_sequence=prohibited_sequence  # TODO: understand what they really wnat
        )
        filters = [
            masked_sequences,
            soft_masked_sequences,
            gc_content,
            melting_temperature,
            secondary_sctructure,
            homopolymeric_runs,
            prohibited_sequences,
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
            file_database = oligo_database.save_database(filename_out="oligo_database_property_filter")
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
            genomic_regions: List[str],
            cross_hybridization_balstn_search_parameters: dict,
            cross_hybridization_balstn_hit_parameters: dict,
            hybridization_probability_balstn_search_parameters: dict,
            hybridization_probability_balstn_hit_parameters: dict,
            hybridization_probability_threshold: float,
            block_size: int,
            n_jobs: int = 1,
            
        ):

        ##### log parameters #####
        logging.info("Parameters Specificty Filters:")
        args, _, _, values = inspect.getargvalues(inspect.currentframe())
        parameters = {i: values[i] for i in args}
        self._log_parameters(parameters)

        num_genes_before, num_oligos_before = self._get_oligo_database_info(oligo_database.database)
        
        # define reference database
        fasta_files = self._parse_genomic_regions(genomic_regions=genomic_regions, block_size=block_size)
        reference_database = ReferenceDatabase()
        for fasta_file in fasta_files:
            reference_database.load_sequences_from_fasta(file_fasta=fasta_file, database_overwrite=False)

        # specificity filters
        exact_matches = ExactMatchFilter()
        cross_hybridization_alignment_method = BlastNFilter(
            blast_search_parameters=cross_hybridization_balstn_search_parameters,
            blast_hit_parameters=cross_hybridization_balstn_hit_parameters,
            dir_output=self.dir_output,
        )
        cross_hybridization_policy = RemoveByDegreePolicy()
        cross_hybridization = CrossHybridizationFilter(
            policy=cross_hybridization_policy, 
            specificity_filter=cross_hybridization_alignment_method,
        )
        hybridization_probability_alignment_method = BlastNFilter(
            blast_search_parameters=hybridization_probability_balstn_search_parameters,
            blast_hit_parameters=hybridization_probability_balstn_hit_parameters,
            dir_output=self.dir_output,
        )
        hybridization_probability = HybridizationProbabilityFilter(
            alignment_method=hybridization_probability_alignment_method,
            hybridization_probability_threshold=hybridization_probability_threshold,
            dir_output=self.dir_output,
        )

        filters = [exact_matches, cross_hybridization, hybridization_probability] #TODO add hybrid. prob
        specificity_filter = SpecificityFilter(filters=filters)
        oligo_database = specificity_filter.apply(
            sequence_type="oligo",
            oligo_database=oligo_database,
            n_jobs=n_jobs,
            reference_database=reference_database,
        )

        # write the intermediate result in a file
        if self.write_intermediate_steps:
            file_database = oligo_database.save_database(filename_out="oligo_database_specificty_filter")
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
        self._log_parameters(parameters)

        num_genes_before, num_oligos_before = self._get_oligo_database_info(oligo_database.database)

        oligos_scoring = PadlockOligoScoring(
            Tm_min=Tm_min,
            Tm_opt=Tm_opt,
            Tm_max=Tm_max,
            GC_content_min=GC_content_min,
            GC_content_opt=GC_content_opt,
            GC_content_max=GC_content_max,
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
            n_sets=n_sets,
            n_jobs=n_jobs,
        )

        # write the intermediate result in a file
        if self.write_intermediate_steps:
            file_database = oligo_database.save_database(filename_out="oligo_database_specificty_filter")
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
    config = initialize_parameters(parser, exp_name="oligo_seq")
    dir_output = os.path.abspath(config["output"])
    Path(dir_output).mkdir(parents=True, exist_ok=True)

    ##### Initialize OligoDesigner Class #####
    oligo_designer = OligoSeq(dir_output=dir_output, log_name="oligo_seq")

    ##### load annotations #####
    oligo_designer.load_annotations(source=config["source"], source_params=config["source_params"])

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

    ##### create oligo database #####
    oligo_database, file_database = oligo_designer.create_oligo_database(
        regions=genes,
        genomic_regions=config["genomic_regions"],
        oligo_length_min=config["oligo_length_min"],
        oligo_length_max=config["oligo_length_max"],
        min_oligos_per_region=config["min_oligos_per_gene"],
        n_jobs=config["n_jobs"],
    )
    print(oligo_database.database["AARS1"]["AARS1::10"])

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
        prohibited_sequence=config["prohibited_sequence"],
        Tm_parameters=config["Tm_parameters"],
        Tm_chem_correction_parameters=config["Tm_chem_correction_parameters"],
        n_jobs=config["n_jobs"],
    )

    ##### filter oligos by specificity #####
    oligo_database, file_database = oligo_designer.filter_by_specificity(
        oligo_database,
        genomic_regions = config["genomic_regions"],
        cross_hybridization_balstn_search_parameters=config["cross_hybridization_balstn_search_parameters"],
        cross_hybridization_balstn_hit_parameters=config["cross_hybridization_balstn_hit_parameters"],
        hybridization_probability_balstn_search_parameters=config["hybridization_probability_balstn_search_parameters"],
        hybridization_probability_balstn_hit_parameters=config["hybridization_probability_balstn_hit_parameters"],
        hybridization_probability_threshold = config['hybridization_probability_threshold'],
        block_size = config["oligo_length_max"] - 1, # only required in case we generate exon_exon_junctions
        n_jobs=config["n_jobs"],
    )

    ##### create oligo sets #####
    oligo_database, file_database, dir_oligosets = oligo_designer.create_oligo_sets(
        oligo_database,
        Tm_min=config["Tm_min"],
        Tm_max=config["Tm_max"],
        Tm_opt=config["Tm_opt"],
        GC_content_min=config["GC_content_min"],
        GC_content_max=config["GC_content_max"],
        GC_content_opt=config["GC_content_opt"],
        oligoset_size=config["oligoset_size"],
        min_oligoset_size=config["min_oligoset_size"],
        max_oligos=config["max_graph_size"],
        n_sets = config["n_sets"],
        n_jobs=config["n_jobs"],
    )

    logging.info(f"Oligo sets were saved in {dir_oligosets}")
    logging.info("##### End of the pipeline. #####")

    # ##### create final padlock sequence #####
    # oligo_designer.create_final_sequences(
    #     oligo_database,
    #     detect_oligo_length_min=config["detect_oligo_length_min"],
    #     detect_oligo_length_max=config["detect_oligo_length_max"],
    #     detect_oligo_Tm_opt=config["detect_oligo_Tm_opt"],
    #     Tm_parameters_detection_oligo=config["Tm_parameters_detection_oligo"],
    #     Tm_chem_correction_param_detection_oligo=config["Tm_chem_correction_param_detection_oligo"],
    # )


if __name__ == "__main__":
    main()
