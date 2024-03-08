import inspect
import logging
import os
import warnings
from argparse import ArgumentParser, RawDescriptionHelpFormatter
from pathlib import Path

from Bio.SeqUtils import MeltingTemp as mt

from oligo_designer_toolsuite.database import OligoDatabase
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
from oligo_designer_toolsuite.pipelines import BaseOligoDesigner
from oligo_designer_toolsuite.pipelines._utils import initialize_parameters


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
        prohibited_sequence: str = "",
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
            oligo_database=oligo_database,
            n_jobs=n_jobs,
        )

        # write the intermediate result in a file
        if self.write_intermediate_steps:
            file_database = oligo_database.save_database(filename_out="oligo_database_property_filter.txt")
        else:
            file_database = ""

        ##### loggig database information #####
        num_genes_after, num_oligos_after = self._get_oligo_database_info(oligo_database.database)
        logging.info(
            f"Step - Filter Oligos by Sequence Property: the database contains {num_oligos_after} oligos from {num_genes_after} genes, while {num_oligos_before - num_oligos_after} oligos and {num_genes_before - num_genes_after} genes have been deleted in this step."
        )

        return oligo_database, file_database

    def filter_by_specificity():
        pass

    def create_oligo_sets():
        pass

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
    print("hi")
    config = initialize_parameters(parser, exp_name="oligo_seq")
    print("hi")
    dir_output = os.path.abspath(config["output"])
    Path(dir_output).mkdir(parents=True, exist_ok=True)

    ##### Initialize OligoDesigner Class #####
    oligo_designer = OligoSeq(dir_output=dir_output, log_name="oligo_seq")

    ##### load annotations #####
    oligo_designer.load_annotations(source=config["source"], source_params=config["source_params"])

    ##### read the genes file #####
    if config["file_genes"] is None:
        warnings.warn(
            "No gene list file was provided! All genes from fasta file are used to generate the oligos. This chioce can use a lot of resources."
        )
        genes = None
    else:
        with open(config["file_genes"]) as handle:
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
        Tm_chem_correction_parameters=config["Tm_chem_correction_param"],
        n_jobs=config["n_jobs"],
    )

    ##### filter oligos by specificity #####
    oligo_database, file_database = oligo_designer.filter_by_specificity(
        oligo_database,
        ligation_region_size=config["ligation_region_size"],
        blast_word_size=config["blast_word_size"],
        blast_percent_identity=config["blast_percent_identity"],
        blast_coverage=config["blast_coverage"],
        n_jobs=config["n_jobs"],
    )

    ##### create oligo sets #####
    oligo_database, file_database, dir_oligosets = oligo_designer.create_oligo_sets(
        oligo_database,
        oligoset_size_opt=config["oligoset_size_opt"],
        oligoset_size_min=config["oligoset_size_min"],
        n_sets=config["n_sets"],
        Tm_min=config["Tm_min"],
        Tm_max=config["Tm_max"],
        Tm_opt=config["Tm_opt"],
        Tm_weight=config["Tm_weight"],
        GC_content_min=config["GC_content_min"],
        GC_content_max=config["GC_content_max"],
        GC_content_opt=config["GC_content_opt"],
        GC_weight=config["GC_weight"],
        n_jobs=config["n_jobs"],
        max_oligos=config["max_graph_size"],
    )

    ##### create final padlock sequence #####
    oligo_designer.create_final_sequences(
        oligo_database,
        detect_oligo_length_min=config["detect_oligo_length_min"],
        detect_oligo_length_max=config["detect_oligo_length_max"],
        detect_oligo_Tm_opt=config["detect_oligo_Tm_opt"],
        Tm_parameters_detection_oligo=config["Tm_parameters_detection_oligo"],
        Tm_chem_correction_param_detection_oligo=config["Tm_chem_correction_param_detection_oligo"],
    )


if __name__ == "__main__":
    main()
