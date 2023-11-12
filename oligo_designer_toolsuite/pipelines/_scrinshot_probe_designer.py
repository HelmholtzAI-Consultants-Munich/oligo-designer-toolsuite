############################################
# imports
############################################

import os
import sys
import yaml
import shutil
import logging
import inspect
import warnings
import gc
import psutil

from pathlib import Path
from datetime import datetime
from subprocess import Popen
from argparse import ArgumentParser, RawDescriptionHelpFormatter
from pympler import summary, muppy, tracker



from Bio.SeqUtils import MeltingTemp as mt

from ._utils import initialize_parameters

from oligo_designer_toolsuite.database import (
    CustomGenomicRegionGenerator,
    NcbiGenomicRegionGenerator,
    EnsemblGenomicRegionGenerator,
    OligoDatabase,
    ReferenceDatabase,
)
from oligo_designer_toolsuite.sequence_design import PadlockSequence
from oligo_designer_toolsuite.oligo_efficiency_filter import (
    PadlockOligoScoring,
    PadlockSetScoring,
)
from oligo_designer_toolsuite.oligo_property_filter import (
    GCContent,
    MaskedSequences,
    MeltingTemperatureNN,
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

from ._base_probe_designer import BaseProbeDesigner


############################################
# Scrinshot probe design class
############################################


class ScrinshotProbeDesigner(BaseProbeDesigner):
    """This class generates all padlock probes for a SCRINSHOT experiment from a transcriptome or custom file for a user-defined set of genes.
    The probe design is done in five steps:
    1. Creating all possible probes for a provided annotation and set of genes and store them in a oligo database
    2. Filter probes by a list of property filters, e.g. CG content filter
    3. Filter probes by specificity against a reference database (e.g. transcriptome)
    4. Select sets of best scoring, non-overlappign oligos for each gene
    5. Create the final ready-to-order padlock sequence

    The user can save the oligo database after each processing step and resume the pipeline later by loading an existing database into this class.
    A logger is automatically created at <output_dir>/log_scrinshot_probe_designer_<timestamp>.txt and logs parameters as well as number of genes/oligos after each step.

    :param dir_output: Output directory, defaults to 'output'.
    :type dir_output: str, optional
    :param write_removed_genes: write removed regions to file ``regions_with_insufficient_oligos.txt``, defaults to True
    :type write_removed_genes: bool, optional
    :param write_intermediate_steps: save oligo database after each processing step, defaults to True
    :type write_intermediate_steps: bool, optional
    """

    def filter_probes_by_property(
        self,
        probe_database,
        GC_content_min: int = 40,
        GC_content_max: int = 60,
        Tm_min: int = 52,
        Tm_max: int = 67,
        min_arm_length: int = 10,
        max_arm_Tm_dif: int = 2,
        arm_Tm_min: int = 38,
        arm_Tm_max: int = 49,
        Tm_parameters_probe: dict = {
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
            "Na": 1.25,
            "K": 75,
            "Tris": 20,
            "Mg": 10,
        },
        Tm_chem_correction_param_probe: dict = {
            "DMSO": 0,
            "DMSOfactor": 0.75,
            "fmdfactor": 0.65,
            "fmdmethod": 1,
            "GC": None,
            "fmd": 20,
        },
        n_jobs: int = 1,
    ):
        ##### log parameters #####
        logging.info("Parameters Property Filters:")
        args, _, _, values = inspect.getargvalues(inspect.currentframe())
        parameters = {i: values[i] for i in args}
        self._log_parameters(parameters)

        num_genes_before, num_probes_before = self._get_probe_database_info(
            probe_database.database
        )

        ##### preprocess melting temperature params #####
        Tm_parameters_probe["nn_table"] = getattr(mt, Tm_parameters_probe["nn_table"])
        Tm_parameters_probe["tmm_table"] = getattr(mt, Tm_parameters_probe["tmm_table"])
        Tm_parameters_probe["imm_table"] = getattr(mt, Tm_parameters_probe["imm_table"])
        Tm_parameters_probe["de_table"] = getattr(mt, Tm_parameters_probe["de_table"])

        ##### initialize the filters classes #####
        masked_sequences = MaskedSequences()
        gc_content = GCContent(
            GC_content_min=GC_content_min, GC_content_max=GC_content_max
        )
        melting_temperature = MeltingTemperatureNN(
            Tm_min=Tm_min,
            Tm_max=Tm_max,
            Tm_parameters=Tm_parameters_probe,
            Tm_chem_correction_parameters=Tm_chem_correction_param_probe,
        )
        padlock_arms = PadlockArms(
            min_arm_length=min_arm_length,
            max_arm_Tm_dif=max_arm_Tm_dif,
            arm_Tm_min=arm_Tm_min,
            arm_Tm_max=arm_Tm_max,
            Tm_parameters=Tm_parameters_probe,
            Tm_chem_correction_parameters=Tm_chem_correction_param_probe,
        )

        ##### apply property filter to the database #####
        filters = [masked_sequences, gc_content, melting_temperature, padlock_arms]
        property_filter = PropertyFilter(filters=filters)
        probe_database = property_filter.apply(
            oligo_database=probe_database, n_jobs=n_jobs
        )

        ##### save database #####
        if self.write_intermediate_steps:
            file_database = probe_database.write_database(
                filename="probe_database_property_filter.txt"
            )
        else:
            file_database = ""

        ##### loggig database information #####
        num_genes_after, num_probes_after = self._get_probe_database_info(
            probe_database.database
        )
        logging.info(
            f"Step - Filter Probes by Sequence Property: the database contains {num_probes_after} probes from {num_genes_after} genes, while {num_probes_before - num_probes_after} probes and {num_genes_before - num_genes_after} genes have been deleted in this step."
        )

        return probe_database, file_database

    def filter_probes_by_specificity(
        self,
        probe_database,
        ligation_region_size: int = 5,
        blast_word_size: int = 10,
        blast_percent_identity: int = 80,
        blast_coverage: int = 50,
        n_jobs: int = 1,
    ):
        dir_specificity = os.path.join(self.dir_output, "specificity_temporary")

        ##### log parameters #####
        logging.info("Parameters Specificity Filters:")
        args, _, _, values = inspect.getargvalues(inspect.currentframe())
        parameters = {i: values[i] for i in args}
        self._log_parameters(parameters)

        num_genes_before, num_probes_before = self._get_probe_database_info(
            probe_database.database
        )

        ##### generate transcriptome for reference #####
        # length of exon_junction_size is longer than probe length to cover bulges in alignments
        (
            probe_length_min,
            probe_length_max,
        ) = self._get_probe_length_min_max_from_database(probe_database.database)
        file_transcriptome = (
            self.region_generator.generate_transcript_reduced_representation(
                include_exon_junctions=True, exon_junction_size=2 * probe_length_max
            )
        )
        reference_database = ReferenceDatabase(
            file_fasta=file_transcriptome,
            metadata=self.metadata,
            dir_output=self.dir_output,
        )

        ##### intialize the filter classes #####
        exact_mathces = ExactMatches(dir_specificity=dir_specificity)
        seed_ligation = LigationRegionCreation(
            ligation_region_size=ligation_region_size
        )
        seed_region = BowtieSeedRegion(
            dir_specificity=dir_specificity,
            seed_region_creation=seed_ligation,
            strand="plus",
        )
        blastn = Blastn(
            dir_specificity=dir_specificity,
            word_size=blast_word_size,
            percent_identity=blast_percent_identity,
            coverage=blast_coverage,
            strand="plus",
        )

        ##### apply specificity filter to the database #####
        filters = [exact_mathces, blastn]  # seed_region
        specificity_filter = SpecificityFilter(filters=filters)
        probe_database = specificity_filter.apply(
            oligo_database=probe_database,
            reference_database=reference_database,
            n_jobs=n_jobs,
        )

        ##### save database #####
        if self.write_intermediate_steps:
            file_database = probe_database.write_database(
                filename="probe_database_specificity_filters.txt"
            )
        else:
            file_database = ""

        ##### loggig database information #####
        num_genes_after, num_probes_after = self._get_probe_database_info(
            probe_database.database
        )
        logging.info(
            f"Step - Filter Probes by Specificity: the database contains {num_probes_after} probes from {num_genes_after} genes, while {num_probes_before - num_probes_after} probes and {num_genes_before - num_genes_after} genes have been deleted in this step."
        )

        shutil.rmtree(dir_specificity)

        return probe_database, file_database

    def create_probe_sets(
        self,
        probe_database,
        probeset_size_opt: int = 5,
        probeset_size_min: int = 2,
        n_sets: int = 100,
        Tm_min: int = 52,
        Tm_max: int = 67,
        Tm_opt: int = 60,
        Tm_weight: int = 1,
        GC_content_min: int = 40,
        GC_content_max: int = 50,
        GC_content_opt: int = 60,
        GC_weight: int = 1,
        n_jobs: int = 1,
        max_oligos: int = 5000,
    ):
        ##### log parameters #####
        logging.info("Parameters Probesets:")
        args, _, _, values = inspect.getargvalues(inspect.currentframe())
        parameters = {i: values[i] for i in args}
        self._log_parameters(parameters)

        num_genes_before, num_probes_before = self._get_probe_database_info(
            probe_database.database
        )

        ##### initialize the scoring and oligoset generator classes #####
        set_scoring = PadlockSetScoring()
        oligos_scoring = PadlockOligoScoring(
            Tm_min=Tm_min,
            Tm_opt=Tm_opt,
            Tm_max=Tm_max,
            GC_content_min=GC_content_min,
            GC_content_opt=GC_content_opt,
            GC_content_max=GC_content_max,
            Tm_weight=Tm_weight,
            GC_weight=GC_weight,
        )
        oligoset_generator = OligosetGenerator(
            oligoset_size=probeset_size_opt,
            min_oligoset_size=probeset_size_min,
            oligos_scoring=oligos_scoring,
            set_scoring=set_scoring,
            heurustic_selection=padlock_heuristic_selection,
            max_oligos=max_oligos,
        )

        ##### generate the oligoset #####
        probe_database = oligoset_generator.apply(
            oligo_database=probe_database, n_sets=n_sets, n_jobs=n_jobs
        )

        ##### save database #####
        if self.write_intermediate_steps:
            dir_oligosets = probe_database.write_oligosets(folder="oligosets")
            file_database = probe_database.write_database(
                filename="probe_database_oligosets.txt"
            )
        else:
            dir_oligosets = ""
            file_database = ""

        ##### loggig database information #####
        num_genes_after, num_probes_after = self._get_probe_database_info(
            probe_database.database
        )
        logging.info(
            f"Step - Generate Oligosets: the database contains {num_probes_after} probes from {num_genes_after} genes, while {num_probes_before - num_probes_after} probes and {num_genes_before - num_genes_after} genes have been deleted in this step."
        )

        return probe_database, file_database, dir_oligosets

    def create_final_sequences(
        self,
        probe_database,
        detect_oligo_length_min: int = 18,
        detect_oligo_length_max: int = 25,
        detect_oligo_Tm_opt: int = 32,
        Tm_parameters_detection_oligo: dict = {
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
            "Na": 39,
            "K": 0,
            "Tris": 0,
            "Mg": 0,
        },
        Tm_chem_correction_param_detection_oligo: dict = {
            "DMSO": 0,
            "DMSOfactor": 0.75,
            "fmdfactor": 0.65,
            "fmdmethod": 1,
            "GC": None,
            "fmd": 30,
        },
    ):
        """Generates the padlock sequences for a OligoDataset class for which oligosets have been already computed."""

        ##### log parameters #####
        logging.info("Parameters Final Sequence Design:")
        args, _, _, values = inspect.getargvalues(inspect.currentframe())
        parameters = {i: values[i] for i in args}
        self._log_parameters(parameters)

        ##### preprocessing of the melting temperature parameters #####
        Tm_parameters_detection_oligo["nn_table"] = getattr(
            mt, Tm_parameters_detection_oligo["nn_table"]
        )
        Tm_parameters_detection_oligo["tmm_table"] = getattr(
            mt, Tm_parameters_detection_oligo["tmm_table"]
        )
        Tm_parameters_detection_oligo["imm_table"] = getattr(
            mt, Tm_parameters_detection_oligo["imm_table"]
        )
        Tm_parameters_detection_oligo["de_table"] = getattr(
            mt, Tm_parameters_detection_oligo["de_table"]
        )

        ##### initilize the padlock sequence designer class #####
        padlock_sequence = PadlockSequence(
            detect_oligo_length_min=detect_oligo_length_min,
            detect_oligo_length_max=detect_oligo_length_max,
            detect_oligo_Tm_opt=detect_oligo_Tm_opt,
            Tm_parameters=Tm_parameters_detection_oligo,
            Tm_chem_correction_parameters=Tm_chem_correction_param_detection_oligo,
            dir_output=self.dir_output,
        )

        ##### generate the final padlock sequence #####
        padlock_sequence.design_final_padlock_sequence(oligo_database=probe_database)
        logging.info(
            f"Step - Design Final Padlock Sequences: padlock sequences are stored in '{os.path.join(padlock_sequence.dir_output, 'padlock_sequences')}' directory."
        )


def main():
    """Command line tool to run a pipeline to design Padlock Probes for SCRINSHOT experiments. To run the tool use the command: ``scrinshot_probe_designer [options]``.

    The program supports two ways to recieve the input parameters:

    - command line input
    - configuration file (recommended)

    A standard configuration file can be generated using the following command ``scrinshot_probe_designer_config [options]``

    Since the number of input parameters requested is too high to be handled only through the command line,
    the progam will use as baseline the configuration file (it will be automatically generated if it wasn't provided)
    and the parameters specified in the command line will be overwritten with the values given.

    REMARK: melting temperature parameters can be given only through the configuration file.
    """
    # memory_tracker = tracker.SummaryTracker()
    # get comman line arguments
    parser = ArgumentParser(
        prog="SCRINSHOT Probe Designer",
        usage="scrinshot_probe_designer [options]",
        description=__doc__,
        formatter_class=RawDescriptionHelpFormatter,
    )

    config = initialize_parameters(parser, exp_name="scrinshot")

    dir_output = os.path.abspath(config["output"])
    Path(dir_output).mkdir(parents=True, exist_ok=True)

    ##### Initialize ProbeDesigner Class #####
    probe_designer = ScrinshotProbeDesigner(dir_output=dir_output)

    ##### load annotations #####
    probe_designer.load_annotations(
        source=config["source"], source_params=config["source_params"]
    )
    # print('\n\n')
    # memory_tracker.print_diff()
    # print(f"annotaion class : {pympler.asizeof.asizeof(probe_designer.region_generator)/1000000000}")
    # print('\n\n')

    ##### read the genes file #####
    if config["file_genes"] is None:
        warnings.warn(
            "No gene list file was provided! All genes from fasta file are used to generate the probes. This chioce can use a lot of resources."
        )
        genes = None
    else:
        with open(config["file_genes"]) as handle:
            lines = handle.readlines()
            genes = [line.rstrip() for line in lines]

    ##### create probe database #####
    probe_database, file_database = probe_designer.create_probe_database(
        genes=genes,
        probe_length_min=config["probe_length_min"],
        region=config["region"],
        probe_length_max=config["probe_length_max"],
        min_probes_per_gene=config["min_probes_per_gene"],
        n_jobs=config["n_jobs"],
    )
    # print('\n\n')
    # memory_tracker.print_diff()
    # print(f"annotaion class : {pympler.asizeof.asizeof(probe_designer.region_generator)/1000000000}")
    # print(f"database class : {pympler.asizeof.asizeof(probe_database)/1000000000}")
    # print('\n\n')


    ##### filter probes by property #####
    probe_database, file_database = probe_designer.filter_probes_by_property(
        probe_database,
        GC_content_min=config["target_GC_content_min"],
        GC_content_max=config["target_GC_content_max"],
        Tm_min=config["target_Tm_min"],
        Tm_max=config["target_Tm_max"],
        min_arm_length=config["min_arm_length"],
        max_arm_Tm_dif=config["max_arm_Tm_dif"],
        arm_Tm_min=config["arm_Tm_min"],
        arm_Tm_max=config["arm_Tm_max"],
        Tm_parameters_probe=config["Tm_parameters_probe"],
        Tm_chem_correction_param_probe=config["Tm_chem_correction_param_probe"],
        n_jobs=config["n_jobs"],
    )

    # print('\n\n')
    # memory_tracker.print_diff()      
    # print(f"annotaion class : {pympler.asizeof.asizeof(probe_designer.region_generator)/1000000000}")
    # print(f"database class : {pympler.asizeof.asizeof(probe_database)/1000000000}")
    # print('\n\n')

    ##### filter probes by specificity #####
    probe_database, file_database = probe_designer.filter_probes_by_specificity(
        probe_database,
        ligation_region_size=config["ligation_region_size"],
        blast_word_size=config["blast_word_size"],
        blast_percent_identity=config["blast_percent_identity"],
        blast_coverage=config["blast_coverage"],
        n_jobs=config["n_jobs"],
    )
    # print('\n\n')
    # memory_tracker.print_diff()
    # print('\n\nRAM Used (GB) (specificity):', psutil.virtual_memory()[3]/1000000000,  '\n\n')
    # print(f"annotaion class : {pympler.asizeof.asizeof(probe_designer.region_generator)/1000000000}")
    # print(f"database class : {pympler.asizeof.asizeof(probe_database)/1000000000}")
    # print('\n\n')


    ##### create probe sets #####
    probe_database, file_database, dir_oligosets = probe_designer.create_probe_sets(
        probe_database,
        probeset_size_opt=config["probeset_size_opt"],
        probeset_size_min=config["probeset_size_min"],
        n_sets=config["n_sets"],
        Tm_min=config["target_Tm_min"],
        Tm_max=config["target_Tm_max"],
        Tm_opt=config["Tm_opt"],
        Tm_weight=config["Tm_weight"],
        GC_content_min=config["target_GC_content_min"],
        GC_content_max=config["target_GC_content_max"],
        GC_content_opt=config["GC_content_opt"],
        GC_weight=config["GC_weight"],
        n_jobs=config["n_jobs"],
        max_oligos = config["max_graph_size"],
    )
    # print('\n\n')
    # memory_tracker.print_diff()
    # print('\n\nRAM Used (GB) (probesets):', psutil.virtual_memory()[3]/1000000000,  '\n\n')
    # print(f"annotaion class : {pympler.asizeof.asizeof(probe_designer.region_generator)/1000000000}")
    # print(f"database class : {pympler.asizeof.asizeof(probe_database)/1000000000}")
    # print('\n\n')

    ##### create final padlock sequence #####
    probe_designer.create_final_sequences(
        probe_database,
        detect_oligo_length_min=config["detect_oligo_length_min"],
        detect_oligo_length_max=config["detect_oligo_length_max"],
        detect_oligo_Tm_opt=config["detect_oligo_Tm_opt"],
        Tm_parameters_detection_oligo=config["Tm_parameters_detection_oligo"],
        Tm_chem_correction_param_detection_oligo=config[
            "Tm_chem_correction_param_detection_oligo"
        ],
    )


if __name__ == "__main__":
    main()
