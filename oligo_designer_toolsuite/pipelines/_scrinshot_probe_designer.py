############################################
# imports
############################################

import os
import yaml
import random
import itertools
import inspect
import logging
import warnings
import pandas as pd

from pathlib import Path
from datetime import datetime

from typing import List

from Bio.SeqUtils import MeltingTemp as mt

from oligo_designer_toolsuite.pipelines._utils import (
    log_parameters,
    base_parser,
    get_oligo_database_info,
    generation_step,
    filtering_step,
)
from oligo_designer_toolsuite.sequence_generator import OligoSequenceGenerator
from oligo_designer_toolsuite.database import OligoDatabase, ReferenceDatabase, OligoAttributes
from oligo_designer_toolsuite.oligo_property_filter import (
    GCContentFilter,
    HardMaskedSequenceFilter,
    HomopolymericRunsFilter,
    MeltingTemperatureNNFilter,
    PropertyFilter,
    SoftMaskedSequenceFilter,
    PadlockArmsFilter,
)
from oligo_designer_toolsuite.oligo_specificity_filter import (
    ExactMatchFilter,
    BlastNFilter,
    BlastNSeedregionLigationsiteFilter,
    CrossHybridizationFilter,
    RemoveByLargerRegionPolicy,
    SpecificityFilter,
)
from oligo_designer_toolsuite.oligo_selection import (
    OligosetGeneratorIndependentSet,
    heuristic_selection_independent_set,
)
from oligo_designer_toolsuite.oligo_efficiency_filter import WeightedIsoformTmGCOligoScoring, LowestSetScoring


############################################
# Oligo-seq Designer Functions
############################################


class ScrinshotProbeDesigner:

    def __init__(
        self, file_regions: list, write_intermediate_steps: bool, dir_output: str, n_jobs: int
    ) -> None:
        """Constructor for the ScrinshotProbeDesigner class."""
        ##### read the genes file #####
        if file_regions is None:
            warnings.warn(
                "No gene list file was provided! All genes from fasta file are used to generate the probes. This chioce can use a lot of resources."
            )
            self.gene_ids = None
        else:
            with open(file_regions) as handle:
                lines = handle.readlines()
                self.gene_ids = [line.rstrip() for line in lines]

        self.write_intermediate_steps = write_intermediate_steps

        ##### create the output folder #####
        self.dir_output = os.path.abspath(dir_output)
        Path(self.dir_output).mkdir(parents=True, exist_ok=True)

        self.subdir_db_probes = "db_probes"
        self.subdir_db_reference = "db_reference"

        self.n_jobs = n_jobs

        self.probe_attributes = OligoAttributes()

        ##### setup logger #####
        timestamp = datetime.now()
        file_logger = os.path.join(
            self.dir_output,
            f"log_scrinshot_probe_designer_{timestamp.year}-{timestamp.month}-{timestamp.day}-{timestamp.hour}-{timestamp.minute}.txt",
        )
        logging.getLogger("log_name")
        logging.basicConfig(
            format="%(asctime)s [%(levelname)s] %(message)s",
            level=logging.NOTSET,
            handlers=[logging.FileHandler(file_logger), logging.StreamHandler()],
        )
        logging.captureWarnings(True)

    @generation_step(step_name="Create Database")
    def create_probe_database(
        self,
        probe_length_min: int,
        probe_length_max: int,
        files_fasta_oligo_database: list[str],
        probeset_size_min: int,
    ):
        ##### creating the probe sequences #####
        probe_sequences = OligoSequenceGenerator(dir_output=self.dir_output)
        probe_fasta_file = probe_sequences.create_sequences_sliding_window(
            filename_out="probe_sequences",
            files_fasta_in=files_fasta_oligo_database,
            length_interval_sequences=(probe_length_min, probe_length_max),
            region_ids=self.gene_ids,
            n_jobs=self.n_jobs,
        )

        ##### creating the probe database #####
        oligo_database = OligoDatabase(
            min_oligos_per_region=probeset_size_min,
            write_regions_with_insufficient_oligos=True,
            database_name=self.subdir_db_probes,
            dir_output=self.dir_output,
        )
        oligo_database.load_sequences_from_fasta(
            files_fasta=[probe_fasta_file],
            sequence_type="target",
            region_ids=self.gene_ids,
        )

        ##### save database #####
        if self.write_intermediate_steps:
            file_database = oligo_database.save_database(filename="1_db_probes_initial")
        else:
            file_database = ""

        return oligo_database, file_database

    @filtering_step(step_name="Property Filters")
    def filter_by_property(
        self,
        oligo_database: OligoDatabase,
        probe_GC_content_min: float,
        probe_GC_content_max: float,
        probe_Tm_min: float,
        probe_Tm_max: float,
        arm_Tm_dif_max: int,
        arm_length_min: int,
        arm_Tm_min: float,
        arm_Tm_max: float,
        homopolymeric_base_n: str,
        Tm_parameters_probe: dict,
        Tm_chem_correction_param_probe: dict,
    ):
        num_genes_before, num_probes_before = get_oligo_database_info(oligo_database.database)

        # define the filters
        hard_masked_sequences = HardMaskedSequenceFilter()
        soft_masked_sequences = SoftMaskedSequenceFilter()
        gc_content = GCContentFilter(GC_content_min=probe_GC_content_min, GC_content_max=probe_GC_content_max)
        melting_temperature = MeltingTemperatureNNFilter(
            Tm_min=probe_Tm_min,
            Tm_max=probe_Tm_max,
            Tm_parameters=Tm_parameters_probe,
            Tm_chem_correction_parameters=Tm_chem_correction_param_probe,
        )
        homopolymeric_runs = HomopolymericRunsFilter(
            base_n=homopolymeric_base_n,
        )
        padlock_arms_filter = PadlockArmsFilter(
            arm_length_min=arm_length_min,
            arm_Tm_dif_max=arm_Tm_dif_max,
            arm_Tm_min=arm_Tm_min,
            arm_Tm_max=arm_Tm_max,
            Tm_parameters=Tm_parameters_probe,
            Tm_chem_correction_parameters=Tm_chem_correction_param_probe,
        )

        filters = [
            hard_masked_sequences,
            soft_masked_sequences,
            gc_content,
            melting_temperature,
            homopolymeric_runs,
            padlock_arms_filter,
        ]

        # initialize the preoperty filter class
        property_filter = PropertyFilter(filters=filters)

        # filter the database
        oligo_database = property_filter.apply(
            sequence_type="oligo",
            oligo_database=oligo_database,
            n_jobs=self.n_jobs,
        )

        # write the intermediate result in a file
        if self.write_intermediate_steps:
            file_database = oligo_database.save_database(filename="2_db_probes_property_filter")
        else:
            file_database = ""

        return oligo_database, file_database

    @filtering_step(step_name="Specificty Filters")
    def filter_by_specificity(
        self,
        oligo_database: OligoDatabase,
        files_fasta_reference_database: List[str],
        specificity_blastn_search_parameters: dict,
        specificity_blastn_hit_parameters: dict,
        cross_hybridization_blastn_search_parameters: dict,
        cross_hybridization_blastn_hit_parameters: dict,
        ligation_region_size: int,
        arm_Tm_dif_max: int,
        arm_length_min: int,
        arm_Tm_min: float,
        arm_Tm_max: float,
        Tm_parameters_probe: dict,
        Tm_chem_correction_param_probe: dict,
    ):
        ##### define reference database #####
        reference_database = ReferenceDatabase(
            database_name=self.subdir_db_reference, dir_output=self.dir_output
        )
        reference_database.load_sequences_from_fasta(
            files_fasta=files_fasta_reference_database, database_overwrite=False
        )

        ##### calculate required probe attributes #####
        oligo_database = self.probe_attributes.calculate_padlock_arms(
            oligo_database=oligo_database,
            sequence_type="oligo",
            arm_length_min=arm_length_min,
            arm_Tm_dif_max=arm_Tm_dif_max,
            arm_Tm_min=arm_Tm_min,
            arm_Tm_max=arm_Tm_max,
            Tm_parameters=Tm_parameters_probe,
            Tm_chem_correction_parameters=Tm_chem_correction_param_probe,
        )

        ##### specificity filters #####
        # removing duplicated probes from the region with the most probes
        exact_matches_policy = RemoveByLargerRegionPolicy()
        exact_matches = ExactMatchFilter(policy=exact_matches_policy)

        cross_hybridization_aligner = BlastNFilter(
            search_parameters=cross_hybridization_blastn_search_parameters,
            hit_parameters=cross_hybridization_blastn_hit_parameters,
            filter_name="blastn_crosshybridization",
            dir_output=self.dir_output,
        )
        cross_hybridization_policy = RemoveByLargerRegionPolicy()
        cross_hybridization = CrossHybridizationFilter(
            policy=cross_hybridization_policy,
            alignment_method=cross_hybridization_aligner,
            database_name_reference=self.subdir_db_reference,
            dir_output=self.dir_output,
        )

        if ligation_region_size > 0:
            specificity = BlastNSeedregionLigationsiteFilter(
                seedregion_size=ligation_region_size,
                search_parameters=specificity_blastn_search_parameters,
                hit_parameters=specificity_blastn_hit_parameters,
                filter_name="blastn_specificity",
                dir_output=self.dir_output,
            )
        else:
            specificity = BlastNFilter(
                search_parameters=specificity_blastn_search_parameters,
                hit_parameters=specificity_blastn_hit_parameters,
                filter_name="blastn_specificity",
                dir_output=self.dir_output,
            )

        filters = [exact_matches, cross_hybridization, specificity]
        specificity_filter = SpecificityFilter(filters=filters)
        oligo_database = specificity_filter.apply(
            sequence_type="oligo",
            oligo_database=oligo_database,
            reference_database=reference_database,
            n_jobs=self.n_jobs,
        )

        # write the intermediate result in a file
        if self.write_intermediate_steps:
            file_database = oligo_database.save_database(filename="3_db_probes_specificty_filter")
        else:
            file_database = ""

        dir = reference_database.dir_output
        os.rmdir(dir) if os.path.exists(dir) else None

        dir = cross_hybridization_aligner.dir_output
        os.rmdir(dir) if os.path.exists(dir) else None

        dir = cross_hybridization.dir_output
        os.rmdir(dir) if os.path.exists(dir) else None

        dir = specificity.dir_output
        os.rmdir(dir) if os.path.exists(dir) else None

        return oligo_database, file_database

    def compute_probe_attributes(
        self,
        oligo_database: OligoDatabase,
        Tm_parameters_probe: dict,
        Tm_chem_correction_param_probe: dict,
    ):
        """
        Computes various attributes for probes and stores them in the provided probe database.

        :param oligo_database: The database containing probes for which attributes are to be computed.
        :type oligo_database: OligoDatabase
        :param secondary_structures_T: Temperature at which secondary structure delta G is calculated.
        :type secondary_structures_T: float
        :param Tm_parameters: Parameters to calculate melting temperature.
        :type Tm_parameters: dict
        :param Tm_chem_correction_param_probe: Chemical correction parameters for melting temperature calculation.
        :type Tm_chem_correction_param_probe: dict
        :return: Updated probe database with new attributes.
        :rtype: OligoDatabase
        """

        oligo_database = self.probe_attributes.calculate_oligo_length(oligo_database=oligo_database)
        oligo_database = self.probe_attributes.calculate_GC_content(
            oligo_database=oligo_database, sequence_type="oligo"
        )
        oligo_database = self.probe_attributes.calculate_TmNN(
            oligo_database=oligo_database,
            sequence_type="oligo",
            Tm_parameters=Tm_parameters_probe,
            Tm_chem_correction_parameters=Tm_chem_correction_param_probe,
        )
        oligo_database = self.probe_attributes.calculate_num_targeted_transcripts(
            oligo_database=oligo_database
        )
        oligo_database = self.probe_attributes.calculate_isoform_consensus(oligo_database=oligo_database)

        return oligo_database

    @filtering_step(step_name="Set Selection")
    def create_probe_sets(
        self,
        oligo_database: OligoDatabase,
        probe_isoform_weight: float,
        probe_Tm_weight: float,
        probe_Tm_min: float,
        probe_Tm_opt: float,
        probe_Tm_max: float,
        Tm_parameters_probe: dict,
        Tm_chem_correction_param_probe: dict,
        probe_GC_weight: float,
        probe_GC_content_min: float,
        probe_GC_content_opt: float,
        probe_GC_content_max: float,
        probeset_size_opt: int,
        probeset_size_min: int,
        max_graph_size: int,
        n_sets: int,
        distance_between_probes: int,
    ):
        probes_scoring = WeightedIsoformTmGCOligoScoring(
            Tm_min=probe_Tm_min,
            Tm_opt=probe_Tm_opt,
            Tm_max=probe_Tm_max,
            GC_content_min=probe_GC_content_min,
            GC_content_opt=probe_GC_content_opt,
            GC_content_max=probe_GC_content_max,
            Tm_parameters=Tm_parameters_probe,
            Tm_chem_correction_parameters=Tm_chem_correction_param_probe,
            isoform_weight=probe_isoform_weight,
            Tm_weight=probe_Tm_weight,
            GC_weight=probe_GC_weight,
        )
        set_scoring = LowestSetScoring()
        probeset_generator = OligosetGeneratorIndependentSet(
            opt_oligoset_size=probeset_size_opt,
            min_oligoset_size=probeset_size_min,
            oligos_scoring=probes_scoring,
            set_scoring=set_scoring,
            heurustic_selection=heuristic_selection_independent_set,
            max_oligos=max_graph_size,
            distance_between_oligos=distance_between_probes,
        )
        oligo_database = probeset_generator.apply(
            oligo_database=oligo_database,
            sequence_type="oligo",
            n_sets=n_sets,
            n_jobs=self.n_jobs,
        )

        # write the intermediate result in a file
        if self.write_intermediate_steps:
            file_database = oligo_database.save_database(filename="4_db_probes_probesets")
            file_probesets = oligo_database.write_oligosets()
        else:
            file_database = ""
            file_probesets = ""

        return oligo_database, file_database, file_probesets


############################################
# Oligo-seq Designer Pipeline
############################################


def main():

    args = base_parser()

    ##### read the config file #####
    with open(args["config"], "r") as handle:
        config = yaml.safe_load(handle)

    ##### initialize probe designer pipeline #####
    pipeline = ScrinshotProbeDesigner(
        file_regions=config["file_regions"],
        write_intermediate_steps=config["write_intermediate_steps"],
        dir_output=config["dir_output"],
        n_jobs=config["n_jobs"],
    )

    ##### create probe database #####
    probe_database, file_database = pipeline.create_probe_database(
        probe_length_min=config["probe_length_min"],
        probe_length_max=config["probe_length_max"],
        files_fasta_oligo_database=config["files_fasta_probe_database"],
        # we should have at least "min_probeset_size" probes per gene to create one set
        probeset_size_min=config["probeset_size_min"],
    )

    ##### preprocess melting temperature params #####
    Tm_parameters_probe = config["Tm_parameters_probe"]
    Tm_parameters_probe["nn_table"] = getattr(mt, Tm_parameters_probe["nn_table"])
    Tm_parameters_probe["tmm_table"] = getattr(mt, Tm_parameters_probe["tmm_table"])
    Tm_parameters_probe["imm_table"] = getattr(mt, Tm_parameters_probe["imm_table"])
    Tm_parameters_probe["de_table"] = getattr(mt, Tm_parameters_probe["de_table"])

    Tm_parameters_detection_oligo = config["Tm_parameters_detection_oligo"]
    Tm_parameters_detection_oligo["nn_table"] = getattr(mt, Tm_parameters_detection_oligo["nn_table"])
    Tm_parameters_detection_oligo["tmm_table"] = getattr(mt, Tm_parameters_detection_oligo["tmm_table"])
    Tm_parameters_detection_oligo["imm_table"] = getattr(mt, Tm_parameters_detection_oligo["imm_table"])
    Tm_parameters_detection_oligo["de_table"] = getattr(mt, Tm_parameters_detection_oligo["de_table"])

    ##### filter probes by property #####
    oligo_database, file_database = pipeline.filter_by_property(
        oligo_database=probe_database,
        probe_GC_content_min=config["probe_GC_content_min"],
        probe_GC_content_max=config["probe_GC_content_max"],
        probe_Tm_min=config["probe_Tm_min"],
        probe_Tm_max=config["probe_Tm_max"],
        arm_Tm_dif_max=config["arm_Tm_dif_max"],
        arm_length_min=config["arm_length_min"],
        arm_Tm_min=config["arm_Tm_min"],
        arm_Tm_max=config["arm_Tm_max"],
        homopolymeric_base_n=config["homopolymeric_base_n"],
        Tm_parameters_probe=Tm_parameters_probe,
        Tm_chem_correction_param_probe=config["Tm_chem_correction_param_probe"],
    )

    # ##### filter probes by specificity #####
    probe_database, file_database = pipeline.filter_by_specificity(
        oligo_database=probe_database,
        files_fasta_reference_database=config["files_fasta_reference_database"],
        specificity_blastn_search_parameters=config["specificity_blastn_search_parameters"],
        specificity_blastn_hit_parameters=config["specificity_blastn_hit_parameters"],
        cross_hybridization_blastn_search_parameters=config["cross_hybridization_blastn_search_parameters"],
        cross_hybridization_blastn_hit_parameters=config["cross_hybridization_blastn_hit_parameters"],
        ligation_region_size=config["ligation_region_size"],
        arm_Tm_dif_max=config["arm_Tm_dif_max"],
        arm_length_min=config["arm_length_min"],
        arm_Tm_min=config["arm_Tm_min"],
        arm_Tm_max=config["arm_Tm_max"],
        Tm_parameters_probe=Tm_parameters_probe,
        Tm_chem_correction_param_probe=config["Tm_chem_correction_param_probe"],
    )

    ##### compute all required attributes #####
    logging.info("Computing Oligo Attributes")
    probe_database = pipeline.compute_probe_attributes(
        oligo_database=probe_database,
        Tm_parameters_probe=Tm_parameters_probe,
        Tm_chem_correction_param_probe=config["Tm_chem_correction_param_probe"],
    )

    ##### create probe sets #####
    probe_database, file_database, dir_probesets = pipeline.create_probe_sets(
        oligo_database=probe_database,
        probe_isoform_weight=config["probe_isoform_weight"],
        probe_Tm_weight=config["probe_Tm_weight"],
        probe_Tm_min=config["probe_Tm_min"],
        probe_Tm_opt=config["probe_Tm_opt"],
        probe_Tm_max=config["probe_Tm_max"],
        Tm_parameters_probe=Tm_parameters_probe,
        Tm_chem_correction_param_probe=config["Tm_chem_correction_param_probe"],
        probe_GC_weight=config["probe_GC_weight"],
        probe_GC_content_min=config["probe_GC_content_min"],
        probe_GC_content_opt=config["probe_GC_content_opt"],
        probe_GC_content_max=config["probe_GC_content_max"],
        probeset_size_opt=config["probeset_size_opt"],
        probeset_size_min=config["probeset_size_min"],
        max_graph_size=config["max_graph_size"],
        n_sets=config["n_sets"],
        distance_between_probes=config["distance_between_probes"],
    )

    logging.info(f"Oligo sets were saved in {dir_probesets}")
    logging.info("##### End of the pipeline. #####")


if __name__ == "__main__":
    main()
