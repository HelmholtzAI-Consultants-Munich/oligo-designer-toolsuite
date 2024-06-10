############################################
# imports
############################################
import logging
import os
import shutil
import warnings
from datetime import datetime
from pathlib import Path
from typing import List

import yaml
from Bio.SeqUtils import MeltingTemp as mt

from oligo_designer_toolsuite.database import (
    OligoAttributes,
    OligoDatabase,
    ReferenceDatabase,
)
from oligo_designer_toolsuite.oligo_efficiency_filter import (
    LowestSetScoring,
    WeightedIsoformTmGCOligoScoring,
)
from oligo_designer_toolsuite.oligo_property_filter import (
    GCClampFilter,
    GCContentFilter,
    HardMaskedSequenceFilter,
    HomopolymericRunsFilter,
    MeltingTemperatureNNFilter,
    PropertyFilter,
    SecondaryStructureFilter,
    SoftMaskedSequenceFilter,
)
from oligo_designer_toolsuite.oligo_selection import (
    OligosetGeneratorIndependentSet,
    heuristic_selection_independent_set,
)
from oligo_designer_toolsuite.oligo_specificity_filter import (
    BlastNFilter,
    CrossHybridizationFilter,
    ExactMatchFilter,
    RemoveByDegreePolicy,
    RemoveByLargerRegionPolicy,
    SpecificityFilter,
)
from oligo_designer_toolsuite.pipelines._utils import (
    base_parser,
    filtering_step,
    generation_step,
)
from oligo_designer_toolsuite.sequence_generator import OligoSequenceGenerator

############################################
# MERFISH Probe Designer Functions
############################################


class MerfishProbeDesigner:
    def __init__(
        self,
        file_regions: list,
        write_intermediate_steps: bool,
        dir_output: str,
        n_jobs: int,
    ) -> None:
        """Constructor for the MerfishProbeDesigner class."""
        ##### read the genes file #####
        if file_regions is None:
            warnings.warn(
                "No gene list file was provided! All genes from fasta file are used to generate the oligos. This chioce can use a lot of resources."
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

        self.subdir_db_oligos = "db_probes"
        self.subdir_db_reference = "db_reference"

        ##### setup logger #####
        timestamp = datetime.now()
        file_logger = os.path.join(
            self.dir_output,
            f"log_merfish_probe_designer_{timestamp.year}-{timestamp.month}-{timestamp.day}-{timestamp.hour}-{timestamp.minute}.txt",
        )
        logging.getLogger("log_name")
        logging.basicConfig(
            format="%(asctime)s [%(levelname)s] %(message)s",
            level=logging.NOTSET,
            handlers=[logging.FileHandler(file_logger), logging.StreamHandler()],
        )
        logging.captureWarnings(True)

        self.n_jobs = n_jobs
        self.oligo_attributes_calculator = OligoAttributes()

    @generation_step(step_name="Create Database")
    def create_oligo_database(
        self,
        oligo_length: int,
        files_fasta_oligo_database: list[str],
        min_oligos_per_region: int,
    ):
        """
        Creates an oligo database using sequences generated through a sliding window approach, loading them from specified FASTA files.
        It logs the creation parameters, generates oligo sequences, constructs the database, logs the database stats, and optionally saves the database.

        :param oligo_length: Length of oligos to generate.
        :type oligo_length: int
        :param oligo_length_max: Maximum length of oligos to generate.
        :type oligo_length_max: int
        :param files_fasta_oligo_database: List of paths to FASTA files to read oligo sequences from.
        :type files_fasta_oligo_database: list[str]
        :param min_oligos_per_region: Minimum number of oligos required per gene region for inclusion in the database.
        :type min_oligos_per_region: int
        :return: A tuple containing the oligo database object and path to the saved database file.
        :rtype: (OligoDatabase, str)
        """

        ##### creating the oligo sequences #####
        oligo_sequences = OligoSequenceGenerator(dir_output=self.dir_output)
        oligo_fasta_file = oligo_sequences.create_sequences_sliding_window(
            filename_out="oligo_sequences",
            files_fasta_in=files_fasta_oligo_database,
            length_interval_sequences=(oligo_length, oligo_length),
            region_ids=self.gene_ids,
            n_jobs=self.n_jobs,
        )

        ##### creating the oligo database #####
        oligo_database = OligoDatabase(
            min_oligos_per_region=min_oligos_per_region,
            write_regions_with_insufficient_oligos=True,
            database_name=self.subdir_db_oligos,
            dir_output=self.dir_output,
        )
        oligo_database.load_sequences_from_fasta(
            files_fasta=[oligo_fasta_file],
            sequence_type="target",
            region_ids=self.gene_ids,
        )

        ##### save database #####
        if self.write_intermediate_steps:
            file_database = oligo_database.save_database(filename="1_db_initial")
        else:
            file_database = ""

        return oligo_database, file_database

    @filtering_step(step_name="Property Filters")
    def filter_by_property(
        self,
        oligo_database: OligoDatabase,
        GC_content_min: int,
        GC_content_max: int,
        Tm_min: int,
        Tm_max: int,
        Tm_parameters: dict,
        secondary_structures_T: float,
        secondary_structures_threshold_deltaG: float,
        homopolymeric_base_n: int,
    ):
        """
        Applies multiple sequence property filters to an oligonucleotide database to refine the selection of oligos based on specified criteria.
        The function logs parameters, applies the filters, logs the reduction in oligos, and optionally saves the filtered database.

        :param oligo_database: Oligo database object containing sequences to filter.
        :type oligo_database: OligoDatabase
        :param GC_content_min: Minimum GC content for oligos.
        :type GC_content_min: int
        :param GC_content_max: Maximum GC content for oligos.
        :type GC_content_max: int
        :param Tm_min: Minimum melting temperature for oligos.
        :type Tm_min: int
        :param Tm_max: Maximum melting temperature for oligos.
        :type Tm_max: int
        :param Tm_parameters: Parameters for calculating melting temperatures.
        :type Tm_parameters: dict
        :param secondary_structures_T: Temperature for secondary structure calculations.
        :type secondary_structures_T: float
        :param secondary_structures_threshold_deltaG: Threshold for secondary structure calculations.
        :type secondary_structures_threshold_deltaG: float
        :param homopolymeric_base_n: Number of repeated bases to filter.
        :type homopolymeric_base_n: int
        :return: A tuple containing the filtered oligo database object and path to the saved database file.
        :rtype: (OligoDatabase, str)


        """
        # define the filters
        hard_masked_sequences = HardMaskedSequenceFilter()
        soft_masked_sequences = SoftMaskedSequenceFilter()
        gc_content = GCContentFilter(GC_content_min=GC_content_min, GC_content_max=GC_content_max)
        melting_temperature = MeltingTemperatureNNFilter(
            Tm_min=Tm_min,
            Tm_max=Tm_max,
            Tm_parameters=Tm_parameters,
        )
        secondary_sctructure = SecondaryStructureFilter(
            T=secondary_structures_T,
            thr_DG=secondary_structures_threshold_deltaG,
        )
        homopolymeric_runs = HomopolymericRunsFilter(
            base_n=homopolymeric_base_n,
        )

        filters = [
            hard_masked_sequences,
            soft_masked_sequences,
            gc_content,
            melting_temperature,
            secondary_sctructure,
            homopolymeric_runs,
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
            file_database = oligo_database.save_database(filename="2_db_property_filter")
        else:
            file_database = ""

        return oligo_database, file_database

    @filtering_step(step_name="Specificty Filters")
    def filter_by_specificity(
        self,
        oligo_database: OligoDatabase,
        files_fasta_reference_database: List[str],
        blastn_search_parameters: dict,
        blastn_hit_parameters: dict,
        cross_hybridization_blastn_search_parameters: dict,
        cross_hybridization_blastn_hit_parameters: dict,
    ):
        """
        Applies a series of filters to an oligonucleotide database to ensure the specificity of the oligos
        based on cross-hybridization and hybridization probability criteria. The function initializes specific alignment
        methods for each filtering criterion and logs the process and results.

        :param oligo_database: Oligo database object containing sequences to filter.
        :type oligo_database: OligoDatabase
        :param files_fasta_reference_database: List of paths to FASTA files to read reference sequences from.
        :type files_fasta_reference_database: List[str]
        :param blastn_search_parameters: Parameters for the BLASTN search.
        :type blastn_search_parameters: dict
        :param blastn_hit_parameters: Parameters for the BLASTN hit.
        :type blastn_hit_parameters: dict
        :param cross_hybridization_blastn_search_parameters: Parameters for the BLASTN search for cross-hybridization.
        :type cross_hybridization_blastn_search_parameters: dict
        :param cross_hybridization_blastn_hit_parameters: Parameters for the BLASTN hit for cross-hybridization.
        :type cross_hybridization_blastn_hit_parameters: dict
        :return: A tuple containing the filtered oligo database object and path to the saved database file.
        :rtype: (OligoDatabase, str)
        """

        reference_database = ReferenceDatabase(
            database_name=self.subdir_db_reference, dir_output=self.dir_output
        )
        reference_database.load_sequences_from_fasta(
            files_fasta=files_fasta_reference_database, database_overwrite=False
        )

        ##### specificity filters #####
        # removing duplicated oligos from the region with the most oligos
        exact_matches = ExactMatchFilter(policy=RemoveByLargerRegionPolicy())
        # BlastN filter
        blastn_filter = BlastNFilter(
            search_parameters=blastn_search_parameters,
            hit_parameters=blastn_hit_parameters,
        )

        # Cross-hybridization filter with BlastN
        cross_hybridization = CrossHybridizationFilter(
            policy=RemoveByLargerRegionPolicy(),
            alignment_method=blastn_filter,
            database_name_reference=self.subdir_db_reference,
            filter_name="cross_hybridization",
            dir_output=self.dir_output,
        )

        filters = [exact_matches, blastn_filter, cross_hybridization]
        specificity_filter = SpecificityFilter(filters=filters)
        oligo_database = specificity_filter.apply(
            sequence_type="oligo",
            oligo_database=oligo_database,
            reference_database=reference_database,
            n_jobs=self.n_jobs,
        )

        # write the intermediate result in a file
        if self.write_intermediate_steps:
            file_database = oligo_database.save_database(filename="3_db_specificty_filter")
        else:
            file_database = ""

        for directory in [
            reference_database.dir_output,
            blastn_filter.dir_output,
            cross_hybridization.dir_output,
        ]:
            if os.path.exists(directory):
                shutil.rmtree(directory)

        return oligo_database, file_database

    def compute_oligo_attributes(
        self,
        oligo_database: OligoDatabase,
        Tm_parameters_probe: dict,
    ):
        oligo_database = self.oligo_attributes_calculator.calculate_oligo_length(
            oligo_database=oligo_database
        )
        oligo_database = self.oligo_attributes_calculator.calculate_GC_content(
            oligo_database=oligo_database, sequence_type="oligo"
        )
        oligo_database = self.oligo_attributes_calculator.calculate_TmNN(
            oligo_database=oligo_database,
            sequence_type="oligo",
            Tm_parameters=Tm_parameters_probe,
        )
        oligo_database = self.oligo_attributes_calculator.calculate_num_targeted_transcripts(
            oligo_database=oligo_database
        )
        oligo_database = self.oligo_attributes_calculator.calculate_isoform_consensus(
            oligo_database=oligo_database
        )

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
            isoform_weight=probe_isoform_weight,
            Tm_weight=probe_Tm_weight,
            GC_weight=probe_GC_weight,
        )
        set_scoring = LowestSetScoring(ascending=True)
        probeset_generator = OligosetGeneratorIndependentSet(
            opt_oligoset_size=probeset_size_opt,
            min_oligoset_size=probeset_size_min,
            oligos_scoring=probes_scoring,
            set_scoring=set_scoring,
            heuristic_selection=heuristic_selection_independent_set,
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
            file_probesets = oligo_database.write_oligosets_to_table()
        else:
            file_database = ""
            file_probesets = ""

        return oligo_database, file_database, file_probesets


############################################
# MERFISH Readout Probe Designer Functions
############################################


class MerfishReadoutProbeDesigner:
    def __init__(
        self,
        write_intermediate_steps: bool,
        dir_output: str,
        n_jobs: int,
    ) -> None:
        """Constructor for the MerfishReadoutProbeDesigner class."""

        self.write_intermediate_steps = write_intermediate_steps

        ##### create the output folder #####
        self.dir_output = os.path.abspath(dir_output)
        Path(self.dir_output).mkdir(parents=True, exist_ok=True)

        self.subdir_db_oligos = "db_readout_probes"
        self.subdir_db_reference = "db_readout_reference"

        ##### setup logger #####
        timestamp = datetime.now()
        file_logger = os.path.join(
            self.dir_output,
            f"log_merfish_readout_probe_designer_{timestamp.year}-{timestamp.month}-{timestamp.day}-{timestamp.hour}-{timestamp.minute}.txt",
        )
        logging.getLogger("log_name")
        logging.basicConfig(
            format="%(asctime)s [%(levelname)s] %(message)s",
            level=logging.NOTSET,
            handlers=[logging.FileHandler(file_logger), logging.StreamHandler()],
        )
        logging.captureWarnings(True)

        self.n_jobs = n_jobs

    @generation_step(step_name="Create Readout Probe Database")
    def create_oligo_database(
        self,
        oligo_length: int,
        readout_sequence_probs: dict,
        initial_num_sequences: int = 100000,
    ):
        """
        Creates an oligo database using sequences generated through a sliding window approach, loading them from specified FASTA files.
        It logs the creation parameters, generates oligo sequences, constructs the database, logs the database stats, and optionally saves the database.

        :param oligo_length: Length of oligos to generate.
        :type oligo_length: int
        :param oligo_length_max: Maximum length of oligos to generate.
        :type oligo_length_max: int
        :param readout_sequence_probs: Dictionary containing the base alphabet and probability for generating readout probe sequences.
        :type readout_sequence_probs: dict
        :param initial_num_sequences: Number of initial sequences to generate. Default is 100,000.
        :type initial_num_sequences: int
        :return: A tuple containing the oligo database object and path to the saved database file.
        :rtype: (OligoDatabase, str)
        """

        ##### creating the oligo sequences #####
        oligo_sequences = OligoSequenceGenerator(dir_output=self.dir_output)
        oligo_fasta_file = oligo_sequences.create_sequences_random(
            filename_out="readout_probes_sequences",
            length_sequences=oligo_length,
            num_sequences=initial_num_sequences,
            name_sequences="readout_probes",
            base_alphabet_with_probability=readout_sequence_probs,
        )

        ##### creating the oligo database #####
        oligo_database = OligoDatabase(
            write_regions_with_insufficient_oligos=True,
            database_name=self.subdir_db_oligos,
            dir_output=self.dir_output,
        )
        oligo_database.load_sequences_from_fasta(
            files_fasta=[oligo_fasta_file],
            sequence_type="target",
        )

        ##### save database #####
        if self.write_intermediate_steps:
            file_database = oligo_database.save_database(filename="1_db_initial")
        else:
            file_database = ""

        return oligo_database, file_database

    @filtering_step(step_name="Property Filters")
    def filter_by_property(
        self,
        oligo_database: OligoDatabase,
        GC_content_min: int,
        GC_content_max: int,
        homopolymeric_base_n: int,
    ):
        """
        Applies multiple sequence property filters to an oligonucleotide database to refine the selection of oligos based on specified criteria.
        The function logs parameters, applies the filters, logs the reduction in oligos, and optionally saves the filtered database.

        :param oligo_database: Oligo database object containing sequences to filter.
        :type oligo_database: OligoDatabase
        :param GC_content_min: Minimum GC content for oligos.
        :type GC_content_min: int
        :param GC_content_max: Maximum GC content for oligos.
        :type GC_content_max: int
        :param homopolymeric_base_n: Number of repeated bases to filter.
        :type homopolymeric_base_n: int
        :return: A tuple containing the filtered oligo database object and path to the saved database file.
        :rtype: (OligoDatabase, str)


        """
        # define the filters
        gc_content = GCContentFilter(GC_content_min=GC_content_min, GC_content_max=GC_content_max)
        homopolymeric_runs = HomopolymericRunsFilter(
            base_n=homopolymeric_base_n,
        )

        filters = [
            gc_content,
            homopolymeric_runs,
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
            file_database = oligo_database.save_database(filename="2_db_property_filter")
        else:
            file_database = ""

        return oligo_database, file_database

    @filtering_step(step_name="Specificty Filters")
    def filter_by_specificity(
        self,
        oligo_database: OligoDatabase,
        files_fasta_reference_database: List[str],
        blastn_search_parameters: dict,
        blastn_hit_parameters: dict,
        cross_hybridization_blastn_search_parameters: dict,
        cross_hybridization_blastn_hit_parameters: dict,
    ):
        """
        Applies a series of filters to an oligonucleotide database to ensure the specificity of the oligos
        based on cross-hybridization and hybridization probability criteria. The function initializes specific alignment
        methods for each filtering criterion and logs the process and results.

        :param oligo_database: Oligo database object containing sequences to filter.
        :type oligo_database: OligoDatabase
        :param files_fasta_reference_database: List of paths to FASTA files to read reference sequences from.
        :type files_fasta_reference_database: List[str]
        :param blastn_search_parameters: Parameters for the BLASTN search.
        :type blastn_search_parameters: dict
        :param blastn_hit_parameters: Parameters for the BLASTN hit.
        :type blastn_hit_parameters: dict
        :param cross_hybridization_blastn_search_parameters: Parameters for the BLASTN search for cross-hybridization.
        :type cross_hybridization_blastn_search_parameters: dict
        :param cross_hybridization_blastn_hit_parameters: Parameters for the BLASTN hit for cross-hybridization.
        :type cross_hybridization_blastn_hit_parameters: dict
        :return: A tuple containing the filtered oligo database object and path to the saved database file.
        :rtype: (OligoDatabase, str)
        """

        reference_database = ReferenceDatabase(
            database_name=self.subdir_db_reference, dir_output=self.dir_output
        )
        reference_database.load_sequences_from_fasta(
            files_fasta=files_fasta_reference_database, database_overwrite=False
        )

        ##### specificity filters #####
        # removing duplicated oligos from the region with the most oligos
        exact_matches = ExactMatchFilter(policy=RemoveByLargerRegionPolicy())
        # BlastN filter
        blastn_filter = BlastNFilter(
            search_parameters=blastn_search_parameters,
            hit_parameters=blastn_hit_parameters,
        )

        # Cross-hybridization filter with BlastN
        cross_hybridization = CrossHybridizationFilter(
            policy=RemoveByDegreePolicy(),
            alignment_method=blastn_filter,
            database_name_reference=self.subdir_db_reference,
            filter_name="cross_hybridization",
            dir_output=self.dir_output,
        )

        filters = [exact_matches, blastn_filter, cross_hybridization]
        specificity_filter = SpecificityFilter(filters=filters)
        oligo_database = specificity_filter.apply(
            sequence_type="oligo",
            oligo_database=oligo_database,
            reference_database=reference_database,
            n_jobs=self.n_jobs,
        )

        # write the intermediate result in a file
        if self.write_intermediate_steps:
            file_database = oligo_database.save_database(filename="3_db_specificty_filter")
        else:
            file_database = ""

        for directory in [
            reference_database.dir_output,
            blastn_filter.dir_output,
            cross_hybridization.dir_output,
        ]:
            if os.path.exists(directory):
                shutil.rmtree(directory)

        return oligo_database, file_database


############################################
# MERFISH Primer Designer Functions
############################################


class MerfishPrimerDesigner:
    @filtering_step(step_name="Property Filters")
    def filter_by_property(
        self,
        oligo_database: OligoDatabase,
        GC_content_min: int,
        GC_content_max: int,
        GC_clamp_n_bases: int,
        GC_clamp_n_GC: int,
        Tm_min: int,
        Tm_max: int,
        Tm_parameters: dict,
        secondary_structures_T: float,
        secondary_structures_threshold_deltaG: float,
    ):

        gc_clamp = GCClampFilter(n_bases=GC_clamp_n_bases, n_GC=GC_clamp_n_GC)  # 1 and 2


############################################
# MERFISH Probe Designer Pipeline
############################################


def main():
    """Main function to execute the oligo probe design pipeline based on user configurations.
    It parses command line arguments, initializes the pipeline with specified parameters, creates an oligo database,
    applies various filtering criteria, computes necessary oligo attributes, and finally generates oligo sets.
    """

    args = base_parser()

    ##### read the config file #####
    with open(args["config"], "r") as handle:
        config = yaml.safe_load(handle)

    ##### initialize probe designer pipeline #####
    pipeline = MerfishProbeDesigner(
        file_regions=config["file_regions"],
        write_intermediate_steps=config["write_intermediate_steps"],
        dir_output=config["dir_output"],
        n_jobs=config["n_jobs"],
    )

    ##### create oligo database #####
    oligo_database, file_database = pipeline.create_oligo_database(
        oligo_length=config["probe_length"],
        files_fasta_oligo_database=config["files_fasta_probe_database"],
        # we should have at least "min_oligoset_size" oligos per gene to create one set
        min_oligos_per_region=config["probeset_size"],
    )

    ##### preprocess melting temperature params #####
    probe_Tm_parameters = config["probe_Tm_parameters"]
    probe_Tm_parameters["nn_table"] = getattr(mt, probe_Tm_parameters["nn_table"])
    probe_Tm_parameters["tmm_table"] = getattr(mt, probe_Tm_parameters["tmm_table"])
    probe_Tm_parameters["imm_table"] = getattr(mt, probe_Tm_parameters["imm_table"])
    probe_Tm_parameters["de_table"] = getattr(mt, probe_Tm_parameters["de_table"])

    ##### filter oligos by property #####
    oligo_database, file_database = pipeline.filter_by_property(
        oligo_database,
        GC_content_min=config["probe_GC_content_min"],
        GC_content_max=config["probe_GC_content_max"],
        Tm_min=config["probe_Tm_min"],
        Tm_max=config["probe_Tm_max"],
        Tm_parameters=probe_Tm_parameters,
        secondary_structures_T=config["probe_T_secondary_structure"],
        secondary_structures_threshold_deltaG=config["probe_thr_deltaG_secondary_structure"],
        homopolymeric_base_n=config["probe_homopolymeric_base_n"],
    )

    # ##### filter oligos by specificity #####
    oligo_database, file_database = pipeline.filter_by_specificity(
        oligo_database=oligo_database,
        files_fasta_reference_database=config["probe_files_fasta_reference_database"],
        blastn_search_parameters=config["specificity_blastn_search_parameters"],
        blastn_hit_parameters=config["specificity_blastn_hit_parameters"],
        cross_hybridization_blastn_search_parameters=config["cross_hybridization_blastn_search_parameters"],
        cross_hybridization_blastn_hit_parameters=config["cross_hybridization_blastn_hit_parameters"],
    )

    ##### compute all required attributes #####
    logging.info("Computing Oligo Attributes")
    oligo_database = pipeline.compute_oligo_attributes(
        oligo_database=oligo_database,
        Tm_parameters_probe=probe_Tm_parameters,
    )

    ##### create probe sets #####
    oligo_database, file_database, dir_oligosets = pipeline.create_probe_sets(
        oligo_database=oligo_database,
        probe_isoform_weight=config["probe_isoform_weight"],
        probe_Tm_weight=config["probe_Tm_weight"],
        probe_Tm_min=config["probe_Tm_min"],
        probe_Tm_opt=config["probe_Tm_opt"],
        probe_Tm_max=config["probe_Tm_max"],
        Tm_parameters_probe=probe_Tm_parameters,
        probe_GC_weight=config["probe_GC_weight"],
        probe_GC_content_min=config["probe_GC_content_min"],
        probe_GC_content_opt=config["probe_GC_content_opt"],
        probe_GC_content_max=config["probe_GC_content_max"],
        probeset_size_opt=config["probeset_size"],
        probeset_size_min=config["probeset_size"],
        max_graph_size=config["max_graph_size"],
        n_sets=config["n_sets"],
        distance_between_probes=config["distance_between_probes"],
    )

    logging.info(f"Oligo sets were saved in {dir_oligosets}")
    logging.info("##### End of the pipeline. #####")


if __name__ == "__main__":
    main()
