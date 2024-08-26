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
    AverageSetScoring,
    WeightedTmGCOligoScoring,
)
from oligo_designer_toolsuite.oligo_property_filter import (
    GCContentFilter,
    HardMaskedSequenceFilter,
    HomopolymericRunsFilter,
    MeltingTemperatureNNFilter,
    PropertyFilter,
    SecondaryStructureFilter,
    SelfComplementFilter,
    SoftMaskedSequenceFilter,
)
from oligo_designer_toolsuite.oligo_selection import (
    GraphBasedSelectionPolicy,
    OligosetGeneratorIndependentSet,
)
from oligo_designer_toolsuite.oligo_specificity_filter import (
    BlastNFilter,
    BowtieFilter,
    CrossHybridizationFilter,
    ExactMatchFilter,
    HybridizationProbabilityFilter,
    RemoveByLargerRegionPolicy,
    SpecificityFilter,
)
from oligo_designer_toolsuite.pipelines._utils import base_parser, pipeline_step_basic
from oligo_designer_toolsuite.sequence_generator import OligoSequenceGenerator

############################################
# Oligo-seq Designer Functions
############################################


class OligoSeqProbeDesigner:
    """This class handles the setup and execution of the oligo-seq probe design pipeline, including reading gene region files,
    managing output directories, and configuring logging. It supports parallel processing to optimize the probe design tasks
    and allows intermediate steps to be saved for review. The class is designed to work with various configuration settings to
    customize the probe generation process according to user-defined parameters.

    :param file_regions: List of file paths containing regions of interest. If None, all regions from the fasta file are used.
    :type file_regions: list
    :param write_intermediate_steps: Flag indicating whether intermediate results should be saved.
    :type write_intermediate_steps: bool
    :param dir_output: Directory path where output and logs will be saved.
    :type dir_output: str
    :param n_jobs: Number of parallel jobs to run.
    :type n_jobs: int
    """

    def __init__(
        self,
        file_regions: list,
        write_intermediate_steps: bool,
        dir_output: str,
        n_jobs: int,
    ) -> None:
        """Constructor for the OligoSeqProbeDesigner class."""

        self.write_intermediate_steps = write_intermediate_steps
        self.n_jobs = n_jobs
        self.oligo_attributes = OligoAttributes()

        ##### create the output folder #####
        self.dir_output = os.path.abspath(dir_output)
        Path(self.dir_output).mkdir(parents=True, exist_ok=True)

        self.subdir_db_oligos = "db_probes"
        self.subdir_db_reference = "db_reference"

        ##### setup logger #####
        timestamp = datetime.now()
        file_logger = os.path.join(
            self.dir_output,
            f"log_oligo_seq_probe_designer_{timestamp.year}-{timestamp.month}-{timestamp.day}-{timestamp.hour}-{timestamp.minute}.txt",
        )
        logging.getLogger("log_name")
        logging.basicConfig(
            format="%(asctime)s [%(levelname)s] %(message)s",
            level=logging.NOTSET,
            handlers=[logging.FileHandler(file_logger)],
        )
        logging.captureWarnings(True)

    @pipeline_step_basic(step_name="Create Database")
    def create_oligo_database(
        self,
        gene_ids: list,
        oligo_length_min: int,
        oligo_length_max: int,
        files_fasta_oligo_database: list[str],
        min_oligos_per_region: int,
    ):
        """
        Creates an oligo database using sequences generated through a sliding window approach, loading them from specified FASTA files.
        It logs the creation parameters, generates oligo sequences, constructs the database, logs the database stats, and optionally saves the database.

        :param oligo_length_min: Minimum length of oligos to generate.
        :type oligo_length_min: int
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
            files_fasta_in=files_fasta_oligo_database,
            length_interval_sequences=(oligo_length_min, oligo_length_max),
            region_ids=gene_ids,
            n_jobs=self.n_jobs,
        )

        ##### creating the oligo database #####
        oligo_database = OligoDatabase(
            min_oligos_per_region=min_oligos_per_region,
            write_regions_with_insufficient_oligos=True,
            lru_db_max_in_memory=self.n_jobs * 2 + 2,
            database_name=self.subdir_db_oligos,
            dir_output=self.dir_output,
            n_jobs=1,
        )
        oligo_database.load_database_from_fasta(
            files_fasta=oligo_fasta_file,
            sequence_type="target",
            region_ids=gene_ids,
        )

        ##### save database #####
        if self.write_intermediate_steps:
            file_database = oligo_database.save_database(dir_database="1_db_initial")
        else:
            file_database = ""

        dir = oligo_sequences.dir_output
        shutil.rmtree(dir) if os.path.exists(dir) else None

        return oligo_database, file_database

    @pipeline_step_basic(step_name="Property Filters")
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
        max_len_selfcomp: int,
        Tm_parameters: dict,
        Tm_chem_correction_parameters: dict,
    ):
        """
        Applies multiple sequence property filters to an oligonucleotide database to refine the selection of oligos based on specified criteria.
        The function logs parameters, applies the filters, logs the reduction in oligos, and optionally saves the filtered database.

        :param oligo_database: The database of oligos to be filtered.
        :type oligo_database: OligoDatabase
        :param GC_content_min: Minimum GC content percentage to allow in an oligo.
        :type GC_content_min: int
        :param GC_content_max: Maximum GC content percentage to allow in an oligo.
        :type GC_content_max: int
        :param Tm_min: Minimum melting temperature (Tm) to allow in an oligo.
        :type Tm_min: int
        :param Tm_max: Maximum melting temperature (Tm) to allow in an oligo.
        :type Tm_max: int
        :param secondary_structures_T: Temperature at which to evaluate secondary structures.
        :type secondary_structures_T: float
        :param secondary_structures_threshold_deltaG: Threshold free energy below which secondary structures are considered problematic.
        :type secondary_structures_threshold_deltaG: float
        :param homopolymeric_base_n: Bases to check for homopolymeric runs.
        :type homopolymeric_base_n: str
        :param max_len_selfcomp: Maximum allowable length of self-complementary sequences.
        :type max_len_selfcomp: int
        :param Tm_parameters: Parameters for melting temperature calculation.
        :type Tm_parameters: dict
        :param Tm_chem_correction_parameters: Parameters for chemical correction of melting temperature.
        :type Tm_chem_correction_parameters: dict
        :return: Tuple containing the filtered oligo database and path to the saved database file.
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
            Tm_chem_correction_parameters=Tm_chem_correction_parameters,
        )
        secondary_sctructure = SecondaryStructureFilter(
            T=secondary_structures_T,
            thr_DG=secondary_structures_threshold_deltaG,
        )
        homopolymeric_runs = HomopolymericRunsFilter(
            base_n=homopolymeric_base_n,
        )
        self_comp = SelfComplementFilter(max_len_selfcomplement=max_len_selfcomp)

        filters = [
            hard_masked_sequences,
            soft_masked_sequences,
            homopolymeric_runs,
            gc_content,
            self_comp,
            melting_temperature,
            secondary_sctructure,
            homopolymeric_runs,
            self_comp,
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
            file_database = oligo_database.save_database(dir_database="2_db_property_filter")
        else:
            file_database = ""

        return oligo_database, file_database

    @pipeline_step_basic(step_name="Specificity Filters")
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
    ):
        """
        Applies a series of filters to an oligonucleotide database to ensure the specificity of the oligos
        based on cross-hybridization and hybridization probability criteria. The function initializes specific alignment
        methods for each filtering criterion and logs the process and results.

        :param oligo_database: Database to filter.
        :type oligo_database: OligoDatabase
        :param files_fasta_reference_database: List of reference database FASTA files.
        :type files_fasta_reference_database: List[str]
        :param cross_hybridization_alignment_method: Alignment method for cross-hybridization.
        :type cross_hybridization_alignment_method: str
        :param cross_hybridization_search_parameters: Search parameters for cross-hybridization alignment.
        :type cross_hybridization_search_parameters: dict
        :param cross_hybridization_hit_parameters: Hit parameters for cross-hybridization alignment.
        :type cross_hybridization_hit_parameters: dict
        :param hybridization_probability_alignment_method: Alignment method for calculating hybridization probability.
        :type hybridization_probability_alignment_method: str
        :param hybridization_probability_search_parameters: Search parameters for hybridization probability.
        :type hybridization_probability_search_parameters: dict
        :param hybridization_probability_hit_parameters: Hit parameters for hybridization probability.
        :type hybridization_probability_hit_parameters: dict
        :param hybridization_probability_threshold: Threshold for deciding significant hybridization probability.
        :type hybridization_probability_threshold: float
        :return: Updated oligo database and path to the saved database file.
        :rtype: (OligoDatabase, str)
        """

        def _get_alignment_method(
            alignment_method,
            search_parameters,
            hit_parameters,
            filter_name_specification: str = "",
        ):
            if alignment_method == "blastn":
                return BlastNFilter(
                    search_parameters=search_parameters,
                    hit_parameters=hit_parameters,
                    dir_output=self.dir_output,
                    filter_name="blast_" + filter_name_specification,
                )
            elif alignment_method == "bowtie":
                return BowtieFilter(
                    search_parameters=search_parameters,
                    hit_parameters=hit_parameters,
                    dir_output=self.dir_output,
                    filter_name="bowtie_" + filter_name_specification,
                )
            else:
                raise ValueError(f"The alignment method {alignment_method} is not supported.")

        ##### define reference database #####
        reference_database = ReferenceDatabase(
            database_name=self.subdir_db_reference, dir_output=self.dir_output
        )
        reference_database.load_database_from_fasta(
            files_fasta=files_fasta_reference_database, database_overwrite=False
        )

        ##### specificity filters #####
        # removing duplicated oligos from the region with the most oligos
        exact_matches = ExactMatchFilter(policy=RemoveByLargerRegionPolicy(), filter_name="exact_match")

        cross_hybridization_aligner = _get_alignment_method(
            alignment_method=cross_hybridization_alignment_method,
            search_parameters=cross_hybridization_search_parameters,
            hit_parameters=cross_hybridization_hit_parameters,
            filter_name_specification="cross_hybridization",
        )
        cross_hybridization = CrossHybridizationFilter(
            policy=RemoveByLargerRegionPolicy(),
            alignment_method=cross_hybridization_aligner,
            database_name_reference=self.subdir_db_reference,
            dir_output=self.dir_output,
        )

        hybridization_probability_aligner = _get_alignment_method(
            alignment_method=hybridization_probability_alignment_method,
            search_parameters=hybridization_probability_search_parameters,
            hit_parameters=hybridization_probability_hit_parameters,
            filter_name_specification="hybridization_probability",
        )
        hybridization_probability = HybridizationProbabilityFilter(
            alignment_method=hybridization_probability_aligner,
            threshold=hybridization_probability_threshold,
            dir_output=self.dir_output,
        )

        filters = [exact_matches, hybridization_probability, cross_hybridization]
        specificity_filter = SpecificityFilter(filters=filters)
        oligo_database = specificity_filter.apply(
            sequence_type="oligo",
            oligo_database=oligo_database,
            reference_database=reference_database,
            n_jobs=self.n_jobs,
        )

        # write the intermediate result in a file
        if self.write_intermediate_steps:
            file_database = oligo_database.save_database(dir_database="3_db_specificity_filter")
        else:
            file_database = ""

        # remove all directories of intermediate steps
        for directory in [
            reference_database.dir_output,
            cross_hybridization_aligner.dir_output,
            cross_hybridization.dir_output,
            hybridization_probability_aligner.dir_output,
            hybridization_probability.dir_output,
        ]:
            if os.path.exists(directory):
                shutil.rmtree(directory)

        return oligo_database, file_database

    @pipeline_step_basic(step_name="Oligo Selection")
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
        opt_oligoset_size: int,
        min_oligoset_size: int,
        max_oligos: int,
        n_sets: int,
        distance_between_oligos: int,
    ):
        """
        Initializes and applies the oligo set generation process with specified melting temperature and GC content parameters.
        This function configures scoring mechanisms for oligos and sets, generates oligo sets, and logs the process.
        It also saves intermediate results if specified.

        :param oligo_database: Database containing oligo sequences and related data.
        :type oligo_database: OligoDatabase
        :param Tm_min: Minimum melting temperature for oligos.
        :type Tm_min: float
        :param Tm_opt: Optimal melting temperature for oligos.
        :type Tm_opt: float
        :param Tm_max: Maximum melting temperature for oligos.
        :type Tm_max: float
        :param Tm_parameters: Parameters for Tm calculations.
        :type Tm_parameters: dict
        :param Tm_chem_correction_parameters: Chemical correction parameters for Tm calculations.
        :type Tm_chem_correction_parameters: dict
        :param GC_content_min: Minimum acceptable GC content percentage.
        :type GC_content_min: float
        :param GC_content_opt: Optimal GC content percentage.
        :type GC_content_opt: float
        :param GC_content_max: Maximum acceptable GC content percentage.
        :type GC_content_max: float
        :param opt_oligoset_size: Optimal size of oligo sets to generate.
        :type opt_oligoset_size: int
        :param min_oligoset_size: Minimum acceptable size of oligo sets.
        :type min_oligoset_size: int
        :param max_oligos: Maximum number of oligos to consider in the generation process.
        :type max_oligos: int
        :param n_sets: Number of sets to generate.
        :type n_sets: int
        :param distance_between_oligos: Minimum distance required between oligos within a set.
        :type distance_between_oligos: int
        :return: The updated oligo database, file path of the database, file path of oligo sets.
        :rtype: tuple(OligoDatabase, str, str)
        """
        oligos_scoring = WeightedTmGCOligoScoring(
            Tm_min=Tm_min,
            Tm_opt=Tm_opt,
            Tm_max=Tm_max,
            GC_content_min=GC_content_min,
            GC_content_opt=GC_content_opt,
            GC_content_max=GC_content_max,
            Tm_parameters=Tm_parameters,
            Tm_chem_correction_parameters=Tm_chem_correction_parameters,
        )
        set_scoring = AverageSetScoring(ascending=True)

        selection_policy = GraphBasedSelectionPolicy(
            set_size_opt=opt_oligoset_size,
            set_size_min=min_oligoset_size,
            n_sets=n_sets,
            ascending=True,
            set_scoring=set_scoring,
        )
        oligoset_generator = OligosetGeneratorIndependentSet(
            opt_oligoset_size=opt_oligoset_size,
            min_oligoset_size=min_oligoset_size,
            oligos_scoring=oligos_scoring,
            set_scoring=set_scoring,
            heuristic_selection=selection_policy,
            max_oligos=max_oligos,
            distance_between_oligos=distance_between_oligos,
        )
        oligo_database = oligoset_generator.apply(
            oligo_database=oligo_database,
            sequence_type="oligo",
            pre_filter=False,
            n_jobs=self.n_jobs,
        )

        # write the intermediate result in a file
        if self.write_intermediate_steps:
            file_database = oligo_database.save_database(dir_database="4_db_oligosets")
            file_oligosets = oligo_database.write_oligosets_to_table()
        else:
            file_database = ""
            file_oligosets = ""

        return oligo_database, file_database, file_oligosets

    def compute_oligo_attributes(
        self,
        oligo_database: OligoDatabase,
        secondary_structures_T: float,
        Tm_parameters: dict,
        Tm_chem_correction_parameters: dict,
    ):
        """
        Computes various attributes for oligos and stores them in the provided oligo database.

        :param oligo_database: The database containing oligos for which attributes are to be computed.
        :type oligo_database: OligoDatabase
        :param secondary_structures_T: Temperature at which secondary structure delta G is calculated.
        :type secondary_structures_T: float
        :param Tm_parameters: Parameters to calculate melting temperature.
        :type Tm_parameters: dict
        :param Tm_chem_correction_parameters: Chemical correction parameters for melting temperature calculation.
        :type Tm_chem_correction_parameters: dict
        :return: Updated oligo database with new attributes.
        :rtype: OligoDatabase
        """
        oligo_database = self.oligo_attributes.calculate_oligo_length(oligo_database=oligo_database)
        oligo_database = self.oligo_attributes.calculate_GC_content(
            oligo_database=oligo_database, sequence_type="oligo"
        )
        oligo_database = self.oligo_attributes.calculate_TmNN(
            oligo_database=oligo_database,
            sequence_type="oligo",
            Tm_parameters=Tm_parameters,
            Tm_chem_correction_parameters=Tm_chem_correction_parameters,
        )
        oligo_database = self.oligo_attributes.calculate_num_targeted_transcripts(
            oligo_database=oligo_database
        )
        oligo_database = self.oligo_attributes.calculate_isoform_consensus(oligo_database=oligo_database)
        oligo_database = self.oligo_attributes.calculate_length_selfcomplement(
            oligo_database=oligo_database, sequence_type="oligo"
        )
        oligo_database = self.oligo_attributes.calculate_DG_secondary_structure(
            oligo_database=oligo_database, sequence_type="oligo", T=secondary_structures_T
        )

        return oligo_database

    def generate_output(self, oligo_database: OligoDatabase, top_n_sets: int):

        attributes = [
            "chromosome",
            "start",
            "end",
            "strand",
            "oligo",
            "target",
            "length",
            "GC_content",
            "TmNN",
            "num_targeted_transcripts",
            "isoform_consensus",
            "length_selfcomplement",
            "DG_secondary_structure" "source",
            "species",
            "annotation_release",
            "genome_assembly",
            "regiontype",
            "gene_id",
            "transcript_id",
            "exon_number",
        ]
        oligo_database.write_oligosets_to_yaml(
            attributes=attributes, top_n_sets=top_n_sets, ascending=True, filename="oligo_seq_probes.yml"
        )

    def compute_oligo_attributes(
        self,
        oligo_database: OligoDatabase,
        secondary_structures_T: float,
        Tm_parameters: dict,
        Tm_chem_correction_parameters: dict,
    ):
        """
        Computes various attributes for oligos and stores them in the provided oligo database.

        :param oligo_database: The database containing oligos for which attributes are to be computed.
        :type oligo_database: OligoDatabase
        :param secondary_structures_T: Temperature at which secondary structure delta G is calculated.
        :type secondary_structures_T: float
        :param Tm_parameters: Parameters to calculate melting temperature.
        :type Tm_parameters: dict
        :param Tm_chem_correction_parameters: Chemical correction parameters for melting temperature calculation.
        :type Tm_chem_correction_parameters: dict
        :return: Updated oligo database with new attributes.
        :rtype: OligoDatabase
        """
        oligo_database = self.oligo_attributes.calculate_oligo_length(oligo_database=oligo_database)
        oligo_database = self.oligo_attributes.calculate_GC_content(
            oligo_database=oligo_database, sequence_type="oligo"
        )
        oligo_database = self.oligo_attributes.calculate_TmNN(
            oligo_database=oligo_database,
            sequence_type="oligo",
            Tm_parameters=Tm_parameters,
            Tm_chem_correction_parameters=Tm_chem_correction_parameters,
        )
        oligo_database = self.oligo_attributes.calculate_num_targeted_transcripts(
            oligo_database=oligo_database
        )
        oligo_database = self.oligo_attributes.calculate_isoform_consensus(oligo_database=oligo_database)
        oligo_database = self.oligo_attributes.calculate_length_selfcomplement(
            oligo_database=oligo_database, sequence_type="oligo"
        )
        oligo_database = self.oligo_attributes.calculate_DG_secondary_structure(
            oligo_database=oligo_database, sequence_type="oligo", T=secondary_structures_T
        )

        return oligo_database

    def generate_output(self, oligo_database: OligoDatabase, top_n_sets: int):

        attributes = [
            "chromosome",
            "start",
            "end",
            "strand",
            "oligo",
            "target",
            "length",
            "GC_content",
            "TmNN",
            "num_targeted_transcripts",
            "isoform_consensus",
            "length_selfcomplement",
            "DG_secondary_structure" "source",
            "species",
            "annotation_release",
            "genome_assembly",
            "regiontype",
            "gene_id",
            "transcript_id",
            "exon_number",
        ]
        oligo_database.write_oligosets_to_yaml(
            attributes=attributes, top_n_sets=top_n_sets, ascending=True, filename="oligo_seq_probes.yml"
        )


############################################
# Oligo-seq Designer Pipeline
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

    ##### process parameters #####
    # preprocess melting temperature params
    Tm_parameters = config["Tm_parameters"]
    Tm_parameters["nn_table"] = getattr(mt, Tm_parameters["nn_table"])
    Tm_parameters["tmm_table"] = getattr(mt, Tm_parameters["tmm_table"])
    Tm_parameters["imm_table"] = getattr(mt, Tm_parameters["imm_table"])
    Tm_parameters["de_table"] = getattr(mt, Tm_parameters["de_table"])

    ##### read the genes file #####
    if config["file_regions"] is None:
        warnings.warn(
            "No gene list file was provided! All genes from fasta file are used to generate the probes. This chioce can use a lot of resources."
        )
        gene_ids = None
    else:
        with open(config["file_regions"]) as handle:
            lines = handle.readlines()
            # ensure that the list contains unique gene ids
            gene_ids = list(set([line.rstrip() for line in lines]))

    ##### initialize probe designer pipeline #####
    pipeline = OligoSeqProbeDesigner(
        file_regions=config["file_regions"],
        write_intermediate_steps=config["write_intermediate_steps"],
        dir_output=config["dir_output"],
        n_jobs=config["n_jobs"],
    )

    ##### create oligo database #####
    oligo_database, file_database = pipeline.create_oligo_database(
        gene_ids=gene_ids,
        oligo_length_min=config["oligo_length_min"],
        oligo_length_max=config["oligo_length_max"],
        files_fasta_oligo_database=config["files_fasta_oligo_database"],
        # we should have at least "min_oligoset_size" oligos per gene to create one set
        min_oligos_per_region=config["min_oligoset_size"],
    )

    ##### filter oligos by property #####
    oligo_database, file_database = pipeline.filter_by_property(
        oligo_database,
        config["GC_content_min"],
        config["GC_content_max"],
        config["Tm_min"],
        Tm_max=config["Tm_max"],
        secondary_structures_T=config["secondary_structures_T"],
        secondary_structures_threshold_deltaG=config["secondary_structures_threshold_deltaG"],
        homopolymeric_base_n=config["homopolymeric_base_n"],
        max_len_selfcomp=config["max_len_selfcomp"],
        Tm_parameters=Tm_parameters,
        Tm_chem_correction_parameters=config["Tm_chem_correction_parameters"],
    )

    # ##### filter oligos by specificity #####
    oligo_database, file_database = pipeline.filter_by_specificity(
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
    )

    ##### create oligo sets #####
    oligo_database, file_database, dir_oligosets = pipeline.create_oligo_sets(
        oligo_database=oligo_database,
        Tm_min=config["Tm_min"],
        Tm_opt=config["Tm_opt"],
        Tm_max=config["Tm_max"],
        Tm_parameters=Tm_parameters,
        Tm_chem_correction_parameters=config["Tm_chem_correction_parameters"],
        GC_content_min=config["GC_content_min"],
        GC_content_opt=config["GC_content_opt"],
        GC_content_max=config["GC_content_max"],
        opt_oligoset_size=config["opt_oligoset_size"],
        min_oligoset_size=config["min_oligoset_size"],
        max_oligos=config["max_graph_size"],
        n_sets=config["n_sets"],
        distance_between_oligos=config["distance_between_oligos"],
    )

    ##### compute all required attributes #####
    logging.info("Computing Oligo Attributes")
    oligo_database = pipeline.compute_oligo_attributes(
        oligo_database=oligo_database,
        secondary_structures_T=config["secondary_structures_T"],
        Tm_parameters=Tm_parameters,
        Tm_chem_correction_parameters=config["Tm_chem_correction_parameters"],
    )

    pipeline.generate_output(oligo_database=oligo_database, top_n_sets=config["top_n_sets"])

    logging.info(f"Oligo sets were saved in {dir_oligosets}")
    logging.info("##### End of the pipeline. #####")


if __name__ == "__main__":
    main()
