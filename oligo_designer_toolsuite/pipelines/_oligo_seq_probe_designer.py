############################################
# imports
############################################

import logging
import os
import shutil
import warnings
from datetime import datetime
from pathlib import Path
from typing import List, Tuple

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
    GreedySelectionPolicy,
    OligosetGeneratorIndependentSet,
)
from oligo_designer_toolsuite.oligo_specificity_filter import (
    BlastNFilter,
    BowtieFilter,
    CrossHybridizationFilter,
    ExactMatchFilter,
    HybridizationProbabilityFilter,
    RemoveByLargerRegionPolicy,
    RemoveAllPolicy,
    SpecificityFilter,
)
from oligo_designer_toolsuite.pipelines._utils import (
    base_parser,
    base_log_parameters,
    pipeline_step_basic,
    check_content_oligo_database,
)
from oligo_designer_toolsuite.sequence_generator import OligoSequenceGenerator

############################################
# Oligo-Seq Probe Designer
############################################


class OligoSeqProbeDesigner:
    """
    A class for designing hybridization probes for the Oligo-Seq experiments.

    An oligo-seq probe is an oligo hybridization probe, which is optimized for
    probe-based targeted sequencing to measure RNA expression.

    :param write_intermediate_steps: Whether to save intermediate results during the probe design pipeline.
    :type write_intermediate_steps: bool
    :param dir_output: Directory path where output files and logs will be saved.
    :type dir_output: str
    :param n_jobs: Number of parallel jobs to use for computationally intensive tasks.
    :type n_jobs: int
    """

    def __init__(self, write_intermediate_steps: bool, dir_output: str, n_jobs: int) -> None:
        """Constructor for the OligoSeqProbeDesigner class."""

        ##### create the output folder #####
        self.dir_output = os.path.abspath(dir_output)
        Path(self.dir_output).mkdir(parents=True, exist_ok=True)

        ##### setup logger #####
        timestamp = datetime.now()
        file_logger = os.path.join(
            self.dir_output,
            f"log_oligoseq_probe_designer_{timestamp.year}-{timestamp.month}-{timestamp.day}-{timestamp.hour}-{timestamp.minute}.txt",
        )
        logging.getLogger("log_name")
        logging.basicConfig(
            format="%(asctime)s [%(levelname)s] %(message)s",
            level=logging.NOTSET,
            handlers=[logging.FileHandler(file_logger)],
        )
        logging.captureWarnings(True)
        logging.info("--------------START PIPELINE--------------")

        ##### set class parameters #####
        self.write_intermediate_steps = write_intermediate_steps
        self.n_jobs = n_jobs
        self.oligo_attributes_calculator = OligoAttributes()
        self.set_developer_parameters()

    def set_developer_parameters(
        self,
        target_probe_hybridization_probability_alignment_method: str = "blastn",
        target_probe_hybridization_probability_blastn_search_parameters: dict = {
            "perc_identity": 80,
            "strand": "minus",
            "word_size": 10,
        },
        target_probe_hybridization_probability_blastn_hit_parameters: dict = {"coverage": 50},
        target_probe_hybridization_probability_bowtie_search_parameters: dict = {
            "-v": 3,
            "--nofw": "",
        },
        target_probe_hybridization_probability_bowtie_hit_parameters: dict = None,
        target_probe_cross_hybridization_alignment_method: str = "blastn",
        target_probe_cross_hybridization_blastn_search_parameters: dict = {
            "perc_identity": 80,
            "strand": "minus",
            "word_size": 10,
        },
        target_probe_cross_hybridization_blastn_hit_parameters: dict = {"coverage": 50},
        target_probe_cross_hybridization_bowtie_search_parameters: dict = {
            "-v": 3,
            "--nofw": "",
        },
        target_probe_cross_hybridization_bowtie_hit_parameters: dict = None,
        max_graph_size: int = 5000,
        n_attempts: int = 100000,
        heuristic: bool = True,
        heuristic_n_attempts: int = 100,
        target_probe_Tm_parameters: dict = {
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
            "saltcorr": 7,
            "Na": 1000,
            "K": 0,
            "Tris": 0,
            "Mg": 0,
            "dNTPs": 0,
        },
        target_probe_Tm_chem_correction_parameters: dict = {
            "DMSO": 0,
            "fmd": 20,
            "DMSOfactor": 0.75,
            "fmdfactor": 0.65,
            "fmdmethod": 1,
            "GC": None,
        },
        target_probe_Tm_salt_correction_parameters: dict = None,
    ):
        """
        Set developer-specific parameters for Oligo-Seq probe designer pipeline.
        These parameters can be used to customize and fine-tune the pipeline.

        :param target_probe_hybridization_probability_alignment_method: Alignment method for computing hybridization probabilities,
            either 'blastn' or 'bowtie'. Defaults to 'blastn'.
        :type target_probe_hybridization_probability_alignment_method: str
        :param target_probe_hybridization_probability_blastn_search_parameters: Parameters for BLASTN search.
        :type target_probe_hybridization_probability_blastn_search_parameters: dict
        :param target_probe_hybridization_probability_blastn_hit_parameters: Parameters for BLASTN hits.
        :type target_probe_hybridization_probability_blastn_hit_parameters: dict
        :param target_probe_hybridization_probability_bowtie_search_parameters: Parameters for Bowtie search.
        :type target_probe_hybridization_probability_bowtie_search_parameters: dict
        :param target_probe_hybridization_probability_bowtie_hit_parameters: Parameters for Bowtie hits.
        :type target_probe_hybridization_probability_bowtie_hit_parameters: dict
        :param target_probe_cross_hybridization_alignment_method: Alignment method for cross-hybridization filtering,
            either 'blastn' or 'bowtie'. Defaults to 'blastn'.
        :type target_probe_cross_hybridization_alignment_method: str
        :param target_probe_cross_hybridization_blastn_search_parameters: Parameters for BLASTN cross-hybridization search.
        :type target_probe_cross_hybridization_blastn_search_parameters: dict
        :param target_probe_cross_hybridization_blastn_hit_parameters: Parameters for BLASTN cross-hybridization hits.
        :type target_probe_cross_hybridization_blastn_hit_parameters: dict
        :param target_probe_cross_hybridization_bowtie_search_parameters: Parameters for Bowtie cross-hybridization search.
        :type target_probe_cross_hybridization_bowtie_search_parameters: dict
        :param target_probe_cross_hybridization_bowtie_hit_parameters: Parameters for Bowtie cross-hybridization hits.
        :type target_probe_cross_hybridization_bowtie_hit_parameters: dict
        :param max_graph_size: Maximum size of the graph used in set selection, defaults to 5000.
        :type max_graph_size: int
        :param n_attempts: Maximum number of attempts for selecting oligo sets, defaults to 100000.
        :type n_attempts: int
        :param heuristic: Whether to apply heuristic methods in oligo set selection, defaults to True.
        :type heuristic: bool
        :param heuristic_n_attempts: Maximum number of attempts for heuristic selecting oligo sets, defaults to 100.
        :type heuristic_n_attempts: int
        :param target_probe_Tm_parameters: Parameters for melting temperature computation.
            For using Bio.SeqUtils.MeltingTemp default parameters set to ``{}``. For more information on parameters,
            see: https://biopython.org/docs/1.75/api/Bio.SeqUtils.MeltingTemp.html#Bio.SeqUtils.MeltingTemp.Tm_NN
        :type target_probe_Tm_parameters: dict
        :param target_probe_Tm_chem_correction_parameters: Chemical correction parameters for melting temperature.
            For using Bio.SeqUtils.MeltingTemp default parameters set to ``{}``. For more information on parameters,
            see: https://biopython.org/docs/1.75/api/Bio.SeqUtils.MeltingTemp.html#Bio.SeqUtils.MeltingTemp.salt_correction
        :type target_probe_Tm_chem_correction_parameters: dict
        :param target_probe_Tm_salt_correction_parameters: Salt correction parameters for melting temperature.
            For using Bio.SeqUtils.MeltingTemp default parameters set to ``{}``. For more information on parameters,
            see: https://biopython.org/docs/1.75/api/Bio.SeqUtils.MeltingTemp.html#Bio.SeqUtils.MeltingTemp.chem_correction
        :type target_probe_Tm_salt_correction_parameters: dict
        """
        ### Parameters for the specificity filters
        self.target_probe_hybridization_probability_alignment_method = (
            target_probe_hybridization_probability_alignment_method
        )
        if target_probe_hybridization_probability_alignment_method == "blastn":
            # Specificity filter with BlastN
            self.target_probe_hybridization_probability_search_parameters = (
                target_probe_hybridization_probability_blastn_search_parameters
            )
            self.target_probe_hybridization_probability_hit_parameters = (
                target_probe_hybridization_probability_blastn_hit_parameters
            )
        elif target_probe_hybridization_probability_alignment_method == "bowtie":
            # Specificity filter with Bowtie
            self.target_probe_hybridization_probability_search_parameters = (
                target_probe_hybridization_probability_bowtie_search_parameters
            )
            self.target_probe_hybridization_probability_hit_parameters = (
                target_probe_hybridization_probability_bowtie_hit_parameters
            )

        self.target_probe_cross_hybridization_alignment_method = (
            target_probe_cross_hybridization_alignment_method
        )
        if target_probe_cross_hybridization_alignment_method == "blastn":
            # Crosshybridization filter with BlastN
            self.target_probe_cross_hybridization_search_parameters = (
                target_probe_cross_hybridization_blastn_search_parameters
            )
            self.target_probe_cross_hybridization_hit_parameters = (
                target_probe_cross_hybridization_blastn_hit_parameters
            )
        elif target_probe_cross_hybridization_alignment_method == "bowtie":
            # Crosshybridization filter with Bowtie
            self.target_probe_cross_hybridization_search_parameters = (
                target_probe_cross_hybridization_bowtie_search_parameters
            )
            self.target_probe_cross_hybridization_hit_parameters = (
                target_probe_cross_hybridization_bowtie_hit_parameters
            )

        ### Parameters for the Oligo set selection
        self.max_graph_size = max_graph_size
        self.heuristic = heuristic
        self.n_attempts = n_attempts
        self.heuristic_n_attempts = heuristic_n_attempts

        ### Parameters for Melting Temperature
        # The melting temperature is used in 2 different stages (property filters and padlock detection probe design), where a few parameters are shared and the others differ.
        # parameters for melting temperature -> for more information on parameters, see: https://biopython.org/docs/1.75/api/Bio.SeqUtils.MeltingTemp.html#Bio.SeqUtils.MeltingTemp.Tm_NN

        # preprocess melting temperature params
        target_probe_Tm_parameters["nn_table"] = getattr(mt, target_probe_Tm_parameters["nn_table"])
        target_probe_Tm_parameters["tmm_table"] = getattr(mt, target_probe_Tm_parameters["tmm_table"])
        target_probe_Tm_parameters["imm_table"] = getattr(mt, target_probe_Tm_parameters["imm_table"])
        target_probe_Tm_parameters["de_table"] = getattr(mt, target_probe_Tm_parameters["de_table"])

        ## target probe
        self.target_probe_Tm_parameters = target_probe_Tm_parameters
        self.target_probe_Tm_chem_correction_parameters = target_probe_Tm_chem_correction_parameters
        self.target_probe_Tm_salt_correction_parameters = target_probe_Tm_salt_correction_parameters

    def design_target_probes(
        self,
        files_fasta_target_probe_database: list,
        files_fasta_reference_database_targe_probe: list,
        gene_ids: list = None,
        target_probe_length_min: int = 26,
        target_probe_length_max: int = 30,
        target_probe_split_region: int = 4,
        target_probe_targeted_exons: list = ["1", "2", "3"],
        target_probe_isoform_consensus: float = 0,
        target_probe_GC_content_min: float = 45,
        target_probe_GC_content_opt: float = 55,
        target_probe_GC_content_max: float = 65,
        target_probe_Tm_min: float = 50,
        target_probe_Tm_opt: float = 60,
        target_probe_Tm_max: float = 70,
        target_probe_secondary_structures_T: int = 37,
        target_probe_secondary_structures_threshold_deltaG: int = 0,
        target_probe_homopolymeric_base_n: dict = {"A": 6, "T": 6, "C": 6, "G": 6},
        target_probe_max_len_selfcomplement: int = 10,
        target_probe_hybridization_probability_threshold: float = 0.001,
        target_probe_read_length_bias: int = 20,
        target_probe_GC_weight: float = 1,
        target_probe_Tm_weight: float = 1,
        set_size_min: int = 3,
        set_size_opt: int = 5,
        distance_between_target_probes: int = 0,
        n_sets: int = 100,
    ):
        """
        Design target probes based on specified parameters, including property and specificity filters.
        The designed probes are organized into sets based on customizable constraints.

        :param files_fasta_target_probe_database: FASTA files containing the target probe database.
        :type files_fasta_target_probe_database: list
        :param files_fasta_reference_database_targe_probe: FASTA files containing the reference database.
        :type files_fasta_reference_database_targe_probe: list
        :param gene_ids: List of gene IDs to target, or None to target all genes.
        :type gene_ids: list, optional
        :param target_probe_length_min: Minimum length of target probes, defaults to 26.
        :type target_probe_length_min: int, optional
        :param target_probe_length_max: Maximum length of target probes, defaults to 30.
        :type target_probe_length_max: int, optional
        :param target_probe_split_region: The number of bases required on each side of a split sequence (e.g. exon junctions) to include it, defaults to 4.
        :type target_probe_split_region: int, optional
        :param target_probe_targeted_exons: Exons to target, defaults to ["1", "2", "3"].
        :type target_probe_targeted_exons: list, optional
        :param target_probe_isoform_consensus: Minimum isoform consensus, defaults to 0.
        :type target_probe_isoform_consensus: float, optional
        :param target_probe_GC_content_min: Minimum GC content for probes, defaults to 45.
        :type target_probe_GC_content_min: float, optional
        :param target_probe_GC_content_opt: Optimal GC content for probes, defaults to 55.
        :type target_probe_GC_content_opt: float, optional
        :param target_probe_GC_content_max: Maximum GC content for probes, defaults to 65.
        :type target_probe_GC_content_max: float, optional
        :param target_probe_Tm_min: Minimum melting temperature for probes, defaults to 50.
        :type target_probe_Tm_min: float, optional
        :param target_probe_Tm_opt: Optimal melting temperature for probes, defaults to 60.
        :type target_probe_Tm_opt: float, optional
        :param target_probe_Tm_max: Maximum melting temperature for probes, defaults to 70.
        :type target_probe_Tm_max: float, optional
        :param target_probe_secondary_structures_T: Temperature for secondary structure evaluation, defaults to 37.
        :type target_probe_secondary_structures_T: int, optional
        :param target_probe_secondary_structures_threshold_deltaG: Threshold delta G for secondary structures, defaults to 0.
        :type target_probe_secondary_structures_threshold_deltaG: int, optional
        :param target_probe_homopolymeric_base_n: Maximum number of consecutive homopolymeric bases, defaults to {"A": 6, "T": 6, "C": 6, "G": 6}.
        :type target_probe_homopolymeric_base_n: dict, optional
        :param target_probe_max_len_selfcomplement: Maximum length of self-complementary sequences, defaults to 10.
        :type target_probe_max_len_selfcomplement: int, optional
        :param target_probe_hybridization_probability_threshold: Threshold for hybridization probability, defaults to 0.001.
        :type target_probe_hybridization_probability_threshold: float, optional
        :param target_probe_read_length_bias: Minimum length of sequencing reads to account for shorter sequencing reads,
            i.e. first <target_probe_read_length_bias> bases of an oligo should not match the <target_probe_read_length_bias> bases of another oligo
        :type target_probe_read_length_bias: int
        :param target_probe_GC_weight: Weight for GC content in set scoring, defaults to 1.
        :type target_probe_GC_weight: float, optional
        :param target_probe_Tm_weight: Weight for melting temperature in set scoring, defaults to 1.
        :type target_probe_Tm_weight: float, optional
        :param set_size_min: Minimum size of probe sets, defaults to 3.
        :type set_size_min: int, optional
        :param set_size_opt: Optimal size of probe sets, defaults to 5.
        :type set_size_opt: int, optional
        :param distance_between_target_probes: Minimum genomic distance between probes in a set, defaults to 0.
        :type distance_between_target_probes: int, optional
        :param n_sets: Number of probe sets to generate, defaults to 100.
        :type n_sets: int, optional
        :return: The designed oligo database with probe sets.
        :rtype: OligoDatabase
        """
        target_probe_designer = TargetProbeDesigner(self.dir_output, self.n_jobs)

        oligo_database = target_probe_designer.create_oligo_database(
            gene_ids=gene_ids,
            oligo_length_min=target_probe_length_min,
            oligo_length_max=target_probe_length_max,
            split_region=target_probe_split_region,
            files_fasta_oligo_database=files_fasta_target_probe_database,
            min_oligos_per_gene=set_size_min,
            isoform_consensus=target_probe_isoform_consensus,
            targeted_exons=target_probe_targeted_exons,
        )
        check_content_oligo_database(oligo_database)

        if self.write_intermediate_steps:
            dir_database = oligo_database.save_database(dir_database="1_db_probes_initial")
            print(f"Saved probe database for step 1 (Create Database) in directory {dir_database}")

        oligo_database = target_probe_designer.filter_by_property(
            oligo_database=oligo_database,
            GC_content_min=target_probe_GC_content_min,
            GC_content_max=target_probe_GC_content_max,
            Tm_min=target_probe_Tm_min,
            Tm_max=target_probe_Tm_max,
            secondary_structures_T=target_probe_secondary_structures_T,
            secondary_structures_threshold_deltaG=target_probe_secondary_structures_threshold_deltaG,
            homopolymeric_base_n=target_probe_homopolymeric_base_n,
            max_len_selfcomplement=target_probe_max_len_selfcomplement,
            Tm_parameters=self.target_probe_Tm_parameters,
            Tm_chem_correction_parameters=self.target_probe_Tm_chem_correction_parameters,
            Tm_salt_correction_parameters=self.target_probe_Tm_salt_correction_parameters,
        )
        check_content_oligo_database(oligo_database)

        if self.write_intermediate_steps:
            dir_database = oligo_database.save_database(dir_database="2_db_probes_property_filter")
            print(f"Saved probe database for step 2 (Property Filters) in directory {dir_database}")

        oligo_database = target_probe_designer.filter_by_specificity(
            oligo_database=oligo_database,
            files_fasta_reference_database=files_fasta_reference_database_targe_probe,
            cross_hybridization_alignment_method=self.target_probe_cross_hybridization_alignment_method,
            cross_hybridization_search_parameters=self.target_probe_cross_hybridization_search_parameters,
            cross_hybridization_hit_parameters=self.target_probe_cross_hybridization_hit_parameters,
            hybridization_probability_alignment_method=self.target_probe_hybridization_probability_alignment_method,
            hybridization_probability_search_parameters=self.target_probe_hybridization_probability_search_parameters,
            hybridization_probability_hit_parameters=self.target_probe_hybridization_probability_hit_parameters,
            hybridization_probability_threshold=target_probe_hybridization_probability_threshold,
            target_probe_read_length_bias=target_probe_read_length_bias,
        )
        check_content_oligo_database(oligo_database)

        if self.write_intermediate_steps:
            dir_database = oligo_database.save_database(dir_database="3_db_probes_specificity_filter")
            print(f"Saved probe database for step 3 (Specificity Filters) in directory {dir_database}")

        oligo_database = target_probe_designer.create_oligo_sets(
            oligo_database=oligo_database,
            GC_content_min=target_probe_GC_content_min,
            GC_content_opt=target_probe_GC_content_opt,
            GC_content_max=target_probe_GC_content_max,
            GC_weight=target_probe_GC_weight,
            Tm_min=target_probe_Tm_min,
            Tm_opt=target_probe_Tm_opt,
            Tm_max=target_probe_Tm_max,
            Tm_weight=target_probe_Tm_weight,
            Tm_parameters=self.target_probe_Tm_parameters,
            Tm_chem_correction_parameters=self.target_probe_Tm_chem_correction_parameters,
            Tm_salt_correction_parameters=self.target_probe_Tm_salt_correction_parameters,
            set_size_min=set_size_min,
            set_size_opt=set_size_opt,
            distance_between_oligos=distance_between_target_probes,
            n_sets=n_sets,
            max_graph_size=self.max_graph_size,
            n_attempts=self.n_attempts,
            heuristic=self.heuristic,
            heuristic_n_attempts=self.heuristic_n_attempts,
        )
        check_content_oligo_database(oligo_database)

        if self.write_intermediate_steps:
            dir_database = oligo_database.save_database(dir_database="4_db_probes_probesets")
            dir_probesets = oligo_database.write_oligosets_to_table()
            print(
                f"Saved probe database for step 4 (Specificity Filters) in directory {dir_database} and probeset table in directory {dir_probesets}"
            )

        return oligo_database

    def generate_output(
        self,
        oligo_database: OligoDatabase,
        top_n_sets: int = 3,
        attributes=[
            "source",
            "species",
            "annotation_release",
            "genome_assembly",
            "regiontype",
            "gene_id",
            "transcript_id",
            "exon_number",
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
            "number_total_transcripts",
            "isoform_consensus",
            "length_selfcomplement",
            "DG_secondary_structure",
        ],
    ) -> None:
        """
        Generate the final output files for the Oligo-Seq probe design pipeline.

        :param oligo_database: The oligo database containing final designed probes and attributes.
        :type oligo_database: OligoDatabase
        :param top_n_sets: Number of top probe sets to include in the output, defaults to 3.
        :type top_n_sets: int
        :param attributes: List of attributes to include in the output files, defaults to a predefined list of probe attributes.
        :type attributes: list

        :return: None
        """
        oligo_database = self.oligo_attributes_calculator.calculate_oligo_length(
            oligo_database=oligo_database
        )
        oligo_database = self.oligo_attributes_calculator.calculate_GC_content(
            oligo_database=oligo_database, sequence_type="oligo"
        )
        oligo_database = self.oligo_attributes_calculator.calculate_TmNN(
            oligo_database=oligo_database,
            sequence_type="oligo",
            Tm_parameters=self.target_probe_Tm_parameters,
            Tm_chem_correction_parameters=self.target_probe_Tm_chem_correction_parameters,
            Tm_salt_correction_parameters=self.target_probe_Tm_salt_correction_parameters,
        )
        oligo_database = self.oligo_attributes_calculator.calculate_num_targeted_transcripts(
            oligo_database=oligo_database
        )
        oligo_database = self.oligo_attributes_calculator.calculate_isoform_consensus(
            oligo_database=oligo_database
        )
        oligo_database = self.oligo_attributes_calculator.calculate_length_selfcomplement(
            oligo_database=oligo_database, sequence_type="oligo"
        )

        oligo_database.write_oligosets_to_yaml(
            attributes=attributes, top_n_sets=top_n_sets, ascending=True, filename="oligo_seq_probesets"
        )
        oligo_database.write_database_to_table(
            attributes=attributes, flatten_attribute=True, filename="oligo_seq_probes"
        )

        logging.info("--------------END PIPELINE--------------")


############################################
# Oligo-Seq Target Probe Designer
############################################


class TargetProbeDesigner:
    """
    A class for designing target probes for Oligo-Seq experiments.
    This class provides methods for creating, filtering, and scoring oligos based
    on specific properties and designing oligo sets for targeted probes.

    :param dir_output: Directory path where output files and intermediate results will be saved.
    :type dir_output: str
    :param n_jobs: Number of parallel jobs to use for computationally intensive tasks.
    :type n_jobs: int
    """

    def __init__(self, dir_output: str, n_jobs: int) -> None:
        """Constructor for the TargetProbeDesigner class."""

        ##### create the output folder #####
        self.dir_output = os.path.abspath(dir_output)
        self.subdir_db_probes = "db_probes"
        self.subdir_db_reference = "db_reference"

        self.n_jobs = n_jobs
        self.oligo_attributes_calculator = OligoAttributes()

    @pipeline_step_basic(step_name="Create Database")
    def create_oligo_database(
        self,
        gene_ids: list,
        oligo_length_min: int,
        oligo_length_max: int,
        split_region: int,
        files_fasta_oligo_database: list[str],
        min_oligos_per_gene: int,
        isoform_consensus: float,
        targeted_exons: list[str],
    ) -> OligoDatabase:
        """
        Creates an oligo database by generating sequences using a sliding window approach
        and filtering based on specified criteria.

        :param gene_ids: List of gene identifiers for which oligos should be generated.
                        If None, all genes in the input fasta file are used.
        :type gene_ids: list
        :param oligo_length_min: Minimum length of oligos to generate.
        :type oligo_length_min: int
        :param oligo_length_max: Maximum length of oligos to generate.
        :type oligo_length_max: int
        :param split_region: The number of bases required on each side of a split sequence (e.g. exon junctions) to include it.
        :type split_region: int
        :param files_fasta_oligo_database: List of FASTA files to use for sequence generation.
        :type files_fasta_oligo_database: list[str]
        :param min_oligos_per_gene: Minimum number of oligos required per gene.
        :type min_oligos_per_gene: int
        :param isoform_consensus: Threshold for isoform consensus filtering.
        :type isoform_consensus: float
        :param targeted_exons: List of exon numbers to target.
        :type targeted_exons: list[str]
        :return: The generated oligo database.
        :rtype: OligoDatabase
        """
        ##### creating the oligo sequences #####
        oligo_sequences = OligoSequenceGenerator(dir_output=self.dir_output)
        oligo_fasta_file = oligo_sequences.create_sequences_sliding_window(
            files_fasta_in=files_fasta_oligo_database,
            length_interval_sequences=(oligo_length_min, oligo_length_max),
            split_region=split_region,
            region_ids=gene_ids,
            n_jobs=self.n_jobs,
        )

        ##### creating the oligo database #####
        oligo_database = OligoDatabase(
            min_oligos_per_region=min_oligos_per_gene,
            write_regions_with_insufficient_oligos=True,
            lru_db_max_in_memory=self.n_jobs * 2 + 2,
            database_name=self.subdir_db_probes,
            dir_output=self.dir_output,
            n_jobs=1,
        )
        oligo_database.load_database_from_fasta(
            files_fasta=oligo_fasta_file,
            sequence_type="target",
            database_overwrite=True,
            region_ids=gene_ids,
        )
        oligo_database = self.oligo_attributes_calculator.calculate_reverse_complement_sequence(
            oligo_database=oligo_database, sequence_type="target", sequence_type_reverse_complement="oligo"
        )

        ##### pre-filter oligo database for certain attributes #####
        oligo_database = self.oligo_attributes_calculator.calculate_isoform_consensus(
            oligo_database=oligo_database
        )
        oligo_database.filter_database_by_attribute_threshold(
            attribute_name="isoform_consensus",
            attribute_thr=isoform_consensus,
            remove_if_smaller_threshold=True,
        )
        oligo_database.filter_database_by_attribute_category(
            attribute_name="exon_number",
            attribute_category=targeted_exons,
            remove_if_equals_category=False,
        )
        oligo_database.remove_regions_with_insufficient_oligos(pipeline_step="Pre-Filters")

        dir = oligo_sequences.dir_output
        shutil.rmtree(dir) if os.path.exists(dir) else None

        return oligo_database

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
        max_len_selfcomplement: int,
        Tm_parameters: dict,
        Tm_chem_correction_parameters: dict,
        Tm_salt_correction_parameters: dict,
    ) -> OligoDatabase:
        """
        Filter the oligo database based on various sequence properties.

        :param oligo_database: The oligo database to filter.
        :type oligo_database: OligoDatabase
        :param GC_content_min: Minimum GC content for filtering.
        :type GC_content_min: int
        :param GC_content_max: Maximum GC content for filtering.
        :type GC_content_max: int
        :param Tm_min: Minimum melting temperature for filtering.
        :type Tm_min: int
        :param Tm_max: Maximum melting temperature for filtering.
        :type Tm_max: int
        :param secondary_structures_T: Temperature for secondary structure analysis.
        :type secondary_structures_T: float
        :param secondary_structures_threshold_deltaG: Threshold for secondary structure deltaG.
        :type secondary_structures_threshold_deltaG: float
        :param homopolymeric_base_n: Threshold for homopolymeric runs for each base.
        :type homopolymeric_base_n: str
        :param max_len_selfcomplement: Maximum self-complementary length allowed.
        :type max_len_selfcomplement: int
        :param Tm_parameters: Parameters for melting temperature calculations.
        :type Tm_parameters: dict
        :param Tm_chem_correction_parameters: Parameters for chemical corrections in Tm calculations.
        :type Tm_chem_correction_parameters: dict
        :param Tm_salt_correction_parameters: Parameters for salt corrections in Tm calculations.
        :type Tm_salt_correction_parameters: dict
        :return: The filtered oligo database.
        :rtype: OligoDatabase
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
            Tm_salt_correction_parameters=Tm_salt_correction_parameters,
        )
        secondary_sctructure = SecondaryStructureFilter(
            T=secondary_structures_T,
            thr_DG=secondary_structures_threshold_deltaG,
        )
        homopolymeric_runs = HomopolymericRunsFilter(
            base_n=homopolymeric_base_n,
        )
        self_comp = SelfComplementFilter(
            max_len_selfcomplement=max_len_selfcomplement,
        )

        filters = [
            hard_masked_sequences,
            soft_masked_sequences,
            homopolymeric_runs,
            gc_content,
            melting_temperature,
            self_comp,
            secondary_sctructure,
            homopolymeric_runs,
        ]

        # initialize the preoperty filter class
        property_filter = PropertyFilter(filters=filters)

        # filter the database
        oligo_database = property_filter.apply(
            oligo_database=oligo_database,
            sequence_type="oligo",
            n_jobs=self.n_jobs,
        )

        return oligo_database

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
        target_probe_read_length_bias: int,
    ) -> OligoDatabase:
        """
        Filter the oligo database based on sequence specificity to remove sequences that
        cross-hybridize to other oligos or hybridization to other genomic regions.

        :param oligo_database: The oligo database to filter.
        :type oligo_database: OligoDatabase
        :param files_fasta_reference_database: List of FASTA files for the reference database.
        :type files_fasta_reference_database: List[str]
        :param cross_hybridization_alignment_method: Alignment method for cross-hybridization analysis.
        :type cross_hybridization_alignment_method: str
        :param cross_hybridization_search_parameters: Search parameters for cross-hybridization analysis.
        :type cross_hybridization_search_parameters: dict
        :param cross_hybridization_hit_parameters: Hit parameters for cross-hybridization analysis.
        :type cross_hybridization_hit_parameters: dict
        :param hybridization_probability_alignment_method: Alignment method for hybridization probability analysis.
        :type hybridization_probability_alignment_method: str
        :param hybridization_probability_search_parameters: Search parameters for hybridization probability analysis.
        :type hybridization_probability_search_parameters: dict
        :param hybridization_probability_hit_parameters: Hit parameters for hybridization probability analysis.
        :type hybridization_probability_hit_parameters: dict
        :param hybridization_probability_threshold: Threshold for hybridization probability filtering.
        :type hybridization_probability_threshold: float
        :param target_probe_read_length_bias: Minimum length of sequencing reads to account for shorter sequencing reads,
            i.e. first <target_probe_read_length_bias> bases of an oligo should not match the <target_probe_read_length_bias> bases of another oligo
        :type target_probe_read_length_bias: int
        :return: The filtered oligo database.
        :rtype: OligoDatabase
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
        # remove sequences that could cause read length biases because the first
        # <target_probe_read_length_bias> bases of both sequences match
        oligo_database = self.oligo_attributes_calculator.calculate_shortened_sequence(
            oligo_database=oligo_database,
            sequence_length=target_probe_read_length_bias,
            sequence_type="oligo",
            reverse=False,
        )

        exact_matches = ExactMatchFilter(policy=RemoveAllPolicy(), filter_name="exact_match_read_length_bias")
        specificity_filter = SpecificityFilter(filters=[exact_matches])
        oligo_database = specificity_filter.apply(
            sequence_type="oligo_short",
            oligo_database=oligo_database,
            n_jobs=self.n_jobs,
        )

        # removing duplicated oligos from the region with the most oligos
        # this step can be redundant with the hybridization probability filter
        # but improves runtuíme as it pre-filters such sequences
        exact_matches = ExactMatchFilter(policy=RemoveAllPolicy(), filter_name="exact_match")

        # remove oligos that potentially cross-hybridize with other oligos in the set
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

        # uses an alignment method to check for off-target hits wrt to reference sequences
        # and refines those hits with a hybridization probability model
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

        # run all filters specified above
        filters = [exact_matches, cross_hybridization, hybridization_probability]
        specificity_filter = SpecificityFilter(filters=filters)
        oligo_database = specificity_filter.apply(
            sequence_type="oligo",
            oligo_database=oligo_database,
            reference_database=reference_database,
            n_jobs=self.n_jobs,
        )

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

        return oligo_database

    @pipeline_step_basic(step_name="Oligo Selection")
    def create_oligo_sets(
        self,
        oligo_database: OligoDatabase,
        GC_content_min: float,
        GC_content_opt: float,
        GC_content_max: float,
        GC_weight: float,
        Tm_min: float,
        Tm_opt: float,
        Tm_max: float,
        Tm_weight: float,
        Tm_parameters: dict,
        Tm_chem_correction_parameters: dict,
        Tm_salt_correction_parameters: dict,
        set_size_opt: int,
        set_size_min: int,
        distance_between_oligos: int,
        n_sets: int,
        max_graph_size: int,
        n_attempts: int,
        heuristic: bool,
        heuristic_n_attempts: int,
    ) -> OligoDatabase:
        """
        Create optimal oligo sets based on weighted scoring criteria, distance constraints and selection policies.

        :param oligo_database: The oligo database to process.
        :type oligo_database: OligoDatabase
        :param GC_content_min: Minimum GC content for scoring.
        :type GC_content_min: float
        :param GC_content_opt: Optimal GC content for scoring.
        :type GC_content_opt: float
        :param GC_content_max: Maximum GC content for scoring.
        :type GC_content_max: float
        :param GC_weight: Weight for GC content scoring.
        :type GC_weight: float
        :param Tm_min: Minimum melting temperature for scoring.
        :type Tm_min: float
        :param Tm_opt: Optimal melting temperature for scoring.
        :type Tm_opt: float
        :param Tm_max: Maximum melting temperature for scoring.
        :type Tm_max: float
        :param Tm_weight: Weight for melting temperature scoring.
        :type Tm_weight: float
        :param Tm_parameters: Parameters for melting temperature calculations.
        :type Tm_parameters: dict
        :param Tm_chem_correction_parameters: Parameters for chemical corrections in Tm calculations.
        :type Tm_chem_correction_parameters: dict
        :param Tm_salt_correction_parameters: Parameters for salt corrections in Tm calculations.
        :type Tm_salt_correction_parameters: dict
        :param set_size_opt: Optimal size for oligo sets.
        :type set_size_opt: int
        :param set_size_min: Minimum size for oligo sets.
        :type set_size_min: int
        :param distance_between_oligos: Minimum genomic distance between oligos in a set.
        :type distance_between_oligos: int
        :param n_sets: Number of oligo sets to generate.
        :type n_sets: int
        :param max_graph_size: Maximum size of the graph used in set selection.
        :type max_graph_size: int
        :param n_attempts: Maximum number of attempts for selecting oligo sets.
        :type n_attempts: int
        :param heuristic: Whether to apply heuristic methods in oligo set selection.
        :type heuristic: bool
        :param heuristic_n_attempts: Maximum number of attempts for heuristic selecting oligo sets.
        :type heuristic_n_attempts: int
        :return: The updated oligo database.
        :rtype: OligoDatabase
        """
        set_scoring = AverageSetScoring(ascending=True)

        # We change the processing dependent on the required number of probes in the probe sets
        # For small sets, we don't pre-filter and find the initial set by iterating
        # through all possible generated sets, which is faster than the max clique approximation.
        if set_size_min < 15:
            pre_filter = False
            clique_init_approximation = False
            selection_policy = GraphBasedSelectionPolicy(
                set_scoring=set_scoring,
                pre_filter=pre_filter,
                n_attempts=n_attempts,
                heuristic=heuristic,
                heuristic_n_attempts=heuristic_n_attempts,
                clique_init_approximation=clique_init_approximation,
            )
            base_log_parameters(
                {
                    "pre_filter": pre_filter,
                    "clique_init_approximation": clique_init_approximation,
                    "selection_policy": "Graph-Based",
                }
            )

        # For medium sized sets, we don't pre-filter but we apply the max clique approximation
        # to find an initial probe set faster.
        if set_size_min > 15:
            pre_filter = False
            clique_init_approximation = True
            selection_policy = GraphBasedSelectionPolicy(
                set_scoring=set_scoring,
                pre_filter=pre_filter,
                n_attempts=n_attempts,
                heuristic=heuristic,
                heuristic_n_attempts=heuristic_n_attempts,
                clique_init_approximation=clique_init_approximation,
            )
            base_log_parameters(
                {
                    "pre_filter": pre_filter,
                    "clique_init_approximation": clique_init_approximation,
                    "selection_policy": "Graph-Based",
                }
            )

        # For large sets, we apply the pre-filter which removes all probes from the
        # graph that are only part of cliques which are smaller than the minimum set size
        # and we apply the Greedy Selection Policy istead of the graph-based selection policy.
        if set_size_min > 30:
            pre_filter = True
            selection_policy = GreedySelectionPolicy(
                set_scoring=set_scoring,
                score_criteria=set_scoring.score_1,
                pre_filter=pre_filter,
                penalty=0.01,
                n_attempts=n_attempts,
            )
            base_log_parameters({"pre_filter": pre_filter, "selection_policy": "Greedy"})

        oligos_scoring = WeightedTmGCOligoScoring(
            Tm_min=Tm_min,
            Tm_opt=Tm_opt,
            Tm_max=Tm_max,
            GC_content_min=GC_content_min,
            GC_content_opt=GC_content_opt,
            GC_content_max=GC_content_max,
            Tm_parameters=Tm_parameters,
            Tm_chem_correction_parameters=Tm_chem_correction_parameters,
            Tm_salt_correction_parameters=Tm_salt_correction_parameters,
            Tm_weight=Tm_weight,
            GC_weight=GC_weight,
        )
        oligoset_generator = OligosetGeneratorIndependentSet(
            selection_policy=selection_policy,
            oligos_scoring=oligos_scoring,
            set_scoring=set_scoring,
            max_oligos=max_graph_size,
            distance_between_oligos=distance_between_oligos,
        )
        oligo_database = oligoset_generator.apply(
            oligo_database=oligo_database,
            sequence_type="oligo",
            set_size_opt=set_size_opt,
            set_size_min=set_size_min,
            n_sets=n_sets,
            n_jobs=self.n_jobs,
        )

        return oligo_database


############################################
# Oligo-seq Designer Pipeline
############################################


def main():
    """
    Main function for running the OligoSeqProbeDesigner pipeline. This function reads the configuration file,
    processes gene IDs, initializes the probe designer, sets developer parameters, and executes probe design
    and output generation steps.

    :param args: Command-line arguments parsed using the base parser. The arguments include:
        - config: Path to the configuration YAML file containing parameters for the pipeline.
    :type args: dict
    """
    print("--------------START PIPELINE--------------")

    args = base_parser()

    ##### read the config file #####
    with open(args["config"], "r") as handle:
        config = yaml.safe_load(handle)

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
        write_intermediate_steps=config["write_intermediate_steps"],
        dir_output=config["dir_output"],
        n_jobs=config["n_jobs"],
    )

    ##### set custom developer parameters #####
    pipeline.set_developer_parameters(
        target_probe_hybridization_probability_alignment_method=config[
            "target_probe_hybridization_probability_alignment_method"
        ],
        target_probe_hybridization_probability_blastn_search_parameters=config[
            "target_probe_hybridization_probability_blastn_search_parameters"
        ],
        target_probe_hybridization_probability_blastn_hit_parameters=config[
            "target_probe_hybridization_probability_blastn_hit_parameters"
        ],
        target_probe_hybridization_probability_bowtie_search_parameters=config[
            "target_probe_hybridization_probability_bowtie_search_parameters"
        ],
        target_probe_hybridization_probability_bowtie_hit_parameters=config[
            "target_probe_hybridization_probability_bowtie_hit_parameters"
        ],
        target_probe_cross_hybridization_alignment_method=config[
            "target_probe_cross_hybridization_alignment_method"
        ],
        target_probe_cross_hybridization_blastn_search_parameters=config[
            "target_probe_cross_hybridization_blastn_search_parameters"
        ],
        target_probe_cross_hybridization_blastn_hit_parameters=config[
            "target_probe_cross_hybridization_blastn_hit_parameters"
        ],
        target_probe_cross_hybridization_bowtie_search_parameters=config[
            "target_probe_cross_hybridization_bowtie_search_parameters"
        ],
        target_probe_cross_hybridization_bowtie_hit_parameters=config[
            "target_probe_cross_hybridization_bowtie_hit_parameters"
        ],
        max_graph_size=config["max_graph_size"],
        n_attempts=config["n_attempts"],
        heuristic=config["heuristic"],
        heuristic_n_attempts=config["heuristic_n_attempts"],
        target_probe_Tm_parameters=config["target_probe_Tm_parameters"],
        target_probe_Tm_chem_correction_parameters=config["target_probe_Tm_chem_correction_parameters"],
        target_probe_Tm_salt_correction_parameters=config["target_probe_Tm_salt_correction_parameters"],
    )

    ##### design probes #####
    oligo_database = pipeline.design_target_probes(
        files_fasta_target_probe_database=config["files_fasta_target_probe_database"],
        files_fasta_reference_database_targe_probe=config["files_fasta_reference_database_targe_probe"],
        gene_ids=gene_ids,
        target_probe_length_min=config["target_probe_length_min"],
        target_probe_length_max=config["target_probe_length_max"],
        target_probe_split_region=config["target_probe_split_region"],
        target_probe_targeted_exons=config["target_probe_targeted_exons"],
        target_probe_isoform_consensus=config["target_probe_isoform_consensus"],
        target_probe_GC_content_min=config["target_probe_GC_content_min"],
        target_probe_GC_content_opt=config["target_probe_GC_content_opt"],
        target_probe_GC_content_max=config["target_probe_GC_content_max"],
        target_probe_Tm_min=config["target_probe_Tm_min"],
        target_probe_Tm_opt=config["target_probe_Tm_opt"],
        target_probe_Tm_max=config["target_probe_Tm_max"],
        target_probe_secondary_structures_T=config["target_probe_secondary_structures_T"],
        target_probe_secondary_structures_threshold_deltaG=config[
            "target_probe_secondary_structures_threshold_deltaG"
        ],
        target_probe_homopolymeric_base_n=config["target_probe_homopolymeric_base_n"],
        target_probe_max_len_selfcomplement=config["target_probe_max_len_selfcomplement"],
        target_probe_hybridization_probability_threshold=config[
            "target_probe_hybridization_probability_threshold"
        ],
        target_probe_read_length_bias=config["target_probe_read_length_bias"],
        target_probe_GC_weight=config["target_probe_GC_weight"],
        target_probe_Tm_weight=config["target_probe_Tm_weight"],
        set_size_min=config["set_size_min"],
        set_size_opt=config["set_size_opt"],
        distance_between_target_probes=config["distance_between_target_probes"],
        n_sets=config["n_sets"],
    )

    pipeline.generate_output(oligo_database=oligo_database, top_n_sets=config["top_n_sets"])

    print("--------------END PIPELINE--------------")


if __name__ == "__main__":
    main()
