############################################
# imports
############################################

import os
import sys
import yaml
import shutil
import logging
import warnings

import numpy as np
import pandas as pd

from typing import List, Tuple
from pathlib import Path
from datetime import datetime
from itertools import product

from Bio.SeqUtils import Seq
from Bio.SeqUtils import MeltingTemp as mt

from oligo_designer_toolsuite.database import (
    OligoAttributes,
    OligoDatabase,
    ReferenceDatabase,
)
from oligo_designer_toolsuite.oligo_efficiency_filter import (
    LowestSetScoring,
    WeightedGCUtrScoring,
)
from oligo_designer_toolsuite.oligo_property_filter import (
    SoftMaskedSequenceFilter,
    HardMaskedSequenceFilter,
    HomopolymericRunsFilter,
    GCContentFilter,
    GCClampFilter,
    MeltingTemperatureNNFilter,
    SelfComplementFilter,
    ComplementFilter,
    SecondaryStructureFilter,
    PropertyFilter,
)
from oligo_designer_toolsuite.oligo_selection import (
    OligosetGeneratorIndependentSet,
    GraphBasedSelectionPolicy,
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
    pipeline_step_basic,
    check_content_oligo_database,
)
from oligo_designer_toolsuite.utils import append_nucleotide_to_sequences
from oligo_designer_toolsuite.sequence_generator import OligoSequenceGenerator

############################################
# SeqFISH Plus Probe Designer
############################################


class SeqFishPlusProbeDesigner:
    """
    Class for designing probes for SeqFISH Plus experiments.

    SeqFISH Plus probes are used to encode spatial transcriptomic information. This designer handles the design 
    of target probes, readout probes, and their integration into encoding probes with optional primer sequences.
    """
    def __init__(self, write_intermediate_steps: bool, dir_output: str, n_jobs: int) -> None:
        
        """
        Constructor for the SeqFishPlusProbeDesigner class.

        :param write_intermediate_steps: Whether to save intermediate results during pipeline execution.
        :type write_intermediate_steps: bool
        :param dir_output: Directory path for output files.
        :type dir_output: str
        :param n_jobs: Number of parallel jobs to use for computations.
        :type n_jobs: int
        """

        ##### create the output folder #####
        self.dir_output = os.path.abspath(dir_output)
        Path(self.dir_output).mkdir(parents=True, exist_ok=True)

        ##### setup logger #####
        timestamp = datetime.now()
        file_logger = os.path.join(
            self.dir_output,
            f"log_seqfishplus_probe_designer_{timestamp.year}-{timestamp.month}-{timestamp.day}-{timestamp.hour}-{timestamp.minute}.txt",
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
        target_probe_specificity_blastn_search_parameters: dict = {
            "perc_identity": 100,
            "strand": "minus",
            "word_size": 7,
            "dust": "no",
            "soft_masking": "false",
            "max_target_seqs": 10,
            "max_hsps": 1000,
        },
        target_probe_specificity_blastn_hit_parameters: dict = {"min_alignment_length": 15},
        target_probe_cross_hybridization_blastn_search_parameters: dict = {
            "perc_identity": 80,
            "strand": "minus",
            "word_size": 7,
            "dust": "no",
            "soft_masking": "false",
            "max_target_seqs": 10,
        },
        target_probe_cross_hybridization_blastn_hit_parameters: dict = {"min_alignment_length": 17},
        readout_probe_initial_num_sequences: int = 100000,
        readout_probe_specificity_blastn_search_parameters: dict = {
            "perc_identity": 100,
            "strand": "minus",
            "word_size": 7,
            "dust": "no",
            "soft_masking": "false",
            "max_target_seqs": 10,
            "max_hsps": 1000,
        },
        readout_probe_specificity_blastn_hit_parameters: dict = {"min_alignment_length": 10},
        readout_probe_cross_hybridization_blastn_search_parameters: dict = {
            "perc_identity": 100,
            "strand": "minus",
            "word_size": 7,
            "dust": "no",
            "soft_masking": "false",
            "max_target_seqs": 10,
        },
        readout_probe_cross_hybridization_blastn_hit_parameters: dict = {"min_alignment_length": 10},
        primer_initial_num_sequences: int = 1000000,
        primer_specificity_refrence_blastn_search_parameters: dict = {
            "perc_identity": 100,
            "strand": "minus",
            "word_size": 7,
            "dust": "no",
            "soft_masking": "false",
            "max_target_seqs": 10,
            "max_hsps": 1000,
        },
        primer_specificity_refrence_blastn_hit_parameters: dict = {"min_alignment_length": 14},
        primer_specificity_encoding_probes_blastn_search_parameters: dict = {
            "perc_identity": 100,
            "strand": "minus",
            "word_size": 7,
            "dust": "no",
            "soft_masking": "false",
            "max_target_seqs": 10,
            "max_hsps": 1000,
        },
        primer_specificity_encoding_probes_blastn_hit_parameters: dict = {"min_alignment_length": 11},
        primer_Tm_parameters: dict = {
            "check": True,
            "strict": True,
            "c_seq": None,
            "shift": 0,
            "nn_table": "DNA_NN4",
            "tmm_table": "DNA_TMM1",
            "imm_table": "DNA_IMM1",
            "de_table": "DNA_DE1",
            "dnac1": 500,
            "dnac2": 25,
            "selfcomp": False,
            "saltcorr": 5,
            "Na": 300,
            "K": 0,
            "Tris": 0,
            "Mg": 0,
            "dNTPs": 0,
        },
        primer_Tm_chem_correction_parameters: dict = None,
        primer_Tm_salt_correction_parameters: dict = None,
        max_graph_size: int = 5000,
        pre_filter: bool = False,
        n_attempts: int = 100000,
        heuristic: bool = True,
        heuristic_n_attempts: int = 100,
    ):
        """
        Configure default parameters for the probe designer.

        This method sets default parameters for various stages of probe design including specificity filters, 
        melting temperature calculations, and graph-based oligo selection.

        :param target_probe_specificity_blastn_search_parameters: Parameters for BlastN specificity search for target probes.
        :type target_probe_specificity_blastn_search_parameters: dict, optional
        :param target_probe_specificity_blastn_hit_parameters: Parameters for BlastN hit filtering for target probes.
        :type target_probe_specificity_blastn_hit_parameters: dict, optional
        :param target_probe_cross_hybridization_blastn_search_parameters: Parameters for BlastN cross-hybridization search.
        :type target_probe_cross_hybridization_blastn_search_parameters: dict, optional
        :param target_probe_cross_hybridization_blastn_hit_parameters: Parameters for filtering BlastN hits in cross-hybridization.
        :type target_probe_cross_hybridization_blastn_hit_parameters: dict, optional
        :param readout_probe_initial_num_sequences: Initial number of sequences for readout probe generation.
        :type readout_probe_initial_num_sequences: int, optional
        :param readout_probe_specificity_blastn_search_parameters: Parameters for BlastN specificity search for readout probes.
        :type readout_probe_specificity_blastn_search_parameters: dict, optional
        :param readout_probe_specificity_blastn_hit_parameters: Parameters for filtering BlastN hits for readout probes.
        :type readout_probe_specificity_blastn_hit_parameters: dict, optional
        :param primer_initial_num_sequences: Initial number of sequences for primer design.
        :type primer_initial_num_sequences: int, optional
        :param primer_specificity_refrence_blastn_search_parameters: Parameters for BlastN specificity search for primers.
        :type primer_specificity_refrence_blastn_search_parameters: dict, optional
        :param primer_specificity_refrence_blastn_hit_parameters: Parameters for filtering BlastN hits for primers.
        :type primer_specificity_refrence_blastn_hit_parameters: dict, optional
        :param primer_specificity_encoding_probes_blastn_search_parameters: Parameters for BlastN specificity search for primers.
        :type primer_specificity_encoding_probes_blastn_search_parameters: dict, optional
        :param primer_specificity_encoding_probes_blastn_hit_parameters: Parameters for filtering BlastN hits for primers.
        :type primer_specificity_encoding_probes_blastn_hit_parameters: dict, optional
        :param primer_Tm_parameters: Parameters for melting temperature calculations for primers.
        :type primer_Tm_parameters: dict, optional
        :param primer_Tm_chem_correction_parameters: Parameters for chemical corrections in melting temperature calculations.
        :type primer_Tm_chem_correction_parameters: dict, optional
        :param primer_Tm_salt_correction_parameters: Parameters for salt corrections in melting temperature calculations.
        :type primer_Tm_salt_correction_parameters: dict, optional
        :param max_graph_size: Maximum size of the graph for oligo set selection.
        :type max_graph_size: int, optional
        :param pre_filter: Whether to pre-filter sequences before oligo selection.
        :type pre_filter: bool, optional
        :param n_attempts: Number of attempts for oligo selection.
        :type n_attempts: int, optional
        :param heuristic: Whether to use heuristic methods for oligo selection.
        :type heuristic: bool, optional
        :param heuristic_n_attempts: Number of heuristic attempts for oligo selection.
        :type heuristic_n_attempts: int, optional
        """
        ### Parameters for the specificity filters
        # Specificity filter with BlastN
        self.target_probe_specificity_blastn_search_parameters = (
            target_probe_specificity_blastn_search_parameters
        )
        self.target_probe_specificity_blastn_hit_parameters = target_probe_specificity_blastn_hit_parameters

        # Crosshybridization filter with BlastN
        self.target_probe_cross_hybridization_blastn_search_parameters = (
            target_probe_cross_hybridization_blastn_search_parameters
        )
        self.target_probe_cross_hybridization_blastn_hit_parameters = (
            target_probe_cross_hybridization_blastn_hit_parameters
        )

        ## readout probe sequence
        self.readout_probe_initial_num_sequences = readout_probe_initial_num_sequences
        # Specificity filter with BlastN
        self.readout_probe_specificity_blastn_search_parameters = (
            readout_probe_specificity_blastn_search_parameters
        )
        self.readout_probe_specificity_blastn_hit_parameters = readout_probe_specificity_blastn_hit_parameters

        # Crosshybridization filter with BlastN
        self.readout_probe_cross_hybridization_blastn_search_parameters = (
            readout_probe_cross_hybridization_blastn_search_parameters
        )
        self.readout_probe_cross_hybridization_blastn_hit_parameters = (
            readout_probe_cross_hybridization_blastn_hit_parameters
        )

        ##forward primer sequence
        self.primer_initial_num_sequences = primer_initial_num_sequences
        # Specificity filter with BlastN against reference
        self.primer_specificity_refrence_blastn_search_parameters = (
            primer_specificity_refrence_blastn_search_parameters
        )
        self.primer_specificity_refrence_blastn_hit_parameters = (
            primer_specificity_refrence_blastn_hit_parameters
        )

        # Specificity filter with BlastN against encoding rpobes
        self.primer_specificity_encoding_probes_blastn_search_parameters = (
            primer_specificity_encoding_probes_blastn_search_parameters
        )
        self.primer_specificity_encoding_probes_blastn_hit_parameters = (
            primer_specificity_encoding_probes_blastn_hit_parameters
        )

        ### Parameters for Melting Temperature
        # The melting temperature is used in 2 different stages (property filters and padlock detection probe design), where a few parameters are shared and the others differ.
        # parameters for melting temperature -> for more information on parameters, see: https://biopython.org/docs/1.75/api/Bio.SeqUtils.MeltingTemp.html#Bio.SeqUtils.MeltingTemp.Tm_NN

        # preprocess melting temperature parameters
        primer_Tm_parameters["nn_table"] = getattr(mt, primer_Tm_parameters["nn_table"])
        primer_Tm_parameters["tmm_table"] = getattr(mt, primer_Tm_parameters["tmm_table"])
        primer_Tm_parameters["imm_table"] = getattr(mt, primer_Tm_parameters["imm_table"])
        primer_Tm_parameters["de_table"] = getattr(mt, primer_Tm_parameters["de_table"])

        ## primer
        self.primer_Tm_parameters = primer_Tm_parameters
        self.primer_Tm_chem_correction_parameters = primer_Tm_chem_correction_parameters
        self.primer_Tm_salt_correction_parameters = primer_Tm_salt_correction_parameters

        ### Parameters for the Oligo set selection
        self.max_graph_size = max_graph_size
        self.pre_filter = pre_filter
        self.heuristic = heuristic
        self.n_attempts = n_attempts
        self.heuristic_n_attempts = heuristic_n_attempts

    def design_target_probes(
        self,
        files_fasta_target_probe_database: list[str],
        files_fasta_reference_database_targe_probe: List[str],
        gene_ids: list = None,
        target_probe_length_min: int = 28,
        target_probe_length_max: int = 28,
        target_probe_isoform_consensus: float = 100,
        target_probe_GC_content_min: float = 45,
        target_probe_GC_content_opt: float = 55,
        target_probe_GC_content_max: float = 65,
        target_probe_homopolymeric_base_n: dict = {"A": 5, "T": 5, "C": 5, "G": 5},
        target_probe_T_secondary_structure: float = 76,
        target_probe_secondary_structures_threshold_deltaG: float = 0,
        target_probe_GC_weight: float = 1,
        target_probe_UTR_weight: float = 10,
        set_size_opt: int = 24,
        set_size_min: int = 24,
        distance_between_target_probes: int = 2,
        n_sets: int = 100,
    ) -> OligoDatabase:
        """
        Design target probes and run the SeqFISHPlus target probe designer.

        This method creates, filters, and optimizes a database of target probes based on various 
        design criteria: property, specificity and does oligo selection.

        :param files_fasta_target_probe_database: List of FASTA files containing target probe sequences.
        :type files_fasta_target_probe_database: list
        :param files_fasta_reference_database_targe_probe: List of FASTA files for reference database.
        :type files_fasta_reference_database_targe_probe: List[str]
        :param gene_ids: List of gene IDs to target. Defaults to None (target all genes).
        :type gene_ids: list, optional
        :param target_probe_length_min: Minimum probe length.
        :type target_probe_length_min: int
        :param target_probe_isoform_consensus: Minimum isoform consensus percentage.
        :type target_probe_isoform_consensus: float
        :param target_probe_GC_content_min: Minimum GC content percentage.
        :type target_probe_GC_content_min: float
        :param target_probe_GC_content_opt: Optimal GC content percentage.
        :type target_probe_GC_content_opt: float
        :param target_probe_GC_content_max: Maximum GC content percentage.
        :type target_probe_GC_content_max: float
        :param target_probe_homopolymeric_base_n: Maximum allowed homopolymeric bases.
        :type target_probe_homopolymeric_base_n: dict
        :param target_probe_T_secondary_structure: Threshold for secondary structure temperature.
        :type target_probe_T_secondary_structure: float
        :param target_probe_secondary_structures_threshold_deltaG: DeltaG threshold for secondary structures.
        :type target_probe_secondary_structures_threshold_deltaG: float
        :param target_probe_GC_weight: Weight for GC content in scoring.
        :type target_probe_GC_weight: float
        :param target_probe_UTR_weight: Weight for UTR regions in scoring.
        :type target_probe_UTR_weight: float
        :param set_size_opt: Optimal size for probe sets.
        :type set_size_opt: int
        :param set_size_min: Minimum size for probe sets.
        :type set_size_min: int
        :param distance_between_target_probes: Minimum distance between probes in a set.
        :type distance_between_target_probes: int
        :param n_sets: Number of probe sets to generate.
        :type n_sets: int
        :return: A database of designed target probes.
        :rtype: OligoDatabase
        """

        target_probe_designer = SeqFishTargetProbeDesigner(self.dir_output, self.n_jobs)

        oligo_database = target_probe_designer.create_oligo_database(
            gene_ids=gene_ids,
            oligo_length_min=target_probe_length_min,
            oligo_length_max=target_probe_length_max,
            files_fasta_oligo_database=files_fasta_target_probe_database,
            min_oligos_per_gene=set_size_min,
            isoform_consensus=target_probe_isoform_consensus,
        )
        check_content_oligo_database(oligo_database)

        if self.write_intermediate_steps:
            dir_database = oligo_database.save_database(dir_database="1_db_target_probes_initial")
            print(f"Saved target probe database for step 1 (Create Database) in directory {dir_database}")

        oligo_database = target_probe_designer.filter_by_property(
            oligo_database=oligo_database,
            GC_content_min=target_probe_GC_content_min,
            GC_content_max=target_probe_GC_content_max,
            homopolymeric_base_n=target_probe_homopolymeric_base_n,
            T_secondary_structure=target_probe_T_secondary_structure,
            secondary_structures_threshold_deltaG=target_probe_secondary_structures_threshold_deltaG,
        )
        check_content_oligo_database(oligo_database)

        if self.write_intermediate_steps:
            dir_database = oligo_database.save_database(dir_database="2_db_target_probes_property_filter")
            print(f"Saved target probe database for step 2 (Property Filters) in directory {dir_database}")

        oligo_database = target_probe_designer.filter_by_specificity(
            oligo_database=oligo_database,
            files_fasta_reference_database=files_fasta_reference_database_targe_probe,
            specificity_blastn_search_parameters=self.target_probe_specificity_blastn_search_parameters,
            specificity_blastn_hit_parameters=self.target_probe_specificity_blastn_hit_parameters,
            cross_hybridization_blastn_search_parameters=self.target_probe_cross_hybridization_blastn_search_parameters,
            cross_hybridization_blastn_hit_parameters=self.target_probe_cross_hybridization_blastn_hit_parameters,
        )
        check_content_oligo_database(oligo_database)

        if self.write_intermediate_steps:
            dir_database = oligo_database.save_database(dir_database="3_db_target_probes_specificity_filter")
            print(f"Saved target probe database for step 3 (Specificity Filters) in directory {dir_database}")

        oligo_database = target_probe_designer.create_oligo_sets(
            oligo_database=oligo_database,
            GC_weight=target_probe_GC_weight,
            GC_content_opt=target_probe_GC_content_opt,
            UTR_weight=target_probe_UTR_weight,
            set_size_opt=set_size_opt,
            set_size_min=set_size_min,
            max_graph_size=self.max_graph_size,
            n_sets=n_sets,
            n_attempts=self.n_attempts,
            pre_filter=self.pre_filter,
            heuristic=self.heuristic,
            heuristic_n_attempts=self.heuristic_n_attempts,
            distance_between_oligos=distance_between_target_probes,
        )
        check_content_oligo_database(oligo_database)

        if self.write_intermediate_steps:
            dir_database = oligo_database.save_database(dir_database="4_db_target_probes_sets")
            dir_oligosets = oligo_database.write_oligosets_to_table()
            print(
                f"Saved target probe database for step 4 (Specificity Filters) in directory {dir_database} and sets table in directory {dir_oligosets}"
            )

        return oligo_database

    def design_readout_probes(
        self,
        n_genes: int,
        files_fasta_reference_database_readout_probe: List[str],
        readout_probe_length: int = 15,
        readout_probe_base_probabilities: dict = {"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25},
        readout_probe_GC_content_min: float = 40,
        readout_probe_GC_content_max: float = 60,
        readout_probe_homopolymeric_base_n: dict = {"G": 3},
        n_barcode_rounds: int = 4,
        n_pseudocolors: int = 20,
        channels_ids: list = ["Alexa488", "Cy3b", "Alexa647"],
    ) -> Tuple[pd.DataFrame, pd.DataFrame]:
        """
        Design readout probes, and run the SeqFISHPlus readout probe designer pipeline.

        This method creates, filters, and optimizes a database of readout probes based on various
        design criteria: property, specificity and does oligo selection.

        :param n_genes: Number of genes to target.
        :type n_genes: int
        :param files_fasta_reference_database_readout_probe: List of FASTA files for reference database.
        :type files_fasta_reference_database_readout_probe: List[str]
        :param readout_probe_length: Length of readout probes.
        :type readout_probe_length: int
        :param readout_probe_base_probabilities: Base probabilities for readout probes.
        :type readout_probe_base_probabilities: dict
        :param readout_probe_GC_content_min: Minimum GC content percentage.
        :type readout_probe_GC_content_min: float
        :param readout_probe_GC_content_max: Maximum GC content percentage.
        :type readout_probe_GC_content_max: float
        :param readout_probe_homopolymeric_base_n: Maximum allowed homopolymeric bases.
        :type readout_probe_homopolymeric_base_n: dict
        :param n_barcode_rounds: Number of barcode rounds.
        :type n_barcode_rounds: int
        :param n_pseudocolors: Number of pseudocolors.
        :type n_pseudocolors: int
        :param channels_ids: List of channel IDs.
        :type channels_ids: list
        :return: A codebook and a table of readout probes.
        :rtype: Tuple[pd.DataFrame, pd.DataFrame]
        """



        readout_probe_designer = SeqFishPlusReadoutProbeDesigner(
            dir_output=self.dir_output,
            n_jobs=self.n_jobs,
        )
        oligo_database = readout_probe_designer.create_oligo_database(
            oligo_length=readout_probe_length,
            oligo_base_probabilities=readout_probe_base_probabilities,
            initial_num_sequences=self.readout_probe_initial_num_sequences,
        )
        check_content_oligo_database(oligo_database)

        if self.write_intermediate_steps:
            dir_database = oligo_database.save_database(dir_database="1_db_readout_probes_initial")
            print(f"Saved readout probe database for step 1 (Create Database) in directory {dir_database}")

        oligo_database = readout_probe_designer.filter_by_property(
            oligo_database=oligo_database,
            GC_content_min=readout_probe_GC_content_min,
            GC_content_max=readout_probe_GC_content_max,
            homopolymeric_base_n=readout_probe_homopolymeric_base_n,
        )
        check_content_oligo_database(oligo_database)

        if self.write_intermediate_steps:
            dir_database = oligo_database.save_database(dir_database="2_db_readout_probes_property_filter")
            print(f"Saved readout probe database for step 2 (Property Filters) in directory {dir_database}")

        oligo_database = readout_probe_designer.filter_by_specificity(
            oligo_database=oligo_database,
            files_fasta_reference_database=files_fasta_reference_database_readout_probe,
            specificity_blastn_search_parameters=self.readout_probe_specificity_blastn_search_parameters,
            specificity_blastn_hit_parameters=self.readout_probe_specificity_blastn_hit_parameters,
            cross_hybridization_blastn_search_parameters=self.readout_probe_cross_hybridization_blastn_search_parameters,
            cross_hybridization_blastn_hit_parameters=self.readout_probe_cross_hybridization_blastn_hit_parameters,
        )
        check_content_oligo_database(oligo_database)

        if self.write_intermediate_steps:
            dir_database = oligo_database.save_database(dir_database="3_db_readout_probes_specificty_filter")
            print(
                f"Saved readout probe database for step 3 (Specificity Filters) in directory {dir_database}"
            )

        codebook = readout_probe_designer.generate_codebook(
            n_regions=n_genes,
            n_barcode_rounds=n_barcode_rounds,
            n_pseudocolors=n_pseudocolors,
            n_channels=len(channels_ids),
        )

        readout_probe_table = readout_probe_designer.create_readout_probe_table(
            readout_probe_database=oligo_database,
            channels_ids=channels_ids,
            n_barcode_rounds=n_barcode_rounds,
            n_pseudocolors=n_pseudocolors,
        )

        return codebook, readout_probe_table

    def design_encoding_probe(
        self,
        target_probe_database: OligoDatabase,
        codebook: pd.DataFrame,
        readout_probe_table: pd.DataFrame,
    ):
        """
        Design encoding probes by combining target probes with readout probe sequences based on the codebook.

        This method generates encoding probes for each region in the target probe database by combining
        sequences from the readout probes and their corresponding barcodes defined in the codebook.

        :param target_probe_database: Database of target probes containing sequence and attribute information.
        :type target_probe_database: OligoDatabase
        :param codebook: A DataFrame containing barcodes for each region. Each row corresponds to a region,
            with columns representing bits in the barcode.
        :type codebook: pd.DataFrame
        :param readout_probe_table: A DataFrame containing readout probe sequences and their associated bit
            identifiers.
        :type readout_probe_table: pd.DataFrame

        :return: Updated `target_probe_database` with attributes for encoding probes, including sequences
            for target probes, readout probes, and the full encoding probe sequence.
        :rtype: OligoDatabase
        """

        region_ids = list(target_probe_database.database.keys())

        codebook.index = region_ids + [
            f"unassigned_barcode_{i+1}" for i in range(len(codebook.index) - len(region_ids))
        ]
        # codebook = codebook.iloc[: len(region_ids)]

        codebook.to_csv(os.path.join(self.dir_output, "codebook.tsv"), sep="\t")
        readout_probe_table.to_csv(os.path.join(self.dir_output, "readout_probes.tsv"), sep="\t")

        for region_id in region_ids:
            barcode = codebook.loc[region_id]

            bits = barcode[barcode == 1].index
            readout_probe_sequences = readout_probe_table.loc[bits, "readout_probe_sequence"]

            new_probe_attributes_encoding_probe = {}

            for probe_id in target_probe_database.database[region_id].keys():

                sequence_readout_probe_1 = readout_probe_sequences[0]
                sequence_readout_probe_2 = readout_probe_sequences[1]
                sequence_readout_probe_3 = readout_probe_sequences[2]
                sequence_readout_probe_4 = readout_probe_sequences[3]

                new_probe_attributes_encoding_probe[probe_id] = {
                    "barcode": barcode,
                    "sequence_target": target_probe_database.get_oligo_attribute_value(
                        attribute="target", region_id=region_id, oligo_id=probe_id, flatten=True
                    ),
                    "sequence_target_probe": target_probe_database.get_oligo_attribute_value(
                        attribute="oligo", region_id=region_id, oligo_id=probe_id, flatten=True
                    ),
                    "sequence_readout_probe_1": sequence_readout_probe_1,
                    "sequence_readout_probe_2": sequence_readout_probe_2,
                    "sequence_readout_probe_3": sequence_readout_probe_3,
                    "sequence_readout_probe_4": sequence_readout_probe_4,
                    "sequence_encoding_probe": (
                        str(Seq(sequence_readout_probe_1).reverse_complement())
                        + str(Seq(sequence_readout_probe_2).reverse_complement())
                        + "T"
                        + target_probe_database.get_oligo_attribute_value(
                            attribute="oligo", region_id=region_id, oligo_id=probe_id, flatten=True
                        )
                        + "T"
                        + str(Seq(sequence_readout_probe_3).reverse_complement())
                        + str(Seq(sequence_readout_probe_4).reverse_complement())
                    ),
                }

            target_probe_database.update_oligo_attributes(new_probe_attributes_encoding_probe)

        return target_probe_database

    def design_primers(
        self,
        encoding_probe_database: str,
        files_fasta_reference_database_primer: List[str],
        reverse_primer_sequence: str = "CCCTATAGTGAGTCGTATTA",
        primer_length: int = 20,
        primer_base_probabilities: dict = {"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25},
        primer_GC_content_min: float = 50,
        primer_GC_content_max: float = 65,
        primer_number_GC_GCclamp: int = 1,
        primer_number_three_prime_base_GCclamp: int = 2,
        primer_homopolymeric_base_n: int = {"A": 4, "T": 4, "C": 4, "G": 4},
        primer_max_len_selfcomplement: int = 6,
        primer_max_len_complement_reverse_primer: int = 5,
        primer_Tm_min: float = 60,
        primer_Tm_max: float = 75,
        primer_T_secondary_structure: float = 76,
        primer_secondary_structures_threshold_deltaG: float = 0,
    ):
        """
        Design forward and reverse primers for the encoding probe database.

        This method runs the SeqFishPlus primer designer pipeline to create, filter, and optimize a database. Then the
        best forward primer is selected based on the melting temperature (Tm) of the reverse primer.

        :param encoding_probe_database: The encoding probe database containing sequences and regions.
        :type encoding_probe_database: str
        :param files_fasta_reference_database_primer: List of FASTA files containing reference sequences 
            for primer specificity filtering.
        :type files_fasta_reference_database_primer: list[str]
        :param reverse_primer_sequence: Sequence of the reverse primer.
        :type reverse_primer_sequence: str, optional
        :param primer_length: Length of the forward primers to design.
        :type primer_length: int, optional
        :param primer_base_probabilities: Dictionary specifying base probabilities for random primer 
            sequence generation.
        :type primer_base_probabilities: dict, optional
        :param primer_GC_content_min: Minimum acceptable GC content for primers.
        :type primer_GC_content_min: float, optional
        :param primer_GC_content_max: Maximum acceptable GC content for primers.
        :type primer_GC_content_max: float, optional
        :param primer_number_GC_GCclamp: Minimum number of GC bases required at the primer's 3' end.
        :type primer_number_GC_GCclamp: int, optional
        :param primer_number_three_prime_base_GCclamp: Minimum number of GC bases required within 
            the last three bases at the primer's 3' end.
        :type primer_number_three_prime_base_GCclamp: int, optional
        :param primer_homopolymeric_base_n: Dictionary specifying maximum allowed length of homopolymeric 
            stretches.
        :type primer_homopolymeric_base_n: dict, optional
        :param primer_max_len_selfcomplement: Maximum length of self-complementary sequences in the primer.
        :type primer_max_len_selfcomplement: int, optional
        :param primer_max_len_complement_reverse_primer: Maximum length of complementary sequences between 
            the forward and reverse primers.
        :type primer_max_len_complement_reverse_primer: int, optional
        :param primer_Tm_min: Minimum acceptable melting temperature (Tm) for primers.
        :type primer_Tm_min: float, optional
        :param primer_Tm_max: Maximum acceptable melting temperature (Tm) for primers.
        :type primer_Tm_max: float, optional
        :param primer_T_secondary_structure: Maximum allowable melting temperature for secondary structures.
        :type primer_T_secondary_structure: float, optional
        :param primer_secondary_structures_threshold_deltaG: Threshold for the free energy (deltaG) 
            of secondary structures.
        :type primer_secondary_structures_threshold_deltaG: float, optional

        :return: The selected reverse primer sequence and the best forward primer sequence.
        :rtype: Tuple[str, str]
        """

        file_fasta_encoding_probes_database = encoding_probe_database.write_database_to_fasta(
            filename=f"db_reference_encoding_probes",
            save_description=False,
            region_ids=None,
            sequence_type="sequence_encoding_probe",
        )

        primer_designer = SeqFishPlusPrimerDesigner(
            dir_output=self.dir_output,
            n_jobs=self.n_jobs,
        )
        oligo_database = primer_designer.create_oligo_database(
            oligo_length=primer_length,
            oligo_base_probabilities=primer_base_probabilities,
            initial_num_sequences=self.primer_initial_num_sequences,
        )
        check_content_oligo_database(oligo_database)

        if self.write_intermediate_steps:
            dir_database = oligo_database.save_database(dir_database="1_db_primers_initial")
            print(f"Saved primer database for step 1 (Create Database) in directory {dir_database}")

        oligo_database = primer_designer.filter_by_property(
            oligo_database=oligo_database,
            GC_content_min=primer_GC_content_min,
            GC_content_max=primer_GC_content_max,
            number_GC_GCclamp=primer_number_GC_GCclamp,
            number_three_prime_base_GCclamp=primer_number_three_prime_base_GCclamp,
            homopolymeric_base_n=primer_homopolymeric_base_n,
            max_len_selfcomplement=primer_max_len_selfcomplement,
            reverse_primer_sequence=reverse_primer_sequence,
            max_len_complement=primer_max_len_complement_reverse_primer,
            Tm_min=primer_Tm_min,
            Tm_max=primer_Tm_max,
            Tm_parameters=self.primer_Tm_parameters,
            Tm_chem_correction_parameters=self.primer_Tm_chem_correction_parameters,
            Tm_salt_correction_parameters=self.primer_Tm_salt_correction_parameters,
            T_secondary_structure=primer_T_secondary_structure,
            secondary_structures_threshold_deltaG=primer_secondary_structures_threshold_deltaG,
        )
        check_content_oligo_database(oligo_database)

        if self.write_intermediate_steps:
            dir_database = oligo_database.save_database(dir_database="2_db_primer_property_filter")
            print(f"Saved primer database for step 2 (Property Filters) in directory {dir_database}")

        oligo_database = primer_designer.filter_by_specificity(
            oligo_database=oligo_database,
            files_fasta_reference_database=files_fasta_reference_database_primer,
            specificity_refrence_blastn_search_parameters=self.primer_specificity_refrence_blastn_search_parameters,
            specificity_refrence_blastn_hit_parameters=self.primer_specificity_refrence_blastn_hit_parameters,
            file_fasta_encoding_probes_database=file_fasta_encoding_probes_database,
            specificity_encoding_probes_blastn_search_parameters=self.primer_specificity_encoding_probes_blastn_search_parameters,
            specificity_encoding_probes_blastn_hit_parameters=self.primer_specificity_encoding_probes_blastn_hit_parameters,
        )
        check_content_oligo_database(oligo_database)

        if self.write_intermediate_steps:
            dir_database = oligo_database.save_database(dir_database="3_db_primer_specificty_filter")
            print(f"Saved primer database for step 3 (Specificity Filters) in directory {dir_database}")

        min_dif_Tm = 100

        # calculate Tm for the reverse primer
        Tm_reverse_primer = OligoAttributes._calc_TmNN(
            sequence=reverse_primer_sequence,
            Tm_parameters=self.primer_Tm_parameters,
            Tm_chem_correction_parameters=self.primer_Tm_chem_correction_parameters,
            Tm_salt_correction_parameters=self.primer_Tm_salt_correction_parameters,
        )

        # iterate over all primers in the database to find the one with Tm closest to the reverse primer Tm
        for database_region in oligo_database.database.values():
            for primer_attributes in database_region.values():
                Tm_forward_primer = OligoAttributes._calc_TmNN(
                    sequence=primer_attributes["oligo"],
                    Tm_parameters=self.primer_Tm_parameters,
                    Tm_chem_correction_parameters=self.primer_Tm_chem_correction_parameters,
                    Tm_salt_correction_parameters=self.primer_Tm_salt_correction_parameters,
                )
                dif_Tm = abs(Tm_forward_primer - Tm_reverse_primer)
                if dif_Tm < min_dif_Tm:
                    min_dif_Tm = dif_Tm
                    forward_primer_sequence = primer_attributes["oligo"]

        os.remove(file_fasta_encoding_probes_database)

        return reverse_primer_sequence, forward_primer_sequence

    def generate_output(
        self,
        encoding_probe_database: OligoDatabase,
        reverse_primer_sequence: str,
        forward_primer_sequence: str,
        top_n_sets: int,
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
            "sequence_seqfish_plus_probe",
            "sequence_encoding_probe",
            "sequence_readout_probe_1",
            "sequence_readout_probe_2",
            "sequence_readout_probe_3",
            "sequence_readout_probe_4",
            "sequence_forward_primer",
            "sequence_reverse_primer",
            "sequence_target",
            "sequence_target_probe",
            "length",
            "GC_content",
            "target_probe_isoform_consensus",
        ],
    ) -> None:
        """
        Generate the final output files for the SeqFishPlus probe design pipeline.

        This method updates the encoding probe database with primer sequences, computes
        additional attributes, and writes the results to YAML files, including a file
        for probe order information.

        :param encoding_probe_database: Database of encoding probes with associated attributes and sequences.
        :type encoding_probe_database: OligoDatabase
        :param reverse_primer_sequence: Sequence of the reverse primer.
        :type reverse_primer_sequence: str
        :param forward_primer_sequence: Sequence of the forward primer.
        :type forward_primer_sequence: str
        :param top_n_sets: Number of top oligo sets to include in the output.
        :type top_n_sets: int
        :param attributes: List of attributes to include in the final output YAML file.
        :type attributes: list[str], optional

        :return: None
        """
        new_probe_attributes_primer = {}

        for region_id in encoding_probe_database.database.keys():
            for probe_id in encoding_probe_database.database[region_id].keys():
                new_probe_attributes_primer[probe_id] = {
                    "sequence_reverse_primer": reverse_primer_sequence,
                    "sequence_forward_primer": forward_primer_sequence,
                    "sequence_seqfish_plus_probe": forward_primer_sequence
                    + encoding_probe_database.get_oligo_attribute_value(
                        attribute="sequence_encoding_probe",
                        region_id=region_id,
                        oligo_id=probe_id,
                        flatten=True,
                    )
                    + reverse_primer_sequence,
                }
        encoding_probe_database.update_oligo_attributes(new_probe_attributes_primer)

        encoding_probe_database = self.oligo_attributes_calculator.calculate_oligo_length(
            oligo_database=encoding_probe_database
        )
        encoding_probe_database = self.oligo_attributes_calculator.calculate_GC_content(
            oligo_database=encoding_probe_database, sequence_type="oligo"
        )
        encoding_probe_database = self.oligo_attributes_calculator.calculate_num_targeted_transcripts(
            oligo_database=encoding_probe_database
        )
        encoding_probe_database = self.oligo_attributes_calculator.calculate_isoform_consensus(
            oligo_database=encoding_probe_database
        )

        encoding_probe_database.write_oligosets_to_yaml(
            attributes=attributes,
            top_n_sets=top_n_sets,
            ascending=True,
            filename="seqfish_plus_probes.yml",
        )

        # write a second file that only contains order information
        yaml_dict_order = {}

        for region_id in encoding_probe_database.database.keys():
            yaml_dict_order[region_id] = {}
            oligosets_region = encoding_probe_database.oligosets[region_id]
            oligosets_oligo_columns = [col for col in oligosets_region.columns if col.startswith("oligo_")]
            oligosets_score_columns = [col for col in oligosets_region.columns if col.startswith("score_")]

            oligosets_region.sort_values(by=oligosets_score_columns, ascending=True)
            oligosets_region = oligosets_region.head(top_n_sets)[oligosets_oligo_columns]
            oligosets_region.reset_index(inplace=True, drop=True)

            # iterate through all oligo sets
            for oligoset_idx, oligoset in oligosets_region.iterrows():
                oligoset_id = f"oligoset_{oligoset_idx + 1}"
                yaml_dict_order[region_id][oligoset_id] = {}
                for oligo_id in oligoset:
                    yaml_dict_order[region_id][oligoset_id][oligo_id] = {
                        "sequence_seqfish_plus_probe": encoding_probe_database.get_oligo_attribute_value(
                            attribute="sequence_seqfish_plus_probe",
                            region_id=region_id,
                            oligo_id=oligo_id,
                            flatten=True,
                        ),
                        "sequence_readout_probe_1": encoding_probe_database.get_oligo_attribute_value(
                            attribute="sequence_readout_probe_1",
                            region_id=region_id,
                            oligo_id=oligo_id,
                            flatten=True,
                        ),
                        "sequence_readout_probe_2": encoding_probe_database.get_oligo_attribute_value(
                            attribute="sequence_readout_probe_2",
                            region_id=region_id,
                            oligo_id=oligo_id,
                            flatten=True,
                        ),
                        "sequence_readout_probe_3": encoding_probe_database.get_oligo_attribute_value(
                            attribute="sequence_readout_probe_3",
                            region_id=region_id,
                            oligo_id=oligo_id,
                            flatten=True,
                        ),
                        "sequence_readout_probe_4": encoding_probe_database.get_oligo_attribute_value(
                            attribute="sequence_readout_probe_4",
                            region_id=region_id,
                            oligo_id=oligo_id,
                            flatten=True,
                        ),
                    }

        with open(os.path.join(self.dir_output, "seqfish_plus_probes_order.yml"), "w") as outfile:
            yaml.dump(yaml_dict_order, outfile, default_flow_style=False, sort_keys=False)

        logging.info("--------------END PIPELINE--------------")


############################################
# SeqFish Plus Target Probe Designer
############################################


class SeqFishTargetProbeDesigner:
    """
    Class for designing target probes for SeqFishPlus experiments.
    """
    def __init__(self, dir_output: str, n_jobs: int) -> None:
        """Constructor for the SeqFishPlusProbeDesigner class."""

        ##### create the output folder #####
        self.dir_output = os.path.abspath(dir_output)
        self.subdir_db_oligos = "db_probes"
        self.subdir_db_reference = "db_reference"

        self.n_jobs = n_jobs
        self.oligo_attributes_calculator = OligoAttributes()

    @pipeline_step_basic(step_name="Target Probe Generation - Create Database")
    def create_oligo_database(
        self,
        gene_ids: list,
        oligo_length_min: int,
        oligo_length_max: int,
        files_fasta_oligo_database: list[str],
        min_oligos_per_gene: int,
        isoform_consensus: float,
    ) -> Tuple[OligoDatabase, str]:
        """
        Create an oligo database by generating sequences and applying pre-filters.

        This method uses a sliding window approach to generate oligo sequences from input FASTA files,
        creates a database of oligos, and applies pre-filters based on isoform consensus and minimum
        oligos per gene.

        :param gene_ids: List of gene IDs to include in the oligo database.
        :type gene_ids: list
        :param oligo_length_min: Minimum length of the oligos to generate.
        :type oligo_length_min: int
        :param oligo_length_max: Maximum length of the oligos to generate.
        :type oligo_length_max: int
        :param files_fasta_oligo_database: List of FASTA files containing sequences for oligo generation.
        :type files_fasta_oligo_database: list[str]
        :param min_oligos_per_gene: Minimum number of oligos required per gene in the database.
        :type min_oligos_per_gene: int
        :param isoform_consensus: Minimum isoform consensus threshold for filtering oligos.
        :type isoform_consensus: float

        :return: A tuple containing:
            - `oligo_database`: The created and filtered oligo database.
            - `dir`: Directory containing intermediate output from the sequence generation step.
        :rtype: Tuple[OligoDatabase, str]
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
            min_oligos_per_region=min_oligos_per_gene,
            write_regions_with_insufficient_oligos=True,
            lru_db_max_in_memory=self.n_jobs * 2 + 2,
            database_name=self.subdir_db_oligos,
            dir_output=self.dir_output,
            n_jobs=1,
        )
        oligo_database.load_database_from_fasta(
            files_fasta=oligo_fasta_file,
            database_overwrite=True,
            sequence_type="target",
            region_ids=gene_ids,
        )

        ##### pre-filter oligo database for certain attributes #####
        oligo_database = self.oligo_attributes_calculator.calculate_isoform_consensus(
            oligo_database=oligo_database
        )
        oligo_database.filter_database_by_attribute_threshold(
            attribute_name="target_probe_isoform_consensus",
            attribute_thr=isoform_consensus,
            remove_if_smaller_threshold=True,
        )
        oligo_database.remove_regions_with_insufficient_oligos(pipeline_step="Pre-Filters")

        dir = oligo_sequences.dir_output
        shutil.rmtree(dir) if os.path.exists(dir) else None

        return oligo_database

    @pipeline_step_basic(step_name="Target Probe Generation - Property Filters")
    def filter_by_property(
        self,
        oligo_database: OligoDatabase,
        GC_content_min: float,
        GC_content_max: float,
        homopolymeric_base_n: str,
        T_secondary_structure: float,
        secondary_structures_threshold_deltaG: float,
    ) -> Tuple[OligoDatabase, str]:
        """
        Filter an oligo database based on various sequence properties.

        This method applies filters to remove oligos that do not meet specified criteria, including
        GC content, homopolymeric base runs, and secondary structure thresholds.

        :param oligo_database: The oligo database to be filtered.
        :type oligo_database: OligoDatabase
        :param GC_content_min: Minimum acceptable GC content for oligos.
        :type GC_content_min: float
        :param GC_content_max: Maximum acceptable GC content for oligos.
        :type GC_content_max: float
        :param homopolymeric_base_n: Maximum allowable length of homopolymeric base runs.
        :type homopolymeric_base_n: str
        :param T_secondary_structure: Maximum allowable melting temperature for secondary structures.
        :type T_secondary_structure: float
        :param secondary_structures_threshold_deltaG: Threshold for the free energy (deltaG) of secondary structures.
        :type secondary_structures_threshold_deltaG: float

        :return: The filtered oligo database.
        :rtype: OligoDatabase
        """
        # define the filters
        hard_masked_sequences = HardMaskedSequenceFilter()
        soft_masked_sequences = SoftMaskedSequenceFilter()
        gc_content = GCContentFilter(GC_content_min=GC_content_min, GC_content_max=GC_content_max)
        homopolymeric_runs = HomopolymericRunsFilter(
            base_n=homopolymeric_base_n,
        )
        secondary_sctructure = SecondaryStructureFilter(
            T=T_secondary_structure,
            thr_DG=secondary_structures_threshold_deltaG,
        )

        filters = [
            hard_masked_sequences,
            soft_masked_sequences,
            homopolymeric_runs,
            gc_content,
            secondary_sctructure,
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

    @pipeline_step_basic(step_name="Target Probe Generation - Specificity Filters")
    def filter_by_specificity(
        self,
        oligo_database: OligoDatabase,
        files_fasta_reference_database: List[str],
        specificity_blastn_search_parameters: dict,
        specificity_blastn_hit_parameters: dict,
        cross_hybridization_blastn_search_parameters: dict,
        cross_hybridization_blastn_hit_parameters: dict,
    ) -> Tuple[OligoDatabase, str]:
        """
        Filter an oligo database based on sequence specificity.

        This method applies specificity filters, including exact match, BLASTN specificity, and
        cross-hybridization filters, to ensure oligos target only the intended sequences.

        :param oligo_database: The oligo database to be filtered.
        :type oligo_database: OligoDatabase
        :param files_fasta_reference_database: List of FASTA files containing reference sequences for specificity filtering.
        :type files_fasta_reference_database: list[str]
        :param specificity_blastn_search_parameters: Parameters for BLASTN specificity search.
        :type specificity_blastn_search_parameters: dict
        :param specificity_blastn_hit_parameters: Parameters for filtering BLASTN specificity hits.
        :type specificity_blastn_hit_parameters: dict
        :param cross_hybridization_blastn_search_parameters: Parameters for BLASTN cross-hybridization search.
        :type cross_hybridization_blastn_search_parameters: dict
        :param cross_hybridization_blastn_hit_parameters: Parameters for filtering BLASTN cross-hybridization hits.
        :type cross_hybridization_blastn_hit_parameters: dict

        :return: The filtered oligo database.
        :rtype: OligoDatabase
        """
        ##### define reference database #####
        reference_database = ReferenceDatabase(
            database_name=self.subdir_db_reference, dir_output=self.dir_output
        )
        reference_database.load_database_from_fasta(
            files_fasta=files_fasta_reference_database, database_overwrite=False
        )

        ##### exact match filter #####
        exact_matches = ExactMatchFilter(policy=RemoveByLargerRegionPolicy(), filter_name="oligo_exact_match")

        ##### specificity filters #####
        specificity = BlastNFilter(
            search_parameters=specificity_blastn_search_parameters,
            hit_parameters=specificity_blastn_hit_parameters,
            filter_name="oligo_blastn_specificity",
            dir_output=self.dir_output,
        )

        cross_hybridization_aligner = BlastNFilter(
            search_parameters=cross_hybridization_blastn_search_parameters,
            hit_parameters=cross_hybridization_blastn_hit_parameters,
            filter_name="oligo_blastn_crosshybridization",
            dir_output=self.dir_output,
        )
        cross_hybridization = CrossHybridizationFilter(
            policy=RemoveByLargerRegionPolicy(),
            alignment_method=cross_hybridization_aligner,
            database_name_reference=self.subdir_db_reference,
            dir_output=self.dir_output,
        )

        filters = [exact_matches, specificity, cross_hybridization]
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
            specificity.dir_output,
        ]:
            if os.path.exists(directory):
                shutil.rmtree(directory)

        return oligo_database

    @pipeline_step_basic(step_name="Target Probe Generation - Set Selection")
    def create_oligo_sets(
        self,
        oligo_database: OligoDatabase,
        GC_weight: float,
        GC_content_opt: float,
        UTR_weight: float,
        set_size_opt: int,
        set_size_min: int,
        max_graph_size: int,
        n_sets: int,
        n_attempts: int,
        pre_filter: bool,
        heuristic: bool,
        heuristic_n_attempts: int,
        distance_between_oligos: int,
    ) -> Tuple[OligoDatabase, str, str]:
        """
        Create optimal oligo sets based on weighted scoring criteria.
        
        This method generates oligo sets based on a weighted scoring function that considers GC content,
        untranslated region (UTR) properties, and oligo set size. The sets are optimized using a graph-based
        selection policy and a heuristic approach.

        :param oligo_database: The oligo database from which sets will be created.
        :type oligo_database: OligoDatabase
        :param GC_weight: Weight assigned to GC content in the scoring function.
        :type GC_weight: float
        :param GC_content_opt: Optimal GC content for scoring.
        :type GC_content_opt: float
        :param UTR_weight: Weight assigned to untranslated region (UTR) properties in the scoring function.
        :type UTR_weight: float
        :param set_size_opt: Optimal size of each oligo set.
        :type set_size_opt: int
        :param set_size_min: Minimum size of each oligo set.
        :type set_size_min: int
        :param max_graph_size: Maximum size of the graph used in the selection policy.
        :type max_graph_size: int
        :param n_sets: Number of oligo sets to generate.
        :type n_sets: int
        :param n_attempts: Number of attempts for set optimization.
        :type n_attempts: int
        :param pre_filter: Whether to apply pre-filtering during the selection process.
        :type pre_filter: bool
        :param heuristic: Whether to use a heuristic approach for set optimization.
        :type heuristic: bool
        :param heuristic_n_attempts: Number of attempts for heuristic optimization.
        :type heuristic_n_attempts: int
        :param distance_between_oligos: Minimum distance between oligos in a set.
        :type distance_between_oligos: int

        :return: The updated oligo database with the generated oligo sets.
        :rtype: OligoDatabase
        """
        oligos_scoring = WeightedGCUtrScoring(
            GC_content_opt=GC_content_opt, GC_weight=GC_weight, UTR_weight=UTR_weight
        )
        set_scoring = LowestSetScoring(ascending=True)

        selection_policy = GraphBasedSelectionPolicy(
            set_scoring=set_scoring,
            pre_filter=pre_filter,
            n_attempts=n_attempts,
            heuristic=heuristic,
            heuristic_n_attempts=heuristic_n_attempts,
        )
        probeset_generator = OligosetGeneratorIndependentSet(
            selection_policy=selection_policy,
            oligos_scoring=oligos_scoring,
            set_scoring=set_scoring,
            max_oligos=max_graph_size,
            distance_between_oligos=distance_between_oligos,
        )
        oligo_database = probeset_generator.apply(
            oligo_database=oligo_database,
            sequence_type="oligo",
            set_size_opt=set_size_opt,
            set_size_min=set_size_min,
            n_sets=n_sets,
            n_jobs=self.n_jobs,
        )

        return oligo_database


############################################
# SeqFish Plus Readout Probe Designer
############################################


class SeqFishPlusReadoutProbeDesigner:
    def __init__(
        self,
        dir_output: str,
        n_jobs: int,
    ) -> None:
        """Constructor for the SeqFishPlusReadoutProbeDesigner class."""

        ##### create the output folder #####
        self.dir_output = os.path.abspath(dir_output)
        self.subdir_db_oligos = "db_readout_probes"
        self.subdir_db_reference = "db_reference"

        self.n_jobs = n_jobs

    @pipeline_step_basic(step_name="Readout Probe Generation - Create Oligo Database")
    def create_oligo_database(
        self,
        oligo_length: int,
        oligo_base_probabilities: dict,
        initial_num_sequences: int,
    ) -> Tuple[OligoDatabase, str]:
        ##### creating the oligo sequences #####
        oligo_sequences = OligoSequenceGenerator(dir_output=self.dir_output)
        oligo_fasta_file = oligo_sequences.create_sequences_random(
            filename_out="readout_probe_sequences",
            length_sequences=oligo_length,
            num_sequences=initial_num_sequences,
            name_sequences="readout_probe",
            base_alphabet_with_probability=oligo_base_probabilities,
        )

        ##### creating the oligo database #####
        oligo_database = OligoDatabase(
            min_oligos_per_region=0,
            write_regions_with_insufficient_oligos=False,
            lru_db_max_in_memory=self.n_jobs * 2 + 2,
            database_name=self.subdir_db_oligos,
            dir_output=self.dir_output,
            n_jobs=1,
        )
        oligo_database.load_database_from_fasta(
            files_fasta=oligo_fasta_file,
            database_overwrite=True,
            sequence_type="oligo",
            region_ids=None,
        )

        dir = oligo_sequences.dir_output
        shutil.rmtree(dir) if os.path.exists(dir) else None

        return oligo_database

    @pipeline_step_basic(step_name="Readout Probe Generation - Property Filters")
    def filter_by_property(
        self,
        oligo_database: OligoDatabase,
        GC_content_min: float,
        GC_content_max: float,
        homopolymeric_base_n: int,
    ) -> Tuple[OligoDatabase, str]:
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
            oligo_database=oligo_database,
            sequence_type="oligo",
            n_jobs=self.n_jobs,
        )

        return oligo_database

    @pipeline_step_basic(step_name="Readout Probe Generation - Specificity Filters")
    def filter_by_specificity(
        self,
        oligo_database: OligoDatabase,
        files_fasta_reference_database: list,
        specificity_blastn_search_parameters: dict,
        specificity_blastn_hit_parameters: dict,
        cross_hybridization_blastn_search_parameters: dict,
        cross_hybridization_blastn_hit_parameters: dict,
    ) -> Tuple[OligoDatabase, str]:
        reference_database = ReferenceDatabase(
            database_name=self.subdir_db_reference, dir_output=self.dir_output
        )
        reference_database.load_database_from_fasta(
            files_fasta=files_fasta_reference_database, database_overwrite=False
        )

        ##### specificity filters #####
        # removing duplicated oligos from the region with the most oligos
        exact_matches = ExactMatchFilter(
            policy=RemoveByDegreePolicy(), filter_name="readout_probes_exact_match"
        )

        # BlastN Filter
        specificity = BlastNFilter(
            search_parameters=specificity_blastn_search_parameters,
            hit_parameters=specificity_blastn_hit_parameters,
            filter_name="readout_probes_blastn_specificity",
            dir_output=self.dir_output,
        )

        # Cross-Hybridization Filter
        cross_hybridization_aligner = BlastNFilter(
            search_parameters=cross_hybridization_blastn_search_parameters,
            hit_parameters=cross_hybridization_blastn_hit_parameters,
            filter_name="readout_probes_blastn_crosshybridization",
            dir_output=self.dir_output,
        )
        cross_hybridization = CrossHybridizationFilter(
            policy=RemoveByDegreePolicy(),
            alignment_method=cross_hybridization_aligner,
            database_name_reference=self.subdir_db_reference,
            dir_output=self.dir_output,
        )

        filters = [exact_matches, specificity, cross_hybridization]
        specificity_filter = SpecificityFilter(filters=filters)
        oligo_database = specificity_filter.apply(
            sequence_type="oligo",
            oligo_database=oligo_database,
            reference_database=reference_database,
            n_jobs=self.n_jobs,
        )

        for directory in [
            reference_database.dir_output,
            specificity.dir_output,
            cross_hybridization_aligner.dir_output,
            cross_hybridization.dir_output,
        ]:
            if os.path.exists(directory):
                shutil.rmtree(directory)

        return oligo_database

    def generate_codebook(
        self, n_regions: int, n_barcode_rounds: int, n_pseudocolors: int, n_channels: int
    ) -> pd.DataFrame:
        def _generate_barcode(pseudocolors: list, channel: int, n_pseudocolors: int, n_channels: int) -> list:
            pseudocolors = pseudocolors + [sum(pseudocolors) % n_pseudocolors]
            assert n_pseudocolors > max(
                pseudocolors
            ), f"The number of pseudocolor is {n_pseudocolors}, while the barcode contains {max(pseudocolors)} pseudocolors."
            assert (
                n_channels > channel
            ), f"The number of channles is {n_channels}, while the barcode contains {channel} channels."
            n_barcode_rounds = len(pseudocolors)
            barcode = np.zeros(n_channels * n_pseudocolors * n_barcode_rounds, dtype=np.int8)
            for i, pseudocolor in enumerate(pseudocolors):
                barcode[i * n_pseudocolors * n_channels + n_channels * pseudocolor + channel] = 1
            return barcode

        codebook = []
        codebook_size = n_channels * (n_pseudocolors ** (n_barcode_rounds - 1))
        barcode_size = n_pseudocolors * n_barcode_rounds * n_channels
        if codebook_size < n_regions:
            raise ValueError(
                f"The number of valid barcodes ({codebook_size}) is lower than the number of regions ({n_regions}). Consider increasing the number of psudocolors or barcoding rounds."
            )
        for pseudocolors in product(range(n_pseudocolors), repeat=n_barcode_rounds - 1):
            pseudocolors = list(pseudocolors)
            for channel in range(n_channels):
                barcode = _generate_barcode(
                    pseudocolors=pseudocolors,
                    channel=channel,
                    n_pseudocolors=n_pseudocolors,
                    n_channels=n_channels,
                )
                codebook.append(barcode)

        codebook = pd.DataFrame(codebook, columns=[f"bit_{i+1}" for i in range(barcode_size)])

        return codebook

    def create_readout_probe_table(
        self,
        readout_probe_database: OligoDatabase,
        channels_ids: list,
        n_barcode_rounds: int,
        n_pseudocolors: int,
    ) -> pd.DataFrame:
        n_channels = len(channels_ids)
        n_bits = n_barcode_rounds * n_pseudocolors * n_channels
        readout_probes = readout_probe_database.get_oligoid_sequence_mapping(
            sequence_type="oligo", sequence_to_upper=False
        )
        assert (
            len(readout_probes) >= n_bits
        ), f"There are less readout probes ({len(readout_probes)}) than bits ({n_bits})."
        readout_probe_table = pd.DataFrame(
            columns=[
                "bit",
                "barcode_round",
                "pseudocolor",
                "channel",
                "readout_probe_sequence",
            ],
            index=list(range(n_bits)),
        )
        barcode_round = 0
        pseudocolor = 0
        channel = 0
        for i, readout_probe_sequence in enumerate(readout_probes.values()):
            readout_probe_table.iloc[i] = [
                f"bit_{i + 1}",
                barcode_round + 1,
                pseudocolor + 1,
                channels_ids[channel],
                readout_probe_sequence,
            ]
            channel = (channel + 1) % n_channels
            if channel == 0:
                pseudocolor = (pseudocolor + 1) % n_pseudocolors
                if pseudocolor == 0:
                    barcode_round = (barcode_round + 1) % n_barcode_rounds
            if i >= n_bits - 1:
                break
        readout_probe_table.set_index("bit", inplace=True)
        return readout_probe_table


############################################
# SeqFish Plus Primer Designer
############################################


class SeqFishPlusPrimerDesigner:
    def __init__(
        self,
        dir_output: str,
        n_jobs: int,
    ) -> None:
        """Constructor for the SeqFishPlusPrimerDesigner class."""

        ##### create the output folder #####
        self.dir_output = os.path.abspath(dir_output)
        self.subdir_db_oligos = "db_primer"
        self.subdir_db_reference = "db_reference"

        self.n_jobs = n_jobs

    @pipeline_step_basic(step_name="Primer Generation - Create Oligo Database")
    def create_oligo_database(
        self,
        oligo_length: int,
        oligo_base_probabilities: dict,
        initial_num_sequences: int,
    ) -> Tuple[OligoDatabase, str]:
        ##### creating the primer sequences #####
        # random forward primer
        forward_primer_sequences = OligoSequenceGenerator(dir_output=self.dir_output)
        forward_primer_fasta_file = forward_primer_sequences.create_sequences_random(
            filename_out="forward_primer_sequences",
            length_sequences=oligo_length - 1,
            num_sequences=initial_num_sequences,
            name_sequences="forward_primer",
            base_alphabet_with_probability=oligo_base_probabilities,
        )
        # we want to keep primer which end with a specific nucleotide, i.e. "T"
        forward_primer_fasta_file = append_nucleotide_to_sequences(forward_primer_fasta_file, nucleotide="T")

        ##### creating the primer database #####
        oligo_database = OligoDatabase(
            min_oligos_per_region=0,
            write_regions_with_insufficient_oligos=False,
            lru_db_max_in_memory=self.n_jobs * 2 + 2,
            database_name=self.subdir_db_oligos,
            dir_output=self.dir_output,
            n_jobs=1,
        )
        oligo_database.load_database_from_fasta(
            files_fasta=forward_primer_fasta_file,
            database_overwrite=True,
            sequence_type="oligo",
            region_ids=None,
        )

        dir = forward_primer_sequences.dir_output
        shutil.rmtree(dir) if os.path.exists(dir) else None

        return oligo_database

    @pipeline_step_basic(step_name="Primer Generation - Property Filters")
    def filter_by_property(
        self,
        oligo_database: OligoDatabase,
        GC_content_min: float,
        GC_content_max: float,
        number_GC_GCclamp: int,
        number_three_prime_base_GCclamp: int,
        homopolymeric_base_n: int,
        max_len_selfcomplement: int,
        reverse_primer_sequence: str,
        max_len_complement: int,
        Tm_min: float,
        Tm_max: float,
        Tm_parameters: dict,
        Tm_chem_correction_parameters: dict,
        Tm_salt_correction_parameters: dict,
        T_secondary_structure: float,
        secondary_structures_threshold_deltaG: float,
    ) -> Tuple[OligoDatabase, str]:
        # define the filters
        gc_content = GCContentFilter(GC_content_min=GC_content_min, GC_content_max=GC_content_max)
        gc_clamp = GCClampFilter(n_bases=number_three_prime_base_GCclamp, n_GC=number_GC_GCclamp)
        homopolymeric_runs = HomopolymericRunsFilter(
            base_n=homopolymeric_base_n,
        )
        self_complement = SelfComplementFilter(max_len_selfcomplement=max_len_selfcomplement)
        complement = ComplementFilter(
            comparison_sequence=reverse_primer_sequence, max_len_complement=max_len_complement
        )
        melting_temperature = MeltingTemperatureNNFilter(
            Tm_min=Tm_min,
            Tm_max=Tm_max,
            Tm_parameters=Tm_parameters,
            Tm_chem_correction_parameters=Tm_chem_correction_parameters,
            Tm_salt_correction_parameters=Tm_salt_correction_parameters,
        )
        secondary_sctructure = SecondaryStructureFilter(
            T=T_secondary_structure,
            thr_DG=secondary_structures_threshold_deltaG,
        )

        filters = [
            gc_content,
            gc_clamp,
            homopolymeric_runs,
            self_complement,
            complement,
            melting_temperature,
            secondary_sctructure,
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

    @pipeline_step_basic(step_name="Primer Generation - Specificity Filters")
    def filter_by_specificity(
        self,
        oligo_database: OligoDatabase,
        files_fasta_reference_database: List[str],
        specificity_refrence_blastn_search_parameters: dict,
        specificity_refrence_blastn_hit_parameters: dict,
        file_fasta_encoding_probes_database: str,
        specificity_encoding_probes_blastn_search_parameters: dict,
        specificity_encoding_probes_blastn_hit_parameters: dict,
    ) -> Tuple[OligoDatabase, str]:
        ##### specificity filters against reference #####
        reference_database = ReferenceDatabase(
            database_name=self.subdir_db_reference, dir_output=self.dir_output
        )
        reference_database.load_database_from_fasta(
            files_fasta=files_fasta_reference_database, database_overwrite=True
        )
        # BlastN Filter
        specificity_refrence = BlastNFilter(
            search_parameters=specificity_refrence_blastn_search_parameters,
            hit_parameters=specificity_refrence_blastn_hit_parameters,
            filter_name="primer_blastn_specificity_reference",
            dir_output=self.dir_output,
        )

        specificity_filter_reference = SpecificityFilter(filters=[specificity_refrence])
        oligo_database = specificity_filter_reference.apply(
            sequence_type="oligo",
            oligo_database=oligo_database,
            reference_database=reference_database,
            n_jobs=self.n_jobs,
        )

        ##### specificity filters against encoding probes #####
        encoding_probes_database = ReferenceDatabase(
            database_name=self.subdir_db_reference, dir_output=self.dir_output
        )
        encoding_probes_database.load_database_from_fasta(
            files_fasta=file_fasta_encoding_probes_database, database_overwrite=True
        )
        # BlastN Filter
        specificity_encoding_probes = BlastNFilter(
            search_parameters=specificity_encoding_probes_blastn_search_parameters,
            hit_parameters=specificity_encoding_probes_blastn_hit_parameters,
            filter_name="primer_blastn_specificity_encoding_probes",
            dir_output=self.dir_output,
        )

        specificity_filter_encoding_probes = SpecificityFilter(filters=[specificity_encoding_probes])
        oligo_database = specificity_filter_encoding_probes.apply(
            sequence_type="oligo",
            oligo_database=oligo_database,
            reference_database=encoding_probes_database,
            n_jobs=self.n_jobs,
        )

        for directory in [
            reference_database.dir_output,
            encoding_probes_database.dir_output,
            specificity_refrence.dir_output,
            specificity_encoding_probes.dir_output,
        ]:
            if os.path.exists(directory):
                shutil.rmtree(directory)

        return oligo_database


############################################
# SeqFish Plus Probe Designer Pipeline
############################################


def main():
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
    pipeline = SeqFishPlusProbeDesigner(
        write_intermediate_steps=config["write_intermediate_steps"],
        dir_output=config["dir_output"],
        n_jobs=config["n_jobs"],
    )

    ##### set custom developer parameters #####
    pipeline.set_developer_parameters(
        target_probe_specificity_blastn_search_parameters=config[
            "target_probe_specificity_blastn_search_parameters"
        ],
        target_probe_specificity_blastn_hit_parameters=config[
            "target_probe_specificity_blastn_hit_parameters"
        ],
        target_probe_cross_hybridization_blastn_search_parameters=config[
            "target_probe_cross_hybridization_blastn_search_parameters"
        ],
        target_probe_cross_hybridization_blastn_hit_parameters=config[
            "target_probe_cross_hybridization_blastn_hit_parameters"
        ],
        readout_probe_initial_num_sequences=config["readout_probe_initial_num_sequences"],
        readout_probe_specificity_blastn_search_parameters=config[
            "readout_probe_specificity_blastn_search_parameters"
        ],
        readout_probe_specificity_blastn_hit_parameters=config[
            "readout_probe_specificity_blastn_hit_parameters"
        ],
        readout_probe_cross_hybridization_blastn_search_parameters=config[
            "readout_probe_cross_hybridization_blastn_search_parameters"
        ],
        readout_probe_cross_hybridization_blastn_hit_parameters=config[
            "readout_probe_cross_hybridization_blastn_hit_parameters"
        ],
        primer_initial_num_sequences=config["primer_initial_num_sequences"],
        primer_specificity_refrence_blastn_search_parameters=config[
            "primer_specificity_refrence_blastn_search_parameters"
        ],
        primer_specificity_refrence_blastn_hit_parameters=config[
            "primer_specificity_refrence_blastn_hit_parameters"
        ],
        primer_specificity_encoding_probes_blastn_search_parameters=config[
            "primer_specificity_encoding_probes_blastn_search_parameters"
        ],
        primer_specificity_encoding_probes_blastn_hit_parameters=config[
            "primer_specificity_encoding_probes_blastn_hit_parameters"
        ],
        primer_Tm_parameters=config["primer_Tm_parameters"],
        primer_Tm_chem_correction_parameters=config["primer_Tm_chem_correction_parameters"],
        primer_Tm_salt_correction_parameters=config["primer_Tm_salt_correction_parameters"],
        max_graph_size=config["max_graph_size"],
        pre_filter=config["pre_filter"],
        n_attempts=config["n_attempts"],
        heuristic=config["heuristic"],
        heuristic_n_attempts=config["heuristic_n_attempts"],
    )

    ##### design probes #####
    target_probe_database = pipeline.design_target_probes(
        files_fasta_target_probe_database=config["files_fasta_target_probe_database"],
        files_fasta_reference_database_targe_probe=config["files_fasta_reference_database_targe_probe"],
        gene_ids=gene_ids,
        target_probe_length_min=config["target_probe_length_min"],
        target_probe_length_max=config["target_probe_length_max"],
        target_probe_isoform_consensus=config["target_probe_isoform_consensus"],
        target_probe_GC_content_min=config["target_probe_GC_content_min"],
        target_probe_GC_content_opt=config["target_probe_GC_content_opt"],
        target_probe_GC_content_max=config["target_probe_GC_content_max"],
        target_probe_homopolymeric_base_n=config["target_probe_homopolymeric_base_n"],
        target_probe_T_secondary_structure=config["target_probe_T_secondary_structure"],
        target_probe_secondary_structures_threshold_deltaG=config[
            "target_probe_secondary_structures_threshold_deltaG"
        ],
        target_probe_GC_weight=config["target_probe_GC_weight"],
        target_probe_UTR_weight=config["target_probe_UTR_weight"],
        set_size_opt=config["set_size_opt"],
        set_size_min=config["set_size_min"],
        distance_between_target_probes=config["distance_between_target_probes"],
        n_sets=config["n_sets"],
    )

    codebook, readout_probe_table = pipeline.design_readout_probes(
        n_genes=len(target_probe_database.database),
        files_fasta_reference_database_readout_probe=config["files_fasta_reference_database_readout_probe"],
        readout_probe_length=config["readout_probe_length"],
        readout_probe_base_probabilities=config["readout_probe_base_probabilities"],
        readout_probe_GC_content_min=config["readout_probe_GC_content_min"],
        readout_probe_GC_content_max=config["readout_probe_GC_content_max"],
        readout_probe_homopolymeric_base_n=config["readout_probe_homopolymeric_base_n"],
        n_barcode_rounds=config["n_barcode_rounds"],
        n_pseudocolors=config["n_pseudocolors"],
        channels_ids=config["channels_ids"],
    )

    encoding_probe_database = pipeline.design_encoding_probe(
        target_probe_database=target_probe_database,
        codebook=codebook,
        readout_probe_table=readout_probe_table,
    )

    reverse_primer_sequence, forward_primer_sequence = pipeline.design_primers(
        encoding_probe_database=encoding_probe_database,
        files_fasta_reference_database_primer=config["files_fasta_reference_database_primer"],
        reverse_primer_sequence=config["reverse_primer_sequence"],
        primer_length=config["primer_length"],
        primer_base_probabilities=config["primer_base_probabilities"],
        primer_GC_content_min=config["primer_GC_content_min"],
        primer_GC_content_max=config["primer_GC_content_max"],
        primer_number_GC_GCclamp=config["primer_number_GC_GCclamp"],
        primer_number_three_prime_base_GCclamp=config["primer_number_three_prime_base_GCclamp"],
        primer_homopolymeric_base_n=config["primer_homopolymeric_base_n"],
        primer_max_len_selfcomplement=config["primer_max_len_selfcomplement"],
        primer_max_len_complement_reverse_primer=config["primer_max_len_complement_reverse_primer"],
        primer_Tm_min=config["primer_Tm_min"],
        primer_Tm_max=config["primer_Tm_max"],
        primer_T_secondary_structure=config["primer_T_secondary_structure"],
        primer_secondary_structures_threshold_deltaG=config["primer_secondary_structures_threshold_deltaG"],
    )

    pipeline.generate_output(
        encoding_probe_database=encoding_probe_database,
        reverse_primer_sequence=reverse_primer_sequence,
        forward_primer_sequence=forward_primer_sequence,
        top_n_sets=config["top_n_sets"],
    )

    print("--------------END PIPELINE--------------")


if __name__ == "__main__":
    main()
