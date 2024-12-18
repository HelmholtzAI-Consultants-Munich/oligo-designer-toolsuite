############################################
# imports
############################################

import os
import yaml
import shutil
import logging
import warnings

import numpy as np
import pandas as pd

from typing import List, Tuple
from pathlib import Path
from datetime import datetime
from itertools import combinations
from scipy.spatial.distance import hamming

from Bio.SeqUtils import Seq
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
    HomogeneousPropertyOligoSetGenerator,
    OligosetGeneratorIndependentSet,
    GraphBasedSelectionPolicy,
    GreedySelectionPolicy,
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
# Merfish Probe Designer
############################################


class MerfishProbeDesigner:
    def __init__(self, write_intermediate_steps: bool, dir_output: str, n_jobs: int) -> None:
        """Constructor for the MerfishProbeDesigner class."""

        ##### create the output folder #####
        self.dir_output = os.path.abspath(dir_output)
        Path(self.dir_output).mkdir(parents=True, exist_ok=True)

        ##### setup logger #####
        timestamp = datetime.now()
        file_logger = os.path.join(
            self.dir_output,
            f"log_Merfish_probe_designer_{timestamp.year}-{timestamp.month}-{timestamp.day}-{timestamp.hour}-{timestamp.minute}.txt",
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
            "perc_identity": 80,
            "strand": "minus",
            "word_size": 10,
            "dust": "no",
            "soft_masking": "false",
            "max_target_seqs": 10,
            "max_hsps": 1000,
        },
        target_probe_specificity_blastn_hit_parameters: dict = {"min_alignment_length": 17},
        target_probe_cross_hybridization_blastn_search_parameters: dict = {
            "perc_identity": 80,
            "strand": "minus",
            "word_size": 7,
            "dust": "no",
            "soft_masking": "false",
            "max_target_seqs": 10,
        },
        target_probe_cross_hybridization_blastn_hit_parameters: dict = {"min_alignment_length": 17},
        target_probe_Tm_parameters: dict = {
            "check": True,
            "strict": True,
            "c_seq": None,
            "shift": 0,
            "nn_table": "DNA_NN4",
            "tmm_table": "DNA_TMM1",
            "imm_table": "DNA_IMM1",
            "de_table": "DNA_DE1",
            "dnac1": 5,
            "dnac2": 0,
            "selfcomp": False,
            "saltcorr": 5,
            "Na": 300,
            "K": 0,
            "Tris": 0,
            "Mg": 0,
            "dNTPs": 0,
        },
        target_probe_Tm_chem_correction_parameters: dict = None,
        target_probe_Tm_salt_correction_parameters: dict = None,
        readout_probe_initial_num_sequences: int = 100000,
        readout_probe_n_combinations: int = 100000,
        readout_probe_specificity_blastn_search_parameters: dict = {
            "perc_identity": 100,
            "strand": "minus",
            "word_size": 7,
            "dust": "no",
            "soft_masking": "false",
            "max_target_seqs": 10,
            "max_hsps": 1000,
        },
        readout_probe_specificity_blastn_hit_parameters: dict = {"min_alignment_length": 11},
        readout_probe_cross_hybridization_blastn_search_parameters: dict = {
            "perc_identity": 100,
            "strand": "minus",
            "word_size": 7,
            "dust": "no",
            "soft_masking": "false",
            "max_target_seqs": 10,
        },
        readout_probe_cross_hybridization_blastn_hit_parameters: dict = {"min_alignment_length": 11},
        readout_probe_Tm_parameters: dict = {
            "check": True,
            "strict": True,
            "c_seq": None,
            "shift": 0,
            "nn_table": "DNA_NN4",
            "tmm_table": "DNA_TMM1",
            "imm_table": "DNA_IMM1",
            "de_table": "DNA_DE1",
            "dnac1": 25,
            "dnac2": 25,
            "selfcomp": False,
            "saltcorr": 5,
            "Na": 300,
            "K": 0,
            "Tris": 0,
            "Mg": 0,
            "dNTPs": 0,
        },
        readout_probe_Tm_chem_correction_parameters: dict = None,
        readout_probe_Tm_salt_correction_parameters: dict = None,
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
        pre_filter: bool = True,
        n_attempts: int = 100000,
        heuristic: bool = True,
        heuristic_n_attempts: int = 100,
    ):
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
        self.readout_probe_n_combinations = readout_probe_n_combinations
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

        # preprocess melting temperature params
        target_probe_Tm_parameters["nn_table"] = getattr(mt, target_probe_Tm_parameters["nn_table"])
        target_probe_Tm_parameters["tmm_table"] = getattr(mt, target_probe_Tm_parameters["tmm_table"])
        target_probe_Tm_parameters["imm_table"] = getattr(mt, target_probe_Tm_parameters["imm_table"])
        target_probe_Tm_parameters["de_table"] = getattr(mt, target_probe_Tm_parameters["de_table"])

        readout_probe_Tm_parameters["nn_table"] = getattr(mt, readout_probe_Tm_parameters["nn_table"])
        readout_probe_Tm_parameters["tmm_table"] = getattr(mt, readout_probe_Tm_parameters["tmm_table"])
        readout_probe_Tm_parameters["imm_table"] = getattr(mt, readout_probe_Tm_parameters["imm_table"])
        readout_probe_Tm_parameters["de_table"] = getattr(mt, readout_probe_Tm_parameters["de_table"])

        primer_Tm_parameters["nn_table"] = getattr(mt, primer_Tm_parameters["nn_table"])
        primer_Tm_parameters["tmm_table"] = getattr(mt, primer_Tm_parameters["tmm_table"])
        primer_Tm_parameters["imm_table"] = getattr(mt, primer_Tm_parameters["imm_table"])
        primer_Tm_parameters["de_table"] = getattr(mt, primer_Tm_parameters["de_table"])

        ## target probe
        self.target_probe_Tm_parameters = target_probe_Tm_parameters
        self.target_probe_Tm_chem_correction_parameters = target_probe_Tm_chem_correction_parameters
        self.target_probe_Tm_salt_correction_parameters = target_probe_Tm_salt_correction_parameters

        ## readout probe
        self.readout_probe_Tm_parameters = readout_probe_Tm_parameters
        self.readout_probe_Tm_chem_correction_parameters = readout_probe_Tm_chem_correction_parameters
        self.readout_probe_Tm_salt_correction_parameters = readout_probe_Tm_salt_correction_parameters

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
        target_probe_length_min: int = 30,
        target_probe_length_max: int = 30,
        target_probe_isoform_consensus: float = 50,
        target_probe_GC_content_min: float = 43,
        target_probe_GC_content_opt: float = 55,
        target_probe_GC_content_max: float = 63,
        target_probe_Tm_min: float = 66,
        target_probe_Tm_opt: float = 72,
        target_probe_Tm_max: float = 76,
        target_probe_homopolymeric_base_n: dict = {"A": 5, "T": 5, "C": 5, "G": 5},
        target_probe_T_secondary_structure: float = 76,
        target_probe_secondary_structures_threshold_deltaG: float = 0,
        target_probe_GC_weight: float = 1,
        target_probe_Tm_weight: float = 1,
        target_probe_isoform_weight: float = 2,
        set_size_opt: int = 50,
        set_size_min: int = 50,
        distance_between_target_probes: int = 0,
        n_sets: int = 100,
    ) -> OligoDatabase:

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
            Tm_min=target_probe_Tm_min,
            Tm_max=target_probe_Tm_max,
            homopolymeric_base_n=target_probe_homopolymeric_base_n,
            T_secondary_structure=target_probe_T_secondary_structure,
            secondary_structures_threshold_deltaG=target_probe_secondary_structures_threshold_deltaG,
            Tm_parameters=self.target_probe_Tm_parameters,
            Tm_chem_correction_parameters=self.target_probe_Tm_chem_correction_parameters,
            Tm_salt_correction_parameters=self.target_probe_Tm_salt_correction_parameters,
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
            isoform_weight=target_probe_isoform_weight,
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
            set_size_opt=set_size_opt,
            set_size_min=set_size_min,
            distance_between_oligos=distance_between_target_probes,
            n_sets=n_sets,
            max_graph_size=self.max_graph_size,
            n_attempts=self.n_attempts,
            pre_filter=self.pre_filter,
            heuristic=self.heuristic,
            heuristic_n_attempts=self.heuristic_n_attempts,
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
        readout_probe_length: int = 20,
        readout_probe_base_probabilities: dict = {"A": 0.25, "C": 0.00, "G": 0.50, "T": 0.25},
        readout_probe_GC_content_min: float = 40,
        readout_probe_GC_content_max: float = 50,
        readout_probe_homopolymeric_base_n: dict = {"G": 3},
        readout_probe_set_size: int = 16,
        readout_probe_homogeneous_properties_weights: dict = {"TmNN": 1, "GC_content": 1},
        n_bits: int = 16,
        min_hamming_dist: int = 4,
        hamming_weight: int = 4,
        channels_ids: list = ["Alexa488", "Cy3b", "Alexa647"],
    ) -> Tuple[pd.DataFrame, pd.DataFrame]:
        readout_probe_designer = MerfishReadoutProbeDesigner(
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

        oligo_database = readout_probe_designer.create_oligo_sets(
            oligo_database=oligo_database,
            set_size=readout_probe_set_size,
            homogeneous_properties_weights=readout_probe_homogeneous_properties_weights,
            n_combinations=self.readout_probe_n_combinations,
            Tm_parameters=self.readout_probe_Tm_parameters,
            Tm_chem_correction_parameters=self.readout_probe_Tm_chem_correction_parameters,
            Tm_salt_correction_parameters=self.readout_probe_Tm_salt_correction_parameters,
        )
        check_content_oligo_database(oligo_database)

        if self.write_intermediate_steps:
            dir_database = oligo_database.save_database(dir_database="4_db_readout_probes_set_selection")
            print(f"Saved readout probe database for step 4 (Set Selection) in directory {dir_database}")

        codebook = readout_probe_designer.generate_codebook(
            n_regions=n_genes,
            n_bits=n_bits,
            min_hamming_dist=min_hamming_dist,
            hamming_weight=hamming_weight,
        )

        readout_probe_table = readout_probe_designer.create_readout_probe_table(
            readout_probe_database=oligo_database,
            channels_ids=channels_ids,
            n_bits=n_bits,
        )

        return codebook, readout_probe_table

    def design_encoding_probe(
        self,
        target_probe_database: OligoDatabase,
        codebook: pd.DataFrame,
        readout_probe_table: pd.DataFrame,
    ):

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
                    "sequence_encoding_probe": (
                        str(Seq(sequence_readout_probe_1).reverse_complement())
                        + "A"
                        + target_probe_database.get_oligo_attribute_value(
                            attribute="oligo", region_id=region_id, oligo_id=probe_id, flatten=True
                        )
                        + "A"
                        + str(Seq(sequence_readout_probe_2).reverse_complement())
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

        file_fasta_encoding_probes_database = encoding_probe_database.write_database_to_fasta(
            filename=f"db_reference_encoding_probes",
            save_description=False,
            region_ids=None,
            sequence_type="sequence_encoding_probe",
        )

        primer_designer = MerfishPrimerDesigner(
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
            "sequence_merfish_probe",
            "sequence_encoding_probe",
            "sequence_readout_probe_1",
            "sequence_readout_probe_2",
            "sequence_forward_primer",
            "sequence_reverse_primer",
            "sequence_target",
            "sequence_target_probe",
            "length",
            "GC_content",
            "target_probe_isoform_consensus",
        ],
    ) -> None:
        new_probe_attributes_primer = {}

        for region_id in encoding_probe_database.database.keys():
            for probe_id in encoding_probe_database.database[region_id].keys():
                new_probe_attributes_primer[probe_id] = {
                    "sequence_reverse_primer": reverse_primer_sequence,
                    "sequence_forward_primer": forward_primer_sequence,
                    "sequence_merfish_probe": forward_primer_sequence
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
            filename="merfish_probes.yml",
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
                        "sequence_merfish_probe": encoding_probe_database.get_oligo_attribute_value(
                            attribute="sequence_merfish_probe",
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
                    }

        with open(os.path.join(self.dir_output, "merfish_probes_order.yml"), "w") as outfile:
            yaml.dump(yaml_dict_order, outfile, default_flow_style=False, sort_keys=False)

        logging.info("--------------END PIPELINE--------------")


############################################
# SeqFish Plus Target Probe Designer
############################################


class SeqFishTargetProbeDesigner:
    def __init__(self, dir_output: str, n_jobs: int) -> None:
        """Constructor for the MerfishProbeDesigner class."""

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
        Tm_min: float,
        Tm_max: float,
        homopolymeric_base_n: str,
        T_secondary_structure: float,
        secondary_structures_threshold_deltaG: float,
        Tm_parameters: dict,
        Tm_chem_correction_parameters: dict,
        Tm_salt_correction_parameters: dict,
    ) -> Tuple[OligoDatabase, str]:
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
        isoform_weight: float,
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
        pre_filter: bool,
        heuristic: bool,
        heuristic_n_attempts: int,
    ) -> Tuple[OligoDatabase, str, str]:
        oligos_scoring = WeightedIsoformTmGCOligoScoring(
            Tm_min=Tm_min,
            Tm_opt=Tm_opt,
            Tm_max=Tm_max,
            GC_content_min=GC_content_min,
            GC_content_opt=GC_content_opt,
            GC_content_max=GC_content_max,
            Tm_parameters=Tm_parameters,
            Tm_chem_correction_parameters=Tm_chem_correction_parameters,
            Tm_salt_correction_parameters=Tm_salt_correction_parameters,
            isoform_weight=isoform_weight,
            Tm_weight=Tm_weight,
            GC_weight=GC_weight,
        )
        set_scoring = LowestSetScoring(ascending=True)

        # selection_policy = GraphBasedSelectionPolicy(
        #     set_scoring=set_scoring,
        #     pre_filter=pre_filter,
        #     n_attempts=n_attempts,
        #     heuristic=heuristic,
        #     heuristic_n_attempts=heuristic_n_attempts,
        # )
        selection_policy = GreedySelectionPolicy(
            set_scoring=set_scoring,
            score_criteria="set_score_worst",
            pre_filter=pre_filter,
            penalty=0.01,
            n_attempts=n_attempts,
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


class MerfishReadoutProbeDesigner:
    def __init__(
        self,
        dir_output: str,
        n_jobs: int,
    ) -> None:
        """Constructor for the MerfishReadoutProbeDesigner class."""

        ##### create the output folder #####
        self.dir_output = os.path.abspath(dir_output)
        self.subdir_db_oligos = "db_readout_probes"
        self.subdir_db_reference = "db_reference"

        self.n_jobs = n_jobs
        self.oligo_attributes_calculator = OligoAttributes()

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

    @pipeline_step_basic(step_name="Readout Probe Generation - Set Selection")
    def create_oligo_sets(
        self,
        oligo_database: OligoDatabase,
        set_size: int,
        homogeneous_properties_weights: dict,
        n_combinations: int,
        Tm_parameters,
        Tm_chem_correction_parameters,
        Tm_salt_correction_parameters,
    ):
        oligo_database = self.oligo_attributes_calculator.calculate_TmNN(
            oligo_database=oligo_database,
            Tm_parameters=Tm_parameters,
            Tm_chem_correction_parameters=Tm_chem_correction_parameters,
            Tm_salt_correction_parameters=Tm_salt_correction_parameters,
            sequence_type="oligo",
        )
        oligo_database = self.oligo_attributes_calculator.calculate_GC_content(oligo_database, "oligo")

        set_generator = HomogeneousPropertyOligoSetGenerator(
            set_size=set_size, properties=homogeneous_properties_weights
        )
        oligo_database = set_generator.apply(
            oligo_database, n_sets=1, n_combinations=n_combinations, n_jobs=self.n_jobs
        )

        return oligo_database

    def generate_codebook(
        self,
        n_regions: int,
        n_bits: int,
        min_hamming_dist: int,
        hamming_weight: int,
    ) -> pd.DataFrame:
        def _generate_barcode(raw_barcode: list, n_bits: int):
            barcode = np.zeros(n_bits, dtype=np.int8)
            for i in raw_barcode:
                barcode[i] = 1
            return barcode

        codebook = []
        for raw_barcode in combinations(iterable=range(n_bits), r=hamming_weight):
            new_barcode = _generate_barcode(raw_barcode=raw_barcode, n_bits=n_bits)
            # check if the barcode passes the requirements
            add_new_barcode = True
            for barcode in codebook:
                hamming_dist = hamming(new_barcode, barcode) * n_bits
                if hamming_dist < min_hamming_dist:
                    add_new_barcode = False
                    break
            if add_new_barcode:
                codebook.append(new_barcode)
        if len(codebook) < n_regions:
            raise ValueError(
                f"The number of valid barcodes ({len(codebook)}) is lower than the number of regions({n_regions}). Consider increasing the number of bits."
            )

        codebook = pd.DataFrame(codebook, columns=[f"bit_{i+1}" for i in range(n_bits)])
        return codebook

    def create_readout_probe_table(
        self, readout_probe_database: OligoDatabase, channels_ids: list, n_bits: int
    ) -> pd.DataFrame:
        readout_probes = readout_probe_database.get_oligoid_sequence_mapping(
            sequence_type="oligo", sequence_to_upper=False
        )
        assert (
            len(readout_probes) >= n_bits
        ), f"There are less readout probes ({len(readout_probes)}) than bits ({n_bits})."
        readout_probe_table = pd.DataFrame(
            columns=["bit", "channel", "readout_probe_id", "readout_probe_sequence"],
            index=list(range(n_bits)),
        )
        n_channels = len(channels_ids)
        channel = 0
        for i, (readout_probe_id, readout_probe_sequence) in enumerate(readout_probes.items()):
            readout_probe_table.iloc[i] = [
                f"bit_{i+1}",
                channels_ids[channel],
                readout_probe_id,
                readout_probe_sequence,
            ]
            channel = (channel + 1) % n_channels
            if i >= n_bits - 1:
                break

        return readout_probe_table


############################################
# SeqFish Plus Primer Designer
############################################


class MerfishPrimerDesigner:
    def __init__(
        self,
        dir_output: str,
        n_jobs: int,
    ) -> None:
        """Constructor for the MerfishPrimerDesigner class."""

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
    pipeline = MerfishProbeDesigner(
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
        target_probe_Tm_parameters=config["target_probe_Tm_parameters"],
        target_probe_Tm_chem_correction_parameters=config["target_probe_Tm_chem_correction_parameters"],
        target_probe_Tm_salt_correction_parameters=config["target_probe_Tm_salt_correction_parameters"],
        readout_probe_initial_num_sequences=config["readout_probe_initial_num_sequences"],
        readout_probe_n_combinations=config["readout_probe_n_combinations"],
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
        readout_probe_Tm_parameters=config["readout_probe_Tm_parameters"],
        readout_probe_Tm_chem_correction_parameters=config["readout_probe_Tm_chem_correction_parameters"],
        readout_probe_Tm_salt_correction_parameters=config["readout_probe_Tm_salt_correction_parameters"],
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
        target_probe_Tm_min=config["target_probe_Tm_min"],
        target_probe_Tm_opt=config["target_probe_Tm_opt"],
        target_probe_Tm_max=config["target_probe_Tm_max"],
        target_probe_homopolymeric_base_n=config["target_probe_homopolymeric_base_n"],
        target_probe_T_secondary_structure=config["target_probe_T_secondary_structure"],
        target_probe_secondary_structures_threshold_deltaG=config[
            "target_probe_secondary_structures_threshold_deltaG"
        ],
        target_probe_GC_weight=config["target_probe_GC_weight"],
        target_probe_Tm_weight=config["target_probe_Tm_weight"],
        target_probe_isoform_weight=config["target_probe_isoform_weight"],
        set_size_opt=config["set_size_opt"],
        set_size_min=config["set_size_min"],
        distance_between_target_probes=config["distance_between_target_probes"],
        n_sets=config["n_sets"],
    )

    codebook, readout_probe_table = pipeline.design_readout_probes(
        n_genes=3,
        files_fasta_reference_database_readout_probe=config["files_fasta_reference_database_readout_probe"],
        readout_probe_length=config["readout_probe_length"],
        readout_probe_base_probabilities=config["readout_probe_base_probabilities"],
        readout_probe_GC_content_min=config["readout_probe_GC_content_min"],
        readout_probe_GC_content_max=config["readout_probe_GC_content_max"],
        readout_probe_homopolymeric_base_n=config["readout_probe_homopolymeric_base_n"],
        readout_probe_set_size=config["readout_probe_set_size"],
        readout_probe_homogeneous_properties_weights=config["readout_probe_homogeneous_properties_weights"],
        n_bits=config["n_bits"],
        min_hamming_dist=config["min_hamming_dist"],
        hamming_weight=config["hamming_weight"],
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
