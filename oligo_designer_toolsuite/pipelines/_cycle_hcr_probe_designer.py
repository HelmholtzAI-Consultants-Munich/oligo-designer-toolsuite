############################################
# imports
############################################

import os
import yaml
import shutil
import logging
import warnings
import itertools

import pandas as pd
import numpy as np

from typing import List, Tuple
from pathlib import Path
from datetime import datetime

from Bio.SeqUtils import Seq
from Bio.SeqUtils import MeltingTemp as mt

from oligo_designer_toolsuite.database import (
    OligoAttributes,
    OligoDatabase,
    ReferenceDatabase,
)
from oligo_designer_toolsuite.oligo_efficiency_filter import (
    AverageSetScoring,
    WeightedIsoformTmScoring,
)
from oligo_designer_toolsuite.oligo_property_filter import (
    SoftMaskedSequenceFilter,
    HardMaskedSequenceFilter,
    HomopolymericRunsFilter,
    GCContentFilter,
    MeltingTemperatureNNFilter,
    SecondaryStructureFilter,
    PropertyFilter,
)
from oligo_designer_toolsuite.oligo_selection import (
    OligosetGeneratorIndependentSet,
    GraphBasedSelectionPolicy,
    GreedySelectionPolicy,
)
from oligo_designer_toolsuite.oligo_specificity_filter import (
    BlastNFilter,
    BlastNSeedregionSiteFilter,
    CrossHybridizationFilter,
    ExactMatchFilter,
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
# CycleHCR Probe Designer
############################################


class CycleHCRProbeDesigner:

    def __init__(self, write_intermediate_steps: bool, dir_output: str, n_jobs: int) -> None:
        """Constructor for the CycleHCRProbeDesigner class."""

        ##### create the output folder #####
        self.dir_output = os.path.abspath(dir_output)
        Path(self.dir_output).mkdir(parents=True, exist_ok=True)

        ##### setup logger #####
        timestamp = datetime.now()
        file_logger = os.path.join(
            self.dir_output,
            f"log_cyclehcr_probe_designer_{timestamp.year}-{timestamp.month}-{timestamp.day}-{timestamp.hour}-{timestamp.minute}.txt",
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
            "strand": "plus",
            "word_size": 10,
            "dust": "no",
            "soft_masking": "false",
            "max_target_seqs": 10,
            "max_hsps": 1000,
        },
        target_probe_specificity_blastn_hit_parameters: dict = {"coverage": 50},
        target_probe_cross_hybridization_blastn_search_parameters: dict = {
            "perc_identity": 80,
            "strand": "minus",
            "word_size": 7,
            "dust": "no",
            "soft_masking": "false",
            "max_target_seqs": 10,
        },
        target_probe_cross_hybridization_blastn_hit_parameters: dict = {"coverage": 50},
        target_probe_Tm_parameters: dict = {
            "check": True,
            "strict": True,
            "c_seq": None,
            "shift": 0,
            "nn_table": "DNA_NN3",
            "tmm_table": "DNA_TMM1",
            "imm_table": "DNA_IMM1",
            "de_table": "DNA_DE1",
            "dnac1": 25,
            "dnac2": 25,
            "selfcomp": False,
            "saltcorr": 0,
            "Na": 50,
            "K": 0,
            "Tris": 0,
            "Mg": 0,
            "dNTPs": 0,
        },
        target_probe_Tm_chem_correction_parameters: dict = None,
        target_probe_Tm_salt_correction_parameters: dict = None,
        max_graph_size: int = 5000,
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

        ### Parameters for the Oligo set selection
        self.max_graph_size = max_graph_size
        self.heuristic = heuristic
        self.n_attempts = n_attempts
        self.heuristic_n_attempts = heuristic_n_attempts

    def design_target_probes(
        self,
        files_fasta_target_probe_database: list[str],
        files_fasta_reference_database_targe_probe: List[str],
        gene_ids: list = None,
        target_probe_isoform_consensus: float = 0,
        target_probe_L_probe_sequence_length: int = 45,
        target_probe_gap_sequence_length: int = 2,
        target_probe_R_probe_sequence_length: int = 45,
        target_probe_GC_content_min: float = 43,
        target_probe_GC_content_max: float = 63,
        target_probe_Tm_min: float = 66,
        target_probe_Tm_max: float = 76,
        target_probe_homopolymeric_base_n: dict = {"A": 5, "T": 5, "C": 5, "G": 5},
        target_probe_T_secondary_structure: float = 76,
        target_probe_secondary_structures_threshold_deltaG: float = 0,
        target_probe_junction_region_size: int = 13,
        target_probe_Tm_weight: float = 1,
        target_probe_isoform_weight: float = 2,
        set_size_opt: int = 50,
        set_size_min: int = 50,
        distance_between_target_probes: int = 0,
        n_sets: int = 100,
    ) -> OligoDatabase:

        target_probe_designer = TargetProbeDesigner(self.dir_output, self.n_jobs)

        oligo_database = target_probe_designer.create_oligo_database(
            gene_ids=gene_ids,
            target_probe_L_probe_sequence_length=target_probe_L_probe_sequence_length,
            target_probe_gap_sequence_length=target_probe_gap_sequence_length,
            target_probe_R_probe_sequence_length=target_probe_R_probe_sequence_length,
            files_fasta_oligo_database=files_fasta_target_probe_database,
            min_oligos_per_gene=set_size_min,
            isoform_consensus=target_probe_isoform_consensus,
        )
        check_content_oligo_database(oligo_database)

        if self.write_intermediate_steps:
            dir_database = oligo_database.save_database(name_database="1_db_target_probes_initial")
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
            dir_database = oligo_database.save_database(name_database="2_db_target_probes_property_filter")
            print(f"Saved target probe database for step 2 (Property Filters) in directory {dir_database}")

        oligo_database = target_probe_designer.filter_by_specificity(
            oligo_database=oligo_database,
            files_fasta_reference_database=files_fasta_reference_database_targe_probe,
            junction_region_size=target_probe_junction_region_size,
            junction_site=target_probe_L_probe_sequence_length + target_probe_gap_sequence_length // 2,
            specificity_blastn_search_parameters=self.target_probe_specificity_blastn_search_parameters,
            specificity_blastn_hit_parameters=self.target_probe_specificity_blastn_hit_parameters,
            cross_hybridization_blastn_search_parameters=self.target_probe_cross_hybridization_blastn_search_parameters,
            cross_hybridization_blastn_hit_parameters=self.target_probe_cross_hybridization_blastn_hit_parameters,
        )
        check_content_oligo_database(oligo_database)

        if self.write_intermediate_steps:
            dir_database = oligo_database.save_database(name_database="3_db_target_probes_specificity_filter")
            print(f"Saved target probe database for step 3 (Specificity Filters) in directory {dir_database}")

        oligo_database = target_probe_designer.create_oligo_sets(
            oligo_database=oligo_database,
            isoform_weight=target_probe_isoform_weight,
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
            heuristic=self.heuristic,
            heuristic_n_attempts=self.heuristic_n_attempts,
        )
        check_content_oligo_database(oligo_database)

        if self.write_intermediate_steps:
            dir_database = oligo_database.save_database(name_database="4_db_target_probes_sets")
            dir_oligosets = oligo_database.write_oligosets_to_table()
            print(
                f"Saved target probe database for step 4 (Specificity Filters) in directory {dir_database} and sets table in directory {dir_oligosets}"
            )

        return oligo_database

    def design_readout_probes(
        self,
        n_regions: int,
        file_readout_probe_table: str,
        file_codebook: str,
    ) -> Tuple[pd.DataFrame, pd.DataFrame]:

        readout_probe_designer = ReadoutProbeDesigner(
            dir_output=self.dir_output,
            n_jobs=self.n_jobs,
        )
        if file_readout_probe_table:
            readout_probe_table, n_channels, n_readout_probes_LR = (
                readout_probe_designer.load_readout_probe_table(
                    file_readout_probe_table=file_readout_probe_table
                )
            )
            logging.info(
                f"Loaded readout probes table from file and retrieved {n_channels} channels and {n_readout_probes_LR} L and R readout probes."
            )
        else:
            raise ValueError("Generation of readout probe table not implemented!")

        if file_codebook:
            raise ValueError("Loading of codeboook not implemented!")
        else:
            codebook = readout_probe_designer.generate_codebook(
                n_regions=n_regions,
                n_channels=n_channels,
                n_readout_probes_LR=n_readout_probes_LR,
            )

        return codebook, readout_probe_table

    def design_encoding_probe(
        self,
        target_probe_database: OligoDatabase,
        codebook: pd.DataFrame,
        readout_probe_table: pd.DataFrame,
        linker_sequence: str,
    ) -> OligoDatabase:

        region_ids = list(target_probe_database.database.keys())

        codebook.index = region_ids + [
            f"unassigned_barcode_{i+1}" for i in range(len(codebook.index) - len(region_ids))
        ]

        codebook.to_csv(os.path.join(self.dir_output, "codebook.tsv"), sep="\t")
        readout_probe_table.to_csv(os.path.join(self.dir_output, "readout_probes.tsv"), sep="\t")

        for region_id in region_ids:
            barcode = codebook.loc[region_id]
            bits = barcode[barcode == 1].index
            readout_probe_sequences = readout_probe_table.loc[bits, "readout_probe_sequence"]

            new_probe_attributes_encoding_probe = {}

            for probe_id in target_probe_database.database[region_id].keys():

                sequence_readout_probe_L = readout_probe_sequences.iloc[0]
                sequence_readout_probe_R = readout_probe_sequences.iloc[1]

                new_probe_attributes_encoding_probe[probe_id] = {
                    "barcode": barcode,
                    "sequence_target": target_probe_database.get_oligo_attribute_value(
                        attribute="target", region_id=region_id, oligo_id=probe_id, flatten=True
                    ),
                    "sequence_target_probe_L": target_probe_database.get_oligo_attribute_value(
                        attribute="oligo_pair_L", region_id=region_id, oligo_id=probe_id, flatten=True
                    ),
                    "sequence_target_probe_R": target_probe_database.get_oligo_attribute_value(
                        attribute="oligo_pair_R", region_id=region_id, oligo_id=probe_id, flatten=True
                    ),
                    "sequence_readout_probe_L": sequence_readout_probe_L,
                    "sequence_readout_probe_R": sequence_readout_probe_R,
                    "sequence_encoding_probe_L": target_probe_database.get_oligo_attribute_value(
                        attribute="oligo_pair_L", region_id=region_id, oligo_id=probe_id, flatten=True
                    )
                    + linker_sequence
                    + str(Seq(sequence_readout_probe_L).reverse_complement()),
                    "sequence_encoding_probe_R": str(Seq(sequence_readout_probe_R).reverse_complement())
                    + linker_sequence
                    + target_probe_database.get_oligo_attribute_value(
                        attribute="oligo_pair_R", region_id=region_id, oligo_id=probe_id, flatten=True
                    ),
                }

            target_probe_database.update_oligo_attributes(new_probe_attributes_encoding_probe)

        return target_probe_database

    def design_primers(
        self,
        forward_primer_sequence: str,
        reverse_primer_sequence: str,
    ):

        primer_designer = PrimerDesigner(
            dir_output=self.dir_output,
            n_jobs=self.n_jobs,
        )

        if forward_primer_sequence:
            forward_primer_sequence = forward_primer_sequence
        else:
            # generate forward primers
            raise ValueError("Forward primer generation not implemented!")

        if reverse_primer_sequence:
            reverse_primer_sequence = reverse_primer_sequence
        else:
            # generate reverse primers
            raise ValueError("Reverse primer generation not implemented!")

        return reverse_primer_sequence, forward_primer_sequence

    def generate_output(
        self,
        encoding_probe_database: OligoDatabase,
        reverse_primer_sequence: str,
        forward_primer_sequence: str,
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
            "sequence_cyclehcr_probe_L",
            "sequence_cyclehcr_probe_R",
            "sequence_encoding_probe_L",
            "sequence_encoding_probe_L_rc",
            "sequence_encoding_probe_R",
            "sequence_encoding_probe_R_rc",
            "sequence_readout_probe_L",
            "sequence_readout_probe_R",
            "sequence_forward_primer",
            "sequence_reverse_primer",
            "sequence_target",
            "TmNN_oligo_pair_L",
            "TmNN_oligo_pair_R",
        ],
    ) -> None:

        new_probe_attributes_primer = {}

        for region_id in encoding_probe_database.database.keys():
            for probe_id in encoding_probe_database.database[region_id].keys():
                new_probe_attributes_primer[probe_id] = {
                    "sequence_reverse_primer": reverse_primer_sequence,
                    "sequence_forward_primer": forward_primer_sequence,
                    "sequence_cyclehcr_probe_L": forward_primer_sequence
                    + encoding_probe_database.get_oligo_attribute_value(
                        attribute="sequence_encoding_probe_L",
                        region_id=region_id,
                        oligo_id=probe_id,
                        flatten=True,
                    )
                    + reverse_primer_sequence,
                    "sequence_cyclehcr_probe_R": forward_primer_sequence
                    + encoding_probe_database.get_oligo_attribute_value(
                        attribute="sequence_encoding_probe_R",
                        region_id=region_id,
                        oligo_id=probe_id,
                        flatten=True,
                    )
                    + reverse_primer_sequence,
                    "sequence_encoding_probe_L_rc": str(
                        Seq(
                            encoding_probe_database.get_oligo_attribute_value(
                                attribute="sequence_encoding_probe_L",
                                region_id=region_id,
                                oligo_id=probe_id,
                                flatten=True,
                            )
                        ).reverse_complement()
                    ),
                    "sequence_encoding_probe_R_rc": str(
                        Seq(
                            encoding_probe_database.get_oligo_attribute_value(
                                attribute="sequence_encoding_probe_R",
                                region_id=region_id,
                                oligo_id=probe_id,
                                flatten=True,
                            )
                        ).reverse_complement()
                    ),
                }
        encoding_probe_database.update_oligo_attributes(new_probe_attributes_primer)

        encoding_probe_database = self.oligo_attributes_calculator.calculate_TmNN(
            oligo_database=encoding_probe_database,
            Tm_parameters=self.target_probe_Tm_parameters,
            Tm_chem_correction_parameters=self.target_probe_Tm_chem_correction_parameters,
            Tm_salt_correction_parameters=self.target_probe_Tm_salt_correction_parameters,
            sequence_type="oligo_pair_L",
        )
        encoding_probe_database = self.oligo_attributes_calculator.calculate_TmNN(
            oligo_database=encoding_probe_database,
            Tm_parameters=self.target_probe_Tm_parameters,
            Tm_chem_correction_parameters=self.target_probe_Tm_chem_correction_parameters,
            Tm_salt_correction_parameters=self.target_probe_Tm_salt_correction_parameters,
            sequence_type="oligo_pair_R",
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
            filename="cyclehcr_probes",
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
                        "sequence_cyclehcr_probe_L": encoding_probe_database.get_oligo_attribute_value(
                            attribute="sequence_cyclehcr_probe_L",
                            region_id=region_id,
                            oligo_id=oligo_id,
                            flatten=True,
                        ),
                        "sequence_cyclehcr_probe_R": encoding_probe_database.get_oligo_attribute_value(
                            attribute="sequence_cyclehcr_probe_R",
                            region_id=region_id,
                            oligo_id=oligo_id,
                            flatten=True,
                        ),
                        "sequence_readout_probe_L": encoding_probe_database.get_oligo_attribute_value(
                            attribute="sequence_readout_probe_L",
                            region_id=region_id,
                            oligo_id=oligo_id,
                            flatten=True,
                        ),
                        "sequence_readout_probe_R": encoding_probe_database.get_oligo_attribute_value(
                            attribute="sequence_readout_probe_R",
                            region_id=region_id,
                            oligo_id=oligo_id,
                            flatten=True,
                        ),
                    }

        with open(os.path.join(self.dir_output, "cyclehcr_probes_order.yml"), "w") as outfile:
            yaml.dump(yaml_dict_order, outfile, default_flow_style=False, sort_keys=False)

        logging.info("--------------END PIPELINE--------------")


############################################
# CycleHCR Target Probe Designer
############################################


class TargetProbeDesigner:

    def __init__(self, dir_output: str, n_jobs: int) -> None:
        """Constructor for the TargetProbeDesigner class."""

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
        target_probe_L_probe_sequence_length: int,
        target_probe_gap_sequence_length: int,
        target_probe_R_probe_sequence_length: int,
        files_fasta_oligo_database: list[str],
        min_oligos_per_gene: int,
        isoform_consensus: float,
    ) -> OligoDatabase:

        ##### creating the oligo sequences #####
        oligo_length = (
            target_probe_L_probe_sequence_length
            + target_probe_gap_sequence_length
            + target_probe_R_probe_sequence_length
        )
        oligo_sequences = OligoSequenceGenerator(dir_output=self.dir_output)
        oligo_fasta_file = oligo_sequences.create_sequences_sliding_window(
            files_fasta_in=files_fasta_oligo_database,
            length_interval_sequences=(oligo_length, oligo_length),
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
        ##### calculate probe pairs
        split_start_end = [
            (0, target_probe_L_probe_sequence_length),
            (
                target_probe_L_probe_sequence_length,
                target_probe_L_probe_sequence_length + target_probe_gap_sequence_length,
            ),
            (
                target_probe_L_probe_sequence_length + target_probe_gap_sequence_length,
                target_probe_L_probe_sequence_length
                + target_probe_gap_sequence_length
                + target_probe_R_probe_sequence_length,
            ),
        ]
        oligo_database = self.oligo_attributes_calculator.calculate_split_sequence(
            oligo_database=oligo_database,
            split_start_end=split_start_end,
            split_names=["oligo_pair_L", "spacer", "oligo_pair_R"],
            sequence_type="target",
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
    ) -> OligoDatabase:

        # define the filters
        # soft_masked_sequences = SoftMaskedSequenceFilter()
        hard_masked_sequences = HardMaskedSequenceFilter()
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
            # soft_masked_sequences,
            hard_masked_sequences,
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
            sequence_type="oligo_pair_L",
            n_jobs=self.n_jobs,
        )
        oligo_database = property_filter.apply(
            oligo_database=oligo_database,
            sequence_type="oligo_pair_R",
            n_jobs=self.n_jobs,
        )

        return oligo_database

    @pipeline_step_basic(step_name="Target Probe Generation - Specificity Filters")
    def filter_by_specificity(
        self,
        oligo_database: OligoDatabase,
        files_fasta_reference_database: List[str],
        junction_region_size: int,
        junction_site: int,
        specificity_blastn_search_parameters: dict,
        specificity_blastn_hit_parameters: dict,
        cross_hybridization_blastn_search_parameters: dict,
        cross_hybridization_blastn_hit_parameters: dict,
    ) -> OligoDatabase:

        ##### define reference database #####
        reference_database = ReferenceDatabase(
            database_name=self.subdir_db_reference, dir_output=self.dir_output
        )
        reference_database.load_database_from_file(
            files=files_fasta_reference_database, file_type="fasta", database_overwrite=False
        )

        ##### define specificity filters #####
        exact_matches = ExactMatchFilter(policy=RemoveAllPolicy(), filter_name="oligo_exact_match")

        if junction_region_size > 0:
            oligo_ids = oligo_database.get_oligoid_list()
            oligo_database.update_oligo_attributes(
                new_oligo_attribute={oligo_id: {"junction_site": junction_site} for oligo_id in oligo_ids}
            )
            specificity = BlastNSeedregionSiteFilter(
                seedregion_size=junction_region_size,
                seedregion_site_name="junction_site",
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
        specificity.set_reference_database(reference_database=reference_database)

        ##### run specificity filters #####
        specificity_filter = SpecificityFilter(filters=[exact_matches, specificity])
        oligo_database = specificity_filter.apply(
            oligo_database=oligo_database,
            sequence_type="target",
            n_jobs=self.n_jobs,
        )

        ##### define cross hybridization filter #####
        cross_hybridization_aligner_oligo_pair_L = BlastNFilter(
            remove_hits=True,
            search_parameters=cross_hybridization_blastn_search_parameters,
            hit_parameters=cross_hybridization_blastn_hit_parameters,
            filter_name="blastn_crosshybridization",
            dir_output=self.dir_output,
        )
        cross_hybridization_oligo_pair_L = CrossHybridizationFilter(
            policy=RemoveByLargerRegionPolicy(),
            alignment_method=cross_hybridization_aligner_oligo_pair_L,
            sequence_type_reference="oligo_pair_L",
            filter_name="blastn_crosshybridization",
            dir_output=self.dir_output,
        )
        cross_hybridization_aligner_oligo_pair_R = BlastNFilter(
            remove_hits=True,
            search_parameters=cross_hybridization_blastn_search_parameters,
            hit_parameters=cross_hybridization_blastn_hit_parameters,
            filter_name="blastn_crosshybridization",
            dir_output=self.dir_output,
        )
        cross_hybridization_oligo_pair_R = CrossHybridizationFilter(
            policy=RemoveByLargerRegionPolicy(),
            alignment_method=cross_hybridization_aligner_oligo_pair_R,
            sequence_type_reference="oligo_pair_R",
            filter_name="blastn_crosshybridization",
            dir_output=self.dir_output,
        )

        ##### run cross hybridization filter #####
        specificity_filter = SpecificityFilter(
            filters=[cross_hybridization_oligo_pair_L, cross_hybridization_oligo_pair_R]
        )
        oligo_database = specificity_filter.apply(
            oligo_database=oligo_database,
            sequence_type="oligo_pair_L",
            n_jobs=self.n_jobs,
        )
        oligo_database = specificity_filter.apply(
            oligo_database=oligo_database,
            sequence_type="oligo_pair_R",
            n_jobs=self.n_jobs,
        )

        ##### remove all directories of intermediate steps #####
        for directory in [
            cross_hybridization_aligner_oligo_pair_L.dir_output,
            cross_hybridization_aligner_oligo_pair_R.dir_output,
            cross_hybridization_oligo_pair_L.dir_output,
            cross_hybridization_oligo_pair_R.dir_output,
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
    ) -> Tuple[OligoDatabase, str, str]:

        # the higher the score the better, because we want to have on average oligos with high melting temperatures
        set_scoring = AverageSetScoring(ascending=False)

        # We change the processing dependent on the required number of probes in the probe sets
        # For small sets, we don't pre-filter and find the initial set by iterating
        # through all possible generated sets, which is faster than the max clique approximation.
        if set_size_opt < 10:
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
        elif 10 < set_size_opt < 30:
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
        else:
            pre_filter = True
            selection_policy = GreedySelectionPolicy(
                set_scoring=set_scoring,
                score_criteria=set_scoring.score_1,
                pre_filter=pre_filter,
                penalty=0.01,
                n_attempts=n_attempts,
            )
            base_log_parameters({"pre_filter": pre_filter, "selection_policy": "Greedy"})

        oligos_scoring = WeightedIsoformTmScoring(
            Tm_content_opt=Tm_max,
            Tm_weight=Tm_weight,
            isoform_weight=isoform_weight,
            Tm_parameters=Tm_parameters,
            Tm_chem_correction_parameters=Tm_chem_correction_parameters,
            Tm_salt_correction_parameters=Tm_salt_correction_parameters,
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
            sequence_type="target",
            set_size_opt=set_size_opt,
            set_size_min=set_size_min,
            n_sets=n_sets,
            n_jobs=self.n_jobs,
        )

        return oligo_database


############################################
# CycleHCR Readout Probe Designer
############################################


class ReadoutProbeDesigner:

    def __init__(
        self,
        dir_output: str,
        n_jobs: int,
    ) -> None:
        """Constructor for the ReadoutProbeDesigner class."""

        ##### create the output folder #####
        self.dir_output = os.path.abspath(dir_output)

        self.n_jobs = n_jobs

    def generate_codebook(self, n_regions: int, n_channels: int, n_readout_probes_LR: int) -> pd.DataFrame:

        def _generate_barcode(combination: set, codebook_size: int) -> list:
            index1 = ((n_channels * 2) * combination[0]) + (2 * combination[2])
            index2 = ((n_channels * 2) * combination[1]) + (2 * combination[2]) + 1
            barcode = np.zeros(codebook_size, dtype=np.int8)
            barcode[[index1, index2]] = 1
            return barcode

        codebook = []
        codebook_size = n_channels * n_readout_probes_LR * 2

        combinations = list(
            itertools.product(
                list(range(n_readout_probes_LR)), list(range(n_readout_probes_LR)), list(range(n_channels))
            )
        )
        combinations = sorted(combinations, key=lambda t: (0 if t[0] == t[1] else 1, t[1]))
        codebook_size_max = len(combinations)

        if codebook_size_max < (2 * n_regions):
            raise ValueError(
                f"The number of valid barcodes ({codebook_size_max}) is lower than the required number of readout probes ({2 * n_regions}) for {n_regions} regions. Consider increasing the number of L/R readout probes."
            )

        for combination in combinations[:n_regions]:
            barcode = _generate_barcode(
                combination=combination,
                codebook_size=codebook_size,
            )
            codebook.append(barcode)

        codebook = pd.DataFrame(codebook, columns=[f"bit_{i+1}" for i in range(codebook_size)])
        codebook

        return codebook

    def load_readout_probe_table(self, file_readout_probe_table: str):

        required_cols = ["channel", "readout_probe_id", "readout_probe_sequence", "L/R"]

        readout_probe_table = pd.read_csv(file_readout_probe_table, sep=None, engine="python")

        # Check if all required columns exist in readout_probe_table
        cols = set(readout_probe_table.columns)
        if not set(required_cols).issubset(cols):
            missing = set(required_cols) - cols
            raise ValueError(f"Missing columns: {missing}")

        if "bit" not in readout_probe_table.columns:
            readout_probe_table = readout_probe_table.sort_values(by=["readout_probe_id", "channel"])
            readout_probe_table.reset_index(inplace=True, drop=True)
            readout_probe_table["bit"] = "bit_" + (readout_probe_table.index + 1).astype(str)

        readout_probe_table.set_index("bit", inplace=True)
        readout_probe_table = readout_probe_table[required_cols]

        n_channels = len(readout_probe_table["channel"].unique())
        n_readout_probes_R = readout_probe_table["L/R"].value_counts()["R"]
        n_readout_probes_L = readout_probe_table["L/R"].value_counts()["L"]
        n_readout_probes_LR = int(min([n_readout_probes_R, n_readout_probes_L]) / n_channels)

        return readout_probe_table, n_channels, n_readout_probes_LR


############################################
# CycleHCR Primer Designer
############################################


class PrimerDesigner:

    def __init__(
        self,
        dir_output: str,
        n_jobs: int,
    ) -> None:
        """Constructor for the PrimerDesigner class."""

        ##### create the output folder #####
        self.dir_output = os.path.abspath(dir_output)
        self.n_jobs = n_jobs


############################################
# CycleHCR Probe Designer Pipeline
############################################


def main():
    """
    Main function for running the CycleHCRProbeDesigner pipeline. This function reads the configuration file,
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
    pipeline = CycleHCRProbeDesigner(
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
        max_graph_size=config["max_graph_size"],
        n_attempts=config["n_attempts"],
        heuristic=config["heuristic"],
        heuristic_n_attempts=config["heuristic_n_attempts"],
    )

    ##### design probes #####
    target_probe_database = pipeline.design_target_probes(
        files_fasta_target_probe_database=config["files_fasta_target_probe_database"],
        files_fasta_reference_database_targe_probe=config["files_fasta_reference_database_targe_probe"],
        gene_ids=gene_ids,
        target_probe_isoform_consensus=config["target_probe_isoform_consensus"],
        target_probe_L_probe_sequence_length=config["target_probe_L_probe_sequence_length"],
        target_probe_gap_sequence_length=config["target_probe_gap_sequence_length"],
        target_probe_R_probe_sequence_length=config["target_probe_R_probe_sequence_length"],
        target_probe_GC_content_min=config["target_probe_GC_content_min"],
        target_probe_GC_content_max=config["target_probe_GC_content_max"],
        target_probe_Tm_min=config["target_probe_Tm_min"],
        target_probe_Tm_max=config["target_probe_Tm_max"],
        target_probe_homopolymeric_base_n=config["target_probe_homopolymeric_base_n"],
        target_probe_T_secondary_structure=config["target_probe_T_secondary_structure"],
        target_probe_secondary_structures_threshold_deltaG=config[
            "target_probe_secondary_structures_threshold_deltaG"
        ],
        target_probe_junction_region_size=config["target_probe_junction_region_size"],
        target_probe_Tm_weight=config["target_probe_Tm_weight"],
        target_probe_isoform_weight=config["target_probe_isoform_weight"],
        set_size_opt=config["set_size_opt"],
        set_size_min=config["set_size_min"],
        distance_between_target_probes=config["distance_between_target_probes"],
        n_sets=config["n_sets"],
    )

    codebook, readout_probe_table = pipeline.design_readout_probes(
        n_regions=len(target_probe_database.database),
        file_readout_probe_table=config["file_readout_probe_table"],
        file_codebook=config["file_codebook"],
    )

    encoding_probe_database = pipeline.design_encoding_probe(
        target_probe_database=target_probe_database,
        codebook=codebook,
        readout_probe_table=readout_probe_table,
        linker_sequence=config["linker_sequence"],
    )

    reverse_primer_sequence, forward_primer_sequence = pipeline.design_primers(
        forward_primer_sequence=config["forward_primer_sequence"],
        reverse_primer_sequence=config["reverse_primer_sequence"],
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
