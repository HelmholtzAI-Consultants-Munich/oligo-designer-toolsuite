############################################
# imports
############################################

import itertools
import logging
import os
import random
import shutil
import warnings
from datetime import datetime
from pathlib import Path
from typing import List, Tuple

import yaml
from Bio.SeqUtils import MeltingTemp as mt
from joblib import Parallel, delayed
from joblib_progress import joblib_progress

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
    DetectionOligoFilter,
    GCContentFilter,
    HardMaskedSequenceFilter,
    HomopolymericRunsFilter,
    MeltingTemperatureNNFilter,
    PropertyFilter,
    SoftMaskedSequenceFilter,
)
from oligo_designer_toolsuite.oligo_selection import (
    GraphBasedSelectionPolicy,
    OligosetGeneratorIndependentSet,
)
from oligo_designer_toolsuite.oligo_specificity_filter import (
    BlastNFilter,
    BlastNSeedregionLigationsiteFilter,
    CrossHybridizationFilter,
    ExactMatchFilter,
    RemoveByLargerRegionPolicy,
    SpecificityFilter,
)
from oligo_designer_toolsuite.pipelines._utils import (
    base_parser,
    pipeline_step_basic,
    check_content_oligo_database,
)
from oligo_designer_toolsuite.sequence_generator import OligoSequenceGenerator

############################################
# SCRINSHOT Probe Designer Functions
############################################


class ScrinshotProbeDesigner:
    def __init__(self, write_intermediate_steps: bool, dir_output: str, n_jobs: int) -> None:
        """Constructor for the ScrinshotProbeDesigner class."""

        ##### create the output folder #####
        self.dir_output = os.path.abspath(dir_output)
        Path(self.dir_output).mkdir(parents=True, exist_ok=True)

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
        target_probe_specificity_blastn_hit_parameters: dict = {"coverage": 50},
        target_probe_cross_hybridization_blastn_search_parameters: dict = {
            "perc_identity": 80,
            "strand": "minus",
            "word_size": 10,
            "dust": "no",
            "soft_masking": "false",
            "max_target_seqs": 10,
        },
        target_probe_cross_hybridization_blastn_hit_parameters: dict = {"coverage": 80},
        max_graph_size: int = 5000,
        pre_filter: bool = False,
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
            "Na": 39,
            "K": 75,
            "Tris": 20,
            "Mg": 10,
            "dNTPs": 0,
        },
        target_probe_Tm_chem_correction_param: dict = {
            "DMSO": 0,
            "fmd": 20,
            "DMSOfactor": 0.75,
            "fmdfactor": 0.65,
            "fmdmethod": 1,
            "GC": None,
        },
        target_probe_Tm_salt_correction_param: dict = None,
        detection_oligo_Tm_parameters: dict = {
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
            "Na": 39,
            "K": 0,
            "Tris": 0,
            "Mg": 0,
            "dNTPs": 0,
        },
        detection_oligo_Tm_chem_correction_param: dict = {
            "DMSO": 0,
            "fmd": 30,
            "DMSOfactor": 0.75,
            "fmdfactor": 0.65,
            "fmdmethod": 1,
            "GC": None,
        },
        detection_oligo_Tm_salt_correction_param: dict = None,
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

        ### Parameters for the Oligo set selection
        self.max_graph_size = max_graph_size
        self.pre_filter = pre_filter
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

        detection_oligo_Tm_parameters["nn_table"] = getattr(mt, detection_oligo_Tm_parameters["nn_table"])
        detection_oligo_Tm_parameters["tmm_table"] = getattr(mt, detection_oligo_Tm_parameters["tmm_table"])
        detection_oligo_Tm_parameters["imm_table"] = getattr(mt, detection_oligo_Tm_parameters["imm_table"])
        detection_oligo_Tm_parameters["de_table"] = getattr(mt, detection_oligo_Tm_parameters["de_table"])

        ## target probe
        self.target_probe_Tm_parameters = target_probe_Tm_parameters
        self.target_probe_Tm_chem_correction_param = target_probe_Tm_chem_correction_param
        self.target_probe_Tm_salt_correction_param = target_probe_Tm_salt_correction_param

        ## detection oligo
        self.detection_oligo_Tm_parameters = detection_oligo_Tm_parameters
        self.detection_oligo_Tm_chem_correction_param = detection_oligo_Tm_chem_correction_param
        self.detection_oligo_Tm_salt_correction_param = detection_oligo_Tm_salt_correction_param

    def design_target_probes(
        self,
        files_fasta_target_probe_database: list,
        files_fasta_reference_database_targe_probe: list,
        gene_ids: list = None,
        target_probe_length_min: int = 40,
        target_probe_length_max: int = 45,
        target_probe_isoform_consensus: float = 50,
        target_probe_isoform_weight: float = 2,
        target_probe_GC_content_min: float = 40,
        target_probe_GC_content_opt: float = 50,
        target_probe_GC_content_max: float = 60,
        target_probe_GC_weight: float = 1,
        target_probe_Tm_min: float = 65,
        target_probe_Tm_opt: float = 70,
        target_probe_Tm_max: float = 75,
        target_probe_Tm_weight: float = 1,
        target_probe_homopolymeric_base_n: dict = {"A": 5, "T": 5, "C": 5, "G": 5},
        detection_oligo_min_thymines: int = 2,
        detection_oligo_length_min: int = 15,
        detection_oligo_length_max: int = 40,
        target_probe_padlock_arm_length_min: int = 10,
        target_probe_padlock_arm_Tm_dif_max: float = 2,
        target_probe_padlock_arm_Tm_min: float = 50,
        target_probe_padlock_arm_Tm_max: float = 60,
        target_probe_ligation_region_size: int = 5,
        set_size_min: int = 3,
        set_size_opt: int = 5,
        distance_between_target_probes: int = 0,
        n_sets: int = 100,
    ) -> OligoDatabase:

        target_probe_designer = TargetProbeDesigner(self.dir_output, self.n_jobs)

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
            dir_database = oligo_database.save_database(dir_database="1_db_probes_initial")
            print(f"Saved probe database for step 1 (Create Database) in directory {dir_database}")

        oligo_database = target_probe_designer.filter_by_property(
            oligo_database=oligo_database,
            GC_content_min=target_probe_GC_content_min,
            GC_content_max=target_probe_GC_content_max,
            Tm_min=target_probe_Tm_min,
            Tm_max=target_probe_Tm_max,
            detect_oligo_length_min=detection_oligo_length_min,
            detect_oligo_length_max=detection_oligo_length_max,
            min_thymines=detection_oligo_min_thymines,
            arm_length_min=target_probe_padlock_arm_length_min,
            arm_Tm_dif_max=target_probe_padlock_arm_Tm_dif_max,
            arm_Tm_min=target_probe_padlock_arm_Tm_min,
            arm_Tm_max=target_probe_padlock_arm_Tm_max,
            homopolymeric_base_n=target_probe_homopolymeric_base_n,
            Tm_parameters=self.target_probe_Tm_parameters,
            Tm_chem_correction_param=self.target_probe_Tm_chem_correction_param,
            Tm_salt_correction_param=self.target_probe_Tm_salt_correction_param,
        )
        check_content_oligo_database(oligo_database)

        if self.write_intermediate_steps:
            dir_database = oligo_database.save_database(dir_database="2_db_probes_property_filter")
            print(f"Saved probe database for step 2 (Property Filters) in directory {dir_database}")

        oligo_database = target_probe_designer.filter_by_specificity(
            oligo_database=oligo_database,
            files_fasta_reference_database=files_fasta_reference_database_targe_probe,
            specificity_blastn_search_parameters=self.target_probe_specificity_blastn_search_parameters,
            specificity_blastn_hit_parameters=self.target_probe_specificity_blastn_hit_parameters,
            cross_hybridization_blastn_search_parameters=self.target_probe_cross_hybridization_blastn_search_parameters,
            cross_hybridization_blastn_hit_parameters=self.target_probe_cross_hybridization_blastn_hit_parameters,
            ligation_region_size=target_probe_ligation_region_size,
            arm_length_min=target_probe_padlock_arm_length_min,
            arm_Tm_dif_max=target_probe_padlock_arm_Tm_dif_max,
            arm_Tm_min=target_probe_padlock_arm_Tm_min,
            arm_Tm_max=target_probe_padlock_arm_Tm_max,
            Tm_parameters=self.target_probe_Tm_parameters,
            Tm_chem_correction_param=self.target_probe_Tm_chem_correction_param,
            Tm_salt_correction_param=self.target_probe_Tm_salt_correction_param,
        )
        check_content_oligo_database(oligo_database)

        if self.write_intermediate_steps:
            dir_database = oligo_database.save_database(dir_database="3_db_probes_specificity_filter")
            print(f"Saved probe database for step 3 (Specificity Filters) in directory {dir_database}")

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
            Tm_chem_correction_param=self.target_probe_Tm_chem_correction_param,
            Tm_salt_correction_param=self.target_probe_Tm_salt_correction_param,
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
            dir_database = oligo_database.save_database(dir_database="4_db_probes_probesets")
            dir_probesets = oligo_database.write_oligosets_to_table()
            print(
                f"Saved probe database for step 4 (Specificity Filters) in directory {dir_database} and probeset table in directory {dir_probesets}"
            )

        return oligo_database

    def design_detection_oligos(
        self,
        oligo_database: OligoDatabase,
        detection_oligo_length_min: int = 15,
        detection_oligo_length_max: int = 40,
        detection_oligo_min_thymines: int = 2,
        detection_oligo_U_distance: int = 5,
        detection_oligo_Tm_opt: float = 56,
    ) -> OligoDatabase:
        detection_oligo_designer = DetectionOligoDesigner(self.n_jobs)
        oligo_database = detection_oligo_designer.create_detection_oligos(
            oligo_database=oligo_database,
            oligo_length_min=detection_oligo_length_min,
            oligo_length_max=detection_oligo_length_max,
            min_thymines=detection_oligo_min_thymines,
            U_distance=detection_oligo_U_distance,
            Tm_opt=detection_oligo_Tm_opt,
            Tm_parameters=self.detection_oligo_Tm_parameters,
            Tm_chem_correction_parameters=self.detection_oligo_Tm_chem_correction_param,
            Tm_salt_correction_parameters=self.detection_oligo_Tm_salt_correction_param,
        )

        return oligo_database

    def design_padlock_backbone(
        self,
        oligo_database: OligoDatabase,
    ) -> OligoDatabase:

        def _get_barcode(number_regions: int, barcode_length: int, seed: int, choices: list) -> list:

            while len(choices) ** barcode_length < number_regions:
                barcode_length += 1

            barcodes = ["".join(nts) for nts in itertools.product(choices, repeat=barcode_length)]
            random.seed(seed)
            random.shuffle(barcodes)

            return barcodes

        region_ids = list(oligo_database.database.keys())

        barcodes = _get_barcode(len(region_ids), barcode_length=4, seed=0, choices=["A", "C", "T", "G"])

        for region_idx, region_id in enumerate(region_ids):
            oligo_sets_region = oligo_database.oligosets[region_id]
            oligo_sets_oligo_columns = [col for col in oligo_sets_region.columns if col.startswith("oligo_")]

            new_oligo_attributes = {}

            for index in range(len(oligo_sets_region.index)):
                for column in oligo_sets_oligo_columns:
                    oligo_id = str(oligo_sets_region.loc[index, column])
                    barcode = barcodes[region_idx]

                    ligation_site = oligo_database.get_oligo_attribute_value(
                        attribute="ligation_site", region_id=region_id, oligo_id=oligo_id, flatten=True
                    )
                    sequence_oligo = oligo_database.get_oligo_attribute_value(
                        attribute="oligo", region_id=region_id, oligo_id=oligo_id, flatten=True
                    )
                    sequence_padlock_arm1 = sequence_oligo[ligation_site:]
                    sequence_padlock_arm2 = sequence_oligo[:ligation_site]
                    sequence_padlock_accessory1 = "TCCTCTATGATTACTGAC"
                    sequence_padlock_ISS_anchor = "TGCGTCTATTTAGTGGAGCC"
                    sequence_padlock_accessory2 = "CTATCTTCTTT"
                    sequence_padlock_backbone = (
                        sequence_padlock_accessory1
                        + sequence_padlock_ISS_anchor
                        + barcode
                        + sequence_padlock_accessory2
                    )
                    sequence_padlock_probe = (
                        sequence_padlock_arm1 + sequence_padlock_backbone + sequence_padlock_arm2
                    )
                    Tm_arm1 = self.oligo_attributes_calculator._calc_TmNN(
                        sequence=sequence_padlock_arm1,
                        Tm_parameters=self.target_probe_Tm_parameters,
                        Tm_chem_correction_parameters=self.target_probe_Tm_chem_correction_param,
                        Tm_salt_correction_parameters=self.target_probe_Tm_salt_correction_param,
                    )
                    Tm_arm2 = self.oligo_attributes_calculator._calc_TmNN(
                        sequence=sequence_padlock_arm2,
                        Tm_parameters=self.target_probe_Tm_parameters,
                        Tm_chem_correction_parameters=self.target_probe_Tm_chem_correction_param,
                        Tm_salt_correction_parameters=self.target_probe_Tm_salt_correction_param,
                    )

                    new_oligo_attributes[oligo_id] = {
                        "barcode": barcode,
                        "sequence_target": oligo_database.get_oligo_attribute_value(
                            attribute="target", region_id=region_id, oligo_id=oligo_id, flatten=True
                        ),
                        "sequence_target_probe": oligo_database.get_oligo_attribute_value(
                            attribute="oligo", region_id=region_id, oligo_id=oligo_id, flatten=True
                        ),
                        "sequence_padlock_arm1": sequence_padlock_arm1,
                        "sequence_padlock_arm2": sequence_padlock_arm2,
                        "sequence_padlock_accessory1": sequence_padlock_accessory1,
                        "sequence_padlock_ISS_anchor": sequence_padlock_ISS_anchor,
                        "sequence_padlock_accessory2": sequence_padlock_accessory2,
                        "sequence_padlock_backbone": sequence_padlock_backbone,
                        "sequence_padlock_probe": sequence_padlock_probe,
                        "Tm_arm1": Tm_arm1,
                        "Tm_arm2": Tm_arm2,
                        "Tm_diff_arms": round(abs(Tm_arm1 - Tm_arm2), 2),
                    }

            oligo_database.update_oligo_attributes(new_oligo_attributes)

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
            "sequence_padlock_probe",
            "sequence_detection_oligo",
            "sequence_padlock_arm1",
            "sequence_padlock_accessory1",
            "sequence_padlock_ISS_anchor",
            "barcode",
            "sequence_padlock_accessory2",
            "sequence_padlock_arm2",
            "sequence_target",
            "sequence_target_probe",
            "length",
            "ligation_site",
            "Tm_arm1",
            "Tm_arm2",
            "Tm_diff_arms",
            "Tm_detection_oligo",
            "isoform_consensus",
        ],
    ) -> None:
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
            Tm_chem_correction_parameters=self.target_probe_Tm_chem_correction_param,
            Tm_salt_correction_parameters=self.target_probe_Tm_salt_correction_param,
        )
        oligo_database = self.oligo_attributes_calculator.calculate_num_targeted_transcripts(
            oligo_database=oligo_database
        )
        oligo_database = self.oligo_attributes_calculator.calculate_isoform_consensus(
            oligo_database=oligo_database
        )

        oligo_database.write_oligosets_to_yaml(
            attributes=attributes,
            top_n_sets=top_n_sets,
            ascending=True,
            filename="padlock_probes.yml",
        )

        # write a second file that only contains order information
        yaml_dict_order = {}

        for region_id in oligo_database.database.keys():
            yaml_dict_order[region_id] = {}
            oligosets_region = oligo_database.oligosets[region_id]
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
                        "sequence_padlock_probe": oligo_database.get_oligo_attribute_value(
                            attribute="sequence_padlock_probe",
                            region_id=region_id,
                            oligo_id=oligo_id,
                            flatten=True,
                        ),
                        "sequence_detection_oligo": oligo_database.get_oligo_attribute_value(
                            attribute="sequence_detection_oligo",
                            region_id=region_id,
                            oligo_id=oligo_id,
                            flatten=True,
                        ),
                    }

        with open(os.path.join(self.dir_output, "padlock_probes_order.yml"), "w") as outfile:
            yaml.dump(yaml_dict_order, outfile, default_flow_style=False, sort_keys=False)

        logging.info("--------------END PIPELINE--------------")


class TargetProbeDesigner:
    def __init__(self, dir_output: str, n_jobs: int) -> None:
        """Constructor for the ScrinshotProbeDesigner class."""

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
        files_fasta_oligo_database: list[str],
        min_oligos_per_gene: int,
        isoform_consensus: float,
    ) -> Tuple[OligoDatabase, str]:
        ##### creating the probe sequences #####
        oligo_sequences = OligoSequenceGenerator(dir_output=self.dir_output)
        oligo_fasta_file = oligo_sequences.create_sequences_sliding_window(
            files_fasta_in=files_fasta_oligo_database,
            length_interval_sequences=(oligo_length_min, oligo_length_max),
            region_ids=gene_ids,
            n_jobs=self.n_jobs,
        )

        ##### creating the probe database #####
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
            database_overwrite=True,
            sequence_type="target",
            region_ids=gene_ids,
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

    @pipeline_step_basic(step_name="Property Filters")
    def filter_by_property(
        self,
        oligo_database: OligoDatabase,
        GC_content_min: float,
        GC_content_max: float,
        Tm_min: float,
        Tm_max: float,
        detect_oligo_length_min: int,
        detect_oligo_length_max: int,
        min_thymines: int,
        arm_Tm_dif_max: int,
        arm_length_min: int,
        arm_Tm_min: float,
        arm_Tm_max: float,
        homopolymeric_base_n: str,
        Tm_parameters: dict,
        Tm_chem_correction_param: dict,
        Tm_salt_correction_param: dict,
    ) -> Tuple[OligoDatabase, str]:
        # define the filters
        hard_masked_sequences = HardMaskedSequenceFilter()
        soft_masked_sequences = SoftMaskedSequenceFilter()
        gc_content = GCContentFilter(GC_content_min=GC_content_min, GC_content_max=GC_content_max)
        melting_temperature = MeltingTemperatureNNFilter(
            Tm_min=Tm_min,
            Tm_max=Tm_max,
            Tm_parameters=Tm_parameters,
            Tm_chem_correction_parameters=Tm_chem_correction_param,
            Tm_salt_correction_parameters=Tm_salt_correction_param,
        )
        homopolymeric_runs = HomopolymericRunsFilter(
            base_n=homopolymeric_base_n,
        )
        # only need detection oligo filter because it also checks for Padlock arms
        detect_oligo_filter = DetectionOligoFilter(
            detect_oligo_length_min=detect_oligo_length_min,
            detect_oligo_length_max=detect_oligo_length_max,
            min_thymines=min_thymines,
            arm_length_min=arm_length_min,
            arm_Tm_dif_max=arm_Tm_dif_max,
            arm_Tm_min=arm_Tm_min,
            arm_Tm_max=arm_Tm_max,
            Tm_parameters=Tm_parameters,
            Tm_chem_correction_parameters=Tm_chem_correction_param,
            Tm_salt_correction_parameters=Tm_salt_correction_param,
        )

        filters = [
            hard_masked_sequences,
            soft_masked_sequences,
            homopolymeric_runs,
            gc_content,
            melting_temperature,
            detect_oligo_filter,
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
        specificity_blastn_search_parameters: dict,
        specificity_blastn_hit_parameters: dict,
        cross_hybridization_blastn_search_parameters: dict,
        cross_hybridization_blastn_hit_parameters: dict,
        ligation_region_size: int,
        arm_Tm_dif_max: int,
        arm_length_min: int,
        arm_Tm_min: float,
        arm_Tm_max: float,
        Tm_parameters: dict,
        Tm_chem_correction_param: dict,
        Tm_salt_correction_param: dict,
    ) -> Tuple[OligoDatabase, str]:
        ##### define reference database #####
        reference_database = ReferenceDatabase(
            database_name=self.subdir_db_reference, dir_output=self.dir_output
        )
        reference_database.load_database_from_fasta(
            files_fasta=files_fasta_reference_database, database_overwrite=False
        )

        ##### exact match filter #####
        # removing duplicated probes from the region with the most probes
        # exectute seperately before specificity filter to compute ligation side for less oligos
        exact_matches = ExactMatchFilter(policy=RemoveByLargerRegionPolicy(), filter_name="exact_match")
        filters = [exact_matches]
        specificity_filter = SpecificityFilter(filters=filters)
        oligo_database = specificity_filter.apply(
            sequence_type="oligo",
            oligo_database=oligo_database,
            reference_database=reference_database,
            n_jobs=self.n_jobs,
        )

        ##### calculate required probe attributes #####
        oligo_database = self.oligo_attributes_calculator.calculate_padlock_arms(
            oligo_database=oligo_database,
            arm_length_min=arm_length_min,
            arm_Tm_dif_max=arm_Tm_dif_max,
            arm_Tm_min=arm_Tm_min,
            arm_Tm_max=arm_Tm_max,
            Tm_parameters=Tm_parameters,
            Tm_chem_correction_parameters=Tm_chem_correction_param,
            Tm_salt_correction_parameters=Tm_salt_correction_param,
        )

        ##### specificity filters #####
        cross_hybridization_aligner = BlastNFilter(
            search_parameters=cross_hybridization_blastn_search_parameters,
            hit_parameters=cross_hybridization_blastn_hit_parameters,
            filter_name="blastn_crosshybridization",
            dir_output=self.dir_output,
        )
        cross_hybridization = CrossHybridizationFilter(
            policy=RemoveByLargerRegionPolicy(),
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

        filters = [specificity, cross_hybridization]
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

    @pipeline_step_basic(step_name="Set Selection")
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
        Tm_chem_correction_param: dict,
        Tm_salt_correction_param: dict,
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
            Tm_chem_correction_parameters=Tm_chem_correction_param,
            Tm_salt_correction_parameters=Tm_salt_correction_param,
            isoform_weight=isoform_weight,
            Tm_weight=Tm_weight,
            GC_weight=GC_weight,
        )
        set_scoring = LowestSetScoring(ascending=True)

        selection_policy = GraphBasedSelectionPolicy(
            set_scoring=set_scoring,
            pre_filter=pre_filter,
            n_attempts=n_attempts,
            heuristic=heuristic,
            heuristic_n_attempts=heuristic_n_attempts,
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


class DetectionOligoDesigner:
    def __init__(self, n_jobs: int) -> None:
        """Constructor for the ScrinshotProbeDesigner class."""

        ##### create the output folder #####
        self.n_jobs = n_jobs
        self.oligo_attributes_calculator = OligoAttributes()

    def create_detection_oligos(
        self,
        oligo_database: OligoDatabase,
        oligo_length_min: int,
        oligo_length_max: int,
        min_thymines: int,
        U_distance: int,
        Tm_opt: float,
        Tm_parameters: dict,
        Tm_chem_correction_parameters: dict,
        Tm_salt_correction_parameters: dict,
    ) -> dict:

        region_ids = list(oligo_database.database.keys())

        with joblib_progress(description="Design Detection Oligos", total=len(region_ids)):
            Parallel(
                n_jobs=self.n_jobs, prefer="threads", require="sharedmem"
            )(  # there should be an explicit return
                delayed(self._create_detection_oligos_region)(
                    oligo_database,
                    region_id,
                    oligo_length_min,
                    oligo_length_max,
                    min_thymines,
                    U_distance,
                    Tm_opt,
                    Tm_parameters,
                    Tm_chem_correction_parameters,
                    Tm_salt_correction_parameters,
                )
                for region_id in region_ids
            )

        return oligo_database

    def _create_detection_oligos_region(
        self,
        oligo_database: OligoDatabase,
        region_id: str,
        oligo_length_min: int,
        oligo_length_max: int,
        min_thymines: int,
        U_distance: int,
        Tm_opt: float,
        Tm_parameters: dict,
        Tm_chem_correction_parameters: dict,
        Tm_salt_correction_parameters: dict,
    ) -> dict:

        oligosets_region = oligo_database.oligosets[region_id]
        oligosets_oligo_columns = [col for col in oligosets_region.columns if col.startswith("oligo_")]

        new_oligo_attributes = {}

        for index in range(len(oligosets_region.index)):
            for column in oligosets_oligo_columns:
                oligo_id = str(oligosets_region.loc[index, column])

                ligation_site = oligo_database.get_oligo_attribute_value(
                    attribute="ligation_site", region_id=region_id, oligo_id=oligo_id, flatten=True
                )
                sequence_oligo = oligo_database.get_oligo_attribute_value(
                    attribute="oligo", region_id=region_id, oligo_id=oligo_id, flatten=True
                )

                (
                    detect_oligo_even,
                    detect_oligo_long_left,
                    detect_oligo_long_right,
                ) = self.oligo_attributes_calculator._calc_detect_oligo(
                    sequence=sequence_oligo,
                    ligation_site=ligation_site,
                    detect_oligo_length_min=oligo_length_min,
                    detect_oligo_length_max=oligo_length_max,
                    min_thymines=min_thymines,
                )

                # Search for best oligos
                initial_oligos = [
                    detect_oligo
                    for detect_oligo in [
                        detect_oligo_even,
                        detect_oligo_long_left,
                        detect_oligo_long_right,
                    ]
                    if (detect_oligo is not None) and (detect_oligo.count("T") >= min_thymines)
                ]

                # Check which of the three initial detection oligo is the best one
                Tm_dif = [
                    self._get_Tm_dif(
                        detect_oligo,
                        Tm_opt,
                        Tm_parameters,
                        Tm_chem_correction_parameters,
                        Tm_salt_correction_parameters,
                    )
                    for detect_oligo in initial_oligos
                ]
                best_initial_oligo = initial_oligos[Tm_dif.index(min(Tm_dif))]

                # Iterative search through shorter oligos
                oligos_cut_from_right, Tm_dif_cut_from_right = self._find_best_oligo(
                    best_initial_oligo,
                    cut_from_right=True,
                    oligo_length_min=oligo_length_min,
                    min_thymines=min_thymines,
                    Tm_opt=Tm_opt,
                    Tm_parameters=Tm_parameters,
                    Tm_chem_correction_parameters=Tm_chem_correction_parameters,
                    Tm_salt_correction_parameters=Tm_salt_correction_parameters,
                )
                oligos_cut_from_left, Tm_dif_cut_from_left = self._find_best_oligo(
                    best_initial_oligo,
                    cut_from_right=False,
                    oligo_length_min=oligo_length_min,
                    min_thymines=min_thymines,
                    Tm_opt=Tm_opt,
                    Tm_parameters=Tm_parameters,
                    Tm_chem_correction_parameters=Tm_chem_correction_parameters,
                    Tm_salt_correction_parameters=Tm_salt_correction_parameters,
                )
                oligos = oligos_cut_from_right + oligos_cut_from_left
                Tm_dif = Tm_dif_cut_from_right + Tm_dif_cut_from_left
                detection_oligo = oligos[Tm_dif.index(min(Tm_dif))]

                Tm_detection_oligo = self.oligo_attributes_calculator._calc_TmNN(
                    sequence=detection_oligo,
                    Tm_parameters=Tm_parameters,
                    Tm_chem_correction_parameters=Tm_chem_correction_parameters,
                    Tm_salt_correction_parameters=Tm_salt_correction_parameters,
                )

                # exchange T's with U (for enzymatic degradation of oligos)
                detection_oligo = self._exchange_T_with_U(detection_oligo, min_thymines, U_distance)

                new_oligo_attributes[oligo_id] = {
                    "Tm_detection_oligo": Tm_detection_oligo,
                    "sequence_detection_oligo": detection_oligo,
                }

        oligo_database.update_oligo_attributes(new_oligo_attributes)

    def _get_Tm_dif(
        self,
        oligo: str,
        Tm_opt: float,
        Tm_parameters: dict,
        Tm_chem_correction_parameters: dict,
        Tm_salt_correction_parameters: dict,
    ) -> int:
        Tm = self.oligo_attributes_calculator._calc_TmNN(
            sequence=oligo,
            Tm_parameters=Tm_parameters,
            Tm_chem_correction_parameters=Tm_chem_correction_parameters,
            Tm_salt_correction_parameters=Tm_salt_correction_parameters,
        )
        return abs(Tm - Tm_opt)

    def _find_best_oligo(
        self,
        oligo: str,
        cut_from_right: bool,
        oligo_length_min: int,
        min_thymines: int,
        Tm_opt: float,
        Tm_parameters: dict,
        Tm_chem_correction_parameters: dict,
        Tm_salt_correction_parameters: dict,
    ) -> Tuple[list, list]:

        oligos = [oligo]
        Tm_dif = [
            self._get_Tm_dif(
                oligo, Tm_opt, Tm_parameters, Tm_chem_correction_parameters, Tm_salt_correction_parameters
            )
        ]

        # either start cut from left or right and make sure that oligo length is >= oligo_length_min
        for count in range(0, len(oligo) - oligo_length_min):
            if bool(count % 2) * cut_from_right:
                oligo = oligo[1:]
            else:
                oligo = oligo[:-1]

            if oligo.count("T") >= min_thymines:
                oligos.append(oligo)
                Tm_dif.append(
                    self._get_Tm_dif(
                        oligo,
                        Tm_opt,
                        Tm_parameters,
                        Tm_chem_correction_parameters,
                        Tm_salt_correction_parameters,
                    )
                )

        return oligos, Tm_dif

    def _exchange_T_with_U(self, oligo: str, min_thymines: int, U_distance: int) -> str:
        if oligo.find("T") < oligo[::-1].find("T"):
            fluorophor_pos = "left"
        else:
            fluorophor_pos = "right"
            oligo = oligo[::-1]

        pos = 0
        new_pos = 1
        for _ in range(min_thymines):
            while True:
                shift = 0 if (pos == 0 and (new_pos != 0)) else U_distance
                start = min(pos + shift, len(oligo))
                new_pos = oligo[start:].find("T")
                if new_pos == -1:
                    pos = oligo.rfind("T") - U_distance
                else:
                    pos = pos + shift + new_pos
                    oligo = oligo[:pos] + "U" + oligo[pos + 1 :]
                    break

        # Add fluorophore
        if fluorophor_pos == "left":
            oligo = "[fluorophore]" + oligo
        elif fluorophor_pos == "right":
            oligo = oligo[::-1] + "[fluorophore]"

        return oligo


############################################
# SCRINSHOT Probe Designer Pipeline
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
    pipeline = ScrinshotProbeDesigner(
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
        max_graph_size=config["max_graph_size"],
        pre_filter=config["pre_filtering"],
        n_attempts=config["n_attempts"],
        heuristic=config["heuristic"],
        heuristic_n_attempts=config["heuristic_n_attempts"],
        target_probe_Tm_parameters=config["target_probe_Tm_parameters"],
        target_probe_Tm_chem_correction_param=config["target_probe_Tm_chem_correction_param"],
        target_probe_Tm_salt_correction_param=config["target_probe_Tm_salt_correction_param"],
        detection_oligo_Tm_parameters=config["detection_oligo_Tm_parameters"],
        detection_oligo_Tm_chem_correction_param=config["detection_oligo_Tm_chem_correction_param"],
        detection_oligo_Tm_salt_correction_param=config["detection_oligo_Tm_salt_correction_param"],
    )

    ##### design probes #####
    oligo_database = pipeline.design_target_probes(
        gene_ids=gene_ids,
        files_fasta_target_probe_database=config["files_fasta_target_probe_database"],
        files_fasta_reference_database_targe_probe=config["files_fasta_reference_database_targe_probe"],
        target_probe_length_min=config["target_probe_length_min"],
        target_probe_length_max=config["target_probe_length_max"],
        target_probe_isoform_consensus=config["target_probe_isoform_consensus"],
        target_probe_isoform_weight=config["target_probe_isoform_weight"],
        target_probe_GC_content_min=config["target_probe_GC_content_min"],
        target_probe_GC_content_opt=config["target_probe_GC_content_opt"],
        target_probe_GC_content_max=config["target_probe_GC_content_max"],
        target_probe_GC_weight=config["target_probe_GC_weight"],
        target_probe_Tm_min=config["target_probe_Tm_min"],
        target_probe_Tm_opt=config["target_probe_Tm_opt"],
        target_probe_Tm_max=config["target_probe_Tm_max"],
        target_probe_Tm_weight=config["target_probe_Tm_weight"],
        target_probe_homopolymeric_base_n=config["target_probe_homopolymeric_base_n"],
        detection_oligo_min_thymines=config["detection_oligo_min_thymines"],
        detection_oligo_length_min=config["detection_oligo_length_min"],
        detection_oligo_length_max=config["detection_oligo_length_max"],
        target_probe_padlock_arm_length_min=config["target_probe_padlock_arm_length_min"],
        target_probe_padlock_arm_Tm_dif_max=config["target_probe_padlock_arm_Tm_dif_max"],
        target_probe_padlock_arm_Tm_min=config["target_probe_padlock_arm_Tm_min"],
        target_probe_padlock_arm_Tm_max=config["target_probe_padlock_arm_Tm_max"],
        target_probe_ligation_region_size=config["target_probe_ligation_region_size"],
        set_size_min=config["set_size_min"],
        set_size_opt=config["set_size_opt"],
        distance_between_target_probes=config["distance_between_target_probes"],
        n_sets=config["n_sets"],
    )

    oligo_database = pipeline.design_detection_oligos(
        oligo_database=oligo_database,
        detection_oligo_length_min=config["detection_oligo_length_min"],
        detection_oligo_length_max=config["detection_oligo_length_max"],
        detection_oligo_min_thymines=config["detection_oligo_min_thymines"],
        detection_oligo_U_distance=config["detection_oligo_U_distance"],
        detection_oligo_Tm_opt=config["detection_oligo_Tm_opt"],
    )

    oligo_database = pipeline.design_padlock_backbone(oligo_database=oligo_database)

    pipeline.generate_output(oligo_database=oligo_database, top_n_sets=config["top_n_sets"])

    print("--------------END PIPELINE--------------")


if __name__ == "__main__":
    main()
