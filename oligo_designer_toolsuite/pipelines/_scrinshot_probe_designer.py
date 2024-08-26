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
from typing import List

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
from oligo_designer_toolsuite.pipelines._utils import base_parser, pipeline_step_basic
from oligo_designer_toolsuite.sequence_generator import OligoSequenceGenerator

############################################
# SCRINSHOT Probe Designer Functions
############################################


class ScrinshotProbeDesigner:
    def __init__(self, write_intermediate_steps: bool, dir_output: str, n_jobs: int) -> None:
        """Constructor for the ScrinshotProbeDesigner class."""

        self.write_intermediate_steps = write_intermediate_steps
        self.n_jobs = n_jobs
        self.probe_attributes_calculator = OligoAttributes()

        ##### create the output folder #####
        self.dir_output = os.path.abspath(dir_output)
        Path(self.dir_output).mkdir(parents=True, exist_ok=True)

        self.subdir_db_probes = "db_probes"
        self.subdir_db_reference = "db_reference"

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

    @pipeline_step_basic(step_name="Create Database")
    def create_probe_database(
        self,
        gene_ids: list,
        probe_length_min: int,
        probe_length_max: int,
        files_fasta_oligo_database: list[str],
        min_probes_per_gene: int,
    ):
        ##### creating the probe sequences #####
        probe_sequences = OligoSequenceGenerator(dir_output=self.dir_output)
        probe_fasta_file = probe_sequences.create_sequences_sliding_window(
            files_fasta_in=files_fasta_oligo_database,
            length_interval_sequences=(probe_length_min, probe_length_max),
            region_ids=gene_ids,
            n_jobs=self.n_jobs,
        )

        ##### creating the probe database #####
        oligo_database = OligoDatabase(
            min_oligos_per_region=min_probes_per_gene,
            write_regions_with_insufficient_oligos=True,
            lru_db_max_in_memory=self.n_jobs * 2 + 2,
            database_name=self.subdir_db_probes,
            dir_output=self.dir_output,
            n_jobs=1,
        )
        oligo_database.load_database_from_fasta(
            files_fasta=probe_fasta_file,
            sequence_type="target",
            region_ids=gene_ids,
        )

        ##### save database #####
        if self.write_intermediate_steps:
            file_database = oligo_database.save_database(dir_database="1_db_probes_initial")
        else:
            file_database = ""

        dir = probe_sequences.dir_output
        shutil.rmtree(dir) if os.path.exists(dir) else None

        return oligo_database, file_database

    @pipeline_step_basic(step_name="Property Filters")
    def filter_by_property(
        self,
        oligo_database: OligoDatabase,
        probe_GC_content_min: float,
        probe_GC_content_max: float,
        probe_Tm_min: float,
        probe_Tm_max: float,
        detect_oligo_length_min: int,
        detect_oligo_length_max: int,
        min_thymines: int,
        arm_Tm_dif_max: int,
        arm_length_min: int,
        arm_Tm_min: float,
        arm_Tm_max: float,
        homopolymeric_base_n: str,
        Tm_parameters_probe: dict,
        Tm_chem_correction_param_probe: dict,
        Tm_salt_correction_param_probe: dict,
    ):
        # define the filters
        hard_masked_sequences = HardMaskedSequenceFilter()
        soft_masked_sequences = SoftMaskedSequenceFilter()
        gc_content = GCContentFilter(GC_content_min=probe_GC_content_min, GC_content_max=probe_GC_content_max)
        melting_temperature = MeltingTemperatureNNFilter(
            Tm_min=probe_Tm_min,
            Tm_max=probe_Tm_max,
            Tm_parameters=Tm_parameters_probe,
            Tm_chem_correction_parameters=Tm_chem_correction_param_probe,
            Tm_salt_correction_parameters=Tm_salt_correction_param_probe,
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
            Tm_parameters=Tm_parameters_probe,
            Tm_chem_correction_parameters=Tm_chem_correction_param_probe,
            Tm_salt_correction_parameters=Tm_salt_correction_param_probe,
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
            sequence_type="oligo",
            oligo_database=oligo_database,
            n_jobs=self.n_jobs,
        )

        # write the intermediate result in a file
        if self.write_intermediate_steps:
            file_database = oligo_database.save_database(dir_database="2_db_probes_property_filter")
        else:
            file_database = ""

        return oligo_database, file_database

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
        Tm_parameters_probe: dict,
        Tm_chem_correction_param_probe: dict,
        Tm_salt_correction_param_probe: dict,
    ):
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
        oligo_database = self.probe_attributes_calculator.calculate_padlock_arms(
            oligo_database=oligo_database,
            arm_length_min=arm_length_min,
            arm_Tm_dif_max=arm_Tm_dif_max,
            arm_Tm_min=arm_Tm_min,
            arm_Tm_max=arm_Tm_max,
            Tm_parameters=Tm_parameters_probe,
            Tm_chem_correction_parameters=Tm_chem_correction_param_probe,
            Tm_salt_correction_parameters=Tm_salt_correction_param_probe,
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

        # write the intermediate result in a file
        if self.write_intermediate_steps:
            file_database = oligo_database.save_database(dir_database="3_db_probes_specificity_filter")
        else:
            file_database = ""

        # remove all directories of intermediate steps
        for directory in [
            reference_database.dir_output,
            cross_hybridization_aligner.dir_output,
            cross_hybridization.dir_output,
            specificity.dir_output,
        ]:
            if os.path.exists(directory):
                shutil.rmtree(directory)

        return oligo_database, file_database

    @pipeline_step_basic(step_name="Set Selection")
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
        Tm_salt_correction_param_probe: dict,
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
            Tm_salt_correction_parameters=Tm_salt_correction_param_probe,
            isoform_weight=probe_isoform_weight,
            Tm_weight=probe_Tm_weight,
            GC_weight=probe_GC_weight,
        )
        set_scoring = LowestSetScoring(ascending=True)

        selection_policy = GraphBasedSelectionPolicy(
            set_size_opt=probeset_size_opt,
            set_size_min=probeset_size_min,
            n_sets=n_sets,
            ascending=True,
            set_scoring=set_scoring,
        )
        probeset_generator = OligosetGeneratorIndependentSet(
            opt_oligoset_size=probeset_size_opt,
            min_oligoset_size=probeset_size_min,
            oligos_scoring=probes_scoring,
            set_scoring=set_scoring,
            heuristic_selection=selection_policy,
            max_oligos=max_graph_size,
            distance_between_oligos=distance_between_probes,
        )
        oligo_database = probeset_generator.apply(
            oligo_database=oligo_database,
            sequence_type="oligo",
            pre_filter=False,
            n_sets=n_sets,
            n_jobs=self.n_jobs,
        )

        # write the intermediate result in a file
        if self.write_intermediate_steps:
            file_database = oligo_database.save_database(dir_database="4_db_probes_probesets")
            file_probesets = oligo_database.write_oligosets_to_table()
        else:
            file_database = ""
            file_probesets = ""

        return oligo_database, file_database, file_probesets

    def design_final_probe_sequence(
        self,
        oligo_database: OligoDatabase,
        min_thymines: int,
        U_distance: int,
        detect_oligo_length_min: int,
        detect_oligo_length_max: int,
        detect_oligo_Tm_opt: float,
        Tm_parameters_detection_oligo: dict,
        Tm_chem_correction_param_detection_oligo: dict,
        Tm_salt_correction_param_detection_oligo: dict,
        Tm_parameters_probe: dict,
        Tm_chem_correction_param_probe: dict,
        Tm_salt_correction_param_probe: dict,
    ):
        """ """

        def _get_barcode(number_regions: int, barcode_length: int, seed: int, choices: list):

            while len(choices) ** barcode_length < number_regions:
                barcode_length += 1

            barcodes = ["".join(nts) for nts in itertools.product(choices, repeat=barcode_length)]
            random.seed(seed)
            random.shuffle(barcodes)

            return barcodes

        def _get_padlock_probe(probe_attributes: dict):

            ligation_site = probe_attributes["ligation_site"]
            probe_attributes["sequence_padlock_arm1"] = probe_attributes["oligo"][ligation_site:]
            probe_attributes["sequence_padlock_arm2"] = probe_attributes["oligo"][:ligation_site]

            probe_attributes["sequence_padlock_accessory1"] = "TCCTCTATGATTACTGAC"
            probe_attributes["sequence_padlock_ISS_anchor"] = "TGCGTCTATTTAGTGGAGCC"
            probe_attributes["sequence_padlock_accessory2"] = "CTATCTTCTTT"

            probe_attributes["sequence_padlock_backbone"] = (
                probe_attributes["sequence_padlock_accessory1"]
                + probe_attributes["sequence_padlock_ISS_anchor"]
                + probe_attributes["barcode"]
                + probe_attributes["sequence_padlock_accessory2"]
            )

            probe_attributes["sequence_padlock_probe"] = (
                probe_attributes["sequence_padlock_arm1"]
                + probe_attributes["sequence_padlock_backbone"]
                + probe_attributes["sequence_padlock_arm2"]
            )

            probe_attributes["Tm_arm1"] = self.probe_attributes_calculator._calc_TmNN(
                sequence=probe_attributes["sequence_padlock_arm1"],
                Tm_parameters=Tm_parameters_probe,
                Tm_chem_correction_parameters=Tm_chem_correction_param_probe,
                Tm_salt_correction_parameters=Tm_salt_correction_param_probe,
            )
            probe_attributes["Tm_arm2"] = self.probe_attributes_calculator._calc_TmNN(
                sequence=probe_attributes["sequence_padlock_arm2"],
                Tm_parameters=Tm_parameters_probe,
                Tm_chem_correction_parameters=Tm_chem_correction_param_probe,
                Tm_salt_correction_parameters=Tm_salt_correction_param_probe,
            )
            probe_attributes["Tm_diff_arms"] = round(
                abs(probe_attributes["Tm_arm1"] - probe_attributes["Tm_arm2"]), 2
            )

            return probe_attributes

        def _get_detection_oligo(probe_attributes: dict):
            def get_Tm_dif(oligo):
                Tm = self.probe_attributes_calculator._calc_TmNN(
                    sequence=oligo,
                    Tm_parameters=Tm_parameters_detection_oligo,
                    Tm_chem_correction_parameters=Tm_chem_correction_param_detection_oligo,
                    Tm_salt_correction_parameters=Tm_salt_correction_param_detection_oligo,
                )
                return abs(Tm - detect_oligo_Tm_opt)

            def _find_best_oligo(oligo, cut_from_right):

                oligos = [oligo]
                Tm_dif = [get_Tm_dif(oligo)]

                # either start cut from left or right and make sure that oligo length is >= detect_oligo_length_min
                for count in range(0, len(oligo) - detect_oligo_length_min):
                    if bool(count % 2) * cut_from_right:
                        oligo = oligo[1:]
                    else:
                        oligo = oligo[:-1]

                    if oligo.count("T") >= min_thymines:
                        oligos.append(oligo)
                        Tm_dif.append(get_Tm_dif(oligo))

                return oligos, Tm_dif

            def _exchange_T_with_U(oligo):
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

            (
                detect_oligo_even,
                detect_oligo_long_left,
                detect_oligo_long_right,
            ) = self.probe_attributes_calculator._calc_detect_oligo(
                sequence=probe_attributes["oligo"],
                ligation_site=probe_attributes["ligation_site"],
                detect_oligo_length_min=detect_oligo_length_min,
                detect_oligo_length_max=detect_oligo_length_max,
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
            Tm_dif = [get_Tm_dif(detect_oligo) for detect_oligo in initial_oligos]
            best_initial_oligo = initial_oligos[Tm_dif.index(min(Tm_dif))]

            # Iterative search through shorter oligos
            oligos_cut_from_right, Tm_dif_cut_from_right = _find_best_oligo(
                best_initial_oligo, cut_from_right=True
            )
            oligos_cut_from_left, Tm_dif_cut_from_left = _find_best_oligo(
                best_initial_oligo, cut_from_right=False
            )
            oligos = oligos_cut_from_right + oligos_cut_from_left
            Tm_dif = Tm_dif_cut_from_right + Tm_dif_cut_from_left
            detection_oligo = oligos[Tm_dif.index(min(Tm_dif))]

            probe_attributes["Tm_detection_oligo"] = self.probe_attributes_calculator._calc_TmNN(
                sequence=detection_oligo,
                Tm_parameters=Tm_parameters_detection_oligo,
                Tm_chem_correction_parameters=Tm_chem_correction_param_detection_oligo,
                Tm_salt_correction_parameters=Tm_salt_correction_param_detection_oligo,
            )

            # exchange T's with U (for enzymatic degradation of oligos)
            detection_oligo = _exchange_T_with_U(detection_oligo)
            probe_attributes["sequence_detection_oligo"] = detection_oligo

            return probe_attributes

        def _assemble_sequence(oligo_database, region_id, region_idx, barcodes):
            database_region = oligo_database.database[region_id]
            probesets_region = oligo_database.oligosets[region_id]
            probesets_probe_columns = [col for col in probesets_region.columns if col.startswith("oligo_")]

            for index in range(len(probesets_region.index)):
                for column in probesets_probe_columns:
                    probe_id = str(probesets_region.loc[index, column])
                    probe_attributes = database_region[probe_id]

                    probe_attributes["barcode"] = barcodes[region_idx]
                    probe_attributes["sequence_mRNA"] = probe_attributes["target"]
                    probe_attributes["sequence_mRNA_probe"] = probe_attributes["oligo"]

                    probe_attributes = _get_padlock_probe(probe_attributes)
                    probe_attributes = _get_detection_oligo(probe_attributes)

                    oligo_database.database[region_id][probe_id] = probe_attributes

        region_ids = list(oligo_database.database.keys())

        barcodes = _get_barcode(len(region_ids), barcode_length=4, seed=0, choices=["A", "C", "T", "G"])

        with joblib_progress(description="Design Final Padlock Sequence", total=len(region_ids)):
            Parallel(
                n_jobs=self.n_jobs, prefer="threads", require="sharedmem"
            )(  # there should be an explicit return
                delayed(_assemble_sequence)(oligo_database, region_id, region_idx, barcodes)
                for region_idx, region_id in enumerate(region_ids)
            )

        return oligo_database

    def compute_probe_attributes(
        self,
        oligo_database: OligoDatabase,
        Tm_parameters_probe: dict,
        Tm_chem_correction_param_probe: dict,
        Tm_salt_correction_param_probe: dict,
    ):
        oligo_database = self.probe_attributes_calculator.calculate_oligo_length(
            oligo_database=oligo_database
        )
        oligo_database = self.probe_attributes_calculator.calculate_GC_content(
            oligo_database=oligo_database, sequence_type="oligo"
        )
        oligo_database = self.probe_attributes_calculator.calculate_TmNN(
            oligo_database=oligo_database,
            sequence_type="oligo",
            Tm_parameters=Tm_parameters_probe,
            Tm_chem_correction_parameters=Tm_chem_correction_param_probe,
            Tm_salt_correction_parameters=Tm_salt_correction_param_probe,
        )
        oligo_database = self.probe_attributes_calculator.calculate_num_targeted_transcripts(
            oligo_database=oligo_database
        )
        oligo_database = self.probe_attributes_calculator.calculate_isoform_consensus(
            oligo_database=oligo_database
        )

        return oligo_database

    def generate_output(self, oligo_database: OligoDatabase, top_n_sets: int):

        attributes = [
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
            "sequence_mRNA",
            "sequence_mRNA_probe",
            "length",
            "ligation_site",
            "Tm_arm1",
            "Tm_arm2",
            "Tm_diff_arms",
            "Tm_detection_oligo",
            "source",
            "species",
            "annotation_release",
            "genome_assembly",
            "regiontype",
            "gene_id",
            "transcript_id",
            "exon_number",
        ]
        oligo_database.write_oligosets_to_yaml(
            attributes=attributes,
            top_n_sets=top_n_sets,
            ascending=True,
            filename="padlock_probes.yml",
        )

        # write a second file that only contains order information
        yaml_dict_order = {}

        for region_id, database_region in oligo_database.database.items():
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
                        "sequence_padlock_probe": database_region[oligo_id]["sequence_padlock_probe"],
                        "sequence_detection_oligo": database_region[oligo_id]["sequence_detection_oligo"],
                    }

        with open(os.path.join(self.dir_output, "padlock_probes_order.yml"), "w") as outfile:
            yaml.dump(yaml_dict_order, outfile, default_flow_style=False, sort_keys=False)


############################################
# SCRINSHOT Probe Designer Pipeline
############################################


def main():

    args = base_parser()

    ##### read the config file #####
    with open(args["config"], "r") as handle:
        config = yaml.safe_load(handle)

    ##### process parameters #####
    # preprocess melting temperature params
    Tm_parameters_probe = config["Tm_parameters_probe"]
    Tm_parameters_probe["nn_table"] = getattr(mt, Tm_parameters_probe["nn_table"])
    Tm_parameters_probe["tmm_table"] = getattr(mt, Tm_parameters_probe["tmm_table"])
    Tm_parameters_probe["imm_table"] = getattr(mt, Tm_parameters_probe["imm_table"])
    Tm_parameters_probe["de_table"] = getattr(mt, Tm_parameters_probe["de_table"])

    # preprocess melting temperature params
    Tm_parameters_detection_oligo = config["Tm_parameters_detection_oligo"]
    Tm_parameters_detection_oligo["nn_table"] = getattr(mt, Tm_parameters_detection_oligo["nn_table"])
    Tm_parameters_detection_oligo["tmm_table"] = getattr(mt, Tm_parameters_detection_oligo["tmm_table"])
    Tm_parameters_detection_oligo["imm_table"] = getattr(mt, Tm_parameters_detection_oligo["imm_table"])
    Tm_parameters_detection_oligo["de_table"] = getattr(mt, Tm_parameters_detection_oligo["de_table"])

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

    ##### create probe database #####
    probe_database, file_database = pipeline.create_probe_database(
        gene_ids=gene_ids,
        probe_length_min=config["probe_length_min"],
        probe_length_max=config["probe_length_max"],
        files_fasta_oligo_database=config["files_fasta_probe_database"],
        # we should have at least "min_probeset_size" probes per gene to create one set
        min_probes_per_gene=config["probeset_size_min"],
    )

    ##### filter probes by property #####
    probe_database, file_database = pipeline.filter_by_property(
        oligo_database=probe_database,
        probe_GC_content_min=config["probe_GC_content_min"],
        probe_GC_content_max=config["probe_GC_content_max"],
        probe_Tm_min=config["probe_Tm_min"],
        probe_Tm_max=config["probe_Tm_max"],
        detect_oligo_length_min=config["detect_oligo_length_min"],
        detect_oligo_length_max=config["detect_oligo_length_max"],
        min_thymines=config["min_thymines"],
        arm_Tm_dif_max=config["arm_Tm_dif_max"],
        arm_length_min=config["arm_length_min"],
        arm_Tm_min=config["arm_Tm_min"],
        arm_Tm_max=config["arm_Tm_max"],
        homopolymeric_base_n=config["homopolymeric_base_n"],
        Tm_parameters_probe=Tm_parameters_probe,
        Tm_chem_correction_param_probe=config["Tm_chem_correction_param_probe"],
        Tm_salt_correction_param_probe=config["Tm_salt_correction_param_probe"],
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
        Tm_salt_correction_param_probe=config["Tm_salt_correction_param_probe"],
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
        Tm_salt_correction_param_probe=config["Tm_salt_correction_param_probe"],
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

    ##### create final padlock probe sequences #####
    probe_database = pipeline.design_final_probe_sequence(
        oligo_database=probe_database,
        min_thymines=config["min_thymines"],
        U_distance=config["U_distance"],
        detect_oligo_length_min=config["detect_oligo_length_min"],
        detect_oligo_length_max=config["detect_oligo_length_max"],
        detect_oligo_Tm_opt=config["detect_oligo_Tm_opt"],
        Tm_parameters_detection_oligo=Tm_parameters_detection_oligo,
        Tm_chem_correction_param_detection_oligo=config["Tm_chem_correction_param_detection_oligo"],
        Tm_salt_correction_param_detection_oligo=config["Tm_salt_correction_param_detection_oligo"],
        Tm_parameters_probe=Tm_parameters_probe,
        Tm_chem_correction_param_probe=config["Tm_chem_correction_param_probe"],
        Tm_salt_correction_param_probe=config["Tm_salt_correction_param_probe"],
    )

    ##### generate output #####
    # compute all required attributes
    probe_database = pipeline.compute_probe_attributes(
        oligo_database=probe_database,
        Tm_parameters_probe=Tm_parameters_probe,
        Tm_chem_correction_param_probe=config["Tm_chem_correction_param_probe"],
        Tm_salt_correction_param_probe=config["Tm_salt_correction_param_probe"],
    )

    pipeline.generate_output(oligo_database=probe_database, top_n_sets=config["top_n_sets"])

    logging.info(f"Oligo sets were saved in {dir_probesets}")
    logging.info("##### End of the pipeline. #####")


if __name__ == "__main__":
    main()
