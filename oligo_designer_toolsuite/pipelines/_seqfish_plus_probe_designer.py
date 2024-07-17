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

from typing import List
from pathlib import Path
from datetime import datetime
from itertools import product
from joblib import Parallel, delayed
from joblib_progress import joblib_progress

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
from oligo_designer_toolsuite.pipelines._utils import base_parser, pipeline_step_basic
from oligo_designer_toolsuite.sequence_generator import OligoSequenceGenerator

############################################
# SeqFISH Plus Probe Designer
############################################


class SeqFishPlusProbeDesigner:
    def __init__(self, write_intermediate_steps: bool, dir_output: str, n_jobs: int) -> None:
        """Constructor for the SeqFishPlusProbeDesigner class."""

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
            f"log_seqfishplus_probe_designer_{timestamp.year}-{timestamp.month}-{timestamp.day}-{timestamp.hour}-{timestamp.minute}.txt",
        )
        logging.getLogger("log_name")
        logging.basicConfig(
            format="%(asctime)s [%(levelname)s] %(message)s",
            level=logging.NOTSET,
            handlers=[logging.FileHandler(file_logger)],
        )
        logging.captureWarnings(True)

    @pipeline_step_basic(step_name="Target Probe Generation - Create Database")
    def create_oligo_database(
        self,
        gene_ids: list,
        oligo_length_min: int,
        oligo_length_max: int,
        files_fasta_oligo_database: list[str],
        min_probes_per_gene: int,
        isoform_consensus: float,
    ):
        ##### creating the oligo sequences #####
        probe_sequences = OligoSequenceGenerator(dir_output=self.dir_output)
        probe_fasta_file = probe_sequences.create_sequences_sliding_window(
            files_fasta_in=files_fasta_oligo_database,
            length_interval_sequences=(oligo_length_min, oligo_length_max),
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

        oligo_database = self.probe_attributes_calculator.calculate_isoform_consensus(
            oligo_database=oligo_database
        )
        oligo_database.filter_oligo_attribute(
            name_attribute="isoform_consensus",
            thr_attribute=isoform_consensus,
            keep_if_smaller_threshold=False,
        )

        ##### save database #####
        if self.write_intermediate_steps:
            file_database = oligo_database.save_database(dir_database="1_db_probes_initial")
        else:
            file_database = ""

        dir = probe_sequences.dir_output
        shutil.rmtree(dir) if os.path.exists(dir) else None

        return oligo_database, file_database

    @pipeline_step_basic(step_name="Target Probe Generation - Property Filters")
    def filter_by_property(
        self,
        oligo_database: OligoDatabase,
        GC_content_min: float,
        GC_content_max: float,
        homopolymeric_base_n: str,
        T_secondary_structure: float,
        secondary_structures_threshold_deltaG: float,
    ):
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

    @pipeline_step_basic(step_name="Target Probe Generation - Specificity Filters")
    def filter_by_specificity(
        self,
        oligo_database: OligoDatabase,
        files_fasta_reference_database: List[str],
        specificity_blastn_search_params: dict,
        specificity_blastn_hit_params: dict,
        cross_hybridization_blastn_search_params: dict,
        cross_hybridization_blastn_hit_params: dict,
    ):
        ##### define reference database #####
        reference_database = ReferenceDatabase(
            database_name=self.subdir_db_reference, dir_output=self.dir_output
        )
        reference_database.load_database_from_fasta(
            files_fasta=files_fasta_reference_database, database_overwrite=False
        )

        ##### exact match filter #####
        exact_matches = ExactMatchFilter(policy=RemoveByLargerRegionPolicy(), filter_name="probe_exact_match")

        ##### specificity filters #####
        specificity = BlastNFilter(
            search_parameters=specificity_blastn_search_params,
            hit_parameters=specificity_blastn_hit_params,
            filter_name="probe_blastn_specificity",
            dir_output=self.dir_output,
        )

        cross_hybridization_aligner = BlastNFilter(
            search_parameters=cross_hybridization_blastn_search_params,
            hit_parameters=cross_hybridization_blastn_hit_params,
            filter_name="probe_blastn_crosshybridization",
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
        distance_between_oligos: int,
    ):
        probes_scoring = WeightedGCUtrScoring(
            GC_content_opt=GC_content_opt, GC_weight=GC_weight, UTR_weight=UTR_weight
        )
        set_scoring = LowestSetScoring(ascending=True)
        probeset_generator = OligosetGeneratorIndependentSet(
            opt_oligoset_size=set_size_opt,
            min_oligoset_size=set_size_min,
            oligos_scoring=probes_scoring,
            set_scoring=set_scoring,
            heuristic_selection=heuristic_selection_independent_set,
            max_oligos=max_graph_size,
            distance_between_oligos=distance_between_oligos,
        )
        oligo_database = probeset_generator.apply(
            oligo_database=oligo_database,
            sequence_type="oligo",
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
        readout_probe_length: int,
        readout_probe_base_probabilities: dict,
        readout_probe_initial_num_sequences: int,
        readout_probe_GC_content_min: float,
        readout_probe_GC_content_max: float,
        readout_probe_homopolymeric_base_n: dict,
        readout_probe_files_fasta_reference_database: List[str],
        readout_probe_specificity_blastn_search_params: dict,
        readout_probe_specificity_blastn_hit_params: dict,
        readout_probe_cross_hybridization_blastn_search_params: dict,
        readout_probe_cross_hybridization_blastn_hit_params: dict,
        n_barcode_rounds: int,
        n_pseudocolors: int,
        channels_ids: list,
        forward_primer_length: int,
        forward_primer_base_probabilities: dict,
        forward_primer_initial_num_sequences: int,
        reverse_primer_sequence: str,
        forward_primer_GC_content_min: float,
        forward_primer_GC_content_max: float,
        forward_primer_number_GC_GCclamp: int,
        forward_primer_number_three_prime_base_GCclamp: int,
        forward_primer_homopolymeric_base_n: int,
        forward_primer_max_len_selfcomplement: int,
        forward_primer_max_len_complement: int,
        forward_primer_Tm_min: float,
        forward_primer_Tm_max: float,
        forward_primer_Tm_parameters: dict,
        forward_primer_Tm_chem_correction_param: dict,
        forward_primer_Tm_salt_correction_parameters: dict,
        forward_primer_T_secondary_structure: float,
        forward_primer_secondary_structures_threshold_deltaG: float,
        forward_primer_files_fasta_reference_database: List[str],
        forward_primer_specificity_refrence_blastn_search_params: dict,
        forward_primer_specificity_refrence_blastn_hit_params: dict,
        forward_primer_specificity_encoding_probes_blastn_search_parameters: dict,
        forward_primer_specificity_encoding_probes_blastn_hit_parameters: dict,
    ):
        def _generate_readout_probes(
            n_genes: int,
        ):
            pipeline = SeqFishPlusReadoutProbeDesigner(
                write_intermediate_steps=self.write_intermediate_steps,
                dir_output=self.dir_output,
                n_jobs=self.n_jobs,
            )
            readout_probe_database, _ = pipeline.create_oligo_database(
                oligo_length=readout_probe_length,
                oligo_base_probabilities=readout_probe_base_probabilities,
                initial_num_sequences=readout_probe_initial_num_sequences,
            )
            readout_probe_database, _ = pipeline.filter_by_property(
                oligo_database=readout_probe_database,
                GC_content_min=readout_probe_GC_content_min,
                GC_content_max=readout_probe_GC_content_max,
                homopolymeric_base_n=readout_probe_homopolymeric_base_n,
            )
            readout_probe_database, _ = pipeline.filter_by_specificity(
                oligo_database=readout_probe_database,
                files_fasta_reference_database=readout_probe_files_fasta_reference_database,
                specificity_blastn_search_parameters=readout_probe_specificity_blastn_search_params,
                specificity_blastn_hit_parameters=readout_probe_specificity_blastn_hit_params,
                cross_hybridization_blastn_search_parameters=readout_probe_cross_hybridization_blastn_search_params,
                cross_hybridization_blastn_hit_parameters=readout_probe_cross_hybridization_blastn_hit_params,
            )
            codebook = pipeline.generate_codebook(
                n_regions=n_genes,
                n_barcode_rounds=n_barcode_rounds,
                n_pseudocolors=n_pseudocolors,
                n_channels=len(channels_ids),
            )
            readout_probe_table = pipeline.create_readout_probe_table(
                readout_probe_database=readout_probe_database,
                channels_ids=channels_ids,
                n_barcode_rounds=n_barcode_rounds,
                n_pseudocolors=n_pseudocolors,
            )
            return codebook, readout_probe_table

        def _assemble_encoding_sequence(oligo_database, region_id, barcode, readout_probe_table):

            bits = barcode[barcode == 1].index
            readout_probe_sequences = readout_probe_table.loc[bits, "readout_probe_sequence"]

            database_region = oligo_database.database[region_id]
            for probe_id, probe_attributes in database_region.items():
                probe_attributes["barcode"] = barcode
                probe_attributes["sequence_target"] = probe_attributes["target"]
                probe_attributes["sequence_target_probe"] = probe_attributes["oligo"]

                for i, readout_probe_sequence in enumerate(readout_probe_sequences):
                    probe_attributes[f"sequence_readout_probe_{i+1}"] = readout_probe_sequence

                probe_attributes["sequence_encoding_probe"] = (
                    str(Seq(probe_attributes["sequence_readout_probe_1"]).reverse_complement())
                    + str(Seq(probe_attributes["sequence_readout_probe_2"]).reverse_complement())
                    + "T"
                    + probe_attributes["oligo"]
                    + "T"
                    + str(Seq(probe_attributes["sequence_readout_probe_3"]).reverse_complement())
                    + str(Seq(probe_attributes["sequence_readout_probe_4"]).reverse_complement())
                )

                oligo_database.database[region_id][probe_id] = probe_attributes

        def _generate_primers(file_fasta_encoding_probes_database):

            def _get_forward_primer_opt():
                # get the sequence with Tm closest to the Tm of the reverse primer
                forward_primer_sequence = ""
                min_dif_Tm = 100

                # calculate Tm for the reverse primer
                Tm_reverse_primer = OligoAttributes._calc_TmNN(
                    sequence=reverse_primer_sequence,
                    Tm_parameters=forward_primer_Tm_parameters,
                    Tm_chem_correction_parameters=forward_primer_Tm_chem_correction_param,
                    Tm_salt_correction_parameters=forward_primer_Tm_salt_correction_parameters,
                )

                # iterate over all primers in the database to find the one with Tm closest to the reverse primer Tm
                for database_region in primer_database.database.values():
                    for primer_attributes in database_region.values():
                        Tm_forward_primer = OligoAttributes._calc_TmNN(
                            sequence=primer_attributes["oligo"],
                            Tm_parameters=forward_primer_Tm_parameters,
                            Tm_chem_correction_parameters=forward_primer_Tm_chem_correction_param,
                            Tm_salt_correction_parameters=forward_primer_Tm_salt_correction_parameters,
                        )
                        dif_Tm = abs(Tm_forward_primer - Tm_reverse_primer)
                        if dif_Tm < min_dif_Tm:
                            min_dif_Tm = dif_Tm
                            forward_primer_sequence = primer_attributes["oligo"]
                return forward_primer_sequence

            pipeline = SeqFishPlusPrimerDesigner(
                write_intermediate_steps=self.write_intermediate_steps,
                dir_output=self.dir_output,
                n_jobs=self.n_jobs,
            )
            primer_database, _ = pipeline.create_oligo_database(
                oligo_length=forward_primer_length,
                oligo_base_probabilities=forward_primer_base_probabilities,
                initial_num_sequences=forward_primer_initial_num_sequences,
            )

            primer_database, _ = pipeline.filter_by_property(
                oligo_database=primer_database,
                GC_content_min=forward_primer_GC_content_min,
                GC_content_max=forward_primer_GC_content_max,
                number_GC_GCclamp=forward_primer_number_GC_GCclamp,
                number_three_prime_base_GCclamp=forward_primer_number_three_prime_base_GCclamp,
                homopolymeric_base_n=forward_primer_homopolymeric_base_n,
                max_len_selfcomplement=forward_primer_max_len_selfcomplement,
                reverse_primer_sequence=reverse_primer_sequence,
                max_len_complement=forward_primer_max_len_complement,
                Tm_min=forward_primer_Tm_min,
                Tm_max=forward_primer_Tm_max,
                Tm_parameters=forward_primer_Tm_parameters,
                Tm_chem_correction_param=forward_primer_Tm_chem_correction_param,
                Tm_salt_correction_parameters=forward_primer_Tm_salt_correction_parameters,
                T_secondary_structure=forward_primer_T_secondary_structure,
                secondary_structures_threshold_deltaG=forward_primer_secondary_structures_threshold_deltaG,
            )

            primer_database, _ = pipeline.filter_by_specificity(
                oligo_database=primer_database,
                files_fasta_reference_database=forward_primer_files_fasta_reference_database,
                specificity_refrence_blastn_search_parameters=forward_primer_specificity_refrence_blastn_search_params,
                specificity_refrence_blastn_hit_parameters=forward_primer_specificity_refrence_blastn_hit_params,
                file_fasta_encoding_probes_database=file_fasta_encoding_probes_database,
                specificity_encoding_probes_blastn_search_parameters=forward_primer_specificity_encoding_probes_blastn_search_parameters,
                specificity_encoding_probes_blastn_hit_parameters=forward_primer_specificity_encoding_probes_blastn_hit_parameters,
            )

            forward_primer_sequence = _get_forward_primer_opt()

            return reverse_primer_sequence, forward_primer_sequence

        def _add_primer_sequences(oligo_database, reverse_primer_sequence, forward_primer_sequence):
            probe_ids = oligo_database.get_oligoid_list()
            primer_attribute = {
                probe_id: {
                    "sequence_reverse_primer": reverse_primer_sequence,
                    "sequence_forward_primer": forward_primer_sequence,
                }
                for probe_id in probe_ids
            }
            oligo_database.update_oligo_attributes(primer_attribute)

            for region_id, database_region in oligo_database.database.items():
                for probe_id, probe_attributes in database_region.items():
                    oligo_database.database[region_id][probe_id]["sequence_seqfish_plus_probe"] = (
                        probe_attributes["sequence_forward_primer"]
                        + "T"
                        + probe_attributes["sequence_encoding_probe"]
                        + probe_attributes["sequence_reverse_primer"]
                    )
            return oligo_database

        region_ids = list(oligo_database.database.keys())

        codebook, readout_probe_table = _generate_readout_probes(n_genes=len(oligo_database.database))
        # codebook = codebook.iloc[: len(region_ids)]
        codebook.index = region_ids + [
            f"unassigned_barcode_{i+1}" for i in range(len(codebook.index) - len(region_ids))
        ]

        codebook.to_csv(os.path.join(self.dir_output, "codebook.tsv"), sep="\t")
        readout_probe_table.to_csv(os.path.join(self.dir_output, "readout_probes.tsv"), sep="\t")

        with joblib_progress(description="Design Encoding Probe Sequence", total=len(region_ids)):
            Parallel(
                n_jobs=self.n_jobs, prefer="threads", require="sharedmem"
            )(  # there should be an explicit return
                delayed(_assemble_encoding_sequence)(
                    oligo_database, region_id, codebook.loc[region_id], readout_probe_table
                )
                for region_id in region_ids
            )

        file_fasta_encoding_probes_database = oligo_database.write_database_to_fasta(
            filename=f"db_reference_encoding_probes",
            region_ids=None,
            sequence_type="sequence_encoding_probe",
        )

        reverse_primer_sequence, forward_primer_sequence = _generate_primers(
            file_fasta_encoding_probes_database=file_fasta_encoding_probes_database
        )

        os.remove(file_fasta_encoding_probes_database)

        oligo_database = _add_primer_sequences(
            oligo_database, reverse_primer_sequence, forward_primer_sequence
        )

        return oligo_database

    def compute_probe_attributes(
        self,
        oligo_database: OligoDatabase,
        T_secondary_structure: float,
    ):
        oligo_database = self.probe_attributes_calculator.calculate_oligo_length(
            oligo_database=oligo_database
        )
        oligo_database = self.probe_attributes_calculator.calculate_GC_content(
            oligo_database=oligo_database, sequence_type="oligo"
        )
        oligo_database = self.probe_attributes_calculator.calculate_num_targeted_transcripts(
            oligo_database=oligo_database
        )
        oligo_database = self.probe_attributes_calculator.calculate_isoform_consensus(
            oligo_database=oligo_database
        )
        oligo_database = self.probe_attributes_calculator.calculate_DG_secondary_structure(
            oligo_database=oligo_database, sequence_type="oligo", T=T_secondary_structure
        )

        return oligo_database

    def generate_output(self, oligo_database: OligoDatabase, top_n_sets: int):

        attributes = [
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
            "isoform_consensus",
            "DG_secondary_structure",
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
            filename="seqfish_plus_probes.yml",
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
                        "sequence_seqfish_plus_probe": database_region[oligo_id][
                            "sequence_seqfish_plus_probe"
                        ],
                        "sequence_readout_probe_1": database_region[oligo_id]["sequence_readout_probe_1"],
                        "sequence_readout_probe_2": database_region[oligo_id]["sequence_readout_probe_2"],
                        "sequence_readout_probe_3": database_region[oligo_id]["sequence_readout_probe_3"],
                        "sequence_readout_probe_4": database_region[oligo_id]["sequence_readout_probe_4"],
                    }

        with open(os.path.join(self.dir_output, "padlock_probes_order.yml"), "w") as outfile:
            yaml.dump(yaml_dict_order, outfile, default_flow_style=False, sort_keys=False)


############################################
# SeqFish Plus Readout Probe Designer
############################################


class SeqFishPlusReadoutProbeDesigner:
    def __init__(
        self,
        write_intermediate_steps: bool,
        dir_output: str,
        n_jobs: int,
    ) -> None:
        """Constructor for the SeqFishPlusReadoutProbeDesigner class."""

        self.write_intermediate_steps = write_intermediate_steps
        self.n_jobs = n_jobs

        ##### create the output folder #####
        self.dir_output = os.path.abspath(dir_output)
        Path(self.dir_output).mkdir(parents=True, exist_ok=True)

        self.subdir_db_probes = "db_readout_probes"
        self.subdir_db_reference = "db_reference"

    @pipeline_step_basic(step_name="Readout Probe Generation - Create Oligo Database")
    def create_oligo_database(
        self,
        oligo_length: int,
        oligo_base_probabilities: dict,
        initial_num_sequences: int,
    ):
        ##### creating the oligo sequences #####
        probe_sequences = OligoSequenceGenerator(dir_output=self.dir_output)
        probe_fasta_file = probe_sequences.create_sequences_random(
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
            database_name=self.subdir_db_probes,
            dir_output=self.dir_output,
            n_jobs=1,
        )
        oligo_database.load_database_from_fasta(
            files_fasta=probe_fasta_file,
            sequence_type="oligo",
        )

        ##### save database #####
        if self.write_intermediate_steps:
            file_database = oligo_database.save_database(dir_database="1_db_readout_probes_initial")
        else:
            file_database = ""

        dir = probe_sequences.dir_output
        shutil.rmtree(dir) if os.path.exists(dir) else None

        return oligo_database, file_database

    @pipeline_step_basic(step_name="Readout Probe Generation - Property Filters")
    def filter_by_property(
        self,
        oligo_database: OligoDatabase,
        GC_content_min: float,
        GC_content_max: float,
        homopolymeric_base_n: int,
    ):
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
            file_database = oligo_database.save_database(dir_database="2_db_readout_probes_property_filter")
        else:
            file_database = ""

        return oligo_database, file_database

    @pipeline_step_basic(step_name="Readout Probe Generation - Specificity Filters")
    def filter_by_specificity(
        self,
        oligo_database: OligoDatabase,
        files_fasta_reference_database: List[str],
        specificity_blastn_search_parameters: dict,
        specificity_blastn_hit_parameters: dict,
        cross_hybridization_blastn_search_parameters: dict,
        cross_hybridization_blastn_hit_parameters: dict,
    ):
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

        # write the intermediate result in a file
        if self.write_intermediate_steps:
            file_database = oligo_database.save_database(dir_database="3_db_readout_probes_specificty_filter")
        else:
            file_database = ""

        for directory in [
            reference_database.dir_output,
            specificity.dir_output,
            cross_hybridization_aligner.dir_output,
            cross_hybridization.dir_output,
        ]:
            if os.path.exists(directory):
                shutil.rmtree(directory)

        return oligo_database, file_database

    def generate_codebook(self, n_regions: int, n_barcode_rounds: int, n_pseudocolors: int, n_channels: int):
        """This function generates the codebook containig a collection of barcodes originated form the given barcode rounds, pseudocolors and channels.
        Specifically, each barcode is generated such that all the readout probes encoding for one region is associated to the same channel.
        To allow for error corrections the barcodes are generated to be robust against a sigle deletion. To do so, the pseudocolor assignoed to the last barcode round
        is deterministically derived from the other pseoudocolors.

        :param n_regions: Number of regions that need to be assigned to a barcode.
        :type n_regions: int
        :param n_barcode_rounds: Number of barcode runs, defaults to 4.
        :type n_barcode_rounds: int, optional
        :param n_pseudocolors: Number of pseudocolors, defaults to 20.
        :type n_pseudocolors: int, optional
        :param n_channels: Number of channels, defaults to 3.
        :type n_channels: int, optional
        :return: Collection of all the valid barcodes generated
        :rtype: list
        """

        def _generate_barcode(pseudocolors: list, channel: int, n_pseudocolors: int, n_channels: int):
            """This function generates a channel-wise binary barcode, i.e. barcodes where each pseudocolor belongs to the same channel.
            To allow for error corrections an additional barcode round is added and the barcodes are generated to be robust against a sigle deletion.
            For example for a barcode with (i,j,k) pseudocolors, to the additional barcode round we will assig the `(i+j+k) mod n_pseudocolors` pseoudocolor.

            :param pseudocolors: List of the pseoudocolors contained in the barcode.
            :type pseudocolors: list
            :param channel: channel to wich the barcode belongs to.
            :type channel: int
            :param n_pseudocolors: Total number of pseudocolors.
            :type n_pseudocolors: int
            :param n_channels: Total number fo channels
            :type n_channels: int
            :return: Binary barcode.
            :rtype: np.Array
            """
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
    ):
        """This function associates the readout probes with the bit of the barcodes and assigns them to a specific barcoding round, pseudocolor and channel.

        :param readout_probes: Oligo databse containing all the readout probes.
        :type readout_probes: OligoDatabase
        :param channels_ids: Names of the available channles.
        :type channels_ids: list
        :param n_barcode_rounds: Number of barcode runs, defaults to 4
        :type n_barcode_rounds: int, optional
        :param n_pseudocolors: Number of pseudocolors, defaults to 20
        :type n_pseudocolors: int, optional
        :return: table containig the information for each readout probe.
        :rtype: pd.DataFrame
        """
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
        write_intermediate_steps: bool,
        dir_output: str,
        n_jobs: int,
    ) -> None:
        """Constructor for the SeqFishPlusPrimerDesigner class."""

        self.write_intermediate_steps = write_intermediate_steps
        self.n_jobs = n_jobs

        ##### create the output folder #####
        self.dir_output = os.path.abspath(dir_output)
        Path(self.dir_output).mkdir(parents=True, exist_ok=True)

        self.subdir_db_primers = "db_primer"
        self.subdir_db_reference = "db_reference"

    @pipeline_step_basic(step_name="Primer Generation - Create Oligo Database")
    def create_oligo_database(
        self,
        oligo_length: int,
        oligo_base_probabilities: dict,
        initial_num_sequences: int,
    ):
        ##### creating the primer sequences #####
        # random forward primer
        forward_primer_sequences = OligoSequenceGenerator(dir_output=self.dir_output)
        forward_primer_fasta_file = forward_primer_sequences.create_sequences_random(
            filename_out="forward_primer_sequences",
            length_sequences=oligo_length,
            num_sequences=initial_num_sequences,
            name_sequences="forward_primer",
            base_alphabet_with_probability=oligo_base_probabilities,
        )

        ##### creating the primer database #####
        oligo_database = OligoDatabase(
            min_oligos_per_region=0,
            write_regions_with_insufficient_oligos=False,
            lru_db_max_in_memory=self.n_jobs * 2 + 2,
            database_name=self.subdir_db_primers,
            dir_output=self.dir_output,
            n_jobs=1,
        )
        oligo_database.load_database_from_fasta(
            files_fasta=forward_primer_fasta_file,
            sequence_type="oligo",
        )

        ##### save database #####
        if self.write_intermediate_steps:
            file_database = oligo_database.save_database(dir_database="1_db_primers_initial")
        else:
            file_database = ""

        dir = forward_primer_sequences.dir_output
        shutil.rmtree(dir) if os.path.exists(dir) else None

        return oligo_database, file_database

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
        Tm_chem_correction_param: dict,
        Tm_salt_correction_parameters: dict,
        T_secondary_structure: float,
        secondary_structures_threshold_deltaG: float,
    ):
        # define the filters
        # we want to keep primer which end with a specific nucleotide, i.e. "T"
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
            Tm_chem_correction_parameters=Tm_chem_correction_param,
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
            sequence_type="oligo",
            oligo_database=oligo_database,
            n_jobs=self.n_jobs,
        )

        # write the intermediate result in a file
        if self.write_intermediate_steps:
            file_database = oligo_database.save_database(dir_database="2_db_primer_property_filter")
        else:
            file_database = ""

        return oligo_database, file_database

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
    ):
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

        # write the intermediate result in a file
        if self.write_intermediate_steps:
            file_database = oligo_database.save_database(dir_database="3_db_primer_specificty_filter")
        else:
            file_database = ""

        for directory in [
            reference_database.dir_output,
            encoding_probes_database.dir_output,
            specificity_refrence.dir_output,
            specificity_encoding_probes.dir_output,
        ]:
            if os.path.exists(directory):
                shutil.rmtree(directory)

        return oligo_database, file_database


############################################
# SeqFish Plus Probe Designer Pipeline
############################################


def main():
    args = base_parser()

    ##### read the config file #####
    with open(args["config"], "r") as handle:
        config = yaml.safe_load(handle)

    ##### process parameters #####
    # preprocess melting temperature params
    Tm_parameters_primer = config["Tm_parameters_primer"]
    Tm_parameters_primer["nn_table"] = getattr(mt, Tm_parameters_primer["nn_table"])
    Tm_parameters_primer["tmm_table"] = getattr(mt, Tm_parameters_primer["tmm_table"])
    Tm_parameters_primer["imm_table"] = getattr(mt, Tm_parameters_primer["imm_table"])
    Tm_parameters_primer["de_table"] = getattr(mt, Tm_parameters_primer["de_table"])

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

    ##### create probe database #####
    probe_database, file_database = pipeline.create_oligo_database(
        gene_ids=gene_ids,
        oligo_length_min=config["probe_length_min"],
        oligo_length_max=config["probe_length_max"],
        files_fasta_oligo_database=config["files_fasta_probe_database"],
        # we should have at least "min_probeset_size" probes per gene to create one set
        min_probes_per_gene=config["probeset_size_min"],
        isoform_consensus=config["probe_isoform_consensus"],
    )

    ##### filter oligos by property #####
    probe_database, file_database = pipeline.filter_by_property(
        oligo_database=probe_database,
        GC_content_min=config["probe_GC_content_min"],
        GC_content_max=config["probe_GC_content_max"],
        homopolymeric_base_n=config["probe_homopolymeric_base_n"],
        T_secondary_structure=config["probe_T_secondary_structure"],
        secondary_structures_threshold_deltaG=config["probe_secondary_structures_threshold_deltaG"],
    )

    # ##### filter probes by specificity #####
    probe_database, file_database = pipeline.filter_by_specificity(
        oligo_database=probe_database,
        files_fasta_reference_database=config["probe_files_fasta_reference_database"],
        specificity_blastn_search_params=config["probe_specificity_blastn_search_params"],
        specificity_blastn_hit_params=config["probe_specificity_blastn_hit_params"],
        cross_hybridization_blastn_search_params=config["probe_cross_hybridization_blastn_search_params"],
        cross_hybridization_blastn_hit_params=config["probe_cross_hybridization_blastn_hit_params"],
    )

    ##### create probe sets #####
    probe_database, file_database, dir_probesets = pipeline.create_oligo_sets(
        oligo_database=probe_database,
        GC_weight=config["probe_GC_weight"],
        GC_content_opt=config["probe_GC_content_opt"],
        UTR_weight=config["probe_UTR_weight"],
        set_size_opt=config["probeset_size_opt"],
        set_size_min=config["probeset_size_min"],
        max_graph_size=config["max_graph_size"],
        n_sets=config["n_sets"],
        distance_between_oligos=config["distance_between_probes"],
    )

    probe_database = pipeline.design_final_probe_sequence(
        oligo_database=probe_database,
        readout_probe_length=config["readout_probe_length"],
        readout_probe_base_probabilities=config["readout_probe_base_probabilities"],
        readout_probe_initial_num_sequences=config["readout_probe_initial_num_sequences"],
        readout_probe_GC_content_min=config["readout_probe_GC_content_min"],
        readout_probe_GC_content_max=config["readout_probe_GC_content_max"],
        readout_probe_homopolymeric_base_n=config["readout_probe_homopolymeric_base_n"],
        readout_probe_files_fasta_reference_database=config["readout_probe_files_fasta_reference_database"],
        readout_probe_specificity_blastn_search_params=config[
            "readout_probe_specificity_blastn_search_params"
        ],
        readout_probe_specificity_blastn_hit_params=config["readout_probe_specificity_blastn_hit_params"],
        readout_probe_cross_hybridization_blastn_search_params=config[
            "readout_probe_cross_hybridization_blastn_search_params"
        ],
        readout_probe_cross_hybridization_blastn_hit_params=config[
            "readout_probe_cross_hybridization_blastn_hit_params"
        ],
        n_barcode_rounds=config["n_barcode_rounds"],
        n_pseudocolors=config["n_pseudocolors"],
        channels_ids=config["channels_ids"],
        forward_primer_length=config["forward_primer_length"],
        forward_primer_base_probabilities=config["forward_primer_base_probabilities"],
        forward_primer_initial_num_sequences=config["forward_primer_initial_num_sequences"],
        reverse_primer_sequence=config["reverse_primer_sequence"],
        forward_primer_GC_content_min=config["primer_GC_content_min"],
        forward_primer_GC_content_max=config["primer_GC_content_max"],
        forward_primer_number_GC_GCclamp=config["primer_number_GC_GCclamp"],
        forward_primer_number_three_prime_base_GCclamp=config["primer_number_three_prime_base_GCclamp"],
        forward_primer_homopolymeric_base_n=config["primer_homopolymeric_base_n"],
        forward_primer_max_len_selfcomplement=config["primer_max_len_selfcomplement"],
        forward_primer_max_len_complement=config["primer_max_len_complement_reverse_primer"],
        forward_primer_Tm_min=config["primer_Tm_min"],
        forward_primer_Tm_max=config["primer_Tm_max"],
        forward_primer_Tm_parameters=Tm_parameters_primer,
        forward_primer_Tm_chem_correction_param=config["Tm_chem_correction_param_primer"],
        forward_primer_Tm_salt_correction_parameters=config["Tm_salt_correction_param_primer"],
        forward_primer_T_secondary_structure=config["primer_T_secondary_structure"],
        forward_primer_secondary_structures_threshold_deltaG=config[
            "primer_secondary_structures_threshold_deltaG"
        ],
        forward_primer_files_fasta_reference_database=config["forward_primer_files_fasta_reference_database"],
        forward_primer_specificity_refrence_blastn_search_params=config[
            "forward_primer_specificity_refrence_blastn_search_params"
        ],
        forward_primer_specificity_refrence_blastn_hit_params=config[
            "forward_primer_specificity_refrence_blastn_hit_params"
        ],
        forward_primer_specificity_encoding_probes_blastn_search_parameters=config[
            "forward_primer_specificity_encoding_probes_blastn_search_parameters"
        ],
        forward_primer_specificity_encoding_probes_blastn_hit_parameters=config[
            "forward_primer_specificity_encoding_probes_blastn_hit_parameters"
        ],
    )

    ##### generate output #####
    # compute all required attributes
    probe_database = pipeline.compute_probe_attributes(
        oligo_database=probe_database,
        T_secondary_structure=config["probe_T_secondary_structure"],
    )

    pipeline.generate_output(oligo_database=probe_database, top_n_sets=config["top_n_sets"])

    logging.info(f"Oligo sets were saved in {dir_probesets}")
    logging.info("##### End of the pipeline. #####")
