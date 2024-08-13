############################################
# imports
############################################
import logging
import os
import shutil
import warnings
from datetime import datetime
from itertools import combinations
from pathlib import Path
from typing import List

import numpy as np
import pandas as pd
import yaml
from Bio.SeqUtils import MeltingTemp as mt
from scipy.spatial.distance import hamming

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
from oligo_designer_toolsuite.pipelines._utils import base_parser, pipeline_step_basic
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
        if not file_regions:
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

    @pipeline_step_basic(step_name="Create Database")
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
        oligo_database.load_database_from_fasta(
            files_fasta=oligo_fasta_file,
            sequence_type="target",
            region_ids=self.gene_ids,
        )

        ##### save database #####
        if self.write_intermediate_steps:
            file_database = oligo_database.save_database(dir_database="1_db_initial")
        else:
            file_database = ""

        return oligo_database, file_database

    @pipeline_step_basic(step_name="Property Filters")
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
            file_database = oligo_database.save_database(dir_database="2_db_property_filter")
        else:
            file_database = ""

        return oligo_database, file_database

    @pipeline_step_basic(step_name="Specificty Filters")
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
        reference_database.load_database_from_fasta(
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
            file_database = oligo_database.save_database(dir_database="3_db_specificty_filter")
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
            file_database = oligo_database.save_database(dir_database="4_db_probes_probesets")
            file_probesets = oligo_database.write_oligosets_to_table()
        else:
            file_database = ""
            file_probesets = ""

        return oligo_database, file_database, file_probesets

    def design_final_probe_sequences(
        self,
        oligo_database: OligoDatabase,
        codebook,
        readout_probe_table,
        forward_primer_sequence,
        reverse_primer_sequence,
    ):
        region_ids = oligo_database.database.keys()

        with joblib_progress(description="Design Encoding Probe Sequence", total=len(region_ids)):
            Parallel(
                n_jobs=self.n_jobs, prefer="threads", require="sharedmem"
            )(  # there should be an explicit return
                delayed(_assemble_final_sequence)(
                    oligo_database, region_id, codebook.loc[region_id], readout_probe_table
                )
                for region_id in region_ids
            )

        def _assemble_final_sequence(
            oligo_database: OligoDatabase,
            region_id: str,
            barcode: pd.Series,
            readout_probe_table: pd.DataFrame,
        ):
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
                    + "A"
                    + probe_attributes["oligo"]
                    + "A"
                    + str(Seq(probe_attributes["sequence_readout_probe_2"]).reverse_complement())
                )

                probe_attributes["sequence_forward_primer"] = forward_primer_sequence
                probe_attributes["sequence_reverse_primer"] = reverse_primer_sequence
                probe_attributes["sequence_final_probe"] = (
                    probe_attributes["sequence_forward_primer"]
                    + "T"
                    + probe_attributes["sequence_encoding_probe"]
                    + probe_attributes["sequence_reverse_primer"]
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
            filename="merfish_probes.yml",
        )

        # write a second file that only contains order information
        yaml_dict_order = {}

        for region_id, database_region in oligo_database.database.items():
            yaml_dict_order[region_id] = {}
            oligosets_region = oligo_database.oligosets[region_id]
            oligosets_oligo_columns = [col for col in oligosets_region.columns if col.startswith("oligo_")]
            oligosets_score_columns = [col for col in oligosets_region.columns if col.startswith("score_")]

            oligosets_region = oligosets_region.sort_values(by=oligosets_score_columns, ascending=True)
            oligosets_region = oligosets_region.head(top_n_sets)[oligosets_oligo_columns]
            oligosets_region.reset_index(inplace=True, drop=True)

            # iterate through all oligo sets
            for oligoset_idx, oligoset in oligosets_region.iterrows():
                oligoset_id = f"oligoset_{oligoset_idx + 1}"
                yaml_dict_order[region_id][oligoset_id] = {}
                for oligo_id in oligoset:
                    yaml_dict_order[region_id][oligoset_id][oligo_id] = {
                        "sequence_merfish_probe": database_region[oligo_id]["sequence_merfish_probe"],
                        "sequence_readout_probe_1": database_region[oligo_id]["sequence_readout_probe_1"],
                        "sequence_readout_probe_2": database_region[oligo_id]["sequence_readout_probe_2"],
                    }

        with open(os.path.join(self.dir_output, "merfish_probes_order.yml"), "w") as outfile:
            yaml.dump(yaml_dict_order, outfile, default_flow_style=False, sort_keys=False)


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

    @pipeline_step_basic(step_name="Create Readout Probe Database")
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
        oligo_database.load_database_from_fasta(
            files_fasta=[oligo_fasta_file],
            sequence_type="target",
        )

        ##### save database #####
        if self.write_intermediate_steps:
            file_database = oligo_database.save_database(dir_database="1_db_initial")
        else:
            file_database = ""

        return oligo_database, file_database

    @pipeline_step_basic(step_name="Property Filters")
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
            file_database = oligo_database.save_database(dir_database="2_db_property_filter")
        else:
            file_database = ""

        return oligo_database, file_database

    @pipeline_step_basic(step_name="Specificty Filters")
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
        reference_database.load_database_from_fasta(
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
            file_database = oligo_database.save_database(dir_database="3_db_specificty_filter")
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

    def make_barcode(self, raw_barcode: list, n_bits: int):
        """Creates the actual barcode from a list containing the indices of the "one" bits.
        :param raw_barcode: list of the "one" bits in the barcode.
        :type raw_barcode: list
        :param n_bits: Number of bits contained in each barcode.
        :type n_bits: int, optional
        """
        barcode = np.zeros(n_bits, dtype=np.int8)
        for i in raw_barcode:
            barcode[i] = 1
        return barcode

    def generate_codebook(
        self,
        n_regions: int,
        n_bits: int = 16,
        min_hamming_dist: int = 4,
        hamming_weight: int = 4,
    ):
        """This function generates the codebook containing a collection of valid barcodes that fulfill the following constraints:
        - The barcodes have all the same `hamming_weight`, i.e. the number of bits containing 1.
        - The Hamming distance between every pair of valid barcodes is higher than the given threshold `min_hamming_dist`.
        The default parameters follow the MHD4 standard.

        :param n_regions: Number of regions that need to be assigned to a barcode.
        :type n_regions: int
        :param n_bits: Number of bits contained in each barcode, defaults to 16
        :type n_bits: int, optional
        :param min_hamming_dist: Minimum distance between two valid barcodes, defaults to 4
        :type min_hamming_dist: int, optional
        :param hamming_weight: Number of bits containing one in each barcode, defaults to 4
        :type hamming_weight: int, optional
        :return: Collection of all the valid barcodes generated
        :rtype: list
        """
        codebook = []
        for raw_barcode in combinations(iterable=range(n_bits), r=hamming_weight):
            new_barcode = self.make_barcode(raw_barcode=raw_barcode, n_bits=n_bits)
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

    def create_readout_probes_table(
        self, readout_probes_database: OligoDatabase, channels_ids: list, n_bits: int
    ):

        readout_probes = readout_probes_database.get_oligoid_sequence_mapping(
            sequence_type="oligo", sequence_to_upper=False
        )
        assert (
            len(readout_probes.database) >= n_bits
        ), f"There are less readout probes ({len(readout_probes.database)}) than bits ({n_bits})."
        table = pd.DataFrame(
            columns=["bit", "channel", "readout_probe_id", "readout_probe_sequence"],
            index=list(range(n_bits)),
        )
        n_channels = len(channels_ids)
        channel = 0
        for i, (readout_probe_id, readout_probe_features) in enumerate(readout_probes.database.items()):
            table.iloc[i] = [
                i,
                channels_ids[channel],
                readout_probe_id,
                readout_probe_features["oligo"],
            ]
            channel = (channel + 1) % n_channels
            if i >= n_bits - 1:
                break

        return table

    # Higher level function to run the whole pipeline, the steps could also be run individually
    def get_readout_probes(
        self,
        oligo_length: int,
        readout_sequence_probs: dict,
        initial_num_sequences: int,
        GC_content_min: int,
        GC_content_max: int,
        homopolymeric_base_n: int,
        files_fasta_reference_database: List[str],
        blastn_search_parameters: dict,
        blastn_hit_parameters: dict,
        cross_hybridization_blastn_search_parameters: dict,
        cross_hybridization_blastn_hit_parameters: dict,
        region_ids: List[str],
        n_bits: int,
        min_hamming_dist: int,
        hamming_weight: int,
        channels_ids: list,
    ):
        readout_probes, file_database = self.create_oligo_database(
            oligo_length=oligo_length,
            readout_sequence_probs=readout_sequence_probs,
            initial_num_sequences=initial_num_sequences,
        )

        readout_probes, file_database = self.filter_by_property(
            oligo_database=readout_probes,
            GC_content_min=GC_content_min,
            GC_content_max=GC_content_max,
            homopolymeric_base_n=homopolymeric_base_n,
        )

        readout_probes, file_database = self.filter_by_specificity(
            oligo_database=readout_probes,
            files_fasta_reference_database=files_fasta_reference_database,
            blastn_search_parameters=blastn_search_parameters,
            blastn_hit_parameters=blastn_hit_parameters,
            cross_hybridization_blastn_search_parameters=cross_hybridization_blastn_search_parameters,
            cross_hybridization_blastn_hit_parameters=cross_hybridization_blastn_hit_parameters,
        )

        codebook = self.generate_codebook(
            n_regions=len(region_ids),
            n_bits=n_bits,
            min_hamming_dist=min_hamming_dist,
            hamming_weight=hamming_weight,
        )

        readout_probe_table = self.create_readout_probes_table(
            readout_probes=readout_probes,
            channels_ids=channels_ids,
            n_bits=n_bits,
        )

        codebook.index = region_ids + [
            f"unassigned_barcode_{i+1}" for i in range(len(codebook.index) - len(region_ids))
        ]

        codebook.to_csv(os.path.join(self.dir_output, "codebook.tsv"), sep="\t")
        readout_probe_table.to_csv(os.path.join(self.dir_output, "readout_probes.tsv"), sep="\t")

        return codebook, readout_probe_table


############################################
# MERFISH Primer Designer Functions
############################################


class MerfishPrimerDesigner:

    # Pipeline for the forward primer (as in the table above)
    def __init__(
        self,
        write_intermediate_steps: bool,
        dir_output: str,
        n_jobs: int,
    ) -> None:
        """Constructor for the MerfishPrimerDesigner class."""

        self.write_intermediate_steps = write_intermediate_steps
        self.n_jobs = n_jobs

        self.dir_output = os.path.abspath(dir_output)
        Path(self.dir_output).mkdir(parents=True, exist_ok=True)

        self.subdir_db_primers = "db_primer"
        self.subdir_db_reference = "db_reference"

        self.oligo_attributes_calculator = OligoAttributes()

    # Higher level function to run the whole pipeline, the steps could also be run individually
    def get_forward_primer(
        self,
        oligo_length: int,
        oligo_base_probabilities: dict,
        initial_num_sequences: int,
        gc_content_min: int,
        gc_content_max: int,
        gc_clamp_n_bases: int,
        gc_clamp_nGC: int,
        tm_min: int,
        tm_max: int,
        tm_parameters: dict,
        homopolymeric_base_n: int,
        files_fasta_reference_database_transcriptome: List[str],
        blastn_search_parameters_transcriptome: dict,
        blastn_hit_parameters_transcriptome: dict,
        files_fasta_reference_database_encoding_probe: List[str],
        blastn_search_parameters_encoding_probe: dict,
        blastn_hit_parameters_encoding_probe: dict,
        tm_reverse_primer: int,
    ):

        while True:
            oligo_database, file_database = self.create_oligo_database(
                oligo_length=oligo_length,
                oligo_base_probabilities=oligo_base_probabilities,
                initial_num_sequences=initial_num_sequences,
            )

            oligo_database, file_database = self.filter_by_property(
                oligo_database=oligo_database,
                GC_content_min=gc_content_min,
                GC_content_max=gc_content_max,
                GC_clamp_n_bases=gc_clamp_n_bases,
                GC_clamp_nGC=gc_clamp_nGC,
                Tm_min=tm_min,
                Tm_max=tm_max,
                Tm_parameters=tm_parameters,
                homopolymeric_base_n=homopolymeric_base_n,
            )

            oligo_database, file_database = self.filter_by_specificity(
                oligo_database=oligo_database,
                files_fasta_reference_database_transcriptome=files_fasta_reference_database_transcriptome,
                blastn_search_parameters_transcriptome=blastn_search_parameters_transcriptome,
                blastn_hit_parameters_transcriptome=blastn_hit_parameters_transcriptome,
                files_fasta_reference_database_encoding_probe=files_fasta_reference_database_encoding_probe,
                blastn_search_parameters_encoding_probe=blastn_search_parameters_encoding_probe,
                blastn_hit_parameters_encoding_probe=blastn_hit_parameters_encoding_probe,
            )
            if len(oligo_database.database) > 0:
                first_region = next(oligo_database.database.values())
                first_oligo = list(first_region.values())[0]
                best_oligo = first_oligo
                for region in oligo_database.database.values():
                    for oligo in region.values():
                        if abs(oligo["TmNN"] - tm_reverse_primer) < abs(
                            best_oligo["TmNN"] - tm_reverse_primer
                        ):
                            best_oligo = oligo

                return best_oligo["oligo"]

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
        GC_content_min: int,
        GC_content_max: int,
        GC_clamp_n_bases: int,
        GC_clamp_nGC: int,
        Tm_min: int,
        Tm_max: int,
        Tm_parameters: dict,
        homopolymeric_base_n: int,
    ):

        gc_content = GCContentFilter(GC_content_min=GC_content_min, GC_content_max=GC_content_max)
        gc_clamp = GCClampFilter(n_bases=GC_clamp_n_bases, n_GC=GC_clamp_nGC)
        melting_temperature = MeltingTemperatureNNFilter(
            Tm_min=Tm_min,
            Tm_max=Tm_max,
            Tm_parameters=Tm_parameters,
        )
        homopolymeric_runs = HomopolymericRunsFilter(
            base_n=homopolymeric_base_n,
        )

        filters = [
            gc_content,
            gc_clamp,
            melting_temperature,
            homopolymeric_runs,
        ]

        property_filter = PropertyFilter(filters=filters)

        oligo_database = property_filter.apply(
            sequence_type="oligo",
            oligo_database=oligo_database,
            n_jobs=self.n_jobs,
        )

        if self.write_intermediate_steps:
            file_database = oligo_database.save_database(dir_database="2_db_primers_property_filter")
        else:
            file_database = ""

        return oligo_database, file_database

    @pipeline_step_basic(step_name="Primer Generation - Specificty Filters")
    def filter_by_specificity(
        self,
        oligo_database: OligoDatabase,
        files_fasta_reference_database_transcriptome: List[str],
        blastn_search_parameters_transcriptome: dict,
        blastn_hit_parameters_transcriptome: dict,
        files_fasta_reference_database_encoding_probe: List[str],
        blastn_search_parameters_encoding_probe: dict,
        blastn_hit_parameters_encoding_probe: dict,
    ):
        reference_database_transcriptome = ReferenceDatabase(
            database_name=self.subdir_db_reference, dir_output=self.dir_output
        )
        reference_database_transcriptome.load_database_from_fasta(
            files_fasta=files_fasta_reference_database_transcriptome,
            database_overwrite=False,
        )

        # BlastN filter
        blastn_filter_transcriptome = BlastNFilter(
            search_parameters=blastn_search_parameters_transcriptome,
            hit_parameters=blastn_hit_parameters_transcriptome,
        )

        reference_database_encoding_probe = ReferenceDatabase(
            database_name=self.subdir_db_reference, dir_output=self.dir_output
        )
        reference_database_encoding_probe.load_database_from_fasta(
            files_fasta=files_fasta_reference_database_encoding_probe,
            database_overwrite=False,
        )

        # BlastN filter
        blastn_filter_encoding_probe = BlastNFilter(
            search_parameters=blastn_search_parameters_encoding_probe,
            hit_parameters=blastn_hit_parameters_encoding_probe,
        )

        filters = [blastn_filter_transcriptome, blastn_filter_encoding_probe]
        specificity_filter = SpecificityFilter(filters=filters)
        oligo_database = specificity_filter.apply(
            sequence_type="oligo",
            oligo_database=oligo_database,
            reference_database=reference_database_transcriptome,
            n_jobs=self.n_jobs,
        )

        if self.write_intermediate_steps:
            file_database = oligo_database.save_database(dir_database="3_db_primers_specificty_filter")
        else:
            file_database = ""

        for directory in [
            reference_database_transcriptome.dir_output,
            blastn_filter_transcriptome.dir_output,
        ]:
            if os.path.exists(directory):
                shutil.rmtree(directory)

        return oligo_database, file_database


############################################
# MERFISH Probe Designer Pipeline
############################################


def main():
    """Main function to execute the oligo probe design pipeline based on user configurations.
    It parses command line arguments, initializes the pipeline with specified parameters, creates an oligo database,
    applies various filtering criteria, computes necessary oligo attributes, and finally generates oligo sets.
    """

    args = base_parser()
    # read config file
    with open(args["config"], "r") as handle:
        config = yaml.safe_load(handle)

    #### Target sequences pipeline ####
    def run_target_probe_pipeline(config):
        pipeline = MerfishProbeDesigner(
            file_regions=config["file_regions"],
            write_intermediate_steps=config["write_intermediate_steps"],
            dir_output=config["dir_output"],
            n_jobs=config["n_jobs"],
        )

        # create the oligo database
        oligo_database, file_database = pipeline.create_oligo_database(
            oligo_length=config["probe_length"],
            files_fasta_oligo_database=config["files_fasta_probe_database"],
            # we should have at least "min_oligoset_size" oligos per gene to create one set
            min_oligos_per_region=config["probeset_size"],
        )

        # preprocess the Tm parameters
        probe_Tm_parameters = config["probe_Tm_parameters"]
        probe_Tm_parameters["nn_table"] = getattr(mt, probe_Tm_parameters["nn_table"])
        probe_Tm_parameters["tmm_table"] = getattr(mt, probe_Tm_parameters["tmm_table"])
        probe_Tm_parameters["imm_table"] = getattr(mt, probe_Tm_parameters["imm_table"])
        probe_Tm_parameters["de_table"] = getattr(mt, probe_Tm_parameters["de_table"])

        # filter oligos by properties
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

        # filter oligos by specificity
        oligo_database, file_database = pipeline.filter_by_specificity(
            oligo_database=oligo_database,
            files_fasta_reference_database=config["probe_files_fasta_reference_database"],
            blastn_search_parameters=config["specificity_blastn_search_parameters"],
            blastn_hit_parameters=config["specificity_blastn_hit_parameters"],
            cross_hybridization_blastn_search_parameters=config[
                "cross_hybridization_blastn_search_parameters"
            ],
            cross_hybridization_blastn_hit_parameters=config["cross_hybridization_blastn_hit_parameters"],
        )

        # compute oligo attributes
        logging.info("Computing Oligo Attributes")
        oligo_database = pipeline.compute_oligo_attributes(
            oligo_database=oligo_database,
            Tm_parameters_probe=probe_Tm_parameters,
        )

        # create the probe sets
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

    #### Readout probes pipeline ####
    def run_readout_probes_pipeline(config):
        # read config file
        readout_probes_config_file = config["readout_probe_designer"]
        # make readout_probes_config_file an absolute path wrt the config file
        readout_probes_config_file = os.path.join(os.path.dirname(args["config"]), readout_probes_config_file)
        with open(readout_probes_config_file, "r") as handle:
            config_readout = yaml.safe_load(handle)

        pipeline_readout = MerfishReadoutProbeDesigner(
            write_intermediate_steps=config_readout["write_intermediate_steps"],
            dir_output=config_readout["dir_output"],
            n_jobs=config_readout["n_jobs"],
        )

        # get the readout probes
        # Alternatively, we could use the readout probes from the MERFISH paper
        codebook, readout_probe_table = pipeline_readout.get_readout_probes(
            oligo_length=config_readout["readout_length"],
            readout_sequence_probs=config_readout["readout_sequence_probs"],
            initial_num_sequences=config_readout["initial_num_sequences"],
            GC_content_min=config_readout["GC_content_min"],
            GC_content_max=config_readout["GC_content_max"],
            homopolymeric_base_n=config_readout["homopolymeric_base_n"],
            files_fasta_reference_database=config_readout["files_fasta_reference_database"],
            blastn_search_parameters=config_readout["blastn_search_parameters"],
            blastn_hit_parameters=config_readout["blastn_hit_parameters"],
            cross_hybridization_blastn_search_parameters=config_readout[
                "cross_hybridization_blastn_search_parameters"
            ],
            cross_hybridization_blastn_hit_parameters=config_readout[
                "cross_hybridization_blastn_hit_parameters"
            ],
            # region_ids=list(oligo_database.database.keys()),
            region_ids=["AARS1", "DECR2", "PRR35"],
            n_bits=config_readout["n_bits"],
            min_hamming_dist=config_readout["min_hamming_distance"],
            hamming_weight=config_readout["hamming_weight"],
            channels_ids=config_readout["channels_ids"],
        )

    #### Primer pipeline ####
    def run_primer_pipeline(config):
        # read config file
        primer_config_file = config["primer_designer"]
        # make primer_config_file an absolute path wrt the config file
        primer_config_file = os.path.join(os.path.dirname(args["config"]), primer_config_file)
        with open(primer_config_file, "r") as handle:
            config_primer = yaml.safe_load(handle)

        pipeline_primer = MerfishPrimerDesigner(
            write_intermediate_steps=config_primer["write_intermediate_steps"],
            dir_output=config_primer["dir_output"],
            n_jobs=config_primer["n_jobs"],
        )

        # get the reverse primer
        reverse_primer = config_primer["reverse_primer_sequence"]
        tm_reverse_primer = OligoAttributes()._calc_TmNN(
            sequence=reverse_primer,
            Tm_parameters=config_primer["tm_parameters"],
        )

        # get the forward primer
        forward_primer = pipeline_primer.get_forward_primer(
            oligo_length=config_primer["oligo_length"],
            oligo_base_probabilities=config_primer["oligo_base_probabilities"],
            initial_num_sequences=config_primer["initial_num_sequences"],
            gc_content_min=config_primer["gc_content_min"],
            gc_content_max=config_primer["gc_content_max"],
            gc_clamp_n_bases=config_primer["gc_clamp_n_bases"],
            gc_clamp_nGC=config_primer["gc_clamp_nGC"],
            tm_min=config_primer["tm_min"],
            tm_max=config_primer["tm_max"],
            tm_parameters=config_primer["tm_parameters"],
            homopolymeric_base_n=config_primer["homopolymeric_base_n"],
            files_fasta_reference_database_transcriptome=config_primer[
                "files_fasta_reference_database_transcriptome"
            ],
            blastn_search_parameters_transcriptome=config_primer["blastn_search_parameters_transcriptome"],
            blastn_hit_parameters_transcriptome=config_primer["blastn_hit_parameters_transcriptome"],
            files_fasta_reference_database_encoding_probe=config_primer[
                "files_fasta_reference_database_encoding_probe"
            ],
            blastn_search_parameters_encoding_probe=config_primer["blastn_search_parameters_encoding_probe"],
            blastn_hit_parameters_encoding_probe=config_primer["blastn_hit_parameters_encoding_probe"],
            tm_reverse_primer=tm_reverse_primer,
        )

    #### Assembling the final output ####
    # oligo_database = pipeline.design_final_probe_sequences(
    #     oligo_database=oligo_database,
    #     forward_primer_sequence=forward_primer,
    #     reverse_primer_sequence=reverse_primer,
    #     codebook=codebook,
    #     readout_probe_table=readout_probe_table,
    # )

    # pipeline.generate_output(oligo_database=oligo_database, top_n_sets=config["n_sets"])

    #### Final output ####
    run_readout_probes_pipeline(config)

    logging.info(f"Oligo sets were saved in {dir_oligosets}")
    logging.info("##### End of the pipeline. #####")


if __name__ == "__main__":
    main()
