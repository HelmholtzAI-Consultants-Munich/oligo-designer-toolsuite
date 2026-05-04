############################################
# imports
############################################

import logging
import os
import shutil
import warnings
from pathlib import Path

import pandas as pd
import yaml
from Bio.SeqUtils import MeltingTemp as mt

from oligo_designer_toolsuite._exceptions import (
    FeatureNotImplementedError,
    FileFormatError,
)
from oligo_designer_toolsuite.database import OligoDatabase, ReferenceDatabase
from oligo_designer_toolsuite.oligo_efficiency_filter import (
    AverageSetScoring,
    IsoformConsensusScorer,
    OligoScoring,
)
from oligo_designer_toolsuite.oligo_property_calculator import (
    BaseProperty,
    IsoformConsensusProperty,
    PropertyCalculator,
    ReverseComplementSequenceProperty,
    SplitSequenceProperty,
    TmNNProperty,
)
from oligo_designer_toolsuite.oligo_property_filter import (
    GCContentFilter,
    HardMaskedSequenceFilter,
    HomopolymericRunsFilter,
    MeltingTemperatureNNFilter,
    PropertyFilter,
    SecondaryStructureFilter,
)
from oligo_designer_toolsuite.oligo_selection import (
    GraphBasedSelectionPolicy,
    GreedySelectionPolicy,
    OligoSelectionPolicy,
    OligosetGeneratorIndependentSet,
)
from oligo_designer_toolsuite.oligo_specificity_filter import (
    AlignmentSpecificityFilter,
    BlastNFilter,
    BlastNSeedregionSiteFilter,
    CrossHybridizationFilter,
    ExactMatchFilter,
    RemoveAllFilterPolicy,
    RemoveByLargerRegionFilterPolicy,
    SpecificityFilter,
)
from oligo_designer_toolsuite.pipelines._utils import (
    base_log_parameters,
    base_parser,
    check_content_oligo_database,
    format_sequence,
    pipeline_step_basic,
    setup_logging,
)
from oligo_designer_toolsuite.sequence_generator import OligoSequenceGenerator

############################################
# HCR Probe Designer
############################################


class HcrProbeDesigner:
    def __init__(
        self,
        write_intermediate_steps: bool,
        dir_output: str,
        n_jobs: int,
    ) -> None:
        """Constructor for the HcrProbeDesigner class."""

        # create the output folder
        self.dir_output = os.path.abspath(dir_output)
        Path(self.dir_output).mkdir(parents=True, exist_ok=True)

        # setup logger
        setup_logging(
            dir_output=self.dir_output,
            pipeline_name="hcr_probe_designer",
            log_start_message=True,
        )

        ##### set class parameters #####
        self.write_intermediate_steps = write_intermediate_steps
        self.n_jobs = n_jobs

    def design_target_probes(
        self,
        files_fasta_target_probe_database: list[str],
        files_fasta_reference_database_target_probe: list[str],
        region_ids: list[str] | None,
        target_probe_isoform_consensus: float,
        target_probe_L_probe_sequence_length: int,
        target_probe_gap_sequence_length: int,
        target_probe_R_probe_sequence_length: int,
        target_probe_GC_content_min: float,
        target_probe_GC_content_max: float,
        target_probe_Tm_min: float,
        target_probe_Tm_max: float,
        target_probe_Tm_parameters: dict,
        target_probe_Tm_chem_correction_parameters: dict | None,
        target_probe_Tm_salt_correction_parameters: dict | None,
        target_probe_homopolymeric_base_n: dict,
        target_probe_T_secondary_structure: float,
        target_probe_secondary_structures_threshold_deltaG: float,
        target_probe_specificity_blastn_search_parameters: dict,
        target_probe_specificity_blastn_hit_parameters: dict,
        target_probe_cross_hybridization_blastn_search_parameters: dict,
        target_probe_cross_hybridization_blastn_hit_parameters: dict,
        target_probe_junction_region_size: int,
        target_probe_isoform_weight: float,
        set_size_opt: int,
        set_size_min: int,
        distance_between_target_probes: int,
        n_sets: int,
        max_graph_size: int,
        n_attempts: int,
        heuristic: bool,
        heuristic_n_attempts: int,
    ) -> OligoDatabase:

        target_probe_designer = TargetProbeDesigner(self.dir_output, self.n_jobs)

        oligo_database: OligoDatabase = target_probe_designer.create_oligo_database(
            region_ids=region_ids,
            target_probe_L_probe_sequence_length=target_probe_L_probe_sequence_length,
            target_probe_gap_sequence_length=target_probe_gap_sequence_length,
            target_probe_R_probe_sequence_length=target_probe_R_probe_sequence_length,
            files_fasta_oligo_database=files_fasta_target_probe_database,
            min_oligos_per_gene=set_size_min,
            isoform_consensus=target_probe_isoform_consensus,
        )

        if self.write_intermediate_steps:
            dir_database = oligo_database.save_database(name_database="1_db_target_probes_initial")
            logging.info(
                f"Saved target probe database for step 1 (Create Database) in directory {dir_database}"
            )

        oligo_database = target_probe_designer.filter_by_property(
            oligo_database=oligo_database,
            GC_content_min=target_probe_GC_content_min,
            GC_content_max=target_probe_GC_content_max,
            Tm_min=target_probe_Tm_min,
            Tm_max=target_probe_Tm_max,
            homopolymeric_base_n=target_probe_homopolymeric_base_n,
            T_secondary_structure=target_probe_T_secondary_structure,
            secondary_structures_threshold_deltaG=target_probe_secondary_structures_threshold_deltaG,
            Tm_parameters=target_probe_Tm_parameters,
            Tm_chem_correction_parameters=target_probe_Tm_chem_correction_parameters,
            Tm_salt_correction_parameters=target_probe_Tm_salt_correction_parameters,
        )

        if self.write_intermediate_steps:
            dir_database = oligo_database.save_database(name_database="2_db_target_probes_property_filter")
            logging.info(
                f"Saved target probe database for step 2 (Property Filters) in directory {dir_database}"
            )

        oligo_database = target_probe_designer.filter_by_specificity(
            oligo_database=oligo_database,
            files_fasta_reference_database=files_fasta_reference_database_target_probe,
            junction_region_size=target_probe_junction_region_size,
            junction_site=target_probe_L_probe_sequence_length + target_probe_gap_sequence_length // 2,
            specificity_blastn_search_parameters=target_probe_specificity_blastn_search_parameters,
            specificity_blastn_hit_parameters=target_probe_specificity_blastn_hit_parameters,
            cross_hybridization_blastn_search_parameters=target_probe_cross_hybridization_blastn_search_parameters,
            cross_hybridization_blastn_hit_parameters=target_probe_cross_hybridization_blastn_hit_parameters,
        )

        if self.write_intermediate_steps:
            dir_database = oligo_database.save_database(name_database="3_db_target_probes_specificity_filter")
            logging.info(
                f"Saved target probe database for step 3 (Specificity Filters) in directory {dir_database}"
            )

        oligo_database = target_probe_designer.create_oligo_sets(
            oligo_database=oligo_database,
            isoform_weight=target_probe_isoform_weight,
            set_size_opt=set_size_opt,
            set_size_min=set_size_min,
            distance_between_oligos=distance_between_target_probes,
            n_sets=n_sets,
            max_graph_size=max_graph_size,
            n_attempts=n_attempts,
            heuristic=heuristic,
            heuristic_n_attempts=heuristic_n_attempts,
        )

        tm_nn_property: BaseProperty = TmNNProperty(
            Tm_parameters=target_probe_Tm_parameters,
            Tm_chem_correction_parameters=target_probe_Tm_chem_correction_parameters,
            Tm_salt_correction_parameters=target_probe_Tm_salt_correction_parameters,
        )

        calculator = PropertyCalculator(properties=[tm_nn_property])
        oligo_database = calculator.apply(
            oligo_database=oligo_database, sequence_type="oligo_L", n_jobs=self.n_jobs
        )
        oligo_database = calculator.apply(
            oligo_database=oligo_database, sequence_type="oligo_R", n_jobs=self.n_jobs
        )

        if self.write_intermediate_steps:
            dir_database = oligo_database.save_database(name_database="4_db_target_probes_sets")
            logging.info(
                f"Saved target probe database for step 4 (Specificity Filters) in directory {dir_database}."
            )

        return oligo_database

    def design_initiators(
        self,
        region_ids: list[str],
        file_initiator_table: str,
        file_codebook: str,
    ) -> tuple[pd.DataFrame, pd.DataFrame]:

        initiator_designer = InitiatorDesigner(
            dir_output=self.dir_output,
            n_jobs=self.n_jobs,
        )
        if file_initiator_table:
            initiator_table = initiator_designer.load_initiator_table(
                file_initiator_table=file_initiator_table
            )
            logging.info(f"Loaded initiator table from file and retrieved {len(initiator_table)} initiators.")
        else:
            raise FeatureNotImplementedError(
                "Generation of initiator table is not yet implemented. "
                "Please provide a file_initiator_table parameter."
            )

        if file_codebook:
            codebook = initiator_designer.load_codebook(file_codebook=file_codebook)
        else:
            raise FeatureNotImplementedError(
                "Generation of codebook is not yet implemented. " "Please provide a file_codebook parameter."
            )

        # Check if all region_ids are in the codebook
        missing_region_ids = set(region_ids) - set(codebook.index)
        if len(missing_region_ids) > 0:
            raise FileFormatError(
                f"Codebook is missing the following region IDs: {sorted(missing_region_ids)}. "
                f"Codebook contains {len(codebook)} regions: {sorted(codebook.index.tolist())}"
            )

        return codebook, initiator_table

    def assemble_hybridization_probes(
        self,
        target_probe_database: OligoDatabase,
        codebook: pd.DataFrame,
        initiator_table: pd.DataFrame,
        linker_sequence: str,
    ) -> OligoDatabase:

        region_ids = list(target_probe_database.database.keys())

        target_probe_database.set_database_sequence_types(
            [
                "sequence_target",
                "sequence_oligo_L",
                "sequence_oligo_R",
                "sequence_linker",
                "sequence_initiator_L",
                "sequence_initiator_R",
                "sequence_hybridization_probe_L",
                "sequence_hybridization_probe_R",
            ]
        )

        for region_id in region_ids:
            barcode = codebook.loc[region_id]
            bits = barcode[barcode == 1].index
            sequence_initiator_L = initiator_table.loc[bits, "initiator_L_sequence"].iloc[0]
            sequence_initiator_R = initiator_table.loc[bits, "initiator_R_sequence"].iloc[0]

            probe_ids = list(target_probe_database.database[region_id].keys())
            new_properties: dict[str, dict[str, str]] = {probe_id: {} for probe_id in probe_ids}

            for probe_id in probe_ids:
                new_properties[probe_id]["sequence_target"] = format_sequence(
                    database=target_probe_database,
                    property="target",
                    region_id=region_id,
                    oligo_id=probe_id,
                )
                new_properties[probe_id]["sequence_oligo_L"] = format_sequence(
                    database=target_probe_database,
                    property="oligo_L",
                    region_id=region_id,
                    oligo_id=probe_id,
                )
                new_properties[probe_id]["sequence_oligo_R"] = format_sequence(
                    database=target_probe_database,
                    property="oligo_R",
                    region_id=region_id,
                    oligo_id=probe_id,
                )
                new_properties[probe_id]["sequence_linker"] = linker_sequence
                new_properties[probe_id]["sequence_initiator_L"] = sequence_initiator_L
                new_properties[probe_id]["sequence_initiator_R"] = sequence_initiator_R

                new_properties[probe_id]["sequence_hybridization_probe_L"] = (
                    sequence_initiator_L
                    + linker_sequence
                    + format_sequence(
                        database=target_probe_database,
                        property="oligo_L",
                        region_id=region_id,
                        oligo_id=probe_id,
                    )
                )

                new_properties[probe_id]["sequence_hybridization_probe_R"] = (
                    format_sequence(
                        database=target_probe_database,
                        property="oligo_R",
                        region_id=region_id,
                        oligo_id=probe_id,
                    )
                    + linker_sequence
                    + sequence_initiator_R
                )

            target_probe_database.update_oligo_properties(new_properties)

        return target_probe_database

    def generate_output(
        self,
        probe_database: OligoDatabase,
        codebook: pd.DataFrame,
        initiator_table: pd.DataFrame,
        top_n_sets: int = 3,
        output_properties: list[str] | None = None,
    ) -> None:

        if output_properties is None:
            output_properties = [
                "source",
                "species",
                "annotation_release",
                "genome_assembly",
                "gene_id",
                "chromosome",
                "start",
                "end",
                "strand",
                "regiontype",
                "transcript_id",
                "exon_number",
                "sequence_target",
                "sequence_linker",
                "sequence_initiator_L",
                "sequence_initiator_R",
                "sequence_hybridization_probe_L",
                "sequence_hybridization_probe_R",
                "TmNN_oligo_L",
                "TmNN_oligo_R",
                "isoform_consensus",
            ]

        # write codebook and readout probe table
        codebook.to_csv(os.path.join(self.dir_output, "codebook.tsv"), sep="\t", index_label="region_id")
        initiator_table.to_csv(os.path.join(self.dir_output, "initiators.tsv"), sep="\t")

        probe_database.write_oligosets_to_yaml(
            properties=output_properties,
            top_n_sets=top_n_sets,
            ascending=True,
            filename="hcr_probes",
        )

        probe_database.write_ready_to_order_yaml(
            properties=[
                "sequence_hybridization_probe_L",
                "sequence_hybridization_probe_R",
                "sequence_initiator_L",
                "sequence_initiator_R",
            ],
            top_n_sets=top_n_sets,
            ascending=True,
            filename="hcr_probes_order",
        )

        probe_database.write_oligosets_to_table(
            properties=[
                "region_id",
                "oligoset_id",
                "oligo_id",
                "sequence_initiator_L",
                "sequence_linker",
                "sequence_oligo_L",
                "sequence_oligo_R",
                "sequence_linker",
                "sequence_initiator_R",
                "sequence_hybridization_probe_L",
                "sequence_hybridization_probe_R",
            ],
            top_n_sets=top_n_sets,
            ascending=True,
            filename="hcr_probes",
        )


############################################
# HCR Target Probe Designer
############################################


class TargetProbeDesigner:

    def __init__(self, dir_output: str, n_jobs: int) -> None:
        """Constructor for the TargetProbeDesigner class."""

        ##### create the output folder #####
        self.dir_output = os.path.abspath(dir_output)
        self.subdir_db_oligos = "db_target_probes"
        self.subdir_db_reference = "db_reference"

        self.n_jobs = n_jobs

    @pipeline_step_basic(step_name="Target Probe Generation - Create Database")
    def create_oligo_database(
        self,
        region_ids: list[str] | None,
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
            region_ids=region_ids,
            n_jobs=self.n_jobs,
        )

        ##### creating the oligo database #####
        oligo_database = OligoDatabase(
            min_oligos_per_region=min_oligos_per_gene,
            write_regions_with_insufficient_oligos=True,
            max_entries_in_memory=self.n_jobs * 2 + 2,
            database_name=self.subdir_db_oligos,
            dir_output=self.dir_output,
            n_jobs=1,
        )
        oligo_database.load_database_from_fasta(
            files_fasta=oligo_fasta_file,
            database_overwrite=True,
            sequence_type="target",
            region_ids=region_ids,
        )
        # Set all sequence types that will be used in this pipeline
        oligo_database.set_database_sequence_types(["target", "oligo", "oligo_L", "oligo_R"])

        ##### pre-filter oligo database for certain properties #####
        isoform_consensus_property: BaseProperty = IsoformConsensusProperty()
        reverse_complement_sequence_property: BaseProperty = ReverseComplementSequenceProperty(
            sequence_type_reverse_complement="oligo"
        )

        calculator = PropertyCalculator(
            properties=[isoform_consensus_property, reverse_complement_sequence_property]
        )
        oligo_database = calculator.apply(
            oligo_database=oligo_database, sequence_type="target", n_jobs=self.n_jobs
        )
        oligo_database.filter_database_by_property_threshold(
            property_name="isoform_consensus",
            property_thr=isoform_consensus,
            remove_if_smaller_threshold=True,
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

        # Calculate split sequence using new PropertyCalculator pattern
        # first right then left sequence because we are splitting the oligo not the target sequence
        split_sequence_property: BaseProperty = SplitSequenceProperty(
            split_start_end=split_start_end,
            split_names=["oligo_R", "spacer", "oligo_L"],
        )

        calculator = PropertyCalculator(properties=[split_sequence_property])
        oligo_database = calculator.apply(
            oligo_database=oligo_database, sequence_type="oligo", n_jobs=self.n_jobs
        )

        dir = oligo_sequences.dir_output
        shutil.rmtree(dir) if os.path.exists(dir) else None

        oligo_database.remove_regions_with_insufficient_oligos(pipeline_step="Pre-Filters")
        check_content_oligo_database(oligo_database)

        return oligo_database

    @pipeline_step_basic(step_name="Target Probe Generation - Property Filters")
    def filter_by_property(
        self,
        oligo_database: OligoDatabase,
        GC_content_min: float,
        GC_content_max: float,
        Tm_min: float,
        Tm_max: float,
        homopolymeric_base_n: dict,
        T_secondary_structure: float,
        secondary_structures_threshold_deltaG: float,
        Tm_parameters: dict,
        Tm_chem_correction_parameters: dict | None,
        Tm_salt_correction_parameters: dict | None,
    ) -> OligoDatabase:

        # define the filters
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
            sequence_type="oligo_L",
            n_jobs=self.n_jobs,
        )
        oligo_database.remove_regions_with_insufficient_oligos(pipeline_step="Property Filters")
        check_content_oligo_database(oligo_database)

        oligo_database = property_filter.apply(
            oligo_database=oligo_database,
            sequence_type="oligo_R",
            n_jobs=self.n_jobs,
        )
        oligo_database.remove_regions_with_insufficient_oligos(pipeline_step="Property Filters")
        check_content_oligo_database(oligo_database)

        return oligo_database

    @pipeline_step_basic(step_name="Target Probe Generation - Specificity Filters")
    def filter_by_specificity(
        self,
        oligo_database: OligoDatabase,
        files_fasta_reference_database: list[str],
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
        exact_matches = ExactMatchFilter(policy=RemoveAllFilterPolicy(), filter_name="oligo_exact_match")

        specificity: AlignmentSpecificityFilter
        if junction_region_size > 0:
            oligo_ids = oligo_database.get_oligoid_list()
            oligo_database.update_oligo_properties(
                new_oligo_property={oligo_id: {"junction_site": junction_site} for oligo_id in oligo_ids}
            )
            specificity = BlastNSeedregionSiteFilter(
                seedregion_size=junction_region_size,
                seedregion_site_name="junction_site",
                search_parameters=specificity_blastn_search_parameters,
                hit_parameters=specificity_blastn_hit_parameters,
                filter_name="oligo_blastn_specificity",
                dir_output=self.dir_output,
            )
        else:
            specificity = BlastNFilter(
                search_parameters=specificity_blastn_search_parameters,
                hit_parameters=specificity_blastn_hit_parameters,
                filter_name="oligo_blastn_specificity",
                dir_output=self.dir_output,
            )
        specificity.set_reference_database(reference_database=reference_database)

        ##### run specificity filters #####
        specificity_filter = SpecificityFilter(filters=[exact_matches, specificity])
        oligo_database = specificity_filter.apply(
            oligo_database=oligo_database,
            sequence_type="oligo",
            n_jobs=self.n_jobs,
        )
        oligo_database.remove_regions_with_insufficient_oligos(pipeline_step="Specificity Filters")
        check_content_oligo_database(oligo_database)

        ##### define cross hybridization filter #####
        cross_hybridization_aligner_oligo_pair_L = BlastNFilter(
            remove_hits=True,
            search_parameters=cross_hybridization_blastn_search_parameters,
            hit_parameters=cross_hybridization_blastn_hit_parameters,
            filter_name="oligo_L_R_blastn_crosshybridization",
            dir_output=self.dir_output,
        )
        cross_hybridization_oligo_pair_L = CrossHybridizationFilter(
            policy=RemoveByLargerRegionFilterPolicy(),
            alignment_method=cross_hybridization_aligner_oligo_pair_L,
            sequence_type_reference="oligo_L",
            filter_name="oligo_L_R_blastn_crosshybridization",
            dir_output=self.dir_output,
        )
        cross_hybridization_aligner_oligo_pair_R = BlastNFilter(
            remove_hits=True,
            search_parameters=cross_hybridization_blastn_search_parameters,
            hit_parameters=cross_hybridization_blastn_hit_parameters,
            filter_name="oligo_L_R_blastn_crosshybridization",
            dir_output=self.dir_output,
        )
        cross_hybridization_oligo_pair_R = CrossHybridizationFilter(
            policy=RemoveByLargerRegionFilterPolicy(),
            alignment_method=cross_hybridization_aligner_oligo_pair_R,
            sequence_type_reference="oligo_R",
            filter_name="oligo_L_R_blastn_crosshybridization",
            dir_output=self.dir_output,
        )

        ##### run cross hybridization filter #####
        specificity_filter = SpecificityFilter(
            filters=[cross_hybridization_oligo_pair_L, cross_hybridization_oligo_pair_R]
        )
        oligo_database = specificity_filter.apply(
            oligo_database=oligo_database,
            sequence_type="oligo_L",
            n_jobs=self.n_jobs,
        )
        oligo_database.remove_regions_with_insufficient_oligos(pipeline_step="Specificity Filters")
        check_content_oligo_database(oligo_database)

        oligo_database = specificity_filter.apply(
            oligo_database=oligo_database,
            sequence_type="oligo_R",
            n_jobs=self.n_jobs,
        )
        oligo_database.remove_regions_with_insufficient_oligos(pipeline_step="Specificity Filters")
        check_content_oligo_database(oligo_database)

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
        set_size_opt: int,
        set_size_min: int,
        distance_between_oligos: int,
        n_sets: int,
        max_graph_size: int,
        n_attempts: int,
        heuristic: bool,
        heuristic_n_attempts: int,
    ) -> OligoDatabase:

        # Define all scorers
        isoform_consensus_scorer = IsoformConsensusScorer(normalize=True, score_weight=isoform_weight)
        oligos_scoring = OligoScoring(scorers=[isoform_consensus_scorer])
        set_scoring = AverageSetScoring(ascending=False)

        # We change the processing dependent on the required number of probes in the probe sets
        # For small sets, we don't pre-filter and find the initial set by iterating
        # through all possible generated sets, which is faster than the max clique approximation.
        selection_policy: OligoSelectionPolicy
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

        oligo_database.remove_regions_with_insufficient_oligos(pipeline_step="Oligo Selection")
        check_content_oligo_database(oligo_database)

        return oligo_database


############################################
# HCR Initiator Designer
############################################


class InitiatorDesigner:

    def __init__(
        self,
        dir_output: str,
        n_jobs: int,
    ) -> None:
        """Constructor for the InitiatorDesigner class."""

        ##### create the output folder #####
        self.dir_output = os.path.abspath(dir_output)
        self.n_jobs = n_jobs

    def load_codebook(self, file_codebook: str) -> pd.DataFrame:

        codebook = pd.read_csv(file_codebook, sep=None, engine="python", index_col="region_id")

        # Check for at least one column
        if len(codebook.columns) == 0:
            raise FileFormatError(f"Codebook file '{file_codebook}' must contain at least one column.")

        # Check that all columns start with "bit_"
        non_bit_columns = [col for col in codebook.columns if not str(col).startswith("bit_")]
        if len(non_bit_columns) > 0:
            raise FileFormatError(
                f"Codebook file '{file_codebook}' must have all columns named with 'bit_*'. "
                f"Found columns that don't match: {non_bit_columns}"
            )

        # remove rows with any NaN
        codebook_clean = codebook[codebook.notna().all(axis=1)]
        # keep only rows with 0/1 values
        codebook_clean = codebook_clean[codebook_clean.isin([0, 1]).all(axis=1)]
        # keep only valid one-hot rows (exactly one "1")
        codebook_clean = codebook_clean[(codebook_clean == 1).sum(axis=1) == 1]
        if len(codebook_clean) == 0:
            raise FileFormatError(f"Codebook file '{file_codebook}' must contain at least one row with data.")

        return codebook

    def load_initiator_table(self, file_initiator_table: str) -> pd.DataFrame:

        required_cols = ["bit", "initiator_L_sequence", "initiator_R_sequence"]

        initiator_table = pd.read_csv(file_initiator_table, sep=None, engine="python")

        # Check if all required columns exist in readout_probe_table
        cols = set(initiator_table.columns)
        if not set(required_cols).issubset(cols):
            missing = set(required_cols) - cols
            raise FileFormatError(
                f"Initiator table is missing required columns: {missing}. "
                f"Required columns are: {required_cols}."
            )

        initiator_table = initiator_table[required_cols]
        initiator_table = initiator_table.set_index("bit")

        return initiator_table


############################################
# HCR Probe Designer Pipeline
############################################


def main() -> None:

    logging.info("--------------START PIPELINE--------------")

    args = base_parser()

    ##### read the config file #####
    with open(args["config"], "r") as handle:
        config = yaml.safe_load(handle)

    ##### read the genes file #####
    if config["file_regions"] is None:
        warnings.warn(
            "No gene list file was provided! All genes from fasta file are used to generate the probes. This choice can use a lot of resources."
        )
        gene_ids = None
    else:
        with open(config["file_regions"]) as handle:
            lines = handle.readlines()
            # ensure that the list contains unique gene ids
            gene_ids = list(set([line.rstrip() for line in lines]))

    # preprocess melting temperature params
    target_probe_Tm_parameters = config["target_probe_Tm_parameters"]
    target_probe_Tm_parameters["nn_table"] = getattr(mt, target_probe_Tm_parameters["nn_table"])
    target_probe_Tm_parameters["tmm_table"] = getattr(mt, target_probe_Tm_parameters["tmm_table"])
    target_probe_Tm_parameters["imm_table"] = getattr(mt, target_probe_Tm_parameters["imm_table"])
    target_probe_Tm_parameters["de_table"] = getattr(mt, target_probe_Tm_parameters["de_table"])

    ##### initialize probe designer pipeline #####
    pipeline = HcrProbeDesigner(
        write_intermediate_steps=config["write_intermediate_steps"],
        dir_output=config["dir_output"],
        n_jobs=config["n_jobs"],
    )

    ##### design probes #####
    target_probe_database = pipeline.design_target_probes(
        region_ids=gene_ids,
        files_fasta_target_probe_database=config["files_fasta_target_probe_database"],
        files_fasta_reference_database_target_probe=config["files_fasta_reference_database_target_probe"],
        target_probe_isoform_consensus=config["target_probe_isoform_consensus"],
        target_probe_L_probe_sequence_length=config["target_probe_L_probe_sequence_length"],
        target_probe_gap_sequence_length=config["target_probe_gap_sequence_length"],
        target_probe_R_probe_sequence_length=config["target_probe_R_probe_sequence_length"],
        target_probe_GC_content_min=config["target_probe_GC_content_min"],
        target_probe_GC_content_max=config["target_probe_GC_content_max"],
        target_probe_Tm_min=config["target_probe_Tm_min"],
        target_probe_Tm_max=config["target_probe_Tm_max"],
        target_probe_Tm_parameters=target_probe_Tm_parameters,
        target_probe_Tm_chem_correction_parameters=config["target_probe_Tm_chem_correction_parameters"],
        target_probe_Tm_salt_correction_parameters=config["target_probe_Tm_salt_correction_parameters"],
        target_probe_homopolymeric_base_n=config["target_probe_homopolymeric_base_n"],
        target_probe_T_secondary_structure=config["target_probe_T_secondary_structure"],
        target_probe_secondary_structures_threshold_deltaG=config[
            "target_probe_secondary_structures_threshold_deltaG"
        ],
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
        n_attempts=config["n_attempts"],
        heuristic=config["heuristic"],
        heuristic_n_attempts=config["heuristic_n_attempts"],
        target_probe_junction_region_size=config["target_probe_junction_region_size"],
        target_probe_isoform_weight=config["target_probe_isoform_weight"],
        set_size_opt=config["set_size_opt"],
        set_size_min=config["set_size_min"],
        distance_between_target_probes=config["distance_between_target_probes"],
        n_sets=config["n_sets"],
    )

    codebook, initiator_table = pipeline.design_initiators(
        region_ids=list(target_probe_database.database.keys()),
        file_initiator_table=config["file_initiator_table"],
        file_codebook=config["file_codebook"],
    )

    hybridization_probe_database = pipeline.assemble_hybridization_probes(
        target_probe_database=target_probe_database,
        codebook=codebook,
        initiator_table=initiator_table,
        linker_sequence=config["linker_sequence"],
    )

    pipeline.generate_output(
        probe_database=hybridization_probe_database,
        codebook=codebook,
        initiator_table=initiator_table,
        top_n_sets=config["top_n_sets"],
    )

    logging.info("--------------END PIPELINE--------------")


if __name__ == "__main__":
    main()
