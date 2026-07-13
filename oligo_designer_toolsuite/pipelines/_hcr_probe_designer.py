############################################
# imports
############################################

import os
import shutil
from pathlib import Path
from typing import Any

import pandas as pd
import yaml

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
    BasePropertyFilter,
    GCContentFilter,
    HardMaskedSequenceFilter,
    HomopolymericRunsFilter,
    MeltingTemperatureNNFilter,
    PropertyFilter,
    SecondaryStructureFilter,
    SoftMaskedSequenceFilter,
)
from oligo_designer_toolsuite.oligo_selection import IndependentSetsOligoSelection
from oligo_designer_toolsuite.oligo_specificity_filter import (
    AlignmentSpecificityFilter,
    BaseSpecificityFilter,
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
    preprocess_tm_parameters,
    validate_bit_mapping_table,
    validate_codebook,
)
from oligo_designer_toolsuite.sequence_generator import OligoSequenceGenerator
from oligo_designer_toolsuite.utils import configure_root_logger, logger

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

        ##### set class parameters #####
        self.write_intermediate_steps = write_intermediate_steps
        self.n_jobs = n_jobs

    def design_target_probes(
        self,
        target_probes_parameters: dict,
    ) -> OligoDatabase:
        """
        Design target probes for HCR experiments.

        Thin wrapper around :py:meth:`TargetProbeDesigner.generate_target_probes` — instantiates
        the inner designer and delegates the full multi-step workflow to it. Kept as a public API
        so callers can drive the target-probe stage without touching the inner class directly.

        :param target_probes_parameters: ``target_probes`` block from the pipeline config. Must
            contain ``oligo_generation``, ``property_filters``, ``specificity_filters``, and
            ``probe_set_selection`` sub-blocks. Tm parameters are expected to have been inlined
            into ``property_filters.Tm_filter`` by :func:`_preprocess_config`.
        :type target_probes_parameters: dict
        :return: An `OligoDatabase` containing the designed target probes organised into sets,
            with per-arm Tm properties annotated on the L/R half-probes.
        :rtype: OligoDatabase
        """
        target_probe_designer = TargetProbeDesigner(self.dir_output, self.n_jobs)
        target_probes_database = target_probe_designer.generate_target_probes(
            target_probes_parameters=target_probes_parameters,
            write_intermediate_steps=self.write_intermediate_steps,
        )

        # Compute per-arm Tm on both L and R halves. Tm parameters were inlined into ``Tm_filter``
        # by ``_preprocess_config``; the same values are reused here.
        tm_nn_property: BaseProperty = TmNNProperty(
            Tm_parameters=target_probes_parameters["property_filters"]["Tm_filter"]["Tm_parameters"],
            Tm_chem_correction_parameters=target_probes_parameters["property_filters"]["Tm_filter"][
                "Tm_chem_correction_parameters"
            ],
            Tm_salt_correction_parameters=target_probes_parameters["property_filters"]["Tm_filter"][
                "Tm_salt_correction_parameters"
            ],
        )
        calculator = PropertyCalculator(properties=[tm_nn_property])
        target_probes_database = calculator.apply(
            oligo_database=target_probes_database, sequence_type="oligo_L", n_jobs=self.n_jobs
        )
        target_probes_database = calculator.apply(
            oligo_database=target_probes_database, sequence_type="oligo_R", n_jobs=self.n_jobs
        )

        return target_probes_database

    def design_initiators(
        self,
        region_ids: list[str],
        initiator_probe_parameters: dict,
    ) -> tuple[pd.DataFrame, pd.DataFrame]:

        initiator_designer = InitiatorDesigner(
            dir_output=self.dir_output,
            n_jobs=self.n_jobs,
        )

        if initiator_probe_parameters["initiator_table"]["source"] == "load":
            initiator_table = initiator_designer.load_initiator_table(
                file_initiator_table=initiator_probe_parameters["initiator_table"]["file"]
            )
            initiator_table_source = initiator_probe_parameters["initiator_table"]["file"]
            logger.info(f"Loaded initiator table from file and retrieved {len(initiator_table)} initiators.")
        else:
            initiator_table = initiator_designer.generate_initiator_table()
            initiator_table_source = initiator_probe_parameters["initiator_table"]["source"]

        if initiator_probe_parameters["codebook"]["source"] == "load":
            codebook = initiator_designer.load_codebook(
                file_codebook=initiator_probe_parameters["codebook"]["file"]
            )
            codebook_source = initiator_probe_parameters["codebook"]["file"]
        else:
            codebook = initiator_designer.generate_codebook(region_ids=region_ids)
            codebook_source = initiator_probe_parameters["codebook"]["source"]

        initiator_designer.validate(
            codebook=codebook,
            initiator_table=initiator_table,
            region_ids=region_ids,
            codebook_source=codebook_source,
            initiator_table_source=initiator_table_source,
        )

        return codebook, initiator_table

    def assemble_hybridization_probes(
        self,
        oligo_database: OligoDatabase,
        hybridization_probe_parameters: dict,
    ) -> OligoDatabase:

        linker_sequence = hybridization_probe_parameters["linker_sequence"]
        codebook = hybridization_probe_parameters["codebook"]
        initiator_table = hybridization_probe_parameters["initiator_table"]

        region_ids = list(oligo_database.database.keys())

        oligo_database.set_database_sequence_types(
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

            probe_ids = list(oligo_database.database[region_id].keys())
            new_properties: dict[str, dict[str, str]] = {probe_id: {} for probe_id in probe_ids}

            for probe_id in probe_ids:
                new_properties[probe_id]["sequence_target"] = format_sequence(
                    database=oligo_database,
                    property="target",
                    region_id=region_id,
                    oligo_id=probe_id,
                )
                new_properties[probe_id]["sequence_oligo_L"] = format_sequence(
                    database=oligo_database,
                    property="oligo_L",
                    region_id=region_id,
                    oligo_id=probe_id,
                )
                new_properties[probe_id]["sequence_oligo_R"] = format_sequence(
                    database=oligo_database,
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
                        database=oligo_database,
                        property="oligo_L",
                        region_id=region_id,
                        oligo_id=probe_id,
                    )
                )

                new_properties[probe_id]["sequence_hybridization_probe_R"] = (
                    format_sequence(
                        database=oligo_database,
                        property="oligo_R",
                        region_id=region_id,
                        oligo_id=probe_id,
                    )
                    + linker_sequence
                    + sequence_initiator_R
                )

            oligo_database.update_oligo_properties(new_properties)

        return oligo_database

    def generate_output(
        self,
        oligo_database: OligoDatabase,
        codebook: pd.DataFrame,
        initiator_table: pd.DataFrame,
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
        codebook.to_csv(os.path.join(self.dir_output, "codebook.tsv"), sep="\t", index_label="gene_name")
        initiator_table.to_csv(os.path.join(self.dir_output, "initiators.tsv"), sep="\t")

        oligo_database.write_oligosets_to_yaml(
            properties=output_properties,
            ascending=True,
            filename="hcr_probes",
        )

        oligo_database.write_ready_to_order_yaml(
            properties=[
                "sequence_hybridization_probe_L",
                "sequence_hybridization_probe_R",
                "sequence_initiator_L",
                "sequence_initiator_R",
            ],
            ascending=True,
            filename="hcr_probes_order",
        )

        oligo_database.write_oligosets_to_table(
            properties=output_properties,
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

    def generate_target_probes(
        self,
        target_probes_parameters: dict,
        write_intermediate_steps: bool = False,
    ) -> OligoDatabase:
        """
        Generate HCR target probes by running the full multi-step target-probe design pipeline.

        HCR probes are split into two halves — an ``oligo_L`` (left) arm and an ``oligo_R`` (right)
        arm — separated by a short gap that carries the ligation site. Internally this method
        orchestrates the existing decorated steps :py:meth:`_create_oligo_database` →
        :py:meth:`_filter_by_property` → :py:meth:`_filter_by_specificity` →
        :py:meth:`_create_oligo_sets`, and then computes the per-arm melting temperature
        (``TmNNProperty``) on both ``oligo_L`` and ``oligo_R`` so downstream initiator assembly
        and reporting have Tm values available for each half-probe.

        :param target_probes_parameters: ``target_probes`` block. Must contain ``oligo_generation``,
            ``property_filters``, ``specificity_filters``, ``probe_set_selection`` sub-blocks.
            Tm parameters are expected to have been inlined into ``property_filters.Tm_filter``
            by :func:`_preprocess_config`.
        :type target_probes_parameters: dict
        :param write_intermediate_steps: If True, save the per-step target-probe databases for
            debugging.
        :type write_intermediate_steps: bool
        :return: An `OligoDatabase` containing the designed target probes organised into sets,
            with per-arm Tm properties annotated on ``oligo_L`` and ``oligo_R`` sequences.
        :rtype: OligoDatabase
        """
        oligo_generation_parameters = target_probes_parameters["oligo_generation"]
        property_filters_parameters = target_probes_parameters["property_filters"]
        specificity_filters_parameters = target_probes_parameters["specificity_filters"]
        probe_set_selection_parameters = target_probes_parameters["probe_set_selection"]

        oligo_database: OligoDatabase = self._create_oligo_database(
            region_ids=oligo_generation_parameters["region_ids"],
            oligo_length=oligo_generation_parameters["oligo_length"],
            L_probe_sequence_length=oligo_generation_parameters["L_probe_sequence_length"],
            gap_sequence_length=oligo_generation_parameters["gap_sequence_length"],
            R_probe_sequence_length=oligo_generation_parameters["R_probe_sequence_length"],
            files_fasta_oligo_database=oligo_generation_parameters["files_fasta_probe_database"],
            min_oligos_per_gene=probe_set_selection_parameters["independent_set_selection"]["set_size_min"],
        )

        if write_intermediate_steps:
            dir_database = oligo_database.save_database(name_database="1_db_target_probes_initial")
            logger.info(
                f"Saved target probe database for step 1 (Create Database) in directory {dir_database}"
            )

        oligo_database = self._filter_by_property(
            oligo_database=oligo_database,
            isoform_consensus_filter=property_filters_parameters["isoform_consensus_filter"],
            hard_masked_sequences_filter=property_filters_parameters["hard_masked_sequences_filter"],
            soft_masked_sequences_filter=property_filters_parameters["soft_masked_sequences_filter"],
            homopolymeric_runs_filter=property_filters_parameters["homopolymeric_runs_filter"],
            GC_content_filter=property_filters_parameters["GC_content_filter"],
            Tm_filter=property_filters_parameters["Tm_filter"],
            secondary_structure_filter=property_filters_parameters["secondary_structure_filter"],
        )

        if write_intermediate_steps:
            dir_database = oligo_database.save_database(name_database="2_db_target_probes_property_filter")
            logger.info(
                f"Saved target probe database for step 2 (Property Filters) in directory {dir_database}"
            )

        oligo_database = self._filter_by_specificity(
            oligo_database=oligo_database,
            specificity_blastn_filter=specificity_filters_parameters["specificity_blastn_filter"],
            cross_hybridization_blastn_filter=specificity_filters_parameters[
                "cross_hybridization_blastn_filter"
            ],
        )

        if write_intermediate_steps:
            dir_database = oligo_database.save_database(name_database="3_db_target_probes_specificity_filter")
            logger.info(
                f"Saved target probe database for step 3 (Specificity Filters) in directory {dir_database}"
            )

        oligo_database = self._create_oligo_sets(
            oligo_database=oligo_database,
            independent_set_selection=probe_set_selection_parameters["independent_set_selection"],
            isoform_consensus_score=probe_set_selection_parameters["isoform_consensus_score"],
        )

        if write_intermediate_steps:
            dir_database = oligo_database.save_database(name_database="4_db_target_probes_sets")
            logger.info(
                f"Saved target probe database for step 4 (Set Selection) in directory {dir_database}."
            )

        return oligo_database

    @pipeline_step_basic(step_name="Target Probe Generation - Create Database")
    def _create_oligo_database(
        self,
        region_ids: list[str] | None,
        oligo_length: int,
        L_probe_sequence_length: int,
        gap_sequence_length: int,
        R_probe_sequence_length: int,
        files_fasta_oligo_database: list[str],
        min_oligos_per_gene: int,
    ) -> OligoDatabase:

        ##### creating the oligo sequences #####
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

        ##### compute reverse complement (always on) #####
        reverse_complement_sequence_property: BaseProperty = ReverseComplementSequenceProperty(
            sequence_type_reverse_complement="oligo"
        )
        calculator = PropertyCalculator(properties=[reverse_complement_sequence_property])
        oligo_database = calculator.apply(
            oligo_database=oligo_database, sequence_type="target", n_jobs=self.n_jobs
        )

        ##### split the oligo sequence into L arm + spacer + R arm (always on) #####
        split_start_end = [
            (0, L_probe_sequence_length),
            (L_probe_sequence_length, L_probe_sequence_length + gap_sequence_length),
            (
                L_probe_sequence_length + gap_sequence_length,
                L_probe_sequence_length + gap_sequence_length + R_probe_sequence_length,
            ),
        ]
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

        oligo_database.remove_regions_with_insufficient_oligos(pipeline_step="Database Creation")
        check_content_oligo_database(oligo_database)

        return oligo_database

    @pipeline_step_basic(step_name="Target Probe Generation - Property Filters")
    def _filter_by_property(
        self,
        oligo_database: OligoDatabase,
        isoform_consensus_filter: dict,
        soft_masked_sequences_filter: dict,
        hard_masked_sequences_filter: dict,
        homopolymeric_runs_filter: dict,
        GC_content_filter: dict,
        Tm_filter: dict,
        secondary_structure_filter: dict,
    ) -> OligoDatabase:

        # Pre-filter by isoform consensus (cheap property lookup before sequence filters)
        if isoform_consensus_filter["enabled"]:
            isoform_consensus_property = IsoformConsensusProperty()
            calculator = PropertyCalculator(properties=[isoform_consensus_property])
            oligo_database = calculator.apply(
                oligo_database=oligo_database, sequence_type="target", n_jobs=self.n_jobs
            )
            oligo_database.filter_database_by_property_threshold(
                property_name="isoform_consensus",
                property_thr=isoform_consensus_filter["isoform_consensus"],
                remove_if_smaller_threshold=True,
            )

        # Build sequence-based property filter list, gating each filter on its own ``enabled`` flag.
        filters: list[BasePropertyFilter] = []
        if hard_masked_sequences_filter["enabled"]:
            hard_masked_sequences = HardMaskedSequenceFilter()
            filters.append(hard_masked_sequences)

        if soft_masked_sequences_filter["enabled"]:
            soft_masked_sequences = SoftMaskedSequenceFilter()
            filters.append(soft_masked_sequences)

        # Composition: homopolymeric runs, GC range, prohibited motifs
        if homopolymeric_runs_filter["enabled"]:
            homopolymeric_runs = HomopolymericRunsFilter(
                base_n=homopolymeric_runs_filter["homopolymeric_base_n"],
            )
            filters.append(homopolymeric_runs)

        if GC_content_filter["enabled"]:
            gc_content = GCContentFilter(
                GC_content_min=GC_content_filter["GC_content_min"],
                GC_content_max=GC_content_filter["GC_content_max"],
            )
            filters.append(gc_content)

        # Thermodynamics: self-complementarity (hairpins), Tm range, secondary structure (ΔG)
        if Tm_filter["enabled"]:
            melting_temperature = MeltingTemperatureNNFilter(
                Tm_min=Tm_filter["Tm_min"],
                Tm_max=Tm_filter["Tm_max"],
                Tm_parameters=Tm_filter["Tm_parameters"],
                Tm_chem_correction_parameters=Tm_filter["Tm_chem_correction_parameters"],
                Tm_salt_correction_parameters=Tm_filter["Tm_salt_correction_parameters"],
            )
            filters.append(melting_temperature)

        if secondary_structure_filter["enabled"]:
            secondary_structure = SecondaryStructureFilter(
                T=secondary_structure_filter["T"],
                thr_DG=secondary_structure_filter["thr_DG"],
            )
            filters.append(secondary_structure)

        # initialize the property filter class
        property_filter = PropertyFilter(filters=filters)

        # filter the database — same filters applied to both L and R arms
        oligo_database = property_filter.apply(
            oligo_database=oligo_database,
            sequence_type="oligo_L",
            n_jobs=self.n_jobs,
        )
        check_content_oligo_database(oligo_database)

        oligo_database = property_filter.apply(
            oligo_database=oligo_database,
            sequence_type="oligo_R",
            n_jobs=self.n_jobs,
        )
        check_content_oligo_database(oligo_database)

        return oligo_database

    @pipeline_step_basic(step_name="Target Probe Generation - Specificity Filters")
    def _filter_by_specificity(
        self,
        oligo_database: OligoDatabase,
        specificity_blastn_filter: dict,
        cross_hybridization_blastn_filter: dict,
    ) -> OligoDatabase:

        ##### exact match + specificity filter (exact_matches always on, specificity gated on enabled) #####
        exact_matches = ExactMatchFilter(policy=RemoveAllFilterPolicy(), filter_name="exact_match")
        filters: list[BaseSpecificityFilter] = [exact_matches]
        directories = []

        if specificity_blastn_filter["enabled"]:
            reference_database = ReferenceDatabase(
                database_name=f"{self.subdir_db_reference}_sequences", dir_output=self.dir_output
            )
            reference_database.load_database_from_file(
                files=specificity_blastn_filter["files_fasta_reference_database"],
                file_type="fasta",
                database_overwrite=True,
            )
            specificity: AlignmentSpecificityFilter
            if specificity_blastn_filter["junction_region_size"] > 0:
                oligo_ids = oligo_database.get_oligoid_list()
                junction_site = specificity_blastn_filter["junction_site"]
                oligo_database.update_oligo_properties(
                    new_oligo_property={oligo_id: {"junction_site": junction_site} for oligo_id in oligo_ids}
                )
                specificity = BlastNSeedregionSiteFilter(
                    seedregion_size=specificity_blastn_filter["junction_region_size"],
                    seedregion_site_name="junction_site",
                    search_parameters=specificity_blastn_filter["search_parameters"],
                    hit_parameters=specificity_blastn_filter["hit_parameters"],
                    filter_name="specificity_blastn_filter",
                    dir_output=self.dir_output,
                )
            else:
                specificity = BlastNFilter(
                    search_parameters=specificity_blastn_filter["search_parameters"],
                    hit_parameters=specificity_blastn_filter["hit_parameters"],
                    filter_name="specificity_blastn_filter",
                    dir_output=self.dir_output,
                )
            specificity.set_reference_database(reference_database=reference_database)
            filters.append(specificity)
            directories.append(specificity.dir_output)

        specificity_filter = SpecificityFilter(filters=filters)
        oligo_database = specificity_filter.apply(
            oligo_database=oligo_database,
            sequence_type="oligo",
            n_jobs=self.n_jobs,
        )
        check_content_oligo_database(oligo_database)

        ##### cross hybridization filter (gated on enabled) #####
        # Built and applied separately because it operates on the "oligo_L" and "oligo_R"
        # sequence types rather than on the joined "oligo" sequence.
        if cross_hybridization_blastn_filter["enabled"]:
            cross_hybridization_aligner_L = BlastNFilter(
                remove_hits=True,
                search_parameters=cross_hybridization_blastn_filter["search_parameters"],
                hit_parameters=cross_hybridization_blastn_filter["hit_parameters"],
                filter_name="cross_hybridization_blastn_filter",
                dir_output=self.dir_output,
            )
            cross_hybridization_L = CrossHybridizationFilter(
                policy=RemoveByLargerRegionFilterPolicy(),
                alignment_method=cross_hybridization_aligner_L,
                sequence_type_reference="oligo_L",
                filter_name="cross_hybridization_blastn_filter",
                dir_output=self.dir_output,
            )
            cross_hybridization_aligner_R = BlastNFilter(
                remove_hits=True,
                search_parameters=cross_hybridization_blastn_filter["search_parameters"],
                hit_parameters=cross_hybridization_blastn_filter["hit_parameters"],
                filter_name="cross_hybridization_blastn_filter",
                dir_output=self.dir_output,
            )
            cross_hybridization_R = CrossHybridizationFilter(
                policy=RemoveByLargerRegionFilterPolicy(),
                alignment_method=cross_hybridization_aligner_R,
                sequence_type_reference="oligo_R",
                filter_name="cross_hybridization_blastn_filter",
                dir_output=self.dir_output,
            )

            directories.extend(
                [
                    cross_hybridization_aligner_L.dir_output,
                    cross_hybridization_aligner_R.dir_output,
                    cross_hybridization_L.dir_output,
                    cross_hybridization_R.dir_output,
                ]
            )

            cross_hybridization_filter = SpecificityFilter(
                filters=[cross_hybridization_L, cross_hybridization_R]
            )
            oligo_database = cross_hybridization_filter.apply(
                oligo_database=oligo_database,
                sequence_type="oligo_L",
                n_jobs=self.n_jobs,
            )
            check_content_oligo_database(oligo_database)

            oligo_database = cross_hybridization_filter.apply(
                oligo_database=oligo_database,
                sequence_type="oligo_R",
                n_jobs=self.n_jobs,
            )
            check_content_oligo_database(oligo_database)

        ##### remove all directories of intermediate steps #####
        for directory in directories:
            if os.path.exists(directory):
                shutil.rmtree(directory)

        return oligo_database

    @pipeline_step_basic(step_name="Target Probe Generation - Set Selection")
    def _create_oligo_sets(
        self,
        oligo_database: OligoDatabase,
        independent_set_selection: dict,
        isoform_consensus_score: dict,
    ) -> OligoDatabase:

        # Define all scorers
        isoform_consensus_scorer = IsoformConsensusScorer(score_weight=isoform_consensus_score["weight"])
        oligos_scoring = OligoScoring(scorers=[isoform_consensus_scorer])
        set_scoring = AverageSetScoring(ascending=False)

        base_log_parameters({"Set Selection": "Independent Sets"})
        oligoset_generator = IndependentSetsOligoSelection(
            oligos_scoring=oligos_scoring,
            set_scoring=set_scoring,
            set_size_opt=independent_set_selection["set_size_opt"],
            set_size_min=independent_set_selection["set_size_min"],
            distance_between_oligos=independent_set_selection["distance_between_probes"],
            n_attempts_graph=independent_set_selection["n_attempts_graph"],
            n_attempts_clique_enum=independent_set_selection["n_attempts_clique_enum"],
            diversification_fraction=independent_set_selection["diversification_fraction"],
            jaccard_opt=independent_set_selection["jaccard_opt"],
            jaccard_step=independent_set_selection["jaccard_step"],
        )
        oligo_database = oligoset_generator.apply(
            oligo_database=oligo_database,
            sequence_type="oligo",
            n_sets=independent_set_selection["n_sets"],
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
        return pd.read_csv(file_codebook, sep=None, engine="python", index_col="gene_name")

    def generate_codebook(self, region_ids: list[str]) -> pd.DataFrame:
        """
        Generate a one-hot HCR codebook for the given regions.

        Placeholder for a future implementation. Once implemented, the output is
        expected to satisfy the same contract as a loaded codebook: ``gene_name``
        index, ``bit_*`` columns, exactly one bit set per row.
        """
        raise FeatureNotImplementedError(
            "Generation of codebook is not yet implemented. "
            "Please provide a codebook.file parameter and set codebook.source to 'load'."
        )

    def load_initiator_table(self, file_initiator_table: str) -> pd.DataFrame:
        initiator_table = pd.read_csv(file_initiator_table, sep=None, engine="python")
        if "bit" not in initiator_table.columns:
            raise FileFormatError(f"Initiator table '{file_initiator_table}' must contain a 'bit' column.")
        return initiator_table.set_index("bit")

    def generate_initiator_table(self) -> pd.DataFrame:
        """
        Generate an HCR initiator table with orthogonal L/R initiator sequences.

        Placeholder for a future implementation. Once implemented, the output is
        expected to satisfy the same contract as a loaded initiator table: ``bit``
        index, ``initiator_L_sequence`` / ``initiator_R_sequence`` columns
        containing DNA sequences.
        """
        raise FeatureNotImplementedError(
            "Generation of initiator table is not yet implemented. "
            "Please provide an initiator_table.file parameter and set initiator_table.source to 'load'."
        )

    def validate(
        self,
        codebook: pd.DataFrame,
        initiator_table: pd.DataFrame,
        region_ids: list[str],
        *,
        codebook_source: str,
        initiator_table_source: str,
    ) -> None:
        """
        Validate that a (codebook, initiator_table) pair forms a valid HCR initiator setup.

        Centralizes the HCR-specific validation contract (one-hot codebook indexed by
        ``gene_name``; bit-indexed initiator table with ``initiator_L_sequence`` /
        ``initiator_R_sequence`` DNA columns; codebook bits covered by the table) so
        that all paths producing these tables — loading from file today, generating
        programmatically in the future — share a single validation gate.

        :param codebook: Codebook DataFrame to validate.
        :type codebook: pd.DataFrame
        :param initiator_table: Initiator table DataFrame to validate.
        :type initiator_table: pd.DataFrame
        :param region_ids: Region IDs required to be present in the codebook index.
        :type region_ids: list[str]
        :param codebook_source: Source identifier (file path or marker) for the codebook.
        :type codebook_source: str
        :param initiator_table_source: Source identifier (file path or marker) for the initiator table.
        :type initiator_table_source: str
        :raises FileFormatError: If either input fails validation.
        """
        validate_codebook(
            codebook=codebook,
            region_ids=region_ids,
            source=codebook_source,
            expected_hamming_weight=1,
            index_name="gene_name",
        )
        validate_bit_mapping_table(
            table=initiator_table,
            codebook=codebook,
            source=initiator_table_source,
            required_columns=["initiator_L_sequence", "initiator_R_sequence"],
            sequence_columns=["initiator_L_sequence", "initiator_R_sequence"],
        )


############################################
# HCR Probe Designer Pipeline
############################################


def _preprocess_config(config: dict[str, Any]) -> dict[str, Any]:
    """
    Preprocess an HCR pipeline configuration dict in place.

    - Resolves the ``nn_table``/``tmm_table``/``imm_table``/``de_table`` strings in
      ``target_probes.global_parameters.Tm_parameters`` to their ``Bio.SeqUtils.MeltingTemp`` objects.
    - For every Tm chem/salt correction block: if ``enabled`` is ``False`` sets ``parameters`` to
      ``None`` so downstream filters receive a clean ``None``.
    - Inlines Tm parameters and chem/salt corrections into every block that consumes them
      (``Tm_filter``) so designer methods don't have to thread ``global_parameters`` through the
      call chain.
    - Computes the derived ``junction_site`` from the L arm + half the gap length and injects it
      into the ``specificity_blastn_filter`` block.
    - Expands ``target_probes.oligo_generation.file_region_ids`` to a sorted unique list under
      ``target_probes.oligo_generation.region_ids`` (or ``None`` if no file was provided).
    """

    # Preprocess Tm tables and set Tm_chem/salt_correction_parameters to None if the correction is disabled
    for section in ["target_probes"]:
        config[section]["global_parameters"]["Tm_parameters"] = preprocess_tm_parameters(
            config[section]["global_parameters"]["Tm_parameters"]
        )
        for correction in ["Tm_chem_correction_parameters", "Tm_salt_correction_parameters"]:
            correction_cfg = config[section]["global_parameters"][correction]
            if not correction_cfg["enabled"]:
                correction_cfg["parameters"] = None

    target_probe_Tm_parameters = config["target_probes"]["global_parameters"]["Tm_parameters"]
    target_probe_Tm_chem_correction_parameters = config["target_probes"]["global_parameters"][
        "Tm_chem_correction_parameters"
    ]["parameters"]
    target_probe_Tm_salt_correction_parameters = config["target_probes"]["global_parameters"][
        "Tm_salt_correction_parameters"
    ]["parameters"]

    # Inline Tm parameters into the Tm_filter block
    config["target_probes"]["property_filters"]["Tm_filter"]["Tm_parameters"] = target_probe_Tm_parameters
    config["target_probes"]["property_filters"]["Tm_filter"][
        "Tm_chem_correction_parameters"
    ] = target_probe_Tm_chem_correction_parameters
    config["target_probes"]["property_filters"]["Tm_filter"][
        "Tm_salt_correction_parameters"
    ] = target_probe_Tm_salt_correction_parameters

    # Compute the derived junction site (in the joined oligo coordinate space) and inject it
    # into the specificity BLASTN filter block for the seed-region variant.
    L_probe_sequence_length = config["target_probes"]["oligo_generation"]["L_probe_sequence_length"]
    gap_sequence_length = config["target_probes"]["oligo_generation"]["gap_sequence_length"]
    R_probe_sequence_length = config["target_probes"]["oligo_generation"]["R_probe_sequence_length"]
    config["target_probes"]["oligo_generation"]["oligo_length"] = (
        L_probe_sequence_length + gap_sequence_length + R_probe_sequence_length
    )
    config["target_probes"]["specificity_filters"]["specificity_blastn_filter"]["junction_site"] = (
        L_probe_sequence_length + gap_sequence_length // 2
    )

    ##### read the genes file #####
    file_region_ids = config["target_probes"]["oligo_generation"]["file_region_ids"]
    if file_region_ids is None:
        logger.warning(
            "No gene list file was provided! All genes from fasta file are used to generate the probes. "
            "This choice can use a lot of resources."
        )
        config["target_probes"]["oligo_generation"]["region_ids"] = None
    else:
        with open(file_region_ids) as f:
            config["target_probes"]["oligo_generation"]["region_ids"] = sorted({line.rstrip() for line in f})

    return config


def hcr_probe_designer(config: dict[str, Any]) -> None:
    """
    Execute the HCR probe design pipeline from a (raw) configuration dict.

    The dict is expected to follow the nested layout of ``data/configs/hcr_probe_designer.yaml``
    (``general``, ``target_probes.*``, ``hybridization_probe``). The caller is responsible for
    configuring the library logger before invoking this function (see :func:`main`).

    :param config: Pipeline configuration loaded via ``yaml.safe_load``.
    :type config: dict
    """

    ##### preprocess the config file #####
    config_dict = _preprocess_config(config)

    ##### initialize probe designer pipeline #####
    pipeline = HcrProbeDesigner(
        write_intermediate_steps=config_dict["general"]["write_intermediate_steps"],
        dir_output=config_dict["general"]["dir_output"],
        n_jobs=config_dict["general"]["n_jobs"],
    )

    ##### design target probes #####
    target_probe_database = pipeline.design_target_probes(
        target_probes_parameters=config_dict["target_probes"],
    )

    ##### design initiators (codebook + initiator table) #####
    codebook, initiator_table = pipeline.design_initiators(
        region_ids=list(target_probe_database.database.keys()),
        initiator_probe_parameters=config_dict["initiator_probes"],
    )

    ##### assemble hybridization probes #####
    config_dict["hybridization_probes"]["codebook"] = codebook
    config_dict["hybridization_probes"]["initiator_table"] = initiator_table
    hybridization_probe_database = pipeline.assemble_hybridization_probes(
        oligo_database=target_probe_database,
        hybridization_probe_parameters=config_dict["hybridization_probes"],
    )

    ##### write outputs #####
    pipeline.generate_output(
        oligo_database=hybridization_probe_database,
        codebook=codebook,
        initiator_table=initiator_table,
    )


def main() -> None:

    print("--------------START PIPELINE--------------")

    args = base_parser()

    ##### read the config file #####
    with open(args["config"], "r") as handle:
        config = yaml.safe_load(handle)

    # setup logger now that we know the output directory
    configure_root_logger(
        dir_output=config["general"]["dir_output"],
        pipeline_name="hcr_probe_designer",
    )

    hcr_probe_designer(config)

    print("--------------END PIPELINE--------------")


if __name__ == "__main__":
    main()
