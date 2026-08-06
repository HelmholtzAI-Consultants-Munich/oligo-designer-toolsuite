"""
Oligo-seq probe designer pipeline.

Oligo-seq is a targeted RNA detection method that uses short DNA probes to bind
selected transcript regions. The bound probes are read out by sequencing, making
the method suitable for focused gene panels and low-input RNA measurements.

See :class:`OligoSeqProbeDesigner` for the full pipeline description and probe
structure. See :func:`oligo_seq_probe_designer` for the config-driven workflow.
"""

############################################
# imports
############################################

import os
import shutil
from pathlib import Path
from typing import Any

import yaml
from pydantic import ValidationError

from oligo_designer_toolsuite.config.pipelines.oligo_seq_probe_designer import OligoSeqProbeDesignerConfig
from oligo_designer_toolsuite.database import OligoDatabase, ReferenceDatabase
from oligo_designer_toolsuite.oligo_efficiency_filter import (
    AverageSetScoring,
    IsoformConsensusScorer,
    NormalizedDeviationFromOptimalGCContentScorer,
    NormalizedDeviationFromOptimalTmScorer,
    OligoScoring,
    OverlapTargetedExonsScorer,
    UniformDistanceScorer,
)
from oligo_designer_toolsuite.oligo_property_calculator import (
    GCContentProperty,
    IsoformConsensusProperty,
    LengthProperty,
    LengthSelfComplementProperty,
    NumTargetedTranscriptsProperty,
    PropertyCalculator,
    ReverseComplementSequenceProperty,
    ShortenedSequenceProperty,
    TmNNProperty,
)
from oligo_designer_toolsuite.oligo_property_filter import (
    BasePropertyFilter,
    GCContentFilter,
    HardMaskedSequenceFilter,
    HomopolymericRunsFilter,
    MeltingTemperatureNNFilter,
    ProhibitedSequenceFilter,
    PropertyFilter,
    SecondaryStructureFilter,
    SelfComplementFilter,
    SoftMaskedSequenceFilter,
)
from oligo_designer_toolsuite.oligo_selection import IndependentSetsOligoSelection
from oligo_designer_toolsuite.oligo_specificity_filter import (
    BaseSpecificityFilter,
    BlastNFilter,
    CrossHybridizationFilter,
    ExactMatchFilter,
    RemoveAllFilterPolicy,
    RemoveByLargerRegionFilterPolicy,
    SpecificityFilter,
    VariantsFilter,
)
from oligo_designer_toolsuite.pipelines._utils import (
    base_log_parameters,
    base_parser,
    check_content_oligo_database,
    get_highly_abundant_kmer_sequences,
    pipeline_step_basic,
    preprocess_tm_parameters,
)
from oligo_designer_toolsuite.sequence_generator import OligoSequenceGenerator
from oligo_designer_toolsuite.utils import configure_root_logger, logger

############################################
# Oligo-Seq Probe Designer
############################################


class OligoSeqProbeDesigner:
    """
    Design probe sets for Oligo-seq RNA detection assays.

    This class runs the design workflow for Oligo-seq probes. The output is a
    set of short single-stranded DNA oligos that bind to selected RNA regions
    and can be used for sequencing-based RNA detection.

    Overview
    --------
    Oligo-seq was developed as part of RNA-GAM, a method that adds RNA detection
    to Genome Architecture Mapping. GAM measures 3D genome architecture from
    ultrathin cryosections of single cells. RNA-GAM extends this idea by adding
    RNA measurements from the same type of material.

    In Oligo-seq, designed single-stranded DNA oligos hybridize to target RNA
    regions. The bound oligos are then read out by sequencing. This makes it
    possible to measure selected transcripts from very small samples, including
    ultrathin cryosections used in the GAM workflow.

    The assay is targeted rather than transcriptome-wide. The user chooses a
    gene panel, and the pipeline designs probes for those genes. Probes can be
    placed in exonic regions to detect mature transcripts and, when configured,
    across exon-intron junctions to capture information related to nascent RNA.

    In the RNA-GAM study, Oligo-seq was incorporated into the GAM pipeline and
    applied to mouse embryonic stem cells, extra-embryonic endoderm cells, and
    liver cells. The combined workflow was used to recover gene expression
    information that can help identify cell type and cell state.

    Probe Structure
    ---------------
    **Oligo-seq Probes**

    Oligo-seq probes are single-stranded DNA oligos that bind directly to RNA.
    Each probe contains one target-binding sequence. The designed probe itself is
    the sequence used for RNA targeting.

    A simplified probe layout is::

        5'--[target-binding sequence]--3'

    The target-binding sequence is usually short, for example around 26 to
    30 nucleotides, depending on the design settings. Several probes are usually
    selected for each gene. Using multiple probes makes the measurement less
    dependent on one local sequence region and helps produce a more stable signal
    for the target gene.

    Probe Library Preparation
    -------------------------
    The designed probes can be ordered as single-stranded DNA oligos or as a
    pooled oligo library, depending on panel size and laboratory workflow.

    In an Oligo-seq experiment, the probes are hybridized to RNA in the sample.
    After hybridization and washing, the oligo-derived signal is read by
    sequencing. In the RNA-GAM setting, this RNA readout is combined with the GAM
    workflow so that gene expression and genome architecture can be studied from
    related ultrathin cryosection material.

    The pipeline focuses on the design and selection of probe sequences. Wet-lab
    processing after synthesis depends on the specific Oligo-seq or RNA-GAM
    protocol used in the experiment.

    Pipeline Overview
    -----------------
    The pipeline performs the main steps needed to design an Oligo-seq probe
    set:

    1. **Target probe design**

       Design gene-specific target-binding sequences that hybridize to RNA
       transcripts.

    2. **Output generation**

       Write the selected probe sequences and their annotations to output files.
       The output includes detailed design files and a simplified file that can
       be used for oligo ordering.

    References
    ----------
    Sparks, T. (2024). Incorporating RNA Detection into Genome Architecture
    Mapping. Doctoral dissertation, Technische Universität Berlin.
    https://doi.org/10.14279/depositonce-19414

    :param dir_output: Directory where output files and intermediate results are
        saved. The directory is created if it does not exist.
    :type dir_output: str
    :param write_intermediate_steps: If ``True``, save intermediate probe
        databases after pipeline steps. This can help with checking a design run
        or finding where probes were removed.
    :type write_intermediate_steps: bool
    :param n_jobs: Number of worker processes used for steps that can run in
        parallel. Use ``1`` to run without parallel processing.
    :type n_jobs: int
    """

    def __init__(self, write_intermediate_steps: bool, dir_output: str, n_jobs: int) -> None:
        """Constructor for the OligoSeqProbeDesigner class."""
        self.dir_output = os.path.abspath(dir_output)
        Path(self.dir_output).mkdir(parents=True, exist_ok=True)

        self.write_intermediate_steps = write_intermediate_steps
        self.n_jobs = n_jobs

    def design_target_probes(
        self,
        target_probes_parameters: dict,
    ) -> OligoDatabase:
        """
        Design the RNA-binding probes for Oligo-seq.

        This step designs single-stranded DNA oligos that bind directly to selected
        RNA regions. Candidate probes are generated from the input sequences,
        filtered for sequence quality and specificity, and then selected into final
        probe sets for each target gene.

        After the probe sets are selected, this method adds useful reporting
        properties to each probe. These include probe length, GC content, melting
        temperature, targeted-transcript count, isoform consensus, and
        self-complementarity. The added values make the final output easier to
        inspect and compare.

        :param target_probes_parameters: Settings from the ``target_probes`` section
            of the pipeline config. This includes candidate generation, sequence
            filters, specificity filters, and probe set selection.
        :type target_probes_parameters: dict
        :return: Database containing the selected Oligo-seq probes for each target
            gene.
        :rtype: OligoDatabase
        """
        target_probe_designer = TargetProbeDesigner(self.dir_output, self.n_jobs)
        target_probes_database = target_probe_designer.generate_target_probes(
            target_probes_parameters=target_probes_parameters,
            write_intermediate_steps=self.write_intermediate_steps,
        )

        # Reporting properties for the output tables (not used as design filters here).
        length_property = LengthProperty()
        gc_content_property = GCContentProperty()
        TmNN_property = TmNNProperty(
            Tm_parameters=target_probes_parameters["property_filters"]["Tm_filter"]["Tm_parameters"],
            Tm_chem_correction_parameters=target_probes_parameters["property_filters"]["Tm_filter"][
                "Tm_chem_correction_parameters"
            ],
            Tm_salt_correction_parameters=target_probes_parameters["property_filters"]["Tm_filter"][
                "Tm_salt_correction_parameters"
            ],
        )
        num_targeted_transcripts_property = NumTargetedTranscriptsProperty()
        isoform_consensus_property = IsoformConsensusProperty()
        length_self_complement_property = LengthSelfComplementProperty()
        calculator = PropertyCalculator(
            properties=[
                length_property,
                gc_content_property,
                TmNN_property,
                num_targeted_transcripts_property,
                isoform_consensus_property,
                length_self_complement_property,
            ]
        )
        target_probes_database = calculator.apply(
            oligo_database=target_probes_database, sequence_type="oligo", n_jobs=self.n_jobs
        )

        return target_probes_database

    def generate_output(
        self,
        oligo_database: OligoDatabase,
        output_properties: list[str] | None = None,
    ) -> None:
        """
        Write the completed Oligo-seq probe design to files.

        This step saves the final probe database with probe sequences and selected
        annotations. It also writes an order-ready file that contains only the probe
        sequences needed for oligo synthesis.

        If no output properties are provided, a default set of annotation fields,
        sequence fields, and probe-quality measurements is written.

        :param oligo_database: Database returned by
            :py:meth:`design_target_probes`.
        :type oligo_database: OligoDatabase
        :param output_properties: Probe properties to include in the detailed output
            files. If ``None``, a default set of annotations, sequences, and probe
            quality values is used. The order file always contains only the
            ``oligo`` sequence.
        :type output_properties: list[str] | None
        :return: None
        :rtype: None
        """

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
                "oligo",
                "target",
                "length_oligo",
                "GC_content_oligo",
                "TmNN_oligo",
                "variants_filter",
                "isoform_consensus",
                "length_selfcomplement_oligo",
            ]

        oligo_database.write_oligosets_to_yaml(
            properties=output_properties,
            ascending=True,
            filename="oligo_seq_probes",
        )

        oligo_database.write_oligosets_to_table(
            properties=output_properties,
            ascending=True,
            filename="oligo_seq_probes",
        )

        oligo_database.write_ready_to_order_yaml(
            properties=["oligo"],
            ascending=True,
            filename="oligo_seq_probes_order",
        )


############################################
# Oligo-Seq Target Probe Designer
############################################


class TargetProbeDesigner:
    """
    Design the RNA-binding probes used in Oligo-seq.

    This class designs the target probes for an Oligo-seq assay. In this
    pipeline, the target probe is the final probe that will bind to RNA. There
    are no readout probes, barcode sequences, or primer-flanked template probes
    in this design step.

    The workflow has four main steps:

    1. **Candidate generation**

       Build candidate probes from transcript or target-region FASTA files. The
       candidates are generated by sliding a window across each input sequence.

    2. **Sequence filtering**

       Remove candidates with unsuitable sequence properties, such as poor GC
       content, long single-base runs, prohibited motifs, strong
       self-complementarity, masked sequence, unsuitable melting temperature, or
       stable secondary structure.

    3. **Specificity filtering**

       Remove candidates that may bind to unintended targets, cross-hybridize
       with other probes, overlap known variants, or be difficult to distinguish
       by short-read sequencing.

    4. **Probe set selection**

       Select final probe sets for each target gene. The selection can favor
       well-spaced probes, target requested exons, isoform coverage, suitable GC
       content, and suitable melting temperature.

    :param dir_output: Directory where output files and intermediate results are
        saved.
    :type dir_output: str
    :param n_jobs: Number of worker processes used for steps that can run in
        parallel. Use ``1`` to run without parallel processing.
    :type n_jobs: int
    """

    def __init__(self, dir_output: str, n_jobs: int) -> None:
        """Constructor for the TargetProbeDesigner class."""
        self.dir_output = os.path.abspath(dir_output)
        self.subdir_db_probes = "db_probes"
        self.subdir_db_reference = "db_reference"

        self.n_jobs = n_jobs

    def generate_target_probes(
        self,
        target_probes_parameters: dict,
        write_intermediate_steps: bool = False,
    ) -> OligoDatabase:
        """
        Run the full Oligo-seq target-probe design workflow.

        This method designs single-stranded DNA probes that bind directly to selected
        RNA regions. It starts from transcript or target-region sequences, creates
        candidate probes, filters them, checks their specificity, and selects final
        probe sets for each target gene.

        If the prohibited-sequence filter is enabled, the method can also add highly
        abundant k-mers from the input reference to the list of blocked motifs. This
        helps avoid short sequence motifs that occur too often and may increase
        background binding.

        :param target_probes_parameters: Settings from the ``target_probes`` section
            of the pipeline config. This includes candidate generation, sequence
            filters, specificity filters, and probe set selection.
        :type target_probes_parameters: dict
        :param write_intermediate_steps: If ``True``, save intermediate probe
            databases after each main step. This can help when checking where probes
            were removed.
        :type write_intermediate_steps: bool
        :return: Database containing selected Oligo-seq probe sets for each target
            gene.
        :rtype: OligoDatabase
        """
        oligo_generation_parameters = target_probes_parameters["oligo_generation"]
        property_filters_parameters = target_probes_parameters["property_filters"]
        specificity_filters_parameters = target_probes_parameters["specificity_filters"]
        probe_set_selection_parameters = target_probes_parameters["probe_set_selection"]

        oligo_database: OligoDatabase = self._create_oligo_database(
            region_ids=oligo_generation_parameters["region_ids"],
            oligo_length_min=oligo_generation_parameters["probe_length_min"],
            oligo_length_max=oligo_generation_parameters["probe_length_max"],
            split_region=oligo_generation_parameters["probe_split_region"],
            files_fasta_oligo_database=oligo_generation_parameters["files_fasta_probe_database"],
            min_oligos_per_gene=probe_set_selection_parameters["independent_set_selection"]["set_size_min"],
        )

        if write_intermediate_steps:
            dir_database = oligo_database.save_database(name_database="1_db_probes_initial")
            logger.info(f"Saved probe database for step 1 (Create Database) in directory {dir_database}")

        # Abundant k-mers bind everywhere in the reference; blacklist them with any
        # user-supplied motifs before property filtering.
        if property_filters_parameters["prohibited_sequences_filter"]["enabled"]:
            prohibited_sequences = []
            if property_filters_parameters["prohibited_sequences_filter"]["kmer_abundance_threshold"]:
                prohibited_sequences = get_highly_abundant_kmer_sequences(
                    files_fasta=oligo_generation_parameters["files_fasta_probe_database"],
                    kmer_abundance_threshold=property_filters_parameters["prohibited_sequences_filter"][
                        "kmer_abundance_threshold"
                    ],
                )
            property_filters_parameters["prohibited_sequences_filter"]["prohibited_sequences"] = list(
                set(
                    property_filters_parameters["prohibited_sequences_filter"]["prohibited_sequences"]
                    + prohibited_sequences
                )
            )

        oligo_database = self._filter_by_property(
            oligo_database=oligo_database,
            isoform_consensus_filter=property_filters_parameters["isoform_consensus_filter"],
            targeted_exons_filter=property_filters_parameters["targeted_exons_filter"],
            hard_masked_sequences_filter=property_filters_parameters["hard_masked_sequences_filter"],
            soft_masked_sequences_filter=property_filters_parameters["soft_masked_sequences_filter"],
            homopolymeric_runs_filter=property_filters_parameters["homopolymeric_runs_filter"],
            GC_content_filter=property_filters_parameters["GC_content_filter"],
            prohibited_sequences_filter=property_filters_parameters["prohibited_sequences_filter"],
            self_complementarity_filter=property_filters_parameters["self_complementarity_filter"],
            Tm_filter=property_filters_parameters["Tm_filter"],
            secondary_structure_filter=property_filters_parameters["secondary_structure_filter"],
        )

        if write_intermediate_steps:
            dir_database = oligo_database.save_database(name_database="2_db_probes_property_filter")
            logger.info(f"Saved probe database for step 2 (Property Filters) in directory {dir_database}")

        oligo_database = self._filter_by_specificity(
            oligo_database=oligo_database,
            read_length_bias_filter=specificity_filters_parameters["read_length_bias_filter"],
            cross_hybridization_blastn_filter=specificity_filters_parameters[
                "cross_hybridization_blastn_filter"
            ],
            specificity_blastn_filter=specificity_filters_parameters["specificity_blastn_filter"],
            variant_filter=specificity_filters_parameters["variant_filter"],
        )

        if write_intermediate_steps:
            dir_database = oligo_database.save_database(name_database="3_db_probes_specificity_filter")
            logger.info(f"Saved probe database for step 3 (Specificity Filters) in directory {dir_database}")

        oligo_database = self._create_oligo_sets(
            oligo_database=oligo_database,
            independent_set_selection=probe_set_selection_parameters["independent_set_selection"],
            uniform_distance_score=probe_set_selection_parameters["uniform_distance_score"],
            isoform_consensus_score=probe_set_selection_parameters["isoform_consensus_score"],
            targeted_exons_score=probe_set_selection_parameters["targeted_exons_score"],
            GC_content_score=probe_set_selection_parameters["GC_content_score"],
            Tm_score=probe_set_selection_parameters["Tm_score"],
        )

        if write_intermediate_steps:
            dir_database = oligo_database.save_database(name_database="4_db_probes_probesets")
            logger.info(f"Saved probe database for step 4 (Set Selection) in directory {dir_database}")

        return oligo_database

    @pipeline_step_basic(step_name="Create Database")
    def _create_oligo_database(
        self,
        region_ids: list | None,
        oligo_length_min: int,
        oligo_length_max: int,
        split_region: int,
        files_fasta_oligo_database: list[str],
        min_oligos_per_gene: int,
    ) -> OligoDatabase:
        """
        Create the first database of candidate Oligo-seq probes.

        Candidate probes are generated by sliding windows across the input sequences.
        All probe lengths between the minimum and maximum length are considered. For
        each candidate, the transcript-facing sequence is stored, and the reverse
        complement is stored as the DNA probe sequence that will bind to the RNA.

        Very long input regions can be split into smaller chunks before candidates
        are generated. This helps keep memory use manageable. Regions with too few
        candidate probes are removed at this stage.

        :param region_ids: Target regions to design probes for, usually gene names
            or gene IDs. If ``None``, all regions in the input FASTA files are used.
        :type region_ids: list[str] | None
        :param oligo_length_min: Minimum candidate probe length in bases.
        :type oligo_length_min: int
        :param oligo_length_max: Maximum candidate probe length in bases.
        :type oligo_length_max: int
        :param split_region: Maximum input-region length before a long region is
            split into smaller chunks.
        :type split_region: int
        :param files_fasta_oligo_database: FASTA files containing the transcript or
            target-region sequences used for probe design.
        :type files_fasta_oligo_database: list[str]
        :param min_oligos_per_gene: Minimum number of candidate probes a region must
            have to remain in the database.
        :type min_oligos_per_gene: int
        :return: Database containing candidate probes with the target sequence, the
            DNA probe sequence, and a reserved short-sequence field used by the
            read-length-bias filter.
        :rtype: OligoDatabase
        """
        oligo_sequences = OligoSequenceGenerator(dir_output=self.dir_output)
        oligo_fasta_file = oligo_sequences.create_sequences_sliding_window(
            files_fasta_in=files_fasta_oligo_database,
            length_interval_sequences=(oligo_length_min, oligo_length_max),
            split_region=split_region,
            region_ids=region_ids,
            n_jobs=self.n_jobs,
        )

        oligo_database = OligoDatabase(
            min_oligos_per_region=min_oligos_per_gene,
            write_regions_with_insufficient_oligos=True,
            max_entries_in_memory=self.n_jobs * 2 + 2,
            database_name=self.subdir_db_probes,
            dir_output=self.dir_output,
            n_jobs=1,
        )
        oligo_database.load_database_from_fasta(
            files_fasta=oligo_fasta_file,
            sequence_type="target",
            database_overwrite=True,
            region_ids=region_ids,
        )
        oligo_database.set_database_sequence_types(["target", "oligo", "oligo_short"])

        # Register oligo_short now so the read-length-bias filter can fill it later.
        # Probe strand is the reverse complement of the transcript ("target") window.
        rc_sequence_property = ReverseComplementSequenceProperty(sequence_type_reverse_complement="oligo")
        calculator = PropertyCalculator(properties=[rc_sequence_property])
        oligo_database = calculator.apply(
            oligo_database=oligo_database, sequence_type="target", n_jobs=self.n_jobs
        )

        dir = oligo_sequences.dir_output
        shutil.rmtree(dir) if os.path.exists(dir) else None

        oligo_database.remove_regions_with_insufficient_oligos(pipeline_step="Database Creation")
        check_content_oligo_database(oligo_database)

        return oligo_database

    @pipeline_step_basic(step_name="Property Filters")
    def _filter_by_property(
        self,
        oligo_database: OligoDatabase,
        isoform_consensus_filter: dict,
        targeted_exons_filter: dict,
        hard_masked_sequences_filter: dict,
        soft_masked_sequences_filter: dict,
        homopolymeric_runs_filter: dict,
        GC_content_filter: dict,
        prohibited_sequences_filter: dict,
        self_complementarity_filter: dict,
        Tm_filter: dict,
        secondary_structure_filter: dict,
    ) -> OligoDatabase:
        """
        Remove candidate probes with unsuitable sequence properties.

        This step checks whether each candidate probe is likely to work well in the
        assay. It can remove probes that overlap masked sequence, contain long
        single-base runs, contain blocked sequence motifs, have unsuitable GC
        content, have a melting temperature outside the chosen range, are strongly
        self-complementary, or are predicted to fold into stable structures.

        Annotation-based checks are applied to the target sequence. These include
        isoform consensus and, when implemented, targeted-exon selection. Sequence
        filters are applied to the DNA probe sequence that will hybridize to RNA.

        :param oligo_database: Candidate probe database returned by
            :py:meth:`_create_oligo_database`. This database is updated by the
            filtering step.
        :type oligo_database: OligoDatabase
        :param isoform_consensus_filter: Settings for keeping probes that target a
            sufficient fraction of annotated isoforms.
        :type isoform_consensus_filter: dict
        :param targeted_exons_filter: Settings for keeping probes that overlap
            selected exons. This option is reserved for targeted-exon filtering and
            is not active in the current implementation.
        :type targeted_exons_filter: dict
        :param hard_masked_sequences_filter: Settings for removing probes that
            overlap hard-masked bases, such as ``N`` bases.
        :type hard_masked_sequences_filter: dict
        :param soft_masked_sequences_filter: Settings for removing probes that
            overlap soft-masked sequence, often used for repetitive or low-complexity
            regions.
        :type soft_masked_sequences_filter: dict
        :param homopolymeric_runs_filter: Settings for removing probes with long
            runs of the same base.
        :type homopolymeric_runs_filter: dict
        :param GC_content_filter: Settings for the allowed GC-content range.
        :type GC_content_filter: dict
        :param prohibited_sequences_filter: Settings for removing probes that
            contain blocked sequence motifs. The motif list can include user-supplied
            sequences and highly abundant k-mers found in the input reference.
        :type prohibited_sequences_filter: dict
        :param self_complementarity_filter: Settings for removing probes with long
            internal self-complementary stretches.
        :type self_complementarity_filter: dict
        :param Tm_filter: Settings for the allowed melting-temperature range and the
            conditions used for the calculation.
        :type Tm_filter: dict
        :param secondary_structure_filter: Settings for removing probes predicted to
            form stable secondary structures.
        :type secondary_structure_filter: dict
        :return: Filtered database containing probes that passed all enabled
            sequence-property checks.
        :rtype: OligoDatabase
        """
        # Cheap property lookup first; drop weak isoform coverage before sequence work.
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

        # targeted_exons_filter is not implemented yet; see the method docstring.
        # if targeted_exons_filter["enabled"]:
        #     targeted_exons_property = TargetedExonsProperty()
        #     calculator = PropertyCalculator(properties=[targeted_exons_property])
        #     oligo_database = calculator.apply(
        #         oligo_database=oligo_database, sequence_type="target", n_jobs=self.n_jobs
        #     )
        #     oligo_database.filter_database_by_property_category(
        #         property_name="targeted_exons",
        #         property_category=targeted_exons_filter["targeted_exons"],
        #         remove_if_equals_category=False,
        #     )

        # Cheapest filters first so failing probes exit before thermodynamics.
        filters: list[BasePropertyFilter] = []
        if hard_masked_sequences_filter["enabled"]:
            hard_masked_sequences = HardMaskedSequenceFilter()
            filters.append(hard_masked_sequences)

        if soft_masked_sequences_filter["enabled"]:
            soft_masked_sequences = SoftMaskedSequenceFilter()
            filters.append(soft_masked_sequences)

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

        if prohibited_sequences_filter["enabled"]:
            prohibited_sequence_filter = ProhibitedSequenceFilter(
                prohibited_sequences=prohibited_sequences_filter["prohibited_sequences"],
            )
            filters.append(prohibited_sequence_filter)

        if self_complementarity_filter["enabled"]:
            self_comp = SelfComplementFilter(
                max_len_selfcomplement=self_complementarity_filter["max_len_selfcomplement"],
            )
            filters.append(self_comp)

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

        property_filter = PropertyFilter(filters=filters)
        oligo_database = property_filter.apply(
            oligo_database=oligo_database,
            sequence_type="oligo",
            n_jobs=self.n_jobs,
        )
        check_content_oligo_database(oligo_database)

        return oligo_database

    @pipeline_step_basic(step_name="Specificity Filters")
    def _filter_by_specificity(
        self,
        oligo_database: OligoDatabase,
        read_length_bias_filter: dict,
        cross_hybridization_blastn_filter: dict,
        specificity_blastn_filter: dict,
        variant_filter: dict,
    ) -> OligoDatabase:
        """
        Remove probes that may give ambiguous or off-target signal.

        This step checks whether candidate probes are specific enough for Oligo-seq.
        It removes exact duplicate probe sequences and, when enabled, probes that
        share the same short 5' sequence.

        The short 5' sequence check is used for sequencing-based readout. It helps
        make sure probes can still be told apart when reads do not cover the full
        probe sequence and the last bases are not sequenced to the end.

        When enabled, BLASTN-based filters remove probes that may bind to unintended
        reference sequences or cross-hybridize with other probes in the same panel.
        A variant filter can also check whether probes overlap known variants from a
        VCF file. Depending on the settings, variant-overlapping probes are either
        removed or marked in the output.

        :param oligo_database: Probe database returned by
            :py:meth:`_filter_by_property`. This database is updated by the
            specificity filters.
        :type oligo_database: OligoDatabase
        :param read_length_bias_filter: Settings for removing probes that share the
            same short 5' sequence. This helps distinguish probes even when
            sequencing reads do not cover the complete probe.
        :type read_length_bias_filter: dict
        :param cross_hybridization_blastn_filter: Settings for checking whether
            probes in the same panel may bind to each other or to unintended probe
            targets.
        :type cross_hybridization_blastn_filter: dict
        :param specificity_blastn_filter: Settings for checking probe specificity
            against reference sequences. This includes the reference FASTA files and
            BLASTN search settings.
        :type specificity_blastn_filter: dict
        :param variant_filter: Settings for checking whether probes overlap known
            variants from VCF files. Depending on the selected action, overlapping
            probes are removed or annotated.
        :type variant_filter: dict
        :return: Filtered database containing probes that passed the enabled
            specificity checks.
        :rtype: OligoDatabase
        """
        # Oligo-seq: short-read prefixes that match cannot be told apart, so drop
        # exact matches on the 5' window stored as oligo_short.
        if read_length_bias_filter["enabled"]:
            shortened_sequence_property = ShortenedSequenceProperty(
                sequence_length=read_length_bias_filter["read_length_bias"], reverse=False
            )
            calculator = PropertyCalculator(properties=[shortened_sequence_property])
            oligo_database = calculator.apply(
                oligo_database=oligo_database, sequence_type="oligo", n_jobs=self.n_jobs
            )
            exact_matches_short = ExactMatchFilter(
                policy=RemoveAllFilterPolicy(), filter_name="exact_match_read_length_bias"
            )
            specificity_filter = SpecificityFilter(filters=[exact_matches_short])
            oligo_database = specificity_filter.apply(
                oligo_database=oligo_database,
                sequence_type="oligo_short",
                n_jobs=self.n_jobs,
            )

        exact_matches = ExactMatchFilter(policy=RemoveAllFilterPolicy(), filter_name="exact_match")
        filters: list[BaseSpecificityFilter] = [exact_matches]
        directories = []

        if specificity_blastn_filter["enabled"]:
            reference_database_alignment = ReferenceDatabase(
                database_name=f"{self.subdir_db_reference}_sequences", dir_output=self.dir_output
            )
            reference_database_alignment.load_database_from_file(
                files=specificity_blastn_filter["files_fasta_reference_database"],
                file_type="fasta",
                database_overwrite=True,
            )
            specificity = BlastNFilter(
                remove_hits=True,
                search_parameters=specificity_blastn_filter["search_parameters"],
                hit_parameters=specificity_blastn_filter["hit_parameters"],
                filter_name="specificity_blastn_filter",
                dir_output=self.dir_output,
            )
            specificity.set_reference_database(reference_database=reference_database_alignment)
            filters.append(specificity)
            directories.append(specificity.dir_output)

        if cross_hybridization_blastn_filter["enabled"]:
            cross_hybridization_aligner = BlastNFilter(
                remove_hits=True,
                search_parameters=cross_hybridization_blastn_filter["search_parameters"],
                hit_parameters=cross_hybridization_blastn_filter["hit_parameters"],
                filter_name="cross_hybridization_blastn_filter",
                dir_output=self.dir_output,
            )
            cross_hybridization = CrossHybridizationFilter(
                policy=RemoveByLargerRegionFilterPolicy(),
                alignment_method=cross_hybridization_aligner,
                filter_name="cross_hybridization_blastn_filter",
                dir_output=self.dir_output,
            )
            filters.append(cross_hybridization)
            directories.append(cross_hybridization_aligner.dir_output)
            directories.append(cross_hybridization.dir_output)

        if variant_filter["enabled"]:
            # action "filter" drops hits; other actions keep probes but annotate them.
            remove_hits = variant_filter["action"] == "filter"
            reference_database_variants = ReferenceDatabase(
                database_name=f"{self.subdir_db_reference}_variants", dir_output=self.dir_output
            )
            reference_database_variants.load_database_from_file(
                files=variant_filter["files_vcf_reference_database"],
                file_type="vcf",
                database_overwrite=True,
            )
            variants = VariantsFilter(
                remove_hits=remove_hits, filter_name="variants_filter", dir_output=self.dir_output
            )
            variants.set_reference_database(reference_database=reference_database_variants)
            filters.append(variants)
            directories.append(variants.dir_output)

        specificity_filter = SpecificityFilter(filters=filters)
        oligo_database = specificity_filter.apply(
            oligo_database=oligo_database,
            sequence_type="oligo",
            n_jobs=self.n_jobs,
        )

        for directory in directories:
            if os.path.exists(directory):
                shutil.rmtree(directory)

        check_content_oligo_database(oligo_database)

        return oligo_database

    @pipeline_step_basic(step_name="Oligo Selection")
    def _create_oligo_sets(
        self,
        oligo_database: OligoDatabase,
        independent_set_selection: dict,
        uniform_distance_score: dict,
        isoform_consensus_score: dict,
        targeted_exons_score: dict,
        GC_content_score: dict,
        Tm_score: dict,
    ) -> OligoDatabase:
        """
        Select final Oligo-seq probe sets for each target gene.

        This step chooses groups of probes from the filtered candidates. The selected
        probes should cover the target region well and meet the requested number of
        probes per gene.

        Probe sets can be scored by several criteria. These include how evenly probes
        are spaced, whether they overlap requested exons, how well they represent
        annotated isoforms, how close their GC content is to the desired value, and
        how close their melting temperature is to the desired value.

        The method can keep more than one possible probe set per gene, which gives
        users alternatives when several good designs are available. Regions without
        enough suitable probes are removed.

        :param oligo_database: Filtered probe database returned by
            :py:meth:`_filter_by_specificity`. This database is updated with the
            selected probe sets.
        :type oligo_database: OligoDatabase
        :param independent_set_selection: Settings that control how many probe sets
            are selected, how many probes each set should contain, and how far apart
            selected probes should be placed.
        :type independent_set_selection: dict
        :param uniform_distance_score: Settings for scoring how evenly probes are
            spaced across the target region.
        :type uniform_distance_score: dict
        :param isoform_consensus_score: Settings for scoring probes by how well they
            represent the annotated isoforms of a target gene.
        :type isoform_consensus_score: dict
        :param targeted_exons_score: Settings for scoring probes that overlap
            selected exons.
        :type targeted_exons_score: dict
        :param GC_content_score: Settings for scoring probes by how close their GC
            content is to the desired value.
        :type GC_content_score: dict
        :param Tm_score: Settings for scoring probes by how close their melting
            temperature is to the desired value.
        :type Tm_score: dict
        :return: Database with selected probe sets attached to each remaining target
            gene.
        :rtype: OligoDatabase
        """
        oligo_lengths = [
            len(sequence) for sequence in oligo_database.get_sequence_list(sequence_type="oligo")
        ]
        # Uniform-distance scoring is in units of mean probe length so sets stay
        # roughly length-invariant when oligo lengths vary.
        average_oligo_length = sum(oligo_lengths) / len(oligo_lengths)

        uniform_distance_scorer = UniformDistanceScorer(
            average_oligo_length=average_oligo_length,
            score_weight=uniform_distance_score["weight"],
        )
        exon_scorer = OverlapTargetedExonsScorer(
            targeted_exons=targeted_exons_score["targeted_exons"],
            score_weight=targeted_exons_score["weight"],
            property_name="exon_number",
        )
        isoform_scorer = IsoformConsensusScorer(score_weight=isoform_consensus_score["weight"])
        Tm_scorer = NormalizedDeviationFromOptimalTmScorer(
            Tm_min=Tm_score["Tm_min"],
            Tm_opt=Tm_score["Tm_opt"],
            Tm_max=Tm_score["Tm_max"],
            Tm_parameters=Tm_score["Tm_parameters"],
            Tm_chem_correction_parameters=Tm_score["Tm_chem_correction_parameters"],
            Tm_salt_correction_parameters=Tm_score["Tm_salt_correction_parameters"],
            score_weight=Tm_score["weight"],
        )
        GC_scorer = NormalizedDeviationFromOptimalGCContentScorer(
            GC_content_min=GC_content_score["GC_content_min"],
            GC_content_opt=GC_content_score["GC_content_opt"],
            GC_content_max=GC_content_score["GC_content_max"],
            score_weight=GC_content_score["weight"],
        )

        oligos_scoring = OligoScoring(
            scorers=[exon_scorer, isoform_scorer, Tm_scorer, GC_scorer, uniform_distance_scorer]
        )
        set_scoring = AverageSetScoring(ascending=True)
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
# Oligo-seq Designer Pipeline
############################################


def _preprocess_config(config_validated: OligoSeqProbeDesignerConfig) -> dict[str, Any]:
    """
    Prepare the Oligo-seq config before the pipeline runs.

    This step converts the validated pydantic config to a plain dict and updates it
    so later design stages can read ready-to-use settings. It resolves
    melting-temperature tables, turns off unused temperature corrections, and copies
    the shared temperature settings into the filters and scoring steps that need
    them.

    It also expands an optional gene-list file into a concrete list of target
    regions. If no gene list is provided, all regions in the input FASTA files are
    used.

    :param config_validated: Pipeline configuration that has already passed pydantic
        validation.
    :type config_validated: OligoSeqProbeDesignerConfig
    :return: A config dict updated with the prepared settings.
    :rtype: dict
    """
    config = config_validated.model_dump()

    # Resolve Tm table names and blank disabled chem/salt corrections to None so
    # downstream filters treat None as "no correction" without checking the flag.
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

    # Inline shared Tm settings into the blocks that consume them.
    config["target_probes"]["property_filters"]["Tm_filter"]["Tm_parameters"] = target_probe_Tm_parameters
    config["target_probes"]["property_filters"]["Tm_filter"][
        "Tm_chem_correction_parameters"
    ] = target_probe_Tm_chem_correction_parameters
    config["target_probes"]["property_filters"]["Tm_filter"][
        "Tm_salt_correction_parameters"
    ] = target_probe_Tm_salt_correction_parameters

    config["target_probes"]["probe_set_selection"]["Tm_score"]["Tm_parameters"] = target_probe_Tm_parameters
    config["target_probes"]["probe_set_selection"]["Tm_score"][
        "Tm_chem_correction_parameters"
    ] = target_probe_Tm_chem_correction_parameters
    config["target_probes"]["probe_set_selection"]["Tm_score"][
        "Tm_salt_correction_parameters"
    ] = target_probe_Tm_salt_correction_parameters

    file_region_ids = config["target_probes"]["oligo_generation"]["file_region_ids"]
    if file_region_ids is None:
        logger.warning(
            "No gene list file was provided! All genes from fasta file are used to generate the probes. This choice can use a lot of resources."
        )
        config["target_probes"]["oligo_generation"]["region_ids"] = None
    else:
        with open(file_region_ids) as f:
            config["target_probes"]["oligo_generation"]["region_ids"] = sorted({line.rstrip() for line in f})

    return config


def oligo_seq_probe_designer(config: OligoSeqProbeDesignerConfig) -> None:
    """
    Run the Oligo-seq probe design pipeline from a config dict.

    This function prepares the config with :func:`_preprocess_config`, then runs
    :class:`OligoSeqProbeDesigner` end to end. It designs target probes and writes
    the final files under ``config.general.dir_output``. The caller should
    configure the library logger before calling this function (see :func:`main`).

    The config should follow ``data/configs/oligo_seq_probe_designer.yaml``.

    Top-level config sections:

    - ``general``: output directory, intermediate-step writing, and worker count.
    - ``target_probes``: candidate generation, sequence filters, specificity filters,
      and probe set selection.

    Files written under ``dir_output``:

    - ``oligo_seq_probes.yml``: full probe records.
    - ``oligo_seq_probes_order.yml``: sequences ready for synthesis.
    - ``oligo_seq_probes.tsv`` / ``oligo_seq_probes.xlsx``: probe sets as tables.

    Intermediate probe databases are also written when
    ``general.write_intermediate_steps`` is ``True``.

    See :class:`OligoSeqProbeDesigner` for the pipeline description and probe
    structure.

    :param config: Validated pipeline configuration. It is converted and prepared by
        :func:`_preprocess_config` before the pipeline runs.
    :type config: OligoSeqProbeDesignerConfig
    :return: None
    :rtype: None
    """
    config_dict = _preprocess_config(config)

    pipeline = OligoSeqProbeDesigner(
        write_intermediate_steps=config_dict["general"]["write_intermediate_steps"],
        dir_output=config_dict["general"]["dir_output"],
        n_jobs=config_dict["general"]["n_jobs"],
    )

    target_probe_database = pipeline.design_target_probes(
        target_probes_parameters=config_dict["target_probes"],
    )

    pipeline.generate_output(oligo_database=target_probe_database)


def main() -> None:
    """
    Run the Oligo-seq probe design pipeline from the command line.

    Parses the required ``-c``/``--config`` argument, loads the YAML configuration
    file, validates it with :class:`OligoSeqProbeDesignerConfig`, and configures
    the library logger to write under the configured output directory. It then
    calls :func:`oligo_seq_probe_designer`.

    :return: None
    :rtype: None
    :raises pydantic.ValidationError: If the configuration file fails schema
        validation.
    """
    print("--------------START PIPELINE--------------")

    args = base_parser(
        prog="Oligo-seq Probe Designer",
        usage="oligo_seq_probe_designer [options]",
        description=__doc__,
    )

    with open(args["config"], "r") as handle:
        config_raw = yaml.safe_load(handle)

    try:
        config_validated = OligoSeqProbeDesignerConfig.model_validate(config_raw)
    except ValidationError as e:
        print("Invalid configuration file:\n%s", e)
        raise

    # Configure logging only after dir_output is known so the log file lands there.
    configure_root_logger(
        dir_output=config_validated.general.dir_output,
        pipeline_name="oligoseq_probe_designer",
    )

    oligo_seq_probe_designer(config_validated)

    print("--------------END PIPELINE--------------")


if __name__ == "__main__":
    main()
