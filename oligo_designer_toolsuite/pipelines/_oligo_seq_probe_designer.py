############################################
# imports
############################################

import logging
import os
import shutil
import warnings
from pathlib import Path
from typing import Any

import yaml

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
    setup_logging,
)
from oligo_designer_toolsuite.sequence_generator import OligoSequenceGenerator

############################################
# Oligo-Seq Probe Designer
############################################


class OligoSeqProbeDesigner:
    """
    A class for designing hybridization probes for Oligo-seq RNA detection assays.

    This class implements a complete pipeline for designing oligonucleotide probes compatible with the **Oligo-seq** method,
    a targeted RNA detection and sequencing approach that combines hybridization-based capture with next-generation sequencing (NGS)
    to profile gene expression at single-cell or subcellular resolution.

    **Oligo-seq Pipeline Overview:**
    1. **Target Probe Design**: Design gene-specific targeting sequences (26-30 nt) that bind to RNA transcripts
    2. **Output Generation**: Generate output files in multiple formats (TSV, YAML)

    Overview
    --------
    **Oligo-seq** (Oligonucleotide sequencing hybridization) is a novel RNA detection tool that merges *in situ hybridization* and *sequencing-based readout*.
    It allows multiplexed, quantitative RNA detection from extremely low input material (as few as ~50 cells) with high reproducibility and specificity. By
    focusing on *targeted exonic regions* and *exon–intron junctions*, Oligo-seq achieves robust capture of nascent and mature transcripts, enabling fine-grained
    resolution of gene expression states in different cell types (e.g., ESCs, XEN, and liver cells).

    Probe Structure
    ---------------
    **Oligo-seq Probes**
    - Single-stranded DNA oligonucleotides designed to hybridize directly to target RNA sequences.
    - Each probe consists of:
        - A **26-30nt targeting sequence** complementary to the RNA transcript.
    - A standard library includes ~92,000 probes covering ~1,800 genes, with an average of ~50 probes per gene.

    References
    ----------

    :ivar dir_output: Directory where probe design and library files will be written.
    :type dir_output: str
    :ivar write_intermediate_steps: Whether to save intermediate design and validation files (default: False).
    :type write_intermediate_steps: bool
    :ivar n_jobs: Number of parallel threads to use for probe design and computational validation.
    :type n_jobs: int
    """

    def __init__(self, write_intermediate_steps: bool, dir_output: str, n_jobs: int) -> None:
        """Constructor for the OligoSeqProbeDesigner class."""

        # create the output folder
        self.dir_output = os.path.abspath(dir_output)
        Path(self.dir_output).mkdir(parents=True, exist_ok=True)

        # setup logger
        setup_logging(
            dir_output=self.dir_output,
            pipeline_name="oligoseq_probe_designer",
            log_start_message=True,
        )

        ##### set class parameters #####
        self.write_intermediate_steps = write_intermediate_steps
        self.n_jobs = n_jobs

    def design_target_probes(
        self,
        oligo_generation_parameters: dict,
        property_filters_parameters: dict,
        specificity_filters_parameters: dict,
        probe_set_selection_parameters: dict,
        global_parameters: dict,
    ) -> OligoDatabase:
        """
        Design target probes for Oligo-seq experiments using a modular multi-step pipeline.

        This method orchestrates the complete target probe design workflow, consisting of:

        1. **Oligo generation**
        - Constructs an initial oligo database from input FASTA files
        - Uses a sliding-window approach with optional region splitting (e.g. exon–exon junction spanning)

        2. **Property-based filtering**
        - Filters candidate oligos based on intrinsic sequence properties, including:
            GC content, melting temperature (Tm), homopolymeric runs, secondary structure,
            self-complementarity, and sequence constraints

        3. **Specificity filtering**
        - Removes oligos with potential off-target binding using sequence alignment (e.g. BLASTN)
        - Filters cross-hybridizing probes and probes overlapping known variants
        - Applies read-length bias constraints where applicable

        4. **Probe set selection**
        - Organizes filtered oligos into optimal probe sets per target
        - Uses graph-based selection with scoring criteria such as:
            isoform consensus, exon targeting, GC content, melting temperature, and spacing constraints

        5. **Property annotation**
        - Computes final probe properties (e.g. GC content, Tm, transcript coverage)
        - Uses shared thermodynamic model parameters provided via `global_parameters`

        The resulting probes are gene-specific targeting sequences (typically 26-30 nt) that bind to
        RNA transcripts. These probes are used directly for hybridization-based capture in Oligo-seq
        experiments, which combine in situ hybridization with next-generation sequencing readout.

        :param oligo_generation_parameters:
            Parameters controlling oligo generation (input sequences, probe lengths, region splitting).
        :type oligo_generation_parameters: dict

        :param property_filters_parameters:
            Parameters for property-based filtering of oligos (e.g. GC content, Tm, sequence constraints).
        :type property_filters_parameters: dict

        :param specificity_filters_parameters:
            Parameters for specificity filtering (e.g. BLAST settings, variant filtering).
        :type specificity_filters_parameters: dict

        :param probe_set_selection_parameters:
            Parameters controlling probe set construction and scoring.
        :type probe_set_selection_parameters: dict

        :param global_parameters:
            Global/shared parameters used across the pipeline, including thermodynamic model settings
            (e.g. melting temperature model, salt and chemical corrections).
        :type global_parameters: dict

        :return:
            An `OligoDatabase` containing the designed target probes organized into sets,
            including computed properties and metadata for each probe.
        :rtype: OligoDatabase
        """

        target_probe_designer = TargetProbeDesigner(self.dir_output, self.n_jobs)

        oligo_database: OligoDatabase = target_probe_designer.create_oligo_database(
            region_ids=oligo_generation_parameters["region_ids"],
            oligo_length_min=oligo_generation_parameters["probe_length_min"],
            oligo_length_max=oligo_generation_parameters["probe_length_max"],
            split_region=oligo_generation_parameters["probe_split_region"],
            files_fasta_oligo_database=oligo_generation_parameters["files_fasta_probe_database"],
            min_oligos_per_gene=probe_set_selection_parameters["independent_set_selection"]["set_size_min"],
        )

        if self.write_intermediate_steps:
            dir_database = oligo_database.save_database(name_database="1_db_probes_initial")
            logging.info(f"Saved probe database for step 1 (Create Database) in directory {dir_database}")

        # Add highly abundant k-mers to prohibited sequences
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

        oligo_database = target_probe_designer.filter_by_property(
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
            global_parameters=global_parameters,
        )

        if self.write_intermediate_steps:
            dir_database = oligo_database.save_database(name_database="2_db_probes_property_filter")
            logging.info(f"Saved probe database for step 2 (Property Filters) in directory {dir_database}")

        oligo_database = target_probe_designer.filter_by_specificity(
            oligo_database=oligo_database,
            read_length_bias_filter=specificity_filters_parameters["read_length_bias_filter"],
            cross_hybridization_blastn_filter=specificity_filters_parameters[
                "cross_hybridization_blastn_filter"
            ],
            specificity_blastn_filter=specificity_filters_parameters["specificity_blastn_filter"],
            variant_filter=specificity_filters_parameters["variant_filter"],
        )

        if self.write_intermediate_steps:
            dir_database = oligo_database.save_database(name_database="3_db_probes_specificity_filter")
            logging.info(f"Saved probe database for step 3 (Specificity Filters) in directory {dir_database}")

        oligo_database = target_probe_designer.create_oligo_sets(
            oligo_database=oligo_database,
            independent_set_selection=probe_set_selection_parameters["independent_set_selection"],
            uniform_distance_score=probe_set_selection_parameters["uniform_distance_score"],
            isoform_consensus_score=probe_set_selection_parameters["isoform_consensus_score"],
            targeted_exons_score=probe_set_selection_parameters["targeted_exons_score"],
            GC_content_score=probe_set_selection_parameters["GC_content_score"],
            Tm_score=probe_set_selection_parameters["Tm_score"],
            global_parameters=global_parameters,
        )

        # Calculate oligo length, GC content, Tm, num targeted transcripts, isoform consensus, and length self complement
        length_property = LengthProperty()
        gc_content_property = GCContentProperty()
        TmNN_property = TmNNProperty(
            Tm_parameters=global_parameters["Tm_parameters"],
            Tm_chem_correction_parameters=global_parameters["Tm_chem_correction_parameters"],
            Tm_salt_correction_parameters=global_parameters["Tm_salt_correction_parameters"],
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
        oligo_database = calculator.apply(
            oligo_database=oligo_database, sequence_type="oligo", n_jobs=self.n_jobs
        )

        if self.write_intermediate_steps:
            dir_database = oligo_database.save_database(name_database="4_db_probes_probesets")
            logging.info(f"Saved probe database for step 4 (Specificity Filters) in directory {dir_database}")

        return oligo_database

    def generate_output(
        self,
        oligo_database: OligoDatabase,
        output_properties: list[str] | None = None,
    ) -> None:
        """
        Generate and export final probe design results in multiple output formats.

        This method serializes the designed probe sets from the provided `OligoDatabase`
        into user-facing output files suitable for downstream analysis and probe ordering.
        All files are written to the pipeline's configured output directory.

        **Generated Output Files**

        1. **oligo_seq_probes.yml**
        - Full probe set information in YAML format
        - Includes all selected properties for each probe

        2. **oligo_seq_probes.tsv / .xlsx**
        - Tabular representation of probe sets
        - Contains the same information as the YAML output
        - Excel output is organized with one sheet per region

        3. **oligo_seq_probes_order.yml**
        - Minimal output for probe ordering
        - Contains only essential sequence information (e.g. oligo sequences)

        :param oligo_database:
            Final `OligoDatabase` containing probe sequences, annotations, and computed properties.
        :type oligo_database: OligoDatabase

        :param output_properties:
            List of properties to include in the output files. If `None`, a default set of
            commonly used properties is written.

            Examples include:
            - metadata: 'source', 'species', 'gene_id', 'chromosome'
            - genomic location: 'start', 'end', 'strand', 'regiontype'
            - sequence: 'oligo', 'target'
            - computed properties: 'length_oligo', 'GC_content_oligo', 'TmNN_oligo',
            'isoform_consensus', 'length_selfcomplement_oligo'
        :type output_properties: list[str] | None

        :return:
            None. Output files are written to disk.
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
    A class for designing target probes for Oligo-seq experiments through a multi-step pipeline.

    This class provides methods for the complete target probe design process, which includes:
    1. Creating an initial oligo database from input FASTA files using a sliding window approach
       with region splitting for memory efficiency
    2. Filtering probes based on sequence properties (GC content, melting temperature, homopolymeric
       runs, self-complementarity, secondary structure)
    3. Filtering probes based on specificity to remove off-target binding, cross-hybridization,
       variants, and off-target binding using BLASTN, and read length bias filtering
    4. Organizing filtered probes into optimal sets based on weighted scoring criteria (targeted
       exons, isoform consensus, GC content, melting temperature) and distance constraints

    The resulting probes are gene-specific targeting sequences (typically 26-30 nt) that bind to
    RNA transcripts. These probes are used directly for hybridization-based capture in Oligo-seq
    experiments, which combine in situ hybridization with next-generation sequencing readout. The
    probes target exonic regions and exon-intron junctions to capture both nascent and mature
    transcripts.

    :param dir_output: Directory path where output files will be saved.
    :type dir_output: str
    :param n_jobs: Number of parallel jobs to use for processing. Set to 1 for serial processing or higher
        values for parallel processing. This affects the parallelization of filtering, property calculation,
        and set generation operations.
    :type n_jobs: int
    """

    def __init__(self, dir_output: str, n_jobs: int) -> None:
        """Constructor for the TargetProbeDesigner class."""

        ##### create the output folder #####
        self.dir_output = os.path.abspath(dir_output)
        self.subdir_db_probes = "db_probes"
        self.subdir_db_reference = "db_reference"

        self.n_jobs = n_jobs

    @pipeline_step_basic(step_name="Create Database")
    def create_oligo_database(
        self,
        region_ids: list | None,
        oligo_length_min: int,
        oligo_length_max: int,
        split_region: int,
        files_fasta_oligo_database: list[str],
        min_oligos_per_gene: int,
    ) -> OligoDatabase:
        """
        Generate an initial oligo database from input sequences using a sliding-window approach.

        This method represents the first step of the target probe design pipeline. It constructs
        a candidate oligo database by extracting sequences from input FASTA files and preparing
        them for downstream filtering and selection.

        **Workflow**

        1. **Oligo generation**
        - Generates candidate sequences using a sliding window across input regions
        - Supports variable oligo lengths within a specified interval
        - Splits large regions into smaller chunks to control memory usage

        2. **Database initialization**
        - Loads generated sequences into an `OligoDatabase`
        - Assigns sequence type `"target"` (original sequence)

        3. **Sequence augmentation**
        - Computes reverse complement sequences (`"oligo"`)
        - Prepares additional sequence types (e.g. `"oligo_short"`) for downstream steps

        4. **Initial filtering**
        - Removes regions with insufficient oligo counts based on `min_oligos_per_gene`
        - Validates database integrity

        The resulting database contains candidate oligos prior to property and specificity filtering.

        :param region_ids:
            Optional list of region or gene identifiers to restrict probe design.
            If `None`, all regions in the input FASTA files are processed.
        :type region_ids: list[str] | None

        :param oligo_length_min:
            Minimum oligo length (in nucleotides).
        :type oligo_length_min: int

        :param oligo_length_max:
            Maximum oligo length (in nucleotides).
        :type oligo_length_max: int

        :param split_region:
            Maximum size (in nucleotides) for splitting large input regions during sequence generation.
            This helps reduce memory usage for long sequences.
        :type split_region: int

        :param files_fasta_oligo_database:
            List of FASTA files containing input sequences (e.g. exons, junctions)
            from which candidate oligos are generated.
        :type files_fasta_oligo_database: list[str]

        :param min_oligos_per_gene:
            Minimum number of oligos required per region after initial generation.
            Regions below this threshold are removed from the database.
        :type min_oligos_per_gene: int

        :return:
            An `OligoDatabase` containing generated candidate oligos with initialized
            sequence types ("target", "oligo", "oligo_short"). The database is filtered
            to include only regions meeting the minimum oligo requirement.
        :rtype: OligoDatabase
        """

        # generate candidate oligo sequences (sliding window over FASTA regions)
        oligo_sequences = OligoSequenceGenerator(dir_output=self.dir_output)
        oligo_fasta_file = oligo_sequences.create_sequences_sliding_window(
            files_fasta_in=files_fasta_oligo_database,
            length_interval_sequences=(oligo_length_min, oligo_length_max),
            split_region=split_region,
            region_ids=region_ids,
            n_jobs=self.n_jobs,
        )

        # load sequences into oligo database
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
        # Set all sequence types that will be used in this pipeline
        oligo_database.set_database_sequence_types(["target", "oligo", "oligo_short"])
        # Calculate reverse complement using new PropertyCalculator pattern
        # compute reverse complement (oligo) and isoform consensus per entry
        rc_sequence_property = ReverseComplementSequenceProperty(sequence_type_reverse_complement="oligo")
        calculator = PropertyCalculator(properties=[rc_sequence_property])
        oligo_database = calculator.apply(
            oligo_database=oligo_database, sequence_type="target", n_jobs=self.n_jobs
        )

        # remove temporary sliding-window output directory
        dir = oligo_sequences.dir_output
        shutil.rmtree(dir) if os.path.exists(dir) else None

        # drop regions with too few oligos and validate database
        oligo_database.remove_regions_with_insufficient_oligos(pipeline_step="Database Creation")
        check_content_oligo_database(oligo_database)

        return oligo_database

    @pipeline_step_basic(step_name="Property Filters")
    def filter_by_property(
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
        global_parameters: dict,
    ) -> OligoDatabase:
        """
        Filter candidate oligos based on intrinsic sequence properties.

        This method applies a sequence of property-based filters to remove oligos with
        undesirable characteristics. Filters are applied in a cost-aware order
        (cheapest first) to reject invalid oligos as early as possible.

        **Workflow**

        1. **Pre-filtering (cheap property-based)**
        - Isoform consensus filtering (if enabled)
        - (Optional) targeted exon filtering (currently not implemented)

        2. **Sequence-based filtering**
        - Masking: removes oligos containing hard-masked (N) or soft-masked (lowercase) bases
        - Composition: filters based on homopolymeric runs, GC content, and prohibited motifs
        - Thermodynamics:
            - Self-complementarity (hairpin potential)
            - Melting temperature (Tm) range
            - Secondary structure stability (ΔG threshold)

        Oligos failing any filter are removed. Regions with insufficient remaining oligos
        are removed from the database.

        :param oligo_database:
            Input `OligoDatabase` containing candidate oligos with initialized sequence types.
        :type oligo_database: OligoDatabase

        :param isoform_consensus_filter:
            Dict with `enabled`, `isoform_consensus`.
        :type isoform_consensus_filter: dict

        :param targeted_exons_filter:
            Dict with `enabled`, `targeted_exons`.
        :type targeted_exons_filter: dict

        :param hard_masked_sequences_filter:
            Dict with `enabled`.
        :type hard_masked_sequences_filter: dict

        :param soft_masked_sequences_filter:
            Dict with `enabled`.
        :type soft_masked_sequences_filter: dict

        :param homopolymeric_runs_filter:
            Dict with `enabled`, `homopolymeric_base_n`.
        :type homopolymeric_runs_filter: dict

        :param GC_content_filter:
            Dict with `enabled`, `GC_content_min`, `GC_content_max`.
        :type GC_content_filter: dict

        :param prohibited_sequences_filter:
            Dict with `enabled`, `prohibited_sequences`, optionally `kmer_abundance_threshold`.
        :type prohibited_sequences_filter: dict

        :param self_complementarity_filter:
            Dict with `enabled`, `max_len_selfcomplement`.
        :type self_complementarity_filter: dict

        :param Tm_filter:
            Dict with `enabled`, `Tm_min`, `Tm_max`.
            Thermodynamic model parameters are provided via `global_parameters`.
        :type Tm_filter: dict

        :param secondary_structure_filter:
            Dict with `enabled`, `T`, `thr_DG`.
        :type secondary_structure_filter: dict

        :param global_parameters:
            Shared global parameters (e.g. thermodynamic model settings).
        :type global_parameters: dict

        :return:
            Filtered `OligoDatabase` containing only oligos passing all property filters.
        :rtype: OligoDatabase
        """
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
        # placeholder for yet missing targeted exon filter
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

        # Instantiate sequence-based property filters
        # Masking: drop oligos with ambiguous or low-complexity bases
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

        if prohibited_sequences_filter["enabled"]:
            prohibited_sequence_filter = ProhibitedSequenceFilter(
                prohibited_sequences=prohibited_sequences_filter["prohibited_sequences"],
            )
            filters.append(prohibited_sequence_filter)

        # Thermodynamics: self-complementarity (hairpins), Tm range, secondary structure (ΔG)
        if self_complementarity_filter["enabled"]:
            self_comp = SelfComplementFilter(
                max_len_selfcomplement=self_complementarity_filter["max_len_selfcomplement"],
            )
            filters.append(self_comp)

        if Tm_filter["enabled"]:
            melting_temperature = MeltingTemperatureNNFilter(
                Tm_min=Tm_filter["Tm_min"],
                Tm_max=Tm_filter["Tm_max"],
                Tm_parameters=global_parameters["Tm_parameters"],
                Tm_chem_correction_parameters=global_parameters["Tm_chem_correction_parameters"],
                Tm_salt_correction_parameters=global_parameters["Tm_salt_correction_parameters"],
            )
            filters.append(melting_temperature)

        if secondary_structure_filter["enabled"]:
            secondary_structure = SecondaryStructureFilter(
                T=secondary_structure_filter["T"],
                thr_DG=secondary_structure_filter["thr_DG"],
            )
            filters.append(secondary_structure)

        # Apply filters in order of cost (cheapest first) so failing oligos are rejected early.
        property_filter = PropertyFilter(filters=filters)
        oligo_database = property_filter.apply(
            oligo_database=oligo_database,
            sequence_type="oligo",
            n_jobs=self.n_jobs,
        )
        check_content_oligo_database(oligo_database)

        return oligo_database

    @pipeline_step_basic(step_name="Specificity Filters")
    def filter_by_specificity(
        self,
        oligo_database: OligoDatabase,
        read_length_bias_filter: dict,
        cross_hybridization_blastn_filter: dict,
        specificity_blastn_filter: dict,
        variant_filter: dict,
    ) -> OligoDatabase:
        """
        Filter candidate oligos based on sequence specificity.

        This method removes or flags oligos that are likely to bind non-specifically,
        cross-hybridize, or overlap with genomic variants.

        **Workflow**

        1. **Read-length bias filtering**
        - Removes oligos sharing identical prefixes of length N
        - Prevents ambiguity in sequencing readout

        2. **Exact match filtering**
        - Removes oligos with identical sequences across different regions

        3. **Cross-hybridization filtering**
        - Identifies oligos that can hybridize to each other
        - Removes oligos from the larger genomic region when conflicts occur

        4. **Reference specificity filtering**
        - Removes oligos with significant similarity to unintended targets
        - Uses alignment-based methods (e.g. BLASTN)

        5. **Variant filtering**
        - Flags or removes oligos overlapping known variants (VCF input)

        All filters operate on the `"oligo"` sequence type unless stated otherwise.
        Intermediate alignment files are removed after filtering.

        :param oligo_database:
            Input `OligoDatabase` containing oligos after property filtering.
        :type oligo_database: OligoDatabase

        :param read_length_bias_filter:
            Dict with `enabled`, `read_length_bias`.
        :type read_length_bias_filter: dict

        :param cross_hybridization_blastn_filter:
            Dict with `enabled`, `search_parameters`, `hit_parameters`.
        :type cross_hybridization_blastn_filter: dict

        :param specificity_blastn_filter:
            Dict with `enabled`, `search_parameters`, `hit_parameters`,
            `files_fasta_reference_database`.
        :type specificity_blastn_filter: dict

        :param variant_filter:
            Dict with `enabled`, `files_vcf_reference_database`, `action`
            ("flag" or "remove").
        :type variant_filter: dict

        :return:
            Filtered `OligoDatabase` containing oligos passing specificity constraints.
            Variant-overlapping oligos are flagged or removed depending on configuration.
        :rtype: OligoDatabase
        """

        # remove sequences that could cause read length biases because the first
        # <target_probe_read_length_bias> bases of both sequences match
        # Calculate shortened sequence using new PropertyCalculator pattern
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

        if variant_filter["enabled"]:
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

        # run all filters specified above
        specificity_filter = SpecificityFilter(filters=filters)
        oligo_database = specificity_filter.apply(
            oligo_database=oligo_database,
            sequence_type="oligo",
            n_jobs=self.n_jobs,
        )

        # remove all directories of intermediate steps
        for directory in directories:
            if os.path.exists(directory):
                shutil.rmtree(directory)

        check_content_oligo_database(oligo_database)

        return oligo_database

    @pipeline_step_basic(step_name="Oligo Selection")
    def create_oligo_sets(
        self,
        oligo_database: OligoDatabase,
        independent_set_selection: dict,
        uniform_distance_score: dict,
        isoform_consensus_score: dict,
        targeted_exons_score: dict,
        GC_content_score: dict,
        Tm_score: dict,
        global_parameters: dict,
    ) -> OligoDatabase:
        """
        Construct optimized oligo sets using scoring and graph-based selection.

        This method groups filtered oligos into sets per region, optimizing for multiple
        criteria such as sequence properties, spatial distribution, and diversity.

        **Workflow**

        1. **Oligo scoring**
        - Computes per-oligo scores based on:
            - Isoform consensus
            - Targeted exon overlap
            - GC content deviation from optimal
            - Melting temperature (Tm) deviation from optimal
            - Uniform spacing (distance-based scoring)

        2. **Set generation**
        - Builds a compatibility graph based on distance constraints
        - Identifies candidate sets using clique-based selection
        - Generates multiple sets per region with diversification

        3. **Set scoring and selection**
        - Scores sets using average oligo score
        - Selects top sets based on ranking (ascending = better)

        4. **Region filtering**
        - Removes regions that cannot produce valid sets
            (based on `set_size_min` constraint)

        The algorithm targets an optimal set size (`set_size_opt`) but may return
        smaller sets down to `set_size_min` if constraints are restrictive.

        :param oligo_database:
            Input `OligoDatabase` containing oligos after all filtering steps.
        :type oligo_database: OligoDatabase

        :param independent_set_selection:
            Dict controlling set generation (size, distance, attempts, diversity, Jaccard thresholds).
        :type independent_set_selection: dict

        :param uniform_distance_score:
            Dict with `weight`.
        :type uniform_distance_score: dict

        :param isoform_consensus_score:
            Dict with `weight`.
        :type isoform_consensus_score: dict

        :param targeted_exons_score:
            Dict with `weight`, `targeted_exons`.
        :type targeted_exons_score: dict

        :param GC_content_score:
            Dict with `GC_content_min`, `GC_content_opt`, `GC_content_max`, `weight`.
        :type GC_content_score: dict

        :param Tm_score:
            Dict with `Tm_min`, `Tm_opt`, `Tm_max`, `weight`.
        :type Tm_score: dict

        :param global_parameters:
            Shared global parameters (e.g. thermodynamic model settings).
        :type global_parameters: dict

        :return:
            Updated `OligoDatabase` containing selected oligo sets per region.
            Each region includes up to `n_sets` sets satisfying size and scoring constraints.
        :rtype: OligoDatabase
        """
        oligo_lengths = [
            len(sequence) for sequence in oligo_database.get_sequence_list(sequence_type="oligo")
        ]
        average_oligo_length = sum(oligo_lengths) / len(oligo_lengths)

        # Define all scorers
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
            Tm_parameters=global_parameters["Tm_parameters"],
            Tm_chem_correction_parameters=global_parameters["Tm_chem_correction_parameters"],
            Tm_salt_correction_parameters=global_parameters["Tm_salt_correction_parameters"],
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
            distance_between_oligos=independent_set_selection["distance_between_target_probes"],
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


def _load_config(config_file: str) -> dict[str, Any]:

    with open(config_file, "r") as handle:
        config: dict[str, Any] = yaml.safe_load(handle)

    config["target_probe"]["global"]["Tm_parameters"] = preprocess_tm_parameters(
        config["target_probe"]["global"]["Tm_parameters"]
    )

    ##### read the genes file #####
    file_region_ids = config["target_probe"]["oligo_generation"]["file_region_ids"]
    if file_region_ids is None:
        warnings.warn(
            "No gene list file was provided! All genes from fasta file are used to generate the probes. This chioce can use a lot of resources."
        )
        config["target_probe"]["oligo_generation"]["region_ids"] = None
    else:
        with open(file_region_ids) as f:
            config["target_probe"]["oligo_generation"]["region_ids"] = sorted({line.rstrip() for line in f})

    return config


def main() -> None:
    """
    Main entry point for running the Oligo-seq probe design pipeline.

    This function orchestrates the complete Oligo-seq probe design workflow:
    1. Parses command-line arguments using the base parser
    2. Reads the configuration YAML file containing all pipeline parameters
    3. Reads the gene IDs file (if provided) or uses all genes from FASTA files
    4. Preprocesses melting temperature parameters for target probes
    5. Preprocesses alignment method parameters for hybridization probability and cross-hybridization
       filtering (BLASTN or Bowtie)
    6. Initializes the OligoSeqProbeDesigner pipeline
    7. Designs target probes for specified genes
    8. Generates output files (YAML, TSV, Excel, order file)

    The function is typically called from the command line:
    ``oligo_seq_probe_designer --config <path_to_config.yaml>``

    Command-line arguments are parsed using `base_parser()`, which expects:
    - `config`: Path to the YAML configuration file containing all pipeline parameters
    """
    logging.info("--------------START PIPELINE--------------")

    args = base_parser()

    ##### read the config file #####
    config = _load_config(args["config"])

    ##### initialize probe designer pipeline #####
    pipeline = OligoSeqProbeDesigner(
        write_intermediate_steps=config["general"]["write_intermediate_steps"],
        dir_output=config["general"]["dir_output"],
        n_jobs=config["general"]["n_jobs"],
    )

    ##### design probes #####
    oligo_database = pipeline.design_target_probes(
        oligo_generation_parameters=config["target_probe"]["oligo_generation"],
        property_filters_parameters=config["target_probe"]["property_filters"],
        specificity_filters_parameters=config["target_probe"]["specificity_filters"],
        probe_set_selection_parameters=config["target_probe"]["probe_set_selection"],
        global_parameters=config["target_probe"]["global"],
    )

    pipeline.generate_output(oligo_database=oligo_database)

    logging.info("--------------END PIPELINE--------------")


if __name__ == "__main__":
    main()
