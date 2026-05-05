############################################
# imports
############################################

import itertools
import logging
import os
import shutil
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from Bio.SeqUtils import Seq

from oligo_designer_toolsuite._exceptions import (
    ConfigurationError,
    FeatureNotImplementedError,
    FileFormatError,
)
from oligo_designer_toolsuite.database import OligoDatabase, ReferenceDatabase
from oligo_designer_toolsuite.oligo_efficiency_filter import (
    AverageSetScoring,
    DeviationFromOptimalTmScorer,
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
from oligo_designer_toolsuite.oligo_selection import IndependentSetsOligoSelection
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
    preprocess_tm_parameters,
    setup_logging,
)
from oligo_designer_toolsuite.sequence_generator import OligoSequenceGenerator

############################################
# CycleHCR Probe Designer
############################################


class CycleHCRProbeDesigner:
    """
    A class for designing hybridization probes for CycleHCR (Cyclic Hybridization Chain Reaction) experiments.

    This class provides a complete pipeline for designing CycleHCR hybridization probes, which are
    fluorescent probes used for multiplexed RNA detection.

    **CycleHCR Pipeline Overview:**
    1. **Target Probe Design**: Design split left/right (L/R) gene-specific targeting sequences that bind
       to adjacent regions on RNA transcripts, separated by a gap. Probes are designed with high melting
       temperature to remain bound during stripping cycles.
    2. **Readout Probe Assignment**: Load or generate readout probe sequences and create a codebook that
       assigns unique barcode pairs to each target region for multiplexed detection.
    3. **Hybridization Probe Assembly**: Combine target probe L/R halves with readout probe barcodes
       based on the codebook, creating complete hybridization probes with linker sequences.
    4. **Primer Design**: Load and validate PCR primers for amplifying DNA template probes. The forward
       primer is selected to match the reverse primer's melting temperature.
    5. **DNA Template Probe Assembly**: Assemble final DNA template probes by combining forward primers,
       target probe sequences, linker sequences, readout probe sequences, and reverse primers.
    6. **Output Generation**: Generate output files in multiple formats (TSV, YAML, Excel) containing
       probe sequences, properties, and codebook information.

    Overview
    --------
    cycleHCR (cyclic Hybridization Chain Reaction) is a multiplexed imaging method that combines
    split-initiator HCR and DNA barcoding to detect RNA and protein targets through multiple
    rounds of hybridization, imaging, and stripping. Each round reads out a subset of targets
    while leaving high–melting-temperature (Tm) hybridization probes bound, enabling hundreds of
    genes or proteins to be assayed sequentially within the same specimen.

    Probe Structure
    ---------------
    **Hybridization (primary) Probes**
    - Each RNA hybridization probe is divided into two ~45 nt halves: Left (L) and Right (R).
      They hybridize to adjacent regions on the target transcript, separated by a 2 nt gap.
    - The complete hybridization probe contains a 92-nt targeting sequence (divided into 45-nt segments
      for the left and right probe pairs, separated by a 2-nt gap), which directs binding to specific
      RNA targets.
    - L and R are designed for high Tm (>80–90 °C RNA:DNA), allowing them to remain bound after
      stringent stripping between imaging rounds.
    - Junctions are screened genome-wide to ensure unique binding; off-target activation requires
      both halves to bind adjacently, minimizing false positives.
    - Each hybridization probe also contains two 14-nt barcode sequences, TT-nucleotide spacers between
      readout and gene-specific regions.

    **Readout Probes**
    - Each target carries two short barcode sequences (L-barcode, R-barcode) within the hybridization probe.
    - In each imaging cycle, a pair of 14 nt readout oligos hybridize to these barcodes. Each readout
      oligo carries half of an 18 bp HCR initiator.
    - When both readouts bind adjacent barcodes, the initiator is reconstituted and triggers
      polymerization of fluorescent HCR hairpins (e.g., B2/B3/B4 hairpin sets).
    - The specific readout sequences contained by a hybridization probe are determined by the binary
      barcode assigned to that RNA target, enabling multiplexed detection of multiple RNA species.

    **DNA Template Probes**
    - The DNA template probe is assembled with a forward primer at the 5' end, followed by the target L/R
      sequences, a linker sequence, the reverse complement of the readout oligo sequences and a reverse
      primer at the 3' end, enabling PCR amplification.


    Probe Library Preparation
    -------------------------
    Hybridization probe libraries are DNA-synthesized, PCR-amplified with a T7 promoter,
    transcribed to RNA, reverse-transcribed, and then processed by USER II digestion
    and alkaline hydrolysis to yield single-stranded DNA probes.

    During imaging cycles, the HCR amplification produces localized fluorescent signal for each
    detected molecule. Readout oligos and hairpins are stripped with 80% formamide, while
    hybridization probes remain bound. New readout pairs are then applied in the next cycle.

    References
    ----------
    Gandin, V., Kim, J., Yang, L. Z., Lian, Y., Kawase, T., Hu, A., ... & Liu, Z. J. (2024).
    Deep-tissue spatial omics: imaging whole-embryo transcriptomics and subcellular structures
    at high spatial resolution. bioRxiv, 2024-05.

    :param dir_output: Directory path where output files will be saved. The directory will be created
        if it does not exist.
    :type dir_output: str
    :param write_intermediate_steps: Whether to save intermediate results during the probe design pipeline.
        If True, intermediate databases and results will be saved at each pipeline step, which is useful
        for debugging and analysis but increases disk usage.
    :type write_intermediate_steps: bool
    :param n_jobs: Number of parallel jobs to use for processing. Set to 1 for serial processing or higher
        values for parallel processing.
    :type n_jobs: int
    """

    def __init__(
        self,
        dir_output: str,
        write_intermediate_steps: bool,
        n_jobs: int,
    ) -> None:
        """Constructor for the CycleHCRProbeDesigner class."""

        # create the output folder
        self.dir_output = os.path.abspath(dir_output)
        Path(self.dir_output).mkdir(parents=True, exist_ok=True)

        # setup logger
        setup_logging(
            dir_output=self.dir_output,
            pipeline_name="cyclehcr_probe_designer",
            log_start_message=True,
        )
        self.write_intermediate_steps = write_intermediate_steps
        self.n_jobs = n_jobs

    def design_target_probes(
        self,
        region_ids: list[str] | None,
        files_fasta_target_probe_database: list[str],
        files_fasta_reference_database_target_probe: list[str],
        # Target Probe Design
        target_probe_isoform_consensus: float,
        target_probe_L_probe_sequence_length: int,
        target_probe_gap_sequence_length: int,
        target_probe_R_probe_sequence_length: int,
        # Property Filter Parameters
        target_probe_GC_content_min: float,
        target_probe_GC_content_max: float,
        target_probe_Tm_min: float,
        target_probe_Tm_max: float,
        target_probe_homopolymeric_base_n: dict,
        target_probe_T_secondary_structure: float,
        target_probe_secondary_structures_threshold_deltaG: float,
        # Melting Temperature Calculation Parameters
        target_probe_Tm_parameters: dict,
        target_probe_Tm_chem_correction_parameters: dict | None,
        target_probe_Tm_salt_correction_parameters: dict | None,
        # Specificity Filter Parameters
        target_probe_junction_region_size: int,
        target_probe_specificity_blastn_search_parameters: dict,
        target_probe_specificity_blastn_hit_parameters: dict,
        target_probe_cross_hybridization_blastn_search_parameters: dict,
        target_probe_cross_hybridization_blastn_hit_parameters: dict,
        # Probe Scoring and Set Selection Parameters
        target_probe_Tm_weight: float,
        target_probe_isoform_weight: float,
        set_size_opt: int,
        set_size_min: int,
        distance_between_target_probes: int,
        n_sets: int,
        n_attempts_graph: int,
        n_attempts_clique_enum: int,
        diversification_fraction: float,
        jaccard_opt: float,
        jaccard_step: float,
    ) -> OligoDatabase:
        """
        Design target probes for CycleHCR experiments through a multi-step pipeline.

        This method performs the complete target probe design process, which includes:
        1. Creating an initial oligo database from input FASTA files using a sliding window approach
        2. Filtering probes based on sequence properties (GC content, melting temperature, homopolymeric
           runs, secondary structure)
        3. Filtering probes based on specificity to remove off-target binding and cross-hybridization
           using BLASTN searches
        4. Organizing filtered probes into optimal sets based on scoring criteria and distance constraints

        The resulting probes are split into left (L) and right (R) halves that hybridize to adjacent
        regions on the target transcript, separated by a gap. These probes are designed with high
        melting temperatures to remain bound during stripping cycles.

        :param region_ids: List of region IDs (e.g., gene IDs) to target for probe design. If None,
            all regions present in the input FASTA files will be used.
        :type region_ids: list[str] | None
        :param files_fasta_target_probe_database: List of paths to FASTA files containing sequences
            from which target probes will be generated. These files should contain genomic regions
            of interest (e.g., exons, exon-exon junctions).
        :type files_fasta_target_probe_database: list[str]
        :param files_fasta_reference_database_target_probe: List of paths to FASTA files containing
            reference sequences used for specificity filtering. These files are used to identify
            off-target binding sites and potential cross-hybridization events (e.g., whole gene sequences).
        :type files_fasta_reference_database_target_probe: list[str]

        **Target Probe Design:**
        :param target_probe_isoform_consensus: Isoform consensus threshold for filtering target probes.
            Probes with isoform consensus values below this threshold will be filtered out. This parameter
            ensures that selected probes target sequences that are conserved across multiple transcript isoforms.
            Value should be between 0.0 and 1.0, where 1.0 indicates perfect consensus across all isoforms.
        :type target_probe_isoform_consensus: float
        :param target_probe_L_probe_sequence_length: Length of the left probe sequence in nucleotides.
            This is the 5' portion of the target probe that binds to the RNA.
        :type target_probe_L_probe_sequence_length: int
        :param target_probe_gap_sequence_length: Length of the gap sequence between left and right probes in nucleotides.
            This gap is not included in the probe sequences but represents the spacing between the two probe halves.
        :type target_probe_gap_sequence_length: int
        :param target_probe_R_probe_sequence_length: Length of the right probe sequence in nucleotides.
            This is the 3' portion of the target probe that binds to the RNA.
        :type target_probe_R_probe_sequence_length: int

        **Property Filter Parameters:**
        :param target_probe_GC_content_min: Minimum GC content (as a fraction between 0.0 and 1.0) for target probes.
            Probes with GC content below this value will be filtered out.
        :type target_probe_GC_content_min: float
        :param target_probe_GC_content_max: Maximum GC content (as a fraction between 0.0 and 1.0) for target probes.
            Probes with GC content above this value will be filtered out.
        :type target_probe_GC_content_max: float
        :param target_probe_Tm_min: Minimum melting temperature (Tm) in degrees Celsius for target probes.
            Probes with calculated Tm below this value will be filtered out.
        :type target_probe_Tm_min: float
        :param target_probe_Tm_max: Maximum melting temperature (Tm) in degrees Celsius for target probes.
            Probes with calculated Tm above this value will be filtered out. This value is also used as
            the optimal Tm target in probe scoring.
        :type target_probe_Tm_max: float
        :param target_probe_homopolymeric_base_n: Dictionary specifying the maximum allowed length of homopolymeric
            runs for each nucleotide base. Keys should be 'A', 'T', 'G', 'C' and values are the maximum run length.
            For example: {'A': 3, 'T': 3, 'G': 3, 'C': 3} allows up to 3 consecutive identical bases.
        :type target_probe_homopolymeric_base_n: dict[str, int]
        :param target_probe_T_secondary_structure: Temperature in degrees Celsius at which to evaluate secondary
            structure formation. Secondary structures that form at this temperature can interfere with probe binding.
        :type target_probe_T_secondary_structure: float
        :param target_probe_secondary_structures_threshold_deltaG: DeltaG threshold (in kcal/mol) for secondary
            structure stability. Probes with secondary structures having deltaG values more negative (more stable)
            than this threshold will be filtered out.
        :type target_probe_secondary_structures_threshold_deltaG: float

        **Specificity Filter Parameters:**
        :param target_probe_junction_region_size: Size of the junction region (in nucleotides) used for seed-based
            specificity filtering. If set to 0, full-length specificity filtering is used instead of seed-based filtering.
            When seed-based filtering is enabled, any probe with a BLASTN hit covering the junction region between
            the left and right probe halves will be removed, regardless of the alignment coverage percentage.
        :type target_probe_junction_region_size: int
        :param target_probe_specificity_blastn_search_parameters: Dictionary of BLASTN search parameters for specificity
            filtering. These parameters control how BLASTN searches are performed to identify off-target binding sites.
            Common parameters include: '-perc_identity', '-strand', '-word_size', '-dust', '-soft_masking',
            '-max_target_seqs', '-max_hsps'.
        :type target_probe_specificity_blastn_search_parameters: dict
        :param target_probe_specificity_blastn_hit_parameters: Dictionary of parameters for filtering BLASTN hits
            during specificity analysis. Common parameters include: 'min_alignment_length', 'coverage', etc.
            Probes with hits matching these criteria will be considered non-specific and filtered out.
        :type target_probe_specificity_blastn_hit_parameters: dict
        :param target_probe_cross_hybridization_blastn_search_parameters: Dictionary of BLASTN search parameters
            for cross-hybridization filtering. These parameters control how BLASTN searches are performed to identify
            potential cross-hybridization between left and right probe pairs within the same set.
        :type target_probe_cross_hybridization_blastn_search_parameters: dict
        :param target_probe_cross_hybridization_blastn_hit_parameters: Dictionary of parameters for filtering BLASTN
            hits during cross-hybridization analysis. Probes with cross-hybridization hits matching these criteria
            will be filtered out to prevent interference between probes in the same set.
        :type target_probe_cross_hybridization_blastn_hit_parameters: dict

        **Melting Temperature Calculation Parameters:**
        :param target_probe_Tm_parameters: Dictionary of parameters for calculating melting temperature (Tm) of target
            probes using the nearest-neighbor method. For using Bio.SeqUtils.MeltingTemp default parameters, set to ``{}``.
            Common parameters include: 'nn_table', 'tmm_table', 'imm_table', 'de_table', 'dnac1', 'dnac2', 'Na', 'K',
            'Tris', 'Mg', 'dNTPs', 'saltcorr', etc. For more information on parameters, see:
            https://biopython.org/docs/1.75/api/Bio.SeqUtils.MeltingTemp.html#Bio.SeqUtils.MeltingTemp.Tm_NN
        :type target_probe_Tm_parameters: dict
        :param target_probe_Tm_chem_correction_parameters: Dictionary of chemical correction parameters for Tm calculation.
            These parameters account for the effects of chemical additives (e.g., DMSO, formamide) on melting temperature.
            Set to ``None`` to disable chemical correction, or set to ``{}`` to use Bio.SeqUtils.MeltingTemp default parameters.
            For more information, see:
            https://biopython.org/docs/1.75/api/Bio.SeqUtils.MeltingTemp.html#Bio.SeqUtils.MeltingTemp.chem_correction
        :type target_probe_Tm_chem_correction_parameters: dict | None
        :param target_probe_Tm_salt_correction_parameters: Dictionary of salt correction parameters for Tm calculation.
            These parameters account for the effects of salt concentration on melting temperature. Set to ``None`` to disable
            salt correction, or set to ``{}`` to use Bio.SeqUtils.MeltingTemp default parameters. For more information, see:
            https://biopython.org/docs/1.75/api/Bio.SeqUtils.MeltingTemp.html#Bio.SeqUtils.MeltingTemp.salt_correction
        :type target_probe_Tm_salt_correction_parameters: dict | None

        **Probe Scoring and Set Selection Parameters:**
        :param target_probe_Tm_weight: Weight assigned to melting temperature (Tm) in the probe scoring function.
            Higher values prioritize probes with Tm closer to the optimal value (target_probe_Tm_max). This weight
            is used in combination with isoform_weight to calculate a composite score for each probe.
        :type target_probe_Tm_weight: float
        :param target_probe_isoform_weight: Weight assigned to isoform consensus in the probe scoring function.
            Higher values prioritize probes with higher isoform consensus values. This weight is used in combination
            with Tm_weight to calculate a composite score for each probe.
        :type target_probe_isoform_weight: float
        :param set_size_opt: Optimal size (number of probes) for each oligo set. The set selection algorithm will
            attempt to generate sets of this size, but may produce sets with fewer probes if constraints cannot be met.
        :type set_size_opt: int
        :param set_size_min: Minimum size (number of probes) required for each oligo set. Sets with fewer probes than
            this value will be rejected, and regions that cannot generate sets meeting this minimum will be removed.
        :type set_size_min: int
        :param distance_between_target_probes: Minimum genomic distance (in nucleotides) required between probes
            within the same set. This spacing constraint prevents probes from binding too close together, which could
            lead to reduced hybridization efficiency.
        :type distance_between_target_probes: int
        :param n_sets: Number of oligo sets to generate per region. Multiple sets allow for redundancy and selection
            of the best-performing set based on scoring criteria.
        :type n_sets: int
        :param n_attempts_graph: Number of randomized graph attempts. In each attempt, a fraction of nodes is randomly
            removed from the compatibility graph to create diversity; more attempts increase diversity at the cost of runtime.
        :type n_attempts_graph: int
        :param n_attempts_clique_enum: Maximum number of cliques enumerated per graph attempt. Limits how many cliques
            are explored before stopping enumeration for the current graph.
        :type n_attempts_clique_enum: int
        :param diversification_fraction: Fraction of oligos to remove from the graph per attempt to create diversity
            in the set selection.
        :type diversification_fraction: float
        :param jaccard_opt: Optimal maximum Jaccard overlap allowed between selected sets. Lower values enforce
            more diversity between sets.
        :type jaccard_opt: float
        :param jaccard_step: Step size used to relax the Jaccard constraint when not enough sets are found.
        :type jaccard_step: float

        :return: An `OligoDatabase` object containing the designed target probes organized into sets.
            The database includes probe sequences, properties, and set assignments for each target region.
        :rtype: OligoDatabase
        """
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
            Tm_max=target_probe_Tm_max,
            Tm_weight=target_probe_Tm_weight,
            Tm_parameters=target_probe_Tm_parameters,
            Tm_chem_correction_parameters=target_probe_Tm_chem_correction_parameters,
            Tm_salt_correction_parameters=target_probe_Tm_salt_correction_parameters,
            set_size_opt=set_size_opt,
            set_size_min=set_size_min,
            distance_between_oligos=distance_between_target_probes,
            n_sets=n_sets,
            n_attempts_graph=n_attempts_graph,
            n_attempts_clique_enum=n_attempts_clique_enum,
            diversification_fraction=diversification_fraction,
            jaccard_opt=jaccard_opt,
            jaccard_step=jaccard_step,
        )

        # Caculate all required properties for output
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

    def design_readout_probes(
        self,
        region_ids: list[str],
        file_readout_probe_table: str | None,
        file_codebook: str | None,
    ) -> tuple[pd.DataFrame, pd.DataFrame]:
        """
        Load or generate readout probes and codebook for CycleHCR experiments.

        This method handles the readout probe assignment process, which involves:
        1. Loading a readout probe table from a file (containing readout probe sequences and their
           associated bit identifiers, channels, and L/R designations)
        2. Loading an existing codebook from a file, or generating a new codebook that assigns
           unique binary barcodes to each target region

        The codebook is a binary matrix where each row corresponds to a region and each column
        represents a bit in the barcode. A value of 1 indicates that the corresponding readout
        probe should be used for that region. Each region is assigned a unique barcode consisting
        of two active bits (one for the left readout probe and one for the right readout probe).

        :param region_ids: List of region IDs (e.g., gene IDs) for which readout probes need to be
            assigned. This is used when generating a new codebook to ensure each region receives
            a unique barcode assignment.
        :type region_ids: list[str]
        :param file_readout_probe_table: Path to a CSV/TSV file containing the readout probe table,
            or None. The file must include columns: 'channel', 'readout_probe_id', 'L/R', and
            'readout_probe_sequence'. If a 'bit' column is not present, it will be automatically
            assigned. Cannot be None (generation of readout probe tables is not yet implemented).
        :type file_readout_probe_table: str | None
        :param file_codebook: Path to a CSV/TSV file containing an existing codebook, or None to
            generate a new codebook. If provided, the codebook must have region IDs as the index
            and bit columns named 'bit_1', 'bit_2', etc. If None, a codebook will be automatically
            generated based on the number of regions and available readout probes.
        :type file_codebook: str | None
        :return: A tuple containing (codebook, readout_probe_table), where:
            - codebook: A pandas DataFrame with region IDs as index and bit columns, where each
              row represents the barcode assignment for a region
            - readout_probe_table: A pandas DataFrame containing readout probe information with
              bit identifiers as index
        :rtype: tuple[pd.DataFrame, pd.DataFrame]
        """
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
            raise FeatureNotImplementedError(
                "Generation of readout probe table is not yet implemented. "
                "Please provide a file_readout_probe_table parameter."
            )

        if file_codebook:
            codebook = readout_probe_designer.load_codebook(file_codebook=file_codebook)
        else:
            codebook = readout_probe_designer.generate_codebook(
                region_ids=region_ids,
                n_channels=n_channels,
                n_readout_probes_LR=n_readout_probes_LR,
            )

        return codebook, readout_probe_table

    def assemble_hybridization_probes(
        self,
        target_probe_database: OligoDatabase,
        codebook: pd.DataFrame,
        readout_probe_table: pd.DataFrame,
        linker_sequence: str,
    ) -> OligoDatabase:
        """
        Assemble hybridization probes by combining target probes with readout probe sequences based on the codebook.

        This method creates the complete hybridization probe sequences by:
        1. Looking up the barcode assignment for each region in the codebook
        2. Identifying the active bits (value = 1) in the barcode, which correspond to the readout probes
           to be used for that region
        3. Retrieving the corresponding readout probe sequences (L and R) from the readout probe table
        4. Assembling the hybridization probe sequences:
           - Left probe: readout_probe_L + reverse_complement(linker) + reverse_complement(targets_sequence_L)
           - Right probe: reverse_complement(targets_sequence_R) + reverse_complement(linker) + readout_probe_R

        The assembled hybridization probes are stored in the database along with all component sequences
        (target sequences, oligo L/R sequences, readout probe sequences, and complete hybridization probes).

        :param target_probe_database: Database of target probes containing sequence and property information.
            This database should contain the designed target probes with their L and R oligo sequences
            organized by region and probe ID.
        :type target_probe_database: OligoDatabase
        :param codebook: A pandas DataFrame containing binary barcodes for each region. Each row corresponds
            to a region ID (index), and each column represents a bit in the barcode (named 'bit_1', 'bit_2', etc.).
            A value of 1 indicates that the corresponding readout probe should be used for that region.
            Each region should have exactly two active bits (one for L and one for R readout probes).
        :type codebook: pd.DataFrame
        :param readout_probe_table: A pandas DataFrame containing readout probe sequences and their associated
            bit identifiers. The DataFrame should have bit identifiers as the index and include a column
            'readout_probe_sequence' containing the probe sequences. The table should also include 'L/R'
            column to distinguish left and right readout probes.
        :type readout_probe_table: pd.DataFrame
        :param linker_sequence: DNA sequence used to link target probes and readout probes in the hybridization probe.
            This sequence is inserted between the target probe sequence and the readout probe sequence during assembly.
            Typically a short spacer sequence (e.g., "TT").
        :type linker_sequence: str
        :return: An updated `OligoDatabase` object containing the assembled hybridization probes with all
            component sequences (target, oligo L/R, readout probe L/R, and complete hybridization probe L/R)
            stored as properties for each probe.
        :rtype: OligoDatabase
        """
        region_ids = list(target_probe_database.database.keys())

        target_probe_database.set_database_sequence_types(
            [
                "sequence_target",
                "sequence_oligo_L",
                "sequence_oligo_R",
                "sequence_readout_probe_L",
                "sequence_readout_probe_R",
                "sequence_hybridization_probe_L",
                "sequence_hybridization_probe_R",
            ]
        )

        for region_id in region_ids:
            barcode = codebook.loc[region_id]
            bits = barcode[barcode == 1].index
            readout_probe_sequences = readout_probe_table.loc[bits, "readout_probe_sequence"]
            sequence_readout_probe_L = readout_probe_sequences.iloc[0]
            sequence_readout_probe_R = readout_probe_sequences.iloc[1]

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
                new_properties[probe_id]["sequence_readout_probe_L"] = sequence_readout_probe_L
                new_properties[probe_id]["sequence_readout_probe_R"] = sequence_readout_probe_R

                new_properties[probe_id]["sequence_hybridization_probe_L"] = (
                    sequence_readout_probe_L
                    + str(Seq(linker_sequence).reverse_complement())
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
                    + str(Seq(linker_sequence).reverse_complement())
                    + sequence_readout_probe_R
                )

            target_probe_database.update_oligo_properties(new_properties)

        return target_probe_database

    def design_primers(
        self,
        forward_primer_sequence: str,
        reverse_primer_sequence: str,
    ) -> tuple[str, str]:
        """
        Load and validate forward and reverse primer sequences for DNA template probe assembly.

        This method processes primer sequences that will be used for PCR amplification of the DNA
        template probes. The primers are incorporated into the DNA template probe structure during
        the assembly step, with the forward primer at the 5' end and the reverse primer at the 3' end.

        Currently, primer generation is not implemented, so both primer sequences must be provided
        as input. The method validates and strips whitespace from the sequences.

        :param forward_primer_sequence: DNA sequence of the forward primer. This primer will be
            placed at the 5' end of the DNA template probe during assembly. Must be a non-empty
            string (primer generation is not yet implemented).
        :type forward_primer_sequence: str
        :param reverse_primer_sequence: DNA sequence of the reverse primer. This primer will be
            placed at the 3' end of the DNA template probe during assembly. Must be a non-empty
            string (primer generation is not yet implemented).
        :type reverse_primer_sequence: str
        :return: A tuple containing (reverse_primer_sequence, forward_primer_sequence) in that order.
            Both sequences have been stripped of leading and trailing whitespace.
        :rtype: tuple[str, str]
        :raises FeatureNotImplementedError: If either primer sequence is empty or None, since primer
            generation is not yet implemented.
        """
        primer_designer = PrimerDesigner(
            dir_output=self.dir_output,
            n_jobs=self.n_jobs,
        )

        if forward_primer_sequence:
            forward_primer_sequence = primer_designer.load_forward_primer(
                forward_primer_sequence=forward_primer_sequence
            )
        else:
            # generate forward primers
            raise FeatureNotImplementedError(
                "Forward primer generation is not yet implemented. "
                "Please provide a forward_primer_sequence parameter."
            )

        if reverse_primer_sequence:
            reverse_primer_sequence = primer_designer.load_reverse_primer(
                reverse_primer_sequence=reverse_primer_sequence
            )
        else:
            # generate reverse primers
            raise FeatureNotImplementedError(
                "Reverse primer generation is not yet implemented. "
                "Please provide a reverse_primer_sequence parameter."
            )

        return reverse_primer_sequence, forward_primer_sequence

    def assemble_dna_template_probes(
        self,
        hybridization_probe_database: OligoDatabase,
        forward_primer_sequence: str,
        reverse_primer_sequence: str,
        linker_sequence: str,
    ) -> OligoDatabase:
        """
        Assemble DNA template probes by combining hybridization probes with forward and reverse primers.

        This method creates the final DNA template probe sequences that are used for PCR amplification
        and subsequent transcription. The DNA template probes are assembled by combining:
        - Forward primer at the 5' end
        - Target sequence
        - Linker sequence
        - Reverse complement of readout probe sequences
        - Reverse primer at the 3' end

        The assembly structure for each probe is:
        - Left DNA template probe: forward_primer + target_L + linker + reverse_complement(readout_probe_L) + reverse_primer
        - Right DNA template probe: forward_primer + reverse_complement(readout_probe_R) + linker + target_R + reverse_primer

        The assembled sequences are stored in the database along with the primer sequences for each probe.

        :param hybridization_probe_database: Database of hybridization probes containing sequence and
            property information. This database should contain the assembled hybridization probes with
            their component sequences (oligo L/R and readout probe L/R sequences).
        :type hybridization_probe_database: OligoDatabase
        :param forward_primer_sequence: DNA sequence of the forward primer that will be placed at the
            5' end of all DNA template probes.
        :type forward_primer_sequence: str
        :param reverse_primer_sequence: DNA sequence of the reverse primer that will be placed at the
            3' end of all DNA template probes.
        :type reverse_primer_sequence: str
        :param linker_sequence: DNA sequence used to link target probes and readout probes in the hybridization probe.
            This sequence is inserted between the target probe sequence and the readout probe sequence during assembly.
            Typically a short spacer sequence (e.g., "TT").
        :type linker_sequence: str
        :return: An updated `OligoDatabase` object containing the assembled DNA template probes with
            all sequences stored as properties, including: sequence_forward_primer, sequence_reverse_primer,
            sequence_dna_template_probe_L, and sequence_dna_template_probe_R for each probe.
        :rtype: OligoDatabase
        """
        region_ids = list(hybridization_probe_database.database.keys())
        hybridization_probe_database.set_database_sequence_types(
            [
                "sequence_reverse_primer",
                "sequence_forward_primer",
                "sequence_dna_template_probe_L",
                "sequence_dna_template_probe_R",
            ]
        )

        for region_id in region_ids:
            probe_ids = list(hybridization_probe_database.database[region_id].keys())
            new_properties: dict[str, dict[str, str]] = {probe_id: {} for probe_id in probe_ids}

            for probe_id in probe_ids:
                new_properties[probe_id]["sequence_reverse_primer"] = reverse_primer_sequence
                new_properties[probe_id]["sequence_forward_primer"] = forward_primer_sequence

                new_properties[probe_id]["sequence_dna_template_probe_L"] = (
                    forward_primer_sequence
                    + str(
                        Seq(
                            format_sequence(
                                database=hybridization_probe_database,
                                property="sequence_oligo_L",
                                region_id=region_id,
                                oligo_id=probe_id,
                            )
                        ).reverse_complement()
                    )
                    + linker_sequence
                    + str(
                        Seq(
                            format_sequence(
                                database=hybridization_probe_database,
                                property="sequence_readout_probe_L",
                                region_id=region_id,
                                oligo_id=probe_id,
                            )
                        ).reverse_complement()
                    )
                    + reverse_primer_sequence
                )
                new_properties[probe_id]["sequence_dna_template_probe_R"] = (
                    forward_primer_sequence
                    + str(
                        Seq(
                            format_sequence(
                                database=hybridization_probe_database,
                                property="sequence_readout_probe_R",
                                region_id=region_id,
                                oligo_id=probe_id,
                            )
                        ).reverse_complement()
                    )
                    + linker_sequence
                    + str(
                        Seq(
                            format_sequence(
                                database=hybridization_probe_database,
                                property="sequence_oligo_R",
                                region_id=region_id,
                                oligo_id=probe_id,
                            )
                        ).reverse_complement()
                    )
                    + reverse_primer_sequence
                )

            hybridization_probe_database.update_oligo_properties(new_properties)

        return hybridization_probe_database

    def generate_output(
        self,
        probe_database: OligoDatabase,
        codebook: pd.DataFrame,
        readout_probe_table: pd.DataFrame,
        output_properties: list[str] | None = None,
    ) -> None:
        """
        Generate the final output files for the CycleHCR probe design pipeline.

        This method creates all output files needed for the CycleHCR experiment, including:
        1. Codebook and readout probe tables (TSV format)
        2. Region-specific readout probe mapping (TSV format)
        3. Complete probe sets in YAML format with all properties
        4. Ready-to-order probe sequences in YAML format
        5. Probe sets in tabular format (TSV and Excel with one sheet per region)

        All output files are written to the output directory specified during pipeline initialization.

        :param probe_database: Database of DNA template probes with associated properties and sequences.
            This database should contain the final assembled probes with all component sequences
            (target, oligo L/R, readout probe L/R, hybridization probe L/R, DNA template probe L/R).
        :type probe_database: OligoDatabase
        :param codebook: A pandas DataFrame containing binary barcodes for each region. Each row
            corresponds to a region ID (index), and columns represent bits in the barcode. This
            codebook determines which readout probes are assigned to each region.
        :type codebook: pd.DataFrame
        :param readout_probe_table: A pandas DataFrame containing readout probe sequences and their
            associated bit identifiers, channels, and L/R designations. This table maps barcode
            bits to specific readout probe sequences.
        :type readout_probe_table: pd.DataFrame
        :param output_properties: List of property names to include in the output files. If None, a default set of
            properties will be included. Available properties include: 'source', 'species', 'gene_id', 'chromosome',
            'start', 'end', 'strand', 'sequence_target', 'sequence_hybridization_probe_L', 'sequence_hybridization_probe_R',
            'sequence_dna_template_probe_L', 'sequence_dna_template_probe_R', 'TmNN_sequence_target_L', etc.
        :type output_properties: list[str] | None

        :return: None

        Output Files Generated:
        -----------------------
        - ``codebook.tsv``: Binary barcode matrix for all regions
        - ``readout_probes.tsv``: Complete readout probe table with all information
        - ``readout_probes_regions.tsv``: Mapping of readout probes to regions
        - ``cyclehcr_probes.yml``: Complete probe sets with all properties in YAML format
        - ``cyclehcr_probes_order.yml``: Ready-to-order sequences (DNA template and readout probes)
        - ``cyclehcr_probes.tsv``: Probe sets in tabular format
        - ``cyclehcr_probes.xlsx``: Probe sets in Excel format with one sheet per region
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
                "sequence_target",
                "sequence_spacer",
                "sequence_readout_probe_L",
                "sequence_readout_probe_R",
                "sequence_hybridization_probe_L",
                "sequence_hybridization_probe_R",
                "sequence_forward_primer",
                "sequence_reverse_primer",
                "sequence_dna_template_probe_L",
                "sequence_dna_template_probe_R",
                "TmNN_oligo_L",
                "TmNN_oligo_R",
                "isoform_consensus",
            ]

        # write codebook and readout probe table
        codebook.to_csv(os.path.join(self.dir_output, "codebook.tsv"), sep="\t", index_label="region_id")
        readout_probe_table.to_csv(os.path.join(self.dir_output, "readout_probes.tsv"), sep="\t")

        readout_probe_table_regions = []
        for region_id, barcode in codebook.iterrows():
            bits = barcode[barcode == 1].index
            readout_probe_info = readout_probe_table.loc[bits, :]
            readout_probe_info["region_id"] = region_id
            readout_probe_table_regions.append(readout_probe_info)
        readout_probe_table_regions_df = pd.concat(readout_probe_table_regions, axis=0)
        readout_probe_table_regions_df[
            ["region_id", "channel", "readout_probe_id", "L/R", "readout_probe_sequence"]
        ].to_csv(os.path.join(self.dir_output, "readout_probes_regions.tsv"), sep="\t", index=False)

        probe_database.write_oligosets_to_yaml(
            properties=output_properties,
            ascending=True,
            filename="cyclehcr_probes",
        )

        probe_database.write_ready_to_order_yaml(
            properties=[
                "sequence_dna_template_probe_L",
                "sequence_dna_template_probe_R",
                "sequence_readout_probe_L",
                "sequence_readout_probe_R",
            ],
            ascending=True,
            filename="cyclehcr_probes_order",
        )

        probe_database.write_oligosets_to_table(
            properties=output_properties,
            ascending=True,
            filename="cyclehcr_probes",
        )


############################################
# CycleHCR Target Probe Designer
############################################


class TargetProbeDesigner:
    """
    A class for designing target probes for CycleHCR experiments.

    This class provides a comprehensive workflow for designing target probes, which are the
    gene-specific portions of CycleHCR hybridization probes. The design process includes:
    1. Creating an initial oligo database from genomic sequences using a sliding window approach
    2. Filtering probes based on sequence properties (GC content, melting temperature,
       homopolymeric runs, secondary structure)
    3. Filtering probes based on specificity to remove off-target binding and cross-hybridization
    4. Organizing filtered probes into optimal sets based on scoring criteria and constraints

    Target probes are designed as split left/right pairs that hybridize to adjacent regions
    on the target transcript, separated by a gap. This split design enables high specificity
    and allows the probes to remain bound during stripping cycles due to their high melting
    temperatures.

    :param dir_output: Directory path where output files and intermediate databases will be saved.
        The directory will be created if it does not exist.
    :type dir_output: str
    :param n_jobs: Number of parallel jobs to use for processing. This affects the parallelization
        of computationally intensive steps such as BLAST searches, property calculations, and
        filtering operations.
    :type n_jobs: int
    """

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
        """
        Create an initial oligo database by generating target probe sequences and performing pre-filtering.

        This method performs the first step of target probe design by:
        1. Generating candidate oligo sequences from input FASTA files using a sliding window approach
        2. Creating an oligo database with the generated sequences
        3. Calculating isoform consensus and reverse complement sequences for each oligo
        4. Filtering oligos based on isoform consensus threshold
        5. Splitting each oligo sequence into left probe, gap spacer, and right probe components

        The total oligo length is calculated as: L_probe_length + gap_length + R_probe_length.
        Each oligo is split into three components: the right probe (5' end), a spacer (gap region),
        and the left probe (3' end). This split is performed on the reverse complement of the target
        sequence to generate the actual probe sequences that will hybridize to the RNA.

        Regions that do not meet the minimum oligo requirement after filtering are removed from
        the database.

        :param region_ids: List of region identifiers (e.g., gene IDs) for which oligos should be
            generated. If None, all regions present in the input FASTA files will be processed.
        :type region_ids: list[str] | None
        :param target_probe_L_probe_sequence_length: Length of the left probe sequence in nucleotides.
            This is the 3' portion of the target probe that will bind to the RNA.
        :type target_probe_L_probe_sequence_length: int
        :param target_probe_gap_sequence_length: Length of the gap sequence between left and right
            probes in nucleotides. This gap is not included in the probe sequences but represents
            the spacing between the two probe halves on the target transcript.
        :type target_probe_gap_sequence_length: int
        :param target_probe_R_probe_sequence_length: Length of the right probe sequence in nucleotides.
            This is the 5' portion of the target probe that will bind to the RNA.
        :type target_probe_R_probe_sequence_length: int
        :param files_fasta_oligo_database: List of paths to FASTA files containing genomic sequences
            from which target probes will be generated. These files should contain sequences for the
            regions of interest (e.g., exons, exon-exon junctions).
        :type files_fasta_oligo_database: list[str]
        :param min_oligos_per_gene: Minimum number of oligos required per region (gene) after filtering.
            Regions with fewer oligos than this threshold will be removed from the database.
        :type min_oligos_per_gene: int
        :param isoform_consensus: Threshold for isoform consensus filtering (typically between 0.0 and 1.0).
            Probes with isoform consensus values below this threshold will be filtered out. This ensures
            that selected probes target sequences that are conserved across multiple transcript isoforms.
        :type isoform_consensus: float
        :return: An `OligoDatabase` object containing the generated target probe sequences with their
            component sequences (target, oligo, oligo_L, oligo_R, spacer) and calculated properties
            (isoform_consensus). The database is filtered to only include regions that meet the
            minimum oligo requirement.
        :rtype: OligoDatabase
        """
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
        """
        Filter the oligo database based on various sequence properties.

        This method applies multiple property-based filters to remove probes that do not meet
        quality criteria. The filters are applied sequentially to both the left (L) and right (R)
        probe sequences. Probes that fail any filter are removed from the database.

        The following filters are applied:
        1. **Hard masked sequences**: Removes probes containing hard-masked nucleotides (lowercase letters)
        2. **Homopolymeric runs**: Removes probes with homopolymeric runs exceeding the specified limits
        3. **GC content**: Removes probes with GC content outside the specified range
        4. **Melting temperature**: Removes probes with calculated Tm outside the specified range
        5. **Secondary structure**: Removes probes that form stable secondary structures at the
           specified temperature

        Regions that do not meet the minimum oligo requirement after filtering are removed from
        the database.

        :param oligo_database: The `OligoDatabase` instance containing oligonucleotide sequences
            and their associated properties. This database should contain target probes with their
            component sequences (oligo_L and oligo_R) already calculated.
        :type oligo_database: OligoDatabase
        :param GC_content_min: Minimum acceptable GC content for oligos, expressed as a fraction
            between 0.0 and 1.0 (e.g., 0.30 for 30% GC content).
        :type GC_content_min: float
        :param GC_content_max: Maximum acceptable GC content for oligos, expressed as a fraction
            between 0.0 and 1.0 (e.g., 0.90 for 90% GC content).
        :type GC_content_max: float
        :param Tm_min: Minimum acceptable melting temperature (Tm) for oligos in degrees Celsius.
            Probes with calculated Tm below this value will be filtered out.
        :type Tm_min: float
        :param Tm_max: Maximum acceptable melting temperature (Tm) for oligos in degrees Celsius.
            Probes with calculated Tm above this value will be filtered out.
        :type Tm_max: float
        :param homopolymeric_base_n: Dictionary specifying the maximum allowed length of homopolymeric
            runs for each nucleotide base. Keys should be 'A', 'T', 'G', 'C' and values are the maximum
            run length. For example: {'A': 3, 'T': 3, 'G': 3, 'C': 3} allows up to 3 consecutive
            identical bases.
        :type homopolymeric_base_n: dict[str, int]
        :param T_secondary_structure: Temperature in degrees Celsius at which to evaluate secondary
            structure formation. Secondary structures that form at this temperature can interfere
            with probe binding.
        :type T_secondary_structure: float
        :param secondary_structures_threshold_deltaG: DeltaG threshold (in kcal/mol) for secondary
            structure stability. Probes with secondary structures having deltaG values more negative
            (more stable) than this threshold will be filtered out.
        :type secondary_structures_threshold_deltaG: float
        :param Tm_parameters: Dictionary of parameters for calculating melting temperature (Tm) using
            the nearest-neighbor method. Common parameters include: 'nn_table', 'tmm_table', 'imm_table',
            'de_table', 'dnac1', 'dnac2', 'Na', 'K', 'Tris', 'Mg', 'dNTPs', 'saltcorr', etc.
        :type Tm_parameters: dict
        :param Tm_chem_correction_parameters: Dictionary of chemical correction parameters for Tm
            calculation. These parameters account for the effects of chemical additives (e.g., DMSO,
            formamide) on melting temperature. Set to None to disable chemical correction.
        :type Tm_chem_correction_parameters: dict | None
        :param Tm_salt_correction_parameters: Dictionary of salt correction parameters for Tm calculation.
            These parameters account for the effects of salt concentration on melting temperature.
            Set to None to disable salt correction.
        :type Tm_salt_correction_parameters: dict | None
        :return: A filtered `OligoDatabase` object containing only probes that pass all property filters.
            Regions with insufficient oligos after filtering are removed.
        :rtype: OligoDatabase
        """
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
        """
        Filter the oligo database based on sequence specificity to remove probes that bind
        non-specifically or cross-hybridize.

        This method applies two types of specificity filters:

        1. **Specificity filtering**: Removes probes that bind to unintended genomic regions
           - **Exact matches**: Removes all probes with exact sequence matches to probes of other regions.
           - **BLASTN specificity**: Uses BLASTN to search for similar sequences in the reference database.
             If `junction_region_size > 0`, all probes where BLASTN hits cover the junction region are
             removed, independent of the coverage threshold.

        2. **Cross-hybridization filtering**: Removes probes where the left (L) and right (R) halves
           as well as L/L and R/R probes cross-hybridize with each other. This is critical for split
           probes because if the probes can bind to each other, they may form dimers instead of binding
           to the target RNA. Probes from the larger genomic region are removed when cross-hybridization
           is detected.

        The reference database is loaded from the provided FASTA files and used for all BLASTN searches.
        Regions that do not meet the minimum oligo requirement after filtering are removed from
        the database.

        :param oligo_database: The `OligoDatabase` instance containing oligonucleotide sequences
            and their associated properties. This database should contain target probes with their
            component sequences (oligo, oligo_L, oligo_R) already calculated.
        :type oligo_database: OligoDatabase
        :param files_fasta_reference_database: List of paths to FASTA files containing reference
            sequences against which specificity will be evaluated. These typically include the
            entire genome or transcriptome to identify off-target binding sites.
        :type files_fasta_reference_database: list[str]
        :param junction_region_size: Size of the junction region (in nucleotides) for seed-based
            specificity filtering. If > 0, all probes where BLASTN hits cover the junction region
            are removed, independent of the coverage threshold.
        :type junction_region_size: int
        :param junction_site: Position of the junction site within the oligo sequence (0-based index).
            This marks the boundary between the left and right probe halves and is used for seed-based
            filtering when `junction_region_size > 0`.
        :type junction_site: int
        :param specificity_blastn_search_parameters: Dictionary of parameters for BLASTN searches
            used in specificity filtering. Common parameters include: 'task', 'word_size', 'evalue',
            'max_target_seqs', 'num_threads', etc.
        :type specificity_blastn_search_parameters: dict
        :param specificity_blastn_hit_parameters: Dictionary of parameters for filtering BLASTN hits
            in specificity searches. Common parameters include: 'identity_min', 'alignment_length_min',
            'mismatches_max', 'gaps_max', etc. Probes with hits meeting these criteria are removed.
        :type specificity_blastn_hit_parameters: dict
        :param cross_hybridization_blastn_search_parameters: Dictionary of parameters for BLASTN
            searches used in cross-hybridization filtering. These searches check if oligo_L sequences
            align to oligo_R sequences (and vice versa). Common parameters are similar to
            `specificity_blastn_search_parameters`.
        :type cross_hybridization_blastn_search_parameters: dict
        :param cross_hybridization_blastn_hit_parameters: Dictionary of parameters for filtering
            BLASTN hits in cross-hybridization searches. Common parameters are similar to
            `specificity_blastn_hit_parameters`. Probes with cross-hybridization hits meeting these
            criteria are removed from the larger region.
        :type cross_hybridization_blastn_hit_parameters: dict
        :return: A filtered `OligoDatabase` object containing only probes that pass all specificity
            and cross-hybridization filters. Regions with insufficient oligos after filtering are removed.
        :rtype: OligoDatabase
        """
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
        Tm_max: float,
        Tm_weight: float,
        Tm_parameters: dict,
        Tm_chem_correction_parameters: dict | None,
        Tm_salt_correction_parameters: dict | None,
        set_size_opt: int,
        set_size_min: int,
        distance_between_oligos: int,
        n_sets: int,
        n_attempts_graph: int,
        n_attempts_clique_enum: int,
        diversification_fraction: float,
        jaccard_opt: float,
        jaccard_step: float,
    ) -> OligoDatabase:
        """
        Create optimal oligo sets based on weighted scoring criteria, distance constraints, and set selection.

        This method selects optimal sets of target probes for each region by:
        1. Scoring each oligo based on weighted criteria (isoform consensus and melting temperature)
        2. Building a graph where edges represent non-overlapping oligos (based on distance constraints)
        3. Selecting sets of oligos that maximize the average score while respecting distance constraints
        4. Generating multiple diverse sets per region (Jaccard-controlled) to provide alternatives

        Regions that do not meet the minimum oligo requirement after set generation are removed from
        the database.

        :param oligo_database: The `OligoDatabase` instance containing oligonucleotide sequences
            and their associated properties. This database should contain filtered target probes
            ready for set selection.
        :type oligo_database: OligoDatabase
        :param isoform_weight: Weight assigned to isoform consensus in the scoring function.
            Higher values prioritize probes that are conserved across multiple transcript isoforms.
        :type isoform_weight: float
        :param Tm_max: Target melting temperature (Tm) in degrees Celsius. The scoring function
            penalizes deviations from this optimal value.
        :type Tm_max: float
        :param Tm_weight: Weight assigned to melting temperature in the scoring function.
            Higher values prioritize probes with Tm closer to the target value.
        :type Tm_weight: float
        :param Tm_parameters: Dictionary of parameters for calculating melting temperature (Tm) using
            the nearest-neighbor method. Common parameters include: 'nn_table', 'tmm_table', 'imm_table',
            'de_table', 'dnac1', 'dnac2', 'Na', 'K', 'Tris', 'Mg', 'dNTPs', 'saltcorr', etc.
        :type Tm_parameters: dict
        :param Tm_chem_correction_parameters: Dictionary of chemical correction parameters for Tm
            calculation. These parameters account for the effects of chemical additives (e.g., DMSO,
            formamide) on melting temperature. Set to None to disable chemical correction.
        :type Tm_chem_correction_parameters: dict | None
        :param Tm_salt_correction_parameters: Dictionary of salt correction parameters for Tm calculation.
            These parameters account for the effects of salt concentration on melting temperature.
            Set to None to disable salt correction.
        :type Tm_salt_correction_parameters: dict | None
        :param set_size_opt: Optimal size (number of probes) for each oligo set. The algorithm
            will attempt to generate sets of this size, but may produce sets as small as `set_size_min`
            if insufficient probes are available.
        :type set_size_opt: int
        :param set_size_min: Minimum acceptable size (number of probes) for each oligo set.
            Sets smaller than this will not be generated.
        :type set_size_min: int
        :param distance_between_oligos: Minimum genomic distance (in nucleotides) required between
            any two oligos in the same set. This ensures probes are sufficiently spaced along the
            transcript to improve hybridization efficiency.
        :type distance_between_oligos: int
        :param n_sets: Number of oligo sets to generate per region. Multiple sets provide alternatives
            in case some sets perform poorly in experiments.
        :type n_sets: int
        :param n_attempts_graph: Number of randomized graph attempts. In each attempt, a fraction of nodes is randomly
            removed from the compatibility graph to create diversity.
        :type n_attempts_graph: int
        :param n_attempts_clique_enum: Maximum number of cliques enumerated per graph attempt.
        :type n_attempts_clique_enum: int
        :param diversification_fraction: Fraction of oligos to remove at random per attempt to create
            diversity between sets.
        :type diversification_fraction: float
        :param jaccard_opt: Optimal maximum Jaccard overlap between selected sets. Sets with overlap
            above this value are discouraged when selecting multiple sets per region.
        :type jaccard_opt: float
        :param jaccard_step: Step size for relaxing Jaccard overlap when not enough sets are found.
        :type jaccard_step: float
        :return: An updated `OligoDatabase` object containing the generated oligo sets. Each region
            will have up to `n_sets` sets stored, with each set containing between `set_size_min` and
            `set_size_opt` probes. Regions with insufficient oligos are removed.
        :rtype: OligoDatabase
        """
        # Define all scorers
        isoform_consensus_scorer = IsoformConsensusScorer(
            score_weight=isoform_weight,
            property_name_transcript_id="transcript_id",
            property_name_number_total_transcripts="number_total_transcripts",
        )
        Tm_scorer = DeviationFromOptimalTmScorer(
            Tm_opt=Tm_max,
            Tm_parameters=Tm_parameters,
            Tm_salt_correction_parameters=Tm_salt_correction_parameters,
            Tm_chem_correction_parameters=Tm_chem_correction_parameters,
            score_weight=Tm_weight,
        )
        oligos_scoring = OligoScoring(scorers=[isoform_consensus_scorer, Tm_scorer])
        # the higher the score the better, because we want to have on average oligos with high melting temperatures
        set_scoring = AverageSetScoring(ascending=False)

        base_log_parameters({"Set Selection": "Independent Sets"})
        oligoset_generator = IndependentSetsOligoSelection(
            oligos_scoring=oligos_scoring,
            set_scoring=set_scoring,
            set_size_opt=set_size_opt,
            set_size_min=set_size_min,
            distance_between_oligos=distance_between_oligos,
            n_attempts_graph=n_attempts_graph,
            n_attempts_clique_enum=n_attempts_clique_enum,
            diversification_fraction=diversification_fraction,
            jaccard_opt=jaccard_opt,
            jaccard_step=jaccard_step,
        )
        oligo_database = oligoset_generator.apply(
            oligo_database=oligo_database,
            sequence_type="oligo",
            n_sets=n_sets,
            n_jobs=self.n_jobs,
        )

        oligo_database.remove_regions_with_insufficient_oligos(pipeline_step="Oligo Selection")
        check_content_oligo_database(oligo_database)

        return oligo_database


############################################
# CycleHCR Readout Probe Designer
############################################


class ReadoutProbeDesigner:
    """
    A class for managing CycleHCR readout probes and codebooks.

    This class provides methods for generating and loading codebooks (barcode matrices) that
    encode multiple genomic regions using CycleHCR readout probes. The codebook assigns each
    region a unique binary barcode, where each barcode consists of two active bits representing
    a left-right (L/R) readout probe pair in a specific fluorescence channel.

    The class also handles loading and validating readout probe tables that contain the
    sequences and metadata for the readout probes used in the encoding scheme.

    :param dir_output: Directory path where output files will be saved. This directory will
        be created if it does not exist.
    :type dir_output: str
    :param n_jobs: Number of parallel jobs to use for processing. This parameter is currently
        reserved for future parallelization of readout probe operations.
    :type n_jobs: int
    """

    def __init__(self, dir_output: str, n_jobs: int) -> None:
        """Constructor for the ReadoutProbeDesigner class."""

        ##### create the output folder #####
        self.dir_output = os.path.abspath(dir_output)
        self.n_jobs = n_jobs

    def generate_codebook(
        self, region_ids: list[str], n_channels: int, n_readout_probes_LR: int
    ) -> pd.DataFrame:
        """
        Generate a codebook (barcode matrix) for encoding multiple regions using CycleHCR readout probes.

        This method creates a binary barcode matrix where each row represents a genomic region and
        each column represents a bit position. Each region is assigned a unique barcode consisting
        of exactly two active bits (value 1), representing a left-right (L/R) readout probe pair
        in a specific fluorescence channel.

        The encoding scheme works as follows:
        - Each barcode is generated from a combination of (probe_L_id, probe_R_id, channel_id)
        - The two active bits correspond to the left and right readout probes in the specified channel
        - Combinations are prioritized: same probe pairs (L=R) are preferred over different pairs
        - The codebook size is limited by the number of available probe/channel combinations

        The method validates that sufficient barcodes are available to encode all requested regions.
        If not, a `ConfigurationError` is raised with suggestions to increase the number of probes
        or reduce the number of regions.

        :param region_ids: List of region identifiers (e.g., gene IDs) to encode in the codebook.
            Each region will be assigned a unique barcode.
        :type region_ids: list[str]
        :param n_channels: Number of fluorescence channels used in the CycleHCR experiment.
            Each channel can use different readout probe pairs.
        :type n_channels: int
        :param n_readout_probes_LR: Number of left/right readout probe pairs available per channel.
            This determines the maximum number of unique barcodes that can be generated.
        :type n_readout_probes_LR: int
        :return: A pandas DataFrame containing the binary barcode matrix. Rows are indexed by
            `region_ids`, and columns are named `bit_1`, `bit_2`, etc. Only columns with at least
            one active bit are included. Each row has exactly two bits set to 1.
        :rtype: pd.DataFrame
        :raises ConfigurationError: If the number of available barcodes is insufficient to encode
            all requested regions (i.e., `codebook_size_max < 2 * n_regions`).
        """

        def _generate_barcode(combination: tuple[int, int, int], codebook_size: int) -> list:
            index1 = ((n_channels * 2) * combination[0]) + (2 * combination[2])
            index2 = ((n_channels * 2) * combination[1]) + (2 * combination[2]) + 1
            barcode = np.zeros(codebook_size, dtype=np.int8)
            barcode[[index1, index2]] = 1
            return list(barcode)

        n_regions = len(region_ids)
        codebook_size = n_channels * n_readout_probes_LR * 2

        combinations = list(
            itertools.product(
                list(range(n_readout_probes_LR)), list(range(n_readout_probes_LR)), list(range(n_channels))
            )
        )
        combinations = sorted(combinations, key=lambda t: (0 if t[0] == t[1] else 1, t[1]))
        codebook_size_max = len(combinations)

        if codebook_size_max < (2 * n_regions):
            raise ConfigurationError(
                f"The number of valid barcodes ({codebook_size_max}) is lower than the required number of readout probes ({2 * n_regions}) for {n_regions} regions. "
                f"Consider increasing the number of L/R readout probes or reducing the number of regions."
            )

        codebook_list = []
        for combination in combinations[:n_regions]:
            barcode = _generate_barcode(
                combination=combination,
                codebook_size=codebook_size,
            )
            codebook_list.append(barcode)

        codebook: pd.DataFrame = pd.DataFrame(
            codebook_list, index=region_ids, columns=[f"bit_{i+1}" for i in range(codebook_size)]
        )

        # Remove columns where all values are 0
        codebook = codebook.loc[:, (codebook != 0).any(axis=0)]

        return codebook

    def load_codebook(self, file_codebook: str) -> pd.DataFrame:
        """
        Load and validate a codebook from a file.

        This method reads a codebook file (CSV or TSV format) and performs validation to ensure
        it meets the required format. The codebook must have:
        - A `region_id` column (or index) identifying each genomic region
        - One or more columns named with the pattern `bit_*` (e.g., `bit_1`, `bit_2`, etc.)
        - Binary values (0 or 1) in the bit columns
        - At least one row with data

        The codebook is used to assign readout probe pairs to each region based on the binary
        barcode encoding.

        :param file_codebook: Path to the CSV or TSV file containing the codebook. The file should
            have `region_id` as the index column (or a column named `region_id`), and columns named
            `bit_1`, `bit_2`, etc. representing the barcode bits.
        :type file_codebook: str
        :return: A pandas DataFrame containing the codebook with region IDs as the index and
            bit columns as data columns. The DataFrame is filtered to only include bit columns
            that have at least one active bit (value 1).
        :rtype: pd.DataFrame
        :raises FileFormatError: If the codebook file:
            - Does not contain at least one column
            - Contains columns that are not named with the `bit_*` pattern
            - Does not contain at least one row with data (after removing empty rows)
        """
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

        # Check for at least one data row (excluding empty rows)
        codebook_clean = codebook.dropna(how="all")
        if len(codebook_clean) == 0:
            raise FileFormatError(f"Codebook file '{file_codebook}' must contain at least one row with data.")

        return codebook

    def load_readout_probe_table(self, file_readout_probe_table: str) -> tuple[pd.DataFrame, int, int]:
        """
        Load and validate a table containing readout probe information.

        This method reads a readout probe table from a file and validates its structure. The table
        must contain the following required columns:
        - `channel`: Fluorescence channel number (integer)
        - `readout_probe_id`: Unique identifier for each readout probe (within a channel)
        - `L/R`: Probe type, either 'L' (left) or 'R' (right)
        - `readout_probe_sequence`: DNA sequence of the readout probe

        If a `bit` column is not present, the method automatically assigns bit labels (`bit_1`,
        `bit_2`, etc.) based on the sorted order of probes by `readout_probe_id` and `channel`.
        The bit labels are used to map probes to positions in the codebook.

        The method calculates the number of channels and the number of L/R probe pairs per channel
        by analyzing the data. It assumes an equal number of L and R probes per channel.

        :param file_readout_probe_table: Path to the CSV or TSV file containing the readout probe
            data. The file should have columns: `channel`, `readout_probe_id`, `L/R`, and
            `readout_probe_sequence`. An optional `bit` column can be included to specify bit
            labels manually.
        :type file_readout_probe_table: str
        :return: A tuple containing:
            - **DataFrame**: The formatted readout probe table with `bit` as the index and the
              required columns as data columns. If `bit` was not in the original file, it is
              automatically generated.
            - **int**: Number of unique fluorescence channels in the table.
            - **int**: Number of left/right readout probe pairs per channel (calculated as the
              minimum of L and R probes divided by the number of channels).
        :rtype: tuple[pd.DataFrame, int, int]
        :raises FileFormatError: If the readout probe table is missing any of the required columns:
            `channel`, `readout_probe_id`, `L/R`, or `readout_probe_sequence`.
        """
        required_cols = ["channel", "readout_probe_id", "L/R", "readout_probe_sequence"]

        readout_probe_table = pd.read_csv(file_readout_probe_table, sep=None, engine="python")

        # Check if all required columns exist in readout_probe_table
        cols = set(readout_probe_table.columns)
        if not set(required_cols).issubset(cols):
            missing = set(required_cols) - cols
            raise FileFormatError(
                f"Readout probe table is missing required columns: {missing}. "
                f"Required columns are: {required_cols}."
            )

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
    """
    A class for loading and validating CycleHCR PCR primers.

    This class provides methods for loading and validating forward and reverse primer sequences
    that are used in the CycleHCR pipeline for PCR amplification of DNA template probes. The primers
    bind to the 5' and 3' ends of the DNA template probes and enable amplification during the
    experimental workflow.

    Currently, this class only supports loading pre-designed primer sequences. Automatic primer
    design functionality is not yet implemented.

    :param dir_output: Directory path where output files will be saved. This directory will
        be created if it does not exist.
    :type dir_output: str
    :param n_jobs: Number of parallel jobs to use for processing. This parameter is currently
        reserved for future parallelization of primer operations.
    :type n_jobs: int
    """

    def __init__(self, dir_output: str, n_jobs: int) -> None:
        """Constructor for the PrimerDesigner class."""

        ##### create the output folder #####
        self.dir_output = os.path.abspath(dir_output)
        self.n_jobs = n_jobs

    def load_forward_primer(self, forward_primer_sequence: str) -> str:
        """
        Load and validate a forward primer sequence.

        This method takes a forward primer sequence string, validates it by converting to string
        and stripping whitespace, then returns the cleaned sequence. The forward primer binds to
        the 5' end of the DNA template probe and initiates PCR amplification in the forward direction.

        :param forward_primer_sequence: DNA sequence of the forward primer. Should be a string
            containing valid nucleotide characters (A, T, G, C). The sequence will be stripped
            of leading and trailing whitespace.
        :type forward_primer_sequence: str
        :return: The cleaned forward primer sequence with whitespace removed. The sequence is
            ready to be used in the CycleHCR pipeline for DNA template probe assembly.
        :rtype: str
        """
        forward_primer = str(forward_primer_sequence).strip()
        return forward_primer

    def load_reverse_primer(self, reverse_primer_sequence: str) -> str:
        """
        Load and validate a reverse primer sequence.

        This method takes a reverse primer sequence string, validates it by converting to string
        and stripping whitespace, then returns the cleaned sequence. The reverse primer binds to
        the 3' end of the DNA template probe and initiates PCR amplification in the reverse direction.

        :param reverse_primer_sequence: DNA sequence of the reverse primer. Should be a string
            containing valid nucleotide characters (A, T, G, C). The sequence will be stripped
            of leading and trailing whitespace.
        :type reverse_primer_sequence: str
        :return: The cleaned reverse primer sequence with whitespace removed. The sequence is
            ready to be used in the CycleHCR pipeline for DNA template probe assembly.
        :rtype: str
        """
        reverse_primer = str(reverse_primer_sequence).strip()
        return reverse_primer


############################################
# CycleHCR Probe Designer Pipeline
############################################


def main() -> None:
    """
    Main entry point for running the CycleHCR probe design pipeline.

    This function orchestrates the complete CycleHCR probe design workflow:
    1. Parses command-line arguments using the base parser
    2. Reads the configuration YAML file containing all pipeline parameters
    3. Reads the gene IDs file (if provided) or uses all genes from FASTA files
    4. Preprocesses melting temperature parameters for target probes
    5. Initializes the CycleHCRProbeDesigner pipeline
    6. Designs target probes for specified genes
    7. Loads readout probes and generates the codebook
    8. Assembles hybridization probes by combining target probes with readout probe barcodes
    9. Loads/validates forward and reverse primers for PCR amplification
    10. Assembles final DNA template probes with primers
    11. Generates output files (codebook, readout probe table, probe sequences, etc.)


    The function is typically called from the command line:
    ``cycle_hcr_probe_designer --config <path_to_config.yaml>``

    Command-line arguments are parsed using `base_parser()`, which expects:
    - `config`: Path to the YAML configuration file containing all pipeline parameters
    """
    logging.info("--------------START PIPELINE--------------")

    args = base_parser()

    ##### read the config file #####
    with open(args["config"], "r") as handle:
        config = yaml.safe_load(handle)

    ##### read the genes file #####
    if config["file_regions"] is None:
        warnings.warn(
            "No gene list file was provided! All genes from fasta file are used to generate the probes. This chioce can use a lot of resources."
        )
        region_ids = None
    else:
        with open(config["file_regions"]) as handle:
            lines = handle.readlines()
            # ensure that the list contains unique gene ids
            region_ids = list(set([line.rstrip() for line in lines]))

    ##### initialize probe designer pipeline #####
    pipeline = CycleHCRProbeDesigner(
        dir_output=config["dir_output"],
        write_intermediate_steps=config["write_intermediate_steps"],
        n_jobs=config["n_jobs"],
    )

    ##### design probes #####
    target_probe_database = pipeline.design_target_probes(
        region_ids=region_ids,
        files_fasta_target_probe_database=config["files_fasta_target_probe_database"],
        files_fasta_reference_database_target_probe=config["files_fasta_reference_database_target_probe"],
        # Target Probe Design
        target_probe_isoform_consensus=config["target_probe_isoform_consensus"],
        target_probe_L_probe_sequence_length=config["target_probe_L_probe_sequence_length"],
        target_probe_gap_sequence_length=config["target_probe_gap_sequence_length"],
        target_probe_R_probe_sequence_length=config["target_probe_R_probe_sequence_length"],
        # Property Filter Parameters
        target_probe_GC_content_min=config["target_probe_GC_content_min"],
        target_probe_GC_content_max=config["target_probe_GC_content_max"],
        target_probe_Tm_min=config["target_probe_Tm_min"],
        target_probe_Tm_max=config["target_probe_Tm_max"],
        target_probe_homopolymeric_base_n=config["target_probe_homopolymeric_base_n"],
        target_probe_T_secondary_structure=config["target_probe_T_secondary_structure"],
        target_probe_secondary_structures_threshold_deltaG=config[
            "target_probe_secondary_structures_threshold_deltaG"
        ],
        # Melting Temperature Calculation Parameters
        target_probe_Tm_parameters=preprocess_tm_parameters(config["target_probe_Tm_parameters"]),
        target_probe_Tm_chem_correction_parameters=config["target_probe_Tm_chem_correction_parameters"],
        target_probe_Tm_salt_correction_parameters=config["target_probe_Tm_salt_correction_parameters"],
        # Specificity Filter Parameters
        target_probe_junction_region_size=config["target_probe_junction_region_size"],
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
        # Probe Scoring and Set Selection Parameters
        target_probe_Tm_weight=config["target_probe_Tm_weight"],
        target_probe_isoform_weight=config["target_probe_isoform_weight"],
        set_size_opt=config["set_size_opt"],
        set_size_min=config["set_size_min"],
        distance_between_target_probes=config["distance_between_target_probes"],
        n_sets=config["n_sets"],
        n_attempts_graph=config["n_attempts_graph"],
        n_attempts_clique_enum=config["n_attempts_clique_enum"],
        diversification_fraction=config["diversification_fraction"],
        jaccard_opt=config["jaccard_opt"],
        jaccard_step=config["jaccard_step"],
    )

    codebook, readout_probe_table = pipeline.design_readout_probes(
        region_ids=list(target_probe_database.database.keys()),
        file_readout_probe_table=config["file_readout_probe_table"],
        file_codebook=config["file_codebook"],
    )

    hybridization_probe_database = pipeline.assemble_hybridization_probes(
        target_probe_database=target_probe_database,
        codebook=codebook,
        readout_probe_table=readout_probe_table,
        linker_sequence=config["linker_sequence"],
    )

    reverse_primer_sequence, forward_primer_sequence = pipeline.design_primers(
        forward_primer_sequence=config["forward_primer_sequence"],
        reverse_primer_sequence=config["reverse_primer_sequence"],
    )

    final_probe_database = pipeline.assemble_dna_template_probes(
        hybridization_probe_database=hybridization_probe_database,
        forward_primer_sequence=forward_primer_sequence,
        reverse_primer_sequence=reverse_primer_sequence,
        linker_sequence=config["linker_sequence"],
    )

    pipeline.generate_output(
        probe_database=final_probe_database,
        codebook=codebook,
        readout_probe_table=readout_probe_table,
    )

    logging.info("--------------END PIPELINE--------------")


if __name__ == "__main__":
    main()
