############################################
# imports
############################################

import logging
import os
import shutil
import warnings
from itertools import product
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from Bio.SeqUtils import Seq

from oligo_designer_toolsuite._exceptions import ConfigurationError
from oligo_designer_toolsuite.database import OligoDatabase, ReferenceDatabase
from oligo_designer_toolsuite.oligo_efficiency_filter import (
    DeviationFromOptimalGCContentScorer,
    LowestSetScoring,
    OligoScoring,
    OverlapUTRScorer,
)
from oligo_designer_toolsuite.oligo_property_calculator import (
    IsoformConsensusProperty,
    PropertyCalculator,
    ReverseComplementSequenceProperty,
)
from oligo_designer_toolsuite.oligo_property_calculator._property_functions import calc_tm_nn
from oligo_designer_toolsuite.oligo_property_filter import (
    ComplementFilter,
    GCClampFilter,
    GCContentFilter,
    HardMaskedSequenceFilter,
    HomopolymericRunsFilter,
    MeltingTemperatureNNFilter,
    PropertyFilter,
    SecondaryStructureFilter,
    SelfComplementFilter,
    SoftMaskedSequenceFilter,
)
from oligo_designer_toolsuite.oligo_selection import IndependentSetsOligoSelection
from oligo_designer_toolsuite.oligo_specificity_filter import (
    BlastNFilter,
    CrossHybridizationFilter,
    ExactMatchFilter,
    RemoveAllFilterPolicy,
    RemoveByDegreeFilterPolicy,
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
from oligo_designer_toolsuite.utils import append_nucleotide_to_sequences

############################################
# SeqFISH Plus Probe Designer
############################################


class SeqFishPlusProbeDesigner:
    """
    A class for designing hybridization probes for SeqFISH+ (sequential fluorescence in situ hybridization plus) experiments.

    This class provides a complete pipeline for designing SeqFISH+ probes, which enable multiplexed RNA detection
    in single cells through combinatorial barcoding with pseudocolors and sequential imaging rounds.

    **SeqFISH+ Pipeline Overview:**
    1. **Target Probe Design**: Design gene-specific targeting sequences (28 nt) that bind to RNA transcripts
    2. **Readout Probe Design**: Generate readout probe sequences (15 nt) and create a codebook using pseudocolors
       and barcode rounds for hybridization
    3. **Hybridization Probe Assembly**: Combine target probes with four readout probe overhang sequences based on the codebook
    4. **Primer Design**: Design PCR primers for amplifying DNA template probes
    5. **DNA Template Probe Assembly**: Assemble final DNA template probes with forward and reverse primers
    6. **Output Generation**: Generate output files in multiple formats (TSV, YAML)

    Overview
    --------
    seqFISH+ (sequential fluorescence in situ hybridization plus) is a spatial transcriptomics
    method that enables the identification, counting, and spatial mapping of up to 10,000 RNA
    species directly in single cells within their native tissue context. By combining sequential
    rounds of hybridization and imaging, seqFISH+ overcomes the optical crowding that limits
    traditional smFISH approaches and achieves super-resolved, transcriptome-scale RNA imaging
    using a standard confocal microscope.

    Each RNA species is assigned a unique multi-bit barcode that is read out through sequential
    hybridization of fluorophore-labeled readout probes. By expanding the barcode space using
    “pseudocolours” across multiple fluorescence channels and hybridization rounds, seqFISH+
    achieves genome-scale multiplexing.

    Probe Structure
    ---------------
    **Hybridization (Encoding) Probes**
    - Single-stranded DNA oligonucleotides that hybridize directly to target mRNAs.
    - Each probe contains:
        - A **28-nt targeting sequence** complementary to the RNA transcript.
        - Four **15-nt overhang sequences** (readout handles, labeled I–IV) that correspond to
        barcode positions decoded by sequential readout hybridizations.
        - Overhangs serve as binding sites for fluorescent readout probes during barcode readout.
    - A set of 24 hybridization probes is typically designed per gene to ensure detection efficiency
    and signal amplification.
    - The hybridization probe has the structure:
        [Overhang I] + [Overhang II] + [Targeting Sequence] + [Overhang III] + [Overhang IV]

    **Readout Probes**
    - Short (typically 15-nt) dye-labeled DNA oligonucleotides complementary to the overhang
    sequences on the hybridization probes.
    - Each readout probe binds to one overhang and carries a fluorophore (Alexa Fluor 488, Cy3b,
    or Alexa Fluor 647) that reports the “on” state of a barcode bit during a hybridization round.
    - The barcoding scheme uses a combinatorial “pseudocolour” system with 60 hybridization
    rounds divided among three fluorescence channels. Each channel encodes ~8,000 genes,
    resulting in a total of ~24,000 barcodes with built-in error correction.

    **DNA Template Probe**
    - The hybridization probes are synthesized from large oligonucleotide pools containing
    forward and reverse PCR primer regions flanking the probe body.
    - These primer regions allow limited-cycle PCR amplification and in vitro transcription
    to generate RNA intermediates that are reverse-transcribed into single-stranded DNA probes.
    - The forward primer includes a uracil residue, enabling enzymatic cleavage (e.g., by USER enzyme)
    to remove primer sequences and yield the final probe.
    - The **DNA template probe** has the structure:
        [Forward Primer] + [Overhang I] + [Overhang II] + [Targeting Sequence] + [Overhang III] + [Overhang IV] + [Reverse Primer]

    Probe Library Preparation
    -------------------------
    The seqFISH+ probe library is generated through a multi-step molecular synthesis and
    purification workflow. First, target genes are selected and assigned barcode sequences
    within an error-correctable codebook. For each gene, a set of 24 hybridization probes is designed,
    each with a 28-nt target-binding region and four 15-nt readout overhangs. These sequences
    are synthesized in an oligonucleotide pool and amplified by limited PCR cycles. Amplified
    products are transcribed and reverse-transcribed to generate single-stranded DNA probes,
    followed by enzymatic cleavage to remove primer regions and purification.

    References
    ----------
    Eng, C. H. L., Lawson, M., Zhu, Q., Dries, R., Koulena, N., Takei, Y., ... & Cai, L. (2019).
    Transcriptome-scale super-resolved imaging in tissues by RNA seqFISH+.
    Nature, 568(7751), 235-239.

    :param dir_output: Directory path where output files will be saved. This directory will be created
        if it does not exist.
    :type dir_output: str
    :param write_intermediate_steps: Whether to save intermediate results during the probe design pipeline.
        If True, intermediate databases and results will be saved at each pipeline step, which is useful
        for debugging and analysis but increases disk usage.
    :type write_intermediate_steps: bool
    :param n_jobs: Number of parallel jobs to use for processing. Set to 1 for serial processing or higher
        values for parallel processing. This affects the parallelization of filtering, property calculation,
        and set generation operations.
    :type n_jobs: int
    """

    def __init__(
        self,
        write_intermediate_steps: bool,
        dir_output: str,
        n_jobs: int,
    ) -> None:
        """Constructor for the SeqFishPlusProbeDesigner class."""

        # create the output folder
        self.dir_output = os.path.abspath(dir_output)
        Path(self.dir_output).mkdir(parents=True, exist_ok=True)

        # setup logger
        setup_logging(
            dir_output=self.dir_output,
            pipeline_name="seqfishplus_probe_designer",
            log_start_message=True,
        )

        ##### set class parameters #####
        self.write_intermediate_steps = write_intermediate_steps
        self.n_jobs = n_jobs

    def design_target_probes(
        self,
        # Step 1: Create Database Parameters
        region_ids: list[str] | None,
        files_fasta_target_probe_database: list[str],
        target_probe_length_min: int,
        target_probe_length_max: int,
        target_probe_isoform_consensus: float,
        # Step 2: Property Filter Parameters
        target_probe_GC_content_min: float,
        target_probe_GC_content_opt: float,
        target_probe_GC_content_max: float,
        target_probe_homopolymeric_base_n: dict[str, int],
        target_probe_T_secondary_structure: float,
        target_probe_secondary_structures_threshold_deltaG: float,
        # Step 3: Specificity Filter Parameters
        files_fasta_reference_database_target_probe: list[str],
        target_probe_specificity_blastn_search_parameters: dict,
        target_probe_specificity_blastn_hit_parameters: dict,
        target_probe_cross_hybridization_blastn_search_parameters: dict,
        target_probe_cross_hybridization_blastn_hit_parameters: dict,
        # Step 4: Probe Scoring and Set Selection Parameters
        target_probe_GC_weight: float,
        target_probe_UTR_weight: float,
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
        Design target probes for SeqFISH+ experiments through a multi-step pipeline.

        This method performs the complete target probe design process, which includes:
        1. Creating an initial oligo database from input FASTA files using a sliding window approach
        2. Filtering probes based on sequence properties (GC content, homopolymeric runs, secondary structure)
        3. Filtering probes based on specificity to remove off-target binding and cross-hybridization
           using BLASTN searches
        4. Organizing filtered probes into optimal sets based on weighted scoring criteria (GC content,
           UTR targeting) and distance constraints

        The resulting probes are gene-specific targeting sequences (typically 28 nt) that bind to RNA
        transcripts. These probes will later be combined with four readout probe overhang sequences to
        create complete hybridization probes.

        **Step 1: Create Database Parameters**

        :param region_ids: List of gene identifiers (e.g., gene IDs) to target for probe design. If None,
            all genes present in the input FASTA files will be used.
        :type region_ids: list[str] | None
        :param files_fasta_target_probe_database: List of paths to FASTA files containing sequences
            from which target probes will be generated. These files should contain genomic regions
            of interest (e.g., exons, exon-exon junctions).
        :type files_fasta_target_probe_database: list[str]
        :param target_probe_length_min: Minimum length (in nucleotides) for target probe sequences.
        :type target_probe_length_min: int
        :param target_probe_length_max: Maximum length (in nucleotides) for target probe sequences.
        :type target_probe_length_max: int
        :param target_probe_isoform_consensus: Threshold for isoform consensus filtering (typically
            between 0.0 and 1.0). Probes with isoform consensus values below this threshold will be
            filtered out. This ensures that selected probes target sequences that are conserved across
            multiple transcript isoforms.
        :type target_probe_isoform_consensus: float

        **Step 2: Property Filter Parameters**

        :param target_probe_GC_content_min: Minimum acceptable GC content for target probes, expressed
            as a fraction between 0.0 and 1.0.
        :type target_probe_GC_content_min: float
        :param target_probe_GC_content_opt: Optimal GC content for target probes, expressed as a fraction
            between 0.0 and 1.0. Used in scoring to prioritize probes closer to this value.
        :type target_probe_GC_content_opt: float
        :param target_probe_GC_content_max: Maximum acceptable GC content for target probes, expressed
            as a fraction between 0.0 and 1.0.
        :type target_probe_GC_content_max: float
        :param target_probe_homopolymeric_base_n: Dictionary specifying the maximum allowed length of
            homopolymeric runs for each nucleotide base (keys: 'A', 'T', 'G', 'C').
        :type target_probe_homopolymeric_base_n: dict[str, int]
        :param target_probe_T_secondary_structure: Temperature in degrees Celsius at which to evaluate
            secondary structure formation.
        :type target_probe_T_secondary_structure: float
        :param target_probe_secondary_structures_threshold_deltaG: DeltaG threshold (in kcal/mol) for
            secondary structure stability. Probes with secondary structures having deltaG values more
            negative than this threshold will be filtered out.
        :type target_probe_secondary_structures_threshold_deltaG: float

        **Step 3: Specificity Filter Parameters**

        :param files_fasta_reference_database_target_probe: List of paths to FASTA files containing
            reference sequences used for specificity filtering. These files are used to identify
            off-target binding sites (e.g., whole genome or transcriptome sequences).
        :type files_fasta_reference_database_target_probe: list[str]
        :param target_probe_specificity_blastn_search_parameters: Dictionary of parameters for BLASTN
            searches used in specificity filtering.
        :type target_probe_specificity_blastn_search_parameters: dict
        :param target_probe_specificity_blastn_hit_parameters: Dictionary of parameters for filtering
            BLASTN hits in specificity searches.
        :type target_probe_specificity_blastn_hit_parameters: dict
        :param target_probe_cross_hybridization_blastn_search_parameters: Dictionary of parameters for
            BLASTN searches used in cross-hybridization filtering.
        :type target_probe_cross_hybridization_blastn_search_parameters: dict
        :param target_probe_cross_hybridization_blastn_hit_parameters: Dictionary of parameters for
            filtering BLASTN hits in cross-hybridization searches.
        :type target_probe_cross_hybridization_blastn_hit_parameters: dict

        **Step 4: Probe Scoring and Set Selection Parameters**

        :param target_probe_GC_weight: Weight assigned to GC content in the scoring function.
        :type target_probe_GC_weight: float
        :param target_probe_UTR_weight: Weight assigned to untranslated region (UTR) targeting in the
            scoring function. Probes that target UTRs are prioritized.
        :type target_probe_UTR_weight: float
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
        :param diversification_fraction: Fraction of oligos to remove at random per attempt to create diversity
            between sets.
        :type diversification_fraction: float
        :param jaccard_opt: Optimal maximum Jaccard overlap between selected sets. Sets with overlap above this
            value are discouraged when selecting multiple sets per region.
        :type jaccard_opt: float
        :param jaccard_step: Step size for relaxing Jaccard overlap when not enough sets are found.
        :type jaccard_step: float
        :return: An `OligoDatabase` object containing the designed target probes organized into sets.
            The database includes probe sequences, properties, and set assignments for each target gene.
        :rtype: OligoDatabase
        """
        target_probe_designer = TargetProbeDesigner(self.dir_output, self.n_jobs)

        oligo_database: OligoDatabase = target_probe_designer.create_oligo_database(
            region_ids=region_ids,
            oligo_length_min=target_probe_length_min,
            oligo_length_max=target_probe_length_max,
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
            homopolymeric_base_n=target_probe_homopolymeric_base_n,
            T_secondary_structure=target_probe_T_secondary_structure,
            secondary_structures_threshold_deltaG=target_probe_secondary_structures_threshold_deltaG,
        )

        if self.write_intermediate_steps:
            dir_database = oligo_database.save_database(name_database="2_db_target_probes_property_filter")
            logging.info(
                f"Saved target probe database for step 2 (Property Filters) in directory {dir_database}"
            )

        oligo_database = target_probe_designer.filter_by_specificity(
            oligo_database=oligo_database,
            files_fasta_reference_database=files_fasta_reference_database_target_probe,
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
            # Scoring Parameters
            GC_weight=target_probe_GC_weight,
            GC_content_opt=target_probe_GC_content_opt,
            UTR_weight=target_probe_UTR_weight,
            # Set Selection Parameters
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

        if self.write_intermediate_steps:
            dir_database = oligo_database.save_database(name_database="4_db_target_probes_sets")
            logging.info(
                f"Saved target probe database for step 4 (Specificity Filters) in directory {dir_database}."
            )

        return oligo_database

    def design_readout_probes(
        self,
        # Step 1: Create Database Parameters
        region_ids: list[str],
        readout_probe_length: int,
        readout_probe_base_probabilities: dict[str, float],
        readout_probe_initial_num_sequences: int,
        # Step 2: Property Filter Parameters
        readout_probe_GC_content_min: float,
        readout_probe_GC_content_max: float,
        readout_probe_homopolymeric_base_n: dict[str, int],
        # Step 3: Specificity Filter Parameters
        files_fasta_reference_database_readout_probe: list[str],
        readout_probe_specificity_blastn_search_parameters: dict,
        readout_probe_specificity_blastn_hit_parameters: dict,
        readout_probe_cross_hybridization_blastn_search_parameters: dict,
        readout_probe_cross_hybridization_blastn_hit_parameters: dict,
        # Step 4: Codebook and Readout Probe Table Parameters
        n_barcode_rounds: int,
        n_pseudocolors: int,
        channels_ids: list,
    ) -> tuple[pd.DataFrame, pd.DataFrame]:
        """
        Design readout probes and generate a codebook for SeqFISH+ experiments through a multi-step pipeline.

        This method performs the complete readout probe design process, which includes:
        1. Creating an initial oligo database by generating random sequences with specified nucleotide
           probabilities and length
        2. Filtering probes based on sequence properties (GC content, homopolymeric runs)
        3. Filtering probes based on specificity to remove off-target binding and cross-hybridization
           using BLASTN searches
        4. Generating a codebook using pseudocolors and barcode rounds for combinatorial hybridization
        5. Creating a readout probe table that maps codebook bits to channels, barcode rounds,
           pseudocolors, and readout probe sequences

        The codebook assigns each region a unique barcode using a combinatorial pseudocolor system
        across multiple barcode rounds and fluorescence channels. The readout probe table assigns
        readout probes to each bit position and distributes them across channels, rounds, and
        pseudocolors for multiplexed detection.

        **Step 1: Create Database Parameters**

        :param region_ids: List of region identifiers (e.g., gene IDs) for which readout probes and
            codebook entries are to be generated. The number of regions determines the minimum number
            of barcodes required in the codebook.
        :type region_ids: list[str]
        :param readout_probe_length: Length (in nucleotides) of each readout probe sequence.
        :type readout_probe_length: int
        :param readout_probe_base_probabilities: Dictionary specifying the probability of each nucleotide
            base in randomly generated sequences. Keys should be 'A', 'T', 'G', 'C' and values should
            sum to 1.0 (e.g., {"A": 0.25, "T": 0.25, "G": 0.25, "C": 0.25}).
        :type readout_probe_base_probabilities: dict[str, float]
        :param readout_probe_initial_num_sequences: Number of random sequences to generate initially
            before filtering. Higher values provide more candidates but increase computation time.
        :type readout_probe_initial_num_sequences: int

        **Step 2: Property Filter Parameters**

        :param readout_probe_GC_content_min: Minimum acceptable GC content for readout probes, expressed
            as a fraction between 0.0 and 1.0.
        :type readout_probe_GC_content_min: float
        :param readout_probe_GC_content_max: Maximum acceptable GC content for readout probes, expressed
            as a fraction between 0.0 and 1.0.
        :type readout_probe_GC_content_max: float
        :param readout_probe_homopolymeric_base_n: Dictionary specifying the maximum allowed length of
            homopolymeric runs for each nucleotide base (keys: 'A', 'T', 'G', 'C').
        :type readout_probe_homopolymeric_base_n: dict[str, int]

        **Step 3: Specificity Filter Parameters**

        :param files_fasta_reference_database_readout_probe: List of paths to FASTA files containing
            reference sequences used for specificity filtering. These files are used to identify
            off-target binding sites and potential cross-hybridization events.
        :type files_fasta_reference_database_readout_probe: list[str]
        :param readout_probe_specificity_blastn_search_parameters: Dictionary of parameters for BLASTN
            searches used in specificity filtering.
        :type readout_probe_specificity_blastn_search_parameters: dict
        :param readout_probe_specificity_blastn_hit_parameters: Dictionary of parameters for filtering
            BLASTN hits in specificity searches.
        :type readout_probe_specificity_blastn_hit_parameters: dict
        :param readout_probe_cross_hybridization_blastn_search_parameters: Dictionary of parameters for
            BLASTN searches used in cross-hybridization filtering.
        :type readout_probe_cross_hybridization_blastn_search_parameters: dict
        :param readout_probe_cross_hybridization_blastn_hit_parameters: Dictionary of parameters for
            filtering BLASTN hits in cross-hybridization searches.
        :type readout_probe_cross_hybridization_blastn_hit_parameters: dict

        **Step 4: Codebook and Readout Probe Table Parameters**

        :param n_barcode_rounds: Number of barcode rounds for hybridization. Each round contributes to the
            combinatorial barcode space.
        :type n_barcode_rounds: int
        :param n_pseudocolors: Number of pseudocolors to use in the barcoding scheme. Pseudocolors
            expand the barcode space by allowing multiple states per round.
        :type n_pseudocolors: int
        :param channels_ids: List of fluorescence channel identifiers (e.g., ['Alexa488', 'Cy3b', 'Alexa647'])
            to which readout probes will be assigned. Readout probes are distributed across channels
            in a round-robin fashion.
        :type channels_ids: list
        :return: A tuple containing:
            - **codebook** (pd.DataFrame): Binary barcode matrix with region IDs as index and bit columns
              (bit_1, bit_2, etc.) as data. Each row represents a region's barcode assignment.
            - **readout_probe_table** (pd.DataFrame): Table mapping each bit to a readout probe sequence,
              barcode round, pseudocolor, channel assignment, and bit label. Indexed by bit labels
              (bit_1, bit_2, etc.).
        :rtype: tuple[pd.DataFrame, pd.DataFrame]
        """
        readout_probe_designer = ReadoutProbeDesigner(
            dir_output=self.dir_output,
            n_jobs=self.n_jobs,
        )
        oligo_database: OligoDatabase = readout_probe_designer.create_oligo_database(
            oligo_length=readout_probe_length,
            oligo_base_probabilities=readout_probe_base_probabilities,
            initial_num_sequences=readout_probe_initial_num_sequences,
        )

        if self.write_intermediate_steps:
            dir_database = oligo_database.save_database(name_database="1_db_readout_probes_initial")
            logging.info(
                f"Saved readout probe database for step 1 (Create Database) in directory {dir_database}"
            )

        oligo_database = readout_probe_designer.filter_by_property(
            oligo_database=oligo_database,
            GC_content_min=readout_probe_GC_content_min,
            GC_content_max=readout_probe_GC_content_max,
            homopolymeric_base_n=readout_probe_homopolymeric_base_n,
        )

        if self.write_intermediate_steps:
            dir_database = oligo_database.save_database(name_database="2_db_readout_probes_property_filter")
            logging.info(
                f"Saved readout probe database for step 2 (Property Filters) in directory {dir_database}"
            )

        oligo_database = readout_probe_designer.filter_by_specificity(
            oligo_database=oligo_database,
            files_fasta_reference_database=files_fasta_reference_database_readout_probe,
            specificity_blastn_search_parameters=readout_probe_specificity_blastn_search_parameters,
            specificity_blastn_hit_parameters=readout_probe_specificity_blastn_hit_parameters,
            cross_hybridization_blastn_search_parameters=readout_probe_cross_hybridization_blastn_search_parameters,
            cross_hybridization_blastn_hit_parameters=readout_probe_cross_hybridization_blastn_hit_parameters,
        )

        if self.write_intermediate_steps:
            dir_database = oligo_database.save_database(name_database="3_db_readout_probes_specificty_filter")
            logging.info(
                f"Saved readout probe database for step 3 (Specificity Filters) in directory {dir_database}"
            )

        codebook = readout_probe_designer.generate_codebook(
            region_ids=region_ids,
            n_barcode_rounds=n_barcode_rounds,
            n_pseudocolors=n_pseudocolors,
            n_channels=len(channels_ids),
        )

        readout_probe_table = readout_probe_designer.create_readout_probe_table(
            readout_probe_database=oligo_database,
            channels_ids=channels_ids,
            n_barcode_rounds=n_barcode_rounds,
            n_pseudocolors=n_pseudocolors,
        )

        return codebook, readout_probe_table

    def assemble_hybridization_probes(
        self,
        target_probe_database: OligoDatabase,
        codebook: pd.DataFrame,
        readout_probe_table: pd.DataFrame,
    ) -> OligoDatabase:
        """
        Assemble hybridization probes by combining target probes with four readout probe overhang sequences based on the codebook.

        This method creates complete SeqFISH+ hybridization probes by combining gene-specific target probes
        with four readout probe overhang sequences according to the codebook assignment. For each region, the method:
        1. Looks up the region's barcode in the codebook to identify which four readout probes are assigned
        2. Retrieves the corresponding readout probe sequences from the readout probe table
        3. Assembles the hybridization probe sequence with the structure:
           [reverse_complement(readout_probe_1)] + [reverse_complement(readout_probe_2)] + "T" +
           [target_probe] + "T" + [reverse_complement(readout_probe_3)] + [reverse_complement(readout_probe_4)]

        The readout probes are reverse-complemented because they will hybridize to the overhang sequences
        embedded in the hybridization probe. The single "T" nucleotides serve as spacers between the
        readout probe binding sites and the target probe sequence.

        The assembled sequences are stored as properties in the database for each probe, enabling downstream
        primer addition and DNA template probe assembly.

        :param target_probe_database: The `OligoDatabase` instance containing target probes with their
            sequences and properties. This database should contain target probes organized by region IDs,
            with each region having one or more probe sets.
        :type target_probe_database: OligoDatabase
        :param codebook: A pandas DataFrame containing binary barcodes for each region. Rows are indexed
            by region IDs, and columns represent bit positions (bit_1, bit_2, etc.). Each row has exactly
            four bits set to 1, indicating which readout probes are assigned to that region.
        :type codebook: pd.DataFrame
        :param readout_probe_table: A pandas DataFrame containing readout probe sequences and their
            associated bit identifiers. The DataFrame should be indexed by bit labels (bit_1, bit_2, etc.)
            and contain a 'readout_probe_sequence' column with the probe sequences.
        :type readout_probe_table: pd.DataFrame
        :return: An updated `OligoDatabase` object containing the assembled hybridization probes. The
            database includes the following new sequence properties for each probe:
            - `sequence_target`: The gene-specific targeting sequence
            - `sequence_readout_probe_1`: The first readout probe sequence (from the barcode)
            - `sequence_readout_probe_2`: The second readout probe sequence (from the barcode)
            - `sequence_readout_probe_3`: The third readout probe sequence (from the barcode)
            - `sequence_readout_probe_4`: The fourth readout probe sequence (from the barcode)
            - `sequence_hybridization_probe`: The complete assembled hybridization probe sequence
        :rtype: OligoDatabase
        """
        region_ids = list(target_probe_database.database.keys())

        target_probe_database.set_database_sequence_types(
            [
                "sequence_target",
                "sequence_readout_probe_1",
                "sequence_readout_probe_2",
                "sequence_readout_probe_3",
                "sequence_readout_probe_4",
                "sequence_hybridization_probe",
            ]
        )

        for region_id in region_ids:
            barcode = codebook.loc[region_id]
            bits = barcode[barcode == 1].index
            readout_probe_sequences = readout_probe_table.loc[bits, "readout_probe_sequence"]
            sequence_readout_probe_1 = readout_probe_sequences.iloc[0]
            sequence_readout_probe_2 = readout_probe_sequences.iloc[1]
            sequence_readout_probe_3 = readout_probe_sequences.iloc[2]
            sequence_readout_probe_4 = readout_probe_sequences.iloc[3]

            probe_ids = list(target_probe_database.database[region_id].keys())
            new_properties: dict[str, dict[str, str]] = {probe_id: {} for probe_id in probe_ids}

            for probe_id in probe_ids:

                new_properties[probe_id]["sequence_target"] = format_sequence(
                    database=target_probe_database,
                    property="target",
                    region_id=region_id,
                    oligo_id=probe_id,
                )

                new_properties[probe_id]["sequence_readout_probe_1"] = sequence_readout_probe_1
                new_properties[probe_id]["sequence_readout_probe_2"] = sequence_readout_probe_2
                new_properties[probe_id]["sequence_readout_probe_3"] = sequence_readout_probe_3
                new_properties[probe_id]["sequence_readout_probe_4"] = sequence_readout_probe_4

                new_properties[probe_id]["sequence_hybridization_probe"] = (
                    str(Seq(sequence_readout_probe_1).reverse_complement())
                    + str(Seq(sequence_readout_probe_2).reverse_complement())
                    + "T"
                    + format_sequence(
                        database=target_probe_database,
                        property="oligo",
                        region_id=region_id,
                        oligo_id=probe_id,
                    )
                    + "T"
                    + str(Seq(sequence_readout_probe_3).reverse_complement())
                    + str(Seq(sequence_readout_probe_4).reverse_complement())
                )

            target_probe_database.update_oligo_properties(new_properties)

        return target_probe_database

    def design_primers(
        self,
        # Step 1: Create Database Parameters
        reverse_primer_sequence: str,
        primer_length: int,
        primer_base_probabilities: dict[str, float],
        primer_initial_num_sequences: int,
        # Step 2: Property Filter Parameters
        primer_GC_content_min: float,
        primer_GC_content_max: float,
        primer_number_GC_GCclamp: int,
        primer_number_three_prime_base_GCclamp: int,
        primer_homopolymeric_base_n: dict[str, int],
        primer_max_len_selfcomplement: int,
        primer_max_len_complement_reverse_primer: int,
        primer_Tm_min: float,
        primer_Tm_max: float,
        primer_T_secondary_structure: float,
        primer_secondary_structures_threshold_deltaG: float,
        primer_Tm_parameters: dict,
        primer_Tm_chem_correction_parameters: dict | None,
        primer_Tm_salt_correction_parameters: dict | None,
        # Step 3: Specificity Filter Parameters
        hybridization_probe_database: OligoDatabase,
        files_fasta_reference_database_primer: list[str],
        primer_specificity_refrence_blastn_search_parameters: dict,
        primer_specificity_refrence_blastn_hit_parameters: dict,
        primer_specificity_hybridization_probes_blastn_search_parameters: dict,
        primer_specificity_hybridization_probes_blastn_hit_parameters: dict,
    ) -> tuple[str, str]:
        """
        Design forward and reverse primers for SeqFISH+ hybridization probes through a multi-step pipeline.

        This method performs the complete primer design process, which includes:
        1. Creating an initial oligo database by generating random sequences with specified nucleotide
           probabilities and length (all ending with "T" nucleotide)
        2. Filtering primers based on sequence properties (GC content, GC clamp, homopolymeric runs,
           self-complementarity, complementarity to reverse primer, melting temperature, secondary structure)
        3. Filtering primers based on specificity to remove primers that bind to reference sequences
           or to the hybridization probes themselves using BLASTN searches
        4. Selecting the forward primer that has a melting temperature closest to the reverse primer's
           melting temperature to ensure balanced PCR amplification

        The reverse primer sequence is provided as input, and the method designs a forward primer that
        matches its melting temperature for optimal PCR performance.

        **Step 1: Create Database Parameters**

        :param reverse_primer_sequence: DNA sequence of the reverse primer that will be used as a reference.
            The forward primer will be selected to match this primer's melting temperature.
        :type reverse_primer_sequence: str
        :param primer_length: Length (in nucleotides) of each primer sequence to generate.
        :type primer_length: int
        :param primer_base_probabilities: Dictionary specifying the probability of each nucleotide base
            in randomly generated sequences. Keys should be 'A', 'T', 'G', 'C' and values should sum to 1.0
            (e.g., {"A": 0.25, "T": 0.25, "G": 0.25, "C": 0.25}).
        :type primer_base_probabilities: dict[str, float]
        :param primer_initial_num_sequences: Number of random sequences to generate initially before filtering.
            Higher values provide more candidates but increase computation time.
        :type primer_initial_num_sequences: int

        **Step 2: Property Filter Parameters**

        :param primer_GC_content_min: Minimum acceptable GC content for primers, expressed as a fraction
            between 0.0 and 1.0.
        :type primer_GC_content_min: float
        :param primer_GC_content_max: Maximum acceptable GC content for primers, expressed as a fraction
            between 0.0 and 1.0.
        :type primer_GC_content_max: float
        :param primer_number_GC_GCclamp: Minimum number of G or C nucleotides required within the specified
            number of bases at the 3' end (GC clamp). This improves primer binding stability.
        :type primer_number_GC_GCclamp: int
        :param primer_number_three_prime_base_GCclamp: Number of bases from the 3' end to consider for
            the GC clamp requirement.
        :type primer_number_three_prime_base_GCclamp: int
        :param primer_homopolymeric_base_n: Dictionary specifying the maximum allowed length of homopolymeric
            runs for each nucleotide base (keys: 'A', 'T', 'G', 'C').
        :type primer_homopolymeric_base_n: dict[str, int]
        :param primer_max_len_selfcomplement: Maximum allowable length of self-complementary sequences.
            Primers with longer self-complementary regions can form hairpins and reduce PCR efficiency.
        :type primer_max_len_selfcomplement: int
        :param primer_max_len_complement_reverse_primer: Maximum allowable length of complementarity to the
            reverse primer sequence. This prevents the forward and reverse primers from binding to each other.
        :type primer_max_len_complement_reverse_primer: int
        :param primer_Tm_min: Minimum acceptable melting temperature (Tm) for primers in degrees Celsius.
        :type primer_Tm_min: float
        :param primer_Tm_max: Maximum acceptable melting temperature (Tm) for primers in degrees Celsius.
        :type primer_Tm_max: float
        :param primer_T_secondary_structure: Temperature in degrees Celsius at which to evaluate secondary
            structure formation.
        :type primer_T_secondary_structure: float
        :param primer_secondary_structures_threshold_deltaG: DeltaG threshold (in kcal/mol) for secondary
            structure stability. Primers with secondary structures having deltaG values more negative than
            this threshold will be filtered out.
        :type primer_secondary_structures_threshold_deltaG: float
        :param primer_Tm_parameters: Dictionary of parameters for calculating melting temperature (Tm) of primers
            using the nearest-neighbor method. For using Bio.SeqUtils.MeltingTemp default parameters, set to ``{}``.
            Common parameters include: 'nn_table', 'tmm_table', 'imm_table', 'de_table', 'dnac1', 'dnac2', 'Na', 'K',
            'Tris', 'Mg', 'dNTPs', 'saltcorr', etc. For more information on parameters, see:
            https://biopython.org/docs/1.75/api/Bio.SeqUtils.MeltingTemp.html#Bio.SeqUtils.MeltingTemp.Tm_NN
        :type primer_Tm_parameters: dict
        :param primer_Tm_chem_correction_parameters: Dictionary of chemical correction parameters for Tm calculation.
            These parameters account for the effects of chemical additives (e.g., DMSO, formamide) on melting temperature.
            Set to ``None`` to disable chemical correction, or set to ``{}`` to use Bio.SeqUtils.MeltingTemp default parameters.
            For more information, see:
            https://biopython.org/docs/1.75/api/Bio.SeqUtils.MeltingTemp.html#Bio.SeqUtils.MeltingTemp.chem_correction
        :type primer_Tm_chem_correction_parameters: dict | None
        :param primer_Tm_salt_correction_parameters: Dictionary of salt correction parameters for Tm calculation.
            These parameters account for the effects of salt concentration on melting temperature. Set to ``None`` to disable
            salt correction, or set to ``{}`` to use Bio.SeqUtils.MeltingTemp default parameters. For more information, see:
            https://biopython.org/docs/1.75/api/Bio.SeqUtils.MeltingTemp.html#Bio.SeqUtils.MeltingTemp.salt_correction
        :type primer_Tm_salt_correction_parameters: dict | None

        **Step 3: Specificity Filter Parameters**

        :param hybridization_probe_database: The `OligoDatabase` instance containing hybridization probes.
            This database is used to create a reference FASTA file for specificity filtering to ensure primers
            do not bind to the hybridization probes themselves.
        :type hybridization_probe_database: OligoDatabase
        :param files_fasta_reference_database_primer: List of paths to FASTA files containing reference sequences
            used for specificity filtering. These files are used to identify off-target binding sites (e.g.,
            whole genome or transcriptome sequences).
        :type files_fasta_reference_database_primer: list[str]
        :param primer_specificity_refrence_blastn_search_parameters: Dictionary of parameters for BLASTN
            searches used in specificity filtering against the reference database.
        :type primer_specificity_refrence_blastn_search_parameters: dict
        :param primer_specificity_refrence_blastn_hit_parameters: Dictionary of parameters for filtering
            BLASTN hits in specificity searches against the reference database.
            Primers with hits meeting these criteria are removed.
        :type primer_specificity_refrence_blastn_hit_parameters: dict
        :param primer_specificity_hybridization_probes_blastn_search_parameters: Dictionary of parameters for BLASTN
            searches used in specificity filtering against the hybridization probes database.
        :type primer_specificity_hybridization_probes_blastn_search_parameters: dict
        :param primer_specificity_hybridization_probes_blastn_hit_parameters: Dictionary of parameters for filtering
            BLASTN hits in specificity searches against the hybridization probes database. Primers with hits meeting these
            criteria are removed.
        :type primer_specificity_hybridization_probes_blastn_hit_parameters: dict
        :return: A tuple containing:
            - **reverse_primer_sequence** (str): The input reverse primer sequence (unchanged)
            - **forward_primer_sequence** (str): The selected forward primer sequence with melting temperature
              closest to the reverse primer's melting temperature
        :rtype: tuple[str, str]
        """
        file_fasta_hybridization_probes_database = hybridization_probe_database.write_database_to_fasta(
            filename=f"db_reference_hybridization_probes",
            save_description=False,
            region_ids=None,
            sequence_type="sequence_hybridization_probe",
        )

        primer_designer = PrimerDesigner(
            dir_output=self.dir_output,
            n_jobs=self.n_jobs,
        )
        oligo_database: OligoDatabase = primer_designer.create_oligo_database(
            oligo_length=primer_length,
            oligo_base_probabilities=primer_base_probabilities,
            initial_num_sequences=primer_initial_num_sequences,
        )

        if self.write_intermediate_steps:
            dir_database = oligo_database.save_database(name_database="1_db_primers_initial")
            logging.info(f"Saved primer database for step 1 (Create Database) in directory {dir_database}")

        oligo_database = primer_designer.filter_by_property(
            oligo_database=oligo_database,
            GC_content_min=primer_GC_content_min,
            GC_content_max=primer_GC_content_max,
            number_GC_GCclamp=primer_number_GC_GCclamp,
            number_three_prime_base_GCclamp=primer_number_three_prime_base_GCclamp,
            homopolymeric_base_n=primer_homopolymeric_base_n,
            max_len_selfcomplement=primer_max_len_selfcomplement,
            reverse_primer_sequence=reverse_primer_sequence,
            max_len_complement=primer_max_len_complement_reverse_primer,
            Tm_min=primer_Tm_min,
            Tm_max=primer_Tm_max,
            Tm_parameters=primer_Tm_parameters,
            Tm_chem_correction_parameters=primer_Tm_chem_correction_parameters,
            Tm_salt_correction_parameters=primer_Tm_salt_correction_parameters,
            T_secondary_structure=primer_T_secondary_structure,
            secondary_structures_threshold_deltaG=primer_secondary_structures_threshold_deltaG,
        )

        if self.write_intermediate_steps:
            dir_database = oligo_database.save_database(name_database="2_db_primer_property_filter")
            logging.info(f"Saved primer database for step 2 (Property Filters) in directory {dir_database}")

        oligo_database = primer_designer.filter_by_specificity(
            oligo_database=oligo_database,
            files_fasta_reference_database=files_fasta_reference_database_primer,
            specificity_refrence_blastn_search_parameters=primer_specificity_refrence_blastn_search_parameters,
            specificity_refrence_blastn_hit_parameters=primer_specificity_refrence_blastn_hit_parameters,
            file_fasta_hybridization_probes_database=file_fasta_hybridization_probes_database,
            specificity_hybridization_probes_blastn_search_parameters=primer_specificity_hybridization_probes_blastn_search_parameters,
            specificity_hybridization_probes_blastn_hit_parameters=primer_specificity_hybridization_probes_blastn_hit_parameters,
        )

        if self.write_intermediate_steps:
            dir_database = oligo_database.save_database(name_database="3_db_primer_specificty_filter")
            logging.info(
                f"Saved primer database for step 3 (Specificity Filters) in directory {dir_database}"
            )

        # calculate Tm for the reverse primer
        Tm_reverse_primer = calc_tm_nn(
            sequence=reverse_primer_sequence,
            Tm_parameters=primer_Tm_parameters,
            Tm_chem_correction_parameters=primer_Tm_chem_correction_parameters,
            Tm_salt_correction_parameters=primer_Tm_salt_correction_parameters,
        )

        # iterate over all primers in the database to find the one with Tm closest to the reverse primer Tm
        min_dif_Tm = float("inf")
        forward_primer_sequence = ""
        for database_region in oligo_database.database.values():
            for primer_properties in database_region.values():
                Tm_forward_primer = calc_tm_nn(
                    sequence=primer_properties["oligo"],
                    Tm_parameters=primer_Tm_parameters,
                    Tm_chem_correction_parameters=primer_Tm_chem_correction_parameters,
                    Tm_salt_correction_parameters=primer_Tm_salt_correction_parameters,
                )
                dif_Tm = abs(Tm_forward_primer - Tm_reverse_primer)
                if dif_Tm < min_dif_Tm:
                    min_dif_Tm = dif_Tm
                    forward_primer_sequence = primer_properties["oligo"]

        os.remove(file_fasta_hybridization_probes_database)

        return reverse_primer_sequence, forward_primer_sequence

    def assemble_dna_template_probes(
        self,
        hybridization_probe_database: OligoDatabase,
        forward_primer_sequence: str,
        reverse_primer_sequence: str,
    ) -> OligoDatabase:
        """
        Assemble DNA template probes by combining hybridization probes with forward and reverse primers.

        This method creates the final DNA template probes used for PCR amplification by combining
        hybridization probes with PCR primer sequences. For each probe in the database, the method
        assembles the DNA template probe with the structure:
        [Forward Primer] + [Hybridization Probe] + [Reverse Primer]

        The assembled sequences are stored as properties in the database for each probe, ready for
        synthesis and experimental use.

        :param hybridization_probe_database: The `OligoDatabase` instance containing hybridization probes
            with their sequences and properties. This database should contain the `sequence_hybridization_probe`
            property for each probe, which was created by the `assemble_hybridization_probes` method.
        :type hybridization_probe_database: OligoDatabase
        :param forward_primer_sequence: DNA sequence of the forward primer that binds to the 5' end of
            the DNA template probe.
        :type forward_primer_sequence: str
        :param reverse_primer_sequence: DNA sequence of the reverse primer that binds to the 3' end of
            the DNA template probe.
        :type reverse_primer_sequence: str
        :return: An updated `OligoDatabase` object containing the assembled DNA template probes. The
            database includes the following new sequence properties for each probe:
            - `sequence_forward_primer`: The forward primer sequence
            - `sequence_reverse_primer`: The reverse primer sequence
            - `sequence_dna_template_probe`: The complete assembled DNA template probe sequence
        :rtype: OligoDatabase
        """
        region_ids = list(hybridization_probe_database.database.keys())
        hybridization_probe_database.set_database_sequence_types(
            [
                "sequence_forward_primer",
                "sequence_reverse_primer",
                "sequence_dna_template_probe",
            ]
        )

        for region_id in region_ids:
            probe_ids = list(hybridization_probe_database.database[region_id].keys())
            new_properties: dict[str, dict[str, str]] = {probe_id: {} for probe_id in probe_ids}

            for probe_id in probe_ids:
                new_properties[probe_id]["sequence_reverse_primer"] = reverse_primer_sequence
                new_properties[probe_id]["sequence_forward_primer"] = forward_primer_sequence

                new_properties[probe_id]["sequence_dna_template_probe"] = (
                    forward_primer_sequence
                    + format_sequence(
                        database=hybridization_probe_database,
                        property="sequence_hybridization_probe",
                        region_id=region_id,
                        oligo_id=probe_id,
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
        Generate the final output files for the SeqFISH+ probe design pipeline.

        This method writes all output files required for the SeqFISH+ experiment, including codebooks,
        readout probe tables, and probe sequences in multiple formats. The output files are written
        to the pipeline's output directory.

        **Generated Output Files:**

        1. **codebook.tsv**: Binary barcode matrix with region IDs as index and bit columns. Each row
           represents a region's barcode assignment using the pseudocolor and barcode round system.

        2. **readout_probes.tsv**: Table mapping readout probe sequences to bit positions, barcode rounds,
           pseudocolors, and channels. Contains columns: bit, barcode_round, pseudocolor, channel,
           readout_probe_sequence.

        3. **seqfish_plus_probes.yml**: Complete probe information in YAML format, including all specified
           properties for each probe set per region.

        4. **seqfish_plus_probes.tsv**: Complete probe information in TSV format, including all specified
           properties for each probe set per region.

        5. **seqfish_plus_probes.xlsx**: Complete probe information in Excel format with one sheet per region.
           Each sheet contains probe sets for that region with all specified properties.

        6. **seqfish_plus_probes_order.yml**: Simplified YAML file containing only the essential sequences
           needed for ordering probes (DNA template probe and all four readout probe sequences).

        :param probe_database: The `OligoDatabase` instance containing the final DNA template probes
            with all sequences and properties. This should be the result of the `assemble_dna_template_probes`
            method.
        :type probe_database: OligoDatabase
        :param codebook: A pandas DataFrame containing binary barcodes for each region. Rows are indexed
            by region IDs, and columns represent bit positions. This should be the codebook generated
            by the `design_readout_probes` method.
        :type codebook: pd.DataFrame
        :param readout_probe_table: A pandas DataFrame containing readout probe sequences and their
            associated bit identifiers, barcode rounds, pseudocolors, and channel assignments. This should
            be the readout probe table generated by the `design_readout_probes` method.
        :type readout_probe_table: pd.DataFrame
        :param output_properties: List of property names to include in the output files. If None, a default
            set of properties will be included. Available properties include: 'source', 'species', 'gene_id',
            'chromosome', 'start', 'end', 'strand', 'sequence_target', 'sequence_readout_probe_1',
            'sequence_readout_probe_2', 'sequence_readout_probe_3', 'sequence_readout_probe_4',
            'sequence_hybridization_probe', 'sequence_forward_primer', 'sequence_reverse_primer',
            'sequence_dna_template_probe', 'isoform_consensus', etc.
        :type output_properties: list[str] | None
        :param properties: Default list of properties to include if `output_properties` is None. This parameter
            is used internally and typically should not be modified.
        :type properties: list[str]
        :return: None. All output files are written to the pipeline's output directory.
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
                "sequence_readout_probe_1",
                "sequence_readout_probe_2",
                "sequence_readout_probe_3",
                "sequence_readout_probe_4",
                "sequence_hybridization_probe",
                "sequence_forward_primer",
                "sequence_reverse_primer",
                "sequence_dna_template_probe",
                "isoform_consensus",
            ]

        codebook.to_csv(os.path.join(self.dir_output, "codebook.tsv"), sep="\t", index_label="region_id")
        readout_probe_table.to_csv(os.path.join(self.dir_output, "readout_probes.tsv"), sep="\t")

        probe_database.write_oligosets_to_yaml(
            properties=output_properties,
            ascending=True,
            filename="seqfish_plus_probes",
        )

        probe_database.write_oligosets_to_table(
            properties=output_properties,
            ascending=True,
            filename="seqfish_plus_probes",
        )

        probe_database.write_ready_to_order_yaml(
            properties=[
                "sequence_dna_template_probe",
                "sequence_readout_probe_1",
                "sequence_readout_probe_2",
                "sequence_readout_probe_3",
                "sequence_readout_probe_4",
            ],
            ascending=True,
            filename="seqfish_plus_probes_order",
        )


############################################
# SeqFish Plus Target Probe Designer
############################################


class TargetProbeDesigner:
    """
    A class for designing target probes for SeqFISH+ experiments through a multi-step pipeline.

    This class provides methods for the complete target probe design process, which includes:
    1. Creating an initial oligo database from input FASTA files using a sliding window approach
    2. Filtering probes based on sequence properties (GC content, homopolymeric runs, secondary structure)
    3. Filtering probes based on specificity to remove off-target binding and cross-hybridization
       using BLASTN searches
    4. Organizing filtered probes into optimal sets based on weighted scoring criteria (GC content,
       UTR targeting) and distance constraints

    The resulting probes are gene-specific targeting sequences (typically 28 nt) that bind to RNA
    transcripts. These probes will later be combined with four readout probe overhang sequences to
    create complete hybridization probes.

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
        self.subdir_db_oligos = "db_probes"
        self.subdir_db_reference = "db_reference"

        self.n_jobs = n_jobs

    @pipeline_step_basic(step_name="Target Probe Generation - Create Database")
    def create_oligo_database(
        self,
        region_ids: list | None,
        oligo_length_min: int,
        oligo_length_max: int,
        files_fasta_oligo_database: list[str],
        min_oligos_per_gene: int,
        isoform_consensus: float,
    ) -> OligoDatabase:
        """
        Create an initial oligo database by generating sequences using a sliding window approach
        and performing pre-filtering based on isoform consensus.

        This is the first step of target probe design. The method:
        1. Generates candidate oligo sequences from input FASTA files using a sliding window approach
           across the specified length range
        2. Creates an `OligoDatabase` and loads the generated sequences
        3. Calculates reverse complement sequences and isoform consensus properties
        4. Pre-filters oligos based on isoform consensus threshold
        5. Removes regions with insufficient oligos after filtering

        The database stores sequences with sequence types "target" (original sequence) and
        "oligo" (reverse complement).

        :param region_ids: List of gene identifiers (e.g., gene IDs) to target for probe design. If None,
            all genes present in the input FASTA files will be used.
        :type region_ids: list | None
        :param oligo_length_min: Minimum length (in nucleotides) for target probe sequences.
        :type oligo_length_min: int
        :param oligo_length_max: Maximum length (in nucleotides) for target probe sequences.
        :type oligo_length_max: int
        :param files_fasta_oligo_database: List of paths to FASTA files containing sequences from which
            target probes will be generated. These files should contain genomic regions of interest
            (e.g., exons, exon-exon junctions).
        :type files_fasta_oligo_database: list[str]
        :param min_oligos_per_gene: Minimum number of oligos required per region (gene) after filtering.
            Regions with fewer oligos than this threshold will be removed from the database.
        :type min_oligos_per_gene: int
        :param isoform_consensus: Threshold for isoform consensus filtering (typically between 0.0 and 1.0).
            Probes with isoform consensus values below this threshold will be filtered out. This ensures
            that selected probes target sequences that are conserved across multiple transcript isoforms.
        :type isoform_consensus: float
        :return: An `OligoDatabase` object containing the generated target probe sequences with their
            component sequences (target, oligo) and calculated properties (isoform_consensus). The database
            is filtered to only include regions that meet the minimum oligo requirement.
        :rtype: OligoDatabase
        """
        ##### creating the oligo sequences #####
        oligo_sequences = OligoSequenceGenerator(dir_output=self.dir_output)
        oligo_fasta_file = oligo_sequences.create_sequences_sliding_window(
            files_fasta_in=files_fasta_oligo_database,
            length_interval_sequences=(oligo_length_min, oligo_length_max),
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
        oligo_database.set_database_sequence_types(["target", "oligo"])
        # Calculate reverse complement and isoform consensus
        properties = [
            ReverseComplementSequenceProperty(sequence_type_reverse_complement="oligo"),
            IsoformConsensusProperty(),
        ]
        calculator = PropertyCalculator(properties=properties)
        oligo_database = calculator.apply(
            oligo_database=oligo_database, sequence_type="target", n_jobs=self.n_jobs
        )

        ##### pre-filter oligo database for certain properties #####
        oligo_database.filter_database_by_property_threshold(
            property_name="isoform_consensus",
            property_thr=isoform_consensus,
            remove_if_smaller_threshold=True,
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
        homopolymeric_base_n: dict[str, int],
        T_secondary_structure: float,
        secondary_structures_threshold_deltaG: float,
    ) -> OligoDatabase:
        """
        Filter the oligo database based on sequence properties to remove probes with undesirable
        characteristics.

        This method applies sequential filtering using multiple property-based filters:
        1. **Hard masked sequences**: Removes probes containing hard-masked nucleotides (N)
        2. **Soft masked sequences**: Removes probes containing soft-masked nucleotides (lowercase)
        3. **Homopolymeric runs**: Removes probes with homopolymeric runs exceeding specified lengths
        4. **GC content**: Removes probes with GC content outside the specified range
        5. **Secondary structure**: Removes probes that form stable secondary structures at the
           specified temperature

        Probes that fail any filter are removed. Regions with insufficient oligos after filtering
        are removed from the database.

        :param oligo_database: The `OligoDatabase` instance containing oligonucleotide sequences
            and their associated properties. This database should contain target probes with their
            component sequences already calculated.
        :type oligo_database: OligoDatabase
        :param GC_content_min: Minimum acceptable GC content for oligos, expressed as a fraction
            between 0.0 and 1.0.
        :type GC_content_min: float
        :param GC_content_max: Maximum acceptable GC content for oligos, expressed as a fraction
            between 0.0 and 1.0.
        :type GC_content_max: float
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
        :return: A filtered `OligoDatabase` object containing only probes that pass all property filters.
            Regions with insufficient oligos after filtering are removed.
        :rtype: OligoDatabase
        """
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
            oligo_database=oligo_database,
            sequence_type="oligo",
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
           - **Exact matches**: Removes all probes with exact sequence matches to probes of other regions
           - **BLASTN specificity**: Uses BLASTN to search for similar sequences in the reference database.
             Probes with hits meeting the specified criteria are removed

        2. **Cross-hybridization filtering**: Removes probes that cross-hybridize with each other.
           This is critical because if probes can bind to each other, they may form dimers instead
           of binding to the target RNA. Probes from the larger genomic region are removed when
           cross-hybridization is detected.

        The reference database is loaded from the provided FASTA files and used for all BLASTN searches.
        Regions that do not meet the minimum oligo requirement after filtering are removed from
        the database.

        :param oligo_database: The `OligoDatabase` instance containing oligonucleotide sequences
            and their associated properties. This database should contain target probes with their
            component sequences already calculated.
        :type oligo_database: OligoDatabase
        :param files_fasta_reference_database: List of paths to FASTA files containing reference
            sequences against which specificity will be evaluated. These typically include the
            entire genome or transcriptome to identify off-target binding sites.
        :type files_fasta_reference_database: list[str]
        :param specificity_blastn_search_parameters: Dictionary of parameters for BLASTN searches
            used in specificity filtering.
        :type specificity_blastn_search_parameters: dict
        :param specificity_blastn_hit_parameters: Dictionary of parameters for filtering BLASTN hits
            in specificity searches. Probes with hits meeting these criteria are removed.
        :type specificity_blastn_hit_parameters: dict
        :param cross_hybridization_blastn_search_parameters: Dictionary of parameters for BLASTN
            searches used in cross-hybridization filtering. These searches check if probes align to
            each other.
        :type cross_hybridization_blastn_search_parameters: dict
        :param cross_hybridization_blastn_hit_parameters: Dictionary of parameters for filtering
            BLASTN hits in cross-hybridization searches. Probes with cross-hybridization hits meeting these
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

        ##### exact match filter #####
        exact_matches = ExactMatchFilter(policy=RemoveAllFilterPolicy(), filter_name="oligo_exact_match")

        ##### specificity filters #####
        specificity = BlastNFilter(
            remove_hits=True,
            search_parameters=specificity_blastn_search_parameters,
            hit_parameters=specificity_blastn_hit_parameters,
            filter_name="oligo_blastn_specificity",
            dir_output=self.dir_output,
        )
        specificity.set_reference_database(reference_database=reference_database)

        cross_hybridization_aligner = BlastNFilter(
            remove_hits=True,
            search_parameters=cross_hybridization_blastn_search_parameters,
            hit_parameters=cross_hybridization_blastn_hit_parameters,
            filter_name="oligo_blastn_crosshybridization",
            dir_output=self.dir_output,
        )
        cross_hybridization_aligner.set_reference_database(reference_database=reference_database)
        cross_hybridization = CrossHybridizationFilter(
            policy=RemoveByLargerRegionFilterPolicy(),
            alignment_method=cross_hybridization_aligner,
            filter_name="oligo_blastn_crosshybridization",
            dir_output=self.dir_output,
        )

        specificity_filter = SpecificityFilter(filters=[exact_matches, specificity, cross_hybridization])
        oligo_database = specificity_filter.apply(
            oligo_database=oligo_database,
            sequence_type="oligo",
            n_jobs=self.n_jobs,
        )

        # remove all directories of intermediate steps
        for directory in [
            cross_hybridization_aligner.dir_output,
            cross_hybridization.dir_output,
            specificity.dir_output,
        ]:
            if os.path.exists(directory):
                shutil.rmtree(directory)

        oligo_database.remove_regions_with_insufficient_oligos(pipeline_step="Specificity Filters")
        check_content_oligo_database(oligo_database)

        return oligo_database

    @pipeline_step_basic(step_name="Target Probe Generation - Set Selection")
    def create_oligo_sets(
        self,
        oligo_database: OligoDatabase,
        GC_weight: float,
        GC_content_opt: float,
        UTR_weight: float,
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

        This method performs the following steps:
        1. **Scoring**: Calculates scores for each oligo based on weighted criteria (GC content, UTR targeting).
           Higher scores indicate better probes.
        2. **Set generation**: Builds a compatibility graph from distance constraints and selects sets via
           a graph-based (clique) strategy. Generates multiple diverse sets per region, controlling overlap
           between sets using a Jaccard threshold (`jaccard_opt`) with optional relaxation (`jaccard_step`).
        3. **Set scoring**: Evaluates each generated set and selects the best sets based on the lowest
           average score (ascending order).
        4. **Region filtering**: Removes regions that cannot generate sets meeting the minimum size requirement.

        The algorithm attempts to find sets with optimal size (`set_size_opt`) but may produce sets
        as small as `set_size_min` if constraints cannot be met.

        :param oligo_database: The `OligoDatabase` instance containing oligonucleotide sequences
            and their associated properties. This database should contain target probes that have
            passed all previous filtering steps.
        :type oligo_database: OligoDatabase
        :param GC_weight: Weight assigned to GC content in the scoring function.
        :type GC_weight: float
        :param GC_content_opt: Optimal GC content for oligos, expressed as a fraction between 0.0
            and 1.0. Used in scoring to prioritize probes closer to this value.
        :type GC_content_opt: float
        :param UTR_weight: Weight assigned to untranslated region (UTR) targeting in the scoring function.
            Probes that target UTRs are prioritized.
        :type UTR_weight: float
        :param set_size_opt: Optimal size (number of probes) for each oligo set. The set selection algorithm will
            attempt to generate sets of this size, but may produce sets with fewer probes if constraints cannot be met.
        :type set_size_opt: int
        :param set_size_min: Minimum size (number of probes) required for each oligo set. Sets with fewer probes than
            this value will be rejected, and regions that cannot generate sets meeting this minimum will be removed.
        :type set_size_min: int
        :param distance_between_oligos: Minimum genomic distance (in nucleotides) required between probes
            within the same set. This spacing constraint prevents probes from binding too close together, which could
            lead to reduced hybridization efficiency.
        :type distance_between_oligos: int
        :param n_sets: Number of oligo sets to generate per region. Multiple sets allow for redundancy and selection
            of the best-performing set based on scoring criteria.
        :type n_sets: int
        :param n_attempts_graph: Number of randomized graph attempts. In each attempt, a fraction of nodes is randomly
            removed from the compatibility graph to create diversity.
        :type n_attempts_graph: int
        :param n_attempts_clique_enum: Maximum number of cliques enumerated per graph attempt.
        :type n_attempts_clique_enum: int
        :param diversification_fraction: Fraction of oligos to remove at random per attempt to create diversity
            between sets.
        :type diversification_fraction: float
        :param jaccard_opt: Optimal maximum Jaccard overlap between selected sets. Sets with overlap above this
            value are discouraged when selecting multiple sets per region.
        :type jaccard_opt: float
        :param jaccard_step: Step size for relaxing Jaccard overlap when not enough sets are found.
        :type jaccard_step: float
        :return: An updated `OligoDatabase` object containing the generated oligo sets. Each region
            will have up to `n_sets` sets stored, with each set containing between `set_size_min` and
            `set_size_opt` probes. Regions with insufficient oligos are removed.
        :rtype: OligoDatabase
        """
        # Define all scorers
        utr_scorer = OverlapUTRScorer(score_weight=UTR_weight, property_name="regiontype")
        GC_scorer = DeviationFromOptimalGCContentScorer(
            GC_content_opt=GC_content_opt,
            score_weight=GC_weight,
        )
        oligos_scoring = OligoScoring(scorers=[utr_scorer, GC_scorer])
        set_scoring = LowestSetScoring(ascending=True)

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
# SeqFish Plus Readout Probe Designer
############################################


class ReadoutProbeDesigner:
    """
    A class for designing SeqFISH+ readout probes and generating codebooks through a multi-step pipeline.

    This class provides methods for the complete readout probe design process, which includes:
    1. Creating an initial oligo database by generating random sequences with specified nucleotide
       probabilities
    2. Filtering probes based on sequence properties (GC content, homopolymeric runs)
    3. Filtering probes based on specificity to remove off-target binding and cross-hybridization
       using BLASTN searches
    4. Generating a combinatorial pseudocolor codebook with barcode rounds for hybridization regions
    5. Creating a readout probe table that maps codebook bits to barcode rounds, pseudocolors,
       channels, and readout probe sequences

    The resulting readout probes are non-targeting sequences (typically 15 nt overhangs) that bind
    to barcode sequences on the hybridization probes. Each readout probe is assigned to a specific bit
    position in the codebook and distributed across fluorescence channels and barcode rounds for
    multiplexed detection using a combinatorial pseudocolor system.

    :param dir_output: Directory path where output files will be saved.
    :type dir_output: str
    :param n_jobs: Number of parallel jobs to use for processing. Set to 1 for serial processing or higher
        values for parallel processing. This affects the parallelization of filtering, property calculation,
        and set generation operations.
    :type n_jobs: int
    """

    def __init__(
        self,
        dir_output: str,
        n_jobs: int,
    ) -> None:
        """Constructor for the ReadoutProbeDesigner class."""

        ##### create the output folder #####
        self.dir_output = os.path.abspath(dir_output)
        self.subdir_db_oligos = "db_readout_probes"
        self.subdir_db_reference = "db_reference"

        self.n_jobs = n_jobs

    @pipeline_step_basic(step_name="Readout Probe Generation - Create Oligo Database")
    def create_oligo_database(
        self,
        oligo_length: int,
        oligo_base_probabilities: dict,
        initial_num_sequences: int,
    ) -> OligoDatabase:
        """
        Create an initial oligo database by generating random sequences with specified nucleotide probabilities.

        This is the first step of readout probe design. The method generates random DNA sequences
        with user-defined nucleotide probabilities and creates an `OligoDatabase` to store them.
        These sequences will be filtered in subsequent steps.

        :param oligo_length: Length (in nucleotides) of each readout probe sequence to generate.
        :type oligo_length: int
        :param oligo_base_probabilities: Dictionary specifying the probability of each nucleotide
            base in randomly generated sequences. Keys should be 'A', 'T', 'G', 'C' and values should
            sum to 1.0 (e.g., {"A": 0.25, "T": 0.25, "G": 0.25, "C": 0.25}).
        :type oligo_base_probabilities: dict
        :param initial_num_sequences: Number of random sequences to generate initially before filtering.
            Higher values provide more candidates but increase computation time.
        :type initial_num_sequences: int
        :return: An `OligoDatabase` object containing the generated random readout probe sequences.
            The database stores sequences with sequence type "oligo".
        :rtype: OligoDatabase
        """
        ##### creating the oligo sequences #####
        oligo_sequences = OligoSequenceGenerator(dir_output=self.dir_output)
        oligo_fasta_file = oligo_sequences.create_sequences_random(
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
            max_entries_in_memory=self.n_jobs * 2 + 2,
            database_name=self.subdir_db_oligos,
            dir_output=self.dir_output,
            n_jobs=1,
        )
        oligo_database.load_database_from_fasta(
            files_fasta=oligo_fasta_file,
            database_overwrite=True,
            sequence_type="oligo",
            region_ids=None,
        )
        # Set all sequence types that will be used in this pipeline
        oligo_database.set_database_sequence_types(["oligo"])

        dir = oligo_sequences.dir_output
        shutil.rmtree(dir) if os.path.exists(dir) else None

        return oligo_database

    @pipeline_step_basic(step_name="Readout Probe Generation - Property Filters")
    def filter_by_property(
        self,
        oligo_database: OligoDatabase,
        GC_content_min: float,
        GC_content_max: float,
        homopolymeric_base_n: dict[str, int],
    ) -> OligoDatabase:
        """
        Filter the oligo database based on sequence properties to remove probes with undesirable
        characteristics.

        This method applies sequential filtering using property-based filters:
        1. **GC content**: Removes probes with GC content outside the specified range
        2. **Homopolymeric runs**: Removes probes with homopolymeric runs exceeding specified lengths

        Probes that fail any filter are removed. Regions with insufficient oligos after filtering
        are removed from the database.

        :param oligo_database: The `OligoDatabase` instance containing oligonucleotide sequences
            and their associated properties. This database should contain readout probe sequences
            generated in the previous step.
        :type oligo_database: OligoDatabase
        :param GC_content_min: Minimum acceptable GC content for readout probes, expressed as a fraction
            between 0.0 and 1.0.
        :type GC_content_min: float
        :param GC_content_max: Maximum acceptable GC content for readout probes, expressed as a fraction
            between 0.0 and 1.0.
        :type GC_content_max: float
        :param homopolymeric_base_n: Dictionary specifying the maximum allowed length of homopolymeric
            runs for each nucleotide base. Keys should be 'A', 'T', 'G', 'C' and values are the maximum
            run length. For example: {'A': 3, 'T': 3, 'G': 3, 'C': 3} allows up to 3 consecutive
            identical bases.
        :type homopolymeric_base_n: dict[str, int]
        :return: A filtered `OligoDatabase` object containing only probes that pass all property filters.
            Regions with insufficient oligos after filtering are removed.
        :rtype: OligoDatabase
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
            oligo_database=oligo_database,
            sequence_type="oligo",
            n_jobs=self.n_jobs,
        )
        oligo_database.remove_regions_with_insufficient_oligos(pipeline_step="Property Filters")
        check_content_oligo_database(oligo_database)

        return oligo_database

    @pipeline_step_basic(step_name="Readout Probe Generation - Specificity Filters")
    def filter_by_specificity(
        self,
        oligo_database: OligoDatabase,
        files_fasta_reference_database: list,
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
           - **Exact matches**: Removes all probes with exact sequence matches to other probes
           - **BLASTN specificity**: Uses BLASTN to search for similar sequences in the reference database.
             Probes with hits meeting the specified criteria are removed

        2. **Cross-hybridization filtering**: Removes probes that cross-hybridize with each other.
           This is critical because if probes can bind to each other, they may form dimers instead
           of binding to their intended targets. Probes are removed based on their degree of
           cross-hybridization (using `RemoveByDegreeFilterPolicy`).

        The reference database is loaded from the provided FASTA files and used for all BLASTN searches.
        Regions that do not meet the minimum oligo requirement after filtering are removed from
        the database.

        :param oligo_database: The `OligoDatabase` instance containing oligonucleotide sequences
            and their associated properties. This database should contain readout probe sequences
            that have passed property filtering.
        :type oligo_database: OligoDatabase
        :param files_fasta_reference_database: List of paths to FASTA files containing reference
            sequences against which specificity will be evaluated. These typically include the
            entire genome or transcriptome to identify off-target binding sites.
        :type files_fasta_reference_database: list
        :param specificity_blastn_search_parameters: Dictionary of parameters for BLASTN searches
            used in specificity filtering. Common parameters include: 'task', 'word_size', 'evalue',
            'max_target_seqs', 'num_threads', etc.
        :type specificity_blastn_search_parameters: dict
        :param specificity_blastn_hit_parameters: Dictionary of parameters for filtering BLASTN hits
            in specificity searches. Common parameters include: 'identity_min', 'alignment_length_min',
            'mismatches_max', 'gaps_max', etc. Probes with hits meeting these criteria are removed.
        :type specificity_blastn_hit_parameters: dict
        :param cross_hybridization_blastn_search_parameters: Dictionary of parameters for BLASTN
            searches used in cross-hybridization filtering. These searches check if probes align to
            each other. Common parameters are similar to `specificity_blastn_search_parameters`.
        :type cross_hybridization_blastn_search_parameters: dict
        :param cross_hybridization_blastn_hit_parameters: Dictionary of parameters for filtering
            BLASTN hits in cross-hybridization searches. Common parameters are similar to
            `specificity_blastn_hit_parameters`. Probes with cross-hybridization hits meeting these
            criteria are removed based on their degree of cross-hybridization.
        :type cross_hybridization_blastn_hit_parameters: dict
        :return: A filtered `OligoDatabase` object containing only probes that pass all specificity
            and cross-hybridization filters. Regions with insufficient oligos after filtering are removed.
        :rtype: OligoDatabase
        """
        reference_database = ReferenceDatabase(
            database_name=self.subdir_db_reference, dir_output=self.dir_output
        )
        reference_database.load_database_from_file(
            files=files_fasta_reference_database, file_type="fasta", database_overwrite=False
        )

        ##### specificity filters #####
        # removing duplicated oligos
        exact_matches = ExactMatchFilter(
            policy=RemoveAllFilterPolicy(), filter_name="readout_probes_exact_match"
        )

        # BlastN Filter
        specificity = BlastNFilter(
            search_parameters=specificity_blastn_search_parameters,
            hit_parameters=specificity_blastn_hit_parameters,
            filter_name="readout_probes_blastn_specificity",
            dir_output=self.dir_output,
        )
        specificity.set_reference_database(reference_database=reference_database)

        # Cross-Hybridization Filter
        cross_hybridization_aligner = BlastNFilter(
            search_parameters=cross_hybridization_blastn_search_parameters,
            hit_parameters=cross_hybridization_blastn_hit_parameters,
            filter_name="readout_probes_blastn_crosshybridization",
            dir_output=self.dir_output,
        )
        cross_hybridization_aligner.set_reference_database(reference_database=reference_database)

        cross_hybridization = CrossHybridizationFilter(
            policy=RemoveByDegreeFilterPolicy(),
            alignment_method=cross_hybridization_aligner,
            filter_name="readout_probes_blastn_crosshybridization",
            dir_output=self.dir_output,
        )

        specificity_filter = SpecificityFilter(filters=[exact_matches, specificity, cross_hybridization])
        oligo_database = specificity_filter.apply(
            oligo_database=oligo_database,
            sequence_type="oligo",
            n_jobs=self.n_jobs,
        )

        # remove all directories of intermediate steps
        for directory in [
            specificity.dir_output,
            cross_hybridization_aligner.dir_output,
            cross_hybridization.dir_output,
        ]:
            if os.path.exists(directory):
                shutil.rmtree(directory)
        oligo_database.remove_regions_with_insufficient_oligos(pipeline_step="Specificity Filters")
        check_content_oligo_database(oligo_database)
        return oligo_database

    def generate_codebook(
        self, region_ids: list[str], n_barcode_rounds: int, n_pseudocolors: int, n_channels: int
    ) -> pd.DataFrame:
        """
        Generate a combinatorial pseudocolor codebook with barcode rounds for hybridization regions.

        This method generates a codebook that assigns each region a unique barcode using a combinatorial
        pseudocolor system. The codebook ensures that:
        - Each barcode is constructed from combinations of pseudocolors across multiple barcode rounds
        - Each barcode includes a channel assignment for fluorescence detection
        - The number of valid barcodes is sufficient to encode all regions

        The combinatorial pseudocolor system allows for efficient encoding: with `n_pseudocolors` pseudocolors
        and `n_barcode_rounds` rounds, the codebook can encode up to `n_channels * (n_pseudocolors ** (n_barcode_rounds - 1))`
        unique regions. Each barcode is represented as a binary vector where each bit corresponds to a specific
        combination of barcode round, pseudocolor, and channel.

        The algorithm generates all possible barcode combinations by iterating through pseudocolor combinations
        and channels. Columns where all values are 0 (unused bits) are automatically removed from the final codebook.

        :param region_ids: List of region identifiers (e.g., gene IDs) to encode in the codebook.
            Each region will be assigned a unique barcode.
        :type region_ids: list[str]
        :param n_barcode_rounds: Number of barcode rounds in the encoding scheme. Each round contributes
            to the combinatorial encoding capacity.
        :type n_barcode_rounds: int
        :param n_pseudocolors: Number of pseudocolors available for encoding. Each barcode round uses
            one pseudocolor, and the final round includes a checksum pseudocolor.
        :type n_pseudocolors: int
        :param n_channels: Number of fluorescence channels available for detection. Each barcode is
            assigned to one channel.
        :type n_channels: int
        :return: A DataFrame containing the codebook with binary encoded bits. Each row represents a
            region's barcode, with columns named `bit_1`, `bit_2`, etc. Unused bit columns (all zeros)
            are automatically removed.
        :rtype: pd.DataFrame
        :raises ConfigurationError: If the number of valid barcodes (based on pseudocolors, barcode rounds,
            and channels) is insufficient for the number of regions. In this case, consider increasing
            `n_pseudocolors`, `n_barcode_rounds`, or `n_channels`, or reducing the number of regions.
        """

        def _generate_barcode(
            pseudocolors: list, channel: int, n_pseudocolors: int, n_channels: int
        ) -> np.ndarray:
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

        n_regions = len(region_ids)

        codebook_list: list[np.ndarray] = []
        codebook_size = n_channels * (n_pseudocolors ** (n_barcode_rounds - 1))
        barcode_size = n_pseudocolors * n_barcode_rounds * n_channels

        if codebook_size < n_regions:
            raise ConfigurationError(
                f"The number of valid barcodes ({codebook_size}) is lower than the number of regions ({n_regions}). "
                f"Consider increasing the number of pseudocolors or barcoding rounds, or reducing the number of regions."
            )
        for pseudocolors in product(range(n_pseudocolors), repeat=n_barcode_rounds - 1):
            for channel in range(n_channels):
                barcode = _generate_barcode(
                    pseudocolors=list(pseudocolors),
                    channel=channel,
                    n_pseudocolors=n_pseudocolors,
                    n_channels=n_channels,
                )
                codebook_list.append(barcode)

        codebook: pd.DataFrame = pd.DataFrame(
            codebook_list[0:n_regions], index=region_ids, columns=[f"bit_{i+1}" for i in range(barcode_size)]
        )

        # Remove columns where all values are 0
        codebook = codebook.loc[:, (codebook != 0).any(axis=0)]

        return codebook

    def create_readout_probe_table(
        self,
        readout_probe_database: OligoDatabase,
        channels_ids: list,
        n_barcode_rounds: int,
        n_pseudocolors: int,
    ) -> pd.DataFrame:
        """
        Create a readout probe table that maps codebook bits to barcode rounds, pseudocolors, channels, and readout probe sequences.

        This method generates a table where each bit position in the codebook is assigned:
        1. A readout probe sequence from the database
        2. A barcode round identifier (indicating which round of barcoding this bit belongs to)
        3. A pseudocolor identifier (indicating which pseudocolor is used in this round)
        4. A fluorescence channel identifier

        The readout probes are distributed across channels, pseudocolors, and barcode rounds in a systematic
        fashion. For each barcode round, pseudocolors are cycled through, and for each pseudocolor, channels
        are cycled through. This ensures that each bit position has a unique combination of barcode round,
        pseudocolor, and channel.

        The table is indexed by bit labels (bit_1, bit_2, etc.) and contains columns for barcode_round,
        pseudocolor, channel, and readout_probe_sequence. This table is used to assign readout probes
        to each bit position in the codebook for multiplexed imaging using the combinatorial pseudocolor system.

        :param readout_probe_database: The `OligoDatabase` instance containing readout probe sequences
            and their associated properties. This database should contain readout probes that have
            been filtered.
        :type readout_probe_database: OligoDatabase
        :param channels_ids: List of fluorescence channel identifiers (e.g., ['Cy3', 'Cy5', 'Alexa488'])
            to which readout probes will be assigned. Readout probes are distributed across channels
            in a systematic fashion.
        :type channels_ids: list
        :param n_barcode_rounds: Number of barcode rounds in the encoding scheme. This determines how
            many rounds of barcoding are used, and each round requires readout probes.
        :type n_barcode_rounds: int
        :param n_pseudocolors: Number of pseudocolors available for encoding. Each barcode round uses
            one pseudocolor, and readout probes are assigned to pseudocolors systematically.
        :type n_pseudocolors: int
        :return: A DataFrame containing the readout probe table with columns:
            - **barcode_round**: Barcode round identifier (1-indexed) for this bit
            - **pseudocolor**: Pseudocolor identifier (1-indexed) for this bit
            - **channel**: Fluorescence channel identifier assigned to this bit
            - **readout_probe_sequence**: DNA sequence of the readout probe
            The DataFrame is indexed by bit labels (bit_1, bit_2, etc.).
        :rtype: pd.DataFrame
        :raises AssertionError: If the number of available readout probes in the database is less
            than the number of required bits (`n_barcode_rounds * n_pseudocolors * n_channels`). In this
            case, generate more readout probes or reduce the number of barcode rounds, pseudocolors, or channels.
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
            readout_probe_table.iloc[i] = pd.Series(
                [
                    f"bit_{i + 1}",
                    barcode_round + 1,
                    pseudocolor + 1,
                    channels_ids[channel],
                    readout_probe_sequence,
                ],
                index=[
                    "bit",
                    "barcode_round",
                    "pseudocolor",
                    "channel",
                    "readout_probe_sequence",
                ],
            )
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


class PrimerDesigner:
    """
    A class for designing SeqFISH+ primers through a multi-step pipeline.

    This class provides methods for the complete primer design process, which includes:
    1. Creating an initial oligo database by generating random sequences with specified nucleotide
       probabilities and length (all ending with "T" nucleotide)
    2. Filtering primers based on sequence properties (GC content, GC clamp, homopolymeric runs,
       self-complementarity, complementarity to reverse primer, melting temperature, secondary structure)
    3. Filtering primers based on specificity to remove primers that bind to reference sequences
       or to the hybridization probes themselves using BLASTN searches

    The class is used to design forward primers that match a provided reverse primer's melting
    temperature for optimal PCR amplification. The resulting primers are used to amplify the
    hybridization probes during probe library preparation.

    :param dir_output: Directory path where output files will be saved.
    :type dir_output: str
    :param n_jobs: Number of parallel jobs to use for processing. Set to 1 for serial processing or higher
        values for parallel processing. This affects the parallelization of filtering, property calculation,
        and set generation operations.
    :type n_jobs: int
    """

    def __init__(
        self,
        dir_output: str,
        n_jobs: int,
    ) -> None:
        """Constructor for the PrimerDesigner class."""

        ##### create the output folder #####
        self.dir_output = os.path.abspath(dir_output)
        self.subdir_db_oligos = "db_primer"
        self.subdir_db_reference = "db_reference"

        self.n_jobs = n_jobs

    @pipeline_step_basic(step_name="Primer Generation - Create Oligo Database")
    def create_oligo_database(
        self,
        oligo_length: int,
        oligo_base_probabilities: dict,
        initial_num_sequences: int,
    ) -> OligoDatabase:
        """
        Create an initial oligo database by generating random sequences with specified nucleotide
        probabilities, all ending with a "T" nucleotide.

        This is the first step of primer design. The method generates random DNA sequences with
        user-defined nucleotide probabilities. All generated sequences end with a "T" nucleotide,
        which is a common requirement for PCR primers to improve binding stability. The sequences
        are created with length `oligo_length - 1`, then a "T" is appended to each sequence.

        :param oligo_length: Length (in nucleotides) of each primer sequence to generate. Note that
            sequences are generated with length `oligo_length - 1`, then a "T" nucleotide is appended,
            resulting in sequences of exactly `oligo_length`.
        :type oligo_length: int
        :param oligo_base_probabilities: Dictionary specifying the probability of each nucleotide base
            in randomly generated sequences. Keys should be 'A', 'T', 'G', 'C' and values should sum to 1.0
            (e.g., {"A": 0.25, "T": 0.25, "G": 0.25, "C": 0.25}).
        :type oligo_base_probabilities: dict
        :param initial_num_sequences: Number of random sequences to generate initially before filtering.
            Higher values provide more candidates but increase computation time.
        :type initial_num_sequences: int
        :return: An `OligoDatabase` object containing the generated random primer sequences. All sequences
            end with a "T" nucleotide. The database stores sequences with sequence type "oligo".
        :rtype: OligoDatabase
        """
        ##### creating the primer sequences #####
        # random forward primer
        forward_primer_sequences = OligoSequenceGenerator(dir_output=self.dir_output)
        forward_primer_fasta_file = forward_primer_sequences.create_sequences_random(
            filename_out="forward_primer_sequences",
            length_sequences=oligo_length - 1,
            num_sequences=initial_num_sequences,
            name_sequences="forward_primer",
            base_alphabet_with_probability=oligo_base_probabilities,
        )
        # we want to keep primer which end with a specific nucleotide, i.e. "T"
        forward_primer_fasta_file = append_nucleotide_to_sequences(forward_primer_fasta_file, nucleotide="T")

        ##### creating the primer database #####
        oligo_database = OligoDatabase(
            min_oligos_per_region=0,
            write_regions_with_insufficient_oligos=False,
            max_entries_in_memory=self.n_jobs * 2 + 2,
            database_name=self.subdir_db_oligos,
            dir_output=self.dir_output,
            n_jobs=1,
        )
        oligo_database.load_database_from_fasta(
            files_fasta=forward_primer_fasta_file,
            database_overwrite=True,
            sequence_type="oligo",
            region_ids=None,
        )
        # Set all sequence types that will be used in this pipeline
        oligo_database.set_database_sequence_types(["oligo"])

        dir = forward_primer_sequences.dir_output
        shutil.rmtree(dir) if os.path.exists(dir) else None

        return oligo_database

    @pipeline_step_basic(step_name="Primer Generation - Property Filters")
    def filter_by_property(
        self,
        oligo_database: OligoDatabase,
        GC_content_min: float,
        GC_content_max: float,
        number_GC_GCclamp: int,
        number_three_prime_base_GCclamp: int,
        homopolymeric_base_n: dict[str, int],
        max_len_selfcomplement: int,
        reverse_primer_sequence: str,
        max_len_complement: int,
        Tm_min: float,
        Tm_max: float,
        Tm_parameters: dict,
        Tm_chem_correction_parameters: dict | None,
        Tm_salt_correction_parameters: dict | None,
        T_secondary_structure: float,
        secondary_structures_threshold_deltaG: float,
    ) -> OligoDatabase:
        """
        Filter the oligo database based on sequence properties to remove primers with undesirable
        characteristics.

        This method applies sequential filtering using multiple property-based filters:
        1. **GC content**: Removes primers with GC content outside the specified range
        2. **GC clamp**: Removes primers that do not have sufficient G or C nucleotides at the 3' end
        3. **Homopolymeric runs**: Removes primers with homopolymeric runs exceeding specified lengths
        4. **Self-complementarity**: Removes primers with excessive self-complementary regions that can
           form hairpins and reduce PCR efficiency
        5. **Complementarity to reverse primer**: Removes primers with excessive complementarity to the
           reverse primer sequence, which prevents primer-dimer formation
        6. **Melting temperature**: Removes primers with Tm outside the specified range
        7. **Secondary structure**: Removes primers that form stable secondary structures at the
           specified temperature

        Probes that fail any filter are removed. Regions with insufficient oligos after filtering
        are removed from the database.

        :param oligo_database: The `OligoDatabase` instance containing oligonucleotide sequences
            and their associated properties. This database should contain primer sequences generated
            in the previous step.
        :type oligo_database: OligoDatabase
        :param GC_content_min: Minimum acceptable GC content for primers, expressed as a fraction
            between 0.0 and 1.0.
        :type GC_content_min: float
        :param GC_content_max: Maximum acceptable GC content for primers, expressed as a fraction
            between 0.0 and 1.0.
        :type GC_content_max: float
        :param number_GC_GCclamp: Minimum number of G or C nucleotides required within the specified
            number of bases at the 3' end (GC clamp). This improves primer binding stability.
        :type number_GC_GCclamp: int
        :param number_three_prime_base_GCclamp: Number of bases from the 3' end to consider for
            the GC clamp requirement.
        :type number_three_prime_base_GCclamp: int
        :param homopolymeric_base_n: Dictionary specifying the maximum allowed length of homopolymeric
            runs for each nucleotide base. Keys should be 'A', 'T', 'G', 'C' and values are the maximum
            run length. For example: {'A': 3, 'T': 3, 'G': 3, 'C': 3} allows up to 3 consecutive
            identical bases.
        :type homopolymeric_base_n: dict[str, int]
        :param max_len_selfcomplement: Maximum allowable length of self-complementary sequences.
            Primers with longer self-complementary regions can form hairpins and reduce PCR efficiency.
        :type max_len_selfcomplement: int
        :param reverse_primer_sequence: DNA sequence of the reverse primer that will be used for
            complementarity filtering. This prevents the forward and reverse primers from binding to
            each other.
        :type reverse_primer_sequence: str
        :param max_len_complement: Maximum allowable length of complementarity to the reverse primer
            sequence. This prevents the forward and reverse primers from binding to each other.
        :type max_len_complement: int
        :param Tm_min: Minimum acceptable melting temperature (Tm) for primers in degrees Celsius.
        :type Tm_min: float
        :param Tm_max: Maximum acceptable melting temperature (Tm) for primers in degrees Celsius.
        :type Tm_max: float
        :param Tm_parameters: Dictionary of parameters for calculating melting temperature (Tm) of primers
            using the nearest-neighbor method. For using Bio.SeqUtils.MeltingTemp default parameters, set to ``{}``.
            Common parameters include: 'nn_table', 'tmm_table', 'imm_table', 'de_table', 'dnac1', 'dnac2', 'Na', 'K',
            'Tris', 'Mg', 'dNTPs', 'saltcorr', etc. For more information on parameters, see:
            https://biopython.org/docs/1.75/api/Bio.SeqUtils.MeltingTemp.html#Bio.SeqUtils.MeltingTemp.Tm_NN
        :type Tm_parameters: dict
        :param Tm_chem_correction_parameters: Dictionary of chemical correction parameters for Tm calculation.
            These parameters account for the effects of chemical additives (e.g., DMSO, formamide) on melting temperature.
            Set to ``None`` to disable chemical correction, or set to ``{}`` to use Bio.SeqUtils.MeltingTemp default parameters.
            For more information, see:
            https://biopython.org/docs/1.75/api/Bio.SeqUtils.MeltingTemp.html#Bio.SeqUtils.MeltingTemp.chem_correction
        :type Tm_chem_correction_parameters: dict | None
        :param Tm_salt_correction_parameters: Dictionary of salt correction parameters for Tm calculation.
            These parameters account for the effects of salt concentration on melting temperature. Set to ``None`` to disable
            salt correction, or set to ``{}`` to use Bio.SeqUtils.MeltingTemp default parameters. For more information, see:
            https://biopython.org/docs/1.75/api/Bio.SeqUtils.MeltingTemp.html#Bio.SeqUtils.MeltingTemp.salt_correction
        :type Tm_salt_correction_parameters: dict | None
        :param T_secondary_structure: Temperature in degrees Celsius at which to evaluate secondary
            structure formation. Secondary structures that form at this temperature can interfere
            with primer binding.
        :type T_secondary_structure: float
        :param secondary_structures_threshold_deltaG: DeltaG threshold (in kcal/mol) for secondary
            structure stability. Primers with secondary structures having deltaG values more negative
            (more stable) than this threshold will be filtered out.
        :type secondary_structures_threshold_deltaG: float
        :return: A filtered `OligoDatabase` object containing only primers that pass all property filters.
            Regions with insufficient oligos after filtering are removed.
        :rtype: OligoDatabase
        """
        # define the filters
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
            Tm_chem_correction_parameters=Tm_chem_correction_parameters,
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
            oligo_database=oligo_database,
            sequence_type="oligo",
            n_jobs=self.n_jobs,
        )
        oligo_database.remove_regions_with_insufficient_oligos(pipeline_step="Property Filters")
        check_content_oligo_database(oligo_database)

        return oligo_database

    @pipeline_step_basic(step_name="Primer Generation - Specificity Filters")
    def filter_by_specificity(
        self,
        oligo_database: OligoDatabase,
        files_fasta_reference_database: list[str],
        specificity_refrence_blastn_search_parameters: dict,
        specificity_refrence_blastn_hit_parameters: dict,
        file_fasta_hybridization_probes_database: str,
        specificity_hybridization_probes_blastn_search_parameters: dict,
        specificity_hybridization_probes_blastn_hit_parameters: dict,
    ) -> OligoDatabase:
        """
        Filter the oligo database based on sequence specificity to remove primers that bind
        non-specifically to reference sequences or to the hybridization probes themselves.

        This method applies two types of specificity filters:

        1. **Reference database specificity filtering**: Removes primers that bind to unintended
           genomic regions using BLASTN searches against reference sequences (e.g., whole genome
           or transcriptome sequences). This ensures primers do not bind to off-target sites.

        2. **Hybridization probes specificity filtering**: Removes primers that bind to the hybridization
           probes themselves. This is critical because if primers can bind to the hybridization probes,
           they may interfere with probe function or cause unwanted amplification. The hybridization
           probe database is used to create a reference FASTA file for this filtering step.

        Both filters use BLASTN searches with configurable search and hit parameters. Primers with
        hits meeting the specified criteria are removed. Regions that do not meet the minimum oligo
        requirement after filtering are removed from the database.

        :param oligo_database: The `OligoDatabase` instance containing oligonucleotide sequences
            and their associated properties. This database should contain primer sequences that have
            passed property filtering.
        :type oligo_database: OligoDatabase
        :param files_fasta_reference_database: List of paths to FASTA files containing reference
            sequences used for specificity filtering. These files are used to identify off-target
            binding sites (e.g., whole genome or transcriptome sequences).
        :type files_fasta_reference_database: list[str]
        :param specificity_refrence_blastn_search_parameters: Dictionary of parameters for BLASTN
            searches used in specificity filtering against the reference database.
        :type specificity_refrence_blastn_search_parameters: dict
        :param specificity_refrence_blastn_hit_parameters: Dictionary of parameters for filtering
            BLASTN hits in specificity searches against the reference database.
            Primers with hits meeting these criteria are removed.
        :type specificity_refrence_blastn_hit_parameters: dict
        :param file_fasta_hybridization_probes_database: Path to the FASTA file containing
            hybridization probe sequences. This file is used to create a reference database
            for specificity filtering to ensure primers do not bind to the hybridization probes themselves.
        :type file_fasta_hybridization_probes_database: str
        :param specificity_hybridization_probes_blastn_search_parameters: Dictionary of parameters for BLASTN
            searches used in specificity filtering against the hybridization probes database.
        :type specificity_hybridization_probes_blastn_search_parameters: dict
        :param specificity_hybridization_probes_blastn_hit_parameters: Dictionary of parameters for filtering
            BLASTN hits in specificity searches against the hybridization probes database. Primers with hits meeting these
            criteria are removed.
        :type specificity_hybridization_probes_blastn_hit_parameters: dict
        :return: A filtered `OligoDatabase` object containing only primers that pass all specificity
            filters. Regions with insufficient oligos after filtering are removed.
        :rtype: OligoDatabase
        """
        ##### specificity filters against reference #####
        reference_database = ReferenceDatabase(
            database_name=self.subdir_db_reference, dir_output=self.dir_output
        )
        reference_database.load_database_from_file(
            files=files_fasta_reference_database, file_type="fasta", database_overwrite=True
        )
        # BlastN Filter
        specificity_refrence = BlastNFilter(
            search_parameters=specificity_refrence_blastn_search_parameters,
            hit_parameters=specificity_refrence_blastn_hit_parameters,
            filter_name="primer_blastn_specificity_reference",
            dir_output=self.dir_output,
        )
        specificity_refrence.set_reference_database(reference_database=reference_database)

        ##### specificity filters against hybridization probes #####
        hybridization_probes_database = ReferenceDatabase(
            database_name=self.subdir_db_reference, dir_output=self.dir_output
        )
        hybridization_probes_database.load_database_from_file(
            files=file_fasta_hybridization_probes_database, file_type="fasta", database_overwrite=True
        )
        # BlastN Filter
        specificity_hybridization_probes = BlastNFilter(
            search_parameters=specificity_hybridization_probes_blastn_search_parameters,
            hit_parameters=specificity_hybridization_probes_blastn_hit_parameters,
            filter_name="primer_blastn_specificity_hybridization_probes",
            dir_output=self.dir_output,
        )
        specificity_hybridization_probes.set_reference_database(
            reference_database=hybridization_probes_database
        )

        specificity_filter = SpecificityFilter(
            filters=[specificity_refrence, specificity_hybridization_probes]
        )
        oligo_database = specificity_filter.apply(
            oligo_database=oligo_database,
            sequence_type="oligo",
            n_jobs=self.n_jobs,
        )

        # remove all directories of intermediate steps
        for directory in [
            specificity_refrence.dir_output,
            specificity_hybridization_probes.dir_output,
        ]:
            if os.path.exists(directory):
                shutil.rmtree(directory)

        oligo_database.remove_regions_with_insufficient_oligos(pipeline_step="Specificity Filters")
        check_content_oligo_database(oligo_database)

        return oligo_database


############################################
# SeqFish Plus Probe Designer Pipeline
############################################


def main() -> None:
    """
    Main entry point for running the SeqFish Plus probe design pipeline.

    This function orchestrates the complete SeqFish Plus probe design workflow:
    1. Parses command-line arguments using the base parser
    2. Reads the configuration YAML file containing all pipeline parameters
    3. Reads the gene IDs file (if provided) or uses all genes from FASTA files
    4. Preprocesses melting temperature parameters for target probes, readout probes, and primers
    5. Initializes the SeqFishPlusProbeDesigner pipeline
    6. Designs target probes for specified genes
    7. Designs readout probes and generates the codebook
    8. Assembles hybridization probes by combining target probes with readout probe barcodes
    9. Designs forward and reverse primers for PCR amplification
    10. Assembles final DNA template probes with primers
    11. Generates output files (codebook, readout probe table, probe sequences, etc.)

    The function is typically called from the command line:
    ``seqfish_plus_probe_designer --config <path_to_config.yaml>``

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
    pipeline = SeqFishPlusProbeDesigner(
        write_intermediate_steps=config["write_intermediate_steps"],
        dir_output=config["dir_output"],
        n_jobs=config["n_jobs"],
    )

    ##### design probes #####
    target_probe_database = pipeline.design_target_probes(
        # Step 1: Create Database Parameters
        region_ids=region_ids,
        files_fasta_target_probe_database=config["files_fasta_target_probe_database"],
        target_probe_length_min=config["target_probe_length_min"],
        target_probe_length_max=config["target_probe_length_max"],
        target_probe_isoform_consensus=config["target_probe_isoform_consensus"],
        # Step 2: Property Filter Parameters
        target_probe_GC_content_min=config["target_probe_GC_content_min"],
        target_probe_GC_content_opt=config["target_probe_GC_content_opt"],
        target_probe_GC_content_max=config["target_probe_GC_content_max"],
        target_probe_homopolymeric_base_n=config["target_probe_homopolymeric_base_n"],
        target_probe_T_secondary_structure=config["target_probe_T_secondary_structure"],
        target_probe_secondary_structures_threshold_deltaG=config[
            "target_probe_secondary_structures_threshold_deltaG"
        ],
        # Step 3: Specificity Filter Parameters
        files_fasta_reference_database_target_probe=config["files_fasta_reference_database_target_probe"],
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
        # Step 4: Probe Scoring and Set Selection Parameters
        target_probe_GC_weight=config["target_probe_GC_weight"],
        target_probe_UTR_weight=config["target_probe_UTR_weight"],
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
        # Step 1: Create Database Parameters
        region_ids=list(target_probe_database.database.keys()),
        readout_probe_length=config["readout_probe_length"],
        readout_probe_base_probabilities=config["readout_probe_base_probabilities"],
        readout_probe_initial_num_sequences=config["readout_probe_initial_num_sequences"],
        # Step 2: Property Filter Parameters
        readout_probe_GC_content_min=config["readout_probe_GC_content_min"],
        readout_probe_GC_content_max=config["readout_probe_GC_content_max"],
        readout_probe_homopolymeric_base_n=config["readout_probe_homopolymeric_base_n"],
        # Step 3: Specificity Filter Parameters
        files_fasta_reference_database_readout_probe=config["files_fasta_reference_database_readout_probe"],
        readout_probe_specificity_blastn_search_parameters=config[
            "readout_probe_specificity_blastn_search_parameters"
        ],
        readout_probe_specificity_blastn_hit_parameters=config[
            "readout_probe_specificity_blastn_hit_parameters"
        ],
        readout_probe_cross_hybridization_blastn_search_parameters=config[
            "readout_probe_cross_hybridization_blastn_search_parameters"
        ],
        readout_probe_cross_hybridization_blastn_hit_parameters=config[
            "readout_probe_cross_hybridization_blastn_hit_parameters"
        ],
        # Step 4: Codebook and Readout Probe Table Parameters
        n_barcode_rounds=config["n_barcode_rounds"],
        n_pseudocolors=config["n_pseudocolors"],
        channels_ids=config["channels_ids"],
    )

    hybridization_probe_database = pipeline.assemble_hybridization_probes(
        target_probe_database=target_probe_database,
        codebook=codebook,
        readout_probe_table=readout_probe_table,
    )

    reverse_primer_sequence, forward_primer_sequence = pipeline.design_primers(
        # Step 1: Create Database Parameters
        reverse_primer_sequence=config["reverse_primer_sequence"],
        primer_length=config["primer_length"],
        primer_base_probabilities=config["primer_base_probabilities"],
        primer_initial_num_sequences=config["primer_initial_num_sequences"],
        # Step 2: Property Filter Parameters
        primer_GC_content_min=config["primer_GC_content_min"],
        primer_GC_content_max=config["primer_GC_content_max"],
        primer_number_GC_GCclamp=config["primer_number_GC_GCclamp"],
        primer_number_three_prime_base_GCclamp=config["primer_number_three_prime_base_GCclamp"],
        primer_homopolymeric_base_n=config["primer_homopolymeric_base_n"],
        primer_max_len_selfcomplement=config["primer_max_len_selfcomplement"],
        primer_max_len_complement_reverse_primer=config["primer_max_len_complement_reverse_primer"],
        primer_Tm_min=config["primer_Tm_min"],
        primer_Tm_max=config["primer_Tm_max"],
        primer_T_secondary_structure=config["primer_T_secondary_structure"],
        primer_secondary_structures_threshold_deltaG=config["primer_secondary_structures_threshold_deltaG"],
        primer_Tm_parameters=preprocess_tm_parameters(config["primer_Tm_parameters"]),
        primer_Tm_chem_correction_parameters=config["primer_Tm_chem_correction_parameters"],
        primer_Tm_salt_correction_parameters=config["primer_Tm_salt_correction_parameters"],
        # Step 3: Specificity Filter Parameters
        hybridization_probe_database=hybridization_probe_database,
        files_fasta_reference_database_primer=config["files_fasta_reference_database_primer"],
        primer_specificity_refrence_blastn_search_parameters=config[
            "primer_specificity_refrence_blastn_search_parameters"
        ],
        primer_specificity_refrence_blastn_hit_parameters=config[
            "primer_specificity_refrence_blastn_hit_parameters"
        ],
        primer_specificity_hybridization_probes_blastn_search_parameters=config[
            "primer_specificity_hybridization_probes_blastn_search_parameters"
        ],
        primer_specificity_hybridization_probes_blastn_hit_parameters=config[
            "primer_specificity_hybridization_probes_blastn_hit_parameters"
        ],
    )

    probe_database = pipeline.assemble_dna_template_probes(
        hybridization_probe_database=hybridization_probe_database,
        reverse_primer_sequence=reverse_primer_sequence,
        forward_primer_sequence=forward_primer_sequence,
    )

    pipeline.generate_output(
        probe_database=probe_database,
        codebook=codebook,
        readout_probe_table=readout_probe_table,
    )

    logging.info("--------------END PIPELINE--------------")


if __name__ == "__main__":
    main()
