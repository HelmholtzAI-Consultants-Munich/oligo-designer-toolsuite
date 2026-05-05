############################################
# imports
############################################

import itertools
import logging
import os
import random
import shutil
import warnings
from pathlib import Path

import yaml
from joblib import Parallel, delayed
from joblib_progress import joblib_progress

from oligo_designer_toolsuite.database import OligoDatabase, ReferenceDatabase
from oligo_designer_toolsuite.oligo_efficiency_filter import (
    IsoformConsensusScorer,
    LowestSetScoring,
    NormalizedDeviationFromOptimalGCContentScorer,
    NormalizedDeviationFromOptimalTmScorer,
    OligoScoring,
)
from oligo_designer_toolsuite.oligo_property_calculator import (
    GCContentProperty,
    IsoformConsensusProperty,
    LengthProperty,
    PadlockArmsProperty,
    PropertyCalculator,
    ReverseComplementSequenceProperty,
    TmNNProperty,
)
from oligo_designer_toolsuite.oligo_property_calculator._property_functions import (
    calc_detect_oligo,
    calc_tm_nn,
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
    pipeline_step_basic,
    preprocess_tm_parameters,
    setup_logging,
)
from oligo_designer_toolsuite.sequence_generator import OligoSequenceGenerator

############################################
# SCRINSHOT Probe Designer
############################################


class ScrinshotProbeDesigner:
    """
    A class for designing padlock and detection probes for SCRINSHOT (Single-Cell Resolution IN Situ Hybridization On Tissues) experiments.

    This class provides a comprehensive pipeline for designing padlock probes and detection oligonucleotides compatible with the SCRINSHOT
    method, a multiplex, single-cell–resolution RNA mapping approach that enables spatial transcriptomics in fixed tissue sections.

    **SCRINSHOT Pipeline Overview:**
    - **Target Probe Design**: Design gene-specific targeting sequences (~40-45 nt) that bind to RNA transcripts
    - **Detection Oligo Design**: Generate 30–35 nt UNG-cleavable detection oligos centered on the ligation site.
    - **Padlock Backbone Assembly**: Combine 5' arm + constant backbone (53 nt) + 3' arm to form full padlock sequences and record ligation-site coordinates.
    - **Output Generation**: Generate output files in multiple formats (TSV, YAML)

    Overview
    --------
    SCRINSHOT (Single-Cell Resolution IN Situ Hybridization On Tissues) is a targeted spatial transcriptomics approach for multiplex
    detection of RNA molecules in fixed tissue sections with single-cell resolution. It combines **direct padlock probe hybridization on RNA**,
    **SplintR ligase–mediated circularization**, and **rolling circle amplification (RCA)** to generate bright, quantifiable signals from individual transcripts.

    By bypassing reverse transcription and using optimized probe design and stringent hybridization conditions, SCRINSHOT achieves high sensitivity,
    specificity, and quantitative performance comparable to scRNA-seq data, across a wide range of expression levels. The method enables spatial mapping
    of abundant and rare cell types across diverse tissues (e.g., lung, heart, kidney, brain) and is compatible with standard epifluorescence microscopy.

    Probe Structure
    ---------------
    **Padlock Probes**
    - Single-stranded DNA oligonucleotides designed to hybridize directly to target RNA sequences.
    - Each probe is composed of:
        - **Target-specific arms**: Each arm is approximately 20 nucleotides, complementary to adjacent regions of the target mRNA that flank the ligation site (Tm ≈ 50–60 °C).
        - **Composite backbone** providing priming and detection regions:
            - accessory sequence 1 = "TCCTCTATGATTACTGAC"`
            - ISS anchor sequence = "TGCGTCTATTTAGTGGAGCC"`
            - accessory sequence 2 = "CTATCTTCTTT"`
            - backbone sequence = [accessory sequence 1] + [ISS anchor sequence] + [barcode] + [accessory sequence 2]
        - **Full probe assembly**: [padlock arm 1] + [backbone sequence] + [padlock arm 2]
    - The ligation junction lies between both arms, enabling circularization by **SplintR ligase**.
    - After ligation, the circularized probe serves as a template for **rolling circle amplification (RCA)**, producing long concatemeric RCA products.
    - RCA products are detected using complementary fluorophore-labeled detection oligos.

    **Detection Oligos**
    - Short (~30–35 nt) single-stranded DNA probes complementary to the gene-specific region of RCA products.
    - Designed with the ligation site centered within the oligo and a melting temperature around 56 °C.
    - Include 2–3 **uracil (U)** substitutions spaced ≤ 10 nt apart to allow enzymatic cleavage by **Uracil DNA Glycosylase (UNG)**, facilitating sequential hybridization cycles.
    - Labeled at the 3' end with fluorophores (FITC, Cy3, Cy5; optionally Texas Red or Atto740 for extended color sets).

    References
    ----------
    Sountoulidis, A., Liontos, A., Nguyen, H. P., Firsova, A. B., Fysikopoulos, A., Qian, X., et al. (2020).
    SCRINSHOT enables spatial mapping of cell states in tissue sections with single-cell resolution.
    *PLOS Biology*, 18(11): e3000675. https://doi.org/10.1371/journal.pbio.3000675

    :ivar dir_output: Directory path where output probe design files will be saved.
    :type dir_output: str
    :ivar write_intermediate_steps: Whether to save intermediate probe and validation data (default: False).
    :type write_intermediate_steps: bool
    :ivar n_jobs: Number of parallel threads to use for sequence design and BLAST validation.
    :type n_jobs: int
    """

    def __init__(self, write_intermediate_steps: bool, dir_output: str, n_jobs: int) -> None:
        """Constructor for the ScrinshotProbeDesigner class."""

        # create the output folder
        self.dir_output = os.path.abspath(dir_output)
        Path(self.dir_output).mkdir(parents=True, exist_ok=True)

        # setup logger
        setup_logging(
            dir_output=self.dir_output,
            pipeline_name="scrinshot_probe_designer",
            log_start_message=True,
        )

        ##### set class parameters #####
        self.write_intermediate_steps = write_intermediate_steps
        self.n_jobs = n_jobs

    def design_target_probes(
        self,
        # Step 1: Create Database Parameters
        region_ids: list | None,
        files_fasta_target_probe_database: list,
        target_probe_length_min: int,
        target_probe_length_max: int,
        target_probe_isoform_consensus: float,
        # Step 2: Property Filter Parameters
        target_probe_GC_content_min: int,
        target_probe_GC_content_max: int,
        target_probe_Tm_min: int,
        target_probe_Tm_max: int,
        target_probe_homopolymeric_base_n: dict,
        detection_oligo_min_thymines: int,
        detection_oligo_length_min: int,
        detection_oligo_length_max: int,
        target_probe_padlock_arm_length_min: int,
        target_probe_padlock_arm_Tm_dif_max: int,
        target_probe_padlock_arm_Tm_min: int,
        target_probe_padlock_arm_Tm_max: int,
        target_probe_Tm_parameters: dict,
        target_probe_Tm_chem_correction_parameters: dict | None,
        target_probe_Tm_salt_correction_parameters: dict | None,
        # Step 3: Specificity Filter Parameters
        files_fasta_reference_database_target_probe: list,
        target_probe_specificity_blastn_search_parameters: dict,
        target_probe_specificity_blastn_hit_parameters: dict,
        target_probe_cross_hybridization_blastn_search_parameters: dict,
        target_probe_cross_hybridization_blastn_hit_parameters: dict,
        target_probe_ligation_region_size: int,
        # Step 4: Probe Scoring and Set Selection Parameters
        target_probe_isoform_weight: float,
        target_probe_GC_content_opt: int,
        target_probe_GC_weight: float,
        target_probe_Tm_opt: int,
        target_probe_Tm_weight: float,
        set_size_min: int,
        set_size_opt: int,
        distance_between_target_probes: int,
        n_sets: int,
        n_attempts_graph: int,
        n_attempts_clique_enum: int,
        diversification_fraction: float,
        jaccard_opt: float,
        jaccard_step: float,
    ) -> OligoDatabase:
        """
        Design target probes for SCRINSHOT experiments through a multi-step pipeline.

        This method performs the complete target probe design process, which includes:
        1. Creating an initial oligo database from input FASTA files using a sliding window approach
        2. Filtering probes based on sequence properties (GC content, melting temperature, homopolymeric
           runs, detection oligo requirements, padlock arm requirements)
        3. Filtering probes based on specificity to remove off-target binding and cross-hybridization
           using BLASTN searches
        4. Organizing filtered probes into optimal sets based on weighted scoring criteria (isoform
           consensus, GC content, melting temperature) and distance constraints

        The resulting probes are gene-specific targeting sequences (typically 40-45 nt) that bind to RNA
        transcripts. These probes will later be split into padlock arms and combined with a backbone
        sequence to create complete padlock probes.

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
        :type target_probe_GC_content_min: int
        :param target_probe_GC_content_max: Maximum acceptable GC content for target probes, expressed
            as a fraction between 0.0 and 1.0.
        :type target_probe_GC_content_max: int
        :param target_probe_Tm_min: Minimum acceptable melting temperature (Tm) for target probes in
            degrees Celsius.
        :type target_probe_Tm_min: int
        :param target_probe_Tm_max: Maximum acceptable melting temperature (Tm) for target probes in
            degrees Celsius.
        :type target_probe_Tm_max: int
        :param target_probe_homopolymeric_base_n: Dictionary specifying the maximum allowed length of
            homopolymeric runs for each nucleotide base (keys: 'A', 'T', 'G', 'C').
        :type target_probe_homopolymeric_base_n: dict[str, int]
        :param detection_oligo_min_thymines: Minimum number of thymine (T) nucleotides required in the
            detection oligo region. These thymines will be converted to uracils (U) for UNG cleavage.
        :type detection_oligo_min_thymines: int
        :param detection_oligo_length_min: Minimum length (in nucleotides) for detection oligo sequences.
        :type detection_oligo_length_min: int
        :param detection_oligo_length_max: Maximum length (in nucleotides) for detection oligo sequences.
        :type detection_oligo_length_max: int
        :param target_probe_padlock_arm_length_min: Minimum length (in nucleotides) for each padlock
            probe arm. Each arm must meet this requirement for the probe to pass filtering.
        :type target_probe_padlock_arm_length_min: int
        :param target_probe_padlock_arm_Tm_dif_max: Maximum allowed difference in melting temperature
            (in degrees Celsius) between the two padlock arms. This ensures balanced binding of both arms.
        :type target_probe_padlock_arm_Tm_dif_max: int
        :param target_probe_padlock_arm_Tm_min: Minimum acceptable melting temperature (Tm) for padlock
            arms in degrees Celsius.
        :type target_probe_padlock_arm_Tm_min: int
        :param target_probe_padlock_arm_Tm_max: Maximum acceptable melting temperature (Tm) for padlock
            arms in degrees Celsius.
        :type target_probe_padlock_arm_Tm_max: int
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
        :param target_probe_ligation_region_size: Size of the ligation region (in nucleotides) around
            the ligation site. This parameter is used for seed-based specificity filtering around the
            junction region where padlock arms meet.
        :type target_probe_ligation_region_size: int

        **Step 4: Probe Scoring and Set Selection Parameters**

        :param target_probe_isoform_weight: Weight assigned to isoform consensus in the scoring function.
        :type target_probe_isoform_weight: float
        :param target_probe_GC_content_opt: Optimal GC content for target probes, expressed as a fraction
            between 0.0 and 1.0. Used in scoring to prioritize probes closer to this value.
        :type target_probe_GC_content_opt: int
        :param target_probe_GC_weight: Weight assigned to GC content in the scoring function.
        :type target_probe_GC_weight: float
        :param target_probe_Tm_opt: Optimal melting temperature (Tm) for target probes in degrees Celsius.
            Used in scoring to prioritize probes closer to this value.
        :type target_probe_Tm_opt: int
        :param target_probe_Tm_weight: Weight assigned to melting temperature in the scoring function.
        :type target_probe_Tm_weight: float
        :param set_size_min: Minimum size (number of probes) required for each oligo set. Sets with fewer probes than
            this value will be rejected, and regions that cannot generate sets meeting this minimum will be removed.
        :type set_size_min: int
        :param set_size_opt: Optimal size (number of probes) for each oligo set. The set selection algorithm will
            attempt to generate sets of this size, but may produce sets with fewer probes if constraints cannot be met.
        :type set_size_opt: int
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
            dir_database = oligo_database.save_database(name_database="1_db_probes_initial")
            logging.info(f"Saved probe database for step 1 (Create Database) in directory {dir_database}")

        oligo_database = target_probe_designer.filter_by_property(
            oligo_database=oligo_database,
            GC_content_min=target_probe_GC_content_min,
            GC_content_max=target_probe_GC_content_max,
            Tm_min=target_probe_Tm_min,
            Tm_max=target_probe_Tm_max,
            Tm_parameters=target_probe_Tm_parameters,
            Tm_chem_correction_parameters=target_probe_Tm_chem_correction_parameters,
            Tm_salt_correction_parameters=target_probe_Tm_salt_correction_parameters,
            homopolymeric_base_n=target_probe_homopolymeric_base_n,
            detect_oligo_length_min=detection_oligo_length_min,
            detect_oligo_length_max=detection_oligo_length_max,
            min_thymines=detection_oligo_min_thymines,
            arm_length_min=target_probe_padlock_arm_length_min,
            arm_Tm_dif_max=target_probe_padlock_arm_Tm_dif_max,
            arm_Tm_min=target_probe_padlock_arm_Tm_min,
            arm_Tm_max=target_probe_padlock_arm_Tm_max,
        )

        if self.write_intermediate_steps:
            dir_database = oligo_database.save_database(name_database="2_db_probes_property_filter")
            logging.info(f"Saved probe database for step 2 (Property Filters) in directory {dir_database}")

        oligo_database = target_probe_designer.filter_by_specificity(
            oligo_database=oligo_database,
            files_fasta_reference_database=files_fasta_reference_database_target_probe,
            specificity_blastn_search_parameters=target_probe_specificity_blastn_search_parameters,
            specificity_blastn_hit_parameters=target_probe_specificity_blastn_hit_parameters,
            cross_hybridization_blastn_search_parameters=target_probe_cross_hybridization_blastn_search_parameters,
            cross_hybridization_blastn_hit_parameters=target_probe_cross_hybridization_blastn_hit_parameters,
            ligation_region_size=target_probe_ligation_region_size,
            arm_length_min=target_probe_padlock_arm_length_min,
            arm_Tm_dif_max=target_probe_padlock_arm_Tm_dif_max,
            arm_Tm_min=target_probe_padlock_arm_Tm_min,
            arm_Tm_max=target_probe_padlock_arm_Tm_max,
            Tm_parameters=target_probe_Tm_parameters,
            Tm_chem_correction_parameters=target_probe_Tm_chem_correction_parameters,
            Tm_salt_correction_parameters=target_probe_Tm_salt_correction_parameters,
        )

        if self.write_intermediate_steps:
            dir_database = oligo_database.save_database(name_database="3_db_probes_specificity_filter")
            logging.info(f"Saved probe database for step 3 (Specificity Filters) in directory {dir_database}")

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

        # Calculate oligo length, GC content, Tm, and isoform consensus
        length_property = LengthProperty()
        gc_content_property = GCContentProperty()
        TmNN_property = TmNNProperty(
            Tm_parameters=target_probe_Tm_parameters,
            Tm_chem_correction_parameters=target_probe_Tm_chem_correction_parameters,
            Tm_salt_correction_parameters=target_probe_Tm_salt_correction_parameters,
        )
        isoform_consensus_property = IsoformConsensusProperty()
        calculator = PropertyCalculator(
            properties=[length_property, gc_content_property, TmNN_property, isoform_consensus_property]
        )

        oligo_database = calculator.apply(
            oligo_database=oligo_database, sequence_type="oligo", n_jobs=self.n_jobs
        )

        if self.write_intermediate_steps:
            dir_database = oligo_database.save_database(name_database="4_db_probes_probesets")
            logging.info(
                f"Saved probe database for step 4 (Specificity Filters) in directory {dir_database}."
            )

        return oligo_database

    def design_detection_oligos(
        self,
        oligo_database: OligoDatabase,
        detection_oligo_length_min: int,
        detection_oligo_length_max: int,
        detection_oligo_min_thymines: int,
        detection_oligo_U_distance: int,
        detection_oligo_Tm_opt: float,
        detection_oligo_Tm_parameters: dict,
        detection_oligo_Tm_chem_correction_parameters: dict | None,
        detection_oligo_Tm_salt_correction_parameters: dict | None,
    ) -> OligoDatabase:
        """
        Design detection oligonucleotides for SCRINSHOT padlock probes.

        This method generates detection oligos that are complementary to the gene-specific region of
        rolling circle amplification (RCA) products. Detection oligos are designed with:
        1. The ligation site centered within the oligo sequence
        2. A target melting temperature around 56 °C for optimal hybridization
        3. 2-3 uracil (U) substitutions spaced ≤ 10 nt apart to allow enzymatic cleavage by
           Uracil DNA Glycosylase (UNG), facilitating sequential hybridization cycles

        The detection oligos are created by extracting sequences centered on the ligation site from
        the target probes, then converting thymines (T) to uracils (U) at appropriate positions.
        The resulting sequences are stored as properties in the database for each probe.

        :param oligo_database: The `OligoDatabase` instance containing target probes with their
            sequences and properties. This database should contain target probes organized by region IDs,
            with each region having one or more probe sets and ligation site information.
        :type oligo_database: OligoDatabase
        :param detection_oligo_length_min: Minimum length (in nucleotides) for detection oligo sequences.
        :type detection_oligo_length_min: int
        :param detection_oligo_length_max: Maximum length (in nucleotides) for detection oligo sequences.
        :type detection_oligo_length_max: int
        :param detection_oligo_min_thymines: Minimum number of thymine (T) nucleotides required in the
            detection oligo sequence. These thymines will be converted to uracils (U) for UNG cleavage.
        :type detection_oligo_min_thymines: int
        :param detection_oligo_U_distance: Maximum distance (in nucleotides) allowed between uracil
            substitutions in the detection oligo. Uracils must be spaced ≤ this distance apart.
        :type detection_oligo_U_distance: int
        :param detection_oligo_Tm_opt: Optimal melting temperature (Tm) for detection oligos in degrees
            Celsius. The algorithm will attempt to select detection oligos with Tm closest to this value.
        :type detection_oligo_Tm_opt: float
        :param detection_oligo_Tm_parameters: Dictionary of parameters for calculating melting temperature (Tm) of detection
            oligos using the nearest-neighbor method. For using Bio.SeqUtils.MeltingTemp default parameters, set to ``{}``.
            Common parameters include: 'nn_table', 'tmm_table', 'imm_table', 'de_table', 'dnac1', 'dnac2', 'Na', 'K',
            'Tris', 'Mg', 'dNTPs', 'saltcorr', etc. For more information on parameters, see:
            https://biopython.org/docs/1.75/api/Bio.SeqUtils.MeltingTemp.html#Bio.SeqUtils.MeltingTemp.Tm_NN
        :type detection_oligo_Tm_parameters: dict
        :param detection_oligo_Tm_chem_correction_parameters: Dictionary of chemical correction parameters for Tm calculation.
            These parameters account for the effects of chemical additives (e.g., DMSO, formamide) on melting temperature.
            Set to ``None`` to disable chemical correction, or set to ``{}`` to use Bio.SeqUtils.MeltingTemp default parameters.
            For more information, see:
            https://biopython.org/docs/1.75/api/Bio.SeqUtils.MeltingTemp.html#Bio.SeqUtils.MeltingTemp.chem_correction
        :type detection_oligo_Tm_chem_correction_parameters: dict | None
        :param detection_oligo_Tm_salt_correction_parameters: Dictionary of salt correction parameters for Tm calculation.
            These parameters account for the effects of salt concentration on melting temperature. Set to ``None`` to disable
            salt correction, or set to ``{}`` to use Bio.SeqUtils.MeltingTemp default parameters. For more information, see:
            https://biopython.org/docs/1.75/api/Bio.SeqUtils.MeltingTemp.html#Bio.SeqUtils.MeltingTemp.salt_correction
        :type detection_oligo_Tm_salt_correction_parameters: dict | None
        :return: An updated `OligoDatabase` object containing the designed detection oligos. The
            database includes the following new sequence properties for each probe:
            - `sequence_detection_oligo`: The detection oligo sequence with uracil substitutions
            - `Tm_detection_oligo`: The melting temperature of the detection oligo
        :rtype: OligoDatabase
        """

        detection_oligo_designer = DetectionOligoDesigner(self.n_jobs)
        oligo_database = detection_oligo_designer.create_detection_oligos(
            oligo_database=oligo_database,
            oligo_length_min=detection_oligo_length_min,
            oligo_length_max=detection_oligo_length_max,
            min_thymines=detection_oligo_min_thymines,
            U_distance=detection_oligo_U_distance,
            Tm_opt=detection_oligo_Tm_opt,
            Tm_parameters=detection_oligo_Tm_parameters,
            Tm_chem_correction_parameters=detection_oligo_Tm_chem_correction_parameters,
            Tm_salt_correction_parameters=detection_oligo_Tm_salt_correction_parameters,
        )

        return oligo_database

    def assemble_padlock_backbone(
        self,
        oligo_database: OligoDatabase,
        target_probe_Tm_parameters: dict,
        target_probe_Tm_chem_correction_parameters: dict | None,
        target_probe_Tm_salt_correction_parameters: dict | None,
    ) -> OligoDatabase:
        """
        Assemble padlock probes by combining target probe arms with a constant backbone sequence.

        This method creates complete SCRINSHOT padlock probes by:
        1. Splitting each target probe sequence at the ligation site into two arms (5' arm and 3' arm)
        2. Generating a unique barcode for each region (gene)
        3. Constructing the composite backbone sequence with the structure:
           [accessory sequence 1] + [ISS anchor sequence] + [barcode] + [accessory sequence 2]
        4. Assembling the full padlock probe with the structure:
           [padlock arm 1] + [backbone sequence] + [padlock arm 2]
        5. Calculating melting temperatures for both arms to verify balanced binding

        The ligation junction lies between both arms, enabling circularization by SplintR ligase.
        After ligation, the circularized probe serves as a template for rolling circle amplification (RCA).

        :param oligo_database: The `OligoDatabase` instance containing target probes with their
            sequences, ligation sites, and properties. This database should contain target probes
            organized by region IDs, with each region having one or more probe sets.
        :type oligo_database: OligoDatabase
        :param target_probe_Tm_parameters: Dictionary of parameters for calculating melting temperature (Tm) of padlock
            arms using the nearest-neighbor method. For using Bio.SeqUtils.MeltingTemp default parameters, set to ``{}``.
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
        :return: An updated `OligoDatabase` object containing the assembled padlock probes. The
            database includes the following new sequence properties for each probe:
            - `barcode`: The unique barcode sequence assigned to this region
            - `sequence_target`: The original gene-specific targeting sequence
            - `sequence_padlock_arm1`: The 5' padlock arm sequence (from ligation site to end)
            - `sequence_padlock_arm2`: The 3' padlock arm sequence (from start to ligation site)
            - `sequence_padlock_accessory1`: The first accessory sequence ("TCCTCTATGATTACTGAC")
            - `sequence_padlock_ISS_anchor`: The ISS anchor sequence ("TGCGTCTATTTAGTGGAGCC")
            - `sequence_padlock_accessory2`: The second accessory sequence ("CTATCTTCTTT")
            - `sequence_padlock_backbone`: The complete backbone sequence (accessory1 + ISS anchor + barcode + accessory2)
            - `sequence_padlock_probe`: The complete assembled padlock probe sequence
            - `Tm_arm1`: The melting temperature of arm 1
            - `Tm_arm2`: The melting temperature of arm 2
            - `Tm_diff_arms`: The absolute difference in melting temperature between the two arms
        :rtype: OligoDatabase
        """

        def _get_barcode(number_regions: int, barcode_length: int, seed: int, choices: list) -> list[str]:

            while len(choices) ** barcode_length < number_regions:
                barcode_length += 1

            barcodes: list[str] = ["".join(nts) for nts in itertools.product(choices, repeat=barcode_length)]
            random.seed(seed)
            random.shuffle(barcodes)

            return barcodes

        region_ids = list(oligo_database.database.keys())

        barcodes = _get_barcode(len(region_ids), barcode_length=4, seed=0, choices=["A", "C", "T", "G"])

        for region_idx, region_id in enumerate(region_ids):
            oligo_sets_region = oligo_database.oligosets[region_id]
            oligo_sets_oligo_columns = [col for col in oligo_sets_region.columns if col.startswith("oligo_")]

            new_oligo_properties = {}

            for index in range(len(oligo_sets_region.index)):
                for column in oligo_sets_oligo_columns:
                    oligo_id = str(oligo_sets_region.loc[index, column])
                    barcode: str = barcodes[region_idx]

                    ligation_site = oligo_database.get_oligo_property_value(
                        property="ligation_site", region_id=region_id, oligo_id=oligo_id, flatten=True
                    )
                    sequence_oligo = oligo_database.get_oligo_property_value(
                        property="oligo", region_id=region_id, oligo_id=oligo_id, flatten=True
                    )
                    # required for type linting since get_oligo_property_value() could return None
                    if not isinstance(sequence_oligo, str) or not isinstance(ligation_site, int):
                        continue
                    sequence_padlock_arm1: str = sequence_oligo[ligation_site:]
                    sequence_padlock_arm2: str = sequence_oligo[:ligation_site]
                    sequence_padlock_accessory1: str = "TCCTCTATGATTACTGAC"
                    sequence_padlock_ISS_anchor: str = "TGCGTCTATTTAGTGGAGCC"
                    sequence_padlock_accessory2: str = "CTATCTTCTTT"
                    sequence_padlock_backbone: str = (
                        sequence_padlock_accessory1
                        + sequence_padlock_ISS_anchor
                        + barcode
                        + sequence_padlock_accessory2
                    )
                    sequence_padlock_probe: str = (
                        sequence_padlock_arm1 + sequence_padlock_backbone + sequence_padlock_arm2
                    )
                    Tm_arm1 = calc_tm_nn(
                        sequence=sequence_padlock_arm1,
                        Tm_parameters=target_probe_Tm_parameters,
                        Tm_chem_correction_parameters=target_probe_Tm_chem_correction_parameters,
                        Tm_salt_correction_parameters=target_probe_Tm_salt_correction_parameters,
                    )
                    Tm_arm2 = calc_tm_nn(
                        sequence=sequence_padlock_arm2,
                        Tm_parameters=target_probe_Tm_parameters,
                        Tm_chem_correction_parameters=target_probe_Tm_chem_correction_parameters,
                        Tm_salt_correction_parameters=target_probe_Tm_salt_correction_parameters,
                    )

                    new_oligo_properties[oligo_id] = {
                        "barcode": barcode,
                        "sequence_target": oligo_database.get_oligo_property_value(
                            property="target", region_id=region_id, oligo_id=oligo_id, flatten=True
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

            oligo_database.update_oligo_properties(new_oligo_properties)

        return oligo_database

    def generate_output(
        self,
        probe_database: OligoDatabase,
        output_properties: list[str] | None = None,
    ) -> None:
        """
        Generate the final output files for the SCRINSHOT probe design pipeline.

        This method writes all output files required for the SCRINSHOT experiment, including padlock
        probe sequences, detection oligo sequences, and probe properties in multiple formats. The
        output files are written to the pipeline's output directory.

        **Generated Output Files:**

        1. **padlock_probes.yml**: Complete probe information in YAML format, including all specified
           properties for each probe set per region.

        2. **padlock_probes.tsv**: Complete probe information in TSV format, including all specified
           properties for each probe set per region.

        3. **padlock_probes.xlsx**: Complete probe information in Excel format with one sheet per region.
           Each sheet contains probe sets for that region with all specified properties.

        4. **padlock_probes_order.yml**: Simplified YAML file containing only the essential sequences
           needed for ordering probes (padlock probe and detection oligo sequences).

        :param probe_database: The `OligoDatabase` instance containing the final padlock probes
            with all sequences and properties. This should be the result of the `design_padlock_backbone`
            and `design_detection_oligos` methods.
        :type probe_database: OligoDatabase
        :param output_properties: List of property names to include in the output files. If None, a default
            set of properties will be included. Available properties include: 'source', 'species', 'gene_id',
            'chromosome', 'start', 'end', 'strand', 'sequence_target', 'sequence_padlock_arm1',
            'sequence_padlock_arm2', 'sequence_padlock_backbone', 'sequence_padlock_probe',
            'sequence_detection_oligo', 'barcode', 'ligation_site', 'Tm_arm1', 'Tm_arm2', 'Tm_diff_arms',
            'Tm_detection_oligo', 'GC_content_oligo', 'TmNN_oligo', 'isoform_consensus', etc.
        :type output_properties: list[str] | None
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
                "sequence_padlock_probe",
                "sequence_detection_oligo",
                "sequence_padlock_arm1",
                "sequence_padlock_accessory1",
                "sequence_padlock_ISS_anchor",
                "barcode",
                "sequence_padlock_accessory2",
                "sequence_padlock_arm2",
                "sequence_target",
                "GC_content_oligo",
                "TmNN_oligo",
                "ligation_site",
                "Tm_arm1",
                "Tm_arm2",
                "Tm_diff_arms",
                "Tm_detection_oligo",
                "isoform_consensus",
            ]

        probe_database.write_oligosets_to_yaml(
            properties=output_properties,
            ascending=True,
            filename="padlock_probes",
        )

        probe_database.write_oligosets_to_table(
            properties=output_properties,
            ascending=True,
            filename="padlock_probes",
        )

        probe_database.write_ready_to_order_yaml(
            properties=[
                "sequence_padlock_probe",
                "sequence_detection_oligo",
            ],
            ascending=True,
            filename="padlock_probes_order",
        )


############################################
# Scrinshot Target Probe Designer
############################################
class TargetProbeDesigner:
    """
    A class for designing target probes (padlock probe arms) for SCRINSHOT experiments through a multi-step pipeline.

    This class provides methods for the complete target probe design process, which includes:
    1. Creating an initial oligo database from input FASTA files using a sliding window approach
    2. Filtering probes based on sequence properties (GC content, melting temperature, homopolymeric
       runs, detection oligo requirements, padlock arm requirements)
    3. Filtering probes based on specificity to remove off-target binding and cross-hybridization
       using BLASTN searches, with junction-based filtering around the ligation region
    4. Organizing filtered probes into optimal sets based on weighted scoring criteria (isoform
       consensus, GC content, melting temperature) and distance constraints

    The resulting probes are gene-specific targeting sequences (typically 40-45 nt) that will be split
    into padlock probe arms. Each probe is split into two arms (5' and 3') that flank the ligation site,
    and these arms will later be combined with a composite backbone to create complete padlock probes.
    The probes must also support detection oligo design centered on the ligation site with sufficient
    thymines for UNG cleavage in sequential hybridization cycles.

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
        "oligo" (reverse complement). These sequences will later be split into padlock arms
        and used to create complete padlock probes.

        :param region_ids: List of gene identifiers (e.g., gene IDs) to target for probe design. If None,
            all genes present in the input FASTA files will be used.
        :type region_ids: list[str] | None
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

        ##### creating the probe sequences #####
        oligo_sequences = OligoSequenceGenerator(dir_output=self.dir_output)
        oligo_fasta_file = oligo_sequences.create_sequences_sliding_window(
            files_fasta_in=files_fasta_oligo_database,
            length_interval_sequences=(oligo_length_min, oligo_length_max),
            region_ids=region_ids,
            n_jobs=self.n_jobs,
        )

        ##### creating the probe database #####
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

    @pipeline_step_basic(step_name="Property Filters")
    def filter_by_property(
        self,
        oligo_database: OligoDatabase,
        GC_content_min: float,
        GC_content_max: float,
        Tm_min: float,
        Tm_max: float,
        Tm_parameters: dict,
        Tm_chem_correction_parameters: dict | None,
        Tm_salt_correction_parameters: dict | None,
        homopolymeric_base_n: dict[str, int],
        detect_oligo_length_min: int,
        detect_oligo_length_max: int,
        min_thymines: int,
        arm_length_min: int,
        arm_Tm_dif_max: int,
        arm_Tm_min: float,
        arm_Tm_max: float,
    ) -> OligoDatabase:
        """
        Filter the oligo database based on sequence properties to remove probes with undesirable
        characteristics.

        This method applies sequential filtering using multiple property-based filters:
        1. **Hard masked sequences**: Removes probes containing hard-masked nucleotides (N)
        2. **Soft masked sequences**: Removes probes containing soft-masked nucleotides (lowercase)
        3. **Homopolymeric runs**: Removes probes with homopolymeric runs exceeding specified lengths
        4. **GC content**: Removes probes with GC content outside the specified range
        5. **Melting temperature**: Removes probes with Tm outside the specified range
        6. **Detection oligo requirements**: Removes probes that cannot form valid detection oligos
           centered on the ligation site with sufficient thymines for UNG cleavage
        7. **Padlock arm requirements**: Removes probes that cannot be split into valid padlock arms
           with balanced melting temperatures

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
        :param Tm_min: Minimum acceptable melting temperature (Tm) for oligos in degrees Celsius.
            Probes with calculated Tm below this value will be filtered out.
        :type Tm_min: float
        :param Tm_max: Maximum acceptable melting temperature (Tm) for oligos in degrees Celsius.
            Probes with calculated Tm above this value will be filtered out.
        :type Tm_max: float
        :param Tm_parameters: Dictionary of parameters for calculating melting temperature (Tm) using
            the nearest-neighbor method. For using Bio.SeqUtils.MeltingTemp default parameters, set to ``{}``.
            Common parameters include: 'nn_table', 'tmm_table', 'imm_table', 'de_table', 'dnac1', 'dnac2', 'Na', 'K',
            'Tris', 'Mg', 'dNTPs', 'saltcorr', etc. For more information on parameters, see:
            https://biopython.org/docs/1.75/api/Bio.SeqUtils.MeltingTemp.html#Bio.SeqUtils.MeltingTemp.Tm_NN
        :type Tm_parameters: dict
        :param Tm_chem_correction_parameters: Dictionary of chemical correction parameters for Tm
            calculation. These parameters account for the effects of chemical additives (e.g., DMSO,
            formamide) on melting temperature. Set to ``None`` to disable chemical correction, or set to ``{}``
            to use Bio.SeqUtils.MeltingTemp default parameters. For more information, see:
            https://biopython.org/docs/1.75/api/Bio.SeqUtils.MeltingTemp.html#Bio.SeqUtils.MeltingTemp.chem_correction
        :type Tm_chem_correction_parameters: dict | None
        :param Tm_salt_correction_parameters: Dictionary of salt correction parameters for Tm calculation.
            These parameters account for the effects of salt concentration on melting temperature. Set to ``None``
            to disable salt correction, or set to ``{}`` to use Bio.SeqUtils.MeltingTemp default parameters.
            For more information, see:
            https://biopython.org/docs/1.75/api/Bio.SeqUtils.MeltingTemp.html#Bio.SeqUtils.MeltingTemp.salt_correction
        :type Tm_salt_correction_parameters: dict | None
        :param homopolymeric_base_n: Dictionary specifying the maximum allowed length of homopolymeric
            runs for each nucleotide base. Keys should be 'A', 'T', 'G', 'C' and values are the maximum
            run length. For example: {'A': 3, 'T': 3, 'G': 3, 'C': 3} allows up to 3 consecutive
            identical bases.
        :type homopolymeric_base_n: dict[str, int]
        :param detect_oligo_length_min: Minimum length (in nucleotides) for detection oligo sequences
            that will be extracted from the probe, centered on the ligation site.
        :type detect_oligo_length_min: int
        :param detect_oligo_length_max: Maximum length (in nucleotides) for detection oligo sequences
            that will be extracted from the probe, centered on the ligation site.
        :type detect_oligo_length_max: int
        :param min_thymines: Minimum number of thymine (T) nucleotides required in the detection oligo
            region. These thymines will be converted to uracils (U) for UNG cleavage in sequential
            hybridization cycles.
        :type min_thymines: int
        :param arm_length_min: Minimum length (in nucleotides) for each padlock probe arm. Each arm
            must meet this requirement for the probe to pass filtering.
        :type arm_length_min: int
        :param arm_Tm_dif_max: Maximum allowed difference in melting temperature (in degrees Celsius)
            between the two padlock arms. This ensures balanced binding of both arms for efficient
            ligation.
        :type arm_Tm_dif_max: int
        :param arm_Tm_min: Minimum acceptable melting temperature (Tm) for padlock arms in degrees Celsius.
        :type arm_Tm_min: float
        :param arm_Tm_max: Maximum acceptable melting temperature (Tm) for padlock arms in degrees Celsius.
        :type arm_Tm_max: float
        :return: A filtered `OligoDatabase` object containing only probes that pass all property filters.
            Regions with insufficient oligos after filtering are removed.
        :rtype: OligoDatabase
        """

        # define the filters
        hard_masked_sequences = HardMaskedSequenceFilter()
        soft_masked_sequences = SoftMaskedSequenceFilter()
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
            Tm_chem_correction_parameters=Tm_chem_correction_parameters,
            Tm_salt_correction_parameters=Tm_salt_correction_parameters,
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
        oligo_database.remove_regions_with_insufficient_oligos(pipeline_step="Property Filters")
        check_content_oligo_database(oligo_database)

        return oligo_database

    @pipeline_step_basic(step_name="Specificity Filters")
    def filter_by_specificity(
        self,
        oligo_database: OligoDatabase,
        files_fasta_reference_database: list[str],
        specificity_blastn_search_parameters: dict,
        specificity_blastn_hit_parameters: dict,
        cross_hybridization_blastn_search_parameters: dict,
        cross_hybridization_blastn_hit_parameters: dict,
        ligation_region_size: int,
        arm_length_min: int,
        arm_Tm_dif_max: int,
        arm_Tm_min: float,
        arm_Tm_max: float,
        Tm_parameters: dict,
        Tm_chem_correction_parameters: dict | None,
        Tm_salt_correction_parameters: dict | None,
    ) -> OligoDatabase:
        """
        Filter the oligo database based on sequence specificity to remove probes that bind
        non-specifically or cross-hybridize.

        This method applies two types of specificity filters:

        1. **Specificity filtering**: Removes probes that bind to unintended genomic regions
           - **Exact matches**: Removes all probes with exact sequence matches to probes of other regions
           - **BLASTN specificity**: Uses BLASTN to search for similar sequences in the reference database.
             Probes with hits meeting the specified criteria are removed. If `ligation_region_size > 0`,
             seed-based filtering is applied around the ligation site, removing all probes where BLASTN
             hits cover the junction region, regardless of the coverage threshold. If `ligation_region_size == 0`,
             full-length specificity filtering is performed.

        2. **Cross-hybridization filtering**: Removes probes that cross-hybridize with each other.
           This is critical because if probes can bind to each other, they may form dimers instead
           of binding to the target RNA. Probes from the larger genomic region are removed when
           cross-hybridization is detected.

        Before applying specificity filters, the method calculates padlock arm properties (arm sequences,
        ligation site, and arm melting temperatures) for all probes. This information is required for
        seed-based specificity filtering when `ligation_region_size > 0`.

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
        :param ligation_region_size: Size of the ligation region (in nucleotides) around the ligation
            site. If > 0, seed-based specificity filtering is applied: all probes where BLASTN hits
            cover the junction region are removed, regardless of the coverage threshold. If 0,
            full-length specificity filtering is performed. Both modes perform full BLASTN searches.
        :type ligation_region_size: int
        :param arm_length_min: Minimum length (in nucleotides) for each padlock probe arm. Used for
            calculating padlock arm properties before specificity filtering.
        :type arm_length_min: int
        :param arm_Tm_dif_max: Maximum allowed difference in melting temperature (in degrees Celsius)
            between the two padlock arms. Used for calculating padlock arm properties.
        :type arm_Tm_dif_max: int
        :param arm_Tm_min: Minimum acceptable melting temperature (Tm) for padlock arms in degrees Celsius.
            Used for calculating padlock arm properties.
        :type arm_Tm_min: float
        :param arm_Tm_max: Maximum acceptable melting temperature (Tm) for padlock arms in degrees Celsius.
            Used for calculating padlock arm properties.
        :type arm_Tm_max: float
        :param Tm_parameters: Dictionary of parameters for calculating melting temperature (Tm) using
            the nearest-neighbor method. Used for calculating padlock arm Tm values. For using
            Bio.SeqUtils.MeltingTemp default parameters, set to ``{}``. Common parameters include:
            'nn_table', 'tmm_table', 'imm_table', 'de_table', 'dnac1', 'dnac2', 'Na', 'K', 'Tris', 'Mg',
            'dNTPs', 'saltcorr', etc. For more information on parameters, see:
            https://biopython.org/docs/1.75/api/Bio.SeqUtils.MeltingTemp.html#Bio.SeqUtils.MeltingTemp.Tm_NN
        :type Tm_parameters: dict
        :param Tm_chem_correction_parameters: Dictionary of chemical correction parameters for Tm
            calculation. These parameters account for the effects of chemical additives (e.g., DMSO,
            formamide) on melting temperature. Set to ``None`` to disable chemical correction, or set to ``{}``
            to use Bio.SeqUtils.MeltingTemp default parameters. For more information, see:
            https://biopython.org/docs/1.75/api/Bio.SeqUtils.MeltingTemp.html#Bio.SeqUtils.MeltingTemp.chem_correction
        :type Tm_chem_correction_parameters: dict | None
        :param Tm_salt_correction_parameters: Dictionary of salt correction parameters for Tm calculation.
            These parameters account for the effects of salt concentration on melting temperature. Set to ``None``
            to disable salt correction, or set to ``{}`` to use Bio.SeqUtils.MeltingTemp default parameters.
            For more information, see:
            https://biopython.org/docs/1.75/api/Bio.SeqUtils.MeltingTemp.html#Bio.SeqUtils.MeltingTemp.salt_correction
        :type Tm_salt_correction_parameters: dict | None
        :return: A filtered `OligoDatabase` object containing only probes that pass all specificity
            and cross-hybridization filters. The database includes calculated padlock arm properties
            (ligation_site, sequence_padlock_arm1, sequence_padlock_arm2, etc.). Regions with insufficient
            oligos after filtering are removed.
        :rtype: OligoDatabase
        """

        ##### exact match filter #####
        # removing duplicated probes from the region with the most probes
        # exectute seperately before specificity filter to compute ligation side for less oligos
        exact_matches = ExactMatchFilter(policy=RemoveAllFilterPolicy(), filter_name="exact_match")
        specificity_filter = SpecificityFilter(filters=[exact_matches])
        oligo_database = specificity_filter.apply(
            oligo_database=oligo_database,
            sequence_type="oligo",
            n_jobs=self.n_jobs,
        )

        ##### calculate required probe properties #####
        # Calculate padlock arms and detection oligo
        padlock_arms_property = PadlockArmsProperty(
            arm_length_min=arm_length_min,
            arm_Tm_dif_max=arm_Tm_dif_max,
            arm_Tm_min=arm_Tm_min,
            arm_Tm_max=arm_Tm_max,
            Tm_parameters=Tm_parameters,
            Tm_chem_correction_parameters=Tm_chem_correction_parameters,
            Tm_salt_correction_parameters=Tm_salt_correction_parameters,
        )
        calculator = PropertyCalculator(properties=[padlock_arms_property])
        oligo_database = calculator.apply(
            oligo_database=oligo_database, sequence_type="oligo", n_jobs=self.n_jobs
        )

        ##### define reference database #####
        reference_database = ReferenceDatabase(
            database_name=self.subdir_db_reference, dir_output=self.dir_output
        )
        reference_database.load_database_from_file(
            files=files_fasta_reference_database, file_type="fasta", database_overwrite=True
        )

        ##### specificity filters #####
        cross_hybridization_aligner = BlastNFilter(
            search_parameters=cross_hybridization_blastn_search_parameters,
            hit_parameters=cross_hybridization_blastn_hit_parameters,
            filter_name="blastn_crosshybridization",
            dir_output=self.dir_output,
        )
        cross_hybridization_aligner.set_reference_database(reference_database=reference_database)
        cross_hybridization = CrossHybridizationFilter(
            policy=RemoveByLargerRegionFilterPolicy(),
            alignment_method=cross_hybridization_aligner,
            filter_name="blastn_crosshybridization",
            dir_output=self.dir_output,
        )

        specificity: AlignmentSpecificityFilter
        if ligation_region_size > 0:
            specificity = BlastNSeedregionSiteFilter(
                seedregion_size=ligation_region_size,
                seedregion_site_name="ligation_site",
                search_parameters=specificity_blastn_search_parameters,
                hit_parameters=specificity_blastn_hit_parameters,
                filter_name="blastn_specificity",
                dir_output=self.dir_output,
            )
            specificity.set_reference_database(reference_database=reference_database)
        else:
            specificity = BlastNFilter(
                search_parameters=specificity_blastn_search_parameters,
                hit_parameters=specificity_blastn_hit_parameters,
                filter_name="blastn_specificity",
                dir_output=self.dir_output,
            )
            specificity.set_reference_database(reference_database=reference_database)

        specificity_filter = SpecificityFilter(filters=[specificity, cross_hybridization])
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

        This method performs the following steps:
        1. **Scoring**: Calculates scores for each oligo based on weighted criteria (isoform consensus,
           GC content, melting temperature). Higher scores indicate better probes.
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
            passed all previous filtering steps, including padlock arm property calculations.
        :type oligo_database: OligoDatabase
        :param isoform_weight: Weight assigned to isoform consensus in the scoring function.
        :type isoform_weight: float
        :param GC_content_min: Minimum acceptable GC content for oligos, expressed as a fraction
            between 0.0 and 1.0. Used in scoring to penalize probes with GC content below this value.
        :type GC_content_min: float
        :param GC_content_opt: Optimal GC content for oligos, expressed as a fraction between 0.0
            and 1.0. Used in scoring to prioritize probes closer to this value.
        :type GC_content_opt: float
        :param GC_content_max: Maximum acceptable GC content for oligos, expressed as a fraction
            between 0.0 and 1.0. Used in scoring to penalize probes with GC content above this value.
        :type GC_content_max: float
        :param GC_weight: Weight assigned to GC content in the scoring function.
        :type GC_weight: float
        :param Tm_min: Minimum acceptable melting temperature (Tm) for oligos in degrees Celsius.
            Used in scoring to penalize probes with Tm below this value.
        :type Tm_min: float
        :param Tm_opt: Optimal melting temperature (Tm) for oligos in degrees Celsius. Used in scoring
            to prioritize probes closer to this value.
        :type Tm_opt: float
        :param Tm_max: Maximum acceptable melting temperature (Tm) for oligos in degrees Celsius.
            Used in scoring to penalize probes with Tm above this value.
        :type Tm_max: float
        :param Tm_weight: Weight assigned to melting temperature in the scoring function.
        :type Tm_weight: float
        :param Tm_parameters: Dictionary of parameters for calculating melting temperature (Tm) using
            the nearest-neighbor method. For using Bio.SeqUtils.MeltingTemp default parameters, set to ``{}``.
            Common parameters include: 'nn_table', 'tmm_table', 'imm_table', 'de_table', 'dnac1', 'dnac2', 'Na', 'K',
            'Tris', 'Mg', 'dNTPs', 'saltcorr', etc. For more information on parameters, see:
            https://biopython.org/docs/1.75/api/Bio.SeqUtils.MeltingTemp.html#Bio.SeqUtils.MeltingTemp.Tm_NN
        :type Tm_parameters: dict
        :param Tm_chem_correction_parameters: Dictionary of chemical correction parameters for Tm
            calculation. These parameters account for the effects of chemical additives (e.g., DMSO,
            formamide) on melting temperature. Set to ``None`` to disable chemical correction, or set to ``{}``
            to use Bio.SeqUtils.MeltingTemp default parameters. For more information, see:
            https://biopython.org/docs/1.75/api/Bio.SeqUtils.MeltingTemp.html#Bio.SeqUtils.MeltingTemp.chem_correction
        :type Tm_chem_correction_parameters: dict | None
        :param Tm_salt_correction_parameters: Dictionary of salt correction parameters for Tm calculation.
            These parameters account for the effects of salt concentration on melting temperature. Set to ``None``
            to disable salt correction, or set to ``{}`` to use Bio.SeqUtils.MeltingTemp default parameters.
            For more information, see:
            https://biopython.org/docs/1.75/api/Bio.SeqUtils.MeltingTemp.html#Bio.SeqUtils.MeltingTemp.salt_correction
        :type Tm_salt_correction_parameters: dict | None
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
        isoform_consensus_scorer = IsoformConsensusScorer(
            score_weight=isoform_weight,
            property_name_transcript_id="transcript_id",
            property_name_number_total_transcripts="number_total_transcripts",
        )
        Tm_scorer = NormalizedDeviationFromOptimalTmScorer(
            Tm_min=Tm_min,
            Tm_opt=Tm_opt,
            Tm_max=Tm_max,
            Tm_parameters=Tm_parameters,
            Tm_chem_correction_parameters=Tm_chem_correction_parameters,
            Tm_salt_correction_parameters=Tm_salt_correction_parameters,
            score_weight=Tm_weight,
        )
        GC_scorer = NormalizedDeviationFromOptimalGCContentScorer(
            GC_content_min=GC_content_min,
            GC_content_opt=GC_content_opt,
            GC_content_max=GC_content_max,
            score_weight=GC_weight,
        )
        oligos_scoring = OligoScoring(scorers=[isoform_consensus_scorer, Tm_scorer, GC_scorer])
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
# Scrinshot Detection Oligo Designer
############################################
class DetectionOligoDesigner:
    """
    A class for designing detection oligonucleotides for SCRINSHOT padlock probes.

    This class provides methods for generating detection oligos that hybridize to rolling circle
    amplification (RCA) products generated from padlock probes. The design process includes:
    1. Extracting candidate detection oligo sequences centered on the ligation site from target probes
    2. Evaluating multiple candidate sequences (even-length, left-extended, right-extended) to find
       optimal melting temperature
    3. Selecting the candidate with melting temperature closest to the target value
    4. Iteratively optimizing the sequence length to achieve the best Tm match
    5. Converting thymines (T) to uracils (U) at strategic positions for UNG cleavage
    6. Adding fluorophore label position indicators

    Detection oligos are designed with specific requirements:
    - The ligation site must be centered within the oligo sequence (~30-35 nt)
    - Target melting temperature around 56 °C for optimal hybridization to RCA products
    - 2-3 uracil (U) substitutions spaced ≤ 10 nt apart to allow enzymatic cleavage by
      Uracil DNA Glycosylase (UNG), enabling sequential hybridization cycles for multiplexing
    - Minimum number of thymines required to ensure sufficient uracil conversion sites

    The detection oligos are complementary to the gene-specific region of RCA products and are
    labeled at the 3' end with fluorophores (FITC, Cy3, Cy5, etc.) for fluorescence detection.

    :param n_jobs: Number of parallel jobs to use for processing. Set to 1 for serial processing or higher
        values for parallel processing. This affects the parallelization of detection oligo design
        across regions.
    :type n_jobs: int
    """

    def __init__(self, n_jobs: int) -> None:
        """Constructor for the DetectionOligoDesigner class."""

        ##### create the output folder #####
        self.n_jobs = n_jobs

    def create_detection_oligos(
        self,
        oligo_database: OligoDatabase,
        oligo_length_min: int,
        oligo_length_max: int,
        min_thymines: int,
        U_distance: int,
        Tm_opt: float,
        Tm_parameters: dict,
        Tm_chem_correction_parameters: dict | None,
        Tm_salt_correction_parameters: dict | None,
    ) -> OligoDatabase:
        """
        Design detection oligonucleotides for SCRINSHOT padlock probes.

        This method generates detection oligos that are complementary to the gene-specific region of
        rolling circle amplification (RCA) products. For each probe in the database, the method:
        1. Extracts candidate detection oligo sequences centered on the ligation site
        2. Evaluates multiple candidate sequences (even-length, left-extended, right-extended)
        3. Selects the candidate with melting temperature closest to the optimal value
        4. Iteratively shortens the selected candidate from both ends to find the best Tm match
        5. Converts thymines (T) to uracils (U) at appropriate positions for UNG cleavage
        6. Adds fluorophore label position indicator

        Detection oligos are designed with:
        - The ligation site centered within the oligo sequence
        - A target melting temperature around 56 °C for optimal hybridization
        - 2-3 uracil (U) substitutions spaced ≤ specified distance apart to allow enzymatic
          cleavage by Uracil DNA Glycosylase (UNG), facilitating sequential hybridization cycles

        :param oligo_database: The `OligoDatabase` instance containing target probes with their
            sequences, ligation sites, and properties. This database should contain target probes
            organized by region IDs, with each region having one or more probe sets and ligation
            site information calculated from padlock arm properties.
        :type oligo_database: OligoDatabase
        :param oligo_length_min: Minimum length (in nucleotides) for detection oligo sequences.
        :type oligo_length_min: int
        :param oligo_length_max: Maximum length (in nucleotides) for detection oligo sequences.
        :type oligo_length_max: int
        :param min_thymines: Minimum number of thymine (T) nucleotides required in the detection
            oligo sequence. These thymines will be converted to uracils (U) for UNG cleavage.
        :type min_thymines: int
        :param U_distance: Maximum distance (in nucleotides) allowed between uracil substitutions
            in the detection oligo. Uracils must be spaced ≤ this distance apart.
        :type U_distance: int
        :param Tm_opt: Optimal melting temperature (Tm) for detection oligos in degrees Celsius.
            The algorithm will select detection oligos with Tm closest to this value.
        :type Tm_opt: float
        :param Tm_parameters: Dictionary of parameters for calculating melting temperature (Tm) of detection
            oligos using the nearest-neighbor method. For using Bio.SeqUtils.MeltingTemp default parameters, set to ``{}``.
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
        :return: An updated `OligoDatabase` object containing the designed detection oligos. The
            database includes the following new sequence properties for each probe:
            - `sequence_detection_oligo`: The detection oligo sequence with uracil substitutions
              and fluorophore position indicator
            - `Tm_detection_oligo`: The melting temperature of the detection oligo
        :rtype: OligoDatabase
        """

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
        Tm_chem_correction_parameters: dict | None,
        Tm_salt_correction_parameters: dict | None,
    ) -> None:
        """
        Create detection oligos for all probes in a single region.

        This is a helper method called by `create_detection_oligos` to process one region at a time.
        For each probe in the region, it extracts candidate detection oligo sequences centered on
        the ligation site, selects the best candidate based on melting temperature, and converts
        thymines to uracils for UNG cleavage.

        :param oligo_database: The `OligoDatabase` instance containing target probes. This will be
            updated in-place with detection oligo properties.
        :type oligo_database: OligoDatabase
        :param region_id: The identifier of the region to process.
        :type region_id: str
        :param oligo_length_min: Minimum length (in nucleotides) for detection oligo sequences.
        :type oligo_length_min: int
        :param oligo_length_max: Maximum length (in nucleotides) for detection oligo sequences.
        :type oligo_length_max: int
        :param min_thymines: Minimum number of thymine (T) nucleotides required in the detection
            oligo sequence.
        :type min_thymines: int
        :param U_distance: Maximum distance (in nucleotides) allowed between uracil substitutions.
        :type U_distance: int
        :param Tm_opt: Optimal melting temperature (Tm) for detection oligos in degrees Celsius.
        :type Tm_opt: float
        :param Tm_parameters: Dictionary of parameters for calculating melting temperature (Tm).
        :type Tm_parameters: dict
        :param Tm_chem_correction_parameters: Dictionary of chemical correction parameters for Tm
            calculation, or None to disable.
        :type Tm_chem_correction_parameters: dict | None
        :param Tm_salt_correction_parameters: Dictionary of salt correction parameters for Tm
            calculation, or None to disable.
        :type Tm_salt_correction_parameters: dict | None
        :return: None. The oligo_database is updated in-place with detection oligo properties.
        """

        oligosets_region = oligo_database.oligosets[region_id]
        oligosets_oligo_columns = [col for col in oligosets_region.columns if col.startswith("oligo_")]

        new_oligo_properties = {}

        for index in range(len(oligosets_region.index)):
            for column in oligosets_oligo_columns:
                oligo_id = str(oligosets_region.loc[index, column])

                ligation_site = oligo_database.get_oligo_property_value(
                    property="ligation_site", region_id=region_id, oligo_id=oligo_id, flatten=True
                )
                sequence_oligo = oligo_database.get_oligo_property_value(
                    property="oligo", region_id=region_id, oligo_id=oligo_id, flatten=True
                )
                # required for type linting since get_oligo_property_value() could return None
                if not isinstance(sequence_oligo, str) or not isinstance(ligation_site, int):
                    continue

                (
                    detect_oligo_even,
                    detect_oligo_long_left,
                    detect_oligo_long_right,
                ) = calc_detect_oligo(
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

                Tm_detection_oligo = calc_tm_nn(
                    sequence=detection_oligo,
                    Tm_parameters=Tm_parameters,
                    Tm_chem_correction_parameters=Tm_chem_correction_parameters,
                    Tm_salt_correction_parameters=Tm_salt_correction_parameters,
                )

                # exchange T's with U (for enzymatic degradation of oligos)
                detection_oligo = self._exchange_T_with_U(detection_oligo, min_thymines, U_distance)

                new_oligo_properties[oligo_id] = {
                    "Tm_detection_oligo": Tm_detection_oligo,
                    "sequence_detection_oligo": detection_oligo,
                }

        oligo_database.update_oligo_properties(new_oligo_properties)

    def _get_Tm_dif(
        self,
        oligo: str,
        Tm_opt: float,
        Tm_parameters: dict,
        Tm_chem_correction_parameters: dict | None,
        Tm_salt_correction_parameters: dict | None,
    ) -> float:
        """
        Calculate the absolute difference between an oligo's melting temperature and the optimal Tm.

        This helper method is used to evaluate how close a detection oligo's melting temperature
        is to the target optimal value. Lower differences indicate better matches.

        :param oligo: The DNA sequence for which to calculate the Tm difference.
        :type oligo: str
        :param Tm_opt: Optimal melting temperature (Tm) in degrees Celsius.
        :type Tm_opt: float
        :param Tm_parameters: Dictionary of parameters for calculating melting temperature (Tm).
        :type Tm_parameters: dict
        :param Tm_chem_correction_parameters: Dictionary of chemical correction parameters for Tm
            calculation, or None to disable.
        :type Tm_chem_correction_parameters: dict | None
        :param Tm_salt_correction_parameters: Dictionary of salt correction parameters for Tm
            calculation, or None to disable.
        :type Tm_salt_correction_parameters: dict | None
        :return: The absolute difference between the calculated Tm and the optimal Tm, in degrees Celsius.
        :rtype: float
        """

        Tm = calc_tm_nn(
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
        Tm_chem_correction_parameters: dict | None,
        Tm_salt_correction_parameters: dict | None,
    ) -> tuple[list[str], list[float]]:
        """
        Iteratively shorten an oligo sequence to find variants with optimal melting temperature.

        This helper method generates shortened variants of the input oligo by removing nucleotides
        from one end (left or right, depending on `cut_from_right`). It evaluates all variants
        that meet the minimum length and thymine requirements, calculating the Tm difference for
        each to identify the best match to the optimal temperature.

        The method alternates between cutting from the specified end and the opposite end to
        explore a range of sequence lengths while maintaining the ligation site centering.

        :param oligo: The initial DNA sequence to shorten and evaluate.
        :type oligo: str
        :param cut_from_right: If True, start cutting from the right end; if False, start cutting
            from the left end. The method alternates between ends.
        :type cut_from_right: bool
        :param oligo_length_min: Minimum length (in nucleotides) for shortened variants.
            Variants shorter than this will not be generated.
        :type oligo_length_min: int
        :param min_thymines: Minimum number of thymine (T) nucleotides required in each variant.
            Variants with fewer thymines will be skipped.
        :type min_thymines: int
        :param Tm_opt: Optimal melting temperature (Tm) in degrees Celsius.
        :type Tm_opt: float
        :param Tm_parameters: Dictionary of parameters for calculating melting temperature (Tm).
        :type Tm_parameters: dict
        :param Tm_chem_correction_parameters: Dictionary of chemical correction parameters for Tm
            calculation, or None to disable.
        :type Tm_chem_correction_parameters: dict | None
        :param Tm_salt_correction_parameters: Dictionary of salt correction parameters for Tm
            calculation, or None to disable.
        :type Tm_salt_correction_parameters: dict | None
        :return: A tuple containing:
            - **oligos** (list[str]): List of all valid shortened oligo variants
            - **Tm_dif** (list[float]): List of Tm differences (absolute difference from optimal)
              for each variant, in the same order as the oligos list
        :rtype: tuple[list[str], list[float]]
        """

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
        """
        Convert thymine (T) nucleotides to uracil (U) in a detection oligo for UNG cleavage.

        This helper method strategically converts T nucleotides to U to enable enzymatic cleavage
        by Uracil DNA Glycosylase (UNG) in sequential hybridization cycles. The method:
        1. Determines the fluorophore position (left or right) based on T distribution
        2. Converts at least `min_thymines` T nucleotides to U, ensuring they are spaced
           ≤ `U_distance` nucleotides apart
        3. Adds a fluorophore position indicator at the appropriate end

        The uracil substitutions allow the detection oligo to be enzymatically cleaved after each
        imaging round, enabling sequential hybridization cycles for multiplexed detection.

        :param oligo: The detection oligo DNA sequence in which to convert T to U.
        :type oligo: str
        :param min_thymines: Minimum number of thymine (T) nucleotides to convert to uracil (U).
            The method will convert at least this many T nucleotides.
        :type min_thymines: int
        :param U_distance: Maximum distance (in nucleotides) allowed between uracil substitutions.
            Uracils will be spaced ≤ this distance apart to ensure efficient UNG cleavage.
        :type U_distance: int
        :return: The detection oligo sequence with T nucleotides converted to U and a fluorophore
            position indicator added. The indicator "[fluorophore]" is added at the end where
            fewer T nucleotides are present (to preserve more T nucleotides for conversion).
        :rtype: str
        """

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


def main() -> None:
    """
    Main entry point for running the SCRINSHOT probe design pipeline.

    This function orchestrates the complete SCRINSHOT probe design workflow:
    1. Parses command-line arguments using the base parser
    2. Reads the configuration YAML file containing all pipeline parameters
    3. Reads the gene IDs file (if provided) or uses all genes from FASTA files
    4. Preprocesses melting temperature parameters for target probes and detection oligos
    5. Initializes the ScrinshotProbeDesigner pipeline
    6. Designs target probes for specified genes (padlock probe arms)
    7. Designs detection oligos with uracil substitutions for optimal Tm
    8. Assembles padlock probes by combining target probe arms with the composite backbone
    9. Generates output files (YAML, TSV, Excel, order file)

    The function is typically called from the command line:
    ``scrinshot_probe_designer --config <path_to_config.yaml>``

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

    ##### Preprocess Tm parameters #####
    target_probe_Tm_parameters = preprocess_tm_parameters(config["target_probe_Tm_parameters"])

    ##### initialize probe designer pipeline #####
    pipeline = ScrinshotProbeDesigner(
        write_intermediate_steps=config["write_intermediate_steps"],
        dir_output=config["dir_output"],
        n_jobs=config["n_jobs"],
    )

    ##### design probes #####
    oligo_database = pipeline.design_target_probes(
        # Step 1: Create Database Parameters
        region_ids=region_ids,
        files_fasta_target_probe_database=config["files_fasta_target_probe_database"],
        target_probe_length_min=config["target_probe_length_min"],
        target_probe_length_max=config["target_probe_length_max"],
        target_probe_isoform_consensus=config["target_probe_isoform_consensus"],
        # Step 2: Property Filter Parameters
        target_probe_GC_content_min=config["target_probe_GC_content_min"],
        target_probe_GC_content_max=config["target_probe_GC_content_max"],
        target_probe_Tm_min=config["target_probe_Tm_min"],
        target_probe_Tm_max=config["target_probe_Tm_max"],
        target_probe_homopolymeric_base_n=config["target_probe_homopolymeric_base_n"],
        detection_oligo_min_thymines=config["detection_oligo_min_thymines"],
        detection_oligo_length_min=config["detection_oligo_length_min"],
        detection_oligo_length_max=config["detection_oligo_length_max"],
        target_probe_padlock_arm_length_min=config["target_probe_padlock_arm_length_min"],
        target_probe_padlock_arm_Tm_dif_max=config["target_probe_padlock_arm_Tm_dif_max"],
        target_probe_padlock_arm_Tm_min=config["target_probe_padlock_arm_Tm_min"],
        target_probe_padlock_arm_Tm_max=config["target_probe_padlock_arm_Tm_max"],
        target_probe_Tm_parameters=target_probe_Tm_parameters,
        target_probe_Tm_chem_correction_parameters=config["target_probe_Tm_chem_correction_parameters"],
        target_probe_Tm_salt_correction_parameters=config["target_probe_Tm_salt_correction_parameters"],
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
        target_probe_ligation_region_size=config["target_probe_ligation_region_size"],
        # Step 4: Probe Scoring and Set Selection Parameters
        target_probe_isoform_weight=config["target_probe_isoform_weight"],
        target_probe_GC_content_opt=config["target_probe_GC_content_opt"],
        target_probe_GC_weight=config["target_probe_GC_weight"],
        target_probe_Tm_opt=config["target_probe_Tm_opt"],
        target_probe_Tm_weight=config["target_probe_Tm_weight"],
        set_size_min=config["set_size_min"],
        set_size_opt=config["set_size_opt"],
        distance_between_target_probes=config["distance_between_target_probes"],
        n_sets=config["n_sets"],
        n_attempts_graph=config["n_attempts_graph"],
        n_attempts_clique_enum=config["n_attempts_clique_enum"],
        diversification_fraction=config["diversification_fraction"],
        jaccard_opt=config["jaccard_opt"],
        jaccard_step=config["jaccard_step"],
    )

    oligo_database = pipeline.design_detection_oligos(
        oligo_database=oligo_database,
        detection_oligo_length_min=config["detection_oligo_length_min"],
        detection_oligo_length_max=config["detection_oligo_length_max"],
        detection_oligo_min_thymines=config["detection_oligo_min_thymines"],
        detection_oligo_U_distance=config["detection_oligo_U_distance"],
        detection_oligo_Tm_opt=config["detection_oligo_Tm_opt"],
        detection_oligo_Tm_parameters=preprocess_tm_parameters(config["detection_oligo_Tm_parameters"]),
        detection_oligo_Tm_chem_correction_parameters=config["detection_oligo_Tm_chem_correction_parameters"],
        detection_oligo_Tm_salt_correction_parameters=config["detection_oligo_Tm_salt_correction_parameters"],
    )

    probe_database = pipeline.assemble_padlock_backbone(
        oligo_database=oligo_database,
        target_probe_Tm_parameters=target_probe_Tm_parameters,
        target_probe_Tm_chem_correction_parameters=config["target_probe_Tm_chem_correction_parameters"],
        target_probe_Tm_salt_correction_parameters=config["target_probe_Tm_salt_correction_parameters"],
    )

    pipeline.generate_output(probe_database=probe_database)

    logging.info("--------------END PIPELINE--------------")


if __name__ == "__main__":
    main()
