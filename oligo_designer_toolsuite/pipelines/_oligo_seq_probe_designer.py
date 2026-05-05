############################################
# imports
############################################

import logging
import os
import shutil
import warnings
from pathlib import Path

import yaml

from oligo_designer_toolsuite._exceptions import ConfigurationError
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
    BlastNFilter,
    BowtieFilter,
    CrossHybridizationFilter,
    ExactMatchFilter,
    HybridizationProbabilityFilter,
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
        # Step 1: Create Database Parameters
        region_ids: list | None,
        files_fasta_target_probe_database: list,
        target_probe_length_min: int,
        target_probe_length_max: int,
        target_probe_split_region: int,
        target_probe_isoform_consensus: float,
        # Step 2: Property Filter Parameters
        target_probe_GC_content_min: int,
        target_probe_GC_content_max: int,
        target_probe_Tm_min: int,
        target_probe_Tm_max: int,
        target_probe_secondary_structures_T: int,
        target_probe_secondary_structures_threshold_deltaG: int,
        target_probe_homopolymeric_base_n: dict,
        target_probe_prohibited_sequences: list[str],
        target_probe_max_len_selfcomplement: int,
        target_probe_Tm_parameters: dict,
        target_probe_Tm_chem_correction_parameters: dict | None,
        target_probe_Tm_salt_correction_parameters: dict | None,
        # Step 3: Specificity Filter Parameters
        files_fasta_reference_database_target_probe: list,
        files_vcf_reference_database_target_probe: list,
        target_probe_cross_hybridization_alignment_method: str,
        target_probe_cross_hybridization_search_parameters: dict,
        target_probe_cross_hybridization_hit_parameters: dict,
        target_probe_hybridization_probability_alignment_method: str,
        target_probe_hybridization_probability_search_parameters: dict,
        target_probe_hybridization_probability_hit_parameters: dict,
        target_probe_hybridization_probability_threshold: float,
        target_probe_read_length_bias: int,
        # Step 4: Probe Scoring and Set Selection Parameters
        target_probe_targeted_exons: list,
        target_probe_targeted_exons_weight: float,
        target_probe_isoform_weight: float,
        target_probe_GC_content_opt: int,
        target_probe_GC_weight: float,
        target_probe_Tm_opt: int,
        target_probe_Tm_weight: float,
        n_sets: int,
        set_size_min: int,
        set_size_opt: int,
        distance_between_target_probes: int,
        uniform_distance_weight: float,
        n_attempts_graph: int,
        n_attempts_clique_enum: int,
        diversification_fraction: float,
        jaccard_opt: float,
        jaccard_step: float,
    ) -> OligoDatabase:
        """
        Design target probes for Oligo-seq experiments through a multi-step pipeline.

        This method performs the complete target probe design process, which includes:
        1. Creating an initial oligo database from input FASTA files using a sliding window approach
           with region splitting capability
        2. Filtering probes based on sequence properties (GC content, melting temperature, homopolymeric
           runs, secondary structure, self-complementarity)
        3. Filtering probes based on specificity to remove off-target binding, cross-hybridization,
           and variants using BLASTN or Bowtie searches, and hybridization probability filtering
        4. Organizing filtered probes into optimal sets based on weighted scoring criteria (targeted
           exons, isoform consensus, GC content, melting temperature) and distance constraints

        The resulting probes are gene-specific targeting sequences (typically 26-30 nt) that bind to
        RNA transcripts. These probes are used directly for hybridization-based capture in Oligo-seq
        experiments, which combine in situ hybridization with next-generation sequencing readout.

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
        :param target_probe_split_region: Size of regions (in nucleotides) to split large genomic
            regions into when generating probe candidates. This helps manage memory usage for very
            long sequences.
        :type target_probe_split_region: int
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
        :param target_probe_secondary_structures_T: Temperature in degrees Celsius at which to evaluate
            secondary structure formation.
        :type target_probe_secondary_structures_T: int
        :param target_probe_secondary_structures_threshold_deltaG: DeltaG threshold (in kcal/mol) for
            secondary structure stability. Probes with secondary structures having deltaG values more
            negative than this threshold will be filtered out.
        :type target_probe_secondary_structures_threshold_deltaG: int
        :param target_probe_homopolymeric_base_n: Dictionary specifying the maximum allowed length of
            homopolymeric runs for each nucleotide base (keys: 'A', 'T', 'G', 'C').
        :type target_probe_homopolymeric_base_n: dict[str, int]
        :param target_probe_prohibited_sequences: The sequences to prohibit in the oligos. If an oligo contains any of these sequences, it will be filtered out.
        :type target_probe_prohibited_sequences: list[str]
        :param target_probe_max_len_selfcomplement: Maximum allowable length of self-complementary
            sequences. Probes with longer self-complementary regions can form hairpins and reduce
            hybridization efficiency.
        :type target_probe_max_len_selfcomplement: int
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
        :param files_vcf_reference_database_target_probe: List of paths to VCF files containing variant
            information used for filtering probes that overlap with known single nucleotide polymorphisms
            (SNPs) or other variants. Probes overlapping variants may have reduced specificity.
        :type files_vcf_reference_database_target_probe: list[str]
        :param target_probe_cross_hybridization_alignment_method: Alignment method to use for
            cross-hybridization filtering. Must be either "blastn" or "bowtie".
        :type target_probe_cross_hybridization_alignment_method: str
        :param target_probe_cross_hybridization_search_parameters: Dictionary of parameters for alignment
            searches used in cross-hybridization filtering. Parameters depend on the alignment method
            (BLASTN or Bowtie).
        :type target_probe_cross_hybridization_search_parameters: dict
        :param target_probe_cross_hybridization_hit_parameters: Dictionary of parameters for filtering
            alignment hits in cross-hybridization searches. Probes with cross-hybridization hits meeting
            these criteria are removed from the larger region.
        :type target_probe_cross_hybridization_hit_parameters: dict
        :param target_probe_hybridization_probability_alignment_method: Alignment method to use for
            hybridization probability filtering. Must be either "blastn" or "bowtie".
        :type target_probe_hybridization_probability_alignment_method: str
        :param target_probe_hybridization_probability_search_parameters: Dictionary of parameters for
            alignment searches used in hybridization probability filtering. Parameters depend on the
            alignment method (BLASTN or Bowtie).
        :type target_probe_hybridization_probability_search_parameters: dict
        :param target_probe_hybridization_probability_hit_parameters: Dictionary of parameters for
            filtering alignment hits in hybridization probability searches.
        :type target_probe_hybridization_probability_hit_parameters: dict
        :param target_probe_hybridization_probability_threshold: Threshold for hybridization probability
            filtering. Probes with hybridization probabilities above this threshold are removed, as they
            may bind non-specifically.
        :type target_probe_hybridization_probability_threshold: float
        :param target_probe_read_length_bias: Number of nucleotides from the 5' end of probes to check
            for read length bias. Probes where the first N bases match exactly with other probes are
            removed to prevent sequencing read length biases.
        :type target_probe_read_length_bias: int

        **Step 4: Probe Scoring and Set Selection Parameters**

        :param target_probe_targeted_exons: List of exon identifiers that should be preferentially
            targeted by probes. Probes overlapping these exons receive higher scores.
        :type target_probe_targeted_exons: list[str]
        :param target_probe_targeted_exons_weight: Weight assigned to targeted exons overlap in the
            scoring function.
        :type target_probe_targeted_exons_weight: float
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
        :param n_sets: Number of oligo sets to generate per region. Multiple sets allow for redundancy and selection
            of the best-performing set based on scoring criteria.
        :type n_sets: int
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
        :param uniform_distance_weight: Weight assigned to uniform distance in the scoring function.
        :type uniform_distance_weight: float
        :param n_attempts_graph: Number of randomized graph attempts. In each attempt, a fraction of nodes is randomly
            removed from the compatibility graph to create diversity; more attempts increase diversity at the cost of runtime.
        :type n_attempts_graph: int
        :param n_attempts_clique_enum: Maximum number of cliques enumerated per graph attempt. Limits how many cliques
            are explored before stopping enumeration for the current graph.
        :type n_attempts_clique_enum: int
        :param diversification_fraction: Fraction of oligos to remove from the graph to create diversity in the set selection.
        :type diversification_fraction: float
        :param jaccard_opt: Optimal maximum Jaccard overlap between selected sets.
        :type jaccard_opt: float
        :param jaccard_step: Step size used to relax the Jaccard constraint when not enough sets are found.
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
            split_region=target_probe_split_region,
            files_fasta_oligo_database=files_fasta_target_probe_database,
            min_oligos_per_gene=set_size_min,
        )

        if self.write_intermediate_steps:
            dir_database = oligo_database.save_database(name_database="1_db_probes_initial")
            logging.info(f"Saved probe database for step 1 (Create Database) in directory {dir_database}")

        oligo_database = target_probe_designer.filter_by_property(
            oligo_database=oligo_database,
            isoform_consensus=target_probe_isoform_consensus,
            GC_content_min=target_probe_GC_content_min,
            GC_content_max=target_probe_GC_content_max,
            Tm_min=target_probe_Tm_min,
            Tm_max=target_probe_Tm_max,
            Tm_parameters=target_probe_Tm_parameters,
            Tm_chem_correction_parameters=target_probe_Tm_chem_correction_parameters,
            Tm_salt_correction_parameters=target_probe_Tm_salt_correction_parameters,
            homopolymeric_base_n=target_probe_homopolymeric_base_n,
            prohibited_sequences=target_probe_prohibited_sequences,
            max_len_selfcomplement=target_probe_max_len_selfcomplement,
            secondary_structures_T=target_probe_secondary_structures_T,
            secondary_structures_threshold_deltaG=target_probe_secondary_structures_threshold_deltaG,
        )

        if self.write_intermediate_steps:
            dir_database = oligo_database.save_database(name_database="2_db_probes_property_filter")
            logging.info(f"Saved probe database for step 2 (Property Filters) in directory {dir_database}")

        oligo_database = target_probe_designer.filter_by_specificity(
            oligo_database=oligo_database,
            files_fasta_reference_database=files_fasta_reference_database_target_probe,
            files_vcf_reference_database=files_vcf_reference_database_target_probe,
            target_probe_read_length_bias=target_probe_read_length_bias,
            cross_hybridization_alignment_method=target_probe_cross_hybridization_alignment_method,
            cross_hybridization_search_parameters=target_probe_cross_hybridization_search_parameters,
            cross_hybridization_hit_parameters=target_probe_cross_hybridization_hit_parameters,
            hybridization_probability_alignment_method=target_probe_hybridization_probability_alignment_method,
            hybridization_probability_search_parameters=target_probe_hybridization_probability_search_parameters,
            hybridization_probability_hit_parameters=target_probe_hybridization_probability_hit_parameters,
            hybridization_probability_threshold=target_probe_hybridization_probability_threshold,
        )

        if self.write_intermediate_steps:
            dir_database = oligo_database.save_database(name_database="3_db_probes_specificity_filter")
            logging.info(f"Saved probe database for step 3 (Specificity Filters) in directory {dir_database}")

        oligo_database = target_probe_designer.create_oligo_sets(
            oligo_database=oligo_database,
            targeted_exons=target_probe_targeted_exons,
            targeted_exons_weight=target_probe_targeted_exons_weight,
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
            n_sets=n_sets,
            set_size_min=set_size_min,
            set_size_opt=set_size_opt,
            distance_between_oligos=distance_between_target_probes,
            uniform_distance_weight=uniform_distance_weight,
            n_attempts_graph=n_attempts_graph,
            n_attempts_clique_enum=n_attempts_clique_enum,
            diversification_fraction=diversification_fraction,
            jaccard_opt=jaccard_opt,
            jaccard_step=jaccard_step,
        )

        # Calculate oligo length, GC content, Tm, num targeted transcripts, isoform consensus, and length self complement
        length_property = LengthProperty()
        gc_content_property = GCContentProperty()
        TmNN_property = TmNNProperty(
            Tm_parameters=target_probe_Tm_parameters,
            Tm_chem_correction_parameters=target_probe_Tm_chem_correction_parameters,
            Tm_salt_correction_parameters=target_probe_Tm_salt_correction_parameters,
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
        Generate the final output files for the Oligo-seq probe design pipeline.

        This method writes all output files required for the Oligo-seq experiment, including probe
        sequences and properties in multiple formats. The output files are written to the pipeline's
        output directory.

        **Generated Output Files:**

        1. **oligo_seq_probes.yml**: Complete probe information in YAML format, including all specified
           properties for each probe in the top N sets per region.

        2. **oligo_seq_probes.tsv**: Complete probe information in TSV format, including all specified
           properties for each probe in the top N sets per region.

        3. **oligo_seq_probes.xlsx**: Complete probe information in Excel format with one sheet per region.
           Each sheet contains probe sets for that region with all specified properties.

        4. **oligo_seq_probes_order.yml**: Simplified YAML file containing only the essential sequences
           needed for ordering probes (oligo sequences).

        :param oligo_database: The `OligoDatabase` instance containing the final target probes
            with all sequences and properties. This should be the result of the `design_target_probes`
            method.
        :type oligo_database: OligoDatabase
        :param output_properties: List of property names to include in the output files. If None, a default
            set of properties will be included. Available properties include: 'source', 'species', 'gene_id',
            'chromosome', 'start', 'end', 'strand', 'oligo', 'target', 'length_oligo', 'GC_content_oligo',
            'TmNN_oligo', 'SNP_filter', 'isoform_consensus', 'length_selfcomplement_oligo', etc.
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
                "oligo",
                "target",
                "length_oligo",
                "GC_content_oligo",
                "TmNN_oligo",
                "SNP_filter",
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
       variants, and high hybridization probability using BLASTN or Bowtie searches, and read length
       bias filtering
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
        Create an initial oligo database by generating sequences using a sliding window approach
        with region splitting and performing pre-filtering based on isoform consensus.

        This is the first step of target probe design. The method:
        1. Generates candidate oligo sequences from input FASTA files using a sliding window approach
           across the specified length range, with large regions split into smaller chunks for memory efficiency
        2. Creates an `OligoDatabase` and loads the generated sequences
        3. Calculates reverse complement sequences and isoform consensus properties
        4. Pre-filters oligos based on isoform consensus threshold
        5. Removes regions with insufficient oligos after filtering

        The database stores sequences with sequence types "target" (original sequence), "oligo" (reverse
        complement), and "oligo_short" (shortened sequences for read length bias filtering).

        :param region_ids: List of gene identifiers (e.g., gene IDs) to target for probe design. If None,
            all genes present in the input FASTA files will be used.
        :type region_ids: list[str] | None
        :param oligo_length_min: Minimum length (in nucleotides) for target probe sequences.
        :type oligo_length_min: int
        :param oligo_length_max: Maximum length (in nucleotides) for target probe sequences.
        :type oligo_length_max: int
        :param split_region: Size of regions (in nucleotides) to split large genomic regions into when
            generating probe candidates. This helps manage memory usage for very long sequences.
        :type split_region: int
        :param files_fasta_oligo_database: List of paths to FASTA files containing sequences from which
            target probes will be generated. These files should contain genomic regions of interest
            (e.g., exons, exon-exon junctions).
        :type files_fasta_oligo_database: list[str]
        :param min_oligos_per_gene: Minimum number of oligos required per region (gene) after filtering.
            Regions with fewer oligos than this threshold will be removed from the database.
        :type min_oligos_per_gene: int
        :return: An `OligoDatabase` object containing the generated target probe sequences with their
            component sequences (target, oligo, oligo_short) and calculated properties (isoform_consensus).
            The database is filtered to only include regions that meet the minimum oligo requirement.
        :rtype: OligoDatabase
        """

        ##### creating the oligo sequences #####
        oligo_sequences = OligoSequenceGenerator(dir_output=self.dir_output)
        oligo_fasta_file = oligo_sequences.create_sequences_sliding_window(
            files_fasta_in=files_fasta_oligo_database,
            length_interval_sequences=(oligo_length_min, oligo_length_max),
            split_region=split_region,
            region_ids=region_ids,
            n_jobs=self.n_jobs,
        )

        ##### creating the oligo database #####
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
        reverse_complement_sequence_property = ReverseComplementSequenceProperty(
            sequence_type_reverse_complement="oligo"
        )
        calculator = PropertyCalculator(properties=[reverse_complement_sequence_property])
        oligo_database = calculator.apply(
            oligo_database=oligo_database, sequence_type="target", n_jobs=self.n_jobs
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
        isoform_consensus: float,
        GC_content_min: int,
        GC_content_max: int,
        Tm_min: int,
        Tm_max: int,
        Tm_parameters: dict,
        Tm_chem_correction_parameters: dict | None,
        Tm_salt_correction_parameters: dict | None,
        homopolymeric_base_n: dict[str, int],
        prohibited_sequences: list[str],
        max_len_selfcomplement: int,
        secondary_structures_T: float,
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
        5. **Melting temperature**: Removes probes with Tm outside the specified range
        6. **Self-complementarity**: Removes probes with excessive self-complementary regions that can
           form hairpins and reduce hybridization efficiency
        7. **Secondary structure**: Removes probes that form stable secondary structures at the
           specified temperature

        Probes that fail any filter are removed. Regions with insufficient oligos after filtering
        are removed from the database.

        :param oligo_database: The `OligoDatabase` instance containing oligonucleotide sequences
            and their associated properties. This database should contain target probes with their
            component sequences already calculated.
        :type oligo_database: OligoDatabase
        :param isoform_consensus: Threshold for isoform consensus filtering (typically between 0.0 and 1.0).
            Probes with isoform consensus values below this threshold will be filtered out. This ensures
            that selected probes target sequences that are conserved across multiple transcript isoforms.
        :type isoform_consensus: float
        :param GC_content_min: Minimum acceptable GC content for oligos, expressed as a fraction
            between 0.0 and 1.0.
        :type GC_content_min: int
        :param GC_content_max: Maximum acceptable GC content for oligos, expressed as a fraction
            between 0.0 and 1.0.
        :type GC_content_max: int
        :param Tm_min: Minimum acceptable melting temperature (Tm) for oligos in degrees Celsius.
            Probes with calculated Tm below this value will be filtered out.
        :type Tm_min: int
        :param Tm_max: Maximum acceptable melting temperature (Tm) for oligos in degrees Celsius.
            Probes with calculated Tm above this value will be filtered out.
        :type Tm_max: int
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
        :param prohibited_sequences: List of sequences that are prohibited from being used in the probes.
            These sequences will be removed from the database.
        :type prohibited_sequences: list[str]
        :param max_len_selfcomplement: Maximum allowable length of self-complementary sequences.
            Probes with longer self-complementary regions can form hairpins and reduce hybridization
            efficiency.
        :type max_len_selfcomplement: int
        :param secondary_structures_T: Temperature in degrees Celsius at which to evaluate secondary
            structure formation. Secondary structures that form at this temperature can interfere
            with probe binding.
        :type secondary_structures_T: float
        :param secondary_structures_threshold_deltaG: DeltaG threshold (in kcal/mol) for secondary
            structure stability. Probes with secondary structures having deltaG values more negative
            (more stable) than this threshold will be filtered out.
        :type secondary_structures_threshold_deltaG: float
        :return: A filtered `OligoDatabase` object containing only probes that pass all property filters.
            Regions with insufficient oligos after filtering are removed.
        :rtype: OligoDatabase
        """
        ##### pre-filter oligo database for certain properties #####
        isoform_consensus_property = IsoformConsensusProperty()
        calculator = PropertyCalculator(properties=[isoform_consensus_property])
        oligo_database = calculator.apply(
            oligo_database=oligo_database, sequence_type="target", n_jobs=self.n_jobs
        )
        oligo_database.filter_database_by_property_threshold(
            property_name="isoform_consensus",
            property_thr=isoform_consensus,
            remove_if_smaller_threshold=True,
        )

        # define the filters
        hard_masked_sequences = HardMaskedSequenceFilter()
        soft_masked_sequences = SoftMaskedSequenceFilter()
        self_comp = SelfComplementFilter(
            max_len_selfcomplement=max_len_selfcomplement,
        )
        homopolymeric_runs = HomopolymericRunsFilter(
            base_n=homopolymeric_base_n,
        )
        gc_content = GCContentFilter(GC_content_min=GC_content_min, GC_content_max=GC_content_max)
        melting_temperature = MeltingTemperatureNNFilter(
            Tm_min=Tm_min,
            Tm_max=Tm_max,
            Tm_parameters=Tm_parameters,
            Tm_chem_correction_parameters=Tm_chem_correction_parameters,
            Tm_salt_correction_parameters=Tm_salt_correction_parameters,
        )
        secondary_sctructure = SecondaryStructureFilter(
            T=secondary_structures_T,
            thr_DG=secondary_structures_threshold_deltaG,
        )
        prohibited_sequence_filter = ProhibitedSequenceFilter(
            prohibited_sequences=prohibited_sequences,
        )

        filters = [
            hard_masked_sequences,
            soft_masked_sequences,
            prohibited_sequence_filter,
            homopolymeric_runs,
            gc_content,
            self_comp,
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

    @pipeline_step_basic(step_name="Specificity Filters")
    def filter_by_specificity(
        self,
        oligo_database: OligoDatabase,
        files_fasta_reference_database: list[str],
        files_vcf_reference_database: list[str],
        target_probe_read_length_bias: int,
        cross_hybridization_alignment_method: str,
        cross_hybridization_search_parameters: dict,
        cross_hybridization_hit_parameters: dict,
        hybridization_probability_alignment_method: str,
        hybridization_probability_search_parameters: dict,
        hybridization_probability_hit_parameters: dict,
        hybridization_probability_threshold: float,
    ) -> OligoDatabase:
        """
        Filter the oligo database based on sequence specificity to remove probes that bind
        non-specifically, cross-hybridize, overlap with variants, or have high hybridization probability.

        This method applies multiple types of specificity filters:

        1. **Read length bias filtering**: Removes probes where the first N bases match exactly with
           other probes. This prevents sequencing read length biases that could occur when multiple
           probes share identical 5' ends.

        2. **Exact match filtering**: Removes all probes with exact sequence matches to probes of
           other regions.

        3. **Variant filtering**: Marks probes that overlap with known single nucleotide polymorphisms
           (SNPs) or other variants from VCF files. Probes overlapping variants may have reduced
           specificity and are flagged (not removed) for downstream analysis.

        4. **Cross-hybridization filtering**: Removes probes that cross-hybridize with each other.
           This is critical because if probes can bind to each other, they may form dimers instead
           of binding to the target RNA. Probes from the larger genomic region are removed when
           cross-hybridization is detected. The alignment method (BLASTN or Bowtie) can be selected.

        5. **Hybridization probability filtering**: Removes probes with high hybridization probability
           to unintended targets in the reference database. This filter uses alignment-based methods
           (BLASTN or Bowtie) to estimate the probability of non-specific binding and removes probes
           above the specified threshold.

        The reference databases are loaded from the provided FASTA and VCF files. Alignment methods
        (BLASTN or Bowtie) can be selected independently for cross-hybridization and hybridization
        probability filtering. Regions that do not meet the minimum oligo requirement after filtering
        are removed from the database.

        :param oligo_database: The `OligoDatabase` instance containing oligonucleotide sequences
            and their associated properties. This database should contain target probes with their
            component sequences already calculated.
        :type oligo_database: OligoDatabase
        :param files_fasta_reference_database: List of paths to FASTA files containing reference
            sequences against which specificity will be evaluated. These typically include the
            entire genome or transcriptome to identify off-target binding sites.
        :type files_fasta_reference_database: list[str]
        :param files_vcf_reference_database: List of paths to VCF files containing variant
            information used for filtering probes that overlap with known single nucleotide polymorphisms
            (SNPs) or other variants. Probes overlapping variants are flagged (not removed) for
            downstream analysis.
        :type files_vcf_reference_database: list[str]
        :param target_probe_read_length_bias: Number of nucleotides from the 5' end of probes to check
            for read length bias. Probes where the first N bases match exactly with other probes are
            removed to prevent sequencing read length biases.
        :type target_probe_read_length_bias: int
        :param cross_hybridization_alignment_method: Alignment method to use for cross-hybridization
            filtering. Must be either "blastn" or "bowtie".
        :type cross_hybridization_alignment_method: str
        :param cross_hybridization_search_parameters: Dictionary of parameters for alignment searches
            used in cross-hybridization filtering. Parameters depend on the alignment method
            (BLASTN or Bowtie).
        :type cross_hybridization_search_parameters: dict
        :param cross_hybridization_hit_parameters: Dictionary of parameters for filtering alignment
            hits in cross-hybridization searches. Probes with cross-hybridization hits meeting these
            criteria are removed from the larger region.
        :type cross_hybridization_hit_parameters: dict
        :param hybridization_probability_alignment_method: Alignment method to use for hybridization
            probability filtering. Must be either "blastn" or "bowtie".
        :type hybridization_probability_alignment_method: str
        :param hybridization_probability_search_parameters: Dictionary of parameters for alignment
            searches used in hybridization probability filtering. Parameters depend on the alignment
            method (BLASTN or Bowtie).
        :type hybridization_probability_search_parameters: dict
        :param hybridization_probability_hit_parameters: Dictionary of parameters for filtering
            alignment hits in hybridization probability searches.
        :type hybridization_probability_hit_parameters: dict
        :param hybridization_probability_threshold: Threshold for hybridization probability filtering.
            Probes with hybridization probabilities above this threshold are removed, as they may bind
            non-specifically.
        :type hybridization_probability_threshold: float
        :return: A filtered `OligoDatabase` object containing only probes that pass all specificity
            filters. Probes overlapping variants are flagged but not removed. Regions with insufficient
            oligos after filtering are removed.
        :rtype: OligoDatabase
        """

        def _get_alignment_method(
            alignment_method: str,
            search_parameters: dict,
            hit_parameters: dict,
            filter_name: str,
            dir_output: str,
        ) -> BlastNFilter | BowtieFilter:
            if alignment_method == "blastn":
                return BlastNFilter(
                    search_parameters=search_parameters,
                    hit_parameters=hit_parameters,
                    filter_name=filter_name,
                    dir_output=dir_output,
                )
            elif alignment_method == "bowtie":
                return BowtieFilter(
                    search_parameters=search_parameters,
                    filter_name=filter_name,
                    dir_output=dir_output,
                )
            else:
                raise ConfigurationError(
                    f"Alignment method '{alignment_method}' is not supported. "
                    f"Supported methods are: 'blastn' or 'bowtie'."
                )

        ##### define reference database #####
        reference_database_alignment = ReferenceDatabase(
            database_name=f"{self.subdir_db_reference}_sequences", dir_output=self.dir_output
        )
        reference_database_alignment.load_database_from_file(
            files=files_fasta_reference_database, file_type="fasta", database_overwrite=True
        )

        reference_database_variants = ReferenceDatabase(
            database_name=f"{self.subdir_db_reference}_variants", dir_output=self.dir_output
        )
        reference_database_variants.load_database_from_file(
            files=files_vcf_reference_database, file_type="vcf", database_overwrite=True
        )

        ##### define exact match filter #####
        # remove sequences that could cause read length biases because the first
        # <target_probe_read_length_bias> bases of both sequences match
        # Calculate shortened sequence using new PropertyCalculator pattern
        shortened_sequence_property = ShortenedSequenceProperty(
            sequence_length=target_probe_read_length_bias, reverse=False
        )
        calculator = PropertyCalculator(properties=[shortened_sequence_property])
        oligo_database = calculator.apply(
            oligo_database=oligo_database, sequence_type="oligo", n_jobs=self.n_jobs
        )

        exact_matches_short = ExactMatchFilter(
            policy=RemoveAllFilterPolicy(), filter_name="exact_match_read_length_bias"
        )

        ##### run exact match filter #####
        specificity_filter = SpecificityFilter(filters=[exact_matches_short])
        oligo_database = specificity_filter.apply(
            oligo_database=oligo_database,
            sequence_type="oligo_short",
            n_jobs=self.n_jobs,
        )

        ##### define specificity filters #####
        exact_matches = ExactMatchFilter(policy=RemoveAllFilterPolicy(), filter_name="exact_match")

        variants = VariantsFilter(remove_hits=False, filter_name="SNP_filter", dir_output=self.dir_output)
        variants.set_reference_database(reference_database=reference_database_variants)

        cross_hybridization_aligner = _get_alignment_method(
            alignment_method=cross_hybridization_alignment_method,
            search_parameters=cross_hybridization_search_parameters,
            hit_parameters=cross_hybridization_hit_parameters,
            filter_name="cross_hybridization_filter",
            dir_output=self.dir_output,
        )
        cross_hybridization = CrossHybridizationFilter(
            policy=RemoveByLargerRegionFilterPolicy(),
            alignment_method=cross_hybridization_aligner,
            filter_name="cross_hybridization_filter",
            dir_output=self.dir_output,
        )

        hybridization_probability_aligner = _get_alignment_method(
            alignment_method=hybridization_probability_alignment_method,
            search_parameters=hybridization_probability_search_parameters,
            hit_parameters=hybridization_probability_hit_parameters,
            filter_name="hybridization_probability_filter",
            dir_output=self.dir_output,
        )
        hybridization_probability_aligner.set_reference_database(
            reference_database=reference_database_alignment
        )
        hybridization_probability = HybridizationProbabilityFilter(
            alignment_method=hybridization_probability_aligner,
            threshold=hybridization_probability_threshold,
            filter_name="hybridization_probability_filter",
            dir_output=self.dir_output,
        )

        # run all filters specified above
        filters = [
            exact_matches,
            variants,
            cross_hybridization,
            hybridization_probability,
        ]
        specificity_filter = SpecificityFilter(filters=filters)
        oligo_database = specificity_filter.apply(
            oligo_database=oligo_database,
            sequence_type="oligo",
            n_jobs=self.n_jobs,
        )

        # remove all directories of intermediate steps
        for directory in [
            cross_hybridization_aligner.dir_output,
            cross_hybridization.dir_output,
            hybridization_probability_aligner.dir_output,
            hybridization_probability.dir_output,
            variants.dir_output,
        ]:
            if os.path.exists(directory):
                shutil.rmtree(directory)

        oligo_database.remove_regions_with_insufficient_oligos("Specificity Filters")
        check_content_oligo_database(oligo_database)

        return oligo_database

    @pipeline_step_basic(step_name="Oligo Selection")
    def create_oligo_sets(
        self,
        oligo_database: OligoDatabase,
        targeted_exons: list[str],
        targeted_exons_weight: float,
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
        n_sets: int,
        set_size_opt: int,
        set_size_min: int,
        distance_between_oligos: int,
        uniform_distance_weight: float,
        n_attempts_graph: int,
        n_attempts_clique_enum: int,
        diversification_fraction: float,
        jaccard_opt: float,
        jaccard_step: float,
    ) -> OligoDatabase:
        """
        Create optimal oligo sets based on weighted scoring criteria, distance constraints, and set selection.

        This method performs the following steps:
        1. **Scoring**: Calculates scores for each oligo based on weighted criteria (targeted exons overlap,
           isoform consensus, GC content, melting temperature). Higher scores indicate better probes.
        2. **Set generation**: Builds a compatibility graph from distance constraints and selects sets via
           a graph-based (clique) strategy. Generates multiple diverse sets per region, controlling overlap
           between sets using a Jaccard threshold (`jaccard_opt`) with optional relaxation (`jaccard_step`).
        3. **Set scoring**: Evaluates each generated set by average oligo score and selects the best sets
           (ascending order, i.e. lower average score is better for this pipeline).
        4. **Region filtering**: Removes regions that cannot generate sets meeting the minimum size requirement.

        The algorithm attempts to find sets with optimal size (`set_size_opt`) but may produce sets
        as small as `set_size_min` if constraints cannot be met.

        :param oligo_database: The `OligoDatabase` instance containing oligonucleotide sequences
            and their associated properties. This database should contain target probes that have
            passed all previous filtering steps.
        :type oligo_database: OligoDatabase
        :param targeted_exons: List of exon identifiers that should be preferentially targeted by probes.
            Probes overlapping these exons receive higher scores.
        :type targeted_exons: list[str]
        :param targeted_exons_weight: Weight assigned to targeted exons overlap in the scoring function.
        :type targeted_exons_weight: float
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
        :param n_sets: Number of oligo sets to generate per region. Multiple sets allow for redundancy and selection
            of the best-performing set based on scoring criteria.
        :type n_sets: int
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
        :param uniform_distance_weight: Weight assigned to uniform distance in the scoring function.
        :type uniform_distance_weight: float
        :param n_attempts_graph: Number of randomized graph attempts. In each attempt, a fraction of nodes is randomly
            removed from the compatibility graph to create diversity.
        :type n_attempts_graph: int
        :param n_attempts_clique_enum: Maximum number of cliques enumerated per graph attempt.
        :type n_attempts_clique_enum: int
        :param diversification_fraction: Fraction of oligos to remove from the graph to create diversity in the set selection.
        :type diversification_fraction: float
        :param jaccard_opt: Optimal maximum Jaccard overlap between selected sets.
        :type jaccard_opt: float
        :param jaccard_step: Step size used to relax the Jaccard constraint when not enough sets are found.
        :type jaccard_step: float
        :return: An updated `OligoDatabase` object containing the generated oligo sets. Each region
            will have up to `n_sets` sets stored, with each set containing between `set_size_min` and
            `set_size_opt` probes. Regions with insufficient oligos are removed.
        :rtype: OligoDatabase
        """
        oligo_lengths = [
            len(sequence) for sequence in oligo_database.get_sequence_list(sequence_type="oligo")
        ]
        average_oligo_length = sum(oligo_lengths) / len(oligo_lengths)

        # Define all scorers
        uniform_distance_scorer = UniformDistanceScorer(
            average_oligo_length=average_oligo_length, score_weight=uniform_distance_weight
        )
        exon_scorer = OverlapTargetedExonsScorer(
            targeted_exons=targeted_exons,
            score_weight=targeted_exons_weight,
            property_name="exon_number",
        )
        isoform_scorer = IsoformConsensusScorer(score_weight=isoform_weight)
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

        oligos_scoring = OligoScoring(
            scorers=[exon_scorer, isoform_scorer, Tm_scorer, GC_scorer, uniform_distance_scorer]
        )
        set_scoring = AverageSetScoring(ascending=True)

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
# Oligo-seq Designer Pipeline
############################################


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
    pipeline = OligoSeqProbeDesigner(
        write_intermediate_steps=config["write_intermediate_steps"],
        dir_output=config["dir_output"],
        n_jobs=config["n_jobs"],
    )

    ##### Add highly abundant k-mers to prohibited sequences #####
    target_probe_prohibited_sequences = list(
        set(
            config["target_probe_prohibited_sequences"]
            + get_highly_abundant_kmer_sequences(
                files_fasta=config["files_fasta_target_probe_database"],
                kmer_abundance_threshold=config["target_probe_kmer_abundance_threshold"],
            )
        )
    )

    ##### design probes #####
    oligo_database = pipeline.design_target_probes(
        # Step 1: Create Database Parameters
        region_ids=region_ids,
        files_fasta_target_probe_database=config["files_fasta_target_probe_database"],
        target_probe_length_min=config["target_probe_length_min"],
        target_probe_length_max=config["target_probe_length_max"],
        target_probe_split_region=config["target_probe_split_region"],
        target_probe_isoform_consensus=config["target_probe_isoform_consensus"],
        # Step 2: Property Filter Parameters
        target_probe_GC_content_min=config["target_probe_GC_content_min"],
        target_probe_GC_content_max=config["target_probe_GC_content_max"],
        target_probe_Tm_min=config["target_probe_Tm_min"],
        target_probe_Tm_max=config["target_probe_Tm_max"],
        target_probe_secondary_structures_T=config["target_probe_secondary_structures_T"],
        target_probe_secondary_structures_threshold_deltaG=config[
            "target_probe_secondary_structures_threshold_deltaG"
        ],
        target_probe_homopolymeric_base_n=config["target_probe_homopolymeric_base_n"],
        target_probe_prohibited_sequences=target_probe_prohibited_sequences,
        target_probe_max_len_selfcomplement=config["target_probe_max_len_selfcomplement"],
        target_probe_Tm_parameters=preprocess_tm_parameters(config["target_probe_Tm_parameters"]),
        target_probe_Tm_chem_correction_parameters=config["target_probe_Tm_chem_correction_parameters"],
        target_probe_Tm_salt_correction_parameters=config["target_probe_Tm_salt_correction_parameters"],
        # Step 3: Specificity Filter Parameters
        files_fasta_reference_database_target_probe=config["files_fasta_reference_database_target_probe"],
        files_vcf_reference_database_target_probe=config["files_vcf_reference_database_target_probe"],
        target_probe_cross_hybridization_alignment_method=config[
            "target_probe_cross_hybridization_alignment_method"
        ],
        target_probe_cross_hybridization_search_parameters=config[
            "target_probe_cross_hybridization_search_parameters"
        ],
        target_probe_cross_hybridization_hit_parameters=config[
            "target_probe_cross_hybridization_hit_parameters"
        ],
        target_probe_hybridization_probability_alignment_method=config[
            "target_probe_hybridization_probability_alignment_method"
        ],
        target_probe_hybridization_probability_search_parameters=config[
            "target_probe_hybridization_probability_search_parameters"
        ],
        target_probe_hybridization_probability_hit_parameters=config[
            "target_probe_hybridization_probability_hit_parameters"
        ],
        target_probe_hybridization_probability_threshold=config[
            "target_probe_hybridization_probability_threshold"
        ],
        target_probe_read_length_bias=config["target_probe_read_length_bias"],
        # Step 4: Probe Scoring and Set Selection Parameters
        target_probe_targeted_exons=config["target_probe_targeted_exons"],
        target_probe_targeted_exons_weight=config["target_probe_targeted_exons_weight"],
        target_probe_isoform_weight=config["target_probe_isoform_weight"],
        target_probe_GC_content_opt=config["target_probe_GC_content_opt"],
        target_probe_GC_weight=config["target_probe_GC_weight"],
        target_probe_Tm_opt=config["target_probe_Tm_opt"],
        target_probe_Tm_weight=config["target_probe_Tm_weight"],
        n_sets=config["n_sets"],
        set_size_min=config["set_size_min"],
        set_size_opt=config["set_size_opt"],
        distance_between_target_probes=config["distance_between_target_probes"],
        uniform_distance_weight=config["uniform_distance_weight"],
        n_attempts_graph=config["n_attempts_graph"],
        n_attempts_clique_enum=config["n_attempts_clique_enum"],
        diversification_fraction=config["diversification_fraction"],
        jaccard_opt=config["jaccard_opt"],
        jaccard_step=config["jaccard_step"],
    )

    pipeline.generate_output(oligo_database=oligo_database)

    logging.info("--------------END PIPELINE--------------")


if __name__ == "__main__":
    main()
