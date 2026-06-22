############################################
# imports
############################################

import itertools
import os
import random
import shutil
from pathlib import Path
from typing import Any

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
    BasePropertyFilter,
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
    pipeline_step_basic,
    preprocess_tm_parameters,
)
from oligo_designer_toolsuite.sequence_generator import OligoSequenceGenerator
from oligo_designer_toolsuite.utils import configure_root_logger, logger

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

        ##### set class parameters #####
        self.write_intermediate_steps = write_intermediate_steps
        self.n_jobs = n_jobs

    def design_target_probes(
        self,
        oligo_generation_parameters: dict,
        property_filters_parameters: dict,
        specificity_filters_parameters: dict,
        probe_set_selection_parameters: dict,
        padlock_arms_parameters: dict,
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

        :param oligo_generation_parameters: ``target_probe.oligo_generation`` block. Must contain
            ``region_ids`` (populated from ``file_region_ids``), ``probe_length_min``, ``probe_length_max``,
            and ``files_fasta_probe_database``.
        :type oligo_generation_parameters: dict
        :param property_filters_parameters: ``target_probe.property_filters`` block. Each filter sub-dict
            (``isoform_consensus_filter``, ``hard_masked_sequences_filter``, ``soft_masked_sequences_filter``,
            ``homopolymeric_runs_filter``, ``GC_content_filter``, ``Tm_filter``, ``detection_oligo_filter``)
            carries an ``enabled`` flag plus its parameters.
        :type property_filters_parameters: dict
        :param specificity_filters_parameters: ``target_probe.specificity_filters`` block. Contains the
            shared ``files_fasta_reference_database`` and ``ligation_region_size`` plus the
            ``specificity_blastn_filter`` and ``cross_hybridization_blastn_filter`` sub-dicts, each with
            ``enabled``, ``search_parameters``, and ``hit_parameters``.
        :type specificity_filters_parameters: dict
        :param probe_set_selection_parameters: ``target_probe.probe_set_selection`` block. Contains the
            ``independent_set_selection`` scalars and the ``isoform_consensus_score`` / ``GC_content_score`` /
            ``Tm_score`` sub-dicts.
        :type probe_set_selection_parameters: dict
        :param padlock_arms_parameters: ``target_probe.padlock_arms_properties`` block. Used to compute padlock-arm
            sequences/Tm (via ``PadlockArmsProperty``) for downstream filtering and backbone assembly.
        :type padlock_arms_parameters: dict
        :return: An `OligoDatabase` object containing the designed target probes organized into sets.
            The database includes probe sequences, properties, and set assignments for each target gene.
        :rtype: OligoDatabase
        """

        target_probe_designer = TargetProbeDesigner(self.dir_output, self.n_jobs)

        oligo_database: OligoDatabase = target_probe_designer.create_oligo_database(
            region_ids=oligo_generation_parameters["region_ids"],
            oligo_length_min=oligo_generation_parameters["probe_length_min"],
            oligo_length_max=oligo_generation_parameters["probe_length_max"],
            files_fasta_oligo_database=oligo_generation_parameters["files_fasta_probe_database"],
            min_oligos_per_gene=probe_set_selection_parameters["independent_set_selection"]["set_size_min"],
        )

        if self.write_intermediate_steps:
            dir_database = oligo_database.save_database(name_database="1_db_probes_initial")
            logger.info(f"Saved probe database for step 1 (Create Database) in directory {dir_database}")

        oligo_database = target_probe_designer.filter_by_property(
            oligo_database=oligo_database,
            isoform_consensus_filter=property_filters_parameters["isoform_consensus_filter"],
            hard_masked_sequences_filter=property_filters_parameters["hard_masked_sequences_filter"],
            soft_masked_sequences_filter=property_filters_parameters["soft_masked_sequences_filter"],
            homopolymeric_runs_filter=property_filters_parameters["homopolymeric_runs_filter"],
            GC_content_filter=property_filters_parameters["GC_content_filter"],
            Tm_filter=property_filters_parameters["Tm_filter"],
            detection_oligo_filter=property_filters_parameters["detection_oligo_filter"],
        )

        if self.write_intermediate_steps:
            dir_database = oligo_database.save_database(name_database="2_db_probes_property_filter")
            logger.info(f"Saved probe database for step 2 (Property Filters) in directory {dir_database}")

        ##### compute padlock arm properties #####
        # Required by assemble_padlock_backbone and by the seed-region BLASTN filter when ligation_region_size > 0.
        padlock_arms_property = PadlockArmsProperty(
            arm_length_min=padlock_arms_parameters["padlock_arm_length_min"],
            arm_Tm_dif_max=padlock_arms_parameters["padlock_arm_Tm_dif_max"],
            arm_Tm_min=padlock_arms_parameters["padlock_arm_Tm_min"],
            arm_Tm_max=padlock_arms_parameters["padlock_arm_Tm_max"],
            Tm_parameters=padlock_arms_parameters["Tm_parameters"],
            Tm_chem_correction_parameters=padlock_arms_parameters["Tm_chem_correction_parameters"],
            Tm_salt_correction_parameters=padlock_arms_parameters["Tm_salt_correction_parameters"],
        )
        calculator = PropertyCalculator(properties=[padlock_arms_property])
        oligo_database = calculator.apply(
            oligo_database=oligo_database, sequence_type="oligo", n_jobs=self.n_jobs
        )

        oligo_database = target_probe_designer.filter_by_specificity(
            oligo_database=oligo_database,
            specificity_blastn_filter=specificity_filters_parameters["specificity_blastn_filter"],
            cross_hybridization_blastn_filter=specificity_filters_parameters[
                "cross_hybridization_blastn_filter"
            ],
        )

        if self.write_intermediate_steps:
            dir_database = oligo_database.save_database(name_database="3_db_probes_specificity_filter")
            logger.info(f"Saved probe database for step 3 (Specificity Filters) in directory {dir_database}")

        oligo_database = target_probe_designer.create_oligo_sets(
            oligo_database=oligo_database,
            independent_set_selection=probe_set_selection_parameters["independent_set_selection"],
            isoform_consensus_score=probe_set_selection_parameters["isoform_consensus_score"],
            GC_content_score=probe_set_selection_parameters["GC_content_score"],
            Tm_score=probe_set_selection_parameters["Tm_score"],
        )

        # Calculate oligo length, GC content, Tm, and isoform consensus on the selected oligos.
        length_property = LengthProperty()
        gc_content_property = GCContentProperty()
        TmNN_property = TmNNProperty(
            Tm_parameters=property_filters_parameters["Tm_filter"]["Tm_parameters"],
            Tm_chem_correction_parameters=property_filters_parameters["Tm_filter"][
                "Tm_chem_correction_parameters"
            ],
            Tm_salt_correction_parameters=property_filters_parameters["Tm_filter"][
                "Tm_salt_correction_parameters"
            ],
        )
        isoform_consensus_property = IsoformConsensusProperty()
        calculator = PropertyCalculator(
            properties=[
                length_property,
                gc_content_property,
                TmNN_property,
                isoform_consensus_property,
            ]
        )

        oligo_database = calculator.apply(
            oligo_database=oligo_database, sequence_type="oligo", n_jobs=self.n_jobs
        )

        if self.write_intermediate_steps:
            dir_database = oligo_database.save_database(name_database="4_db_probes_probesets")
            logger.info(f"Saved probe database for step 4 (Specificity Filters) in directory {dir_database}.")

        return oligo_database

    def design_detection_oligos(
        self,
        oligo_database: OligoDatabase,
        oligo_generation_parameters: dict,
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
        :param oligo_generation_parameters: ``detection_oligo.oligo_generation`` block. Must contain
            ``oligo_length_min``, ``oligo_length_max``, ``min_thymines``, ``U_distance``, and ``Tm_opt``.
        :type oligo_generation_parameters: dict
        :return: An updated `OligoDatabase` object containing the designed detection oligos. The
            database includes the following new sequence properties for each probe:
            - `sequence_detection_oligo`: The detection oligo sequence with uracil substitutions
            - `Tm_detection_oligo`: The melting temperature of the detection oligo
        :rtype: OligoDatabase
        """

        detection_oligo_designer = DetectionOligoDesigner(self.n_jobs)
        oligo_database = detection_oligo_designer.create_detection_oligos(
            oligo_database=oligo_database,
            oligo_length_min=oligo_generation_parameters["oligo_length_min"],
            oligo_length_max=oligo_generation_parameters["oligo_length_max"],
            min_thymines=oligo_generation_parameters["min_thymines"],
            U_distance=oligo_generation_parameters["U_distance"],
            Tm_opt=oligo_generation_parameters["Tm_opt"],
            Tm_parameters=oligo_generation_parameters["Tm_parameters"],
            Tm_chem_correction_parameters=oligo_generation_parameters["Tm_chem_correction_parameters"],
            Tm_salt_correction_parameters=oligo_generation_parameters["Tm_salt_correction_parameters"],
        )

        return oligo_database

    def assemble_padlock_backbone(
        self,
        oligo_database: OligoDatabase,
        padlock_arms_parameters: dict,
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

        Tm_parameters = padlock_arms_parameters["Tm_parameters"]
        Tm_chem_correction_parameters = padlock_arms_parameters["Tm_chem_correction_parameters"]
        Tm_salt_correction_parameters = padlock_arms_parameters["Tm_salt_correction_parameters"]

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
                        Tm_parameters=Tm_parameters,
                        Tm_chem_correction_parameters=Tm_chem_correction_parameters,
                        Tm_salt_correction_parameters=Tm_salt_correction_parameters,
                    )
                    Tm_arm2 = calc_tm_nn(
                        sequence=sequence_padlock_arm2,
                        Tm_parameters=Tm_parameters,
                        Tm_chem_correction_parameters=Tm_chem_correction_parameters,
                        Tm_salt_correction_parameters=Tm_salt_correction_parameters,
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
    ) -> OligoDatabase:
        """
        Create an initial oligo database by generating sequences using a sliding window approach.

        This is the first step of target probe design. The method:
        1. Generates candidate oligo sequences from input FASTA files using a sliding window approach
           across the specified length range
        2. Creates an `OligoDatabase` and loads the generated sequences
        3. Calculates the reverse complement sequence (``oligo`` sequence type) — always on
        4. Removes regions with insufficient oligos

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
        :return: An `OligoDatabase` object containing the generated target probe sequences with their
            component sequences (target, oligo). The database is filtered to only include regions that
            meet the minimum oligo requirement.
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
        # Compute the reverse complement -> "oligo" sequence type (always on).
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
    def filter_by_property(
        self,
        oligo_database: OligoDatabase,
        isoform_consensus_filter: dict,
        hard_masked_sequences_filter: dict,
        soft_masked_sequences_filter: dict,
        homopolymeric_runs_filter: dict,
        GC_content_filter: dict,
        Tm_filter: dict,
        detection_oligo_filter: dict,
    ) -> OligoDatabase:
        """
        Filter the oligo database based on sequence properties to remove probes with undesirable
        characteristics.

        This method applies sequential filtering using multiple property-based filters, each gated on
        its own ``enabled`` flag:

        1. **Isoform consensus** (cheap pre-filter): computes ``IsoformConsensusProperty`` on the
           ``target`` sequence and removes regions below the configured threshold.
        2. **Hard masked sequences**: Removes probes containing hard-masked nucleotides (N)
        3. **Soft masked sequences**: Removes probes containing soft-masked nucleotides (lowercase)
        4. **Homopolymeric runs**: Removes probes with homopolymeric runs exceeding specified lengths
        5. **GC content**: Removes probes with GC content outside the specified range
        6. **Melting temperature**: Removes probes with Tm outside the specified range
        7. **Padlock arm / detection oligo requirements**: Removes probes that cannot form valid padlock
           arms with balanced melting temperatures, nor valid detection oligos centered on the ligation
           site with sufficient thymines for UNG cleavage.

        Probes that fail any enabled filter are removed. Regions with insufficient oligos after
        filtering are removed from the database.

        :param oligo_database: The `OligoDatabase` instance containing oligonucleotide sequences
            and their associated properties.
        :type oligo_database: OligoDatabase
        :param isoform_consensus_filter: Dict with ``enabled``, ``isoform_consensus``.
        :type isoform_consensus_filter: dict
        :param hard_masked_sequences_filter: Dict with ``enabled``.
        :type hard_masked_sequences_filter: dict
        :param soft_masked_sequences_filter: Dict with ``enabled``.
        :type soft_masked_sequences_filter: dict
        :param homopolymeric_runs_filter: Dict with ``enabled``, ``homopolymeric_base_n`` (mapping
            ``A``/``T``/``C``/``G`` to maximum allowed run lengths).
        :type homopolymeric_runs_filter: dict
        :param GC_content_filter: Dict with ``enabled``, ``GC_content_min``, ``GC_content_max``.
        :type GC_content_filter: dict
        :param Tm_filter: Dict with ``enabled``, ``Tm_min``, ``Tm_max``, plus thermodynamic model parameters
            (``Tm_parameters``, ``Tm_chem_correction_parameters``, ``Tm_salt_correction_parameters``) injected during
            config preprocessing.
        :type Tm_filter: dict
        :param detection_oligo_filter: Dict with detection-oligo and padlock-arm constraints required by the
            composite ``DetectionOligoFilter``.
        :type detection_oligo_filter: dict
        :return: A filtered `OligoDatabase` object containing only probes that pass all enabled property
            filters. Regions with insufficient oligos after filtering are removed.
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

        # Note: detetcion oligo filter already checks if padlock arms are feasible.
        # No need to apply PadlockArmsFilter here.
        detection_oligo = DetectionOligoFilter(
            detect_oligo_length_min=detection_oligo_filter["oligo_length_min"],
            detect_oligo_length_max=detection_oligo_filter["oligo_length_max"],
            min_thymines=detection_oligo_filter["min_thymines"],
            arm_length_min=detection_oligo_filter["padlock_arm_length_min"],
            arm_Tm_dif_max=detection_oligo_filter["padlock_arm_Tm_dif_max"],
            arm_Tm_min=detection_oligo_filter["padlock_arm_Tm_min"],
            arm_Tm_max=detection_oligo_filter["padlock_arm_Tm_max"],
            Tm_parameters=detection_oligo_filter["Tm_parameters"],
            Tm_chem_correction_parameters=detection_oligo_filter["Tm_chem_correction_parameters"],
            Tm_salt_correction_parameters=detection_oligo_filter["Tm_salt_correction_parameters"],
        )
        filters.append(detection_oligo)

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
        specificity_blastn_filter: dict,
        cross_hybridization_blastn_filter: dict,
    ) -> OligoDatabase:
        """
        Filter the oligo database based on sequence specificity to remove probes that bind
        non-specifically or cross-hybridize.

        Before applying the alignment-based filters, the method computes the ``PadlockArmsProperty``
        for every probe (always on) — this is required by the seed-region BLASTN filter when
        ``ligation_region_size > 0`` and by :py:meth:`ScrinshotProbeDesigner.assemble_padlock_backbone`
        downstream.

        The filter list is seeded with an :class:`ExactMatchFilter` (always on) and then conditionally
        extended with the BLASTN-specificity and cross-hybridization filters depending on their
        ``enabled`` flags. All filters are applied in a single :class:`SpecificityFilter` invocation
        so the database is iterated once.

        1. **Exact matches** (always on): Removes all probes with exact sequence matches to probes of
           other regions.
        2. **BLASTN specificity** (gated on ``specificity_blastn_filter['enabled']``): Uses BLASTN to
           search for similar sequences in the reference database. Probes with hits meeting the
           specified criteria are removed. If ``ligation_region_size > 0``, seed-based filtering is
           applied around the ligation site, removing all probes where BLASTN hits cover the junction
           region, regardless of the coverage threshold. If ``ligation_region_size == 0``, full-length
           specificity filtering is performed.
        3. **Cross-hybridization** (gated on ``cross_hybridization_blastn_filter['enabled']``): Removes
           probes that cross-hybridize with each other. This is critical because if probes can bind to
           each other, they may form dimers instead of binding to the target RNA. Probes from the
           larger genomic region are removed when cross-hybridization is detected.

        The reference database is loaded from the provided FASTA files and shared between the BLASTN
        specificity and cross-hybridization filters. Regions that do not meet the minimum oligo
        requirement after filtering are removed from the database.

        :param oligo_database: The `OligoDatabase` instance containing oligonucleotide sequences
            and their associated properties.
        :type oligo_database: OligoDatabase
        :param files_fasta_reference_database: List of paths to FASTA files containing reference
            sequences against which specificity will be evaluated. These typically include the
            entire genome or transcriptome to identify off-target binding sites. The database is
            shared between both BLASTN-based filters.
        :type files_fasta_reference_database: list[str]
        :param ligation_region_size: Size of the ligation region (in nucleotides) around the ligation
            site. If > 0, seed-based specificity filtering is applied: all probes where BLASTN hits
            cover the junction region are removed, regardless of the coverage threshold. If 0,
            full-length specificity filtering is performed. Both modes perform full BLASTN searches.
        :type ligation_region_size: int
        :param specificity_blastn_filter: Dict with ``enabled``, ``search_parameters``,
            ``hit_parameters``.
        :type specificity_blastn_filter: dict
        :param cross_hybridization_blastn_filter: Dict with ``enabled``, ``search_parameters``,
            ``hit_parameters``.
        :type cross_hybridization_blastn_filter: dict
        :param padlock_arm_filter: Dict with ``enabled``, ``arm_length_min``, ``arm_Tm_dif_max``,
            ``arm_Tm_min``, ``arm_Tm_max``. Provides the constraints used to compute the always-on
            ``PadlockArmsProperty``. The ``enabled`` flag here is consulted by
            :py:meth:`filter_by_property`; the property computation in this method runs unconditionally
            since the arm sequences are required downstream.
        :type padlock_arm_filter: dict
        :return: A filtered `OligoDatabase` object containing only probes that pass all enabled
            specificity and cross-hybridization filters. The database includes calculated padlock arm
            properties (``ligation_site``, arm sequences, arm Tm values). Regions with insufficient
            oligos after filtering are removed.
        :rtype: OligoDatabase
        """
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
            reference_database = ReferenceDatabase(
                database_name=f"{self.subdir_db_reference}_sequences", dir_output=self.dir_output
            )
            reference_database.load_database_from_file(
                files=specificity_blastn_filter["files_fasta_reference_database"],
                file_type="fasta",
                database_overwrite=True,
            )
            specificity: AlignmentSpecificityFilter
            if specificity_blastn_filter["ligation_region_size"] > 0:
                specificity = BlastNSeedregionSiteFilter(
                    remove_hits=True,
                    seedregion_size=specificity_blastn_filter["ligation_region_size"],
                    seedregion_site_name="ligation_site",
                    search_parameters=specificity_blastn_filter["search_parameters"],
                    hit_parameters=specificity_blastn_filter["hit_parameters"],
                    filter_name="specificity_blastn_filter",
                    dir_output=self.dir_output,
                )
            else:
                specificity = BlastNFilter(
                    remove_hits=True,
                    search_parameters=specificity_blastn_filter["search_parameters"],
                    hit_parameters=specificity_blastn_filter["hit_parameters"],
                    filter_name="specificity_blastn_filter",
                    dir_output=self.dir_output,
                )
            specificity.set_reference_database(reference_database=reference_database)
            filters.append(specificity)
            directories.append(specificity.dir_output)

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

    @pipeline_step_basic(step_name="Set Selection")
    def create_oligo_sets(
        self,
        oligo_database: OligoDatabase,
        independent_set_selection: dict,
        isoform_consensus_score: dict,
        GC_content_score: dict,
        Tm_score: dict,
    ) -> OligoDatabase:
        """
        Create optimal oligo sets based on weighted scoring criteria, distance constraints, and set selection.

        This method performs the following steps:
        1. **Scoring**: Calculates scores for each oligo based on weighted criteria (isoform consensus,
           GC content, melting temperature). Higher scores indicate better probes.
        2. **Set generation**: Builds a compatibility graph from distance constraints and selects sets via
           a graph-based (clique) strategy. Generates multiple diverse sets per region, controlling overlap
           between sets using a Jaccard threshold (``jaccard_opt``) with optional relaxation (``jaccard_step``).
        3. **Set scoring**: Evaluates each generated set and selects the best sets based on the lowest
           average score (ascending order).
        4. **Region filtering**: Removes regions that cannot generate sets meeting the minimum size requirement.

        The algorithm attempts to find sets with optimal size (``set_size_opt``) but may produce sets
        as small as ``set_size_min`` if constraints cannot be met.

        :param oligo_database: The `OligoDatabase` instance containing oligonucleotide sequences
            and their associated properties. This database should contain target probes that have
            passed all previous filtering steps, including padlock arm property calculations.
        :type oligo_database: OligoDatabase
        :param independent_set_selection: Dict controlling set generation. Must contain ``n_sets``,
            ``set_size_min``, ``set_size_opt``, ``distance_between_probes``, ``n_attempts_graph``,
            ``n_attempts_clique_enum``, ``diversification_fraction``, ``jaccard_opt``, ``jaccard_step``.
        :type independent_set_selection: dict
        :param isoform_consensus_score: Dict with ``weight``.
        :type isoform_consensus_score: dict
        :param GC_content_score: Dict with ``weight``, ``GC_content_min``, ``GC_content_opt``,
            ``GC_content_max``.
        :type GC_content_score: dict
        :param Tm_score: Dict with ``weight``, ``Tm_min``, ``Tm_opt``, ``Tm_max``.
        :type Tm_score: dict
        :return: An updated `OligoDatabase` object containing the generated oligo sets. Each region
            will have up to ``n_sets`` sets stored, with each set containing between ``set_size_min`` and
            ``set_size_opt`` probes. Regions with insufficient oligos are removed.
        :rtype: OligoDatabase
        """

        # Define all scorers
        isoform_consensus_scorer = IsoformConsensusScorer(
            score_weight=isoform_consensus_score["weight"],
            property_name_transcript_id="transcript_id",
            property_name_number_total_transcripts="number_total_transcripts",
        )
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
        oligos_scoring = OligoScoring(scorers=[isoform_consensus_scorer, Tm_scorer, GC_scorer])
        set_scoring = LowestSetScoring(ascending=True)

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
# SCRINSHOT Designer Pipeline
############################################


def _preprocess_config(config: dict[str, Any]) -> dict[str, Any]:

    # Preprocess Tm tables and set Tm_chem/salt_correction_parameters to None if the correction is disabled
    for section in ["target_probe", "detection_oligo"]:
        config[section]["global_parameters"]["Tm_parameters"] = preprocess_tm_parameters(
            config[section]["global_parameters"]["Tm_parameters"]
        )
        for correction in ["Tm_chem_correction_parameters", "Tm_salt_correction_parameters"]:
            correction_cfg = config[section]["global_parameters"][correction]
            if not correction_cfg["enabled"]:
                correction_cfg["parameters"] = None

    target_probe_Tm_parameters = config["target_probe"]["global_parameters"]["Tm_parameters"]
    target_probe_Tm_chem_correction_parameters = config["target_probe"]["global_parameters"][
        "Tm_chem_correction_parameters"
    ]["parameters"]
    target_probe_Tm_salt_correction_parameters = config["target_probe"]["global_parameters"][
        "Tm_salt_correction_parameters"
    ]["parameters"]

    # Note: for the detection oligo filter we have to set the Tm parameters to the ones form the target
    # probe because it is used to check if padlock arms are feasible.
    config["target_probe"]["property_filters"]["detection_oligo_filter"] = {
        "oligo_length_min": config["detection_oligo"]["oligo_generation"]["oligo_length_min"],
        "oligo_length_max": config["detection_oligo"]["oligo_generation"]["oligo_length_max"],
        "min_thymines": config["detection_oligo"]["oligo_generation"]["min_thymines"],
        "padlock_arm_length_min": config["target_probe"]["padlock_arms_properties"]["padlock_arm_length_min"],
        "padlock_arm_Tm_dif_max": config["target_probe"]["padlock_arms_properties"]["padlock_arm_Tm_dif_max"],
        "padlock_arm_Tm_min": config["target_probe"]["padlock_arms_properties"]["padlock_arm_Tm_min"],
        "padlock_arm_Tm_max": config["target_probe"]["padlock_arms_properties"]["padlock_arm_Tm_max"],
        "Tm_parameters": target_probe_Tm_parameters,
        "Tm_chem_correction_parameters": target_probe_Tm_chem_correction_parameters,
        "Tm_salt_correction_parameters": target_probe_Tm_salt_correction_parameters,
    }

    config["target_probe"]["property_filters"]["Tm_filter"]["Tm_parameters"] = target_probe_Tm_parameters
    config["target_probe"]["property_filters"]["Tm_filter"][
        "Tm_chem_correction_parameters"
    ] = target_probe_Tm_chem_correction_parameters
    config["target_probe"]["property_filters"]["Tm_filter"][
        "Tm_salt_correction_parameters"
    ] = target_probe_Tm_salt_correction_parameters

    config["target_probe"]["padlock_arms_properties"]["Tm_parameters"] = target_probe_Tm_parameters
    config["target_probe"]["padlock_arms_properties"][
        "Tm_chem_correction_parameters"
    ] = target_probe_Tm_chem_correction_parameters
    config["target_probe"]["padlock_arms_properties"][
        "Tm_salt_correction_parameters"
    ] = target_probe_Tm_salt_correction_parameters

    config["target_probe"]["probe_set_selection"]["Tm_score"]["Tm_parameters"] = target_probe_Tm_parameters
    config["target_probe"]["probe_set_selection"]["Tm_score"][
        "Tm_chem_correction_parameters"
    ] = target_probe_Tm_chem_correction_parameters
    config["target_probe"]["probe_set_selection"]["Tm_score"][
        "Tm_salt_correction_parameters"
    ] = target_probe_Tm_salt_correction_parameters

    detection_oligo_Tm_parameters = config["detection_oligo"]["global_parameters"]["Tm_parameters"]
    detection_oligo_Tm_chem_correction_parameters = config["detection_oligo"]["global_parameters"][
        "Tm_chem_correction_parameters"
    ]["parameters"]
    detection_oligo_Tm_salt_correction_parameters = config["detection_oligo"]["global_parameters"][
        "Tm_salt_correction_parameters"
    ]["parameters"]

    config["detection_oligo"]["oligo_generation"]["Tm_parameters"] = detection_oligo_Tm_parameters
    config["detection_oligo"]["oligo_generation"][
        "Tm_chem_correction_parameters"
    ] = detection_oligo_Tm_chem_correction_parameters
    config["detection_oligo"]["oligo_generation"][
        "Tm_salt_correction_parameters"
    ] = detection_oligo_Tm_salt_correction_parameters

    ##### read the genes file #####
    file_region_ids = config["target_probe"]["oligo_generation"]["file_region_ids"]
    if file_region_ids is None:
        logger.warning(
            "No gene list file was provided! All genes from fasta file are used to generate the probes. This choice can use a lot of resources."
        )
        config["target_probe"]["oligo_generation"]["region_ids"] = None
    else:
        with open(file_region_ids) as f:
            config["target_probe"]["oligo_generation"]["region_ids"] = sorted({line.rstrip() for line in f})

    return config


def scrinshot_probe_designer(config: dict[str, Any]) -> None:
    """
    Execute the SCRINSHOT probe design pipeline from a (raw) configuration dict.

    The dict is expected to follow the nested layout of ``data/configs/scrinshot_probe_designer.yaml``
    (``general``, ``target_probe.*``, ``detection_oligo.*``). The caller is responsible for configuring
    the library logger before invoking this function (see :func:`main`).

    :param config: Pipeline configuration loaded via ``yaml.safe_load``.
    :type config: dict
    """

    ##### preprocess the config file #####
    config_dict = _preprocess_config(config)

    ##### initialize probe designer pipeline #####
    pipeline = ScrinshotProbeDesigner(
        write_intermediate_steps=config_dict["general"]["write_intermediate_steps"],
        dir_output=config_dict["general"]["dir_output"],
        n_jobs=config_dict["general"]["n_jobs"],
    )

    ##### design target probes #####
    oligo_database = pipeline.design_target_probes(
        oligo_generation_parameters=config_dict["target_probe"]["oligo_generation"],
        property_filters_parameters=config_dict["target_probe"]["property_filters"],
        specificity_filters_parameters=config_dict["target_probe"]["specificity_filters"],
        probe_set_selection_parameters=config_dict["target_probe"]["probe_set_selection"],
        padlock_arms_parameters=config_dict["target_probe"]["padlock_arms_properties"],
    )

    ##### design detection oligos #####
    oligo_database = pipeline.design_detection_oligos(
        oligo_database=oligo_database,
        oligo_generation_parameters=config_dict["detection_oligo"]["oligo_generation"],
    )

    ##### assemble padlock backbone #####
    probe_database = pipeline.assemble_padlock_backbone(
        oligo_database=oligo_database,
        padlock_arms_parameters=config_dict["target_probe"]["padlock_arms_properties"],
    )

    ##### write outputs #####
    pipeline.generate_output(probe_database=probe_database)


def main() -> None:
    """
    Main entry point for running the SCRINSHOT probe design pipeline.

    Parses ``--config``, loads the YAML, configures the library logger to write into the configured
    output directory, then delegates to :func:`scrinshot_probe_designer`.
    """
    print("--------------START PIPELINE--------------")

    args = base_parser()

    ##### read the config file #####
    with open(args["config"], "r") as handle:
        config = yaml.safe_load(handle)

    # setup logger now that we know the output directory
    configure_root_logger(
        dir_output=config["general"]["dir_output"],
        pipeline_name="scrinshot_probe_designer",
    )

    scrinshot_probe_designer(config)

    print("--------------END PIPELINE--------------")


if __name__ == "__main__":
    main()
