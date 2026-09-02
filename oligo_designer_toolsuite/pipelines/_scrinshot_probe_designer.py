"""
SCRINSHOT probe designer pipeline.

SCRINSHOT, or Single-Cell Resolution IN Situ Hybridization On Tissues, is a
targeted spatial transcriptomics method for detecting RNA molecules in fixed
tissue sections. It uses padlock probes that are circularized on the RNA target
and then amplified by rolling circle amplification.

See :class:`ScrinshotProbeDesigner` for the full pipeline description and probe
structure. See :func:`scrinshot_probe_designer` for the config-driven workflow.
"""

############################################
# imports
############################################

import itertools
import os
import random
import shutil
from pathlib import Path
from typing import Any

from joblib import Parallel, delayed
from joblib_progress import joblib_progress

from oligo_designer_toolsuite.config.pipelines.scrinshot_probe_designer import ScrinshotProbeDesignerConfig
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
    apply_required_parameters,
    base_log_parameters,
    base_parser,
    check_content_oligo_database,
    load_config,
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
    Design padlock probes and detection oligos for SCRINSHOT experiments.

    This class runs the design workflow for SCRINSHOT probe libraries. It
    designs target-binding padlock probes, designs detection oligos for the
    rolling-circle products, assembles the padlock backbone, and writes the
    final order-ready files.

    Overview
    --------
    SCRINSHOT, or Single-Cell Resolution IN Situ Hybridization On Tissues, is a
    targeted spatial transcriptomics method for detecting RNA molecules in fixed
    tissue sections.

    The method uses padlock probes that hybridize directly to mRNA. When both
    target-binding arms of a padlock probe bind next to each other on the RNA,
    SplintR ligase can join the probe ends and circularize the padlock. The
    circularized probe is then amplified by rolling circle amplification. This
    creates a local RCA product at the RNA molecule.

    RCA products are detected with short fluorescent detection oligos. Because
    the amplified product contains many repeated copies of the padlock sequence,
    several detection oligos can bind to one RCA product. This gives a bright
    signal that can be counted in tissue sections.

    SCRINSHOT can be used to map marker-gene expression across many cells while
    preserving tissue context. In the original study, marker-gene measurements
    showed high correlation with published scRNA-seq data, and the method was
    used to map abundant and rare cell types in tissue sections.

    Probe Structure
    ---------------
    **Padlock Probes**

    Padlock probes are single-stranded DNA oligos that bind directly to target
    RNA. Each padlock probe contains two target-binding arms and a backbone.

    The two arms bind to adjacent regions of the target RNA and flank the
    ligation site. If both arms bind correctly, SplintR ligase can close the
    padlock into a circle. The circle then serves as the template for rolling
    circle amplification.

    Each padlock probe contains:

    - two target-binding arms, usually around 20 nucleotides each,
    - a constant backbone used for priming and detection,
    - a gene-specific barcode sequence in the backbone.

    In this pipeline, the backbone is assembled from constant accessory
    sequences, an ISS anchor sequence, and a barcode assigned to the target
    region.

    A simplified padlock layout is::

        [target-binding arm 1] + [backbone with barcode] + [target-binding arm 2]

    The ligation site lies between the two target-binding arms after the probe
    has hybridized to the RNA. This placement makes ligation dependent on
    correct binding at the target site.

    **Detection Oligos**

    Detection oligos are short single-stranded DNA probes that bind to the RCA
    product generated from the circularized padlock probe. They carry a
    fluorophore, so each RCA product can be detected as a bright spot during
    imaging.

    Detection oligos are designed around the padlock ligation site. Their
    melting temperature is chosen to fit the imaging conditions. For sequential
    imaging, detection oligos can include uracil bases. These uracils allow the
    detection oligos to be cleaved by Uracil DNA Glycosylase, so one signal can
    be removed before the next detection round.

    The fluorophore depends on the imaging setup. Common choices include FITC,
    Cy3, and Cy5. Other dyes can be used if the microscope and protocol support
    them.

    Probe Library Preparation
    -------------------------
    SCRINSHOT padlock probes and detection oligos are usually ordered as
    synthetic DNA oligos. Padlock probes are hybridized to RNA in fixed tissue
    sections. After ligation, circularized padlocks are amplified by rolling
    circle amplification.

    Fluorescent detection oligos are then added to read out the RCA products.
    For multiplexed experiments, detection can be performed over several imaging
    cycles. After each cycle, the detection signal is removed, and another set
    of detection oligos is applied.

    Pipeline Overview
    -----------------
    The pipeline performs the main steps needed to design a SCRINSHOT probe set:

    1. **Target probe design**

       Design the gene-specific target sequence that will form the two
       target-binding arms of the padlock probe.

    2. **Detection oligo design**

       Design detection oligos for the RCA products. These oligos are centered
       around the ligation site and can include uracils for sequential imaging.

    3. **Padlock backbone assembly**

       Combine target-binding arms, the constant backbone, and the target
       barcode to build the full padlock probe.

    4. **Output generation**

       Write the ready-to-order padlock probes, detection oligos, and probe
       annotations to output files.

    References
    ----------
    Sountoulidis, A., Liontos, A., Nguyen, H. P., Firsova, A. B.,
    Fysikopoulos, A., Qian, X., et al. (2020). SCRINSHOT enables spatial
    mapping of cell states in tissue sections with single-cell resolution.
    PLOS Biology, 18(11), e3000675.
    https://doi.org/10.1371/journal.pbio.3000675

    :param write_intermediate_steps: If ``True``, save intermediate probe
        databases after pipeline steps. This can help with checking a design run
        or finding where probes were removed.
    :type write_intermediate_steps: bool
    :param dir_output: Directory where output files and intermediate results are
        saved. The directory is created if it does not exist.
    :type dir_output: str
    :param n_jobs: Number of worker processes used for steps that can run in
        parallel. Use ``1`` to run without parallel processing.
    :type n_jobs: int
    """

    def __init__(self, write_intermediate_steps: bool, dir_output: str, n_jobs: int) -> None:
        """Constructor for the ScrinshotProbeDesigner class."""
        self.dir_output = os.path.abspath(dir_output)
        Path(self.dir_output).mkdir(parents=True, exist_ok=True)

        self.write_intermediate_steps = write_intermediate_steps
        self.n_jobs = n_jobs

    def design_target_probes(
        self,
        target_probes_parameters: dict,
    ) -> OligoDatabase:
        """
        Design the RNA-binding part of SCRINSHOT padlock probes.

        This step designs target sequences that bind directly to RNA. Each selected
        target sequence is later split into two padlock arms around the ligation
        site. Candidate probes are generated, filtered for sequence quality and
        specificity, and then selected into final probe sets.

        After the probe sets are selected, this method adds useful reporting values
        to each probe, including probe length, GC content, melting temperature, and
        isoform consensus.

        :param target_probes_parameters: Settings from the ``target_probes`` section
            of the pipeline config. This includes candidate generation, sequence
            filters, specificity filters, and probe set selection.
        :type target_probes_parameters: dict
        :return: Database containing the selected SCRINSHOT target probes for each
            target region.
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
        isoform_consensus_property = IsoformConsensusProperty()
        calculator = PropertyCalculator(
            properties=[
                length_property,
                gc_content_property,
                TmNN_property,
                isoform_consensus_property,
            ]
        )
        target_probes_database = calculator.apply(
            oligo_database=target_probes_database, sequence_type="oligo", n_jobs=self.n_jobs
        )

        return target_probes_database

    def design_detection_oligos(
        self,
        oligo_database: OligoDatabase,
        oligo_generation_parameters: dict,
    ) -> OligoDatabase:
        """
        Design detection oligos for SCRINSHOT RCA products.

        Detection oligos bind to the rolling-circle amplification product generated
        from each circularized padlock probe. They are centered around the padlock
        ligation site and are chosen to match the imaging conditions.

        For sequential imaging, detection oligos can include uracil bases. These
        uracils allow the detection oligo to be cleaved enzymatically, so the signal
        can be removed before the next imaging round.

        :param oligo_database: Database returned by :py:meth:`design_target_probes`.
            This database is updated with detection-oligo sequences and related
            properties.
        :type oligo_database: OligoDatabase
        :param oligo_generation_parameters: Settings from the
            ``detection_oligo.oligo_generation`` section of the pipeline config.
            This includes detection-oligo length, uracil placement, and melting
            temperature settings.
        :type oligo_generation_parameters: dict
        :return: Database with detection-oligo properties added to each probe.
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
        Build the full SCRINSHOT padlock probes.

        This step splits each selected target sequence into two padlock arms and
        inserts the padlock backbone between them. The backbone contains constant
        accessory sequences, an ISS anchor sequence, and a short barcode assigned to
        the target region.

        The resulting padlock probe can bind to RNA through its two arms. When the
        arms bind next to each other on the target, the probe can be circularized by
        ligation and then amplified by rolling circle amplification.

        A simplified layout is::

            [padlock arm 1] + [backbone with barcode] + [padlock arm 2]

        The melting temperature of both arms is calculated and stored, so users can
        check whether the two arms are reasonably balanced.

        :param oligo_database: Database returned by
            :py:meth:`design_detection_oligos`. This database is updated with the
            assembled padlock probe sequences and arm properties.
        :type oligo_database: OligoDatabase
        :param padlock_arms_parameters: Settings from the
            ``target_probes.padlock_arms_properties`` section of the pipeline config.
            This includes the conditions used to calculate melting temperatures for
            the padlock arms.
        :type padlock_arms_parameters: dict
        :return: Database with assembled padlock probes and padlock-arm properties
            added.
        :rtype: OligoDatabase
        """

        def _get_barcode(number_regions: int, barcode_length: int, seed: int, choices: list) -> list[str]:
            """
            Generate barcodes for the SCRINSHOT padlock backbone.

            The method creates all possible barcode sequences of the requested length
            from the allowed bases. If that length is not enough to cover all target
            regions, the barcode length is increased until enough unique barcodes are
            available.

            The resulting barcode list is shuffled with a fixed random seed, so the same
            input settings give the same barcode assignment.

            :param number_regions: Number of target regions that need a barcode.
            :type number_regions: int
            :param barcode_length: Starting barcode length in nucleotides. The length is
                increased if more unique barcodes are needed.
            :type barcode_length: int
            :param seed: Random seed used to shuffle the barcode list reproducibly.
            :type seed: int
            :param choices: Allowed bases used to build the barcode sequences.
            :type choices: list
            :return: Shuffled list of barcode sequences. The list contains at least one
                barcode for each target region.
            :rtype: list[str]
            """
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
                    # Skip probes whose ligation_site was never set (failed arm constraints upstream).
                    if not isinstance(sequence_oligo, str) or not isinstance(ligation_site, int):
                        continue

                    # On the probe strand the junction is arm2 then arm1 (5'→3'); assemble as
                    # arm1–backbone–arm2 so the free ends meet for ligation after hybridization.
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
        oligo_database: OligoDatabase,
        output_properties: list[str] | None = None,
    ) -> None:
        """
        Write the completed SCRINSHOT probe design to files.

        This step saves the final probe database with padlock probe sequences,
        detection oligos, and selected annotations. It also writes an order-ready
        file with the padlock probe and detection oligo sequences needed for
        synthesis.

        If no output properties are provided, a default set of annotation fields,
        padlock fields, detection-oligo fields, and melting-temperature values is
        written.

        :param oligo_database: Database returned by
            :py:meth:`assemble_padlock_backbone`.
        :type oligo_database: OligoDatabase
        :param output_properties: Probe properties to include in the detailed output
            files. If ``None``, a default set of annotations, sequences, and melting
            temperature values is used.
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

        oligo_database.write_oligosets_to_yaml(
            properties=output_properties,
            ascending=True,
            filename="padlock_probes",
        )

        oligo_database.write_oligosets_to_table(
            properties=output_properties,
            ascending=True,
            filename="padlock_probes",
        )

        oligo_database.write_ready_to_order_yaml(
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
    Design the RNA-binding part of SCRINSHOT padlock probes.

    This class designs the target sequence that will later become the two
    target-binding arms of a padlock probe. The sequence binds directly to RNA
    and is split around the ligation site during padlock assembly.

    Good SCRINSHOT target probes need to satisfy several checks at the same
    time. The padlock arms should have suitable and balanced melting
    temperatures. The ligation site should be specific to the intended RNA
    target. The probe also needs enough suitable sequence around the ligation
    site to design a detection oligo for the RCA product.

    The workflow has four main steps:

    1. **Candidate generation**

       Build candidate target probes from transcript or target-region FASTA
       files.

    2. **Sequence filtering**

       Remove candidates with unsuitable sequence properties, such as masked
       sequence, long single-base runs, unsuitable GC content, unsuitable melting
       temperature, or no suitable detection oligo region.

    3. **Specificity filtering**

       Remove candidates that may bind to unintended targets or cross-hybridize
       with other probes in the panel.

    4. **Probe set selection**

       Select final probe sets for each target region, using criteria such as
       isoform coverage, GC content, melting temperature, and probe spacing.

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
        Run the full SCRINSHOT target-probe design workflow.

        This method designs the RNA-binding target sequences used to build padlock
        probes. It starts from transcript or target-region sequences, creates
        candidate probes, filters them, checks their specificity, and selects final
        probe sets for each target region.

        Before specificity filtering, the method calculates padlock-arm properties
        for each candidate. These include the ligation site and arm melting
        temperatures. The ligation site is needed for specificity checks around the
        padlock junction, and the same arm information is later used during padlock
        backbone assembly.

        :param target_probes_parameters: Settings from the ``target_probes`` section
            of the pipeline config. This includes candidate generation, sequence
            filters, specificity filters, padlock-arm settings, and probe set
            selection.
        :type target_probes_parameters: dict
        :param write_intermediate_steps: If ``True``, save intermediate probe
            databases after each main step. This can help when checking where probes
            were removed.
        :type write_intermediate_steps: bool
        :return: Database containing selected SCRINSHOT target probes for each
            target region.
        :rtype: OligoDatabase
        """
        oligo_generation_parameters = target_probes_parameters["oligo_generation"]
        property_filters_parameters = target_probes_parameters["property_filters"]
        specificity_filters_parameters = target_probes_parameters["specificity_filters"]
        probe_set_selection_parameters = target_probes_parameters["probe_set_selection"]
        padlock_arms_parameters = target_probes_parameters["padlock_arms_properties"]

        oligo_database: OligoDatabase = self._create_oligo_database(
            region_ids=oligo_generation_parameters["region_ids"],
            oligo_length_min=oligo_generation_parameters["probe_length_min"],
            oligo_length_max=oligo_generation_parameters["probe_length_max"],
            files_fasta_oligo_database=oligo_generation_parameters["files_fasta_probe_database"],
            min_oligos_per_gene=probe_set_selection_parameters["independent_set_selection"]["set_size_min"],
        )

        if write_intermediate_steps:
            dir_database = oligo_database.save_database(name_database="1_db_probes_initial")
            logger.info(f"Saved probe database for step 1 (Create Database) in directory {dir_database}")

        oligo_database = self._filter_by_property(
            oligo_database=oligo_database,
            isoform_consensus_filter=property_filters_parameters["isoform_consensus_filter"],
            hard_masked_sequences_filter=property_filters_parameters["hard_masked_sequences_filter"],
            soft_masked_sequences_filter=property_filters_parameters["soft_masked_sequences_filter"],
            homopolymeric_runs_filter=property_filters_parameters["homopolymeric_runs_filter"],
            GC_content_filter=property_filters_parameters["GC_content_filter"],
            Tm_filter=property_filters_parameters["Tm_filter"],
            detection_oligo_filter=property_filters_parameters["detection_oligo_filter"],
        )

        if write_intermediate_steps:
            dir_database = oligo_database.save_database(name_database="2_db_probes_property_filter")
            logger.info(f"Saved probe database for step 2 (Property Filters) in directory {dir_database}")

        # Arm Tm/length and ligation_site are needed before seed-region BLAST and backbone assembly.
        padlock_arms_property = PadlockArmsProperty(
            arm_length_min=padlock_arms_parameters["length_min"],
            arm_Tm_dif_max=padlock_arms_parameters["Tm_dif_max"],
            arm_Tm_min=padlock_arms_parameters["Tm_min"],
            arm_Tm_max=padlock_arms_parameters["Tm_max"],
            Tm_parameters=padlock_arms_parameters["Tm_parameters"],
            Tm_chem_correction_parameters=padlock_arms_parameters["Tm_chem_correction_parameters"],
            Tm_salt_correction_parameters=padlock_arms_parameters["Tm_salt_correction_parameters"],
        )
        calculator = PropertyCalculator(properties=[padlock_arms_property])
        oligo_database = calculator.apply(
            oligo_database=oligo_database, sequence_type="oligo", n_jobs=self.n_jobs
        )

        oligo_database = self._filter_by_specificity(
            oligo_database=oligo_database,
            specificity_blastn_filter=specificity_filters_parameters["specificity_blastn_filter"],
            cross_hybridization_blastn_filter=specificity_filters_parameters[
                "cross_hybridization_blastn_filter"
            ],
        )

        if write_intermediate_steps:
            dir_database = oligo_database.save_database(name_database="3_db_probes_specificity_filter")
            logger.info(f"Saved probe database for step 3 (Specificity Filters) in directory {dir_database}")

        oligo_database = self._create_oligo_sets(
            oligo_database=oligo_database,
            independent_set_selection=probe_set_selection_parameters["independent_set_selection"],
            isoform_consensus_score=probe_set_selection_parameters["isoform_consensus_score"],
            GC_content_score=probe_set_selection_parameters["GC_content_score"],
            Tm_score=probe_set_selection_parameters["Tm_score"],
        )

        if write_intermediate_steps:
            dir_database = oligo_database.save_database(name_database="4_db_probes_probesets")
            logger.info(f"Saved probe database for step 4 (Set Selection) in directory {dir_database}.")

        return oligo_database

    @pipeline_step_basic(step_name="Create Database")
    def _create_oligo_database(
        self,
        region_ids: list | None,
        oligo_length_min: int,
        oligo_length_max: int,
        files_fasta_oligo_database: list[str],
        min_oligos_per_gene: int,
    ) -> OligoDatabase:
        """
        Create the first database of candidate SCRINSHOT target probes.

        Candidate probes are generated by sliding windows across the input
        sequences. All probe lengths between the minimum and maximum length are
        considered. For each candidate, the transcript-facing sequence is stored,
        and the reverse complement is stored as the DNA probe sequence that will
        bind to the RNA.

        Regions with too few candidate probes are removed at this stage.

        :param region_ids: Target regions to design probes for, usually gene names
            or gene IDs. If ``None``, all regions in the input FASTA files are used.
        :type region_ids: list[str] | None
        :param oligo_length_min: Minimum candidate probe length in bases.
        :type oligo_length_min: int
        :param oligo_length_max: Maximum candidate probe length in bases.
        :type oligo_length_max: int
        :param files_fasta_oligo_database: FASTA files containing the transcript or
            target-region sequences used for probe design.
        :type files_fasta_oligo_database: list[str]
        :param min_oligos_per_gene: Minimum number of candidate probes a region must
            have to remain in the database.
        :type min_oligos_per_gene: int
        :return: Database containing candidate probes with the target sequence and
            the DNA probe sequence.
        :rtype: OligoDatabase
        """

        oligo_sequences = OligoSequenceGenerator(dir_output=self.dir_output)
        oligo_fasta_file = oligo_sequences.create_sequences_sliding_window(
            files_fasta_in=files_fasta_oligo_database,
            length_interval_sequences=(oligo_length_min, oligo_length_max),
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
            database_overwrite=True,
            sequence_type="target",
            region_ids=region_ids,
        )

        oligo_database.set_database_sequence_types(["target", "oligo"])
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
        hard_masked_sequences_filter: dict,
        soft_masked_sequences_filter: dict,
        homopolymeric_runs_filter: dict,
        GC_content_filter: dict,
        Tm_filter: dict,
        detection_oligo_filter: dict,
    ) -> OligoDatabase:
        """
        Remove candidate probes with unsuitable sequence properties.

        This step checks whether each candidate probe is likely to work well as a
        SCRINSHOT padlock target. It can remove probes that overlap masked sequence,
        contain long single-base runs, have unsuitable GC content, or have a melting
        temperature outside the chosen range.

        The detection-oligo filter is also applied here. It checks whether a
        suitable detection oligo can be placed around the future ligation site and
        whether the two padlock arms can meet the requested length and melting
        temperature criteria.

        Isoform consensus filtering, when enabled, checks how well a probe
        represents the annotated isoforms of the target gene.

        :param oligo_database: Candidate probe database returned by
            :py:meth:`_create_oligo_database`. This database is updated by the
            filtering step.
        :type oligo_database: OligoDatabase
        :param isoform_consensus_filter: Settings for keeping probes that target a
            sufficient fraction of annotated isoforms.
        :type isoform_consensus_filter: dict
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
        :param Tm_filter: Settings for the allowed melting-temperature range and the
            conditions used for the calculation.
        :type Tm_filter: dict
        :param detection_oligo_filter: Settings for checking whether a suitable
            detection oligo and valid padlock arms can be designed for the candidate
            probe.
        :type detection_oligo_filter: dict
        :return: Filtered database containing probes that passed the enabled
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

        # Full-length Tm is separate from DetectionOligoFilter, which also validates padlock arms.
        if Tm_filter["enabled"]:
            melting_temperature = MeltingTemperatureNNFilter(
                Tm_min=Tm_filter["Tm_min"],
                Tm_max=Tm_filter["Tm_max"],
                Tm_parameters=Tm_filter["Tm_parameters"],
                Tm_chem_correction_parameters=Tm_filter["Tm_chem_correction_parameters"],
                Tm_salt_correction_parameters=Tm_filter["Tm_salt_correction_parameters"],
            )
            filters.append(melting_temperature)

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
        specificity_blastn_filter: dict,
        cross_hybridization_blastn_filter: dict,
    ) -> OligoDatabase:
        """
        Remove probes that may bind to the wrong place.

        This step checks whether candidate probes are specific to their intended RNA
        target. It removes exact duplicate probe sequences and, when enabled, uses
        BLASTN to find probes that may also bind to unintended reference sequences.

        For SCRINSHOT, the ligation site is especially important. If a probe has an
        off-target hit that spans the ligation site, the padlock could potentially
        be circularized on the wrong target. When the ligation-region check is
        enabled, the BLASTN filter focuses on this junction region.

        The method can also remove probes that may cross-hybridize with other probes
        in the same panel.

        :param oligo_database: Probe database returned by
            :py:meth:`_filter_by_property`. This database is updated by the
            specificity filters.
        :type oligo_database: OligoDatabase
        :param specificity_blastn_filter: Settings for checking probe specificity
            against reference sequences. This includes the reference FASTA files and
            BLASTN search settings. It can also include a ligation-region setting to
            focus the check around the padlock ligation site.
        :type specificity_blastn_filter: dict
        :param cross_hybridization_blastn_filter: Settings for checking whether
            probes in the same panel may bind to each other or to unintended probe
            targets.
        :type cross_hybridization_blastn_filter: dict
        :return: Filtered database containing probes that passed the enabled
            specificity checks.
        :rtype: OligoDatabase
        """
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
            if specificity_blastn_filter["ligation_region_size"] > 0:
                # Prefer off-targets that span the ligation site; hits that miss the
                # junction cannot circularize a padlock and are less harmful.
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

    @pipeline_step_basic(step_name="Set Selection")
    def _create_oligo_sets(
        self,
        oligo_database: OligoDatabase,
        independent_set_selection: dict,
        isoform_consensus_score: dict,
        GC_content_score: dict,
        Tm_score: dict,
    ) -> OligoDatabase:
        """
        Select final SCRINSHOT target-probe sets.

        This step chooses groups of probes from the filtered candidates. The selected
        probes should be well spaced along the target region and should meet the
        requested number of probes per target.

        Probe sets are scored using isoform consensus, GC content, and melting
        temperature. The method can keep more than one possible probe set per target
        region, which gives users alternatives when several good designs are
        available. Regions without enough suitable probes are removed.

        :param oligo_database: Filtered probe database returned by
            :py:meth:`_filter_by_specificity`. This database is updated with the
            selected probe sets.
        :type oligo_database: OligoDatabase
        :param independent_set_selection: Settings that control how many probe sets
            are selected, how many probes each set should contain, and how far apart
            selected probes should be placed.
        :type independent_set_selection: dict
        :param isoform_consensus_score: Settings for scoring probes by how well they
            represent the annotated isoforms of a target gene.
        :type isoform_consensus_score: dict
        :param GC_content_score: Settings for scoring probes by how close their GC
            content is to the desired value.
        :type GC_content_score: dict
        :param Tm_score: Settings for scoring probes by how close their melting
            temperature is to the desired value.
        :type Tm_score: dict
        :return: Database with selected probe sets attached to each remaining target
            region.
        :rtype: OligoDatabase
        """

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
    Design detection oligos for SCRINSHOT RCA products.

    Detection oligos are fluorescent DNA probes that bind to the rolling-circle
    amplification product of a ligated padlock probe. They are used during
    imaging to make each RCA product visible as a bright spot.

    Each detection oligo is designed around the padlock ligation site. This is
    the sequence junction that is created only after correct padlock ligation.
    The oligo length is chosen so its melting temperature is close to the
    imaging conditions.

    For sequential imaging, detection oligos can include uracil bases. These
    uracils allow the oligo to be cleaved by Uracil DNA Glycosylase, so the
    fluorescent signal can be removed before the next imaging round.

    :param n_jobs: Number of worker processes used for steps that can run in
        parallel. Use ``1`` to run without parallel processing.
    :type n_jobs: int
    """

    def __init__(self, n_jobs: int) -> None:
        """Constructor for the DetectionOligoDesigner class."""
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
        Create detection oligos for all SCRINSHOT target probes.

        This method designs one detection oligo for each selected padlock target
        probe. The detection oligo is placed around the ligation site and selected
        from candidate lengths between ``oligo_length_min`` and
        ``oligo_length_max``.

        For each probe, the candidate whose melting temperature is closest to
        ``Tm_opt`` is chosen. Thymines are then replaced with uracils so the
        detection oligo can be cleaved during sequential imaging. A fluorophore
        marker is added to the final sequence.

        :param oligo_database: Database returned by the SCRINSHOT target-probe
            design step. This database is updated with detection-oligo sequences and
            melting temperatures.
        :type oligo_database: OligoDatabase
        :param oligo_length_min: Minimum allowed detection-oligo length in bases.
        :type oligo_length_min: int
        :param oligo_length_max: Maximum allowed detection-oligo length in bases.
        :type oligo_length_max: int
        :param min_thymines: Minimum number of thymines required in a candidate
            detection oligo. These positions can be converted to uracils.
        :type min_thymines: int
        :param U_distance: Maximum spacing in nucleotides between uracil
            substitutions.
        :type U_distance: int
        :param Tm_opt: Desired melting temperature for the detection oligo in °C.
        :type Tm_opt: float
        :param Tm_parameters: Settings used for melting-temperature calculation.
        :type Tm_parameters: dict
        :param Tm_chem_correction_parameters: Optional chemical correction settings
            for melting-temperature calculation. Use ``None`` when no chemical
            correction is applied.
        :type Tm_chem_correction_parameters: dict | None
        :param Tm_salt_correction_parameters: Optional salt correction settings for
            melting-temperature calculation. Use ``None`` when no salt correction is
            applied.
        :type Tm_salt_correction_parameters: dict | None
        :return: The same database, updated with ``sequence_detection_oligo`` and
            ``Tm_detection_oligo`` for each probe.
        :rtype: OligoDatabase
        """

        region_ids = list(oligo_database.database.keys())

        with joblib_progress(description="Design Detection Oligos", total=len(region_ids)):
            # sharedmem: workers mutate oligo_database in place and return nothing.
            Parallel(n_jobs=self.n_jobs, prefer="threads", require="sharedmem")(
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
        Design detection oligos for one target region.

        This method processes all selected probes for one region, usually one gene.
        For each probe, candidate detection oligos are built around the ligation
        site. The method tries a centered window and two slightly shifted windows,
        then keeps candidates with enough thymines for later uracil substitution.

        The best starting candidate is the one with a melting temperature closest to
        ``Tm_opt``. The candidate is then shortened from either side to test nearby
        lengths. The best final length is selected by melting temperature. After
        that, selected thymines are changed to uracils and a fluorophore marker is
        added.

        The final detection oligo sequence and its melting temperature are written
        back to the database.

        :param oligo_database: Database containing the selected SCRINSHOT target
            probes. This database is updated in place.
        :type oligo_database: OligoDatabase
        :param region_id: Target region to process, usually a gene name or gene ID.
        :type region_id: str
        :param oligo_length_min: Minimum allowed detection-oligo length in bases.
        :type oligo_length_min: int
        :param oligo_length_max: Maximum allowed detection-oligo length in bases.
        :type oligo_length_max: int
        :param min_thymines: Minimum number of thymines required in a candidate
            detection oligo.
        :type min_thymines: int
        :param U_distance: Maximum spacing in nucleotides between uracil
            substitutions.
        :type U_distance: int
        :param Tm_opt: Desired melting temperature for the detection oligo in °C.
        :type Tm_opt: float
        :param Tm_parameters: Settings used for melting-temperature calculation.
        :type Tm_parameters: dict
        :param Tm_chem_correction_parameters: Optional chemical correction settings
            for melting-temperature calculation. Use ``None`` when no chemical
            correction is applied.
        :type Tm_chem_correction_parameters: dict | None
        :param Tm_salt_correction_parameters: Optional salt correction settings for
            melting-temperature calculation. Use ``None`` when no salt correction is
            applied.
        :type Tm_salt_correction_parameters: dict | None
        :return: None. The database is updated in place.
        :rtype: None
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
                if not isinstance(sequence_oligo, str) or not isinstance(ligation_site, int):
                    continue

                # Three windows around the ligation site (centered and ±shifted); keep the
                # closest to Tm_opt, then trim from either end for nearby lengths.
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

                initial_oligos = [
                    detect_oligo
                    for detect_oligo in [
                        detect_oligo_even,
                        detect_oligo_long_left,
                        detect_oligo_long_right,
                    ]
                    if (detect_oligo is not None) and (detect_oligo.count("T") >= min_thymines)
                ]

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

                # Score Tm on the T-only DNA; U substitutions for UNG cleavage come after.
                Tm_detection_oligo = calc_tm_nn(
                    sequence=detection_oligo,
                    Tm_parameters=Tm_parameters,
                    Tm_chem_correction_parameters=Tm_chem_correction_parameters,
                    Tm_salt_correction_parameters=Tm_salt_correction_parameters,
                )

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
        Calculate how far a candidate is from the desired melting temperature.

        The returned value is the absolute difference between the candidate melting
        temperature and ``Tm_opt``. Smaller values mean the candidate is closer to
        the requested imaging conditions.

        :param oligo: Candidate detection-oligo sequence.
        :type oligo: str
        :param Tm_opt: Desired melting temperature in °C.
        :type Tm_opt: float
        :param Tm_parameters: Settings used for melting-temperature calculation.
        :type Tm_parameters: dict
        :param Tm_chem_correction_parameters: Optional chemical correction settings
            for melting-temperature calculation. Use ``None`` when no chemical
            correction is applied.
        :type Tm_chem_correction_parameters: dict | None
        :param Tm_salt_correction_parameters: Optional salt correction settings for
            melting-temperature calculation. Use ``None`` when no salt correction is
            applied.
        :type Tm_salt_correction_parameters: dict | None
        :return: Absolute difference between the candidate Tm and ``Tm_opt`` in °C.
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
        Test shorter versions of a detection-oligo candidate.

        Starting from one candidate sequence, this method removes one nucleotide at
        a time and scores each shortened version by melting temperature. Trimming is
        alternated between the two ends so the ligation site stays close to the
        center of the oligo.

        Only variants that are still long enough and contain enough thymines are
        kept. The thymines are needed later for uracil substitution and enzymatic
        cleavage.

        :param oligo: Starting detection-oligo candidate.
        :type oligo: str
        :param cut_from_right: If ``True``, trimming starts from the right end. If
            ``False``, trimming starts from the left end.
        :type cut_from_right: bool
        :param oligo_length_min: Minimum allowed detection-oligo length in bases.
        :type oligo_length_min: int
        :param min_thymines: Minimum number of thymines required for a candidate to
            be kept.
        :type min_thymines: int
        :param Tm_opt: Desired melting temperature for the detection oligo in °C.
        :type Tm_opt: float
        :param Tm_parameters: Settings used for melting-temperature calculation.
        :type Tm_parameters: dict
        :param Tm_chem_correction_parameters: Optional chemical correction settings
            for melting-temperature calculation. Use ``None`` when no chemical
            correction is applied.
        :type Tm_chem_correction_parameters: dict | None
        :param Tm_salt_correction_parameters: Optional salt correction settings for
            melting-temperature calculation. Use ``None`` when no salt correction is
            applied.
        :type Tm_salt_correction_parameters: dict | None
        :return: Candidate oligos and their absolute Tm differences from
            ``Tm_opt``.
        :rtype: tuple[list[str], list[float]]
        """
        oligos = [oligo]
        Tm_dif = [
            self._get_Tm_dif(
                oligo, Tm_opt, Tm_parameters, Tm_chem_correction_parameters, Tm_salt_correction_parameters
            )
        ]

        # Alternate ends each step so the ligation site stays near the center as length shrinks.
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
        Add uracils and a fluorophore marker to a detection oligo.

        This method replaces selected thymines with uracils. The uracils provide
        cleavage sites for Uracil DNA Glycosylase, which allows the detection signal
        to be removed between imaging rounds.

        Uracils are placed along the sequence with the requested maximum spacing. A
        ``[fluorophore]`` marker is added to one end of the oligo. The marker is
        placed on the side with fewer nearby thymines, keeping it away from the
        main cleavage positions when possible.

        :param oligo: Detection-oligo DNA sequence before uracil substitution.
        :type oligo: str
        :param min_thymines: Number of thymines to replace with uracils.
        :type min_thymines: int
        :param U_distance: Maximum spacing in nucleotides between uracil
            substitutions.
        :type U_distance: int
        :return: Detection oligo with uracil substitutions and a ``[fluorophore]``
            marker.
        :rtype: str
        """
        # Attach the dye on the end with fewer nearby T's so uracil cleavage sites stay away from the label.
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

        if fluorophor_pos == "left":
            oligo = "[fluorophore]" + oligo
        elif fluorophor_pos == "right":
            oligo = oligo[::-1] + "[fluorophore]"

        return oligo


############################################
# SCRINSHOT Designer Pipeline
############################################


def _preprocess_config(config_validated: ScrinshotProbeDesignerConfig) -> dict[str, Any]:
    """
    Prepare the SCRINSHOT config before the pipeline runs.

    This step converts the configuration to a dict and updates the config in place
    so later design stages can read ready-to-use settings.
    It resolves melting-temperature tables for both target probes and the
    detection oligo, turns off unused temperature corrections, and copies the shared
    temperature settings into the filters and scoring steps that need them.

    It also builds the detection-oligo filter from the detection-oligo length settings
    and the padlock-arm constraints, and expands an optional gene-list file into a
    concrete list of target regions. If no gene list is provided, all regions in the
    input FASTA files are used.

    Lastly, it inserts the parameters from required_parameters into the correct sections.

    :param config_validated: Validated pipeline configuration (pydantic model).
    :type config_validated: ScrinshotProbeDesignerConfig
    :return: The configuration converted to a dict, updated with the prepared settings.
    :rtype: dict
    """

    config = config_validated.model_dump()

    apply_required_parameters(config)

    # Resolve Tm table names and blank disabled chem/salt corrections to None so
    # downstream filters treat None as "no correction" without checking the flag.
    for section in ["target_probes", "detection_oligo"]:
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

    # DetectionOligoFilter also scores padlock-arm Tm, so it uses target-probe conditions,
    # not the detection-oligo imaging Tm settings.
    config["target_probes"]["property_filters"]["detection_oligo_filter"] = {
        "oligo_length_min": config["detection_oligo"]["oligo_generation"]["oligo_length_min"],
        "oligo_length_max": config["detection_oligo"]["oligo_generation"]["oligo_length_max"],
        "min_thymines": config["detection_oligo"]["oligo_generation"]["min_thymines"],
        "padlock_arm_length_min": config["target_probes"]["padlock_arms_properties"]["length_min"],
        "padlock_arm_Tm_dif_max": config["target_probes"]["padlock_arms_properties"]["Tm_dif_max"],
        "padlock_arm_Tm_min": config["target_probes"]["padlock_arms_properties"]["Tm_min"],
        "padlock_arm_Tm_max": config["target_probes"]["padlock_arms_properties"]["Tm_max"],
        "Tm_parameters": target_probe_Tm_parameters,
        "Tm_chem_correction_parameters": target_probe_Tm_chem_correction_parameters,
        "Tm_salt_correction_parameters": target_probe_Tm_salt_correction_parameters,
    }

    # Inline shared Tm settings into the blocks that consume them.
    config["target_probes"]["property_filters"]["Tm_filter"]["Tm_parameters"] = target_probe_Tm_parameters
    config["target_probes"]["property_filters"]["Tm_filter"][
        "Tm_chem_correction_parameters"
    ] = target_probe_Tm_chem_correction_parameters
    config["target_probes"]["property_filters"]["Tm_filter"][
        "Tm_salt_correction_parameters"
    ] = target_probe_Tm_salt_correction_parameters

    config["target_probes"]["padlock_arms_properties"]["Tm_parameters"] = target_probe_Tm_parameters
    config["target_probes"]["padlock_arms_properties"][
        "Tm_chem_correction_parameters"
    ] = target_probe_Tm_chem_correction_parameters
    config["target_probes"]["padlock_arms_properties"][
        "Tm_salt_correction_parameters"
    ] = target_probe_Tm_salt_correction_parameters

    config["target_probes"]["probe_set_selection"]["Tm_score"]["Tm_parameters"] = target_probe_Tm_parameters
    config["target_probes"]["probe_set_selection"]["Tm_score"][
        "Tm_chem_correction_parameters"
    ] = target_probe_Tm_chem_correction_parameters
    config["target_probes"]["probe_set_selection"]["Tm_score"][
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


def scrinshot_probe_designer(config: ScrinshotProbeDesignerConfig) -> None:
    """
    Run the SCRINSHOT probe design pipeline from a config dict.

    This function prepares the config with :func:`_preprocess_config`, then runs
    :class:`ScrinshotProbeDesigner` end to end. It designs target probes and
    detection oligos, assembles padlock probes, and writes the final files under
    ``config['general']['dir_output']``. The caller should configure the library
    logger before calling this function (see :func:`main`).

    The config should follow ``data/configs/scrinshot_probe_designer.yaml``.

    Top-level config sections:

    - ``general``: output directory, intermediate-step writing, and worker count.
    - ``target_probes``: candidate generation, sequence filters, specificity filters,
      probe set selection, and padlock-arm settings.
    - ``detection_oligo``: detection-oligo generation and temperature settings.

    Files written under ``dir_output``:

    - ``padlock_probes.yml``: full probe records.
    - ``padlock_probes_order.yml``: sequences ready for synthesis.
    - ``padlock_probes.tsv`` / ``padlock_probes.xlsx``: probe sets as tables.

    Intermediate probe databases are also written when
    ``general.write_intermediate_steps`` is ``True``.

    See :class:`ScrinshotProbeDesigner` for the pipeline description and probe
    structure.

    :param config: Validated pipeline configuration. It is converted to a dict and
        updated in place by :func:`_preprocess_config` before the pipeline runs.
    :type config: ScrinshotProbeDesignerConfig
    :return: None
    :rtype: None
    """

    config_dict = _preprocess_config(config)

    pipeline = ScrinshotProbeDesigner(
        write_intermediate_steps=config_dict["general"]["write_intermediate_steps"],
        dir_output=config_dict["general"]["dir_output"],
        n_jobs=config_dict["general"]["n_jobs"],
    )

    target_probe_database = pipeline.design_target_probes(
        target_probes_parameters=config_dict["target_probes"],
    )

    target_probe_database = pipeline.design_detection_oligos(
        oligo_database=target_probe_database,
        oligo_generation_parameters=config_dict["detection_oligo"]["oligo_generation"],
    )

    target_probe_database = pipeline.assemble_padlock_backbone(
        oligo_database=target_probe_database,
        padlock_arms_parameters=config_dict["target_probes"]["padlock_arms_properties"],
    )

    pipeline.generate_output(oligo_database=target_probe_database)


def main() -> None:
    """
    Run the SCRINSHOT probe design pipeline from the command line.

    Parses the required ``-c``/``--config`` argument, loads the YAML configuration
    file, and configures the library logger to write under the configured output
    directory. It then calls :func:`scrinshot_probe_designer`.

    :return: None
    :rtype: None
    """
    print("--------------START PIPELINE--------------")

    args = base_parser(
        prog="SCRINSHOT Probe Designer",
        usage="scrinshot_probe_designer [options]",
        description=__doc__,
    )

    config_validated = load_config(args["config"], ScrinshotProbeDesignerConfig)

    # Configure logging only after dir_output is known so the log file lands there.
    configure_root_logger(
        dir_output=config_validated.general.dir_output,
        pipeline_name="scrinshot_probe_designer",
    )

    scrinshot_probe_designer(config_validated)

    print("--------------END PIPELINE--------------")


if __name__ == "__main__":
    main()
