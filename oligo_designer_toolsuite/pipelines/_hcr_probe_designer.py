"""
HCR probe designer pipeline.

HCR, or Hybridization Chain Reaction, is an RNA-FISH method that uses pairs of
split-initiator probes to trigger fluorescent signal amplification at RNA
targets. Signal is produced by HCR hairpins after both probe halves bind next to
each other on the transcript.

See :class:`HcrProbeDesigner` for the full pipeline description and probe
structure. See :func:`hcr_probe_designer` for the config-driven workflow.
"""

############################################
# imports
############################################

import os
import shutil
from pathlib import Path
from typing import Any

import pandas as pd
import yaml
from pydantic import ValidationError

from oligo_designer_toolsuite._exceptions import (
    FeatureNotImplementedError,
    FileFormatError,
)
from oligo_designer_toolsuite.config.pipelines.hcr_probe_designer import HcrProbeDesignerConfig
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
    """
    Design probe sets for HCR RNA-FISH experiments.

    This class runs the design workflow for split-initiator HCR probes. The
    final output is a set of left and right probe sequences for each target gene,
    together with the assigned HCR amplifier information.

    Overview
    --------
    HCR, or Hybridization Chain Reaction, is an RNA-FISH method that amplifies
    fluorescent signal without enzymes. It is used to detect RNA molecules in
    fixed cells, tissues, embryos, and whole-mount samples.

    In HCR RNA-FISH, each RNA target is detected by pairs of DNA probes. The two
    probes in a pair bind next to each other on the same transcript. Each probe
    carries one half of an HCR initiator. When both probes bind at the correct
    position, the two initiator halves form a complete initiator and start HCR
    amplification.

    The amplification step uses two fluorescent DNA hairpins. Once triggered,
    the hairpins open one after the other and form a fluorescent polymer at the
    RNA target. A single probe that binds somewhere else should not start
    amplification on its own, which helps reduce background signal.

    Probe Structure
    ---------------
    **Hybridization (Target) Probes**

    Each target gene is detected by several probe pairs. A pair consists of a
    left probe and a right probe that bind to nearby positions on the RNA
    transcript, usually with a small gap between them.

    Each probe contains:

    - a target-binding sequence that is complementary to the RNA,
    - one half of an HCR initiator,
    - a short linker between the target-binding sequence and the initiator half.

    The left and right probes are arranged so that the initiator halves face
    each other when both probes bind to the RNA.

    A simplified layout is::

        Left probe:
            [initiator half L] + [linker] + [target-binding sequence L]

        Right probe:
            [target-binding sequence R] + [linker] + [initiator half R]

    Only the full probe pair should trigger HCR amplification. This split design
    is important because it makes accidental signal from single off-target probes
    much less likely.

    **HCR Amplifiers**

    Each HCR amplifier consists of two DNA hairpins, usually called H1 and H2.
    The hairpins are stable until they meet the matching initiator. After
    initiation, H1 and H2 open in turn and build a fluorescent polymer at the
    target site.

    Different amplifiers, such as B1, B2, B3, and others, can be used in the
    same experiment. Each target gene is assigned to one amplifier. The chosen
    amplifier determines which initiator sequences are attached to the probes
    and which fluorescent hairpins are used during imaging.

    **Codebook and Initiator Table**

    The codebook assigns each target gene to an HCR amplifier. In this pipeline,
    the codebook is a one-hot table: each row is a gene, each column is an
    amplifier bit, and each gene has exactly one active bit.

    This means that each gene is detected in one fluorescence channel. Standard
    HCR is therefore not a combinatorial barcoding method. The number of genes
    that can be imaged together depends on the number of available amplifiers
    and fluorophores in the experiment.

    The initiator table links each amplifier bit to the left and right
    half-initiator sequences used to build the probe pairs.

    Probe Library Preparation
    -------------------------
    Standard HCR probe sets are often small enough to order as individual
    single-stranded DNA oligos. A typical target uses several probe pairs, often
    around 10 to 30 pairs per gene, depending on transcript length and design
    choices.

    The output of this pipeline is a sequence table that can be checked, shared,
    and submitted for oligo synthesis. The synthesized probes are then pooled by
    target or by experiment before hybridization.

    During the experiment, the probe set is hybridized to the sample. After
    washing, matching fluorescent HCR hairpins are added. Hairpin amplification
    creates the signal at the RNA target.

    Pipeline Overview
    -----------------
    The pipeline performs the main steps needed to design an HCR probe set:

    1. **Target probe design**

       Design left and right target-binding sequences for each gene. The two
       sequences in a pair are selected so they bind close to each other on the
       transcript.

    2. **Initiator assignment**

       Load or generate a codebook that assigns each gene to one HCR amplifier.
       Load the initiator table that provides the matching left and right
       half-initiator sequences.

    3. **Hybridization probe assembly**

       Combine each target-binding sequence with the assigned initiator half and
       linker sequence to build the final left and right HCR probes.

    4. **Output generation**

       Write the final probe sequences, probe properties, codebook, and
       initiator table to files that can be inspected and used for ordering.

    References
    ----------
    Choi, H. M. T., Schwarzkopf, M., Fornace, M. E., Acharya, A.,
    Artavanis, G., Stegmaier, J., Cunha, A., & Pierce, N. A. (2018).
    Third-generation in situ hybridization chain reaction: multiplexed,
    quantitative, sensitive, versatile, robust. Development, 145(12),
    dev165753. https://doi.org/10.1242/dev.165753

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

    def __init__(
        self,
        write_intermediate_steps: bool,
        dir_output: str,
        n_jobs: int,
    ) -> None:
        """Constructor for the HcrProbeDesigner class."""
        self.dir_output = os.path.abspath(dir_output)
        Path(self.dir_output).mkdir(parents=True, exist_ok=True)

        self.write_intermediate_steps = write_intermediate_steps
        self.n_jobs = n_jobs

    def design_target_probes(
        self,
        target_probes_parameters: dict,
    ) -> OligoDatabase:
        """
        Design the RNA-binding parts of the HCR probes.

        This step designs left and right target-binding sequences for each target
        gene. The two halves are selected so they bind close to each other on the RNA
        transcript. This is needed for split-initiator HCR, because amplification
        should only start when both probe halves bind at the same target site.

        Candidate probes are generated, filtered for sequence quality and
        specificity, and then selected into final probe sets. The melting
        temperature is calculated separately for the left and right halves and saved
        with each probe.

        :param target_probes_parameters: Settings from the ``target_probes`` section
            of the pipeline config. This includes candidate generation, sequence
            filters, specificity filters, and probe set selection.
        :type target_probes_parameters: dict
        :return: Database containing the selected left and right target-binding
            probe halves for each target gene.
        :rtype: OligoDatabase
        """
        target_probe_designer = TargetProbeDesigner(self.dir_output, self.n_jobs)
        target_probes_database = target_probe_designer.generate_target_probes(
            target_probes_parameters=target_probes_parameters,
            write_intermediate_steps=self.write_intermediate_steps,
        )

        # Per-arm Tm: each half is scored independently for the output tables
        # (params from Tm_filter). High Tm is still preferred, though stripping is less critical than CycleHCR.
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
        """
        Load or create the HCR codebook and initiator table.

        The codebook assigns each target gene to one HCR amplifier. In standard HCR,
        this is a one-hot assignment: each gene has one active bit and is detected in
        one fluorescence channel.

        The initiator table links each codebook bit to the left and right
        half-initiator sequences for the matching amplifier. These half-initiators
        are attached to the left and right target probes during probe assembly.

        The codebook and initiator table can either be loaded from files or
        generated from the config. The initiator table is checked as soon as it
        is available. After the codebook is ready, both tables are checked
        together, so missing targets or missing initiator sequences are caught
        before probe assembly.

        :param region_ids: Target regions that must be represented in the codebook,
            usually gene names or gene IDs.
        :type region_ids: list[str]
        :param initiator_probe_parameters: Settings from the ``initiator_probes``
            section of the pipeline config. This includes the codebook settings and
            the initiator table settings.
        :type initiator_probe_parameters: dict
        :return: Codebook and initiator table used to assign HCR amplifiers to
            target genes.
        :rtype: tuple[pd.DataFrame, pd.DataFrame]
        """
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

        # Table-only validation: catch structural issues (index name, columns, DNA validity)
        # as soon as the table is available.
        initiator_designer.validate(
            initiator_table=initiator_table,
            initiator_table_source=initiator_table_source,
        )

        if initiator_probe_parameters["codebook"]["source"] == "load":
            codebook = initiator_designer.load_codebook(
                file_codebook=initiator_probe_parameters["codebook"]["file"]
            )
            codebook_source = initiator_probe_parameters["codebook"]["file"]
        else:
            codebook = initiator_designer.generate_codebook(region_ids=region_ids)
            codebook_source = initiator_probe_parameters["codebook"]["source"]

        # Full pair validation: also checks that every codebook bit is covered by the table.
        initiator_designer.validate(
            initiator_table=initiator_table,
            initiator_table_source=initiator_table_source,
            codebook=codebook,
            region_ids=region_ids,
            codebook_source=codebook_source,
        )

        return codebook, initiator_table

    def assemble_hybridization_probes(
        self,
        oligo_database: OligoDatabase,
        hybridization_probe_parameters: dict,
        codebook: pd.DataFrame,
        initiator_table: pd.DataFrame,
    ) -> OligoDatabase:
        """
        Build the final HCR probe sequences.

        This step combines each target-binding probe half with the assigned
        half-initiator and linker sequence. Each target site receives a left and a
        right HCR probe. The target-binding parts bind to the RNA, while the
        initiator halves face outward and remain available to start HCR
        amplification.

        A simplified layout is::

            Left probe:
                [initiator half L] + [linker] + [target-binding sequence L]

            Right probe:
                [target-binding sequence R] + [linker] + [initiator half R]

        The assembled sequences are added to the existing probe database.

        :param oligo_database: Database returned by :py:meth:`design_target_probes`.
            This database is updated with the assembled HCR probe sequences.
        :type oligo_database: OligoDatabase
        :param hybridization_probe_parameters: Settings from the
            ``hybridization_probes`` section of the pipeline config, including the
            linker sequence.
        :type hybridization_probe_parameters: dict
        :param codebook: Table returned by :py:meth:`design_initiators`. Rows are
            target regions and columns are barcode bits.
        :type codebook: pd.DataFrame
        :param initiator_table: Table returned by :py:meth:`design_initiators`.
            Links each amplifier bit to its left and right half-initiator sequences.
        :type initiator_table: pd.DataFrame
        :return: Database with left and right HCR probe sequences added.
        :rtype: OligoDatabase
        """
        linker_sequence = hybridization_probe_parameters["linker_sequence"]

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
            # Weight-1 (one-hot) codebook: the single active bit supplies both L and R initiator halves.
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
        """
        Write the completed HCR probe design to files.

        This step saves the final probe database, the codebook, and the initiator
        table. It also writes an order-ready file with the left and right HCR probe
        sequences and the initiator sequences needed for synthesis or checking.

        If no output properties are provided, a default set of annotations and
        sequence fields is written.

        :param oligo_database: Database returned by
            :py:meth:`assemble_hybridization_probes`.
        :type oligo_database: OligoDatabase
        :param codebook: Table assigning each target gene to one HCR amplifier bit.
            Rows are target regions and columns are barcode bits.
        :type codebook: pd.DataFrame
        :param initiator_table: Table linking each amplifier bit to its left and
            right half-initiator sequences.
        :type initiator_table: pd.DataFrame
        :param output_properties: Probe properties to include in the detailed output
            files. If ``None``, a default set of annotations and sequences is used.
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

        codebook.to_csv(os.path.join(self.dir_output, "codebook.tsv"), sep="\t", index_label="gene_name")
        initiator_table.to_csv(os.path.join(self.dir_output, "initiators.tsv"), sep="\t")

        oligo_database.write_oligosets_to_yaml(
            properties=output_properties,
            ascending=True,
            filename="hcr_probes",
        )

        oligo_database.write_oligosets_to_table(
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


############################################
# HCR Target Probe Designer
############################################


class TargetProbeDesigner:
    """
    Design the RNA-binding probe pairs used in HCR experiments.

    This class designs the gene-specific part of HCR probes. Each candidate is
    split into a left and a right target-binding arm. The two arms are selected
    so they bind close to each other on the RNA transcript.

    Both arms must have suitable sequence properties and should bind only to the
    intended target. This matters because HCR amplification should start only
    when both probe halves bind next to each other on the same RNA molecule.

    The workflow has four main steps:

    1. **Candidate generation**

       Build candidate probes from transcript FASTA files and split each
       candidate into a left arm, a gap sequence, and a right arm.

    2. **Sequence filtering**

       Remove candidates with unsuitable sequence properties, such as poor GC
       content, long single-base runs, masked sequence, unsuitable melting
       temperature, or strong secondary structure.

    3. **Specificity filtering**

       Remove candidates that are likely to bind to unintended transcripts or to
       other probes in the panel.

    4. **Probe set selection**

       Select suitable probe sets for each target gene, while keeping probes
       well spaced across the transcript.

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
        self.subdir_db_oligos = "db_target_probes"
        self.subdir_db_reference = "db_reference"

        self.n_jobs = n_jobs

    def generate_target_probes(
        self,
        target_probes_parameters: dict,
        write_intermediate_steps: bool = False,
    ) -> OligoDatabase:
        """
        Run the full HCR target-probe design workflow.

        This method designs the left and right RNA-binding arms used in HCR probe
        pairs. It starts from transcript sequences, creates candidate probe pairs,
        filters them, checks their specificity, and selects final probe sets for
        each target gene.

        Each surviving probe contains two target-binding arms that bind close to
        each other on the RNA. Both arms must pass the same quality checks.

        :param target_probes_parameters: Settings from the ``target_probes`` section
            of the pipeline config. This includes candidate generation, sequence
            filters, specificity filters, and probe set selection.
        :type target_probes_parameters: dict
        :param write_intermediate_steps: If ``True``, save intermediate probe
            databases after each main step. This can help when checking where probes
            were removed.
        :type write_intermediate_steps: bool
        :return: Database containing the selected target-probe pairs for each target
            gene.
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
        """
        Create the first database of candidate target probes.

        Candidate probes are generated by sliding a fixed-size window across the
        input transcript sequences. Each candidate covers the full target region for
        one HCR probe pair: the left binding site, the gap, and the right binding
        site.

        The candidate is then converted into the DNA probe strand and split into the
        left arm, spacer, and right arm. Regions with too few candidate probes are
        removed at this stage.

        :param region_ids: Target regions to design probes for, usually gene names
            or gene IDs. If ``None``, all regions in the input FASTA files are used.
        :type region_ids: list[str] | None
        :param oligo_length: Total length of the candidate target window in bases.
            This should match the left arm length, gap length, and right arm length
            combined.
        :type oligo_length: int
        :param L_probe_sequence_length: Length of the left target-binding arm in
            nucleotides.
        :type L_probe_sequence_length: int
        :param gap_sequence_length: Length of the gap between the two arms on the RNA
            target, in nucleotides.
        :type gap_sequence_length: int
        :param R_probe_sequence_length: Length of the right target-binding arm in
            nucleotides.
        :type R_probe_sequence_length: int
        :param files_fasta_oligo_database: FASTA files containing the transcript or
            target-region sequences used for probe design.
        :type files_fasta_oligo_database: list[str]
        :param min_oligos_per_gene: Minimum number of candidate probes a region must
            have to remain in the database.
        :type min_oligos_per_gene: int
        :return: Database containing candidate probes with target, probe, left-arm,
            right-arm, and spacer sequences.
        :rtype: OligoDatabase
        """
        oligo_sequences = OligoSequenceGenerator(dir_output=self.dir_output)
        oligo_fasta_file = oligo_sequences.create_sequences_sliding_window(
            files_fasta_in=files_fasta_oligo_database,
            length_interval_sequences=(oligo_length, oligo_length),
            region_ids=region_ids,
            n_jobs=self.n_jobs,
        )

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
        oligo_database.set_database_sequence_types(["target", "oligo", "oligo_L", "oligo_R"])

        # Probe strand is the reverse complement of the transcript ("target") window.
        reverse_complement_sequence_property: BaseProperty = ReverseComplementSequenceProperty(
            sequence_type_reverse_complement="oligo"
        )
        calculator = PropertyCalculator(properties=[reverse_complement_sequence_property])
        oligo_database = calculator.apply(
            oligo_database=oligo_database, sequence_type="target", n_jobs=self.n_jobs
        )

        # Split the probe strand, not the target: 5'→3' on the oligo is R, spacer, then L.
        split_start_end = [
            (0, L_probe_sequence_length),
            (L_probe_sequence_length, L_probe_sequence_length + gap_sequence_length),
            (
                L_probe_sequence_length + gap_sequence_length,
                L_probe_sequence_length + gap_sequence_length + R_probe_sequence_length,
            ),
        ]
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
        """
        Remove candidate probes with unsuitable sequence properties.

        This step checks whether each candidate probe is likely to behave well in
        the experiment. It can remove probes that overlap masked sequence, contain
        long single-base runs, have unsuitable GC content, have a melting temperature
        outside the chosen range, or are predicted to fold strongly onto themselves.

        The left and right arms are checked separately. A probe pair is kept only if
        both arms pass the enabled filters. If isoform consensus filtering is
        enabled, probes are also checked for how well they represent the annotated
        isoforms of the target gene.

        :param oligo_database: Candidate probe database returned by
            :py:meth:`_create_oligo_database`. This database is updated by the
            filtering step.
        :type oligo_database: OligoDatabase
        :param isoform_consensus_filter: Settings for keeping probes that target a
            sufficient fraction of annotated isoforms.
        :type isoform_consensus_filter: dict
        :param soft_masked_sequences_filter: Settings for removing probes that
            overlap soft-masked sequence, often used for repetitive or low-complexity
            regions.
        :type soft_masked_sequences_filter: dict
        :param hard_masked_sequences_filter: Settings for removing probes that
            overlap hard-masked bases, such as ``N`` bases.
        :type hard_masked_sequences_filter: dict
        :param homopolymeric_runs_filter: Settings for removing probes with long
            runs of the same base.
        :type homopolymeric_runs_filter: dict
        :param GC_content_filter: Settings for the allowed GC-content range.
        :type GC_content_filter: dict
        :param Tm_filter: Settings for the allowed melting-temperature range and the
            conditions used for the calculation.
        :type Tm_filter: dict
        :param secondary_structure_filter: Settings for removing probes predicted to
            form stable self-structures.
        :type secondary_structure_filter: dict
        :return: Filtered database in which both arms of each remaining probe passed
            the enabled sequence checks.
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

        # Filters were queued cheapest-first so failing probes exit before thermodynamics.
        property_filter = PropertyFilter(filters=filters)

        # Both arms must pass independently; a weak L or R half fails the whole pair.
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
        """
        Remove probes that may bind to the wrong place.

        This step checks whether candidate probes are specific to their intended
        target. It removes exact duplicate matches and, when enabled, uses BLASTN to
        find probes that may also bind to other transcript or reference sequences.

        It can also check whether probe arms are likely to bind to other probes in
        the same panel. For HCR, this is done separately for the left and right arms,
        because both arms need to behave well on their own.

        :param oligo_database: Probe database returned by
            :py:meth:`_filter_by_property`. This database is updated by the
            specificity filters.
        :type oligo_database: OligoDatabase
        :param specificity_blastn_filter: Settings for checking probe specificity
            against reference sequences. This includes the reference FASTA files and
            BLASTN search settings.
        :type specificity_blastn_filter: dict
        :param cross_hybridization_blastn_filter: Settings for checking whether
            probes in the same panel may bind to each other or to the wrong target
            probe arm.
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
            if specificity_blastn_filter["junction_region_size"] > 0:
                # Prefer off-targets that span the L/R junction; single-arm hits cannot
                # bring both initiator halves together and are less harmful.
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

        # Cross-hyb is arm-specific (L vs L, R vs R), so it cannot ride on "oligo".
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
        """
        Select final probe sets for each target gene.

        This step chooses groups of probes from the filtered candidates. The selected
        probes should be well spaced along the transcript and should meet the
        requested number of probes per gene.

        Probe sets are scored by how well they represent the annotated isoforms of
        the target gene. The method can keep more than one possible probe set per
        gene, which gives users alternatives when several good designs are
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
        :return: Database with selected probe sets attached to each remaining target
            gene.
        :rtype: OligoDatabase
        """
        isoform_consensus_scorer = IsoformConsensusScorer(score_weight=isoform_consensus_score["weight"])
        oligos_scoring = OligoScoring(scorers=[isoform_consensus_scorer])
        # ascending=False: higher aggregate (isoform) scores win.
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
    """
    Manage the HCR codebook and initiator table.

    In HCR, each target gene is assigned to one HCR amplifier. Each amplifier is
    represented by one bit in the codebook. A gene has exactly one active bit,
    which means it is detected with one amplifier and one fluorescence channel.

    The initiator table links each bit to the two half-initiator sequences used
    for that amplifier. During probe assembly, these half-initiators are added
    to the left and right target-binding probes. When both probes bind next to
    each other on the RNA, the two halves form the initiator that starts HCR
    amplification.

    This class loads, generates, and checks the two tables needed for this
    assignment:

    - the **codebook**, which assigns each gene to one amplifier bit,
    - the **initiator table**, which gives the left and right half-initiator
      sequences for each bit.

    Standard HCR uses a one-hot codebook. It is not a combinatorial barcode and
    does not provide distance-based error correction. The number of genes that
    can be imaged together is limited by the available amplifiers and
    fluorophores.

    :param dir_output: Directory where output files and intermediate results are
        saved.
    :type dir_output: str
    :param n_jobs: Number of worker processes used for steps that can run in
        parallel. This value is kept for consistency with the other designer
        classes.
    :type n_jobs: int
    """

    def __init__(
        self,
        dir_output: str,
        n_jobs: int,
    ) -> None:
        """Constructor for the InitiatorDesigner class."""
        self.dir_output = os.path.abspath(dir_output)
        self.n_jobs = n_jobs

    def load_codebook(self, file_codebook: str) -> pd.DataFrame:
        """
        Load an HCR codebook from a table file.

        The codebook assigns each target gene to one HCR amplifier. The file should
        contain one row per target region and one column per amplifier bit. The
        target region column must be named ``gene_name``. Bit columns are expected to
        be named ``bit_1``, ``bit_2``, and so on.

        The table is loaded here, but the full content check is done later together
        with the initiator table. This makes sure that every active bit in the
        codebook has matching initiator sequences.

        :param file_codebook: Path to the codebook file. The file must contain a
            ``gene_name`` column and one or more ``bit_*`` columns.
        :type file_codebook: str
        :return: Codebook table indexed by ``gene_name``.
        :rtype: pandas.DataFrame
        """
        return pd.read_csv(file_codebook, sep=None, engine="python", index_col="gene_name")

    def generate_codebook(self, region_ids: list[str]) -> pd.DataFrame:
        """
        Generate an HCR codebook.

        Automatic codebook generation is not available yet. For now, the codebook
        must be supplied as an input file.

        :param region_ids: Target regions that need amplifier assignments, usually
            gene names or gene IDs.
        :type region_ids: list[str]
        :return: Codebook table with one active bit per target region.
        :rtype: pandas.DataFrame
        :raises FeatureNotImplementedError: Always raised until codebook generation
            is implemented.
        """
        raise FeatureNotImplementedError(
            "Generation of codebook is not yet implemented. "
            "Please provide a codebook.file parameter and set codebook.source to 'load'."
        )

    def load_initiator_table(self, file_initiator_table: str) -> pd.DataFrame:
        """
        Load an HCR initiator table from a table file.

        The initiator table links each amplifier bit to the two DNA sequences that
        form the split HCR initiator. One sequence is used on the left probe and the
        other on the right probe.

        The file must contain a ``bit`` column, an ``initiator_L_sequence`` column,
        and an ``initiator_R_sequence`` column. Column and sequence checks are done
        later by :py:meth:`validate`.

        :param file_initiator_table: Path to the initiator table file. The file must
            contain ``bit``, ``initiator_L_sequence``, and
            ``initiator_R_sequence`` columns.
        :type file_initiator_table: str
        :return: Initiator table indexed by ``bit``.
        :rtype: pandas.DataFrame
        :raises FileFormatError: If the ``bit`` column is missing.
        """
        initiator_table = pd.read_csv(file_initiator_table, sep=None, engine="python")
        if "bit" not in initiator_table.columns:
            raise FileFormatError(f"Initiator table '{file_initiator_table}' must contain a 'bit' column.")
        return initiator_table.set_index("bit")

    def generate_initiator_table(self) -> pd.DataFrame:
        """
        Generate an HCR initiator table.

        Automatic initiator table generation is not available yet. For now, the
        initiator table must be supplied as an input file.

        :return: Initiator table indexed by ``bit``.
        :rtype: pandas.DataFrame
        :raises FeatureNotImplementedError: Always raised until initiator table
            generation is implemented.
        """
        raise FeatureNotImplementedError(
            "Generation of initiator table is not yet implemented. "
            "Please provide an initiator_table.file parameter and set initiator_table.source to 'load'."
        )

    def validate(
        self,
        initiator_table: pd.DataFrame | None = None,
        *,
        initiator_table_source: str | None = None,
        codebook: pd.DataFrame | None = None,
        region_ids: list[str] | None = None,
        codebook_source: str | None = None,
    ) -> None:
        """
        Check the codebook and/or initiator table.

        This method can check either table on its own, or both together. When a
        codebook is provided, it must contain all requested target regions, use
        ``gene_name`` as the row index, and contain only ``0`` and ``1`` values in
        its bit columns. For HCR, each target region must have exactly one active
        bit.

        When an initiator table is provided, it must include a valid left and right
        half-initiator sequence for each bit. If a codebook is also provided, every
        bit used by the codebook must appear in the table.

        Running these checks before probe assembly helps catch mismatched input
        files early, such as a codebook that refers to an amplifier bit missing
        from the initiator table.

        :param initiator_table: Optional table linking each amplifier bit to its
            left and right half-initiator sequences. If ``None``, the table is not
            checked.
        :type initiator_table: pd.DataFrame | None
        :param initiator_table_source: File path or source label for the initiator
            table. Used in error messages when the table is checked.
        :type initiator_table_source: str | None
        :param codebook: Optional codebook table assigning target regions to
            amplifier bits. Rows are target regions and columns are ``bit_*``
            entries. If ``None``, the codebook is not checked.
        :type codebook: pd.DataFrame | None
        :param region_ids: Target regions that must be present in the codebook.
            Required when a codebook is checked.
        :type region_ids: list[str] | None
        :param codebook_source: File path or source label for the codebook. Used
            in error messages when the codebook is checked.
        :type codebook_source: str | None
        :return: None
        :rtype: None
        :raises ValueError: If a table is checked without its required companion
            arguments (for example ``region_ids`` / ``codebook_source`` with a
            codebook, or ``initiator_table_source`` with an initiator table).
        :raises FileFormatError: If the codebook or initiator table is missing
            required information or contains invalid values.
        """
        if codebook is not None:
            if region_ids is None:
                raise ValueError("region_ids must be provided when validating a codebook.")
            if codebook_source is None:
                raise ValueError("codebook_source must be provided when validating a codebook.")
            validate_codebook(
                codebook=codebook,
                region_ids=region_ids,
                source=codebook_source,
                expected_hamming_weight=1,
                index_name="gene_name",
            )
        if initiator_table is not None:
            if initiator_table_source is None:
                raise ValueError(
                    "initiator_table_source must be provided when validating an initiator table."
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


def _preprocess_config(config_validated: HcrProbeDesignerConfig) -> dict[str, Any]:
    """
    Prepare the HCR config before the pipeline runs.

    This step converts the configuration to a dict and updates the config in place
    so later design stages can read ready-to-use settings.
    It resolves melting-temperature tables, turns off unused temperature
    corrections, and copies the shared temperature settings into the filters that
    need them.

    It also derives the full probe length and the left/right junction position from
    the arm and gap lengths, and expands an optional gene-list file into a concrete
    list of target regions. If no gene list is provided, all regions in the input
    FASTA files are used.

    :param config_validated: Validated pipeline configuration.
    :type config_validated: HcrProbeDesignerConfig
    :return: The configuration converted to a dict, updated with the prepared settings.
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

    L_probe_sequence_length = config["target_probes"]["oligo_generation"]["L_probe_sequence_length"]
    gap_sequence_length = config["target_probes"]["oligo_generation"]["gap_sequence_length"]
    R_probe_sequence_length = config["target_probes"]["oligo_generation"]["R_probe_sequence_length"]
    config["target_probes"]["oligo_generation"]["oligo_length"] = (
        L_probe_sequence_length + gap_sequence_length + R_probe_sequence_length
    )
    # Junction site is mid-gap on the target strand; the seed-region BLAST filter
    # uses it to bias hits toward junction-spanning off-targets.
    config["target_probes"]["specificity_filters"]["specificity_blastn_filter"]["junction_site"] = (
        L_probe_sequence_length + gap_sequence_length // 2
    )

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


def hcr_probe_designer(config: HcrProbeDesignerConfig) -> None:
    """
    Run the HCR probe design pipeline from a validated configuration (pydantic model).

    This function prepares the config with :func:`_preprocess_config`, then runs
    :class:`HcrProbeDesigner` end to end. It designs target probes, loads or creates
    the initiators and codebook, assembles the hybridization probes, and writes the
    final files under ``config['general']['dir_output']``. The caller should
    configure the library logger before calling this function (see :func:`main`).

    The config should follow ``data/configs/hcr_probe_designer.yaml``.

    Top-level config sections:

    - ``general``: output directory, intermediate-step writing, and worker count.
    - ``target_probes``: candidate generation, sequence filters, specificity filters,
      and probe set selection.
    - ``initiator_probes``: codebook and initiator table settings.
    - ``hybridization_probes``: linker sequence used during probe assembly.

    Files written under ``dir_output``:

    - ``codebook.tsv``: barcode assignments for each target gene.
    - ``initiators.tsv``: initiator sequences and related bit information.
    - ``hcr_probes.yml``: full probe records.
    - ``hcr_probes_order.yml``: sequences ready for synthesis.
    - ``hcr_probes.tsv`` / ``hcr_probes.xlsx``: probe sets as tables.

    Intermediate probe databases are also written when
    ``general.write_intermediate_steps`` is ``True``.

    See :class:`HcrProbeDesigner` for the pipeline description and probe
    structure.

    :param config: Pipeline configuration loaded from the YAML config file. It is
        updated in place by :func:`_preprocess_config` before the pipeline runs.
    :type config: dict
    :return: None
    :rtype: None
    """

    config_dict = _preprocess_config(config)

    pipeline = HcrProbeDesigner(
        write_intermediate_steps=config_dict["general"]["write_intermediate_steps"],
        dir_output=config_dict["general"]["dir_output"],
        n_jobs=config_dict["general"]["n_jobs"],
    )

    target_probe_database = pipeline.design_target_probes(
        target_probes_parameters=config_dict["target_probes"],
    )

    codebook, initiator_table = pipeline.design_initiators(
        region_ids=list(target_probe_database.database.keys()),
        initiator_probe_parameters=config_dict["initiator_probes"],
    )

    hybridization_probe_database = pipeline.assemble_hybridization_probes(
        oligo_database=target_probe_database,
        hybridization_probe_parameters=config_dict["hybridization_probes"],
        codebook=codebook,
        initiator_table=initiator_table,
    )

    pipeline.generate_output(
        oligo_database=hybridization_probe_database,
        codebook=codebook,
        initiator_table=initiator_table,
    )


def main() -> None:
    """
    Run the HCR probe design pipeline from the command line.

    Parses the required ``-c``/``--config`` argument, loads the YAML configuration
    file, and configures the library logger to write under the configured output
    directory. It then calls :func:`hcr_probe_designer`.

    :return: None
    :rtype: None
    """
    print("--------------START PIPELINE--------------")

    args = base_parser(
        prog="HCR Probe Designer",
        usage="hcr_probe_designer [options]",
        description=__doc__,
    )

    with open(args["config"], "r") as handle:
        config_raw = yaml.safe_load(handle)

    try:
        config_validated = HcrProbeDesignerConfig.model_validate(config_raw)
    except ValidationError as e:
        print(f"Invalid configuration file:\n{e}")
        raise

    # Configure logging only after dir_output is known so the log file lands there.
    configure_root_logger(
        dir_output=config_validated.general.dir_output,
        pipeline_name="hcr_probe_designer",
    )

    hcr_probe_designer(config_validated)

    print("--------------END PIPELINE--------------")


if __name__ == "__main__":
    main()
