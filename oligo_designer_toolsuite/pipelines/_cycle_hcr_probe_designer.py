"""
CycleHCR probe designer pipeline.

CycleHCR, or cyclic Hybridization Chain Reaction, is a multiplexed RNA imaging
method that combines DNA barcoding with split-initiator HCR amplification.
Stable primary probes remain bound to the target RNA, while readout probes and
fluorescent HCR hairpins are exchanged across imaging cycles.

See :class:`CycleHCRProbeDesigner` for the full pipeline description and probe
structure. See :func:`cycle_hcr_probe_designer` for the config-driven workflow.
"""

############################################
# imports
############################################

import itertools
import os
import shutil
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import yaml
from Bio.SeqUtils import Seq
from pydantic import ValidationError

from oligo_designer_toolsuite._exceptions import (
    ConfigurationError,
    FeatureNotImplementedError,
    FileFormatError,
)
from oligo_designer_toolsuite.config import CycleHcrProbeDesignerConfig
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
    validate_primer_sequence,
)
from oligo_designer_toolsuite.sequence_generator import OligoSequenceGenerator
from oligo_designer_toolsuite.utils import configure_root_logger, logger

############################################
# CycleHCR Probe Designer
############################################


class CycleHCRProbeDesigner:
    """
    Design probe libraries for CycleHCR experiments.

    This class runs the full design workflow for CycleHCR primary probes,
    readout assignments, and DNA template probes. The output can be used to
    prepare probe libraries for multiplexed RNA detection with cyclic HCR
    readout.

    Overview
    --------
    CycleHCR, or cyclic Hybridization Chain Reaction, is an imaging method for
    detecting many RNA or protein targets in the same sample. It combines DNA
    barcoding with HCR signal amplification.

    In RNA CycleHCR, stable primary probes first bind to their target
    transcripts. These probes stay bound during repeated imaging cycles. In each
    cycle, short readout probes bind to barcode sequences on the primary probes.
    The readout probes carry split HCR initiators. When the correct left and
    right readouts bind next to each other, the initiator is brought together and
    starts HCR amplification with fluorescent hairpins.

    After imaging, the readout probes and HCR hairpins are stripped away. The
    primary probes remain on the RNA, so a new set of readout probes can be
    added in the next cycle.

    Probe Structure
    ---------------
    **Hybridization (Primary) Probes**

    CycleHCR uses pairs of primary probes for each target site. One probe binds
    to the left side of the target region and the other binds to the right side.
    The two target-binding regions are placed next to each other on the RNA,
    usually with a small gap between them.

    Each primary probe contains:

    - a gene-specific target-binding sequence,
    - a short barcode sequence used by the readout probes,
    - a linker sequence between the target-binding and barcode parts.

    The left and right probes are designed to bind strongly to the RNA target.
    This is important because the primary probes should remain bound while
    readout probes and HCR hairpins are removed between imaging cycles.

    A simplified layout is::

        Left probe:
            [Readout barcode L] + [linker] + [target-binding sequence L]

        Right probe:
            [target-binding sequence R] + [linker] + [Readout barcode R]

    **Readout Probes**

    Readout probes bind to the barcode sequences carried by the primary probes.
    For one target site, the left and right readout probes bind to matching
    barcode sequences on the primary probe pair.

    Each readout probe carries one half of an HCR initiator. When both readout
    probes bind at the same target site, the two initiator halves are brought
    together and trigger HCR amplification.

    The fluorescent signal comes from the HCR hairpins, not from the primary
    probes themselves.

    **Codebook**

    The codebook assigns each target to a barcode pattern. Rows correspond to
    target genes, and columns correspond to readout bits or imaging-cycle
    positions, depending on the design.

    A value of ``1`` means that the target carries the matching left and right
    readout barcode pair. During an imaging cycle, the corresponding readout
    probes bind to the primary probes and bring together the split HCR initiator.

    This barcode assignment determines in which cycle a target is detected.
    Because the signal is generated through HCR hairpins, the codebook controls
    readout identity, while the fluorescence comes from the amplifier used in
    that cycle.

    **DNA Template Probes**

    DNA template probes are the synthesis-ready sequences used to produce the
    final primary probe library. Each template contains primer binding sites for
    amplification, the target-binding sequence, the linker, and the readout
    barcode sequence.

    A simplified layout is::

        Left template:
            [forward primer] + [rc(target-binding sequence L)] + [linker]
            + [rc(readout barcode L)] + [reverse primer]

        Right template:
            [forward primer] + [rc(readout barcode R)] + [linker]
            + [rc(target-binding sequence R)] + [reverse primer]

    Probe Library Preparation
    -------------------------
    The designed DNA template library can be ordered as a pooled oligo library.
    The pool is PCR-amplified, transcribed, reverse-transcribed, and processed to
    produce single-stranded DNA primary probes.

    During the experiment, these primary probes are hybridized to the sample.
    Imaging is then performed over several cycles. In each cycle, selected
    readout probes and fluorescent HCR hairpins are added, imaged, and stripped
    before the next cycle begins.

    Pipeline Overview
    -----------------
    The pipeline performs the main design steps needed for a CycleHCR probe
    library:

    1. **Target probe design**

       Design left and right target-binding sequences for each RNA target. The
       probe halves are placed next to each other on the transcript and filtered
       for sequence properties and specificity.

    2. **Readout probe assignment**

       Load or generate readout sequences and assign barcode patterns to target
       regions through a codebook.

    3. **Hybridization probe assembly**

       Combine target-binding sequences, linker sequences, and assigned readout
       barcodes to build the left and right primary probes.

    4. **Primer handling**

       Load and validate the forward and reverse primer sequences used for
       template amplification.

    5. **DNA template assembly**

       Build the final DNA template probes from primers, target-binding
       sequences, linker sequences, and readout barcode sequences.

    6. **Output generation**

       Write the designed probes, codebook, readout assignments, and related
       probe information to the output directory.

    References
    ----------
    Gandin, V., Kim, J., Yang, L. Z., Lian, Y., Kawase, T., Hu, A., et al.
    Deep-tissue transcriptomics and subcellular imaging at high spatial
    resolution. Science. 2025. doi: 10.1126/science.adq2084

    :param dir_output: Directory where output files and intermediate results are
        saved. The directory is created if it does not exist.
    :type dir_output: str
    :param write_intermediate_steps: If ``True``, save intermediate probe
        databases after pipeline steps. This can help with checking or
        debugging a design run.
    :type write_intermediate_steps: bool
    :param n_jobs: Number of worker processes used for steps that can run in
        parallel. Use ``1`` to run without parallel processing.
    :type n_jobs: int
    """

    def __init__(
        self,
        dir_output: str,
        write_intermediate_steps: bool,
        n_jobs: int,
    ) -> None:
        """Constructor for the CycleHCRProbeDesigner class."""

        self.dir_output = os.path.abspath(dir_output)
        Path(self.dir_output).mkdir(parents=True, exist_ok=True)

        self.write_intermediate_steps = write_intermediate_steps
        self.n_jobs = n_jobs

    def design_target_probes(
        self,
        target_probes_parameters: dict,
    ) -> OligoDatabase:
        """
        Design the RNA-binding parts of the CycleHCR primary probes.

        This step designs left and right target-binding sequences for each target
        region. The two halves are chosen so they bind close to each other on the RNA
        transcript. The candidate probes are filtered for sequence quality and
        specificity, and a final probe set is selected for each target.

        The melting temperature is also calculated separately for the left and right
        halves. These values are stored with each probe and can be inspected in the
        output files.

        :param target_probes_parameters: Settings for target probe design from the
            ``target_probes`` section of the pipeline config. This includes probe
            generation, sequence filters, specificity filters, and probe set
            selection settings.
        :type target_probes_parameters: dict
        :return: Database containing the selected left and right target-binding
            probe halves for each target region.
        :rtype: OligoDatabase
        """
        target_probe_designer = TargetProbeDesigner(self.dir_output, self.n_jobs)
        target_probes_database = target_probe_designer.generate_target_probes(
            target_probes_parameters=target_probes_parameters,
            write_intermediate_steps=self.write_intermediate_steps,
        )

        # Per-arm Tm: primary probes must stay bound during stripping, so each half is
        # scored independently for the output tables (Tm params come from Tm_filter).
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

    def design_readout_probes(
        self,
        region_ids: list[str],
        readout_probe_parameters: dict,
    ) -> tuple[pd.DataFrame, pd.DataFrame]:
        """
        Load or create the CycleHCR readout probes and codebook.

        The codebook assigns each target region to a pair of readout probes. These
        readout probes later bind to barcode sequences on the primary probes and
        bring together the split HCR initiator during imaging.

        The codebook and readout probe table can either be loaded from files or
        generated from the config. The readout probe table is checked as soon as it
        is available. After the codebook is ready, both tables are checked
        together, so missing targets or missing readout sequences are caught
        early.

        :param region_ids: Target regions that must be represented in the codebook,
            usually gene names or gene IDs.
        :type region_ids: list[str]
        :param readout_probe_parameters: Settings from the ``readout_probes`` section
            of the pipeline config. This includes the codebook settings and the
            readout probe table settings.
        :type readout_probe_parameters: dict
        :return: Codebook and readout probe table used to assign readout probes to
            target regions.
        :rtype: tuple[pd.DataFrame, pd.DataFrame]
        """
        readout_probe_designer = ReadoutProbeDesigner(
            dir_output=self.dir_output,
            n_jobs=self.n_jobs,
        )

        if readout_probe_parameters["readout_probe_table"]["source"] == "load":
            readout_probe_table = readout_probe_designer.load_readout_probe_table(
                file_readout_probe_table=readout_probe_parameters["readout_probe_table"]["file"],
                codebook_source=readout_probe_parameters["codebook"]["source"],
            )
            readout_probe_table_source = readout_probe_parameters["readout_probe_table"]["file"]
            logger.info(
                f"Loaded readout probes table from file and retrieved {len(readout_probe_table)} readout probes."
            )
        else:
            readout_probe_table = readout_probe_designer.generate_readout_probe_table()
            readout_probe_table_source = readout_probe_parameters["readout_probe_table"]["source"]

        # Table-only validation: catch structural issues (index name, columns, DNA validity)
        # before generate_codebook consumes the table.
        readout_probe_designer.validate(
            readout_probe_table=readout_probe_table,
            readout_probe_table_source=readout_probe_table_source,
        )

        if readout_probe_parameters["codebook"]["source"] == "load":
            codebook = readout_probe_designer.load_codebook(
                file_codebook=readout_probe_parameters["codebook"]["file"]
            )
            codebook_source = readout_probe_parameters["codebook"]["file"]
        else:
            codebook = readout_probe_designer.generate_codebook(
                region_ids=region_ids,
                readout_probe_table=readout_probe_table,
                min_hamming_distance=readout_probe_parameters["codebook"]["min_hamming_distance"],
            )
            codebook_source = readout_probe_parameters["codebook"]["source"]

        # Full pair validation: also checks that every codebook bit is covered by the table.
        readout_probe_designer.validate(
            readout_probe_table=readout_probe_table,
            readout_probe_table_source=readout_probe_table_source,
            codebook=codebook,
            region_ids=region_ids,
            codebook_source=codebook_source,
        )

        return codebook, readout_probe_table

    def assemble_hybridization_probes(
        self,
        oligo_database: OligoDatabase,
        hybridization_probe_parameters: dict,
        codebook: pd.DataFrame,
        readout_probe_table: pd.DataFrame,
    ) -> OligoDatabase:
        """
        Build the CycleHCR primary probes.

        This step combines the RNA-binding probe halves with the readout barcode
        sequences assigned by the codebook. Each target receives a left and right
        primary probe. The target-binding part keeps the probe on the RNA, while the
        readout barcode is used during imaging cycles.

        A simplified layout is::

            Left primary probe:
                [readout barcode L] + [rc(linker)] + [target-binding sequence L]

            Right primary probe:
                [target-binding sequence R] + [rc(linker)] + [readout barcode R]

        The assembled sequences are added to the existing probe database.

        :param oligo_database: Database returned by :py:meth:`design_target_probes`.
            This database is updated with the assembled primary probe sequences.
        :type oligo_database: OligoDatabase
        :param hybridization_probe_parameters: Settings from the
            ``hybridization_probes`` section of the pipeline config, including the
            linker sequence.
        :type hybridization_probe_parameters: dict
        :param codebook: Table returned by :py:meth:`design_readout_probes`. Rows
            are target regions and columns are barcode bits.
        :type codebook: pd.DataFrame
        :param readout_probe_table: Table returned by
            :py:meth:`design_readout_probes`. Links each barcode bit to its readout
            probe sequence.
        :type readout_probe_table: pd.DataFrame
        :return: Database with left and right CycleHCR primary probe sequences added.
        :rtype: OligoDatabase
        """
        linker_sequence = hybridization_probe_parameters["linker_sequence"]

        region_ids = list(oligo_database.database.keys())

        oligo_database.set_database_sequence_types(
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
            # Weight-2 codebook: first active bit is L, second is R (validated upstream).
            sequence_readout_probe_L = readout_probe_sequences.iloc[0]
            sequence_readout_probe_R = readout_probe_sequences.iloc[1]

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
                new_properties[probe_id]["sequence_readout_probe_L"] = sequence_readout_probe_L
                new_properties[probe_id]["sequence_readout_probe_R"] = sequence_readout_probe_R

                # Linker is stored on the template strand; RC so it sits correctly on the probe.
                new_properties[probe_id]["sequence_hybridization_probe_L"] = (
                    sequence_readout_probe_L
                    + str(Seq(linker_sequence).reverse_complement())
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
                    + str(Seq(linker_sequence).reverse_complement())
                    + sequence_readout_probe_R
                )

            oligo_database.update_oligo_properties(new_properties)

        return oligo_database

    def design_primers(
        self,
        primer_parameters: dict,
    ) -> tuple[str, str]:
        """
        Load and check the PCR primers for DNA template probes.

        The forward and reverse primers are used to amplify the DNA template probe
        library before the final primary probes are prepared.

        The primers can either be loaded from files or generated from the config.
        Both primer sequences are checked to make sure they are valid DNA sequences.

        :param primer_parameters: Settings from the ``primers`` section of the
            pipeline config. This includes the forward and reverse primer entries
            and their sequences.
        :type primer_parameters: dict
        :return: Reverse primer sequence and forward primer sequence.
        :rtype: tuple[str, str]
        :raises FeatureNotImplementedError: If primer generation is requested.
        :raises FileFormatError: If a primer is empty or contains characters other
            than ``A``, ``C``, ``G``, and ``T``.
        """
        primer_designer = PrimerDesigner(
            dir_output=self.dir_output,
            n_jobs=self.n_jobs,
        )

        if primer_parameters["forward_primer"]["source"] == "load":
            forward_primer_sequence = primer_designer.load_forward_primer(
                sequence=primer_parameters["forward_primer"]["sequence"]
            )
        else:
            forward_primer_sequence = primer_designer.generate_forward_primer()

        if primer_parameters["reverse_primer"]["source"] == "load":
            reverse_primer_sequence = primer_designer.load_reverse_primer(
                sequence=primer_parameters["reverse_primer"]["sequence"]
            )
        else:
            reverse_primer_sequence = primer_designer.generate_reverse_primer()

        primer_designer.validate(
            forward_primer=forward_primer_sequence,
            reverse_primer=reverse_primer_sequence,
        )

        return reverse_primer_sequence, forward_primer_sequence

    def assemble_dna_template_probes(
        self,
        oligo_database: OligoDatabase,
        hybridization_probe_parameters: dict,
        forward_primer_sequence: str,
        reverse_primer_sequence: str,
    ) -> OligoDatabase:
        """
        Build the synthesis-ready DNA template probes.

        This step adds the forward and reverse primer sequences to each left and
        right probe template. The internal probe parts are written as reverse
        complements so that the downstream library preparation produces the intended
        CycleHCR primary probe strand.

        A simplified layout is::

            Left template:
                [forward primer] + [rc(target-binding sequence L)] + [linker]
                + [rc(readout barcode L)] + [reverse primer]

            Right template:
                [forward primer] + [rc(readout barcode R)] + [linker]
                + [rc(target-binding sequence R)] + [reverse primer]

        The resulting DNA template sequences are the sequences used for pooled oligo
        synthesis.

        :param oligo_database: Database returned by
            :py:meth:`assemble_hybridization_probes`. This database is updated with
            the DNA template sequences.
        :type oligo_database: OligoDatabase
        :param hybridization_probe_parameters: Settings from the
            ``hybridization_probes`` section of the pipeline config, including the
            linker sequence.
        :type hybridization_probe_parameters: dict
        :param forward_primer_sequence: Forward primer sequence returned by
            :py:meth:`design_primers`.
        :type forward_primer_sequence: str
        :param reverse_primer_sequence: Reverse primer sequence returned by
            :py:meth:`design_primers`.
        :type reverse_primer_sequence: str
        :return: Database with left and right DNA template probe sequences added.
        :rtype: OligoDatabase
        """
        linker_sequence = hybridization_probe_parameters["linker_sequence"]

        region_ids = list(oligo_database.database.keys())
        oligo_database.set_database_sequence_types(
            [
                "sequence_reverse_primer",
                "sequence_forward_primer",
                "sequence_dna_template_probe_L",
                "sequence_dna_template_probe_R",
            ]
        )

        for region_id in region_ids:
            probe_ids = list(oligo_database.database[region_id].keys())
            new_properties: dict[str, dict[str, str]] = {probe_id: {} for probe_id in probe_ids}

            for probe_id in probe_ids:
                new_properties[probe_id]["sequence_reverse_primer"] = reverse_primer_sequence
                new_properties[probe_id]["sequence_forward_primer"] = forward_primer_sequence

                # Internal parts are RC so reverse transcription yields the probe strand.
                new_properties[probe_id]["sequence_dna_template_probe_L"] = (
                    forward_primer_sequence
                    + str(
                        Seq(
                            format_sequence(
                                database=oligo_database,
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
                                database=oligo_database,
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
                                database=oligo_database,
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
                                database=oligo_database,
                                property="sequence_oligo_R",
                                region_id=region_id,
                                oligo_id=probe_id,
                            )
                        ).reverse_complement()
                    )
                    + reverse_primer_sequence
                )

            oligo_database.update_oligo_properties(new_properties)

        return oligo_database

    def generate_output(
        self,
        oligo_database: OligoDatabase,
        codebook: pd.DataFrame,
        readout_probe_table: pd.DataFrame,
        output_properties: list[str] | None = None,
    ) -> None:
        """
        Write the completed CycleHCR probe design to files.

        This step saves the final probe database, the codebook, and the readout
        probe table. It also writes an order-ready file that contains the DNA
        template probe sequences and readout probe sequences needed for synthesis.

        If no output properties are provided, a default set of probe annotations and
        sequence fields is written.

        :param oligo_database: Database returned by
            :py:meth:`assemble_dna_template_probes`.
        :type oligo_database: OligoDatabase
        :param codebook: Table assigning target regions to readout barcode bits.
            Rows are target regions and columns are barcode bits.
        :type codebook: pd.DataFrame
        :param readout_probe_table: Table linking each barcode bit to its readout
            probe sequence and related readout information.
        :type readout_probe_table: pd.DataFrame
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

        codebook.to_csv(os.path.join(self.dir_output, "codebook.tsv"), sep="\t", index_label="gene_name")
        readout_probe_table.to_csv(os.path.join(self.dir_output, "readout_probes.tsv"), sep="\t")

        oligo_database.write_oligosets_to_yaml(
            properties=output_properties,
            ascending=True,
            filename="cyclehcr_probes",
        )

        oligo_database.write_oligosets_to_table(
            properties=output_properties,
            ascending=True,
            filename="cyclehcr_probes",
        )

        oligo_database.write_ready_to_order_yaml(
            properties=[
                "sequence_dna_template_probe_L",
                "sequence_dna_template_probe_R",
                "sequence_readout_probe_L",
                "sequence_readout_probe_R",
            ],
            ascending=True,
            filename="cyclehcr_probes_order",
        )


############################################
# CycleHCR Target Probe Designer
############################################


class TargetProbeDesigner:
    """
    Design the RNA-binding probe pairs used by CycleHCR.

    This class designs the gene-specific part of the CycleHCR primary probes.
    Each probe candidate is split into a left and a right target-binding arm.
    The two arms are placed close to each other on the RNA target, with a short
    gap between them.

    Both arms must bind well and bind specifically. This is important because
    CycleHCR primary probes should stay on the RNA while readout probes and HCR
    hairpins are removed between imaging cycles. The split design also helps
    reduce background, because both arms need to bind at the right place to give
    a useful readout later.

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

       Select suitable probe sets for each target region, while keeping probes
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
        Run the full target-probe design workflow.

        This method designs the left and right RNA-binding arms used in CycleHCR
        primary probes. It starts from transcript sequences, creates candidate probe
        pairs, filters them, checks their specificity, and selects final probe sets
        for each target region.

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
            region.
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
            Tm_score=probe_set_selection_parameters["Tm_score"],
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
        one CycleHCR probe pair: the left binding site, the gap, and the right
        binding site.

        The candidate is then converted into the strand that will be used as the DNA
        probe and split into the left arm, spacer, and right arm. Regions with too
        few candidate probes are removed at this stage.

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
        hard_masked_sequences_filter: dict,
        soft_masked_sequences_filter: dict,
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
        the same panel. For CycleHCR, this is done separately for the left and right
        arms, because both arms need to behave well on their own.

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
        directories: list[str] = []

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
        Tm_score: dict,
    ) -> OligoDatabase:
        """
        Select final probe sets for each target region.

        This step chooses groups of probes from the filtered candidates. The selected
        probes should be well spaced along the transcript and should meet the
        requested number of probes per target region.

        Probe sets are scored using properties such as isoform coverage and melting
        temperature. The method can keep more than one possible probe set per region,
        which gives users alternatives when several good designs are available.
        Regions without enough suitable probes are removed.

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
        Tm_scorer = DeviationFromOptimalTmScorer(
            Tm_opt=Tm_score["Tm_opt"],
            Tm_parameters=Tm_score["Tm_parameters"],
            Tm_salt_correction_parameters=Tm_score["Tm_salt_correction_parameters"],
            Tm_chem_correction_parameters=Tm_score["Tm_chem_correction_parameters"],
            score_weight=Tm_score["weight"],
        )
        oligos_scoring = OligoScoring(scorers=[isoform_consensus_scorer, Tm_scorer])
        # ascending=False: higher aggregate scores win; Tm proximity to Tm_opt keeps
        # probes bound through cyclic stripping.
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
# CycleHCR Readout Probe Designer
############################################


class ReadoutProbeDesigner:
    """
    Manage the readout probes and codebook used for CycleHCR decoding.

    In CycleHCR, each target is assigned a two-bit barcode. One active bit
    selects the left readout probe, and one active bit selects the right readout
    probe. These two readout probes bind to the barcode sequences carried by the
    primary probe pair. Together, they bring the split HCR initiator into place
    during imaging.

    This class handles the two tables needed for that assignment:

    - the **readout probe table**, which links each bit to a readout probe
      sequence, an L/R role, a channel, and a readout probe ID,
    - the **codebook**, which assigns each target region to exactly two active
      bits.

    Both tables can be loaded from files. The codebook can also be generated
    from an existing readout probe table. Before the tables are used, they are
    checked together to make sure that all requested target regions and barcode
    bits are covered.

    :param dir_output: Directory where output files and intermediate results are
        saved.
    :type dir_output: str
    :param n_jobs: Number of worker processes used for steps that can run in
        parallel. This value is kept for consistency with the other designer
        classes.
    :type n_jobs: int
    """

    def __init__(self, dir_output: str, n_jobs: int) -> None:
        """Constructor for the ReadoutProbeDesigner class."""

        self.dir_output = os.path.abspath(dir_output)
        self.n_jobs = n_jobs

    def load_codebook(self, file_codebook: str) -> pd.DataFrame:
        """
        Load a CycleHCR codebook from a table file.

        The codebook file should contain one row per target region and one column
        per barcode bit. The target region column must be named ``gene_name``. Bit
        columns are expected to be named ``bit_1``, ``bit_2``, and so on.

        The table is loaded here, but the full content check is done later together
        with the readout probe table. This makes sure that the codebook bits and the
        available readout probes match.

        :param file_codebook: Path to the codebook file. The file must contain a
            ``gene_name`` column and one or more ``bit_*`` columns.
        :type file_codebook: str
        :return: Codebook table indexed by ``gene_name``.
        :rtype: pd.DataFrame
        """
        return pd.read_csv(file_codebook, sep=None, engine="python", index_col="gene_name")

    def generate_codebook(
        self, region_ids: list[str], readout_probe_table: pd.DataFrame, min_hamming_distance: int
    ) -> pd.DataFrame:
        """
        Generate a CycleHCR codebook from the readout probe table.

        The generated codebook assigns each target region to one left readout probe
        and one right readout probe. Each row is a target region. Each column is a
        barcode bit. A value of ``1`` means that the target uses the readout probe
        linked to that bit.

        CycleHCR uses two active bits per target. One bit must point to an ``L``
        readout probe, and one bit must point to an ``R`` readout probe. Both bits
        are assigned within the same fluorescence channel, so the two readout probes
        form a matched L/R pair for that channel.

        The available barcode space is taken from the readout probe table:

        - ``channel`` defines the fluorescence channel or readout channel,
        - ``L/R`` defines whether a readout probe is used on the left or right side,
        - ``readout_probe_id`` groups matching left and right readout probes.

        For each channel, the method builds barcode candidates from all possible
        combinations of left and right readout probe IDs. Combinations where the
        left and right readout probes have the same ID are tried first. These matched
        pairs are useful because they keep the barcode design simple and allow
        larger distances between codewords.

        The requested ``min_hamming_distance`` controls how strictly barcodes must
        differ from each other:

        - ``4`` keeps only matched pairs where ``L`` and ``R`` use the same
        readout probe ID. This gives disjoint two-bit barcodes and allows
        single-bit errors to be detected.
        - ``2`` allows all left/right combinations. This gives more barcodes, but
        barcodes may share one bit.
        - ``0`` is accepted for compatibility and uses the same full combination
        space as ``2``.

        Because every barcode has exactly two active bits, no other minimum Hamming
        distances are possible for distinct valid barcodes.

        :param region_ids: Target regions that need barcode assignments, usually
            gene names or gene IDs.
        :type region_ids: list[str]
        :param readout_probe_table: Table linking each bit to a readout probe
            sequence, channel, L/R role, and readout probe ID.
        :type readout_probe_table: pd.DataFrame
        :param min_hamming_distance: Required minimum distance between barcodes.
            Must be ``0``, ``2``, or ``4``.
        :type min_hamming_distance: int
        :return: Codebook table with target regions as rows and barcode bits as
            columns. Each row contains exactly two active bits.
        :rtype: pd.DataFrame
        :raises ConfigurationError: If the requested Hamming distance is not
            possible, or if there are not enough available barcodes for all target
            regions.
        """
        if min_hamming_distance not in (0, 2, 4):
            raise ConfigurationError(
                f"min_hamming_distance must be one of 0, 2, or 4 (got {min_hamming_distance}). "
                f"Each codeword has Hamming weight 2, so no other minimum distances are achievable."
            )

        n_channels = len(readout_probe_table["channel"].unique())
        n_readout_probes_R = readout_probe_table["L/R"].value_counts()["R"]
        n_readout_probes_L = readout_probe_table["L/R"].value_counts()["L"]
        n_readout_probes_LR = int(min([n_readout_probes_R, n_readout_probes_L]) / n_channels)

        def _generate_barcode(combination: tuple[int, int, int], codebook_size: int) -> list:
            """
            Convert one left/right/channel choice into a two-bit barcode.

            The ``combination`` contains three numbers: the left readout probe index,
            the right readout probe index, and the channel index. The method turns this
            choice into a barcode vector with two ``1`` values.

            Bit positions are arranged in repeated channel blocks. Within each channel
            block, the left bit comes first and the right bit comes second. This layout
            keeps the generated codebook aligned with the generated readout probe table.

            :param combination: Tuple containing the left probe index, right probe index,
                and channel index.
            :type combination: tuple[int, int, int]
            :param codebook_size: Total number of bit columns in the codebook.
            :type codebook_size: int
            :return: Barcode vector with exactly two active bits.
            :rtype: list
            """
            # Bit layout matches the generated readout table: per probe-id block,
            # channels alternate L then R (even = L, odd = R).
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
        # Prefer matched L/R IDs first; they give disjoint bit pairs (distance 4).
        combinations = sorted(combinations, key=lambda t: (0 if t[0] == t[1] else 1, t[1]))

        if min_hamming_distance == 4:
            # Only L == R combinations are pairwise distance 4.
            combinations = [c for c in combinations if c[0] == c[1]]
            if len(combinations) < n_regions:
                raise ConfigurationError(
                    f"Only {len(combinations)} barcodes are available at min_hamming_distance=4 "
                    f"(= n_readout_probes_LR * n_channels = {n_readout_probes_LR} * {n_channels}), "
                    f"which is fewer than the {n_regions} requested regions. "
                    f"Consider increasing n_readout_probes_LR or n_channels, reducing the number of "
                    f"regions, or lowering min_hamming_distance to 2."
                )
        else:
            # Distance 0 or 2: any L/R pair is allowed; only the total count matters.
            if len(combinations) < n_regions:
                raise ConfigurationError(
                    f"Only {len(combinations)} barcodes are available "
                    f"(= n_readout_probes_LR**2 * n_channels = {n_readout_probes_LR**2} * {n_channels}), "
                    f"which is fewer than the {n_regions} requested regions. "
                    f"Consider increasing n_readout_probes_LR or n_channels, or reducing the number "
                    f"of regions."
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
        codebook.index.name = "gene_name"

        return codebook

    def load_readout_probe_table(self, file_readout_probe_table: str, codebook_source: str) -> pd.DataFrame:
        """
        Load a CycleHCR readout probe table from a table file.

        The readout probe table links each barcode bit to a readout probe sequence.
        It must contain the readout probe sequence, the fluorescence channel, the
        L/R role, and the readout probe ID.

        Bit handling depends on how the codebook is provided.

        If the codebook is loaded from file, the readout probe table must already
        contain a ``bit`` column. These bit names must match the bit columns in the
        codebook. In this case, the user controls the bit-to-probe mapping.

        If the codebook will be generated, any existing ``bit`` column is ignored.
        The table is sorted by ``readout_probe_id``, ``channel``, and ``L/R``. New
        bit names are then assigned in this order. This creates the bit layout that
        :py:meth:`generate_codebook` expects, with left and right readout probes
        arranged in a predictable order.

        Column and sequence checks are done later by :py:meth:`validate`.

        :param file_readout_probe_table: Path to the readout probe table file. The
            file must contain ``readout_probe_sequence``, ``channel``,
            ``readout_probe_id``, and ``L/R`` columns. If the codebook is loaded
            from file, it must also contain a ``bit`` column.
        :type file_readout_probe_table: str
        :param codebook_source: Source of the codebook. Use ``"load"`` when the
            codebook is read from file. Other values indicate that the codebook will
            be generated.
        :type codebook_source: str
        :return: Readout probe table indexed by ``bit``.
        :rtype: pd.DataFrame
        :raises FileFormatError: If the ``bit`` column is missing when a codebook
            is loaded from file.
        """
        readout_probe_table = pd.read_csv(file_readout_probe_table, sep=None, engine="python")
        if codebook_source == "load":
            if "bit" not in readout_probe_table.columns:
                raise FileFormatError(
                    f"Readout probe table '{file_readout_probe_table}' must contain a 'bit' column "
                    f"when loading a codebook from file. The 'bit' values must match the codebook's "
                    f"bit columns; the user is responsible for that mapping."
                )
        else:
            # Regenerating the codebook: ignore any existing bit column and assign bits
            # in the order generate_codebook expects (probe id, channel, L then R).
            if "bit" in readout_probe_table.columns:
                readout_probe_table = readout_probe_table.drop(columns=["bit"])
            readout_probe_table = readout_probe_table.sort_values(by=["readout_probe_id", "channel", "L/R"])
            readout_probe_table.reset_index(inplace=True, drop=True)
            readout_probe_table["bit"] = "bit_" + (readout_probe_table.index + 1).astype(str)
        return readout_probe_table.set_index("bit")

    def generate_readout_probe_table(self) -> pd.DataFrame:
        """
        Generate a CycleHCR readout probe table.

        This feature is not implemented yet. For now, the readout probe table must
        be supplied as an input file.

        :return: Readout probe table indexed by ``bit``.
        :rtype: pd.DataFrame
        :raises FeatureNotImplementedError: Always raised until readout probe table
            generation is implemented.
        """
        raise FeatureNotImplementedError(
            "Generation of readout probe table is not yet implemented. "
            "Please provide a file_readout_probe_table parameter."
        )

    def validate(
        self,
        readout_probe_table: pd.DataFrame | None = None,
        *,
        readout_probe_table_source: str | None = None,
        codebook: pd.DataFrame | None = None,
        region_ids: list[str] | None = None,
        codebook_source: str | None = None,
    ) -> None:
        """
        Check the codebook and/or readout probe table.

        This method can check either table on its own, or both together. When a
        codebook is provided, it must contain all requested target regions, use
        ``gene_name`` as the row index, and contain only ``0`` and ``1`` values in
        its bit columns. For CycleHCR, each target region must have exactly two
        active bits.

        When a readout probe table is provided, it must include a valid DNA
        sequence, channel, L/R role, and readout probe ID for each bit. If a
        codebook is also provided, every bit used by the codebook must appear in
        the table.

        Running these checks before probe assembly helps catch mismatched files
        early, such as a codebook that refers to a bit missing from the readout
        probe table.

        :param readout_probe_table: Optional table linking each barcode bit to a
            readout probe sequence and its readout information. If ``None``, the
            table is not checked.
        :type readout_probe_table: pd.DataFrame | None
        :param readout_probe_table_source: File path or source label for the
            readout probe table. Used in error messages when the table is checked.
        :type readout_probe_table_source: str | None
        :param codebook: Optional codebook table assigning target regions to
            barcode bits. Rows are target regions and columns are ``bit_*``
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
            codebook, or ``readout_probe_table_source`` with a readout probe table).
        :raises FileFormatError: If the codebook or readout probe table is missing
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
                expected_hamming_weight=2,
                index_name="gene_name",
            )
        if readout_probe_table is not None:
            if readout_probe_table_source is None:
                raise ValueError(
                    "readout_probe_table_source must be provided when validating a readout probe table."
                )
            validate_bit_mapping_table(
                table=readout_probe_table,
                codebook=codebook,
                source=readout_probe_table_source,
                required_columns=["channel", "readout_probe_id", "L/R", "readout_probe_sequence"],
                sequence_columns=["readout_probe_sequence"],
            )


############################################
# CycleHCR Primer Designer
############################################


class PrimerDesigner:
    """
    Load and check the PCR primers used for CycleHCR DNA templates.

    CycleHCR primary probes are prepared from DNA template sequences. These
    templates contain the probe body flanked by a forward primer and a reverse
    primer. The primer pair is used to amplify the pooled template library before
    the final single-stranded probes are prepared.

    This class currently expects primer sequences to be provided in the pipeline
    config. The primer generation methods are placeholders for a future automatic
    primer-design step.

    :param dir_output: Directory where output files and intermediate results are
        saved.
    :type dir_output: str
    :param n_jobs: Number of worker processes used for steps that can run in
        parallel. This value is kept for consistency with the other designer
        classes.
    :type n_jobs: int
    """

    def __init__(self, dir_output: str, n_jobs: int) -> None:
        """Constructor for the PrimerDesigner class."""

        self.dir_output = os.path.abspath(dir_output)
        self.n_jobs = n_jobs

    def load_reverse_primer(self, sequence: str) -> str:
        """
        Load the reverse primer sequence from the config.

        The reverse primer is used together with the forward primer to amplify the
        DNA template probe pool. This method trims surrounding whitespace. The
        sequence itself is checked later by :py:meth:`validate`.

        :param sequence: Reverse primer sequence used for PCR amplification of the
            DNA template probes.
        :type sequence: str
        :return: Reverse primer sequence ready for validation.
        :rtype: str
        """
        reverse_primer = str(sequence).strip()
        return reverse_primer

    def generate_reverse_primer(self) -> str:
        """
        Generate a reverse primer sequence.

        Automatic reverse primer design is not available yet. Provide the reverse
        primer sequence in the pipeline config instead.

        :return: Reverse primer sequence.
        :rtype: str
        :raises FeatureNotImplementedError: Always raised until reverse primer
            generation is implemented.
        """
        raise FeatureNotImplementedError(
            "Generation of reverse primer is not yet implemented. "
            "Please provide a reverse_primer.sequence parameter and set reverse_primer.source to 'load'."
        )

    def load_forward_primer(self, sequence: str) -> str:
        """
        Load the forward primer sequence from the config.

        The forward primer is used together with the reverse primer to amplify the
        DNA template probe pool. This method trims surrounding whitespace. The
        sequence itself is checked later by :py:meth:`validate`.

        :param sequence: Forward primer sequence used for PCR amplification of the
            DNA template probes.
        :type sequence: str
        :return: Forward primer sequence ready for validation.
        :rtype: str
        """
        forward_primer = str(sequence).strip()
        return forward_primer

    def generate_forward_primer(self) -> str:
        """
        Generate a forward primer sequence.

        Automatic forward primer design is not available yet. Provide the forward
        primer sequence in the pipeline config instead.

        :return: Forward primer sequence.
        :rtype: str
        :raises FeatureNotImplementedError: Always raised until forward primer
            generation is implemented.
        """
        raise FeatureNotImplementedError(
            "Generation of forward primer is not yet implemented. "
            "Please provide a forward_primer.sequence parameter and set forward_primer.source to 'load'."
        )

    def validate(self, forward_primer: str, reverse_primer: str) -> None:
        """
        Check that both primer sequences are valid DNA sequences.

        Each primer must be a non-empty sequence containing only ``A``, ``C``,
        ``G``, and ``T``. This check catches missing primers and accidental
        characters before the template probes are assembled.

        :param forward_primer: Forward primer sequence to check.
        :type forward_primer: str
        :param reverse_primer: Reverse primer sequence to check.
        :type reverse_primer: str
        :raises FileFormatError: If either primer is empty or contains characters
            other than ``A``, ``C``, ``G``, and ``T``.
        """
        validate_primer_sequence(forward_primer, source="forward_primer")
        validate_primer_sequence(reverse_primer, source="reverse_primer")


############################################
# CycleHCR Probe Designer Pipeline
############################################


def _preprocess_config(config_validated: CycleHcrProbeDesignerConfig) -> dict[str, Any]:
    """
    Prepare the CycleHCR config before the pipeline runs.

    This step updates the config in place so later design stages can read ready-to-use
    settings. It resolves melting-temperature tables, turns off unused temperature
    corrections, and copies the shared temperature settings into the filters and
    scoring steps that need them.

    It also derives the full probe length and the left/right junction position from
    the arm and gap lengths, and expands an optional gene-list file into a concrete
    list of target regions. If no gene list is provided, all regions in the input
    FASTA files are used.

    :param config: Pipeline configuration loaded from the YAML config file.
    :type config: dict
    :return: The same config dict, updated with the prepared settings.
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


def cycle_hcr_probe_designer(config: CycleHcrProbeDesignerConfig) -> None:
    """
    Run the CycleHCR probe design pipeline from a config dict.

    This function prepares the config with :func:`_preprocess_config`, then runs
    :class:`CycleHCRProbeDesigner` end to end. It designs target probes, loads or
    creates the readout probes and codebook, assembles primary probes and DNA
    templates, and writes the final files under ``config['general']['dir_output']``.
    The caller should configure the library logger before calling this function
    (see :func:`main`).

    The config should follow ``data/configs/cycle_hcr_probe_designer.yaml``.

    Top-level config sections:

    - ``general``: output directory, intermediate-step writing, and worker count.
    - ``target_probes``: candidate generation, sequence filters, specificity filters,
      and probe set selection.
    - ``readout_probes``: codebook and readout probe table settings.
    - ``primers``: forward and reverse primer settings.
    - ``hybridization_probes``: linker sequence used during probe assembly.

    Files written under ``dir_output``:

    - ``codebook.tsv``: barcode assignments for each target gene.
    - ``readout_probes.tsv``: readout probe sequences and related bit information.
    - ``cyclehcr_probes.yml``: full probe records.
    - ``cyclehcr_probes_order.yml``: sequences ready for synthesis.
    - ``cyclehcr_probes.tsv`` / ``cyclehcr_probes.xlsx``: probe sets as tables.

    Intermediate probe databases are also written when
    ``general.write_intermediate_steps`` is ``True``.

    See :class:`CycleHCRProbeDesigner` for the pipeline description and probe
    structure.

    :param config: Validated pipeline configuration. It is converted and prepared by
        :func:`_preprocess_config` before the pipeline runs.
    :type config: CycleHcrProbeDesignerConfig
    :return: None
    :rtype: None
    """

    config_dict = _preprocess_config(config)

    pipeline = CycleHCRProbeDesigner(
        dir_output=config_dict["general"]["dir_output"],
        write_intermediate_steps=config_dict["general"]["write_intermediate_steps"],
        n_jobs=config_dict["general"]["n_jobs"],
    )

    target_probe_database = pipeline.design_target_probes(
        target_probes_parameters=config_dict["target_probes"],
    )

    codebook, readout_probe_table = pipeline.design_readout_probes(
        region_ids=list(target_probe_database.database.keys()),
        readout_probe_parameters=config_dict["readout_probes"],
    )

    hybridization_probe_database = pipeline.assemble_hybridization_probes(
        oligo_database=target_probe_database,
        hybridization_probe_parameters=config_dict["hybridization_probes"],
        codebook=codebook,
        readout_probe_table=readout_probe_table,
    )

    reverse_primer_sequence, forward_primer_sequence = pipeline.design_primers(
        primer_parameters=config_dict["primers"],
    )

    dna_template_probe_database = pipeline.assemble_dna_template_probes(
        oligo_database=hybridization_probe_database,
        hybridization_probe_parameters=config_dict["hybridization_probes"],
        forward_primer_sequence=forward_primer_sequence,
        reverse_primer_sequence=reverse_primer_sequence,
    )

    pipeline.generate_output(
        oligo_database=dna_template_probe_database,
        codebook=codebook,
        readout_probe_table=readout_probe_table,
    )


def main() -> None:
    """
    Run the CycleHCR probe design pipeline from the command line.

    Parses the required ``-c``/``--config`` argument, loads the YAML configuration
    file, and configures the library logger to write under the configured output
    directory. It then calls :func:`cycle_hcr_probe_designer`.

    :return: None
    :rtype: None
    """
    print("--------------START PIPELINE--------------")

    args = base_parser(
        prog="CycleHCR Probe Designer",
        usage="cycle_hcr_probe_designer [options]",
        description=__doc__,
    )

    with open(args["config"], "r") as handle:
        config_raw = yaml.safe_load(handle)

    try:
        config_validated = CycleHcrProbeDesignerConfig.model_validate(config_raw)
    except ValidationError as e:
        print(f"Invalid configuration file:\n{e}")
        raise

    configure_root_logger(
        dir_output=config_validated.general.dir_output,
        pipeline_name="cyclehcr_probe_designer",
    )

    cycle_hcr_probe_designer(config_validated)

    print("--------------END PIPELINE--------------")


if __name__ == "__main__":
    main()
