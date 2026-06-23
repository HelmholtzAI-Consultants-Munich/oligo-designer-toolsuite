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

from oligo_designer_toolsuite._exceptions import (
    ConfigurationError,
    FeatureNotImplementedError,
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

        self.write_intermediate_steps = write_intermediate_steps
        self.n_jobs = n_jobs

    def design_target_probes(
        self,
        oligo_generation_parameters: dict,
        property_filters_parameters: dict,
        specificity_filters_parameters: dict,
        probe_set_selection_parameters: dict,
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

        :param oligo_generation_parameters: ``target_probe.oligo_generation`` block. Contains
            ``region_ids`` (populated from ``file_region_ids``), ``files_fasta_probe_database``,
            ``L_probe_sequence_length``, ``gap_sequence_length``, ``R_probe_sequence_length``, and the
            derived ``oligo_length`` injected by :func:`_preprocess_config`.
        :type oligo_generation_parameters: dict
        :param property_filters_parameters: ``target_probe.property_filters`` block. Each filter
            sub-dict carries an ``enabled`` flag plus its parameters; ``Tm_filter`` additionally
            receives the inlined ``Tm_parameters`` / ``Tm_chem_correction_parameters`` /
            ``Tm_salt_correction_parameters`` from :func:`_preprocess_config`.
        :type property_filters_parameters: dict
        :param specificity_filters_parameters: ``target_probe.specificity_filters`` block.
            ``specificity_blastn_filter`` carries ``files_fasta_reference_database``,
            ``junction_region_size``, and the derived ``junction_site``;
            ``cross_hybridization_blastn_filter`` carries its own search/hit parameters.
        :type specificity_filters_parameters: dict
        :param probe_set_selection_parameters: ``target_probe.probe_set_selection`` block. Contains
            the ``independent_set_selection`` scalars and the ``isoform_consensus_score`` / ``Tm_score``
            sub-dicts (Tm parameters inlined into ``Tm_score`` by :func:`_preprocess_config`).
        :type probe_set_selection_parameters: dict
        :return: An `OligoDatabase` object containing the designed target probes organized into sets.
            The database includes probe sequences, properties, and set assignments for each target region.
        :rtype: OligoDatabase
        """
        target_probe_designer = TargetProbeDesigner(self.dir_output, self.n_jobs)

        oligo_database: OligoDatabase = target_probe_designer.create_oligo_database(
            region_ids=oligo_generation_parameters["region_ids"],
            oligo_length=oligo_generation_parameters["oligo_length"],
            L_probe_sequence_length=oligo_generation_parameters["L_probe_sequence_length"],
            gap_sequence_length=oligo_generation_parameters["gap_sequence_length"],
            R_probe_sequence_length=oligo_generation_parameters["R_probe_sequence_length"],
            files_fasta_oligo_database=oligo_generation_parameters["files_fasta_probe_database"],
            min_oligos_per_gene=probe_set_selection_parameters["independent_set_selection"]["set_size_min"],
        )

        if self.write_intermediate_steps:
            dir_database = oligo_database.save_database(name_database="1_db_target_probes_initial")
            logger.info(
                f"Saved target probe database for step 1 (Create Database) in directory {dir_database}"
            )

        oligo_database = target_probe_designer.filter_by_property(
            oligo_database=oligo_database,
            isoform_consensus_filter=property_filters_parameters["isoform_consensus_filter"],
            hard_masked_sequences_filter=property_filters_parameters["hard_masked_sequences_filter"],
            soft_masked_sequences_filter=property_filters_parameters["soft_masked_sequences_filter"],
            homopolymeric_runs_filter=property_filters_parameters["homopolymeric_runs_filter"],
            GC_content_filter=property_filters_parameters["GC_content_filter"],
            Tm_filter=property_filters_parameters["Tm_filter"],
            secondary_structure_filter=property_filters_parameters["secondary_structure_filter"],
        )

        if self.write_intermediate_steps:
            dir_database = oligo_database.save_database(name_database="2_db_target_probes_property_filter")
            logger.info(
                f"Saved target probe database for step 2 (Property Filters) in directory {dir_database}"
            )

        oligo_database = target_probe_designer.filter_by_specificity(
            oligo_database=oligo_database,
            specificity_blastn_filter=specificity_filters_parameters["specificity_blastn_filter"],
            cross_hybridization_blastn_filter=specificity_filters_parameters[
                "cross_hybridization_blastn_filter"
            ],
        )

        if self.write_intermediate_steps:
            dir_database = oligo_database.save_database(name_database="3_db_target_probes_specificity_filter")
            logger.info(
                f"Saved target probe database for step 3 (Specificity Filters) in directory {dir_database}"
            )

        oligo_database = target_probe_designer.create_oligo_sets(
            oligo_database=oligo_database,
            independent_set_selection=probe_set_selection_parameters["independent_set_selection"],
            isoform_consensus_score=probe_set_selection_parameters["isoform_consensus_score"],
            Tm_score=probe_set_selection_parameters["Tm_score"],
        )

        # Calculate Tm for both arms (always on). Tm parameters were inlined into the ``Tm_filter`` block by ``_preprocess_config``.
        Tm_filter = property_filters_parameters["Tm_filter"]
        tm_nn_property: BaseProperty = TmNNProperty(
            Tm_parameters=Tm_filter["Tm_parameters"],
            Tm_chem_correction_parameters=Tm_filter["Tm_chem_correction_parameters"],
            Tm_salt_correction_parameters=Tm_filter["Tm_salt_correction_parameters"],
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
            logger.info(
                f"Saved target probe database for step 4 (Specificity Filters) in directory {dir_database}."
            )

        return oligo_database

    def design_readout_probes(
        self,
        region_ids: list[str],
        readout_probe_parameters: dict,
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
        :param readout_probe_parameters: ``readout_probes`` block. Must contain
            ``file_readout_probe_table`` (required — generation of readout probe tables is not yet
            implemented) and ``file_codebook`` (optional — if not provided, a codebook is auto-generated).
        :type readout_probe_parameters: dict
        :return: A tuple containing (codebook, readout_probe_table), where:
            - codebook: A pandas DataFrame with region IDs as index and bit columns, where each
              row represents the barcode assignment for a region
            - readout_probe_table: A pandas DataFrame containing readout probe information with
              bit identifiers as index
        :rtype: tuple[pd.DataFrame, pd.DataFrame]
        """
        file_readout_probe_table = readout_probe_parameters["file_readout_probe_table"]
        file_codebook = readout_probe_parameters["file_codebook"]

        readout_probe_designer = ReadoutProbeDesigner(
            dir_output=self.dir_output,
            n_jobs=self.n_jobs,
        )

        if file_readout_probe_table:
            readout_probe_table = readout_probe_designer.load_readout_probe_table(
                file_readout_probe_table=file_readout_probe_table
            )
            readout_probe_table_source = file_readout_probe_table
            logger.info(
                f"Loaded readout probes table from file and retrieved {len(readout_probe_table)} readout probes."
            )
        else:
            readout_probe_table = readout_probe_designer.generate_readout_probe_table()
            readout_probe_table_source = "<generated>"

        if file_codebook:
            codebook = readout_probe_designer.load_codebook(file_codebook=file_codebook)
            codebook_source = file_codebook
        else:
            codebook = readout_probe_designer.generate_codebook(
                region_ids=region_ids,
                readout_probe_table=readout_probe_table,
            )
            codebook_source = "<generated>"

        readout_probe_designer.validate(
            codebook=codebook,
            readout_probe_table=readout_probe_table,
            region_ids=region_ids,
            codebook_source=codebook_source,
            readout_probe_table_source=readout_probe_table_source,
        )

        return codebook, readout_probe_table

    def assemble_hybridization_probes(
        self,
        oligo_database: OligoDatabase,
        hybridization_probe_parameters: dict,
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

        :param oligo_database: Database of target probes containing sequence and property information.
            This database should contain the designed target probes with their L and R oligo sequences
            organized by region and probe ID.
        :type oligo_database: OligoDatabase
        :param hybridization_probe_parameters: ``hybridization_probes`` block (with codebook and
            readout_probe_table injected by the orchestrator at runtime). Must contain
            ``linker_sequence``, ``codebook``, and ``readout_probe_table``.
        :type hybridization_probe_parameters: dict
        :return: An updated `OligoDatabase` object containing the assembled hybridization probes with all
            component sequences (target, oligo L/R, readout probe L/R, and complete hybridization probe L/R)
            stored as properties for each probe.
        :rtype: OligoDatabase
        """
        linker_sequence = hybridization_probe_parameters["linker_sequence"]
        codebook = hybridization_probe_parameters["codebook"]
        readout_probe_table = hybridization_probe_parameters["readout_probe_table"]

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
        Load and validate forward and reverse primer sequences for DNA template probe assembly.

        This method processes primer sequences that will be used for PCR amplification of the DNA
        template probes. The primers are incorporated into the DNA template probe structure during
        the assembly step, with the forward primer at the 5' end and the reverse primer at the 3' end.

        Currently, primer generation is not implemented, so both primer sequences must be provided
        in the configuration.

        :param primer_parameters: ``primers`` block. Must contain ``forward_primer_sequence`` and
            ``reverse_primer_sequence`` (both required — primer generation is not yet implemented).
        :type primer_parameters: dict
        :return: A tuple containing (reverse_primer_sequence, forward_primer_sequence) in that order.
            Both sequences have been stripped of leading and trailing whitespace.
        :rtype: tuple[str, str]
        :raises FeatureNotImplementedError: If either primer sequence is empty or None, since primer
            generation is not yet implemented.
        """
        forward_primer_sequence = primer_parameters["forward_primer_sequence"]
        reverse_primer_sequence = primer_parameters["reverse_primer_sequence"]

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
        oligo_database: OligoDatabase,
        hybridization_probe_parameters: dict,
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

        :param oligo_database: Database of hybridization probes containing sequence and property
            information. This database should contain the assembled hybridization probes with their
            component sequences (oligo L/R and readout probe L/R sequences).
        :type oligo_database: OligoDatabase
        :param hybridization_probe_parameters: ``hybridization_probes`` block. Must contain
            ``linker_sequence`` plus the orchestrator-injected ``forward_primer_sequence`` and
            ``reverse_primer_sequence`` (both produced by :py:meth:`design_primers`).
        :type hybridization_probe_parameters: dict
        :return: An updated `OligoDatabase` object containing the assembled DNA template probes with
            all sequences stored as properties, including: sequence_forward_primer, sequence_reverse_primer,
            sequence_dna_template_probe_L, and sequence_dna_template_probe_R for each probe.
        :rtype: OligoDatabase
        """
        linker_sequence = hybridization_probe_parameters["linker_sequence"]
        forward_primer_sequence = hybridization_probe_parameters["forward_primer_sequence"]
        reverse_primer_sequence = hybridization_probe_parameters["reverse_primer_sequence"]

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
        Generate the final output files for the CycleHCR probe design pipeline.

        This method creates all output files needed for the CycleHCR experiment, including:
        1. Codebook and readout probe tables (TSV format)
        2. Region-specific readout probe mapping (TSV format)
        3. Complete probe sets in YAML format with all properties
        4. Ready-to-order probe sequences in YAML format
        5. Probe sets in tabular format (TSV and Excel with one sheet per region)

        All output files are written to the output directory specified during pipeline initialization.

        :param oligo_database: Database of DNA template probes with associated properties and sequences.
            This database should contain the final assembled probes with all component sequences
            (target, oligo L/R, readout probe L/R, hybridization probe L/R, DNA template probe L/R).
        :type oligo_database: OligoDatabase
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
        codebook.to_csv(os.path.join(self.dir_output, "codebook.tsv"), sep="\t", index_label="gene_name")
        readout_probe_table.to_csv(os.path.join(self.dir_output, "readout_probes.tsv"), sep="\t")

        readout_probe_table_regions = []
        for gene_name, barcode in codebook.iterrows():
            bits = barcode[barcode == 1].index
            readout_probe_info = readout_probe_table.loc[bits, :]
            readout_probe_info["gene_name"] = gene_name
            readout_probe_table_regions.append(readout_probe_info)
        readout_probe_table_regions_df = pd.concat(readout_probe_table_regions, axis=0)
        readout_probe_table_regions_df[
            ["gene_name", "channel", "readout_probe_id", "L/R", "readout_probe_sequence"]
        ].to_csv(os.path.join(self.dir_output, "readout_probes_regions.tsv"), sep="\t", index=False)

        oligo_database.write_oligosets_to_yaml(
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

        oligo_database.write_oligosets_to_table(
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
        oligo_length: int,
        L_probe_sequence_length: int,
        gap_sequence_length: int,
        R_probe_sequence_length: int,
        files_fasta_oligo_database: list[str],
        min_oligos_per_gene: int,
    ) -> OligoDatabase:
        """
        Create an initial oligo database by generating target probe sequences.

        This method performs the first step of target probe design by:
        1. Generating candidate oligo sequences from input FASTA files using a sliding window approach
           of length ``oligo_length`` (= L + gap + R, precomputed by :func:`_preprocess_config`)
        2. Creating an oligo database with the generated sequences
        3. Calculating the reverse complement sequence (``oligo`` sequence type) — always on
        4. Splitting each oligo sequence into the left probe, gap spacer, and right probe components

        Each oligo is split into three components: the right probe (5' end), a spacer (gap region),
        and the left probe (3' end). This split is performed on the reverse complement of the target
        sequence to generate the actual probe sequences that will hybridize to the RNA.

        Isoform-consensus property computation and filtering live in :py:meth:`filter_by_property`
        so all ``filter_database_by_property_threshold`` calls remain in the property-filter phase.

        :param region_ids: List of region identifiers (e.g., gene IDs) for which oligos should be
            generated. If None, all regions present in the input FASTA files will be processed.
        :type region_ids: list[str] | None
        :param oligo_length: Total oligo length (= L + gap + R), precomputed by :func:`_preprocess_config`.
        :type oligo_length: int
        :param L_probe_sequence_length: Length of the left probe sequence in nucleotides.
        :type L_probe_sequence_length: int
        :param gap_sequence_length: Length of the gap sequence between left and right probes
            in nucleotides. This gap is not included in the probe sequences but represents the
            spacing between the two probe halves on the target transcript.
        :type gap_sequence_length: int
        :param R_probe_sequence_length: Length of the right probe sequence in nucleotides.
        :type R_probe_sequence_length: int
        :param files_fasta_oligo_database: List of paths to FASTA files containing genomic sequences
            from which target probes will be generated.
        :type files_fasta_oligo_database: list[str]
        :param min_oligos_per_gene: Minimum number of oligos required per region (gene) after
            generation. Regions with fewer oligos than this threshold will be removed from the database.
        :type min_oligos_per_gene: int
        :return: An `OligoDatabase` object containing the generated target probe sequences with their
            component sequences (target, oligo, oligo_L, oligo_R, spacer). The database is filtered
            to only include regions that meet the minimum oligo requirement.
        :rtype: OligoDatabase
        """
        ##### creating the oligo sequences #####
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

        ##### compute reverse complement (always on) #####
        reverse_complement_sequence_property: BaseProperty = ReverseComplementSequenceProperty(
            sequence_type_reverse_complement="oligo"
        )
        calculator = PropertyCalculator(properties=[reverse_complement_sequence_property])
        oligo_database = calculator.apply(
            oligo_database=oligo_database, sequence_type="target", n_jobs=self.n_jobs
        )

        ##### split the oligo sequence into L arm + spacer + R arm (always on) #####
        split_start_end = [
            (0, L_probe_sequence_length),
            (L_probe_sequence_length, L_probe_sequence_length + gap_sequence_length),
            (
                L_probe_sequence_length + gap_sequence_length,
                L_probe_sequence_length + gap_sequence_length + R_probe_sequence_length,
            ),
        ]
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

        oligo_database.remove_regions_with_insufficient_oligos(pipeline_step="Database Creation")
        check_content_oligo_database(oligo_database)

        return oligo_database

    @pipeline_step_basic(step_name="Target Probe Generation - Property Filters")
    def filter_by_property(
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
        Filter the oligo database based on various sequence properties.

        This method applies multiple property-based filters, each gated on its own ``enabled``
        flag. The isoform-consensus filter is applied first as a cheap pre-filter on the
        ``target`` sequence type; the remaining sequence-based filters are then applied to both
        the left (L) and right (R) probe sequences. Probes that fail any enabled filter are removed
        from the database.

        The following filters are applied (in this order, when enabled):
        1. **Isoform consensus** (cheap pre-filter): computes ``IsoformConsensusProperty`` on the
           ``target`` sequence and removes regions below the configured threshold.
        2. **Hard masked sequences**: Removes probes containing hard-masked nucleotides (N)
        3. **Soft masked sequences**: Removes probes containing soft-masked nucleotides (lowercase)
        4. **Homopolymeric runs**: Removes probes with homopolymeric runs exceeding specified lengths
        5. **GC content**: Removes probes with GC content outside the specified range
        6. **Melting temperature**: Removes probes with calculated Tm outside the specified range
        7. **Secondary structure**: Removes probes that form stable secondary structures at the
           specified temperature

        Regions that do not meet the minimum oligo requirement after filtering are removed from
        the database.

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
        :param Tm_filter: Dict with ``enabled``, ``Tm_min``, ``Tm_max``. The thermodynamic model
            parameters (``Tm_parameters``, ``Tm_chem_correction_parameters``,
            ``Tm_salt_correction_parameters``) are inlined into this dict by :func:`_preprocess_config`.
        :type Tm_filter: dict
        :param secondary_structure_filter: Dict with ``enabled``, ``T``, ``thr_DG``.
        :type secondary_structure_filter: dict
        :return: A filtered `OligoDatabase` object containing only probes that pass all enabled
            property filters. Regions with insufficient oligos after filtering are removed.
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

        # Build sequence-based property filter list, gating each filter on its own ``enabled`` flag.
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

        if secondary_structure_filter["enabled"]:
            secondary_structure = SecondaryStructureFilter(
                T=secondary_structure_filter["T"],
                thr_DG=secondary_structure_filter["thr_DG"],
            )
            filters.append(secondary_structure)

        # initialize the property filter class
        property_filter = PropertyFilter(filters=filters)

        # filter the database — same filters applied to both L and R arms
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
    def filter_by_specificity(
        self,
        oligo_database: OligoDatabase,
        specificity_blastn_filter: dict,
        cross_hybridization_blastn_filter: dict,
    ) -> OligoDatabase:
        """
        Filter the oligo database based on sequence specificity to remove probes that bind
        non-specifically or cross-hybridize.

        The filter list is seeded with an :class:`ExactMatchFilter` (always on) and then conditionally
        extended with the BLASTN-specificity filter depending on its ``enabled`` flag. The
        cross-hybridization filter is built and applied separately because it operates on the
        ``oligo_L`` and ``oligo_R`` sequence types rather than on the joined ``oligo`` sequence.

        1. **Exact matches** (always on): Removes all probes with exact sequence matches to probes
           of other regions.
        2. **BLASTN specificity** (gated on ``specificity_blastn_filter['enabled']``): Uses BLASTN
           to search for similar sequences in the reference database. If ``junction_region_size > 0``,
           all probes where BLASTN hits cover the junction region are removed, independent of the
           coverage threshold.
        3. **Cross-hybridization** (gated on ``cross_hybridization_blastn_filter['enabled']``):
           Removes probes where the left (L) and right (R) halves cross-hybridize with each other.
           This is critical for split probes because if the probes can bind to each other, they may
           form dimers instead of binding to the target RNA. Probes from the larger genomic region
           are removed when cross-hybridization is detected.

        The reference database is loaded from the FASTA file(s) inside
        ``specificity_blastn_filter`` (the only filter that uses it). Regions that do not meet the
        minimum oligo requirement after filtering are removed from the database.

        :param oligo_database: The `OligoDatabase` instance containing oligonucleotide sequences
            and their associated properties.
        :type oligo_database: OligoDatabase
        :param specificity_blastn_filter: Dict with ``enabled``, ``junction_region_size``,
            ``search_parameters``, ``hit_parameters``, ``files_fasta_reference_database``, plus the
            derived ``junction_site`` injected by :func:`_preprocess_config`.
        :type specificity_blastn_filter: dict
        :param cross_hybridization_blastn_filter: Dict with ``enabled``, ``search_parameters``,
            ``hit_parameters``.
        :type cross_hybridization_blastn_filter: dict
        :return: A filtered `OligoDatabase` object containing only probes that pass all enabled
            specificity and cross-hybridization filters. Regions with insufficient oligos after
            filtering are removed.
        :rtype: OligoDatabase
        """
        ##### exact match + specificity filter (exact_matches always on, specificity gated on enabled) #####
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

        ##### cross hybridization filter (gated on enabled) #####
        # Built and applied separately because it operates on the "oligo_L" and "oligo_R"
        # sequence types rather than on the joined "oligo" sequence.
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

        ##### remove all directories of intermediate steps #####
        for directory in directories:
            if os.path.exists(directory):
                shutil.rmtree(directory)

        return oligo_database

    @pipeline_step_basic(step_name="Target Probe Generation - Set Selection")
    def create_oligo_sets(
        self,
        oligo_database: OligoDatabase,
        independent_set_selection: dict,
        isoform_consensus_score: dict,
        Tm_score: dict,
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
        :param independent_set_selection: Dict controlling set generation. Must contain ``n_sets``,
            ``set_size_min``, ``set_size_opt``, ``distance_between_probes``, ``n_attempts_graph``,
            ``n_attempts_clique_enum``, ``diversification_fraction``, ``jaccard_opt``, ``jaccard_step``.
        :type independent_set_selection: dict
        :param isoform_consensus_score: Dict with ``weight``.
        :type isoform_consensus_score: dict
        :param Tm_score: Dict with ``weight`` and ``Tm_opt``. The thermodynamic model parameters
            (``Tm_parameters``, ``Tm_chem_correction_parameters``, ``Tm_salt_correction_parameters``)
            are inlined into this dict by :func:`_preprocess_config`.
        :type Tm_score: dict
        :return: An updated `OligoDatabase` object containing the generated oligo sets. Each region
            will have up to ``n_sets`` sets stored, with each set containing between ``set_size_min``
            and ``set_size_opt`` probes. Regions with insufficient oligos are removed.
        :rtype: OligoDatabase
        """
        # Define all scorers
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
        # the higher the score the better, because we want to have on average oligos with high melting temperatures
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

    def load_codebook(self, file_codebook: str) -> pd.DataFrame:
        return pd.read_csv(file_codebook, sep=None, engine="python", index_col="gene_name")

    def load_readout_probe_table(self, file_readout_probe_table: str) -> pd.DataFrame:
        readout_probe_table = pd.read_csv(file_readout_probe_table, sep=None, engine="python")
        if "bit" not in readout_probe_table.columns:
            readout_probe_table = readout_probe_table.sort_values(by=["readout_probe_id", "channel"])
            readout_probe_table.reset_index(inplace=True, drop=True)
            readout_probe_table["bit"] = "bit_" + (readout_probe_table.index + 1).astype(str)
        return readout_probe_table.set_index("bit")

    def generate_codebook(self, region_ids: list[str], readout_probe_table: pd.DataFrame) -> pd.DataFrame:
        """
        Generate a codebook (barcode matrix) for encoding multiple regions using CycleHCR readout probes.

        This method creates a binary barcode matrix where each row represents a genomic region and
        each column represents a bit position. Each region is assigned a unique barcode consisting
        of exactly two active bits (value 1), representing a left-right (L/R) readout probe pair
        in a specific fluorescence channel.

        The encoding scheme works as follows:
        - The number of channels and the number of L/R probe pairs per channel are derived from
          the ``readout_probe_table`` (``channel`` and ``L/R`` columns).
        - Each barcode is generated from a combination of (probe_L_id, probe_R_id, channel_id).
        - The two active bits correspond to the left and right readout probes in the specified channel.
        - Combinations are prioritized: same probe pairs (L=R) are preferred over different pairs.
        - The codebook size is limited by the number of available probe/channel combinations.

        The method validates that sufficient barcodes are available to encode all requested regions.
        If not, a ``ConfigurationError`` is raised with suggestions to increase the number of probes
        or reduce the number of regions.

        :param region_ids: List of region identifiers (e.g., gene IDs) to encode in the codebook.
            Each region will be assigned a unique barcode.
        :type region_ids: list[str]
        :param readout_probe_table: Bit-indexed table of readout probes (the same table consumed by
            ``assemble_hybridization_probes``). The ``channel`` and ``L/R`` columns are used to
            derive the number of channels and the number of L/R probe pairs per channel.
        :type readout_probe_table: pd.DataFrame
        :return: A pandas DataFrame containing the binary barcode matrix. Rows are indexed by
            ``region_ids``, and columns are named ``bit_1``, ``bit_2``, etc. Only columns with at
            least one active bit are included. Each row has exactly two bits set to 1.
        :rtype: pd.DataFrame
        :raises ConfigurationError: If the number of available barcodes is insufficient to encode
            all requested regions (i.e., ``codebook_size_max < 2 * n_regions``).
        """
        n_channels = len(readout_probe_table["channel"].unique())
        n_readout_probes_R = readout_probe_table["L/R"].value_counts()["R"]
        n_readout_probes_L = readout_probe_table["L/R"].value_counts()["L"]
        n_readout_probes_LR = int(min([n_readout_probes_R, n_readout_probes_L]) / n_channels)

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
        codebook.index.name = "gene_name"

        # Remove columns where all values are 0
        codebook = codebook.loc[:, (codebook != 0).any(axis=0)]

        return codebook

    def generate_readout_probe_table(self) -> pd.DataFrame:
        """
        Generate a CycleHCR readout probe table.

        Placeholder for a future implementation. Once implemented, the output is expected to
        satisfy the same contract as a loaded readout probe table: ``bit`` index, ``channel`` /
        ``readout_probe_id`` / ``L/R`` / ``readout_probe_sequence`` columns.
        """
        raise FeatureNotImplementedError(
            "Generation of readout probe table is not yet implemented. "
            "Please provide a file_readout_probe_table parameter."
        )

    def validate(
        self,
        codebook: pd.DataFrame,
        readout_probe_table: pd.DataFrame,
        region_ids: list[str],
        *,
        codebook_source: str,
        readout_probe_table_source: str,
    ) -> None:
        """
        Validate that a (codebook, readout_probe_table) pair forms a valid CycleHCR readout setup.

        Centralizes the CycleHCR-specific validation contract (two-hot codebook indexed by
        ``gene_name``; bit-indexed readout probe table with ``channel`` / ``readout_probe_id`` /
        ``L/R`` / ``readout_probe_sequence`` columns; codebook bits covered by the table) so that
        all paths producing these tables — loading from file today, generating programmatically in
        the future — share a single validation gate.

        :param codebook: Codebook DataFrame to validate.
        :type codebook: pd.DataFrame
        :param readout_probe_table: Readout probe table DataFrame to validate.
        :type readout_probe_table: pd.DataFrame
        :param region_ids: Region IDs required to be present in the codebook index.
        :type region_ids: list[str]
        :param codebook_source: Source identifier (file path or marker) for the codebook.
        :type codebook_source: str
        :param readout_probe_table_source: Source identifier (file path or marker) for the
            readout probe table.
        :type readout_probe_table_source: str
        :raises FileFormatError: If either input fails validation.
        """
        validate_codebook(
            codebook=codebook,
            region_ids=region_ids,
            source=codebook_source,
            expected_hamming_weight=2,
            index_name="gene_name",
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


def _preprocess_config(config: dict[str, Any]) -> dict[str, Any]:
    """
    Preprocess a CycleHCR pipeline configuration dict in place.

    - Resolves the ``nn_table``/``tmm_table``/``imm_table``/``de_table`` strings in
      ``target_probe.global_parameters.Tm_parameters`` to their ``Bio.SeqUtils.MeltingTemp`` objects.
    - For every Tm chem/salt correction block: if ``enabled`` is ``False`` sets ``parameters`` to
      ``None`` so downstream filters receive a clean ``None``.
    - Inlines Tm parameters and chem/salt corrections into every block that consumes them
      (``Tm_filter`` and ``Tm_score``) so designer methods don't have to thread
      ``global_parameters`` through the call chain.
    - Computes the derived ``oligo_length = L + gap + R`` and ``junction_site = L + gap // 2``
      and injects them into the blocks that consume them (``oligo_generation`` and
      ``specificity_blastn_filter`` respectively).
    - Expands ``target_probe.oligo_generation.file_region_ids`` to a sorted unique list under
      ``target_probe.oligo_generation.region_ids`` (or ``None`` if no file was provided).
    """

    # Preprocess Tm tables and set Tm_chem/salt_correction_parameters to None if the correction is disabled
    for section in ["target_probe"]:
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

    # Inline Tm parameters into Tm_filter (consumed by the property filter and by the final TmNN property calculator)
    config["target_probe"]["property_filters"]["Tm_filter"]["Tm_parameters"] = target_probe_Tm_parameters
    config["target_probe"]["property_filters"]["Tm_filter"][
        "Tm_chem_correction_parameters"
    ] = target_probe_Tm_chem_correction_parameters
    config["target_probe"]["property_filters"]["Tm_filter"][
        "Tm_salt_correction_parameters"
    ] = target_probe_Tm_salt_correction_parameters

    # Inline Tm parameters into Tm_score (consumed by the Tm scorer in create_oligo_sets)
    config["target_probe"]["probe_set_selection"]["Tm_score"]["Tm_parameters"] = target_probe_Tm_parameters
    config["target_probe"]["probe_set_selection"]["Tm_score"][
        "Tm_chem_correction_parameters"
    ] = target_probe_Tm_chem_correction_parameters
    config["target_probe"]["probe_set_selection"]["Tm_score"][
        "Tm_salt_correction_parameters"
    ] = target_probe_Tm_salt_correction_parameters

    # Derive oligo_length and junction_site and inject them into the consuming blocks.
    L_probe_sequence_length = config["target_probe"]["oligo_generation"]["L_probe_sequence_length"]
    gap_sequence_length = config["target_probe"]["oligo_generation"]["gap_sequence_length"]
    R_probe_sequence_length = config["target_probe"]["oligo_generation"]["R_probe_sequence_length"]
    config["target_probe"]["oligo_generation"]["oligo_length"] = (
        L_probe_sequence_length + gap_sequence_length + R_probe_sequence_length
    )
    config["target_probe"]["specificity_filters"]["specificity_blastn_filter"]["junction_site"] = (
        L_probe_sequence_length + gap_sequence_length // 2
    )

    ##### read the genes file #####
    file_region_ids = config["target_probe"]["oligo_generation"]["file_region_ids"]
    if file_region_ids is None:
        logger.warning(
            "No gene list file was provided! All genes from fasta file are used to generate the probes. "
            "This choice can use a lot of resources."
        )
        config["target_probe"]["oligo_generation"]["region_ids"] = None
    else:
        with open(file_region_ids) as f:
            config["target_probe"]["oligo_generation"]["region_ids"] = sorted({line.rstrip() for line in f})

    return config


def cycle_hcr_probe_designer(config: dict[str, Any]) -> None:
    """
    Execute the CycleHCR probe design pipeline from a (raw) configuration dict.

    The dict is expected to follow the nested layout of ``data/configs/cycle_hcr_probe_designer.yaml``
    (``general``, ``target_probe.*``, ``readout_probes``, ``primers``, ``hybridization_probes``).
    The caller is responsible for configuring the library logger before invoking this function
    (see :func:`main`).

    :param config: Pipeline configuration loaded via ``yaml.safe_load``.
    :type config: dict
    """

    ##### preprocess the config file #####
    config_dict = _preprocess_config(config)

    ##### initialize probe designer pipeline #####
    pipeline = CycleHCRProbeDesigner(
        dir_output=config_dict["general"]["dir_output"],
        write_intermediate_steps=config_dict["general"]["write_intermediate_steps"],
        n_jobs=config_dict["general"]["n_jobs"],
    )

    ##### design target probes #####
    target_probe_database = pipeline.design_target_probes(
        oligo_generation_parameters=config_dict["target_probe"]["oligo_generation"],
        property_filters_parameters=config_dict["target_probe"]["property_filters"],
        specificity_filters_parameters=config_dict["target_probe"]["specificity_filters"],
        probe_set_selection_parameters=config_dict["target_probe"]["probe_set_selection"],
    )

    ##### design readout probes (codebook + readout probe table) #####
    codebook, readout_probe_table = pipeline.design_readout_probes(
        region_ids=list(target_probe_database.database.keys()),
        readout_probe_parameters=config_dict["readout_probes"],
    )

    ##### assemble hybridization probes #####
    # Inject runtime-derived dataframes into the parameters dict.
    config_dict["hybridization_probes"]["codebook"] = codebook
    config_dict["hybridization_probes"]["readout_probe_table"] = readout_probe_table
    hybridization_probe_database = pipeline.assemble_hybridization_probes(
        oligo_database=target_probe_database,
        hybridization_probe_parameters=config_dict["hybridization_probes"],
    )

    ##### design primers (load / validate forward + reverse primer sequences) #####
    reverse_primer_sequence, forward_primer_sequence = pipeline.design_primers(
        primer_parameters=config_dict["primers"],
    )

    ##### assemble DNA template probes #####
    # Inject runtime-derived primer sequences into the parameters dict.
    config_dict["hybridization_probes"]["forward_primer_sequence"] = forward_primer_sequence
    config_dict["hybridization_probes"]["reverse_primer_sequence"] = reverse_primer_sequence
    dna_template_probe_database = pipeline.assemble_dna_template_probes(
        oligo_database=hybridization_probe_database,
        hybridization_probe_parameters=config_dict["hybridization_probes"],
    )

    ##### write outputs #####
    pipeline.generate_output(
        oligo_database=dna_template_probe_database,
        codebook=codebook,
        readout_probe_table=readout_probe_table,
    )


def main() -> None:

    print("--------------START PIPELINE--------------")

    args = base_parser()

    ##### read the config file #####
    with open(args["config"], "r") as handle:
        config = yaml.safe_load(handle)

    # setup logger now that we know the output directory
    configure_root_logger(
        dir_output=config["general"]["dir_output"],
        pipeline_name="cyclehcr_probe_designer",
    )

    cycle_hcr_probe_designer(config)

    print("--------------END PIPELINE--------------")


if __name__ == "__main__":
    main()
