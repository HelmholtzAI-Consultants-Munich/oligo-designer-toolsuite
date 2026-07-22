############################################
# imports
############################################

import os
import shutil
from itertools import combinations
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import yaml
from Bio.SeqUtils import Seq
from scipy.spatial.distance import hamming

from oligo_designer_toolsuite._exceptions import (
    ConfigurationError,
    FeatureNotImplementedError,
    FileFormatError,
)
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
    PropertyCalculator,
    ReverseComplementSequenceProperty,
    TmNNProperty,
)
from oligo_designer_toolsuite.oligo_property_calculator._property_functions import calc_tm_nn
from oligo_designer_toolsuite.oligo_property_filter import (
    BasePropertyFilter,
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
from oligo_designer_toolsuite.oligo_selection import (
    HomogeneousPropertyOligoSelection,
    IndependentSetsOligoSelection,
)
from oligo_designer_toolsuite.oligo_specificity_filter import (
    BaseSpecificityFilter,
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
    validate_bit_mapping_table,
    validate_codebook,
    validate_primer_sequence,
)
from oligo_designer_toolsuite.sequence_generator import OligoSequenceGenerator
from oligo_designer_toolsuite.utils import append_nucleotide_to_sequences, configure_root_logger, logger

############################################
# Merfish Probe Designer
############################################


class MerfishProbeDesigner:
    """
    A class for designing hybridization probes for MERFISH (Multiplexed Error-Robust Fluorescence In Situ Hybridization) experiments.

    This class provides a complete pipeline for designing MERFISH probes, which enable multiplexed RNA detection
    in single cells through combinatorial barcoding and sequential imaging rounds.

    **MERFISH Pipeline Overview:**
    1. **Target Probe Design**: Design gene-specific targeting sequences that bind to RNA transcripts
    2. **Readout Probe Design**: Generate readout probe sequences and create a binary codebook for encoding
    3. **Hybridization Probe Assembly**: Combine target probes with readout probes based on the codebook
    4. **Primer Design**: Design PCR primers for amplifying DNA template probes
    5. **DNA Template Probe Assembly**: Assemble final DNA template probes with primers
    6. **Output Generation**: Generate output files in multiple formats (TSV, YAML, Excel)

    Overview
    --------
    MERFISH is an image-based single-cell transcriptomics method that allows
    hundreds to thousands of RNA species to be identified, counted, and spatially
    localized in individual cells while preserving their native context.

    Each RNA species is labeled with a unique binary barcode, which is read out
    through sequential rounds of single-molecule FISH (smFISH) imaging. This
    combinatorial barcoding strategy enables massively multiplexed RNA detection
    with single-molecule precision.

    Probe Structure
    ---------------
    **Hybridization (Encoding) Probes**
    - Single-stranded DNA oligonucleotides that hybridize directly to target RNA transcripts.
    - Each probe contains:
        - A **30-nt targeting sequence** complementary to the target mRNA.
        - Two **20-nt barcode sequences** (readout sequences) that correspond to bits in the RNA’s
        binary barcode and are read out by fluorescently labeled secondary probes.
        - **Single A-nucleotide spacers** separating readout and targeting regions to minimize
        secondary structure formation and nonspecific hybridization.
    - The hybridization probe has the structure:
        [Readout 1] + [Targeting Sequence] + [Readout 2]

    **Readout Probes**
    - Short (typically 20-nt), dye-labeled DNA oligonucleotides that hybridize to the readout
    sequences within the hybridization probes.
    - Each readout probe is complementary to one barcode sequence and carries a fluorophore
    used to report the “on” state of a specific bit during an imaging round.
    - Sequential rounds of hybridization, imaging, and fluorophore cleavage allow decoding
    of the binary barcode for each RNA molecule.
    - The barcoding scheme uses an **error-robust Modified Hamming Distance 4 (MHD4) code**
    with constant Hamming weight (typically four “1” bits). This ensures that each barcode
    differs by at least four bits from all others, enabling detection and correction of
    hybridization or imaging errors during decoding.

    **DNA Template Probe**
    - The hybridization probes are synthesized from oligonucleotide pools containing
    forward and reverse PCR priming regions flanking the probe body.
     - These primer regions enable limited-cycle PCR amplification and in vitro transcription
    to produce RNA intermediates that are later reverse-transcribed into single-stranded
    DNA probes.
    - Primer regions are cleaved after synthesis (e.g., using USER enzyme) to generate
    the final hybridization probe ready for cellular labeling.
    - The **DNA template probe** has the structure:
        [Forward Primer] + [Readout 1] + [Targeting Sequence] + [Readout 2] + [Reverse Primer]

    Probe Library Preparation
    -------------------------
    The MERFISH probe library is generated through a multi-step molecular workflow. First, target
    genes are selected and assigned binary barcodes from an error-robust codebook designed with
    sufficient Hamming distance to enable accurate barcode decoding. For each gene, a set of ~92
    hybridization probes is designed, each containing a 30-nucleotide target-binding region and two
    20-nucleotide readout sequences corresponding to that gene’s barcode. These sequences are
    synthesized as part of a large oligonucleotide pool, flanked by forward and reverse PCR primer
    regions for amplification. The oligo pool is then PCR-amplified under limited cycles,
    transcribed into RNA using a T7 promoter, and reverse-transcribed back into DNA to produce
    single-stranded probes. Following synthesis, the primer regions are enzymatically cleaved
    (e.g., using USER enzyme) to yield the final hybridization probes, approximately 72 nucleotides
    in length. These probes are purified and stored until use.

    References
    ----------
    Wang, G., Moffitt, J. R., & Zhuang, X. (2018).
    "Multiplexed imaging of high-density libraries of RNAs with MERFISH and expansion microscopy."
    *Scientific Reports*, 8, 4847. DOI: 10.1038/s41598-018-22297-7

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
        dir_output: str,
        write_intermediate_steps: bool,
        n_jobs: int,
    ) -> None:
        """Constructor for the MerfishProbeDesigner class."""

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
        Design target probes for MERFISH experiments through a multi-step pipeline.

        This method performs the complete target probe design process, which includes:
        1. Creating an initial oligo database from input FASTA files using a sliding window approach
        2. Filtering probes based on sequence properties (isoform consensus, GC content, Tm,
           homopolymeric runs, secondary structure)
        3. Filtering probes based on specificity to remove off-target binding and cross-hybridization
           using BLASTN searches
        4. Organizing filtered probes into optimal sets based on weighted scoring criteria (isoform
           consensus, GC content, Tm) and distance constraints

        The resulting probes are gene-specific targeting sequences (typically 30 nt) that bind to RNA
        transcripts. These probes will later be combined with readout probe barcodes to create complete
        hybridization probes.

        :param oligo_generation_parameters: ``target_probe.oligo_generation`` block. Contains
            ``region_ids`` (populated from ``file_region_ids`` by :func:`_preprocess_config`),
            ``files_fasta_probe_database``, ``probe_length_min``, ``probe_length_max``.
        :type oligo_generation_parameters: dict
        :param property_filters_parameters: ``target_probe.property_filters`` block. Each filter
            sub-dict carries an ``enabled`` flag plus its parameters; ``Tm_filter`` additionally
            receives the inlined ``Tm_parameters`` / ``Tm_chem_correction_parameters`` /
            ``Tm_salt_correction_parameters`` from :func:`_preprocess_config`.
        :type property_filters_parameters: dict
        :param specificity_filters_parameters: ``target_probe.specificity_filters`` block.
            ``specificity_blastn_filter`` carries ``files_fasta_reference_database`` (shared with
            the cross-hybridization filter).
        :type specificity_filters_parameters: dict
        :param probe_set_selection_parameters: ``target_probe.probe_set_selection`` block. Contains
            the ``independent_set_selection`` scalars and the ``GC_content_score`` / ``Tm_score`` /
            ``isoform_consensus_score`` sub-dicts (Tm parameters + min/max inlined by
            :func:`_preprocess_config`).
        :type probe_set_selection_parameters: dict
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
            GC_content_score=probe_set_selection_parameters["GC_content_score"],
            Tm_score=probe_set_selection_parameters["Tm_score"],
            isoform_consensus_score=probe_set_selection_parameters["isoform_consensus_score"],
        )

        if self.write_intermediate_steps:
            dir_database = oligo_database.save_database(name_database="4_db_target_probes_sets")
            logger.info(
                f"Saved target probe database for step 4 (Set Selection) in directory {dir_database}."
            )

        return oligo_database

    def design_readout_probes(
        self,
        region_ids: list[str],
        readout_probe_parameters: dict,
    ) -> tuple[pd.DataFrame, pd.DataFrame]:
        """
        Design readout probes and generate a codebook for MERFISH experiments through a multi-step pipeline.

        Each artefact (codebook, readout probe table) is loaded from file when
        ``...source == "load"`` or generated programmatically when ``...source == "generate"``.
        The readout probe table always carries all ``n_bits`` rows (the MERFISH standard 16-bit
        readout set), even when the codebook only references a subset — this keeps the full
        readout set available if a caller later loads a different codebook that references other
        bits. Once both are obtained, the pair is validated via
        :py:meth:`ReadoutProbeDesigner.validate` (a warning is emitted for readout bits not
        referenced by the codebook, which is expected).

        The codebook assigns each region a unique binary barcode with a fixed Hamming weight
        (``hamming_weight``) and a minimum Hamming distance (``min_hamming_dist``) that enables
        single-bit error detection/correction. The readout probe table assigns readout probes to
        each bit position and distributes them across fluorescence channels for multiplexed
        detection.

        :param region_ids: List of region identifiers (e.g., gene IDs) for which readout probes and
            codebook entries are to be generated. The number of regions determines the minimum number
            of barcodes required in the codebook.
        :type region_ids: list[str]
        :param readout_probe_parameters: ``readout_probes`` block. Must contain ``codebook`` (with
            ``source``, optional ``file``, ``n_bits``, ``min_hamming_dist``, ``hamming_weight``)
            and ``readout_probe_table`` (with ``source``, optional ``file``, ``channels_ids``, and
            the multi-step generation parameters — ``oligo_generation``, ``property_filters``,
            ``specificity_filters``, ``probe_set_selection``, ``global_parameters`` — when
            ``source == "generate"``).
        :type readout_probe_parameters: dict
        :return: A tuple containing:
            - **codebook** (pd.DataFrame): Binary barcode matrix with ``gene_name`` as index and bit
              columns (bit_1, bit_2, etc.) as data. Each row represents a region's barcode, with
              exactly ``hamming_weight`` bits set to 1.
            - **readout_probe_table** (pd.DataFrame): Table mapping each bit referenced by the
              codebook to a readout probe sequence, channel assignment, and probe ID. Indexed by
              bit labels (bit_1, bit_2, etc.).
        :rtype: tuple[pd.DataFrame, pd.DataFrame]
        """
        readout_probe_designer = ReadoutProbeDesigner(
            dir_output=self.dir_output,
            n_jobs=self.n_jobs,
        )

        ##### codebook: load or generate #####
        if readout_probe_parameters["codebook"]["source"] == "load":
            codebook = readout_probe_designer.load_codebook(
                file_codebook=readout_probe_parameters["codebook"]["file"]
            )
            codebook_source = readout_probe_parameters["codebook"]["file"]
        else:
            codebook = readout_probe_designer.generate_codebook(
                region_ids=region_ids,
                codebook_parameters=readout_probe_parameters["codebook"],
            )
            codebook_source = readout_probe_parameters["codebook"]["source"]

        ##### readout probe table: load or generate via the multi-step pipeline #####
        if readout_probe_parameters["readout_probe_table"]["source"] == "load":
            readout_probe_table = readout_probe_designer.load_readout_probe_table(
                file_readout_probe_table=readout_probe_parameters["readout_probe_table"]["file"]
            )
            readout_probe_table_source = readout_probe_parameters["readout_probe_table"]["file"]
        else:
            readout_probe_table = readout_probe_designer.generate_readout_probe_table(
                readout_probe_parameters=readout_probe_parameters["readout_probe_table"],
                codebook_parameters=readout_probe_parameters["codebook"],
                write_intermediate_steps=self.write_intermediate_steps,
            )
            readout_probe_table_source = readout_probe_parameters["readout_probe_table"]["source"]

        readout_probe_designer.validate(
            codebook=codebook,
            readout_probe_table=readout_probe_table,
            region_ids=region_ids,
            codebook_source=codebook_source,
            readout_probe_table_source=readout_probe_table_source,
            hamming_weight=readout_probe_parameters["codebook"]["hamming_weight"],
        )

        return codebook, readout_probe_table

    def assemble_hybridization_probes(
        self,
        oligo_database: OligoDatabase,
        codebook: pd.DataFrame,
        readout_probe_table: pd.DataFrame,
    ) -> OligoDatabase:
        """
        Assemble hybridization probes by combining target probes with readout probe sequences based on the codebook.

        This method creates complete MERFISH hybridization probes by combining gene-specific target probes
        with readout probe barcodes according to the codebook assignment. For each region, the method:
        1. Looks up the region's barcode in the codebook to identify which two readout probes are assigned
        2. Retrieves the corresponding readout probe sequences from the readout probe table
        3. Assembles the hybridization probe sequence with the structure:
           [reverse_complement(readout_probe_1)] + "A" + [target_probe] + "A" + [reverse_complement(readout_probe_2)]

        The readout probes are reverse-complemented because they will hybridize to the barcode sequences
        embedded in the hybridization probe. The single "A" nucleotides serve as spacers between the
        readout probe binding sites and the target probe sequence.

        The assembled sequences are stored as properties in the database for each probe, enabling downstream
        primer addition and DNA template probe assembly.

        :param oligo_database: The target-probe `OligoDatabase` (returned by :py:meth:`design_target_probes`).
            The method assembles hybridization probes onto this database and returns it.
        :type oligo_database: OligoDatabase
        :param codebook: A pandas DataFrame containing binary barcodes for each region. Rows are indexed
            by ``gene_name``, and columns represent bit positions (bit_1, bit_2, etc.). Each row has exactly
            ``hamming_weight`` bits set to 1, indicating which readout probes are assigned to that region.
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
            - `sequence_hybridization_probe`: The complete assembled hybridization probe sequence
        :rtype: OligoDatabase
        """
        region_ids = list(oligo_database.database.keys())

        oligo_database.set_database_sequence_types(
            [
                "sequence_target",
                "sequence_readout_probe_1",
                "sequence_readout_probe_2",
                "sequence_hybridization_probe",
            ]
        )

        for region_id in region_ids:
            barcode = codebook.loc[region_id]
            bits = barcode[barcode == 1].index
            readout_probe_sequences = readout_probe_table.loc[bits, "readout_probe_sequence"]
            sequence_readout_probe_1 = readout_probe_sequences.iloc[0]
            sequence_readout_probe_2 = readout_probe_sequences.iloc[1]

            probe_ids = list(oligo_database.database[region_id].keys())
            new_properties: dict[str, dict[str, str]] = {probe_id: {} for probe_id in probe_ids}

            for probe_id in probe_ids:

                new_properties[probe_id]["sequence_target"] = format_sequence(
                    database=oligo_database,
                    property="target",
                    region_id=region_id,
                    oligo_id=probe_id,
                )

                new_properties[probe_id]["sequence_readout_probe_1"] = sequence_readout_probe_1
                new_properties[probe_id]["sequence_readout_probe_2"] = sequence_readout_probe_2

                new_properties[probe_id]["sequence_hybridization_probe"] = (
                    str(Seq(sequence_readout_probe_1).reverse_complement())
                    + "A"
                    + format_sequence(
                        database=oligo_database,
                        property="oligo",
                        region_id=region_id,
                        oligo_id=probe_id,
                    )
                    + "A"
                    + str(Seq(sequence_readout_probe_2).reverse_complement())
                )

            oligo_database.update_oligo_properties(new_properties)

        return oligo_database

    def design_primers(
        self,
        oligo_database: OligoDatabase,
        primer_parameters: dict,
    ) -> tuple[str, str]:
        """
        Design forward and reverse primers for MERFISH hybridization probes through a multi-step pipeline.

        Each primer (forward, reverse) is loaded from its config block when
        ``...source == "load"`` or generated programmatically when ``...source == "generate"``.
        :py:meth:`PrimerDesigner.generate_forward_primer` runs the full multi-step primer pipeline
        and selects the candidate whose Tm is closest to the reverse primer's Tm;
        :py:meth:`PrimerDesigner.generate_reverse_primer` is currently a placeholder. After both
        primers are obtained the pair is validated via :py:meth:`PrimerDesigner.validate`.

        :param oligo_database: The hybridization-probe `OligoDatabase` (returned by
            :py:meth:`assemble_hybridization_probes`). Written to a temporary FASTA so it can be
            used as a reference for the second specificity-filter pass (the primer must not bind
            the hybridization probes themselves).
        :type oligo_database: OligoDatabase
        :param primer_parameters: ``primers`` block. Must contain ``forward_primer`` (with
            ``source`` and either ``sequence`` when ``source == "load"`` or the full
            ``oligo_generation`` / ``property_filters`` / ``specificity_filters`` /
            ``global_parameters`` sub-blocks when ``source == "generate"`` — the Tm parameters
            inside ``global_parameters`` are inlined into ``property_filters.Tm_filter`` by
            :func:`_preprocess_config`) and ``reverse_primer`` (with ``source`` and ``sequence``
            when ``source == "load"``).
        :type primer_parameters: dict
        :return: A tuple ``(reverse_primer_sequence, forward_primer_sequence)``.
        :rtype: tuple[str, str]
        """
        # Write the hybridization probes to a FASTA file so they can be used as a reference for
        # the second specificity pass (primer must not bind the hybridization probes themselves).
        file_fasta_hybridization_probes_database = oligo_database.write_database_to_fasta(
            filename="db_reference_hybridization_probes",
            save_description=False,
            region_ids=None,
            sequence_type="sequence_hybridization_probe",
        )

        primer_designer = PrimerDesigner(
            dir_output=self.dir_output,
            n_jobs=self.n_jobs,
        )

        ##### reverse primer: load or generate (placeholder) #####
        if primer_parameters["reverse_primer"]["source"] == "load":
            reverse_primer_sequence = primer_designer.load_reverse_primer(
                sequence=primer_parameters["reverse_primer"]["sequence"]
            )
        else:
            reverse_primer_sequence = primer_designer.generate_reverse_primer()

        ##### forward primer: load or generate via the multi-step pipeline #####
        if primer_parameters["forward_primer"]["source"] == "load":
            forward_primer_sequence = primer_designer.load_forward_primer(
                sequence=primer_parameters["forward_primer"]["sequence"]
            )
        else:
            forward_primer_sequence = primer_designer.generate_forward_primer(
                parameters=primer_parameters["forward_primer"],
                reverse_primer_sequence=reverse_primer_sequence,
                file_fasta_hybridization_probes_database=file_fasta_hybridization_probes_database,
                write_intermediate_steps=self.write_intermediate_steps,
            )

        os.remove(file_fasta_hybridization_probes_database)

        primer_designer.validate(
            forward_primer=forward_primer_sequence,
            reverse_primer=reverse_primer_sequence,
        )

        return reverse_primer_sequence, forward_primer_sequence

    def assemble_dna_template_probes(
        self,
        oligo_database: OligoDatabase,
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

        :param oligo_database: The hybridization-probe `OligoDatabase` (returned by
            :py:meth:`assemble_hybridization_probes`) — must expose the
            ``sequence_hybridization_probe`` property on each probe. The method attaches DNA
            template probe sequences to this database and returns it.
        :type oligo_database: OligoDatabase
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
        region_ids = list(oligo_database.database.keys())
        oligo_database.set_database_sequence_types(
            [
                "sequence_forward_primer",
                "sequence_reverse_primer",
                "sequence_dna_template_probe",
            ]
        )

        for region_id in region_ids:
            probe_ids = list(oligo_database.database[region_id].keys())
            new_properties: dict[str, dict[str, str]] = {probe_id: {} for probe_id in probe_ids}

            for probe_id in probe_ids:
                new_properties[probe_id]["sequence_reverse_primer"] = reverse_primer_sequence
                new_properties[probe_id]["sequence_forward_primer"] = forward_primer_sequence

                new_properties[probe_id]["sequence_dna_template_probe"] = (
                    forward_primer_sequence
                    + format_sequence(
                        database=oligo_database,
                        property="sequence_hybridization_probe",
                        region_id=region_id,
                        oligo_id=probe_id,
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
        Generate the final output files for the MERFISH probe design pipeline.

        This method writes all output files required for the MERFISH experiment, including codebooks,
        readout probe tables, and probe sequences in multiple formats. The output files are written
        to the pipeline's output directory.

        **Generated Output Files:**

        1. **codebook.tsv**: Binary barcode matrix with region IDs as index and bit columns. Each row
           represents a region's barcode assignment.

        2. **readout_probes.tsv**: Table mapping readout probe sequences to bit positions and channels.
           Contains columns: bit, channel, readout_probe_id, readout_probe_sequence.

        3. **merfish_probes.yml**: Complete probe information in YAML format, including all specified
           properties for each probe in the top N sets per region.

        4. **merfish_probes.tsv**: Complete probe information in TSV format, including all specified
           properties for each probe in the top N sets per region.

        5. **merfish_probes.xlsx**: Complete probe information in Excel format with one sheet per region.
           Each sheet contains probe sets for that region with all specified properties.

        6. **merfish_probes_order.yml**: Simplified YAML file containing only the essential sequences
           needed for ordering probes (DNA template probe and readout probe sequences).

        :param oligo_database: The DNA-template-probe `OligoDatabase` (returned by
            :py:meth:`assemble_dna_template_probes`) — the final probe database written to the
            output files.
        :type oligo_database: OligoDatabase
        :param codebook: A pandas DataFrame containing binary barcodes for each region. Rows are indexed
            by region IDs, and columns represent bit positions. This should be the codebook generated
            by the `design_readout_probes` method.
        :type codebook: pd.DataFrame
        :param readout_probe_table: A pandas DataFrame containing readout probe sequences and their
            associated bit identifiers and channel assignments. This should be the readout probe table
            generated by the `design_readout_probes` method.
        :type readout_probe_table: pd.DataFrame
        :param output_properties: List of property names to include in the output files. If None, a default
            set of properties will be included. Available properties include: 'source', 'species', 'gene_id',
            'chromosome', 'start', 'end', 'strand', 'sequence_target', 'sequence_readout_probe_1',
            'sequence_readout_probe_2', 'sequence_hybridization_probe', 'sequence_forward_primer',
            'sequence_reverse_primer', 'sequence_dna_template_probe', 'isoform_consensus', etc.
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
                "sequence_target",
                "sequence_readout_probe_1",
                "sequence_readout_probe_2",
                "sequence_hybridization_probe",
                "sequence_forward_primer",
                "sequence_reverse_primer",
                "sequence_dna_template_probe",
                "isoform_consensus",
            ]

        codebook.to_csv(os.path.join(self.dir_output, "codebook.tsv"), sep="\t", index_label="gene_name")
        readout_probe_table.to_csv(os.path.join(self.dir_output, "readout_probes.tsv"), sep="\t")

        oligo_database.write_oligosets_to_yaml(
            properties=output_properties,
            ascending=True,
            filename="merfish_probes",
        )

        oligo_database.write_oligosets_to_table(
            properties=output_properties,
            ascending=True,
            filename="merfish_probes",
        )

        oligo_database.write_ready_to_order_yaml(
            properties=[
                "sequence_dna_template_probe",
                "sequence_readout_probe_1",
                "sequence_readout_probe_2",
            ],
            ascending=True,
            filename="merfish_probes_order",
        )


############################################
# Merfish Target Probe Designer
############################################


class TargetProbeDesigner:
    """
    A class for designing target probes for MERFISH experiments through a multi-step pipeline.

    This class provides methods for the complete target probe design process, which includes:
    1. Creating an initial oligo database from input FASTA files using a sliding window approach
    2. Filtering probes based on sequence properties (GC content, melting temperature, homopolymeric
       runs, secondary structure)
    3. Filtering probes based on specificity to remove off-target binding and cross-hybridization
       using BLASTN searches
    4. Organizing filtered probes into optimal sets based on weighted scoring criteria (isoform
       consensus, GC content, melting temperature) and distance constraints

    The resulting probes are gene-specific targeting sequences (typically 30 nt) that bind to RNA
    transcripts. These probes will later be combined with readout probe barcodes to create complete
    hybridization probes.

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
        region_ids: list[str] | None,
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
        3. Calculates the reverse complement of each sequence

        The database stores sequences with sequence types "target" (original sequence) and
        "oligo" (reverse complement). Isoform-consensus filtering is applied downstream in
        :py:meth:`filter_by_property`.

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
            component sequences (target, oligo). The database is filtered to only include regions
            that meet the minimum oligo requirement.
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
        oligo_database.set_database_sequence_types(["target", "oligo"])
        oligo_database.load_database_from_fasta(
            files_fasta=oligo_fasta_file,
            database_overwrite=True,
            sequence_type="target",
            region_ids=region_ids,
        )

        ##### compute reverse complement (always on) #####
        reverse_complement_sequence_property = ReverseComplementSequenceProperty(
            sequence_type_reverse_complement="oligo"
        )
        calculator = PropertyCalculator(properties=[reverse_complement_sequence_property])
        oligo_database = calculator.apply(
            oligo_database=oligo_database, sequence_type="target", n_jobs=self.n_jobs
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

        This method applies multiple property-based filters, each gated on its own ``enabled`` flag.
        The isoform-consensus filter is applied first as a cheap pre-filter on the ``target``
        sequence type; the remaining sequence-based filters are then applied to the ``oligo``
        sequence type. Probes that fail any enabled filter are removed from the database.

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

        # initialize the property filter class
        property_filter = PropertyFilter(filters=filters)

        # filter the database
        oligo_database = property_filter.apply(
            oligo_database=oligo_database,
            sequence_type="oligo",
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
        extended with BLASTN-specificity and cross-hybridization filters depending on their
        ``enabled`` flags. All filters are applied in a single :class:`SpecificityFilter` invocation
        so the database is iterated once.

        1. **Exact matches** (always on): Removes all probes with exact sequence matches to probes of
           other regions.
        2. **BLASTN specificity** (gated on ``specificity_blastn_filter['enabled']``): Uses BLASTN to
           search for similar sequences in the reference database. Probes with hits meeting the
           specified criteria are removed.
        3. **Cross-hybridization** (gated on ``cross_hybridization_blastn_filter['enabled']``):
           Removes probes that cross-hybridize with each other. This is critical because if probes
           can bind to each other, they may form dimers instead of binding to the target RNA. Probes
           from the larger genomic region are removed when cross-hybridization is detected.

        The reference database is loaded from the FASTA file(s) inside ``specificity_blastn_filter``
        and is shared by both BLASTN-based filters. Regions that do not meet the minimum oligo
        requirement after filtering are removed from the database.

        :param oligo_database: The `OligoDatabase` instance containing oligonucleotide sequences
            and their associated properties.
        :type oligo_database: OligoDatabase
        :param specificity_blastn_filter: Dict with ``enabled``, ``search_parameters``,
            ``hit_parameters``, ``files_fasta_reference_database``.
        :type specificity_blastn_filter: dict
        :param cross_hybridization_blastn_filter: Dict with ``enabled``, ``search_parameters``,
            ``hit_parameters``.
        :type cross_hybridization_blastn_filter: dict
        :return: A filtered `OligoDatabase` object containing only probes that pass all enabled
            specificity and cross-hybridization filters. Regions with insufficient oligos after
            filtering are removed.
        :rtype: OligoDatabase
        """
        ##### exact match filter (always on); BLASTN filters gated on ``enabled`` #####
        exact_matches = ExactMatchFilter(policy=RemoveAllFilterPolicy(), filter_name="exact_match")
        filters: list[BaseSpecificityFilter] = [exact_matches]
        directories: list[str] = []

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

    @pipeline_step_basic(step_name="Target Probe Generation - Set Selection")
    def create_oligo_sets(
        self,
        oligo_database: OligoDatabase,
        independent_set_selection: dict,
        GC_content_score: dict,
        Tm_score: dict,
        isoform_consensus_score: dict,
    ) -> OligoDatabase:
        """
        Create optimal oligo sets based on weighted scoring criteria, distance constraints, and set selection.

        This method selects optimal sets of target probes for each region by:
        1. Scoring each oligo based on weighted criteria (isoform consensus, GC content, Tm) — the
           MERFISH pipeline uses ``Normalized`` deviation-from-optimum scorers, so each score dict
           needs both the ``opt`` target and the surrounding ``min``/``max`` bounds (inlined by
           :func:`_preprocess_config` from the corresponding filter blocks).
        2. Building a graph where edges represent non-overlapping oligos (based on distance constraints)
        3. Selecting sets of oligos that minimize the average score while respecting distance constraints
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
        :param GC_content_score: Dict with ``weight``, ``GC_content_opt`` and inlined
            ``GC_content_min`` / ``GC_content_max`` (from ``GC_content_filter``).
        :type GC_content_score: dict
        :param Tm_score: Dict with ``weight``, ``Tm_opt``, inlined ``Tm_min`` / ``Tm_max`` (from
            ``Tm_filter``), and the inlined thermodynamic-model parameters (``Tm_parameters``,
            ``Tm_chem_correction_parameters``, ``Tm_salt_correction_parameters``).
        :type Tm_score: dict
        :param isoform_consensus_score: Dict with ``weight``.
        :type isoform_consensus_score: dict
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
# Merfish Readout Probe Designer
############################################


class ReadoutProbeDesigner:
    """
    A class for designing MERFISH readout probes and generating codebooks through a multi-step pipeline.

    This class provides methods for the complete readout probe design process, which includes:
    1. Creating an initial oligo database by generating random sequences with specified nucleotide
       probabilities
    2. Filtering probes based on sequence properties (GC content, homopolymeric runs)
    3. Filtering probes based on specificity to remove off-target binding and cross-hybridization
       using BLASTN searches
    4. Organizing filtered probes into sets with homogeneous properties (GC content and melting
       temperature) for consistent hybridization behavior
    5. Generating a binary barcode codebook with Hamming distance constraints for error correction
    6. Creating a readout probe table that maps codebook bits to channels and readout probe sequences

    The resulting readout probes are non-targeting sequences that bind to barcode sequences on the
    hybridization probes. Each readout probe is assigned to a specific bit position in the codebook
    and distributed across fluorescence channels for multiplexed detection.

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

    def load_codebook(self, file_codebook: str) -> pd.DataFrame:
        """Load a MERFISH codebook from a TSV/CSV file (index column: ``gene_name``)."""
        return pd.read_csv(file_codebook, sep=None, engine="python", index_col="gene_name")

    def generate_codebook(
        self,
        region_ids: list[str],
        codebook_parameters: dict,
    ) -> pd.DataFrame:
        """
        Generate a binary barcode codebook with specified Hamming distance constraints and fixed Hamming weight.

        The codebook assigns each region a unique binary barcode where:
        - each barcode has a fixed Hamming weight (``hamming_weight`` active bits),
        - the minimum Hamming distance between any two barcodes meets ``min_hamming_dist``,
        - the number of valid barcodes is sufficient to encode all ``region_ids``.

        The Hamming distance constraint provides error correction capability: if a barcode is misread
        during imaging it can still be correctly identified as long as the number of bit errors is less
        than half the minimum Hamming distance.

        The algorithm generates all possible barcodes with the specified Hamming weight, then filters
        them to ensure the minimum Hamming distance constraint is met. The returned codebook always
        carries all ``n_bits`` columns — bit positions that no accepted codeword uses appear as
        all-zero columns rather than being dropped, so the codebook shape matches the ``n_bits``
        readout probe set 1:1.

        :param region_ids: List of region identifiers (e.g., gene IDs) to encode in the codebook.
            Each region will be assigned a unique barcode.
        :type region_ids: list[str]
        :param codebook_parameters: ``readout_probes.codebook`` block. Must contain ``n_bits``,
            ``min_hamming_dist``, ``hamming_weight``.
        :type codebook_parameters: dict
        :return: A DataFrame containing the codebook with binary encoded bits. Rows are indexed by
            ``gene_name`` (from ``region_ids``); columns are ``bit_1`` .. ``bit_{n_bits}``. Bits not
            used by any codeword appear as all-zero columns.
        :rtype: pd.DataFrame
        :raises ConfigurationError: If the number of valid barcodes (meeting Hamming distance constraints)
            is insufficient for the number of regions. In this case, consider increasing ``n_bits``,
            reducing ``min_hamming_dist``, or reducing the number of regions.
        """

        def _generate_barcode(raw_barcode: list, n_bits: int) -> np.ndarray:
            barcode = np.zeros(n_bits, dtype=np.int8)
            for i in raw_barcode:
                barcode[i] = 1
            return barcode

        n_bits = codebook_parameters["n_bits"]
        min_hamming_dist = codebook_parameters["min_hamming_dist"]
        hamming_weight = codebook_parameters["hamming_weight"]

        n_regions = len(region_ids)
        codebook_list: list[np.ndarray] = []
        for raw_barcode in combinations(iterable=range(n_bits), r=hamming_weight):
            new_barcode = _generate_barcode(raw_barcode=list(raw_barcode), n_bits=n_bits)
            # check if the barcode passes the requirements
            add_new_barcode = True
            for barcode in codebook_list:
                hamming_dist = hamming(new_barcode, barcode) * n_bits
                if hamming_dist < min_hamming_dist:
                    add_new_barcode = False
                    break
            if add_new_barcode:
                codebook_list.append(new_barcode)
        if len(codebook_list) < n_regions:
            raise ConfigurationError(
                f"The number of valid barcodes ({len(codebook_list)}) is lower than the number of regions ({n_regions}). "
                f"Consider increasing the number of bits or reducing the number of regions."
            )

        codebook: pd.DataFrame = pd.DataFrame(
            codebook_list[0:n_regions], index=region_ids, columns=[f"bit_{i+1}" for i in range(n_bits)]
        )
        codebook.index.name = "gene_name"

        return codebook

    def load_readout_probe_table(self, file_readout_probe_table: str) -> pd.DataFrame:
        """Load a MERFISH readout probe table from a TSV/CSV file. The ``bit`` column must be present."""
        readout_probe_table = pd.read_csv(file_readout_probe_table, sep=None, engine="python")
        if "bit" not in readout_probe_table.columns:
            raise FileFormatError(
                f"Readout probe table '{file_readout_probe_table}' must contain a 'bit' column."
            )
        return readout_probe_table.set_index("bit")

    def generate_readout_probe_table(
        self,
        readout_probe_parameters: dict,
        codebook_parameters: dict,
        write_intermediate_steps: bool = False,
    ) -> pd.DataFrame:
        """
        Generate a MERFISH readout probe table by running the full multi-step readout probe
        design pipeline and assigning the surviving sequences to bit / channel slots.

        Internally orchestrates the existing decorated steps:
        :py:meth:`_create_oligo_database` → :py:meth:`_filter_by_property` →
        :py:meth:`_filter_by_specificity` → :py:meth:`_create_oligo_sets`, then formats the
        selected set into a bit-indexed table via :py:meth:`_format_readout_probe_table`. Each
        decorated step keeps its own ``@pipeline_step_basic`` logging.

        :param readout_probe_parameters: ``readout_probes.readout_probe_table`` block. Must contain
            ``oligo_generation``, ``property_filters``, ``specificity_filters``,
            ``probe_set_selection`` (with the inlined ``Tm_parameters`` /
            ``Tm_chem_correction_parameters`` / ``Tm_salt_correction_parameters`` from
            :func:`_preprocess_config`) and ``channels_ids``.
        :type readout_probe_parameters: dict
        :param codebook_parameters: ``readout_probes.codebook`` block; ``n_bits`` is used to size
            the returned table.
        :type codebook_parameters: dict
        :param write_intermediate_steps: If True, save the per-step readout-probe databases for
            debugging.
        :type write_intermediate_steps: bool
        :return: Bit-indexed readout probe table with ``channel``, ``readout_probe_id``,
            ``readout_probe_sequence`` columns.
        :rtype: pd.DataFrame
        """
        oligo_database: OligoDatabase = self._create_oligo_database(
            oligo_length=readout_probe_parameters["oligo_generation"]["probe_length"],
            oligo_base_probabilities=readout_probe_parameters["oligo_generation"]["base_probabilities"],
            initial_num_sequences=readout_probe_parameters["oligo_generation"]["initial_num_sequences"],
        )

        if write_intermediate_steps:
            dir_database = oligo_database.save_database(name_database="1_db_readout_probes_initial")
            logger.info(
                f"Saved readout probe database for step 1 (Create Database) in directory {dir_database}"
            )

        oligo_database = self._filter_by_property(
            oligo_database=oligo_database,
            GC_content_filter=readout_probe_parameters["property_filters"]["GC_content_filter"],
            homopolymeric_runs_filter=readout_probe_parameters["property_filters"][
                "homopolymeric_runs_filter"
            ],
        )

        if write_intermediate_steps:
            dir_database = oligo_database.save_database(name_database="2_db_readout_probes_property_filter")
            logger.info(
                f"Saved readout probe database for step 2 (Property Filters) in directory {dir_database}"
            )

        oligo_database = self._filter_by_specificity(
            oligo_database=oligo_database,
            specificity_blastn_filter=readout_probe_parameters["specificity_filters"][
                "specificity_blastn_filter"
            ],
            cross_hybridization_blastn_filter=readout_probe_parameters["specificity_filters"][
                "cross_hybridization_blastn_filter"
            ],
        )

        if write_intermediate_steps:
            dir_database = oligo_database.save_database(name_database="3_db_readout_probes_specificty_filter")
            logger.info(
                f"Saved readout probe database for step 3 (Specificity Filters) in directory {dir_database}"
            )

        oligo_database = self._create_oligo_sets(
            oligo_database=oligo_database,
            set_size=readout_probe_parameters["probe_set_selection"]["set_size"],
            homogeneous_properties_weights=readout_probe_parameters["probe_set_selection"][
                "homogeneous_properties_weights"
            ],
            n_combinations=readout_probe_parameters["probe_set_selection"]["n_combinations"],
            Tm_parameters=readout_probe_parameters["probe_set_selection"]["Tm_parameters"],
            Tm_chem_correction_parameters=readout_probe_parameters["probe_set_selection"][
                "Tm_chem_correction_parameters"
            ],
            Tm_salt_correction_parameters=readout_probe_parameters["probe_set_selection"][
                "Tm_salt_correction_parameters"
            ],
        )

        if write_intermediate_steps:
            dir_database = oligo_database.save_database(name_database="4_db_readout_probes_set_selection")
            logger.info(
                f"Saved readout probe database for step 4 (Set Selection) in directory {dir_database}"
            )

        return self._format_readout_probe_table(
            oligo_database=oligo_database,
            channels_ids=readout_probe_parameters["channels_ids"],
            n_bits=codebook_parameters["n_bits"],
        )

    def validate(
        self,
        codebook: pd.DataFrame,
        readout_probe_table: pd.DataFrame,
        region_ids: list[str],
        *,
        codebook_source: str,
        readout_probe_table_source: str,
        hamming_weight: int,
    ) -> None:
        """
        Validate that a (codebook, readout_probe_table) pair forms a valid MERFISH readout setup.

        Delegates to the shared :func:`validate_codebook` and :func:`validate_bit_mapping_table`
        helpers with MERFISH-specific knobs (codebook indexed by ``gene_name``; readout probe
        table with the MERFISH three-column layout; every codeword has exactly ``hamming_weight``
        active bits).
        """
        validate_codebook(
            codebook=codebook,
            region_ids=region_ids,
            source=codebook_source,
            expected_hamming_weight=hamming_weight,
            index_name="gene_name",
        )
        validate_bit_mapping_table(
            table=readout_probe_table,
            codebook=codebook,
            source=readout_probe_table_source,
            required_columns=["channel", "readout_probe_id", "readout_probe_sequence"],
            sequence_columns=["readout_probe_sequence"],
        )

    @pipeline_step_basic(step_name="Readout Probe Generation - Create Oligo Database")
    def _create_oligo_database(
        self,
        oligo_length: int,
        oligo_base_probabilities: dict[str, float],
        initial_num_sequences: int,
    ) -> OligoDatabase:
        """
        Create an initial oligo database by generating random sequences with specified nucleotide probabilities.

        Private helper of :py:meth:`generate_readout_probe_table`. Generates random DNA sequences
        with user-defined nucleotide probabilities and creates an `OligoDatabase` to store them.
        These sequences are then filtered and organized into sets by downstream steps.

        :param oligo_length: Length (in nucleotides) of each readout probe sequence to generate.
        :type oligo_length: int
        :param oligo_base_probabilities: Dictionary specifying the probability of each nucleotide
            base in randomly generated sequences. Keys should be 'A', 'T', 'G', 'C' and values should
            sum to 1.0 (e.g., {"A": 0.25, "T": 0.25, "G": 0.25, "C": 0.25}).
        :type oligo_base_probabilities: dict[str, float]
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
        oligo_database.set_database_sequence_types(["oligo"])
        oligo_database.load_database_from_fasta(
            files_fasta=oligo_fasta_file,
            database_overwrite=True,
            sequence_type="oligo",
            region_ids=None,
        )

        dir = oligo_sequences.dir_output
        shutil.rmtree(dir) if os.path.exists(dir) else None

        return oligo_database

    @pipeline_step_basic(step_name="Readout Probe Generation - Property Filters")
    def _filter_by_property(
        self,
        oligo_database: OligoDatabase,
        GC_content_filter: dict,
        homopolymeric_runs_filter: dict,
    ) -> OligoDatabase:
        """
        Filter the oligo database based on sequence properties.

        Each filter is gated on its own ``enabled`` flag. Probes that fail any enabled filter are
        removed from the database. Private helper of :py:meth:`generate_readout_probe_table`.

        :param oligo_database: The `OligoDatabase` instance containing oligonucleotide sequences
            and their associated properties. This database should contain readout probe sequences
            generated in the previous step.
        :type oligo_database: OligoDatabase
        :param GC_content_filter: Dict with ``enabled``, ``GC_content_min``, ``GC_content_max``.
        :type GC_content_filter: dict
        :param homopolymeric_runs_filter: Dict with ``enabled``, ``homopolymeric_base_n`` (mapping
            ``A``/``T``/``C``/``G`` to maximum allowed run lengths).
        :type homopolymeric_runs_filter: dict
        :return: A filtered `OligoDatabase` object containing only probes that pass all enabled
            property filters. Regions with insufficient oligos after filtering are removed.
        :rtype: OligoDatabase
        """
        # Build property-filter list, gating each filter on its own ``enabled`` flag.
        filters: list[BasePropertyFilter] = []

        if GC_content_filter["enabled"]:
            gc_content = GCContentFilter(
                GC_content_min=GC_content_filter["GC_content_min"],
                GC_content_max=GC_content_filter["GC_content_max"],
            )
            filters.append(gc_content)

        if homopolymeric_runs_filter["enabled"]:
            homopolymeric_runs = HomopolymericRunsFilter(
                base_n=homopolymeric_runs_filter["homopolymeric_base_n"],
            )
            filters.append(homopolymeric_runs)

        # initialize the property filter class
        property_filter = PropertyFilter(filters=filters)

        # filter the database
        oligo_database = property_filter.apply(
            oligo_database=oligo_database,
            sequence_type="oligo",
            n_jobs=self.n_jobs,
        )
        check_content_oligo_database(oligo_database)

        return oligo_database

    @pipeline_step_basic(step_name="Readout Probe Generation - Specificity Filters")
    def _filter_by_specificity(
        self,
        oligo_database: OligoDatabase,
        specificity_blastn_filter: dict,
        cross_hybridization_blastn_filter: dict,
    ) -> OligoDatabase:
        """
        Filter the oligo database based on sequence specificity.

        The filter list is seeded with an :class:`ExactMatchFilter` (always on) and then conditionally
        extended with BLASTN-specificity and cross-hybridization filters depending on their
        ``enabled`` flags. The reference database is loaded from the FASTA file(s) inside
        ``specificity_blastn_filter`` and is shared by both BLASTN-based filters. Private helper
        of :py:meth:`generate_readout_probe_table`.

        :param oligo_database: The `OligoDatabase` instance containing oligonucleotide sequences
            and their associated properties. This database should contain readout probe sequences
            that have passed property filtering.
        :type oligo_database: OligoDatabase
        :param specificity_blastn_filter: Dict with ``enabled``, ``search_parameters``,
            ``hit_parameters``, ``files_fasta_reference_database``.
        :type specificity_blastn_filter: dict
        :param cross_hybridization_blastn_filter: Dict with ``enabled``, ``search_parameters``,
            ``hit_parameters``.
        :type cross_hybridization_blastn_filter: dict
        :return: A filtered `OligoDatabase` object containing only probes that pass all enabled
            specificity and cross-hybridization filters. Regions with insufficient oligos after
            filtering are removed.
        :rtype: OligoDatabase
        """
        ##### exact match filter (always on); BLASTN filters gated on ``enabled`` #####
        exact_matches = ExactMatchFilter(
            policy=RemoveAllFilterPolicy(), filter_name="readout_probes_exact_match"
        )
        filters: list[BaseSpecificityFilter] = [exact_matches]
        directories: list[str] = []

        if specificity_blastn_filter["enabled"]:
            reference_database = ReferenceDatabase(
                database_name=self.subdir_db_reference, dir_output=self.dir_output
            )
            reference_database.load_database_from_file(
                files=specificity_blastn_filter["files_fasta_reference_database"],
                file_type="fasta",
                database_overwrite=False,
            )
            specificity = BlastNFilter(
                search_parameters=specificity_blastn_filter["search_parameters"],
                hit_parameters=specificity_blastn_filter["hit_parameters"],
                filter_name="readout_probes_blastn_specificity",
                dir_output=self.dir_output,
            )
            specificity.set_reference_database(reference_database=reference_database)
            filters.append(specificity)
            directories.append(specificity.dir_output)

        if cross_hybridization_blastn_filter["enabled"]:
            cross_hybridization_aligner = BlastNFilter(
                search_parameters=cross_hybridization_blastn_filter["search_parameters"],
                hit_parameters=cross_hybridization_blastn_filter["hit_parameters"],
                filter_name="readout_probes_blastn_crosshybridization",
                dir_output=self.dir_output,
            )
            cross_hybridization = CrossHybridizationFilter(
                policy=RemoveByDegreeFilterPolicy(),
                alignment_method=cross_hybridization_aligner,
                filter_name="readout_probes_blastn_crosshybridization",
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

        # remove all directories of intermediate steps
        for directory in directories:
            if os.path.exists(directory):
                shutil.rmtree(directory)

        check_content_oligo_database(oligo_database)
        return oligo_database

    @pipeline_step_basic(step_name="Readout Probe Generation - Set Selection")
    def _create_oligo_sets(
        self,
        oligo_database: OligoDatabase,
        set_size: int,
        homogeneous_properties_weights: dict[str, float],
        n_combinations: int,
        Tm_parameters: dict,
        Tm_chem_correction_parameters: dict | None,
        Tm_salt_correction_parameters: dict | None,
    ) -> OligoDatabase:
        """
        Create optimal oligo sets with homogeneous properties (GC content and melting temperature).

        Private helper of :py:meth:`generate_readout_probe_table`. Steps:
        1. **Property calculation**: Calculates GC content and melting temperature (Tm) for each oligo
           using the specified Tm parameters and corrections
        2. **Set generation**: Organizes oligos into sets that have homogeneous properties. The algorithm
           selects sets where all probes have similar GC content and Tm values, which is important for
           consistent hybridization behavior across the readout probe set
        3. **Region filtering**: Removes regions that cannot generate sets meeting the size requirement

        Uses :class:`HomogeneousPropertyOligoSelection` to evaluate multiple combinations and select
        the best set with the most uniform properties.

        :param oligo_database: The `OligoDatabase` instance containing oligonucleotide sequences
            and their associated properties. This database should contain readout probe sequences
            that have passed all previous filtering steps.
        :type oligo_database: OligoDatabase
        :param set_size: Number of readout probes to include in each set. Sets are selected to have
            homogeneous properties (GC content and melting temperature).
        :type set_size: int
        :param homogeneous_properties_weights: Dictionary specifying weights for property homogeneity
            in set selection. Keys should be property names (e.g., 'GC_content', 'TmNN') and values
            are weights that determine the relative importance of each property in the homogeneity score.
        :type homogeneous_properties_weights: dict[str, float]
        :param n_combinations: Number of combinations to evaluate during set creation. Higher values
            may find better sets but increase computation time.
        :type n_combinations: int
        :param Tm_parameters: Dictionary of parameters for calculating melting temperature (Tm) of readout
            probes using the nearest-neighbor method. For using Bio.SeqUtils.MeltingTemp default parameters, set to ``{}``.
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
        :return: An updated `OligoDatabase` object containing the generated oligo sets. Each set contains
            `set_size` probes with homogeneous properties. Regions with insufficient oligos are removed.
        :rtype: OligoDatabase
        """
        # Calculate Tm and GC content
        properties = [
            TmNNProperty(
                Tm_parameters=Tm_parameters,
                Tm_chem_correction_parameters=Tm_chem_correction_parameters,
                Tm_salt_correction_parameters=Tm_salt_correction_parameters,
            ),
            GCContentProperty(),
        ]
        calculator = PropertyCalculator(properties=properties)
        oligo_database = calculator.apply(
            oligo_database=oligo_database, sequence_type="oligo", n_jobs=self.n_jobs
        )

        set_generator = HomogeneousPropertyOligoSelection(
            set_size=set_size,
            properties=homogeneous_properties_weights,
            n_combinations=n_combinations,
        )
        oligo_database = set_generator.apply(
            oligo_database=oligo_database,
            sequence_type="oligo",
            n_sets=1,
            n_jobs=self.n_jobs,
        )

        oligo_database.remove_regions_with_insufficient_oligos(pipeline_step="Oligo Selection")
        check_content_oligo_database(oligo_database)

        return oligo_database

    def _format_readout_probe_table(
        self, oligo_database: OligoDatabase, channels_ids: list[str], n_bits: int
    ) -> pd.DataFrame:
        """
        Format a filtered readout-probe database into a bit-indexed table mapping each bit to a
        channel and readout probe sequence.

        Private helper of :py:meth:`generate_readout_probe_table`. Steps:
        1. A readout probe sequence from the database
        2. A fluorescence channel identifier

        The readout probes are distributed across channels in a round-robin fashion. For example,
        if there are 3 channels (e.g., ['Cy3', 'Cy5', 'Alexa488']), bit_1 gets channel 0, bit_2 gets
        channel 1, bit_3 gets channel 2, bit_4 gets channel 0 again, and so on.

        The table is indexed by bit labels (bit_1, bit_2, etc.) and contains columns for channel,
        readout probe ID, and readout probe sequence. This table is used to assign readout probes
        to each bit position in the codebook for multiplexed imaging.

        :param oligo_database: Filtered ``OligoDatabase`` of readout probe sequences (returned by
            :py:meth:`_create_oligo_sets`).
        :type oligo_database: OligoDatabase
        :param channels_ids: List of fluorescence channel identifiers (e.g., ['Cy3', 'Cy5', 'Alexa488'])
            to which readout probes will be assigned. Readout probes are distributed across channels
            in a round-robin fashion.
        :type channels_ids: list[str]
        :param n_bits: Number of bits in the codebook, representing the number of readout probes needed.
            This should match the number of bit columns in the codebook generated by `generate_codebook`.
        :type n_bits: int
        :return: A DataFrame containing the readout probe table with columns:
            - **channel**: Fluorescence channel identifier assigned to this bit
            - **readout_probe_id**: Unique identifier for the readout probe
            - **readout_probe_sequence**: DNA sequence of the readout probe
            The DataFrame is indexed by bit labels (bit_1, bit_2, etc.).
        :rtype: pd.DataFrame
        :raises AssertionError: If the number of available readout probes in the database is less
            than the number of required bits (`n_bits`). In this case, generate more readout probes
            or reduce the number of bits in the codebook.
        """
        readout_probes = oligo_database.get_oligoid_sequence_mapping(
            sequence_type="oligo", sequence_to_upper=False
        )
        assert (
            len(readout_probes) >= n_bits
        ), f"There are less readout probes ({len(readout_probes)}) than bits ({n_bits})."
        readout_probe_table = pd.DataFrame(
            columns=["bit", "channel", "readout_probe_id", "readout_probe_sequence"],
            index=list(range(n_bits)),
        )
        n_channels = len(channels_ids)
        channel = 0
        for i, (readout_probe_id, readout_probe_sequence) in enumerate(readout_probes.items()):
            readout_probe_table.iloc[i] = pd.Series(
                [
                    f"bit_{i+1}",
                    channels_ids[channel],
                    readout_probe_id,
                    readout_probe_sequence,
                ],
                index=[
                    "bit",
                    "channel",
                    "readout_probe_id",
                    "readout_probe_sequence",
                ],
            )
            channel = (channel + 1) % n_channels
            if i >= n_bits - 1:
                break
        readout_probe_table.set_index("bit", inplace=True)
        return readout_probe_table


############################################
# Merfish Primer Designer
############################################


class PrimerDesigner:
    """
    A class for designing MERFISH primers through a multi-step pipeline.

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

    def load_reverse_primer(self, sequence: str) -> str:
        """Load and validate a reverse primer sequence (whitespace stripped)."""
        return str(sequence).strip()

    def generate_reverse_primer(self) -> str:
        """
        Generate a MERFISH reverse primer.

        Placeholder for a future implementation. MERFISH currently expects the reverse primer
        sequence to be provided via the config (``primers.reverse_primer.source = "load"``); the
        forward primer is selected to match the reverse primer's melting temperature.
        """
        raise FeatureNotImplementedError(
            "Generation of the reverse primer is not yet implemented. "
            "Please provide a reverse_primer.sequence parameter and set reverse_primer.source to 'load'."
        )

    def load_forward_primer(self, sequence: str) -> str:
        """Load and validate a forward primer sequence (whitespace stripped)."""
        return str(sequence).strip()

    def generate_forward_primer(
        self,
        parameters: dict,
        reverse_primer_sequence: str,
        file_fasta_hybridization_probes_database: str,
        write_intermediate_steps: bool = False,
    ) -> str:
        """
        Generate a MERFISH forward primer by running the full multi-step primer design pipeline
        and selecting the candidate whose melting temperature is closest to the reverse primer's Tm.

        Internally orchestrates the existing decorated steps:
        :py:meth:`_create_oligo_database` → :py:meth:`_filter_by_property` →
        :py:meth:`_filter_by_specificity`, then delegates the best-Tm match selection to
        :py:meth:`_get_best_forward_primer`. Each decorated step keeps its own
        ``@pipeline_step_basic`` logging.

        :param parameters: ``primers.forward_primer`` block. Must contain ``oligo_generation``,
            ``property_filters``, ``specificity_filters`` sub-blocks. Tm parameters are expected to
            have been inlined into ``property_filters.Tm_filter`` by :func:`_preprocess_config`
            (Tm tables resolved, Tm chem/salt ``parameters`` normalized to ``None`` when disabled);
            the best-Tm matcher reads them from there.
        :type parameters: dict
        :param reverse_primer_sequence: Reverse primer sequence whose Tm the selected forward primer
            should match.
        :type reverse_primer_sequence: str
        :param file_fasta_hybridization_probes_database: Path to a FASTA file with the assembled
            hybridization probe sequences (created at runtime by :py:meth:`design_primers`). Used as
            reference for the second specificity-filter pass (primer must not bind the hybridization
            probes themselves).
        :type file_fasta_hybridization_probes_database: str
        :param write_intermediate_steps: If True, save the per-step primer databases for debugging.
        :type write_intermediate_steps: bool
        :return: Selected forward primer sequence (Tm closest to the reverse primer).
        :rtype: str
        """
        oligo_database: OligoDatabase = self._create_oligo_database(
            oligo_length=parameters["oligo_generation"]["probe_length"],
            oligo_base_probabilities=parameters["oligo_generation"]["base_probabilities"],
            initial_num_sequences=parameters["oligo_generation"]["initial_num_sequences"],
        )

        if write_intermediate_steps:
            dir_database = oligo_database.save_database(name_database="1_db_primers_initial")
            logger.info(f"Saved primer database for step 1 (Create Database) in directory {dir_database}")

        oligo_database = self._filter_by_property(
            oligo_database=oligo_database,
            GC_content_filter=parameters["property_filters"]["GC_content_filter"],
            GC_clamp_filter=parameters["property_filters"]["GC_clamp_filter"],
            homopolymeric_runs_filter=parameters["property_filters"]["homopolymeric_runs_filter"],
            self_complementarity_filter=parameters["property_filters"]["self_complementarity_filter"],
            complement_reverse_primer_filter=parameters["property_filters"][
                "complement_reverse_primer_filter"
            ],
            Tm_filter=parameters["property_filters"]["Tm_filter"],
            secondary_structure_filter=parameters["property_filters"]["secondary_structure_filter"],
            reverse_primer_sequence=reverse_primer_sequence,
        )

        if write_intermediate_steps:
            dir_database = oligo_database.save_database(name_database="2_db_primer_property_filter")
            logger.info(f"Saved primer database for step 2 (Property Filters) in directory {dir_database}")

        oligo_database = self._filter_by_specificity(
            oligo_database=oligo_database,
            specificity_blastn_filter=parameters["specificity_filters"]["specificity_blastn_filter"],
            hybridization_probes_blastn_filter=parameters["specificity_filters"][
                "hybridization_probes_blastn_filter"
            ],
            file_fasta_hybridization_probes_database=file_fasta_hybridization_probes_database,
        )

        if write_intermediate_steps:
            dir_database = oligo_database.save_database(name_database="3_db_primer_specificty_filter")
            logger.info(f"Saved primer database for step 3 (Specificity Filters) in directory {dir_database}")

        forward_primer_sequence = self._get_best_forward_primer(
            oligo_database=oligo_database,
            reverse_primer_sequence=reverse_primer_sequence,
            Tm_parameters=parameters["property_filters"]["Tm_filter"]["Tm_parameters"],
            Tm_chem_correction_parameters=parameters["property_filters"]["Tm_filter"][
                "Tm_chem_correction_parameters"
            ],
            Tm_salt_correction_parameters=parameters["property_filters"]["Tm_filter"][
                "Tm_salt_correction_parameters"
            ],
        )

        return forward_primer_sequence

    def validate(self, forward_primer: str, reverse_primer: str) -> None:
        """
        Validate a (forward_primer, reverse_primer) pair as valid DNA sequences.

        Delegates to the shared :func:`validate_primer_sequence` helper for each primer; the
        source identifiers are hardcoded as ``"forward_primer"`` / ``"reverse_primer"`` because
        primer roles are fixed.
        """
        validate_primer_sequence(sequence=forward_primer, source="forward_primer")
        validate_primer_sequence(sequence=reverse_primer, source="reverse_primer")

    @pipeline_step_basic(step_name="Primer Generation - Create Oligo Database")
    def _create_oligo_database(
        self,
        oligo_length: int,
        oligo_base_probabilities: dict[str, float],
        initial_num_sequences: int,
    ) -> OligoDatabase:
        """
        Create an initial oligo database by generating random sequences with specified nucleotide
        probabilities, all ending with a "T" nucleotide.

        Private helper of :py:meth:`generate_forward_primer`. Generates random DNA sequences with
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
        :type oligo_base_probabilities: dict[str, float]
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
    def _filter_by_property(
        self,
        oligo_database: OligoDatabase,
        GC_content_filter: dict,
        GC_clamp_filter: dict,
        homopolymeric_runs_filter: dict,
        self_complementarity_filter: dict,
        complement_reverse_primer_filter: dict,
        Tm_filter: dict,
        secondary_structure_filter: dict,
        reverse_primer_sequence: str,
    ) -> OligoDatabase:
        """
        Filter the oligo database based on sequence properties.

        Each filter is gated on its own ``enabled`` flag. Probes that fail any enabled filter are
        removed from the database. Private helper of :py:meth:`generate_forward_primer`.

        :param oligo_database: The `OligoDatabase` instance containing oligonucleotide sequences
            and their associated properties. This database should contain primer sequences generated
            in the previous step.
        :type oligo_database: OligoDatabase
        :param GC_content_filter: Dict with ``enabled``, ``GC_content_min``, ``GC_content_max``.
        :type GC_content_filter: dict
        :param GC_clamp_filter: Dict with ``enabled``, ``number_GC_GCclamp``,
            ``number_three_prime_base_GCclamp``.
        :type GC_clamp_filter: dict
        :param homopolymeric_runs_filter: Dict with ``enabled``, ``homopolymeric_base_n``.
        :type homopolymeric_runs_filter: dict
        :param self_complementarity_filter: Dict with ``enabled``, ``max_len_selfcomplement``.
        :type self_complementarity_filter: dict
        :param complement_reverse_primer_filter: Dict with ``enabled``, ``max_len_complement``.
        :type complement_reverse_primer_filter: dict
        :param Tm_filter: Dict with ``enabled``, ``Tm_min``, ``Tm_max``, ``Tm_parameters``,
            ``Tm_chem_correction_parameters``, ``Tm_salt_correction_parameters`` (Tm parameters
            are inlined into this dict by :func:`_preprocess_config`).
        :type Tm_filter: dict
        :param secondary_structure_filter: Dict with ``enabled``, ``T``, ``thr_DG``.
        :type secondary_structure_filter: dict
        :param reverse_primer_sequence: DNA sequence of the reverse primer used by the
            complement-reverse-primer filter to prevent primer-dimer formation.
        :type reverse_primer_sequence: str
        :return: A filtered `OligoDatabase` object containing only primers that pass all enabled
            property filters. Regions with insufficient oligos after filtering are removed.
        :rtype: OligoDatabase
        """
        # Build property-filter list, gating each filter on its own ``enabled`` flag.
        filters: list[BasePropertyFilter] = []

        if GC_content_filter["enabled"]:
            gc_content = GCContentFilter(
                GC_content_min=GC_content_filter["GC_content_min"],
                GC_content_max=GC_content_filter["GC_content_max"],
            )
            filters.append(gc_content)

        if GC_clamp_filter["enabled"]:
            gc_clamp = GCClampFilter(
                n_bases=GC_clamp_filter["number_three_prime_base_GCclamp"],
                n_GC=GC_clamp_filter["number_GC_GCclamp"],
            )
            filters.append(gc_clamp)

        if homopolymeric_runs_filter["enabled"]:
            homopolymeric_runs = HomopolymericRunsFilter(
                base_n=homopolymeric_runs_filter["homopolymeric_base_n"],
            )
            filters.append(homopolymeric_runs)

        if self_complementarity_filter["enabled"]:
            self_complement = SelfComplementFilter(
                max_len_selfcomplement=self_complementarity_filter["max_len_selfcomplement"],
            )
            filters.append(self_complement)

        if complement_reverse_primer_filter["enabled"]:
            complement = ComplementFilter(
                comparison_sequence=reverse_primer_sequence,
                max_len_complement=complement_reverse_primer_filter["max_len_complement"],
            )
            filters.append(complement)

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

        # filter the database
        oligo_database = property_filter.apply(
            oligo_database=oligo_database,
            sequence_type="oligo",
            n_jobs=self.n_jobs,
        )
        check_content_oligo_database(oligo_database)

        return oligo_database

    @pipeline_step_basic(step_name="Primer Generation - Specificity Filters")
    def _filter_by_specificity(
        self,
        oligo_database: OligoDatabase,
        specificity_blastn_filter: dict,
        hybridization_probes_blastn_filter: dict,
        file_fasta_hybridization_probes_database: str,
    ) -> OligoDatabase:
        """
        Filter the oligo database based on sequence specificity.

        Each BLASTN filter is gated on its own ``enabled`` flag. Primers with hits meeting the
        specified criteria are removed. Private helper of :py:meth:`generate_forward_primer`.

        :param oligo_database: The `OligoDatabase` instance containing oligonucleotide sequences
            and their associated properties. This database should contain primer sequences that have
            passed property filtering.
        :type oligo_database: OligoDatabase
        :param specificity_blastn_filter: Dict with ``enabled``, ``search_parameters``,
            ``hit_parameters``, ``files_fasta_reference_database``.
        :type specificity_blastn_filter: dict
        :param hybridization_probes_blastn_filter: Dict with ``enabled``, ``search_parameters``,
            ``hit_parameters``.
        :type hybridization_probes_blastn_filter: dict
        :param file_fasta_hybridization_probes_database: Path to the FASTA file containing the
            assembled hybridization probe sequences (created at runtime). Used as the reference
            database for the hybridization-probe BLASTN filter so primers do not bind the
            hybridization probes themselves.
        :type file_fasta_hybridization_probes_database: str
        :return: A filtered `OligoDatabase` object containing only primers that pass all enabled
            specificity filters. Regions with insufficient oligos after filtering are removed.
        :rtype: OligoDatabase
        """
        ##### BLASTN filters gated on ``enabled`` #####
        filters: list[BaseSpecificityFilter] = []
        directories: list[str] = []

        if specificity_blastn_filter["enabled"]:
            reference_database = ReferenceDatabase(
                database_name=self.subdir_db_reference, dir_output=self.dir_output
            )
            reference_database.load_database_from_file(
                files=specificity_blastn_filter["files_fasta_reference_database"],
                file_type="fasta",
                database_overwrite=True,
            )
            specificity_refrence = BlastNFilter(
                search_parameters=specificity_blastn_filter["search_parameters"],
                hit_parameters=specificity_blastn_filter["hit_parameters"],
                filter_name="primer_blastn_specificity_reference",
                dir_output=self.dir_output,
            )
            specificity_refrence.set_reference_database(reference_database=reference_database)
            filters.append(specificity_refrence)
            directories.append(specificity_refrence.dir_output)

        if hybridization_probes_blastn_filter["enabled"]:
            hybridization_probes_database = ReferenceDatabase(
                database_name=self.subdir_db_reference, dir_output=self.dir_output
            )
            hybridization_probes_database.load_database_from_file(
                files=file_fasta_hybridization_probes_database,
                file_type="fasta",
                database_overwrite=True,
            )
            specificity_hybridization_probes = BlastNFilter(
                search_parameters=hybridization_probes_blastn_filter["search_parameters"],
                hit_parameters=hybridization_probes_blastn_filter["hit_parameters"],
                filter_name="primer_blastn_specificity_hybridization_probes",
                dir_output=self.dir_output,
            )
            specificity_hybridization_probes.set_reference_database(
                reference_database=hybridization_probes_database
            )
            filters.append(specificity_hybridization_probes)
            directories.append(specificity_hybridization_probes.dir_output)

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

    def _get_best_forward_primer(
        self,
        oligo_database: OligoDatabase,
        reverse_primer_sequence: str,
        Tm_parameters: dict,
        Tm_chem_correction_parameters: dict | None,
        Tm_salt_correction_parameters: dict | None,
    ) -> str:
        """
        Pick the candidate primer in ``oligo_database`` whose melting temperature is closest to
        the reverse primer's Tm, for balanced PCR amplification.

        :param oligo_database: Filtered ``OligoDatabase`` of forward primer candidates.
        :type oligo_database: OligoDatabase
        :param reverse_primer_sequence: Reverse primer sequence whose Tm the selected forward
            primer should match.
        :type reverse_primer_sequence: str
        :param Tm_parameters: Nearest-neighbor Tm parameters.
        :type Tm_parameters: dict
        :param Tm_chem_correction_parameters: Tm chemistry correction parameters (``None`` if disabled).
        :type Tm_chem_correction_parameters: dict | None
        :param Tm_salt_correction_parameters: Tm salt correction parameters (``None`` if disabled).
        :type Tm_salt_correction_parameters: dict | None
        :return: Forward primer sequence with the Tm closest to the reverse primer's.
        :rtype: str
        """
        # Compute Tm of the reverse primer once.
        Tm_reverse_primer = calc_tm_nn(
            sequence=reverse_primer_sequence,
            Tm_parameters=Tm_parameters,
            Tm_chem_correction_parameters=Tm_chem_correction_parameters,
            Tm_salt_correction_parameters=Tm_salt_correction_parameters,
        )

        # Iterate over surviving primer candidates and pick the one whose Tm is closest to the reverse primer's.
        min_dif_Tm = float("inf")
        forward_primer_sequence = ""
        for database_region in oligo_database.database.values():
            for primer_properties in database_region.values():
                Tm_forward_primer = calc_tm_nn(
                    sequence=primer_properties["oligo"],
                    Tm_parameters=Tm_parameters,
                    Tm_chem_correction_parameters=Tm_chem_correction_parameters,
                    Tm_salt_correction_parameters=Tm_salt_correction_parameters,
                )
                dif_Tm = abs(Tm_forward_primer - Tm_reverse_primer)
                if dif_Tm < min_dif_Tm:
                    min_dif_Tm = dif_Tm
                    forward_primer_sequence = primer_properties["oligo"]

        return forward_primer_sequence


############################################
# Merfish Probe Designer Pipeline
############################################


def _preprocess_config(config: dict[str, Any]) -> dict[str, Any]:
    """
    Preprocess a MERFISH pipeline configuration dict in place.

    - Resolves the ``nn_table`` / ``tmm_table`` / ``imm_table`` / ``de_table`` strings in every
      ``global_parameters.Tm_parameters`` block (target_probe, readout_probes.readout_probe_table,
      primers.forward_primer when ``source == "generate"``) to their
      ``Bio.SeqUtils.MeltingTemp`` objects.
    - For every Tm chem/salt correction block, sets ``parameters`` to ``None`` when ``enabled`` is
      ``False`` so downstream filters receive a clean ``None``.
    - Inlines Tm parameters and chem/salt corrections into every block that consumes them
      (``target_probe.property_filters.Tm_filter``, ``target_probe.probe_set_selection.Tm_score``,
      ``readout_probes.readout_probe_table.probe_set_selection``, and
      ``primers.forward_primer.property_filters.Tm_filter``) so designer methods don't have to
      thread ``global_parameters`` through the call chain.
    - Copies ``Tm_min`` / ``Tm_max`` from ``target_probe.property_filters.Tm_filter`` into
      ``target_probe.probe_set_selection.Tm_score`` (the ``NormalizedDeviationFromOptimalTmScorer``
      needs these bounds); same for ``GC_content_min`` / ``GC_content_max`` into
      ``GC_content_score``.
    - Expands ``target_probe.oligo_generation.file_region_ids`` to a sorted unique list under
      ``target_probe.oligo_generation.region_ids`` (or ``None`` if no file was provided).
    """

    ##### Tm preprocessing for target_probe #####
    config["target_probe"]["global_parameters"]["Tm_parameters"] = preprocess_tm_parameters(
        config["target_probe"]["global_parameters"]["Tm_parameters"]
    )
    for correction in ["Tm_chem_correction_parameters", "Tm_salt_correction_parameters"]:
        correction_cfg = config["target_probe"]["global_parameters"][correction]
        if not correction_cfg["enabled"]:
            correction_cfg["parameters"] = None

    target_probe_Tm_parameters = config["target_probe"]["global_parameters"]["Tm_parameters"]
    target_probe_Tm_chem_correction_parameters = config["target_probe"]["global_parameters"][
        "Tm_chem_correction_parameters"
    ]["parameters"]
    target_probe_Tm_salt_correction_parameters = config["target_probe"]["global_parameters"][
        "Tm_salt_correction_parameters"
    ]["parameters"]

    # Inline Tm parameters into Tm_filter (consumed by the property filter).
    config["target_probe"]["property_filters"]["Tm_filter"]["Tm_parameters"] = target_probe_Tm_parameters
    config["target_probe"]["property_filters"]["Tm_filter"][
        "Tm_chem_correction_parameters"
    ] = target_probe_Tm_chem_correction_parameters
    config["target_probe"]["property_filters"]["Tm_filter"][
        "Tm_salt_correction_parameters"
    ] = target_probe_Tm_salt_correction_parameters

    # Inline Tm parameters + Tm_min/max into Tm_score (consumed by NormalizedDeviationFromOptimalTmScorer).
    config["target_probe"]["probe_set_selection"]["Tm_score"]["Tm_min"] = config["target_probe"][
        "property_filters"
    ]["Tm_filter"]["Tm_min"]
    config["target_probe"]["probe_set_selection"]["Tm_score"]["Tm_max"] = config["target_probe"][
        "property_filters"
    ]["Tm_filter"]["Tm_max"]
    config["target_probe"]["probe_set_selection"]["Tm_score"]["Tm_parameters"] = target_probe_Tm_parameters
    config["target_probe"]["probe_set_selection"]["Tm_score"][
        "Tm_chem_correction_parameters"
    ] = target_probe_Tm_chem_correction_parameters
    config["target_probe"]["probe_set_selection"]["Tm_score"][
        "Tm_salt_correction_parameters"
    ] = target_probe_Tm_salt_correction_parameters

    # Inline GC_content_min/max into GC_content_score (needed by NormalizedDeviationFromOptimalGCContentScorer).
    config["target_probe"]["probe_set_selection"]["GC_content_score"]["GC_content_min"] = config[
        "target_probe"
    ]["property_filters"]["GC_content_filter"]["GC_content_min"]
    config["target_probe"]["probe_set_selection"]["GC_content_score"]["GC_content_max"] = config[
        "target_probe"
    ]["property_filters"]["GC_content_filter"]["GC_content_max"]

    ##### Tm preprocessing for readout_probes.readout_probe_table (used by the SET-selection Tm scorer) #####
    if config["readout_probes"]["readout_probe_table"]["source"] == "generate":
        config["readout_probes"]["readout_probe_table"]["global_parameters"]["Tm_parameters"] = (
            preprocess_tm_parameters(
                config["readout_probes"]["readout_probe_table"]["global_parameters"]["Tm_parameters"]
            )
        )
        for correction in ["Tm_chem_correction_parameters", "Tm_salt_correction_parameters"]:
            correction_cfg = config["readout_probes"]["readout_probe_table"]["global_parameters"][correction]
            if not correction_cfg["enabled"]:
                correction_cfg["parameters"] = None

        # Inline Tm parameters into probe_set_selection (consumed by the homogeneous-property set scorer).
        config["readout_probes"]["readout_probe_table"]["probe_set_selection"]["Tm_parameters"] = config[
            "readout_probes"
        ]["readout_probe_table"]["global_parameters"]["Tm_parameters"]
        config["readout_probes"]["readout_probe_table"]["probe_set_selection"][
            "Tm_chem_correction_parameters"
        ] = config["readout_probes"]["readout_probe_table"]["global_parameters"][
            "Tm_chem_correction_parameters"
        ][
            "parameters"
        ]
        config["readout_probes"]["readout_probe_table"]["probe_set_selection"][
            "Tm_salt_correction_parameters"
        ] = config["readout_probes"]["readout_probe_table"]["global_parameters"][
            "Tm_salt_correction_parameters"
        ][
            "parameters"
        ]

    ##### Tm preprocessing for primers.forward_primer #####
    if config["primers"]["forward_primer"]["source"] == "generate":
        config["primers"]["forward_primer"]["global_parameters"]["Tm_parameters"] = preprocess_tm_parameters(
            config["primers"]["forward_primer"]["global_parameters"]["Tm_parameters"]
        )
        for correction in ["Tm_chem_correction_parameters", "Tm_salt_correction_parameters"]:
            correction_cfg = config["primers"]["forward_primer"]["global_parameters"][correction]
            if not correction_cfg["enabled"]:
                correction_cfg["parameters"] = None

        primer_Tm_parameters = config["primers"]["forward_primer"]["global_parameters"]["Tm_parameters"]
        primer_Tm_chem_correction_parameters = config["primers"]["forward_primer"]["global_parameters"][
            "Tm_chem_correction_parameters"
        ]["parameters"]
        primer_Tm_salt_correction_parameters = config["primers"]["forward_primer"]["global_parameters"][
            "Tm_salt_correction_parameters"
        ]["parameters"]

        # Inline Tm parameters into Tm_filter (consumed by the property filter and by the best-Tm matcher).
        config["primers"]["forward_primer"]["property_filters"]["Tm_filter"][
            "Tm_parameters"
        ] = primer_Tm_parameters
        config["primers"]["forward_primer"]["property_filters"]["Tm_filter"][
            "Tm_chem_correction_parameters"
        ] = primer_Tm_chem_correction_parameters
        config["primers"]["forward_primer"]["property_filters"]["Tm_filter"][
            "Tm_salt_correction_parameters"
        ] = primer_Tm_salt_correction_parameters

    ##### region ids #####
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


def merfish_probe_designer(config: dict[str, Any]) -> None:
    """
    Execute the MERFISH probe design pipeline from a (raw) configuration dict.

    The dict is expected to follow the nested layout of ``data/configs/merfish_probe_designer.yaml``
    (``general``, ``target_probe.*``, ``readout_probes``, ``primers``). The caller is responsible
    for configuring the library logger before invoking this function (see :func:`main`).

    :param config: Pipeline configuration loaded via ``yaml.safe_load``.
    :type config: dict
    """

    ##### preprocess the config file #####
    config_dict = _preprocess_config(config)

    ##### initialize probe designer pipeline #####
    pipeline = MerfishProbeDesigner(
        write_intermediate_steps=config_dict["general"]["write_intermediate_steps"],
        dir_output=config_dict["general"]["dir_output"],
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
    hybridization_probe_database = pipeline.assemble_hybridization_probes(
        oligo_database=target_probe_database,
        codebook=codebook,
        readout_probe_table=readout_probe_table,
    )

    ##### design primers (load / generate forward + reverse primer sequences) #####
    reverse_primer_sequence, forward_primer_sequence = pipeline.design_primers(
        oligo_database=hybridization_probe_database,
        primer_parameters=config_dict["primers"],
    )

    ##### assemble DNA template probes #####
    dna_template_probe_database = pipeline.assemble_dna_template_probes(
        oligo_database=hybridization_probe_database,
        reverse_primer_sequence=reverse_primer_sequence,
        forward_primer_sequence=forward_primer_sequence,
    )

    ##### write outputs #####
    pipeline.generate_output(
        oligo_database=dna_template_probe_database,
        codebook=codebook,
        readout_probe_table=readout_probe_table,
    )


def main() -> None:
    """
    Main entry point for running the MERFISH probe design pipeline.

    Parses command-line arguments, loads the YAML config, configures the logger, and delegates
    to :func:`merfish_probe_designer`. The function is typically called from the command line:
    ``merfish_probe_designer --config <path_to_config.yaml>``.
    """
    print("--------------START PIPELINE--------------")

    args = base_parser()

    ##### read the config file #####
    with open(args["config"], "r") as handle:
        config = yaml.safe_load(handle)

    # setup logger now that we know the output directory
    configure_root_logger(
        dir_output=config["general"]["dir_output"],
        pipeline_name="merfish_probe_designer",
    )

    merfish_probe_designer(config)

    print("--------------END PIPELINE--------------")


if __name__ == "__main__":
    main()
